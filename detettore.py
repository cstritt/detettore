#!/usr/bin/env python

""" dettetore - a program to detect and characterize transposable element polymorphisms

Copyright (C) 2021 C. Stritt
License: GNU General Public License v3. See LICENSE.txt for details.
        
"""

import argparse
import logging
import os
import shutil
import time
import gzip
import sys
import warnings

from scripts import strumenti


def get_args():

    parser = argparse.ArgumentParser(
        description='Detect TE polymorphisms from paired-end read data')

    # Parameter groups
    parser_input = parser.add_argument_group('INPUT / OUTPUT')

    parser_settings = parser.add_argument_group('PROGRAM SETTINGS')

    parser_thresholds = parser.add_argument_group('ALIGNMENT AND MAPPING THRESHOLDS')


    # REQUIRED
    parser_input.add_argument(
        '-b', dest="bamfile",
        required=True,
        help='Paired-end reads aligned to reference genome.')

    parser_input.add_argument(
        "-r", dest="reference",
        required=True,
        help='Reference genome in fasta format. Required to calculate read statistics')

    parser_input.add_argument(
        '-t', dest="targets",
        required=True,
        help='TE consensus sequences in fasta format.')

    parser_input.add_argument(
        '-a', dest="annot",
        help='TE annotation in gff or bed format. Required for absence module.')

    parser_settings.add_argument(
        '-o', dest='outname',
        required=True,
        help='Sample name (prefix for output files).')


    # PROGRAM SETTINGS
    parser_settings.add_argument(
        '-m', dest='modus',
        nargs='+',
        required=True,
        choices=['tips','taps'],
        help='Program modus.')

    parser_settings.add_argument(
            '-c', dest='cpus',
            type= int, default=1,
            help='Number of CPUs. [1]')

    parser_settings.add_argument(
        '--region', dest='region',
        type=str,
        help='Restrict TE polymorphism search to specific region <chromosome:start:end>.')

    parser_settings.add_argument(
        '--include_invariant',
        action='store_true',
        help='Include conserved TEs in vcf output.')

    parser_settings.add_argument(
        '--require_split',
        action='store_true',
        help='Discard variant candidates for which there is no splitread evidence.')

    parser_settings.add_argument(
        '--keep',
        action='store_true',
        help='Keep intermediate files.')


    # THRESHOLDS
    parser_thresholds.add_argument(
        '-q', dest="mapq",
        type=int, default=20,
        help='Minimum mapping quality of reference-aligned reads. [20]')

    parser_thresholds.add_argument(
        '-lSR', dest='aln_len_SR',
        type=int, default=15,
        help='Minimum alignment length for splitread target hits. [15]')

    parser_thresholds.add_argument(
        '-lDR', dest='aln_len_DR',
        type=int, default=50,
        help='Minimum alignment length for discordant read-pair target hits. [50]')

    args=parser.parse_args()

    return args


#%% MAIN

def main():

    args = get_args()

    parameters = strumenti.mise_en_place(args)

    try:
        os.mkdir(args.outname + '_tmp')
    except OSError:
        pass
    os.chdir(args.outname + '_tmp')

    # Turn off warnings (switch with python -W on command line)
    if not sys.warnoptions:
        warnings.simplefilter("ignore")


    #%% TIPs
    if "tips" in parameters.modus:

        print('\nSEARCHING TE INSERTION POLYMORPHISMS')

        tips = strumenti.TIPs()

        # Get discordant read pairs and splitreads
        print('\nGetting candidate split reads and discordant read pairs from bam file ...')
        DR_anchors, splitreads, split_positions = strumenti.get_split_and_discordant_reads(parameters)
        print('Analysis begins with %i discordant read pairs and %i splitreads.\n' % \
              (len(DR_anchors), len(splitreads)))

        # Discordant read pairs: find anchor clusters with mates mapping to TEs
        if parameters.paired_end:
            tips.discordant_read_pairs(parameters, DR_anchors)
            print(str(len(tips.DR_clusters)) + ' clusters')

        # Same for splitreads
        tips.splitreads(parameters, splitreads, split_positions)
        print(str(len(tips.SR_clusters)) + ' clusters')

        # Combine split and discordant
        tips.combineDR_SR(parameters)

        # VCF output
        tips_vcf, tips_stats = tips.output_vcf(parameters)

        print('\nTIP search finished successfully')


    #%% TAPs
    if "taps" in parameters.modus:

        print('\n\nSEARCHING TE ABSENCE POLYMORPHISMS')
        print('\nExtracting read pairs around annotated TEs ...')

        taps = strumenti.TAPs(parameters)
        taps_vcf, taps_stats = taps.output_vcf(parameters)

        print('TAP search finished successfully\n')


    #%% Write VCF, stats, and log file

    if "tips" in parameters.modus and "taps" in parameters.modus:
        combined_vcf = tips_vcf + taps_vcf
        combined_vcf = sorted(combined_vcf, key=lambda x: (x[0], int(x[1])))

    elif "tips" in parameters.modus:
        combined_vcf = tips_vcf

    elif "taps" in parameters.modus:
        combined_vcf = taps_vcf

    date = time.strftime("%d/%m/%Y")

    metainfo = [

        '##fileFormat=VCFv4.2',
        '##fileDate=%s' % date,
        '##source==detettore v2.0',
        '##reference=%s' % parameters.reference,
        #'##contig=<ID=%s,length=%i,assembly=%s>',

        '##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
        '##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">',
        '##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">',
        '##INFO=<ID=DPADJ,Number=.,Type=Integer,Description="Read Depth of adjacency">',
        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
        '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants"',
        '##INFO=<ID=DR,Number=2,Type=Integer,Description="Discordant reads">',
        '##INFO=<ID=SR,Number=2,Type=Integer,Description="Split reads">',
        '##INFO=<ID=AL,Number=1,Type=Integer,Description="TE alignment length">',
        '##INFO=<ID=AP,Number=1,Type=Integer,Description="Proportion of TE covered by reads">',
        '##INFO=<ID=BPIQR,Number=1,Type=Integer,Description="Interquartile range of insertion breakpoints. Large for spurious read clusters">'
        
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">',

        # Multiple ALT alleles for INS, one per family!
        '##ALT=<ID=INS:ME,Description=Transposable element insertion polymorphism (TIP)>',
        '##ALT=<ID=DEL:ME,Description=Transposable element absence polymorphism (TAP)>',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % (parameters.outname)
        ]


    os.chdir('..')

    with gzip.open(args.outname + '.detettore.vcf.gz', 'wt') as f:

        f.write('\n'.join(metainfo)+'\n')
        for line in combined_vcf:
            f.write('\t'.join(map(str, line)) + '\n')


    # # Stats file
    if 'tips' in parameters.modus:
        print(
            'TIPs summary:\n0/1:%i\n1/1:%i\n' % (tips_stats['0/1'], tips_stats['1/1'])
            )

    if 'taps' in parameters.modus:
        print(
            'TAPs summary:\n0/0:%i\n0/1:%i\n1/1:%i\n./.:%i\n' % (taps_stats['0/0'],taps_stats['0/1'], taps_stats['1/1'], taps_stats['./.'])
            )

    # with open(args.outname + '.stats.tsv', 'w') as f:

    # Create log file
    logging.basicConfig(
            filename = args.outname +'.log',
            filemode='a',
            format = '%(levelname)-10s %(asctime)s %(message)s',
            level = logging.INFO
            )

    log = logging.getLogger(args.outname +'.log')
    
    readtype = 'PE' if parameters.paired_end else 'SE'
    
    log.info('readlength: %i, readtype: %s, isize_mean: %i, isize_stdev: %i' % (
        parameters.readlength, readtype, parameters.isize_mean, parameters.isize_stdev))

    log.info(args)

    # Clean up
    if not args.keep:
        shutil.rmtree(args.outname + '_tmp')

if __name__ == '__main__':
    main()
