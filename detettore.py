#!/usr/bin/env python

""" dettetore -a program to detect and characterize transposable element polymorphisms

Copyright (C) 2018 C. Stritt
License: GNU General Public License v3. See LICENSE.txt for details.
"""

import argparse
import logging
import os
import shutil

from scripts import nonreference_insertions
from scripts import reference_insertions
from scripts import bamstats


def get_args():

    parser = argparse.ArgumentParser(description='Detect TE polymorphisms from paired-end read data')

    parser.add_argument('-b', dest="bamfile", 
                        required=True,
                        help='bwa mem alignment to reference genome.')
    
    parser.add_argument('-m', dest='modus', 
                        nargs='+', required=True,
                        choices=['tips','taps'],
                        help='Program modus.')
    
    parser.add_argument('-o', dest='outfolder',
                        default='resultati', required=True,
                        help='Name of output folder.')

    parser.add_argument('-t', dest="targets", 
                        required=True,
                        help='Target sequences in fasta format.')

    parser.add_argument("-r", dest="reference",  
                        required=True,
                        help='Reference genome in fasta format. Required to calculate read statistics')

    parser.add_argument('-a', dest="annot",
                        help='TE annotation in gff or bed format. Required for absence module.')

    parser.add_argument('-u', dest="uniq",
                        type=int, default=30,
                        help='Difference between XS and AS for a read to be considered uniquely mapped. [30]')

    parser.add_argument('-l', dest='aln_len',
                        type=int, default=30,
                        help='Minimum alignment length for target hits. [30]')

    parser.add_argument('-id', dest='perc_id',
                        type=int, default=80,
                        help='Minimum percentage identity of target hits. [80]')
    
    parser.add_argument('-ws', dest='word_size',
                        type=int, default=11,
                        help='Word size (minimum length of best perfect match) for blasting splitreads against TEs. [11]')
    
    parser.add_argument('-keep',
                        action='store_true',
                        help='Keep intermediate files.')
    
    parser.add_argument('-reconstruct',
                        action='store_true',
                        help='For the absence module. Return sequences of TEs present in both reference and non-reference. Terribly slow!')

    parser.add_argument('-c', dest='cpus',
                        type= int, default=1,
                        help='Number of CPU used for the blast search. [1]')

    args=parser.parse_args()
    return args


def main():

    args = get_args()

    # Files
    bamfile = os.path.abspath(args.bamfile)
    targets = os.path.abspath(args.targets)
    reference = os.path.abspath(args.reference)
    
    # Alignment and mapping thresholds
    thresholds = [args.uniq, args.aln_len, args.perc_id, args.word_size]

    # Program settings
    modus = args.modus
    cpus = args.cpus
    reconstruct = args.reconstruct
    if 'taps' in modus:
        annotation_file = os.path.abspath(args.annot)

    try:
        os.mkdir(args.outfolder)
    except OSError:
        pass
    os.chdir(args.outfolder)

    # Get readlength, estimate mean coverage and insert size
    readinfo = [f for f in os.listdir('.') if f.endswith('_bamstats.txt')]
    if readinfo:
        print('Using ' + readinfo[0])
        readinfo = readinfo[0]
    else:
        print('Estimating bam read statistics\n...')
        readinfo = bamstats.write_output(bamfile, reference, 0.001)
    
    
    if 'tips' in modus:
        print('\nINCIPIT PRESENTIAE DETECTIO')
        nonreference_insertions.run_module(bamfile,
                                           targets,
                                           thresholds,
                                           readinfo,
                                           cpus)
        
        shutil.rmtree('blastdb')
        if not args.keep:
            for g in ['discordant.fasta', 'discordant.fasta.blast', 
                      'softclipped.fasta', 'softclipped.fasta.blast']:
                os.remove(g)
        
        print('FINIVIT PRESENTIAE DETECTIO\n')


    if 'taps' in modus:
        print('INCIPIT ABSENTIAE DETECTIO')
        reference_insertions.run_module(bamfile,
                                        readinfo,
                                        annotation_file,
                                        reference,
                                        thresholds,
                                        reconstruct,
                                        cpus)
        print('FINIVIT ABSENTIAE DETECTIO\n')
        
        
    # Create log file
    logging.basicConfig(
            filename = 'logfile.txt',
            format = '%(levelname)-10s %(asctime)s %(message)s',
            level = logging.INFO)
    log = logging.getLogger('logfile.txt')
    log.info(args)

if __name__ == '__main__':
    main()
