#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""


Created on Tue Mar  2 11:22:13 2021
@author: cristobal
"""


import os
import pysam
import sys

from scipy.stats import norm
from scripts import bamstats
from scripts import strumenti


import pysam
import re
import statistics

from Bio import SeqIO
from collections import Counter
from joblib import Parallel, delayed
from scipy.stats import norm
from scripts import strumenti




#%% VCF

def write_vcf(tips_out, taps_out, reference):

    """ Write VCF file for single strain/accession.
    """

    import time
    date = time.strftime("%d/%m/%Y")

    metainfo = [

        '##fileFormat=VCFv4.2',
        '##fileDate=(%s)' % date,
        '##source==detettore v0.6',
        '##reference=%s' % reference,
        '##contig=<ID=%s,length=%i,assembly=%s>',

        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
        '##INFO=<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">',
        '##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">',
        '##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">',
        '##INFO=<ID=DPADJ,Number=.,Type=Integer,Description="Read Depth of adjacency">',

        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality according to Li 2001">',
        '##FORMAT=<ID=DP,Number=2,Type=Integer,Description="Supporting reads. For insertions, this is the \
            sum of discordant read pairs and splitreads, for deletions it is the nr of read pairs \
            with deviant insert sizes">',



        '##FORMAT=<ID=DR,Number=2,Type=Integer,Description="Discordant reads">',
        '##FORMAT=<ID=SR,Number=2,Type=Integer,Description="Split reads">',
        '##FORMAT=<ID=BR,Number=2,Type=Integer,Description="<Number of reads bridgin the insertion breakpoint">',
        '##FORMAT=<ID=AL,Number=1,Type=Integer,Description="TE alignment length">',

        '##ALT=<ID=INS:ME,Description=description>',
        '##ALT=<ID=DEL:ME,Description=description>'
        ]






# TIPs
for contig in candidates:
    for site in candidates[contig]:


        """Get the target with the highest cumulative bit score,
            both from the 5' of the insertion site and from the 3'
        """
        name_l, name_r = strumenti.get_highest_scoring_target(site)

        d_anchors_l = site[name_l][0].discordant.anchors
        d_anchors_r = site[name_r][1].discordant.anchors
        s_anchors_l = site[name_l][0].splitreads.anchors
        s_anchors_r = site[name_r][1].splitreads.anchors

        A = (d_anchors_l and d_anchors_r)
        B = (s_anchors_l and s_anchors_r)

        if name_l == name_r:
            target_name = name_l
        else:
            target_name = '%s/%s' % (name_l, name_r)


        """ Strand of the insertion

        """
        if site[name_l][0].strand['plus'] >= site[name_l][0].strand['minus']:
            strand_l = '+'
        else:
            strand_l = '-'

        if site[name_r][1].strand['plus'] >= site[name_r][1].strand['minus']:
            strand_r = '+'
        else:
            strand_r = '-'

        strand = '%s/%s' % (strand_l, strand_r)


        """ Number of supporting reads
        """
        discordant_l, discordant_r = len(d_anchors_l), len(d_anchors_r)
        split_l, split_r = len(s_anchors_l), len(s_anchors_r)

        support = '%i/%i' % (discordant_l + split_l, discordant_r + split_r)
        nr_discordant = '%i/%i' % (discordant_l, discordant_r)
        nr_splitreads = '%i/%i' % (split_l, split_r)



        """ Number of nucleotides that could be aligned to the target
        """
        nr_aln_l = site[name_l][0].aligned_positions
        nr_aln_r = site[name_r][1].aligned_positions
        nr_aln = len(nr_aln_l | nr_aln_r)


        """ Guess insertion site position. If splitreads are present,
        these are used to derive the position and to extract, if present,
        the target site duplication (TSD)"""
        TSD = '-'

        if A:
            position = max(remove_outliers(site[name_l][0].discordant.end))

        if B:
            position = max(remove_outliers(site[name_l][0].splitreads.end))
            position_reverse_sr = min(remove_outliers(site[name_r][1].splitreads.start))

            """ If splitreads overlap, extract the overlapping
            sequence if it is shorter than 10 bp"""
            if position_reverse_sr < position:
                region = [chrmsm, position_reverse_sr + 1, position]
                TSD = consensus_from_bam(region, bamfile, [20, 30])
                if len(TSD) > 15:
                    TSD = '-'


        """ Coordinates of the insertion site region and the sequencing
        coverage in this region (mean and standard deviation"""
        region_start = min([site[name_l][0].region_start, site[name_r][1].region_start])
        region_end = max([site[name_l][0].region_end, site[name_r][1].region_end])
        region = [chrmsm, region_start, region_end]
        region_coverage = bamstats.local_cov(region, bamfile)


        """ Number of reads spanning the insertion site. If more than n
        non-split reads span it, the insertion is considered heterozygous"""
        n_bridge_reads = overlapping_reads(chrmsm, position, bamfile, 20)


        anchor_reads = ','.join(d_anchors_l) + ';' +\
                       ','.join(d_anchors_r) + ';' +\
                       ','.join(s_anchors_l) + ';' +\
                       ','.join(s_anchors_r)






# TAPs

for candidate in taps.candidates:

    if candidate == None:
         continue

    """ If there are deviant read pairs, check if the whole element
    or a part of it is missing
    """

    if candidate.nr_deviant > 0:

        chromosome = annot[candidate.indx].chromosome
        element_start = annot[candidate.indx].start
        element_end = annot[candidate.indx].end

        # Ignoring partial deletions
        if candidate.isize > candidate.length:


            """ Check if there are reads spanning the breakpoints
            """
            start_overlap = strumenti.overlapping_reads(
                chromosome,
                element_start,
                bamfile,
                20)

            end_overlap = strumenti.overlapping_reads(
                chromosome,
                element_end,
                bamfile,
                20)

            if start_overlap + end_overlap == 0 and \
               (candidate.length - 5*isize_stdev) < candidate.isize < (candidate.length + 5*isize_stdev):


                me_info = '%s,%i,%i,%s' % (
                    annot[candidate.indx].id,
                    element_start,
                    element_end,
                    annot[candidate.indx].strand
                    )

                INFO = 'SVTYPE=DEL;MEINFO=%s;END=%i;SVLEN=%i;DPADJ=%i' % \
                    (me_info, candidate.end, candidate.end-candidate.start, candidate.cov_mean)

                REF = contigs[chromosome][candidate.start - 2] # because candidate.start is 1-based

                FORMAT = '%s:%i:%i' % ('1/1', GQ, candidate.nr_deviant)

                vcf_line = [
                    chromosome,
                    candidate.start - 1,
                    '.',
                    REF,
                    'DEL:ME',
                    GQ,
                    'PASS',
                    INFO,
                    'GT:GQ:DP:'
                    FORMAT
                    ]




#%% Tabular output
class rawvariant:

    def __init__(self, line, accession):

        fields = line.strip().split('\t')

        self.source = accession

        self.chr = fields[0]
        self.pos = int(fields[1])
        self.TE = fields[2]
        self.strand = fields[3].split('/')

        self.support = map(int, fields[4].split('/'))
        self.discordant = map(int, fields[5].split('/'))
        self.split = map(int, fields[6].split('/'))
        self.aligned = int(fields[7])
        self.tsd = fields[8]
        self.nbridge = int(fields[9])

        self.region_start = int(fields[10])
        self.region_end = int(fields[12])
        self.coverage = int(fields[12])
        self.coverage_stdev = int(fields[13])


def tip_output_table(combined, bamfile):

    """ Write final output table.

   To do: replace with single vcf output for TIPs and TAPs!
    """

    header = [
        'chromosome',
        'position',
        'TE',
        'strand',
        'nr_supporting_reads',
        'nr_discordant_reads',
        'nr_splitreads',
        'nr_aligned_positions',
        'TSD',
        'nr_bridge_reads',
        'region_start',
        'region_end',
        'region_cov_mean',
        'region_cov_stdev']

    tabula = []

    for  chrmsm in combined:

        for site in combined[chrmsm]:


            """Get the target with the highest cumulative bit score,
            both from the 5' of the insertion site and from the 3'
            """
            name_l, name_r = get_highest_scoring_target(site)

            d_anchors_l = site[name_l][0].discordant.anchors
            d_anchors_r = site[name_r][1].discordant.anchors
            s_anchors_l = site[name_l][0].splitreads.anchors
            s_anchors_r = site[name_r][1].splitreads.anchors

            A = (d_anchors_l and d_anchors_r)
            B = (s_anchors_l and s_anchors_r)

            if name_l == name_r:
                target_name = name_l
            else:
                target_name = '%s/%s' % (name_l, name_r)


            """ Strand of the insertion

            """
            if site[name_l][0].strand['plus'] >= site[name_l][0].strand['minus']:
                strand_l = '+'
            else:
                strand_l = '-'

            if site[name_r][1].strand['plus'] >= site[name_r][1].strand['minus']:
                strand_r = '+'
            else:
                strand_r = '-'

            strand = '%s/%s' % (strand_l, strand_r)


            """ Number of supporting reads
            """
            discordant_l, discordant_r = len(d_anchors_l), len(d_anchors_r)
            split_l, split_r = len(s_anchors_l), len(s_anchors_r)

            support = '%i/%i' % (discordant_l + split_l, discordant_r + split_r)
            nr_discordant = '%i/%i' % (discordant_l, discordant_r)
            nr_splitreads = '%i/%i' % (split_l, split_r)



            """ Number of nucleotides that could be aligned to the target
            """
            nr_aln_l = site[name_l][0].aligned_positions
            nr_aln_r = site[name_r][1].aligned_positions
            nr_aln = len(nr_aln_l | nr_aln_r)


            """ Guess insertion site position. If splitreads are present,
            these are used to derive the position and to extract, if present,
            the target site duplication (TSD)"""
            TSD = '-'

            if A:
                position = max(remove_outliers(site[name_l][0].discordant.end))

            if B:
                position = max(remove_outliers(site[name_l][0].splitreads.end))
                position_reverse_sr = min(remove_outliers(site[name_r][1].splitreads.start))

                """ If splitreads overlap, extract the overlapping
                sequence if it is shorter than 10 bp"""
                if position_reverse_sr < position:
                    region = [chrmsm, position_reverse_sr + 1, position]
                    TSD = consensus_from_bam(region, bamfile, [20, 30])
                    if len(TSD) > 15:
                        TSD = '-'


            """ Coordinates of the insertion site region and the sequencing
            coverage in this region (mean and standard deviation"""
            region_start = min([site[name_l][0].region_start, site[name_r][1].region_start])
            region_end = max([site[name_l][0].region_end, site[name_r][1].region_end])
            region = [chrmsm, region_start, region_end]
            region_coverage = local_cov(region, bamfile)


            """ Number of reads spanning the insertion site. If more than n
            non-split reads span it, the insertion is considered heterozygous"""
            n_bridge_reads = overlapping_reads(chrmsm, position, bamfile, 20)


            anchor_reads = ','.join(d_anchors_l) + ';' +\
                           ','.join(d_anchors_r) + ';' +\
                           ','.join(s_anchors_l) + ';' +\
                           ','.join(s_anchors_r)

            outline = [
                chrmsm,
                position,
                target_name,
                strand,
                support,
                nr_discordant,
                nr_splitreads,
                nr_aln,
                TSD,
                n_bridge_reads,
                region_start,
                region_end,
                region_coverage[0],
                region_coverage[1],
                anchor_reads
                ]

            outline = map(str, outline)
            tabula.append(list(outline))

    tabula = sorted(tabula, key=lambda x: (x[0], int(x[1])))
    with open('tips.tsv', 'w') as f, \
         open('tips.supporting_reads.txt', 'w') as g:

        f.write('\t'.join(header) + '\n')
        g.write('discordant_left;discordant_right;split_left;split_right\n')
        for line in tabula:
            f.write('\t'.join(line[:-1]) + '\n')
            g.write(line[-1] + '\n')



def tap_output_table(out, annot, bamfile):
    """

    """
    header = [
        'chromosome',
        'feature_start',
        'feature_end',
        'feature_id',
        'feature_strand',
        'feature_length',
        'isize',
        'nr_deviant',
        'absence_start',
        'absence_end',
        'cov_mean',
        'cov_stdev'
        ]

    recs = []

    outfile = open('taps.tsv', 'w')
    outfile.write('\t'.join(header) + '\n')

    for candidate in out:

        if candidate == None:
             continue

        """ If there are deviant read pairs, check if the whole element
        or a part of it is missing
        """

        is_tap = False

        if candidate.nr_deviant > 0:

            chromosome = annot[candidate.indx].chromosome
            element_start = annot[candidate.indx].start
            element_end = annot[candidate.indx].end

            # easiest case: the whole element is absent
            if candidate.isize > candidate.length:


                """ Check if there are reads spanning the breakpoints
                """
                start_overlap = overlapping_reads(chromosome,
                                                            element_start,
                                                            bamfile, 20)

                end_overlap = overlapping_reads(chromosome,
                                                          element_end,
                                                          bamfile, 20)

                if start_overlap + end_overlap == 0 and \
                   (candidate.length - 5*isize_stdev) < candidate.isize < (candidate.length + 5*isize_stdev):
                    """ Way too restrictive in previous version (2*isize_stdev)"""

                    outline = candidate.write(annot)
                    outfile.write('\t'.join(outline) + '\n')
                    is_tap = True


        if not is_tap and reconstruct:

            annot_row = annot[candidate.indx]
            region = [annot_row.chromosome, annot_row.start, annot_row.end]
            seq, known = candidate.reconstruct(region, bamfile, filters)

            # How many reads bridge the start and end of the TE?
            five_bridge = overlapping_reads(annot_row.chromosome,
                                                      annot_row.start,
                                                      bamfile, 20)

            three_bridge = overlapping_reads(annot_row.chromosome,
                                                      annot_row.end,
                                                      bamfile, 20)

            infofield = 'covered:%s,five_nbridge:%s,three_nbridge:%s' % (known, five_bridge, three_bridge)

            seqrec = SeqRecord(Seq(seq),
                            id = '_'.join(map(str,region)),
                            name = annot_row.id,
                            description = infofield)
            recs.append(seqrec)

    if reconstruct:
        SeqIO.write(recs, 'omnipresent.fasta', 'fasta')

    outfile.close()



