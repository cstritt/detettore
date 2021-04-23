#!/usr/bin/env python

""" dettetore - a program to detect and characterize transposable element polymorphisms

Toolbox for TE polymorphism detection

Copyright (C) 2018 C. Stritt
License: GNU General Public License v3. See LICENSE.txt for details.
"""

import os
import statistics
import subprocess
import pysam
import sys
import numpy as np
import math

from scripts import bamstats

from Bio import SeqIO
from math import fabs
from random import random
from collections import Counter
from joblib import Parallel, delayed
from numba import prange

from scipy import percentile
from scipy.stats import norm

#%% CLASSI

class mise_en_place:

    """ Object for easy access to program settings

    """

    def __init__(self, args):

        # Files
        self.bamfile = os.path.abspath(args.bamfile)
        self.targets = os.path.abspath(args.targets)
        self.reference = os.path.abspath(args.reference)
        self.annotation_file = os.path.abspath(args.annot)

        # Thresholds
        self.uniq = args.uniq
        self.aln_len_DR = args.aln_len_DR
        self.aln_len_SR = args.aln_len_SR
        self.perc_id = args.perc_id
        self.word_size = args.word_size

        # Program settings
        self.modus = args.modus
        self.cpus = args.cpus

        # Parsed files
        self.annotation = []
        self.ref_contigs = {seq_record.id:
                            seq_record.seq for seq_record in SeqIO.parse(self.reference, "fasta")}


        # Load reference TE annotation
        if 'taps' in self.modus:

            # Check file format
            formt = self.annotation_file.split('.')[-1]

            # Read in annotation
            with open(self.annotation_file) as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')

                    if formt.startswith('bed'):

                        entry = {
                            'chromosome' : fields[0],
                            'start' : int(fields[1]),
                            'end' :  int(fields[2]),
                            'id' : fields[3],
                            'strand' : fields[5]
                            }

                    elif formt.startswith('gff'):

                        entry = {
                            # GFF uses 1-based indexing!
                            'chromosome' : fields[0],
                            'start' : int(fields[3]),
                            'end' : int(fields[4]),
                            'id' : fields[-1].split(';')[0].split('=')[-1],
                            'strand' : fields[6]
                            }

                    else:
                        sys.exit('Check annotation format.')

                    self.annotation.append(entry)

        # Estimate read statistics
        print('Estimating bam read statistics\n...')
        readinfo = bamstats.write_output(self.bamfile, 'dict')

        self.readlength = readinfo['readlength']
        self.isize_mean = readinfo['isize_mean']
        self.isize_stdev = readinfo['isize_stdev']

        for k in readinfo:
            print(k + ' : ' + str(readinfo[k]))

class TIPs:

    """
    Detect TE insertions absent in the reference and present in the
    resequenced accession using discordant read pairs and split reads
    """

    def __init__(self):

        self.DR_clusters = {} # discordant read pairs
        self.SR_clusters = {} # splitreads


    def discordant_read_pairs(self, parameters, DR_anchors):
        """
        Blast anchor mates against TE sequences and find clusters of anchor reads
        """

        print('Aligning discordant reads to target sequences...')
        DR_hits = minimap('discordant.fasta', parameters.targets)
        print('%i anchor mates map to a TE' % (len(DR_hits)))


        print('Clustering discordant reads...')
        # Only keep anchors with mate mapping to TE
        anchor_positions = get_anchor_positions(
            DR_hits,
            DR_anchors,
            parameters.bamfile,
            parameters.isize_mean)

        # Only keep clusters with more than 4 items
        DR_clusters = cluster_reads(
            anchor_positions,
            'discordant',
            overlap_proportion=0.05,
            min_cl_size=4
            )

        # Summarize clusters
        self.DR_clusters = summarize_clusters(
            DR_clusters,
            DR_anchors,
            DR_hits,
            'discordant')


    def splitreads(self, parameters, splitreads, split_positions):
        """
        Find clusters of splitreads and blast clipped parts of
        clustering reads against TE database. This is the opposite order than for
        clustering discordant read pairs, as in this way the time for the blast search
        can be greatly reduced.
        """

        print('\nClustering split reads...')
        SR_clusters = cluster_reads(
            split_positions,
            'splitreads',
            overlap_proportion=0.05,
            min_cl_size=4
            )

        print('Aligning split parts to target sequences...')
        fasta = write_clipped_to_fasta(
            SR_clusters,
            splitreads,
            parameters.aln_len_SR)

        SR_minimapd = minimap(fasta, parameters.targets)
        print('%i soft-clipped read parts map to a TE' % (len(SR_minimapd)))

        self.SR_clusters = summarize_clusters(
            SR_clusters,
            splitreads,
            SR_minimapd,
            'splitreads')


    def combineDR_SR(self, parameters):
        """
        Combine discordant read pair and splitread clusters,
        search for reads briding the suspected insertion site and
        thus supporting the reference allele (ie no TIP)

        """

        overlapping = []

        combined_clusters = []

        x = 0 # variable to store the position in the loop below, to avoid going through the whole thing every time

        # Find overlapping split and discordant clusters
        for i in prange(len(self.DR_clusters)):

            candidate = self.DR_clusters[i]
            overlap = False

            for j in prange(x, len(self.SR_clusters)):

                split_candidate = self.SR_clusters[j]

                if candidate[0] == split_candidate[0]:

                    overlaprange = range(split_candidate[1] - (2*parameters.readlength),
                                         split_candidate[1] + (2*parameters.readlength)
                                         )

                    if candidate[1] in overlaprange:

                        # Use position inferred from split reads
                        candidate[1] = split_candidate[1]

                        # Add split cluster summary
                        candidate.append(split_candidate[2])

                        combined_clusters.append(candidate)

                        overlapping.append((i, j))

                        overlap = True
                        x = j

            if not overlap:
                candidate.append(False)
                combined_clusters.append(candidate)

        for j in prange(len(self.SR_clusters)):

            if j not in [x[1] for x in overlapping]:

                self.SR_clusters[j].insert(2, False)
                combined_clusters.append(self.SR_clusters[j])

        self.candidates = sorted(combined_clusters, key=lambda x: (x[0], int(x[1])))



    def output_vcf(self, parameters):

        """ Positions in vcf are 1-based
        """

        vcf_lines = []

        for site in self.candidates:

            CHROM = site[0]
            POS = site[1]

            REF = parameters.ref_contigs[CHROM][POS - 1]
            ALT = '<INS:ME>'

            discordant = site[2]
            split = site[3]


            """ Which transposable element?
            """
            # Merge DR and SR information
            if discordant and split:
               te_hits = merge_TE_dict(site[2].te_hits, site[3].te_hits)

            elif discordant:
                te_hits = site[2].te_hits

            elif split:
                te_hits = site[3].te_hits


            # Get highest scoring TE
            best = get_highest_scoring_TE(te_hits)

            te = best[0]
            te_start = min(te_hits[te]['aligned_positions'])
            te_end = max(te_hits[te]['aligned_positions'])
            te_strand = '+' if te_hits[te]['strand']['+'] > te_hits[te]['strand']['-'] else '-'

            MEINFO = '%s,%i,%i,%s' % (te, te_start, te_end, te_strand)


            """ Genotype and genotype quality

            Infer genotype and genotype quality from split reads if present, else
            from discordant read pairs

            """
            if split:

                split.get_REF_support(parameters, 'split')

                ref_Q = [split.ref_support[x] for x in split.ref_support]
                alt_Q = [split.te_hits[te]['combined_mapqs'][x] for x in split.te_hits[te]['combined_mapqs']]

                genotype = get_genotype(ref_Q, alt_Q)


            else:

                discordant.get_REF_support(parameters, 'discordant')

                ref_Q = [discordant.ref_support[x] for x in discordant.ref_support]
                alt_Q = [discordant.te_hits[te]['combined_mapqs'][x] for x in discordant.te_hits[te]['combined_mapqs']]

                genotype = get_genotype(ref_Q, alt_Q)

            # if genotype[0] == '0/0':
            #     continue


            """ Nr of supporting discordant and splitreads
            """
            DR = len(discordant.te_hits[te]['hit_mapqs']) if discordant else 0
            SR = len(split.te_hits[te]['hit_mapqs']) if split else 0


            """ Number of TE positions covered
            """
            AL = len(te_hits[te]['aligned_positions'])


            """ Target site duplication:
            If splitreads overlap, extract the overlapping sequence if it is shorter than 15 bp
            """

            if split:
                position_reverse_sr = min(remove_outliers(split.breakpoint[1]))

                if position_reverse_sr < POS:
                    region = [CHROM, position_reverse_sr + 1, POS]
                    HOMSEQ = consensus_from_bam(region, parameters.bamfile, [20, 30])
                    if len(HOMSEQ) > 20:
                        HOMSEQ = ''

            else:
                HOMSEQ = ''

            HOMLEN = len(HOMSEQ)


            """ Regional coverage
            """
            if discordant:
                region = [CHROM, discordant.region_start, discordant.region_end]

            else:
                region = [CHROM, split.region_start, split.region_end]

            region_coverage = bamstats.local_cov(region, parameters.bamfile)
            DPADJ = region_coverage[0]

            """ Confidence interval for imprecise variants
            """
            if not split:

                closest = [
                    max(remove_outliers(site[2].breakpoint[0])),
                    min(remove_outliers(site[2].breakpoint[1]))
                    ]

                CI_lower = min(closest)
                CI_upper = max(closest)


            """ Put together FORMAT and INFO fields
            """

            # Imprecise event if there is no splitread information
            if split:

                INFO = 'MEINFO=%s;SVTYPE=INS;HOMSEQ=%s;HOMLEN=%i;DPADJ=%i;DR=%i;SR=%i;AL=%i' \
                    % (MEINFO, HOMSEQ, HOMLEN, DPADJ, DR, SR, AL)
            else:
                INFO = 'MEINFO=%s;SVTYPE=INS;HOMSEQ=%s;HOMLEN=%i;DPADJ=%i;DR=%i;SR=%i;AL=%i;IMPRECISE;CIPOS=%i,%i' \
                    % (MEINFO, HOMSEQ, HOMLEN, DPADJ, DR, SR, AL, CI_lower, CI_upper)

            FORMAT = 'GT:GQ:DP'
            GT = '%s:%i:%i' % (genotype[0], genotype[1], DR+SR)


            outline = [CHROM, POS, '.', REF, ALT, genotype[1], 'PASS', INFO, FORMAT, GT]
            vcf_lines.append(outline)

        return vcf_lines


class TAPs:

    def __init__(self, parameters):

        isize_stats = (parameters.isize_mean,
                       parameters.isize_stdev,
                       parameters.readlength)

        print('\nExtracting read pairs with deviant insert sizes around annotated features\n...')
        inputs = [(i, x) for i, x  in enumerate(parameters.annotation)]

        self.candidates = Parallel(n_jobs=parameters.cpus)(delayed(self.process_taps)\
                                               (i,
                                                isize_stats,
                                                parameters.bamfile,
                                                parameters.ref_contigs,
                                                parameters.uniq) for i in inputs
                                               )


    def extract_deviant_reads(self, region, isize_stats, bamfile, mapq):

        """
        Returns a dictionary d[read name] = [isize, AlignedSegment objects]
        """

        pybam = pysam.AlignmentFile(bamfile, "rb")

        isize_mean, isize_stdev, readlength = isize_stats
        q_upper = int(norm.ppf(0.99, isize_mean,isize_stdev))

        chromosome, start, end = region

        deviant_read_pairs = {}

        for read in pybam.fetch(chromosome, start, end):

            name = read.query_name

            # bad read pair
            if read.mapping_quality < mapq or not read.is_paired:
                continue

            # discard non-properly oriented reads
            if (read.mate_is_reverse and read.is_reverse) or \
                (not read.is_reverse and not read.mate_is_reverse):
                continue

            # insert size not deviant
            isize = abs(read.template_length) - (2*readlength)
            if not (isize > q_upper and isize < 30000):
                continue


            """ Get neat deviant reads with 5' on forward and 3' on reverse
            """

            if not read.is_reverse and read.mate_is_reverse:
                deviant_read_pairs[name] = [isize, read]

            elif name in deviant_read_pairs and read.is_reverse and not read.mate_is_reverse:
                deviant_read_pairs[name].append(read)

        pybam.close()

        remove = []

        for name in deviant_read_pairs:
            if len(deviant_read_pairs[name]) != 3:
                remove.append(name)

        for name in remove:
            deviant_read_pairs.pop(name)

        return deviant_read_pairs


    class TAP_candidate:

        def __init__(self, fields):

            self.annotated_TE = {
                'indx' : fields[0],
                'length'  : fields[1]
                }

            self.context = {
                'region_start' : fields[2],
                'region_end' : fields[3],
                'mean_coverage' : fields[4]
                }

            self.deviant_reads = {
                'nr_deviant' : int(),
                'deviation' : int(),
                'start' : int(),
                'end' : int(),
                'mapqs' : {}
                }

        def add_deviant_info(self, deviant_summary):

            self.deviant_reads['nr_deviant'] = len(deviant_summary)
            self.deviant_reads['deviation'] = statistics.median(
                [deviant_summary[k][0] for k in deviant_summary]
                )

            try:
                self.deviant_reads['start'] = max(
                    [deviant_summary[k][1].reference_end for k in deviant_summary]
                    )

                self.deviant_reads['end'] = min(
                    [deviant_summary[k][2].reference_start for k in deviant_summary]
                    )

            except ValueError:
                pass

            # Read mapping qualities
            """ Same format as for TIPs: dictionary[read name] = mapq
            """
            for k in deviant_summary:

                mapq_sum = deviant_summary[k][1].mapq + deviant_summary[k][2].mapq
                self.deviant_reads['mapqs'][k] = mapq_sum


        #def add_reference_support(self):


    def process_taps(
            self,
            annot_entry,
            isize_stats,
            bamfile,
            contigs,
            uniq):

        """ Check if there are deviant read-pairs around the annotated feature
        """

        i, feature = annot_entry
        isize_mean, isize_stdev, readlength = isize_stats

        region_start = feature['start'] - (isize_mean + isize_stdev)
        region_end = feature['end'] + (isize_mean + isize_stdev)

        # Reset coordinates if they are 'outside' chromosomes
        if region_start < 0:
            region_start = 0

        if region_end > len(contigs[feature['chromosome']]):
            region_end = len(contigs[feature['chromosome']]) - 1

        region = (feature['chromosome'], region_start, region_end)

        feature_length = feature['end'] - feature['start']


        deviant = self.extract_deviant_reads(region, isize_stats, bamfile, uniq)

        if not deviant:
            return None

        cov_mean, cov_stdev = bamstats.local_cov(region, bamfile)

        summary = self.TAP_candidate(
            [
                i,
                feature_length,
                region_start,
                region_end,
                cov_mean
                ]
            )

        summary.add_deviant_info(deviant)

        return summary


    def output_vcf(self, parameters):


        vcf_lines = []

        for i, site in enumerate(self.candidates):

            te = parameters.annotation[i]

            CHROM = te['chromosome']
            POS = te['start']

            REF = parameters.ref_contigs[CHROM][POS - 1]
            ALT = '<DEL:ME>'

            MEINFO = '%s,%i,%i,%s' % (te['id'], 1, te['end']-te['start'], te['strand'])
            SVLEN = -(te['end']-te['start'])


            if site:

                if site.deviant_reads['deviation'] < site.annotated_TE['length']:
                    continue

                """ Get reads supporting reference state, ie the presence of the annotated TE
                """

                region_start = site.context['region_start']
                region_end = site.context['region_end']
                DPADJ = site.context['mean_coverage']

                ref_support = {}

                pybam = pysam.AlignmentFile(parameters.bamfile, "rb")

                for read in pybam.fetch(CHROM, region_start, region_end):

                    if read.mapq == 0:
                        continue

                    # Ignore reads identified as deviant by detettore but as properly paired by the aligner
                    if read in site.deviant_reads['mapqs']:
                        continue

                    # Exclude clipped reads
                    if 'S' in read.cigarstring or 'H' in read.cigarstring:
                        continue

                    if read.is_proper_pair:

                        if read.query_name in ref_support:
                            ref_support[read.query_name] += read.mapping_quality

                        else:

                            insert_interval = range(read.reference_end, read.next_reference_start)

                            if \
                                te['start'] in insert_interval and not te['end'] in insert_interval or \
                                te['end'] in insert_interval and not te['start'] in insert_interval:

                                ref_support[read.query_name] = read.mapping_quality

                pybam.close()


                """ Calculate genotype likelihood and quality
                """
                ref_Q = [ref_support[k] for k in ref_support]
                alt_Q = [site.deviant_reads['mapqs'][k] for k in site.deviant_reads['mapqs']]


                genotype = get_genotype(ref_Q, alt_Q)

                # if genotype[0] == '0/0':
                #     continue

                INFO = 'MEINFO=%s;SVTYPE=DEL;SVLEN=%i;DPADJ=%i' % (MEINFO, SVLEN, DPADJ)
                FORMAT = 'GT:GQ:DP'
                GT = '%s:%i:%i' % (genotype[0], genotype[1], site.deviant_reads['nr_deviant'])


                outline = [CHROM, POS, '.', REF, ALT, genotype[1], 'PASS', INFO, FORMAT, GT]
                vcf_lines.append(outline)

        return vcf_lines



class read_cluster:

    def __init__(self, reads):

        self.chromosome = str()
        self.position = int()

        # Set with names of clustering reads
        self.reads = reads
        self.te_hits = {}
        self.breakpoint = [ [] , [] ]
        self.region_start = float('inf')
        self.region_end = -float('inf')

        # Mapping quality of reads bridging the supposed TE insertion site,
        # thus providing evidence against an insertion. Used to calculate genotype quality
        self.ref_support = {}


    def combine_hits_and_anchors(self, anchors, hits, modus, weight=1):

        for readname in self.reads:

            if readname not in hits:
                continue

            read = anchors[readname]
            cigar = read.cigartuples

            # 5' or 3' site?
            if modus == 'discordant':
                s = 1 if read.is_reverse else 0

            elif modus == 'splitreads':
                s = 1 if cigar[0][0] == 4 else 0

            # region for local coverage
            if read.reference_start < self.region_start:
                self.region_start = read.reference_start
            if read.reference_end > self.region_end:
                self.region_end = read.reference_end

            # Insertion breakpoint
            break_pos = read.reference_end if s == 0 else read.reference_start
            self.breakpoint[s].append(break_pos)

            # TE hits dictionary
            hit = hits[readname]
            target = hit['target_name']

            if target not in self.te_hits:

                self.te_hits[target] = {
                    'strand' : {'+' : 0, '-' : 0},
                    'hit_mapqs' : {},
                    'anchor_mapqs' : {},
                    'combined_mapqs' : {},
                    'aligned_positions' : set()
                    }


            # Strand of the insertion
            if modus == 'splitreads':
                strand = '+' if hit['strand'] == '+' else '-'


            elif modus == 'discordant':

                # Tedious because the bam files produced by bwa report sequences
                # "on the same strand as the reference": if a read is_reverse, its
                # reverse complement is given in the bam file and mapped by detettore against the TEs
                # minimap: + if query/target on same strand

                # Anchor on forward strand
                if not read.is_reverse:

                    # FF
                    if not read.mate_is_reverse: # mate shoud be reverse but is mapped as forward
                        strand = '-' if hit['strand'] == '+' else '+'

                    # FR
                    elif read.mate_is_reverse:
                        strand = '+' if hit['strand'] == '+' else '-'

                elif read.is_reverse:

                    if not read.mate_is_reverse:
                        strand = '+' if hit['strand'] == '+' else '-'

                    if read.mate_is_reverse:
                        strand = '-' if hit['strand'] == '+' else '+'


            self.te_hits[target]['strand'][strand] += 1

            self.te_hits[target]['hit_mapqs'][readname] = hit['mapping_quality']
            self.te_hits[target]['anchor_mapqs'][readname] = read.mapping_quality
            self.te_hits[target]['combined_mapqs'][readname] = hit['mapping_quality'] + read.mapping_quality

            # Nr aligned positions on the target
            self.te_hits[target]['aligned_positions'].update(
                range(hit['target_start'], hit['target_end'])
                )


    def get_highest_scoring_TE(self):
        """ If reads map to multiple TEs, get the TE with
        the highest cumulative mapping quality

        """
        score = 0
        best = ''

        for te in self.te_hits:
            mapq_score = sum(self.te_hits[te]['hit_mapqs'])
            if mapq_score > score:
                score = mapq_score
                best = te
        return best, score


    def get_REF_support(self, parameters, modus, overlap=20):

        """
        Obtain reads and read pairs supporting the reference allele.
        Returns a list of discordant read pairs bridging the insertion breakpoint and their mapq,
        and a list of single reads spanning the insertion breakpoint, and their mapq.

        """

        pybam = pysam.AlignmentFile(parameters.bamfile, "rb")

        # Single reads bridging breakpoint
        if modus == 'splitreads':
            for read in pybam.fetch(self.chromosome, self.position, self.position+1):

                if read.mapq == 0:
                    continue

                down = [x for x in range(read.reference_start, read.reference_end) if
                        x < self.position]
                up = [x for x in range(read.reference_start, read.reference_end) if
                        x > self.position]

                if len(down) > overlap and len(up) > overlap:

                    self.ref_support[read.query_name] = read.mapping_quality

        # Properly paired reads spanning break point
        if modus == 'discordant':

            for read in pybam.fetch(self.chromosome, self.region_start, self.region_end):

                if read.mapq == 0:
                    continue

                if read.is_proper_pair:

                    if read.query_name in self.ref_support:
                        self.ref_support[read.query_name] += read.mapping_quality


                    elif self.position in range(read.reference_end, read.next_reference_start):
                        self.ref_support[read.query_name] = read.mapping_quality

        pybam.close()


#%% FUNZIONI

def get_split_and_discordant_reads(parameters):

        """
        Traverse bam and extract disordant read pairs and split reads.
        Returns 2 dictionaries 'd[read_name] = pybam read object', and a dictionary 'd[chrom]: pos'
        with the splitread positions.

        A read pair is discordant if one read maps uniquely to the reference,
        i.e. AS - XS > uniqueness_threshold, while its mate does not map properly.

        Write uniquely mapping reads to anchors dicitonary, their mate to fasta file.

        To do: make compatible with minimap2

        """
        pybam = pysam.AlignmentFile(parameters.bamfile, "rb")

        discordant_anchors = dict()
        discordant_mates = dict()
        fasta = open('discordant.fasta', 'w')

        splitreads = dict()
        split_positions = dict()

        for read in pybam.fetch():

            clipped = is_softclipped(read, parameters.aln_len_SR)

            if (read.is_proper_pair and not clipped) or read.is_secondary:
                continue

            name = read.query_name

            """
            While BWA always outputs the XS tag, Bowtie2 only does it when
            there is a secondary alignment
            """
            AS = read.get_tag('AS')
            XS = read.get_tag('XS')


            """
            Splitreads. Only reads are considered where exactly one end
            is clipped off.
            """
            if clipped:

                if AS-XS < parameters.uniq:
                    continue

                clseq = read.query_sequence
                if clseq.count('N')/float(len(clseq)) > 0.1:
                    continue

                cigar = read.cigartuples

                chrmsm = pybam.getrname(read.reference_id)
                strt, end = read.reference_start, read.reference_end

                if read.is_read1:
                    name = name +'/1'
                else:
                    name = name + '/2'

                splitreads[name] = read

                # get positions of the splitreads
                if cigar[0][0] == 4:
                    interval = (end - parameters.readlength, end)

                elif cigar[-1][0] == 4:
                    interval = (strt, strt + parameters.readlength)

                if chrmsm not in split_positions:
                    split_positions[chrmsm] = []

                split_positions[chrmsm].append([name, interval[0], interval[1]])

            """
            Discordant read pairs
            """
            if not read.is_proper_pair:

                if AS-XS < parameters.uniq:
                    uniq = 0
                else:
                    uniq = 1

                if read.is_read1:
                    read_nr = 1
                else:
                    read_nr = 2

                if name in discordant_mates:
                    uniq_mate = discordant_mates[name][2]
                    if uniq + uniq_mate == 1:

                        mate_read = discordant_mates[name][0]

                        if uniq == 1:

                            seq = mate_read.query_sequence
                            if seq.count('N')/float(len(seq)) > 0.1:
                                continue

                            discordant_anchors[name] = read
                            header = '>' + name
                            fasta.write(header + '\n' + seq + '\n')

                        else:
                            seq = read.query_sequence
                            if seq.count('N')/float(len(seq)) > 0.1:
                                continue

                            discordant_anchors[name] = mate_read
                            header = '>' + name
                            fasta.write(header + '\n' + seq + '\n')

                    del discordant_mates[name]

                else:
                    discordant_mates[name] = [read, read_nr, uniq]

        fasta.close()
        pybam.close()

        return [discordant_anchors, splitreads, split_positions]


def read_stats(bamfile, proportion):
    """
    Takes a file containing insert sizes and calculates mean and standard
    deviation of the core distribution. Using the core distribution instead of
    all data mitigates the influence of outliers (source: Piccard toolbox)

    """

    pybam = pysam.AlignmentFile(bamfile, "rb")

    isizes_out = "isizes.txt"
    readinfo_out = "readinfo.txt"

    isizes = list()
    readlength = int()

    for read in pybam.fetch():

        while not readlength:
            readlength = read.infer_query_length()

        if read.is_proper_pair:
            if random() <= proportion:
                isizes.append(abs(read.isize))
            continue

    # calculate median absolute deviation

    median = statistics.median(isizes)
    abs_deviations = [fabs(x - median) for x in isizes]
    mad = statistics.median(abs_deviations)

    # calculate mean and standard deviation of core distribution
    under_limit = median - (10*mad)
    upper_limit = median + (10*mad)

    core_dist = [i for i in isizes if i > under_limit and i < upper_limit]

    mean = int(statistics.mean(core_dist))
    stdev = int(statistics.stdev(core_dist))

    # write to files
    out = ["readlength:"+str(readlength), "isize_mean:"+str(mean),
           "isize_stdev:"+str(stdev)]

    with open(readinfo_out, "w") as f:
        f.write('\n'.join(out))

    with open(isizes_out, "w") as f:
        f.write('\n'.join(str(x) for x in isizes))

    pybam.close()

    return readinfo_out


def remove_outliers(lista):
    # outliers defined as in R's boxplots
    q_1 = percentile(lista, 25)
    q_3 = percentile(lista, 75)

    boxlength = q_3 - q_1

    lower = max(min(lista), q_1 - 1.5*boxlength)
    upper = min(max(lista), q_3 + 1.5*boxlength)

    lista_filt = [x for x in lista if x >= lower and x <= upper]
    return lista_filt


def is_overlapping(a,b):
    a_range = set(range(a[0], a[1]))
    b_range = set(range(b[0], b[1]))

    intersection = a_range & b_range
    union = a_range | b_range

    if len(intersection) == 0:
        return False
    else:
        return len(intersection) / float(len(union))


def is_softclipped(read, min_clip_length):
    """
    Return FALSE if read is not clipped.
    Only reads are used where exactly one end is clipped off and not both;
    and where the clipped part => min_clip_length
    """

    cigar = read.cigartuples
    if cigar:
        first, last = cigar[0], cigar[-1]

        a = (first[0] == 4 and first[1] >= min_clip_length)
        b = (last[0] == 4 and last[1] >= min_clip_length)

        if (a and not b) or (b and not a):
             return True
        else:
            return False
    else:
        return False


def write_clipped_to_fasta(clusters, reads, min_clip_length):

    fasta = 'softclipped.fasta'
    f = open(fasta, 'w')

    for k in clusters:
        for cl in clusters[k]:

            for name in cl:
                header = '>' + name
                read = reads[name]
                cigar = read.cigartuples
                readseq = read.query_sequence
                seq = clip_seq(readseq, cigar)

                if len(seq) < min_clip_length:
                    continue
                outline = header + '\n' + seq + '\n'
                f.write(outline)
    f.close()
    return fasta


def cluster_reads(positions, modus, overlap_proportion=0.05, min_cl_size=4):

    """
    Takes a dictionary with chromosomes as keys and ordered intervals.
    Adjacent intervals are assigned to the same cluster if they overlap at
    least by the specified overlap_proportion. The latter is an important
    parameter, since if defined too loosely, read clusters might include
    reads which map to other close-by TE insertions

    """

    clusters = dict()
    chromosomi = sorted(positions.keys())

    for k in chromosomi:

        cluster_indx = 0
        cluster_switch = 0

        clusters[k] = list()
        positions_chrmsm = positions[k]

        for i in prange(len(positions_chrmsm)):

            line = positions_chrmsm[i]
            name = line[0]

            # add 100 bp to splitread positions. this improves the clustering of
            # splitreads for short DNA transposon insertions, which often come with
            # a mapping gap such that splitreads from the 5' and 3' do not overlap
            # anymore and don't cluster

            if modus == 'splitreads':
                strt, end = line[1]-100, line[2]+100
            else:
               strt, end = line[1], line[2]
            interval = strt, end

            # start of a new chromosome
            if i == 0:
                previous = [name, interval]
                continue

            # main algorithm within chromosome
            overlap = is_overlapping(interval, previous[1])
            if overlap >= overlap_proportion:

                try:  # extend exisiting cluster
                    clusters[k][cluster_indx].add(name)

                except IndexError:  # beginn a new cluster
                    clusters[k].append({name})
                    clusters[k][cluster_indx].add(previous[0])
                    cluster_switch = 1

            else:  # close cluster if there is an open one
                if cluster_switch == 1:
                    cluster_indx += 1
                    cluster_switch = 0

            previous = [name, interval]

    # Finally, remove clusters with less than min_cl_size items
    clusters_filt = dict()
    for k in clusters:
        for cl in clusters[k]:

            if len(cl) < min_cl_size:
                continue
            else:
                if k not in clusters_filt:
                    clusters_filt[k] = [cl]
                else:
                    clusters_filt[k].append(cl)

    return clusters_filt


def create_blastdb(fasta):

    outname = 'blastdb/targets'

    cmd = ['makeblastdb',
           '-in', fasta,
           '-out', outname,
           '-dbtype', 'nucl',
           '-parse_seqids']

    subprocess.call(cmd, stdout=open(os.devnull, 'w'))


def blastn(query, min_perc_id, word_size, cpus):

    outfile = query + '.blast'

    cmd = ['blastn',
           '-query', query,
           '-db', 'blastdb/targets',
           '-outfmt', '6 qseqid length bitscore sstrand sseqid sstart send qseq',
           '-word_size', str(word_size),
           '-perc_identity', str(min_perc_id),
           #'-max_target_seqs', '1', # triggers Warning: [blastn] Examining 5 or more matches is recommended
           '-num_threads', str(cpus)]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output = proc.stdout.read()
    blast_out = output.splitlines()

    with open(outfile, 'w') as f:
        for line in blast_out:
            f.write(line.decode('utf-8') + '\n')

    return blast_out


def minimap(queries, targets):

    cmd = ['minimap2', '-x', 'sr', targets, queries]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    output_raw = proc.stdout.read().splitlines()

    minimapd = {}

    for line in output_raw:

        line = line.decode('utf-8')

        fields = line.strip().split()
        readname = fields[0]

        if readname in minimapd:
            if minimapd[readname]['aln_block_length'] > int(fields[10]):
                continue

        minimapd[readname] = {

            #'query_len' : fields[1],
            #'query_start' : fields[2],
            #'query_end' : fields[3],
            'strand' : fields[4],
            'target_name' : fields[5],
            #'target_length' : fields[6],
            'target_start' : int(fields[7]),
            'target_end' : int(fields[8]),
            #'nr_residue_matches' : fields[9], # number of sequence matches
            'aln_block_length' : int(fields[10]), # total nr of matches, mismatches, indels
            'mapping_quality' : int(fields[11])

            }

    return minimapd


def hit_dictionary(blast_out, min_aln_len):

    """ Dictionary of TE alignments with read IDs as keys. Makes use of the
    'hit' class defined above

    """

    TEhits = dict()

    for line in blast_out:

        fields = line.decode('utf-8').split('\t')
        aln_len = float(fields[1])

        if aln_len < min_aln_len:
            continue

        else:
            name = fields[0]
            TEhit = {
                'read' : fields[0],
                'aln_len' : float(fields[1]),
                'score' : float(fields[2]),
                'target_strand' : fields[3],
                'target_id' : fields[4],
                'target_aln_start' : int(fields[5]),
                'target_aln_end' : int(fields[6]),
                'seq' : fields[7]
                }

            # Only use best hit
            if name in TEhits:

                previous_score = TEhits[name][0]['score']
                new_score = TEhit['score']

                if new_score == previous_score:
                    TEhits[name].append(TEhit)
                elif new_score < previous_score:
                    continue
                elif new_score > previous_score:
                    TEhits[name] = [TEhit]

            else:
                TEhits[name] = [TEhit]

    return TEhits


def clip_seq(seq, cigar):
    """
    Return clipped part of sequence, as indicated in cigar. Assumes that there
    is only one clipped part at the beginning or end of the read. Cigar input
    parameter is a pysam cigartuples object.
    """
    if cigar[0][0] == 4:
        end = cigar[0][1]
        clipseq = seq[:end]
    elif cigar[-1][0] == 4:
        strt = sum([x[1] for x in cigar[:-1]])
        end = strt + cigar[-1][1]
        clipseq = seq[strt:end]
    return clipseq



def consensus_from_bam(region, bamfile, filters):

    """ Create a pileup file for a region and extract consensus sequence
    Quality filtering:
            https://gist.github.com/amblina/5ed9a61ce74ad90668f4e29d62e7eb79
    """

    pybam = pysam.AlignmentFile(bamfile, "rb")

    chrmsm, strt, end = region[0], region[1]-1, region[2]
    baseq, mapq = filters

    pile_dict = dict()
    crap_reads = set()

    for pileupcolumn in pybam.pileup(chrmsm, strt, end, **{"truncate": True}):
        pos = pileupcolumn.pos
        pile_dict[pos] = []

        for pileupread in pileupcolumn.pileups:

            read = pileupread.alignment

            if read.query_name in crap_reads:
                continue

            if pileupread.is_del or pileupread.is_refskip:
                continue

            elif read.mapping_quality < mapq:
                crap_reads.add(read.query_name)
                continue

            elif read.query_qualities[pileupread.query_position] < baseq:
                continue

            base = pileupread.alignment.query_sequence[pileupread.query_position]
            pile_dict[pos].append(base)

    consensus_seq = ''
    for i in range(strt, end):
        try:
            bases = pile_dict[i]
            cov = len(bases)
            if cov == 0:
                consensus_seq += 'N'
                continue
            base_counts = Counter(bases)

            # most common base
            cons_base_list = base_counts.most_common(1)
            consensus_base = max(cons_base_list, key=lambda x: x[1])[0]

            consensus_seq += consensus_base

        except KeyError:
            consensus_seq += 'N'
            continue

    pybam.close()
    return consensus_seq


def get_anchor_positions(hits, anchors, bamfile, isize):

    """ Go through the blast hits and get the corresponding anchor positions
    """

    pybam = pysam.AlignmentFile(bamfile, "rb")
    positions = dict()

    for read in hits:

        anchor = anchors[read]
        chromosome = pybam.getrname(anchor.reference_id)
        strt, end = anchor.reference_start, anchor.reference_end

        if anchor.is_reverse:
            interval = strt-isize+10, strt
        else:
            interval = end-10, end + isize

        try:
            positions[chromosome].append([read, interval[0], interval[1]])
        except KeyError:
            positions[chromosome] = [[read, interval[0], interval[1]]]

    for k in positions.keys():
        positions[k] = sorted(positions[k], key=lambda x: int(x[1]))

    pybam.close()
    return positions


def summarize_clusters(clusters, anchors, hits, modus):

    cl_summaries = []

    for chromosome in clusters:

        for reads in clusters[chromosome]:

            c = read_cluster(reads)

            c.combine_hits_and_anchors(
                anchors,
                hits,
                modus)

            if not c.te_hits or \
                len(c.breakpoint[0]) < 2 or len(c.breakpoint[1]) < 2:
                    continue

            # Breakpoint
            """ Remove positional outliers due to stray reads
            """
            bp = max(remove_outliers(c.breakpoint[0]))
            c.position = bp
            c.chromosome = chromosome

            cl_summaries.append([chromosome, bp, c])

    return cl_summaries


def merge_TE_dict(DR_TE_hits, SR_TE_hits):
    """
    Combine DR and SR evidence on the identity and polarity of the TIP,
    by merging te_hit dictionaries saved in the read_cluster class

    """

    te_hits = DR_TE_hits

    # Shared TEs
    shared = [te for te in te_hits if te in SR_TE_hits]

    # TEs only detected through short reads
    SR_only = [te for te in SR_TE_hits if te not in te_hits]

    for te in shared:

        te_hits[te]['aligned_positions'].update(SR_TE_hits[te]['aligned_positions'])

        te_hits[te]['strand']['+'] += SR_TE_hits[te]['strand']['+']
        te_hits[te]['strand']['-'] += SR_TE_hits[te]['strand']['-']


        # Mapping qualities used for calculating genotype likelihoods
        # -> read_cluster.combibe_hits_and_anchors()
        for k in ['anchor_mapqs', 'hit_mapqs']:

            for read in SR_TE_hits[te][k]:

                # if read in te_hits[te][k]:
                #     print(read)

                te_hits[te][k][read] = SR_TE_hits[te][k][read]

    for te in SR_only:
        te_hits[te] = SR_TE_hits[te]


    return te_hits


def get_highest_scoring_TE(te_hits):
    """ If reads map to multiple TEs, get the TE with
    the highest cumulative mapping quality

    """
    score = 0
    best = ''

    for te in te_hits:
        mapq_score = sum([te_hits[te]['hit_mapqs'][x] for x in te_hits[te]['hit_mapqs']])
        if mapq_score > score:
            score = mapq_score
            best = te
    return best, score


def get_genotype(ref_Q, alt_Q):

    """ Genotype likelihood and quality

    INPUT: two dictionaries of the form d[read name] = read quality

    Phred score Q = -10*log10(P) <->  P = 10^(-Q/10)

    Genotype qualities:
    https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs

    """

    # Convert phred score Q to error probability

    # ... reference-supporting reads
    ref_errP = [pow(10, (-x/10)) for x in ref_Q]

    # ... TE-supporting reads
    alt_errP = [pow(10, (-x/10)) for x in alt_Q]


    # Calculate likelihoods
    gt_lik = []
    for n_ref_alleles in [0, 1, 2]:

        gt_lik.append(GL(n_ref_alleles, ref_errP, alt_errP, 2))

    gt_phred = [-10 * np.log10(x) for x in gt_lik]

    minQ = min(gt_phred)
    gt = gt_phred.index(minQ)

    if gt == 0:
        GT = '1/1'
    elif gt == 1:
        GT = '0/1'
    elif gt == 2:
        GT = '0/0'

    # Genotype quality: difference between lowest and second lowest Q
    Q_norm = [x - minQ for x in gt_phred]
    Q_norm.sort()

    if math.isinf(Q_norm[1]):
        Q_norm[1] = 999

    GQ = round(Q_norm[1] - Q_norm[0])

    return GT, GQ


def GL(g, ref_errP, alt_errP, m):

    """
    Calculate genotype likelihoods

    Formula 2 in Li 2011, doi:10.1093/bioinformatics/btr509

    m : ploidy
    k : nr of reads
    g : nr of ref alleles
    l : nr of reads supporting ref

    """

    l = len(ref_errP)

    k = l + len(alt_errP)

    L = ( 1 / pow(m, k) ) * \
        np.prod(
            [ ( (m - g) * e ) + ( g * (1 - e) ) for e in ref_errP]
            )  * \
        np.prod(
            [( (m -g) * (1 - e) ) + ( g * e ) for e in alt_errP]
            )

    return L


