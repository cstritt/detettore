#!/usr/bin/env python

""" dettetore -a program to detect and characterize transposable element polymorphisms

Detect TE absence polymorphisms (TAPs)

Copyright (C) 2018 C. Stritt
License: GNU General Public License v3. See LICENSE.txt for details.
"""

import pysam
import re
import statistics

from Bio import SeqIO
from collections import Counter
from joblib import Parallel, delayed
from scipy.stats import norm
from scripts import strumenti


class annotation:

    def __init__(self, file_path):

        # Check file format
        self.file_path = file_path
        self.format = file_path.split('.')[-1]
        if not self.format.startswith('gff') or self.format.startswith('bed'):
            print('Check annotation format.')

        # Read in annotation
        self.annot = []
        with open(self.file_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split()

                if self.format.startswith('bed'):

                    entry = {
                        'chromosome' : fields[0],
                        'start' : int(fields[1]),
                        'end' :  int(fields[2]),
                        'id' : fields[3],
                        'strand' : fields[5]
                        }

                elif self.format.startswith('gff'):

                    entry = {
                        # GFF uses 1-based indexing!
                        'chromosome' : fields[0],
                        'start' : int(fields[3]),
                        'end' : int(fields[4]),
                        'id' : fields[-1].split(';')[0].split('=')[-1],
                        'strand' : fields[6]
                        }

                self.annot.append(entry)


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
            'mapq_f' : [],
            'mapq_r' : []
            }

    def add_deviant_info(self, deviant_summary):

        self.nr_deviant = len(deviant_summary['name'])
        self.isize = statistics.median(deviant_summary['isizes'])
        try:
            self.start = max(deviant_summary['f_strand_ends'])
            self.end = min(deviant_summary['r_strand_starts'])
        except ValueError:
            pass

    #def add_reference_support(self):



def extract_deviant_reads(region, isize_stats, bamfile, mapq):

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


def process_taps(
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


    deviant = extract_deviant_reads(region, isize_stats, bamfile, uniq)

    if not deviant['name']:
        return None

    cov_mean, cov_stdev = strumenti.local_cov(region, bamfile)

    summary = TAP_candidate([
        i,
        feature_length,
        region_start,
        region_end,
        cov_mean
        ])


    if deviant['name']:
        summary.add_deviant_info(deviant)

    return summary


#%%

def run_module(
        bamfile,
        readinfo,
        annotation_file,
        reference,
        thresholds,
        cpus):

    """ Read in annotation and read information
    """

    reference_TEs = annotation(annotation_file)

    with open(readinfo) as f:
        readlength = int(next(f).split('\t')[1])
        next(f)
        next(f)
        isize_mean = int(next(f).split('\t')[1])
        isize_stdev = int(next(f).split('\t')[1])

    isize_stats = (isize_mean, isize_stdev, readlength)
    uniq = thresholds[0]

    contigs = {seq_record.id: seq_record.seq for seq_record in SeqIO.parse(reference, "fasta")}

    print('Extracting read pairs with deviant insert sizes around annotated features\n...')
    inputs = [(i, x) for i, x  in enumerate(reference_TEs.annot)]
    TAP_candidates = Parallel(n_jobs=cpus)(delayed(process_taps)\
                                           (i,
                                            isize_stats,
                                            bamfile,
                                            contigs,
                                            uniq) for i in inputs)


    #%% write vcf output

    # Exact position from splitreads: strumenti line 956ff

    for candidate in out:

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









