#!/usr/bin/env python

""" dettetore -a program to detect and characterize transposable element polymorphisms

Get basic bam file statistics:
    - readlength
    - coverage mean (of core distribution)
    - coverage standard deviation
    - insert size mean
    - insert size standard deviation

Copyright (C) 2018 C. Stritt
License: GNU General Public License v3. See LICENSE.txt for details.
"""

import argparse
import pysam
import statistics
import random

from math import fabs
from Bio import SeqIO


def get_args():

    parser = argparse.ArgumentParser(description='Get basic coverage and \
                                     insert size statistics for a bam file.')

    parser.add_argument("bamfile",
                        help='bwa mem alignment to reference genome.')

    parser.add_argument("reference",
                        help='Reference genome in fasta format.')

    parser.add_argument('-prop', default = 0.001,
                        help='Proportion of data to be used for estimating \
                        statistics.')

    args = parser.parse_args()
    return args


def chromosome_length(fasta):
    chromosomes = {}

    for seq_record in SeqIO.parse(fasta, 'fasta'):
        chromosomes[seq_record.id] = len(seq_record)

    return chromosomes


def core_dist_stats(dist):
    """ Get the core of a distributions and calculates its moments. The core
    rather than the full distribution is used in order to mitigate the
    influence if coverage and insert size outliers in the bam file.
    """
    # Median absolute deviation
    median = statistics.median(dist)
    abs_deviations = [fabs(x - median) for x in dist]
    mad = statistics.median(abs_deviations)

    # Mean and standard deviation of core distribution
    under_limit = median - (10*mad)
    upper_limit = median + (10*mad)

    core_dist = [i for i in dist if (i >= under_limit and i < upper_limit)]

    mean = int(statistics.mean(core_dist))
    stdev = int(statistics.stdev(core_dist))
    return mean, stdev


def coverage_stats(bamfile, refgen):

    chromosomes = chromosome_length(refgen)
    genome_size = 0
    for k in chromosomes:
        genome_size += chromosomes[k]

    pybam = pysam.AlignmentFile(bamfile, "rb")

    # Dictionary with coverage as key
    cov_d = {}
    positions_count = 0

    for pileupcolumn in pybam.pileup(**{"truncate":True}):

        cov = pileupcolumn.nsegments
        positions_count += 1

        if not cov in cov_d:
            cov_d[cov] = 0

        cov_d[cov] += 1
    pybam.close()
    
    # 
    dist =  []
    for k in cov_d:
        dist += cov_d[k] * [k]

    dist += (genome_size-len(dist)) * [0]
    
    if len(dist) < 1e6:
        stats = core_dist_stats(dist)
    else:
        stats = core_dist_stats(random.sample(dist, int(1e6)))
    return stats


def isize_stats(bamfile, proportion):
    """
    Takes a file containing insert sizes and calculates mean and standard
    deviation of the core distribution. Using the core distribution instead of
    all data mitigates the influence of outliers (source: Piccard toolbox)

    """

    pybam = pysam.AlignmentFile(bamfile, "rb")

    isizes = list()
    readlength = int()

    for read in pybam.fetch():

        while not readlength:
            readlength = read.infer_query_length()

        if read.is_proper_pair:
            if random.random() <= proportion:
                isizes.append(abs(read.isize))
    pybam.close()

    stats = core_dist_stats(isizes)

    return stats, readlength


def write_output(bamfile, refgen, proportion):

    cov_mean, cov_stdev = coverage_stats(bamfile, refgen)
    isize, readlength = isize_stats(bamfile, proportion)

    outlist = ['readlength\t' + str(readlength),
               'coverage_mean\t' + str(cov_mean),
               'coverage_stdev\t' + str(cov_stdev),
               'isize_mean\t' + str(isize[0]),
               'isize_stdev\t' + str(isize[1])]

    basename = bamfile.split('/')[-1].split('.')[0]
    outfile = basename + '_bamstats.txt'
    with open(outfile, 'w') as f:
        for line in outlist:
            print(line)
            f.write(line + '\n')
    return outfile


def main():

    args = get_args()
    write_output(args.bamfile, args.reference, args.prop)


if __name__ == '__main__':
    main()
