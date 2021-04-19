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


def get_args():

    parser = argparse.ArgumentParser(
        description='Get basic coverage and insert size statistics for a bam file.')

    parser.add_argument("bamfile",
                        help='bwa mem alignment to reference genome.')

    args = parser.parse_args()
    return args


def local_cov(region, bamfile):
    """ Estimates mean coverage and standard deviation in the
    region chromosome:start-end. Attencion: indexing in pybam is 0-based """

    chromosome, start, end = region
    pybam = pysam.AlignmentFile(bamfile, "rb")
    depth = dict()

    for pileupcolumn in pybam.pileup(chromosome, start, end, **{"truncate":True}):

        d = pileupcolumn.nsegments
        pos = pileupcolumn.reference_pos
        depth[pos] = d

    pybam.close()

    coverage = []
    for i in range(start, end):
        try:
            c = depth[i]
        except KeyError:
            c = 0
        coverage.append(c)

    if len(coverage) > 1:
        mean = int(statistics.mean(coverage))
        stdev = int(statistics.stdev(coverage))
    else:
        mean = int(coverage[0])
        stdev = 'NA'

    return mean, stdev


def local_cov_filt(region, bamfile, filters):

    """ Estimate mean coverage and standard deviation in the region of
    interest. Filters for base and mapping quality are applied
    """
    chromosome, start, end = region
    baseq, mapq = filters

    pybam = pysam.AlignmentFile(bamfile, "rb")
    cov_dict = {}

    crap_reads = set()

    for pileupcolumn in pybam.pileup(chromosome, start, end,
                                     **{"truncate": True}):
        pos = pileupcolumn.pos
        cov_dict[pos] = 0

        for pileupread in pileupcolumn.pileups:

            read = pileupread.alignment

            read_id = read.query_name
            if read_id in crap_reads:
                continue

            if pileupread.is_del or pileupread.is_refskip or\
                not read.is_proper_pair or read.mapping_quality < mapq:

                crap_reads.add(read_id)
                continue

            if read.query_qualities[pileupread.query_position] < baseq:
                continue

            cov_dict[pos] += 1

    coverage = []
    for i in range(start, end):
        try:
            c = cov_dict[i]
        except KeyError:
            c = 0
        coverage.append(c)

    if len(coverage) > 1:
        mean = int(statistics.mean(coverage))
        stdev = int(statistics.stdev(coverage))
    else:
        mean = int(coverage[0])
        stdev = 'NA'

    return mean, stdev


def core_dist_stats(dist):
    """
    Get the core of a distributions and calculate its moments. The core
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


def coverage_stats(bamfile):
    """
    Create list with coverage at every single position. If the list is longer than
    1e6, use random sample to estimate coverage statistics.

    To do: only use first chromosome, traverse bam only once for coverage and isize stats
    count_coverage function.

    Better: only use 1e6 bp

    """

    pybam = pysam.AlignmentFile(bamfile, "rb")
    first_chromosome = pybam.get_reference_name(0)

    cov = [pileupcolumn.nsegments for pileupcolumn in pybam.pileup(reference=first_chromosome, **{"truncate":True})]
    pybam.close()

    if len(cov) < 1e6:
        stats = core_dist_stats(cov)
    else:
        stats = core_dist_stats(random.sample(cov, int(1e6)))
    return stats


def isize_stats(bamfile):
    """
    Get insert sizes of properly paired reads, calculate mean and standard
    deviation of the core distribution. Using the core distribution instead of
    all data mitigates the influence of outliers (source: Piccard toolbox).

    A second traversal of the bam file is required because here reads are traversed rather
    than positions, as in the coverage_stats function.

    """

    pybam = pysam.AlignmentFile(bamfile, "rb")
    first_chromosome = pybam.get_reference_name(0)

    isizes = list()
    readlengths = set()

    for read in pybam.fetch(reference=first_chromosome):

        readlengths.add(read.infer_query_length())

        if read.is_proper_pair:
            isizes.append(abs(read.isize))

    pybam.close()

    if len(isizes) < 1e6:
        stats = core_dist_stats(isizes)
    else:
        stats = core_dist_stats(random.sample(isizes, int(1e6)))

    return stats, max([x for x in readlengths if x])


def write_output(bamfile, outfmt):

    # cov_mean, cov_stdev = coverage_stats(bamfile)
    isize, readlength = isize_stats(bamfile)

    outlist = [
        'readlength\t' + str(readlength),
        #'coverage_mean\t' + str(cov_mean),
        #'coverage_stdev\t' + str(cov_stdev),
        'isize_mean\t' + str(isize[0]),
        'isize_stdev\t' + str(isize[1])
        ]

    if outfmt == 'file':
        basename = bamfile.split('/')[-1].split('.')[0]
        outfile = basename + '_bamstats.txt'
        with open(outfile, 'w') as f:

            for line in outlist:
                print(line)
                f.write(line + '\n')
        return outfile

    elif outfmt == 'dict':

        return {
            'readlength' : readlength,
            #'coverage_mean' : cov_mean,
            #'coverage_stdev' : cov_stdev,
            'isize_mean' : isize[0],
            'isize_stdev' : isize[1]
            }


def main():

    args = get_args()
    bamfile = args.bamfile
    write_output(bamfile)

if __name__ == '__main__':
    main()
