#!/usr/bin/env python

""" dettetore -a program to detect and characterize transposable element polymorphisms

Get basic bam file statistics:
    - readlength
    - coverage mean (of core distribution)
    - coverage standard deviation
    - insert size mean
    - insert size standard deviation

Copyright (C) 2021 C. Stritt
License: GNU General Public License v3. See LICENSE.txt for details.
"""

import argparse
import pysam
import statistics
import random
import sys

from math import fabs


def get_args():

    parser = argparse.ArgumentParser(
        description='Get basic coverage and insert size statistics for a bam file.')

    parser.add_argument('bamfile')

    args = parser.parse_args()
    return args


def local_cov(region, bamfile):
    """ 
    Estimates mean coverage and standard deviation in the region 
    chromosome:start-end. Attencion: indexing in pybam is 0-based 
    """

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
        c = depth[i] if i in depth else 0
        coverage.append(c)

    if len(coverage) > 1:
        mean = int(statistics.mean(coverage))
        stdev = int(statistics.stdev(coverage))
    elif len(coverage) == 1:
        mean = int(coverage[0])
        stdev = 'NA'
    else:
        mean = 0
        stdev = 'NA'

    return mean, stdev


def core_dist_stats(dist):
    """
    Get the core of a distributions and calculate its moments. The core
    rather than the full distribution is used in order to mitigate the
    influence if coverage and insert size outliers in the bam file.
    """

    # Overall mean (for comparison...)
    mean_tutto = int(statistics.mean(dist))

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
    
    return mean_tutto, mean, stdev


def coverage_stats(bamfile):
    """
    Create list with coverage at every single position. If the list is longer than
    1e6, use random sample to estimate coverage statistics.

    Only traverses the first chromosome if there are multiple.

    """

    pybam = pysam.AlignmentFile(bamfile, "rb")
    first_chromosome = pybam.get_reference_name(0)

    cov = [pileupcolumn.nsegments for pileupcolumn in pybam.pileup(reference=first_chromosome, **{"truncate":True})]
    pybam.close()

    if len(cov) < 1e6:
        mean, mean_cd, stdev_cd = core_dist_stats(cov)
    else:
        mean, mean_cd, stdev_cd = core_dist_stats(random.sample(cov, int(1e6)))
        
    return {
        'cov_mean' : mean,
        'cov_mean_cd' : mean_cd,
        'cov_stdev_cd' : stdev_cd
        }


def length_stats(bamfile):
    """
    Get insert sizes of properly paired reads, calculate mean and standard
    deviation of the core distribution. 

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
    
    # Readlengths    
    readlengths = [x for x in readlengths if x]
    
    read_d = {
        'readlength_max' : max(readlengths),
        'readlength_mean' : int(statistics.mean(readlengths)),
        'readlength_stdev' : int(statistics.stdev(readlengths))
        } 
        
    if not isizes:
        return read_d
        
    # Insert sizes
    if len(isizes) < 1e6:
        mean, mean_cd, stdev_cd = core_dist_stats(isizes)
    else:
        mean, mean_cd, stdev_cd = core_dist_stats(random.sample(isizes, int(1e6)))
        
    isize_d = {
        'isize_mean' : mean,
        'isize_mean_cd' : mean_cd,
        'isize_stdev_cd' : stdev_cd
        }    
    
    return read_d, isize_d


def return_dict(bamfile):
    """
    Returns dictionary used  by detettore. Coverage stats disabled, as 
    they require second traversal of bam file.
    
    """

    length_statistiks = length_stats(bamfile)
    
    # No paired end reads
    if isinstance(length_statistiks, dict):
        
        return { 
            'readlength' : length_statistiks['readlength_max']
            }
         
    else:
        read_d, isize_d = length_statistiks
    
        return {
            'readlength' : read_d['readlength_max'],
            'isize_mean' : isize_d['isize_mean_cd'],
            'isize_stdev' : isize_d['isize_stdev_cd']
            }
    
    
def main():
    """
    Write to stdout:
        
    cov_mean, cov_mean_cd, cov_stdev_cd, readlength, PE/SE, isize_mean, isize_stdev
    
    """
    
    args = get_args()
    
    cov_d = coverage_stats(args.bamfile)
    length_statistiks = length_stats(args.bamfile)
    
    if isinstance(length_statistiks, dict):
        rtype = 'SE'
        read_d = length_statistiks
        isize_d = {
            'isize_mean_cd' : 'NA',
            'isize_stdev_cd' : 'NA'
            }
        
    else:
        rtype = 'PE'
        read_d, isize_d = length_statistiks
        
    outline = map( str, [
        cov_d['cov_mean'], 
        cov_d['cov_mean_cd'], 
        cov_d['cov_stdev_cd'], 
        read_d['readlength_max'],
        read_d['readlength_mean'],
        read_d['readlength_stdev'],
        rtype, 
        isize_d['isize_mean_cd'],
        isize_d['isize_stdev_cd']
        ] )
    
    
    sys.stdout.write('\t'.join([
        'cov_mean', 'cov_mean_cd', 'cov_stdev_cd', 'readlength_max', 
        'readlength_mean', 'readlength_stdev', 'read_type', 'isize_mean_cd', 'isize_stdev_cd']) +'\n')
    sys.stdout.write('\t'.join(outline)+'\n')
    
    
if __name__ == '__main__':
    main()

