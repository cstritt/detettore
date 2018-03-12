#!/usr/bin/env python

""" dettetore -a program to detect and characterize transposable element polymorphisms

Filter detettore output

Copyright (C) 2018 C. Stritt
License: GNU General Public License v3. See LICENSE.txt for details.
"""

import argparse
import sys
from scripts import strumenti


def get_args():
    
    parser = argparse.ArgumentParser(description='Apply various filters \
                                     to the detettore raw output')
    
    parser.add_argument("-f", dest='candidates',
                        help='TIP candidates file')
    
    # coverage-based filters
    parser.add_argument('-cov_prop',
                        type=float, 
                        help='Discard results for which the number of \
                        supporting reads on both sides is not at least \
                        the specified proportion of the coverage in the \
                        region, as well as results in extremely high or low \
                        coverage regions')
    
    parser.add_argument('-cov_max', type=int, 
                        help='Maximum coverage in the region')
    
    parser.add_argument('-cov_min', type=int, 
                        help='Minimum coverage in the region')
    
    # consistency filters
    parser.add_argument('-strand',
                        action='store_true',
                        help='Discard results with different strands \
                        inferred from the 5 and the 3 prime site')
    
    parser.add_argument('-id',
                        action='store_true',
                        help='Discard results with different elements \
                         inferred from the 5 and the 3 prime site')
    
    # others
    parser.add_argument('-comb',
                        action='store_true',
                        help='Discard results if not supported by both \
                        discordant pairs and splitreads')
   
    parser.add_argument('-aligned',
                        type=int, 
                        help='Minimum number of base pairs that could be \
                        reconstructed in the target sequence')
    
    parser.add_argument('-nbridge',  type=int, 
                        help='Filter out candidates with more than this number \
                        of reads bridging the insertion point')
    
#    parser.add_argument('-exclude',
#                        help='File containing names of annotated genes that \
#                        should be filtered out, e.g. because they overlap with TEs')

    parser.add_argument('-supp_reads',
                        help='Output a new supporting reads file containing \
                        only the reads for the filtered output')
    
    args = parser.parse_args()
    return args
    

def passare(candidate, filtri):
    """ Filter detettore raw output """

    """ Minimum proportion of reads supporting the insertion relative to
    local read depth
    """
    if filtri['cov_prop']:
        f_supp = candidate.support[0]
        r_supp = candidate.support[1]
        local_cov = float(candidate.coverage)

        if (f_supp/local_cov < filtri['cov_prop']) or \
                (r_supp/local_cov < filtri['cov_prop']):
            return False
        
    if filtri['cov_max']:
        if candidate.coverage > filtri['cov_max']:
            return False
        
    if filtri['cov_min']:
        if candidate.coverage < filtri['cov_min']:
            return False
    

    """ The same strand has to be identfied from forward and reverse
    """
    if filtri['strand']:
        if candidate.strand[0] != candidate.strand[1]:
            return False

    """ The same TE has to be identified from forward and reverse
    """
    TEs = candidate.TE.split('/')
    if filtri['id'] and len(TEs) > 1:
        if TEs[0] != TEs[1]:
            return False

    """ Both discordant read pairs and splitreads have to support the insertion
    """
    if filtri['combined']:
        nr_discordant = sum(candidate.discordant)
        nr_split = sum(candidate.split)
        if nr_discordant == 0 or nr_split == 0:
            return False

    """ Minimum number of base pairs of the identifed TE that can be
    reconstructed
    """
    if filtri['aligned']:
        if candidate.aligned < filtri['aligned']:
            return False
        
    """ Reads bridging the insertion breakpoint
    """
    if filtri['nbridge']:
        if candidate.nbridge > filtri['nbridge']:
            return False

    """ Minimum proportion of reads supporting the insertion relative to
    local read depth
    """
    if filtri['cov_prop'] > 0:
        f_supp = candidate.support[0]
        r_supp = candidate.support[1]
        local_cov = float(candidate.coverage)

        if (f_supp/local_cov < filtri['cov_prop']) or \
                (r_supp/local_cov < filtri['cov_prop']):
            return False

    return True


def main():
    
    args = get_args()
    
#    if args.exclude:
#        with open(args.exclude) as f:
#            exclude = {x.strip() for x in f}

    filtri = {'cov_prop': args.cov_prop,
              'cov_max': args.cov_max,
              'cov_min': args.cov_min,
              'aligned': args.aligned,
              'combined': args.comb,
              'nbridge' : args.nbridge,
              'strand': args.strand,
              'id': args.id}

    indices = []
    with open(args.candidates, 'r') as f:
       
        sys.stdout.write(next(f))
        
        for i, line in enumerate(f):
            x = strumenti.rawvariant(line, 'NA')
            if passare(x, filtri):
                
                indices.append(i)
                sys.stdout.write(line)
                
    if args.supp_reads:
        outfile = open('TIPs_supporting_reads_filt.txt', 'w')
        with open(args.supp_reads) as g:
            for j, line in enumerate(g):
                if j in indices:
                    outfile.write(line)

            
if __name__ == '__main__':
    main()
    