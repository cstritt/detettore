#!/usr/bin/env python

""" dettetore -a program to detect and characterize transposable element polymorphisms

Call TE variants and create vcf file

Copyright (C) 2018 C. Stritt
License: GNU General Public License v3. See LICENSE.txt for details.
"""

import argparse
import copy
from scripts import strumenti
import time


def get_args():

    parser = argparse.ArgumentParser(description='Create vcf from multiple TIP candidate files.')

    parser.add_argument('-f', dest='rawfiles',
                        help='File with paths to TIP candidates files.')

    parser.add_argument('-w ', type=int, default=50, dest='wsize',
                        help='Insertions within wsize are considered identical \
                        if a) the TE has no splitread evidence and b) the same \
                        TE is identified on the same strand.')


    # coverage-based filtering
    parser.add_argument('-cov_prop', type=float, default=None,
                        help='Discard results for which the number of \
                        supporting reads on both sides is not at least \
                        the specified proportion of the coverage in the \
                        regions')

    parser.add_argument('-cov_max', type=int,
                        help='Maximum coverage in the region')

    parser.add_argument('-cov_min', type=int,
                        help='Minimum coverage in the region')

    # other filters
    parser.add_argument('-aligned', type=int,
                        help='Minimum number of base pairs that could be \
                        reconstructed in the target sequence')

    parser.add_argument('-het',  type=int, dest='nbridge',
                        help='Minimum number of bridge reads for an \
                        insertion to be considered heterozygous')

    parser.add_argument('-comb', action='store_true',
                        help='Discard results if not supported by both \
                        discordant pairs and splitreads')


    # consistency-based filtering
    parser.add_argument('-strand', action='store_true',
                        help='Discard results with different strands \
                        inferred from the 5 and the 3 prime site')

    parser.add_argument('-id', action='store_true',
                        help='Discard results with different elements \
                         inferred from the 5 and the 3 prime site')

    args = parser.parse_args()
    return args


def passare(candidate, filtri):

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

    """ Region coverage minumum and maximum
    """
    if filtri['cov_max']:
        if candidate.coverage > filtri['cov_max']:
            return False

    if filtri['cov_min']:
        if candidate.coverage < filtri['cov_max']:
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


def add_sample(path, d, accs, filtri):

    """ Create dictionary[chromosome][position] containing all TIPs which pass
    the filtering criteria
    """

    accession = path.split('/')[-2]
    accs.append(accession)
    cnt = 0
    with open(path) as f:

        next(f)
        for line in f:

            x = strumenti.rawvariant(line, accession)

            if passare(x, filtri):

                if x.chr not in d:
                    d[x.chr] = dict()

                if x.pos not in d[x.chr]:
                    d[x.chr][x.pos] = []

                d[x.chr][x.pos].append(x)
                cnt += 1

        print '%i TIPs retained after filtering in %s' % (cnt, accession)
        return d, accs


def cluster_TIP_positions(d, window_size):
    """ Cluster positions
    """
    window_size = 50

    clusters = {}
    for chrmsm in d:
        sorted_positions = sorted(d[chrmsm].keys())
        previous = 0-window_size

        clusters[chrmsm] = []

        for p in sorted_positions:
            if p in range(previous, previous + window_size):
                clusters[chrmsm][-1].append(p)
            else:
                clusters[chrmsm].append([p])
            previous = p

    cluster_sizes = {}
    for chrmsm in clusters:
        for cl in clusters[chrmsm]:
            size = len(cl)
            if size not in cluster_sizes:
                cluster_sizes[size] = 0
            cluster_sizes[size] += 1

    print('\n# of TIPs in a window: # of windows with that many TIPs')
    for k in sorted(cluster_sizes.keys()):
        print str(k) + ': ' + str(cluster_sizes[k])

    return clusters


def get_unique_insertion_sites(clusters, d):

    unique_insertions = {}

    for chrmsm in clusters:

        unique_insertions[chrmsm] = {}

        for cl in clusters[chrmsm]:

            """ Within clusters, get candidates with positions known from
            splitreads and candidates with ambiguous positions
            """
            uniq =  {}
            ambiguous =  {}

            for pos in cl:
                for variant in d[chrmsm][pos]:

                    TE = variant.TE
                    var_id = (pos, TE)

                    if variant.tsd:
                        if var_id not in uniq:
                            uniq[var_id] = [[], []]
                        uniq[var_id][0].append(variant)

                    else:
                        if var_id not in ambiguous:
                            ambiguous[var_id] = []
                        ambiguous[var_id].append(variant)

            """ Go through ambiguous insertions and asign them to the closest
            known position of the same TE family
            """
            for var_id in ambiguous:

                uniq_k = uniq.keys()
                uniq_pos = [x[0] for x in uniq_k]
                uniq_TEs = [x[1] for x in uniq_k]

                pos, TE = var_id[0], var_id[1]

                if TE in uniq_TEs:

                    TE_positions = [uniq_pos[i] for i, x in enumerate(uniq_TEs)
                                    if x == TE]

                    closest = min(TE_positions, key=lambda x:abs(x-pos))
                    closest_i = uniq_pos.index(closest)
                    uniq[uniq_k[closest_i]][1] += ambiguous[var_id]

                #  If there is no identical TE among the known positions, create a new insertion event
                else:
                    uniq[var_id] = [[], ambiguous[var_id]]

            unique_insertions[chrmsm].update(uniq)

    return unique_insertions


def write_vcf(accs, unique_insertions, nbridge):

    """ TE vcf format. Info fields:
        GT: genotype
        TE: TE family name
        STR: strand
        DR: number of supporting discordant reads
        SR: number of supporting split reads
        TSD: target site duplication
        DP: Total depth
        AD: Average depth
        DI: Distance to precise insertion location, if known
    """

    date = time.strftime("%d/%m/%Y")

    accs_d_mc = {a:'' for a in accs}
    metainfo = ['##fileFormat=VCFv4.2',
                '##fileFormat=(%s)' % date,
                '##source==detettore v0.6',

                '##INFO=<ID=TE,Number=,Type=String,Description="TE family name">',

                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=ST,Number=2,Type=String,Description="Strand">',
                '##FORMAT=<ID=DP,Number=2,Type=Integer,Description="Supporting reads">',
                '##FORMAT=<ID=DR,Number=2,Type=Integer,Description="Discordant reads">',
                '##FORMAT=<ID=SR,Number=2,Type=Integer,Description="Split reads">',
                '##FORMAT=<ID=AL,Number=1,Type=Integer,Description="TE alignment length">',
                '##FORMAT=<ID=TSD,Number=1,Type=Character,Description="Target site duplication">',
                '##FORMAT=<ID=DI,Number=1,Type=Integer,Description="Distance from insertion site">']


    header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT',
              'QUAL', 'FILTER', 'INFO', 'FORMAT'] + accs

    vcf = []

    for chrmsm in unique_insertions.keys():
        for site in unique_insertions[chrmsm]:

            vcfline = [chrmsm,
                       str(site[0]),
                       '.',
                       '.',
                       site[1],
                       '.',
                       '.',
                       '.',
                       'GT:ST:DP:DR:SR:AL:TSD:DI']

            accs_d = copy.deepcopy(accs_d_mc)
            pos = site[0]
            uniq, ambiguous = unique_insertions[chrmsm][site]

            for variant in uniq + ambiguous:

                """ Write heterozygous variant if there are more than
                het_threshold reads bridging the insertion breakpoint
                """
                if variant.nbridge > nbridge:
                    gt = '1/0'
                else:
                    gt = '1/1'

                if uniq:
                    di = abs(pos-uniq[0].pos)
                else:
                    di = 0

                if not variant.tsd:
                    variant.tsd = '.'

                acc_info = [gt,
                            ','.join(variant.strand),
                            ','.join(map(str, variant.support)),
                            ','.join(map(str, variant.discordant)),
                            ','.join(map(str, variant.split)),
                            str(variant.aligned),
                            variant.tsd,
                            str(di)]

                accs_d[variant.source] = ':'.join(acc_info)

            for a in accs:
                if not accs_d[a]:
                    vcfline.append('.')
                else:
                    vcfline.append(accs_d[a])

            vcf.append(vcfline)

    vcf = sorted(vcf, key=lambda x: (x[0], int(x[1])))
    with open('TIPs.vcf', 'w') as f:
        f.write('\n'.join(metainfo) + '\n')
        f.write('\t'.join(header) + '\n')

        for line in vcf:
            f.write('\t'.join(line) + '\n')


def main():

    args = get_args()

    with open(args.rawfiles) as f:
        files = [x.strip() for x in f]

    filtri = {'cov_prop': args.cov_prop,
              'cov_max': args.cov_max,
              'cov_min': args.cov_min,
              'aligned': args.aligned,
              'combined': args.comb,
              'strand': args.strand,
              'id': args.id}

    d = {}
    accs = []

    for f in files:
        d, accs = add_sample(f, d, accs, filtri)


    clusters = cluster_TIP_positions(d, args.wsize)
    unique_insertions = get_unique_insertion_sites(clusters, d)
    write_vcf(accs, unique_insertions, args.nbridge)


if __name__ == '__main__':
    main()
