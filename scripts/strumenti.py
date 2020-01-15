#!/usr/bin/env python

""" dettetore -a program to detect and characterize transposable element polymorphisms

Toolbox for TE polymorphism detection

Copyright (C) 2018 C. Stritt
License: GNU General Public License v3. See LICENSE.txt for details.
"""

import os
import scipy
import statistics
import subprocess
import pysam

from math import fabs
from random import random
from collections import Counter


class hit:
    """ Container for the blastn output,
    including the query sequence aligned to the TE
    """

    def __init__(self, fields):
        self.read = fields[0]
        self.aln_len = float(fields[1])
        self.score = float(fields[2])
        self.target_strand = fields[3]
        self.target_id = fields[4]
        self.target_aln_start = int(fields[5])
        self.target_aln_end = int(fields[6])
        self.seq = fields[7]


class cluster_summary:

    """ In this container the clustered anchor reads
    and the mates aligning to TEs are combined.
    """

    def __init__(self):

        self.discordant = self.attributes()
        self.splitreads = self.attributes()
        self.score = int()
        self.strand = {'plus': 0, 'minus': 0}
        self.aligned_positions = set()
        self.region_start = float('inf')
        self.region_end = -float('inf')

    class attributes:

        def __init__(self):

            self.start = []
            self.end = []
            self.anchors = []

    def update(self, hit, anchor, readname, modus, weight):

        self.score += (weight * hit.score)
        self.strand[hit.target_strand] += 1

        aln_start = hit.target_aln_start
        aln_end = hit.target_aln_end

        if aln_start < aln_end:
            self.aligned_positions.update(range(aln_start, aln_end))
        else:
            self.aligned_positions.update(range(aln_end, aln_start))


        if anchor.reference_start < self.region_start:
            self.region_start = anchor.reference_start
        if anchor.reference_end > self.region_end:
            self.region_end = anchor.reference_end

        if modus == 'discordant':
            self.discordant.start.append(anchor.reference_start)
            self.discordant.end.append(anchor.reference_end)
            self.discordant.anchors.append(readname)

        elif modus == 'splitreads':
            self.splitreads.start.append(anchor.reference_start)
            self.splitreads.end.append(anchor.reference_end)
            self.splitreads.anchors.append(readname)


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


""" FUNZIONI """


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
    q_1 = scipy.percentile(lista, 25)
    q_3 = scipy.percentile(lista, 75)

    boxlength = q_3 - q_1

    lower = max(min(lista), q_1 - 1.5*boxlength)
    upper = min(max(lista), q_3 + 1.5*boxlength)

    lista_filt = [x for x in lista if x >= lower and x <= upper]
    return lista_filt


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


def parse_gff(gff_file):
    gff = dict()
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            chrmsm = fields[0]
            try:
                gff[chrmsm].append(fields)
            except KeyError:
                gff[chrmsm] = [fields]
    return gff


def is_overlapping(a,b):
    a_range = set(range(a[0], a[1]))
    b_range = set(range(b[0], b[1]))

    intersection = a_range & b_range
    union = a_range | b_range

    if len(intersection) == 0:
        return False
    else:
        return len(intersection) / float(len(union))


def cluster_reads(positions, overlap_proportion):

    """ Takes a dictionary with chromosomes as keys and ordered intervals.
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

        for i in range(len(positions_chrmsm)):

            line = positions_chrmsm[i]
            name = line[0]
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

    return clusters


def prefilter_clusters(clusters, n):
    """ Clusters with less than n reads are filtered out"""

    clusters_filt = dict()
    for k in clusters:
        for cl in clusters[k]:

            if len(cl) < n:
                continue
            else:
                try:
                    clusters_filt[k].append(cl)
                except KeyError:
                    clusters_filt[k] = [cl]
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
           '-max_target_seqs', '1',
           '-num_threads', str(cpus)]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output = proc.stdout.read()
    blast_out = output.splitlines()

    with open(outfile, 'w') as f:
        for line in blast_out:
            f.write(line.decode('utf-8') + '\n')

    return blast_out


def hit_dictionary(blast_out, min_aln_len, modus):

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
            TEhit = hit(fields)

            if name in TEhits:
                if modus == 'one_hit_only':

                    previous_score = TEhits[name][0].score
                    new_score = TEhit.score

                    if new_score == previous_score:
                        TEhits[name].append(TEhit)
                    elif new_score < previous_score:
                        continue
                    elif new_score > previous_score:
                        TEhits[name] = [TEhit]

                elif modus == 'multiple_hits':
                    TEhits[name].append(TEhit)
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


def summarize_cluster(cluster, read_dictionary, hits, modus, weight):

    summary = {}

    for readname in cluster:

        if readname not in hits:
            continue

        read = read_dictionary[readname]
        cigar = read.cigartuples

        if modus == 'discordant':
            if read.is_reverse:
                s = 1
            else:
                s=0

        elif modus == 'splitreads':
            if cigar[0][0] == 4:
                s = 1
            elif cigar[-1][0] == 4:
                s = 0

        for hit in hits[readname]:

            targetname = hit.target_id

            if targetname not in summary:
                summary[targetname] = [cluster_summary(), cluster_summary()]

            summary[targetname][s].update(hit, read, readname, modus, weight)

    targets = list(summary.keys())
    deletenda = []

    for target in targets:
        if modus == 'discordant':
            nr_obs_l = len(summary[target][0].discordant.anchors)
            nr_obs_r = len(summary[target][1].discordant.anchors)

        elif modus == 'splitreads':
            nr_obs_l = len(summary[target][0].splitreads.anchors)
            nr_obs_r = len(summary[target][1].splitreads.anchors)

        if not (nr_obs_l and nr_obs_r):
            deletenda.append(target)

    for target in deletenda:
        del summary[target]

    return summary


def combine_hits_and_anchors(clusters,
                             anchors,
                             hits):

    """ Takes the clustered anchors of the discordant read pairs, checks if
    their mates map to a TE, and, if yes, creates a summary of the cluster in
    form of two cluster_summary objects: one for the evidence from the 5' side,
    one for the evidence from the 3' side.

    """

    summaries = dict()
    posizioni = dict()

    chromosomi = sorted(clusters.keys())

    for chromosome in chromosomi:

        for cluster in clusters[chromosome]:

            summary = summarize_cluster(cluster, anchors, hits, 'discordant', 1)
            if not summary:
                continue

            region_start = min([summary[x][0].region_start for x in summary])
            region_end = max([summary[x][1].region_end for x in summary])

            try:
                summaries[chromosome].append(summary)
                posizioni[chromosome].append((region_start, region_end))

            except KeyError:
                summaries[chromosome] = [summary]
                posizioni[chromosome] = [(region_start, region_end)]

    return summaries, posizioni


def add_splitread_evidence(cluster_summaries,
                           split_clusters,
                           splitreads,
                           split_hits,
                           discordant_cluster_positions):

    """ Add split read information to existing discordant read summaries
    or create new summaries if the split read cluster does not overlap with
    any discordant-read summary
    """

    chromosomi = sorted(split_clusters.keys())

    for chromosome in chromosomi:
        x = 0 # variable to store the position in the loop below, to avoid going through the whole thing every time
        for cluster in split_clusters[chromosome]:

            summary = summarize_cluster(cluster, splitreads, split_hits, 'splitreads', 2)

            if not summary:
                continue

            else:

                clipping_positions = []
                for y in summary:
                    clipping_positions += summary[y][0].splitreads.end
                position = max(clipping_positions)

                overlap = False

                if chromosome in discordant_cluster_positions:

                    """ Go through discordant read summaries and check
                    if there is positional overlap.
                    """
                    for i in range(x, len(discordant_cluster_positions[chromosome])):

                        region_start = discordant_cluster_positions[chromosome][i][0]
                        region_end = discordant_cluster_positions[chromosome][i][1]

                        if position in range(region_start, region_end):

                            cluster_summaries[chromosome][i] = merge_summaries(
                                    cluster_summaries[chromosome][i],
                                    summary)
                            x = i
                            overlap = True

                    if not overlap:
                        cluster_summaries[chromosome].append(summary)

                else:
                    cluster_summaries[chromosome] = [summary]

    return cluster_summaries


def merge_summaries(existing, new):
    """ Merge discordant read and split read summaries
    """

    for target in new:

        if not target in existing:
            existing[target] = [cluster_summary(), cluster_summary()]

        # 0 and 1 for the forward and the reverse strand, resp.
        for i in [0,1]:

            existing[target][i].splitreads.start = new[target][i].splitreads.start
            existing[target][i].splitreads.end = new[target][i].splitreads.end
            existing[target][i].splitreads.anchors = new[target][i].splitreads.anchors

            existing[target][i].score += new[target][i].score
            existing[target][i].aligned_positions.update(new[target][i].aligned_positions)

            if new[target][i].region_start < existing[target][i].region_start:
                existing[target][i].region_start = new[target][i].region_start

            if new[target][i].region_end > existing[target][i].region_end:
                existing[target][i].region_end = new[target][i].region_end

    return existing


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


def get_highest_scoring_target(site):

    """ Returns the TE family with the highest cumulative blast bit score,
    separately for the 5' and the 3' site of the insertion.
    """

    highest_scoring = []

    for i in [0,1]:

        highscore = 0

        for target in site:

            if site[target][i].score > highscore:
                highscore = site[target][i].score
                best = target

        highest_scoring.append(best)

    return highest_scoring


def overlapping_reads(chromosome, position, bamfile, overlap):

    """ Checks if there are reads bridging the suspected insertion
    break point.
    """

    bridge_reads = []

    pybam = pysam.AlignmentFile(bamfile, "rb")

    for read in pybam.fetch(chromosome, position, position+1):

        if not read.mapq > 0:
            continue

        down = [x for x in range(read.reference_start, read.reference_end) if
                x < position]
        up = [x for x in range(read.reference_start, read.reference_end) if
                x > position]

        if len(down) > overlap and len(up) > overlap:

            bridge_reads.append(read.query_name)

    pybam.close()
    return len(bridge_reads)


def create_table(combined, bamfile):

    header = ['chromosome',
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

            outline = [chrmsm,
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
                       anchor_reads]

            outline = map(str, outline)
            tabula.append(list(outline))

    tabula = sorted(tabula, key=lambda x: (x[0], int(x[1])))
    with open('TIPs_candidates.tsv', 'w') as f, \
         open('TIPs_supporting_reads.txt', 'w') as g:

        f.write('\t'.join(header) + '\n')
        g.write('discordant_left;discordant_right;split_left;split_right\n')
        for line in tabula:
            f.write('\t'.join(line[:-1]) + '\n')
            g.write(line[-1] + '\n')
