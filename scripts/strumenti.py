#!/usr/bin/env python

""" dettetore - a program to detect and characterize transposable element polymorphisms

Toolbox for TE polymorphism detection

Copyright (C) 2021 C. Stritt
License: GNU General Public License v3. See LICENSE.txt for details.


Issues:
    -l. 509: GQ used instead of QUAL; use '.' for the moment
        
    - TE annotation: require bed format, with name in col5 and (optionally)
    family in col6. If family is not provided, it is inferred from the consensus library
    
    

"""

import os
import subprocess
import pysam
import sys
import numpy as np
import math
import re
import scipy
import copy
import multiprocessing as mp

from scripts import bamstats

from Bio import SeqIO
from random import sample
from collections import Counter
from numba import prange
from scipy.stats import ttest_1samp



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
        self.outname = args.outname

        # Thresholds
        self.mapq = args.mapq
        self.aln_len_DR = args.aln_len_DR
        self.aln_len_SR = args.aln_len_SR

        # Program settings
        self.modus = args.modus
        self.cpus = args.cpus
        self.include_invariant = args.include_invariant     
        self.require_split = args.require_split
        
        # Parse fasta files
        self.ref_contigs = {
            seq_record.id:
                seq_record.seq for seq_record in SeqIO.parse(self.reference, "fasta")}
            
        self.telib = {
            seq_record.id:
                seq_record.seq for seq_record in SeqIO.parse(self.targets, "fasta")}
            
        # Region
        self.region = args.region.split(':') if args.region else None
        
        if self.region and len(self.region)==2:
            self.region[1] = int(self.region[1])
            self.region[2] = int(self.region[2])

        elif self.region and len(self.region)==1:
            self.region.append(0)
            self.region.append(len(self.ref_contigs[self.region[0]]))
            
            
        # Load reference TE annotation
        self.annotation = []

        if 'taps' in self.modus:

            # Check file format
            formt = self.annotation_file.split('.')[-1]

            # Maximum TE size
            self.longest = 0

            with open(self.annotation_file) as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    
                    if not formt.startswith('bed') and not formt.startswith('gff'):
                        sys.exit('Check annotation format.')
                        
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
                            'id' : re.search(r'.*?=(.*?)[,;]', fields[-1]).group(1),
                            'strand' : fields[6]
                            }
                        
                    if entry['end'] -  entry['start'] > self.longest:
                        self.longest = entry['end'] -  entry['start']

                    self.annotation.append(entry)

        # Estimate read statistics
        print('Estimating bam read statistics\n...')
        readinfo = bamstats.return_dict(self.bamfile)
        for k in readinfo:
                print(k + ' : ' + str(readinfo[k]))   
        
        # No paired-end reads!
        if len(readinfo) == 1:
            self.paired_end = False
            self.readlength = readinfo['readlength']
            self.isize_mean = 0
            self.isize_stdev = 0            

            print('\nNo paired-end reads discovered! Only TIPs based on split reads will be reported.\nTAPs will be scored as missing data (./.)')
            
        else:
            self.paired_end = True
            self.readlength = readinfo['readlength']
            self.isize_mean = readinfo['isize_mean']
            self.isize_stdev = readinfo['isize_stdev']
    

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
        1) Map anchor mates against TE sequences and 2) find clusters of anchor reads
        """

        print('Aligning discordant reads to target sequences...')
        
        DR_hits = minimap(
            'discordant.fasta', 
            parameters, 
            parameters.aln_len_DR,
            k=15,
            w=11,
            )
        
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
            min_cl_size=2
            )

        # Summarize clusters
        self.DR_clusters = summarize_clusters(
            DR_clusters,
            DR_anchors,
            DR_hits,
            'discordant')


    def splitreads(self, parameters, splitreads, split_positions):
        """
        1) Find clusters of splitreads and 2) map clipped parts of
        clustering reads against TE consensus sequences. This is the opposite 
        order than for clustering discordant read pairs, as in this way the 
        time for the mapping is greatly reduced.
        """

        print('\nClustering split reads...')
        
        SR_clusters = cluster_reads(
            split_positions,
            'splitreads',
            overlap_proportion=0.05,
            min_cl_size=2
            )

        print('Aligning split parts to target sequences...')
        
        write_clipped_to_fasta(
            SR_clusters,
            splitreads,
            parameters.aln_len_SR)

        SR_minimapd = minimap(
            'softclipped.fasta', 
            parameters,
            parameters.aln_len_SR,
            k=9,
            w=5
            )
        
        print('%i soft-clipped read parts map to a TE' % (len(SR_minimapd)))

        self.SR_clusters = summarize_clusters(
            SR_clusters,
            splitreads,
            SR_minimapd,
            'splitreads')


    def combineDR_SR(self, parameters):
        """
        Combine discordant read pair and splitread clusters. This creates a 
        list with entries [chromosome, position, DR cluster, SR cluster], 
        with 'False' if only a DR or SR cluster is present at a site

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

        """ 
        Get all the information required in the vcf, as listed below.  
        
        Positions in vcf are 1-based.
        
        Get rid of some noise: use ibpiqr to filter out spread out clusters
            
        """

        vcf_lines = []
        
        stats = {
            "0/1" : 0,
            "1/1" : 0,
            "TEs" : {}
            }

        for site in self.candidates:

            CHROM = site[0]
            POS = site[1]
            REF = parameters.ref_contigs[CHROM][POS - 1]
            
            discordant = site[2]
            split = site[3]
            
            if parameters.require_split and not split:
                continue


            """ Dispersal of read breakpoints in cluster (interquartile range). 
            Large IR expected for spurious cluster, while for a true TIP reakpoints
            should be centered around the insertion site.
            
            Filter out candidates with a spread larger than twice the readlength.
            
            """
            breakpoints = []
            if discordant:
                breakpoints += discordant.breakpoint[0] + discordant.breakpoint[1]
            if split:
                breakpoints += split.breakpoint[0] + split.breakpoint[1]
            
            q1 = scipy.percentile(breakpoints , 25)
            q3 = scipy.percentile(breakpoints , 75)
            IQR = q3 - q1
            
            if IQR > 2*parameters.readlength:
                continue


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
            
            ALT = '<INS:ME:%s>' % (te)
            
            # At least two supporting reads from either side
            if te_hits[te]['side'][0] < 2 or  te_hits[te]['side'][1] < 2:
                continue
                        
            te_start = min(te_hits[te]['aligned_positions'])
            te_end = max(te_hits[te]['aligned_positions'])
            te_strand = '+' if te_hits[te]['strand']['+'] > te_hits[te]['strand']['-'] else '-'
            
            MEINFO = '%s,%i,%i,%s' % (te, te_start, te_end, te_strand)


            """ Genotype and genotype quality

            Infer genotype and genotype quality from split reads if present, else
            from discordant read pairs.
            
            Issue: low GQ and QD values with splitreads only

            """
            if split and te in split.te_hits:

                split.get_REF_support(parameters, 'split')

                ref_Q = [split.ref_support[x] for x in split.ref_support]
                alt_Q = [split.te_hits[te]['combined_mapqs'][x] for x in split.te_hits[te]['combined_mapqs']]

                gt_sr = get_genotype(ref_Q, alt_Q)
                
            else:
                gt_sr = []


            if discordant and te in discordant.te_hits:

                discordant.get_REF_support(parameters, 'discordant')

                ref_Q = [discordant.ref_support[x] for x in discordant.ref_support]
                alt_Q = [discordant.te_hits[te]['combined_mapqs'][x] for x in discordant.te_hits[te]['combined_mapqs']]
                
                # Subsample to avoid numpy issues with extremely low or high values
                if len(ref_Q) > 200:
                    ref_Q = sample(ref_Q, 200)
                    
                if len(alt_Q) > 200:
                    alt_Q = sample(alt_Q, 200)
                
                gt_dr = get_genotype(ref_Q, alt_Q)
                
            else:
                gt_dr = []
                
            
            # Add phred qualities if split and discordant support same GT;
            # if conflicting, prefer split GT
            
            if (gt_sr and gt_dr):
                
                if (gt_sr[0] == gt_dr[0]):
                    genotype = [
                        gt_sr[0], 
                        gt_sr[1] + gt_dr[1], 
                        [gt_sr[2][0] + gt_dr[2][0], gt_sr[2][1] + gt_dr[2][1], gt_sr[2][2] + gt_dr[2][2]]
                        ]
                
                else:
                    genotype = [gt_sr[0], gt_sr[1], gt_sr[2]]
                
            elif gt_sr and not gt_dr:
                genotype = gt_sr
                
            elif gt_dr and not gt_sr:
                genotype = gt_dr
            
            
            # Or keep for joined genotyping at later stage?
            if genotype[0] == '0/0':
                 continue

            """ Nr of ALT and REF supporting reads
            
            Reads present both as split and discordant reads are counted once.
      
            """
            supp_REF = set()
            if discordant and te in discordant.te_hits:
                
                supp_REF.update(
                    set(discordant.ref_support.keys())
                    )
                
            if split and te in split.te_hits:
                
                supp_REF.update(
                    {x[:-2] for x in split.ref_support.keys()}
                    )
                
            supp_ALT = 0
            if discordant and te in discordant.te_hits:
                supp_ALT += len(discordant.te_hits[te]['hit_mapqs']) 
            if split and te in split.te_hits:
                supp_ALT += len(split.te_hits[te]['hit_mapqs']) 
                
                  
            AD = '%i,%i' % (len(supp_REF), supp_ALT)
            DP = len(supp_REF) + supp_ALT
            
            
            """ Nr of split and discordant reads
            """
            if discordant and te in discordant.te_hits:
                DR =  len(discordant.te_hits[te]['hit_mapqs']) 
            else:
                DR = 0
                
            if split and te in split.te_hits:
                SR = len(split.te_hits[te]['hit_mapqs'])
            else:
                SR = 0
            
 
            """ Number and proportion of TE positions covered
            """
            AL = len(te_hits[te]['aligned_positions'])
            PROP = AL / len(parameters.telib[te])
            
            
            """ Target site duplication:
            If splitreads overlap, extract the overlapping sequence if it is shorter than 15 bp
            """
            if split:
                position_reverse_sr = min(remove_outliers(split.breakpoint[1]))

                if position_reverse_sr < POS:
                    region = [CHROM, position_reverse_sr + 1, POS]
                    TSD = consensus_from_bam(region, parameters.bamfile, [20, 30])
                    if len(TSD) > 20:
                        TSD = ''
                else:
                    TSD = ''

            else:
                TSD = ''

            TSDLEN = len(TSD)


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
                INFO = 'MEINFO=%s;SVTYPE=INS;TSD=%s;TSDLEN=%i;DPADJ=%i;DR=%i;SR=%i;AL=%i;AP=%f;BPIQR=%i' \
                    % (MEINFO, TSD, TSDLEN, DPADJ, DR, SR, AL, round(PROP, 2), IQR)
            else:
                INFO = 'MEINFO=%s;SVTYPE=INS;HOMSEQ=%s;HOMLEN=%i;DPADJ=%i;DR=%i;SR=%i;AL=%i;IMPRECISE;CIPOS=%i,%i;AP=%f;BPIQR=%i' \
                    % (MEINFO, TSD, TSDLEN, DPADJ, DR, SR, AL, CI_lower, CI_upper, round(PROP, 2), IQR)

            FORMAT = 'GT:GQ:AD:DP:PL'
            
            # Max GT likelihood of 1e6 (to avoid inf values)
            PL = [ str(int(x)) if not math.isinf(x) else '1e6' for x in genotype[2] ]
            GT = '%s:%i:%s:%i:%s' % (genotype[0], genotype[1], AD, DP, ','.join(PL))
            
            outline = [CHROM, POS, '.', REF, ALT, '.', 'PASS', INFO, FORMAT, GT]
            
            vcf_lines.append(outline)
            
            """ Stats
            """
            stats[genotype[0]] += 1
            if te not in stats['TEs']:
                stats['TEs'][te] = 0
                
            stats['TEs'][te] +=1 
            
            
        return vcf_lines, stats


class read_cluster:
    
    """
    Container for read clusters supporting non-reference TE insertion. Here
    mapping results are combined with the anchor reads.
    
    """

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
        """
        Invoked in summmarize_cluster. Here the read_cluster attributes defined
        above are filled, most importantly the te_hits dictionary.
        
        """

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
                    'aligned_positions' : set(),
                    'side' : [0,0]
                    }

            # Counting support from down- and upstream
            self.te_hits[target]['side'][s] += 1

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
        
        """ 
        If reads map to multiple TEs, get the TE with the highest cumulative 
        mapping quality

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
        
        
class TAPs:
    
    """
    Detect TE insertions present in the reference and absent in the
    resequenced accession using read pairs with deviant insert sizes
    
    """
    
    def __init__(self, parameters):
              
        inputs = []
        
        for te in parameters.annotation:
 
            if parameters.region:
                
                if te['chromosome'] == parameters.region[0] and \
                    te['start'] in range(parameters.region[1], parameters.region[2]) and \
                        te['end'] in range(parameters.region[1], parameters.region[2]):
                            
                            inputs.append((te, parameters))
                            
            else:
                inputs.append((te, parameters))
        
        with mp.Pool(parameters.cpus) as pool:
            self.candidates = pool.map(self.process_tap, inputs)
            
    
    def process_tap(self, te_params):
        
        te, parameters = te_params
        
        # Region in which to look for reads
        region_start = te['start'] - 2*parameters.readlength
        region_end = te['end'] + 2*parameters.readlength
        
        # Reset coordinates if they are 'outside' chromosomes
        if region_start < 0:
            region_start = 0
        
        if region_end > len(parameters.ref_contigs[te['chromosome']]):
            region_end = len(parameters.ref_contigs[te['chromosome']]) - 1
            
        region = [te['chromosome'], region_start, region_end]
        
        
        # Get REF and ALT supporting reads and their mapping qualities
        deviant = extract_deviant_reads(te, region, parameters)
        
        # Regional coverage
        cov_mean, cov_stdev = bamstats.local_cov(region, parameters.bamfile)
        
        if not parameters.include_invariant and not deviant[0]:
            return None
            
        else:
            summary = self.TAP_candidate(
                te,
                region,
                deviant, 
                cov_mean
                )
            
            return summary
        
                
    class TAP_candidate:

        def __init__(self, te, region, deviant, cov_mean):

            self.te = te

            self.context = {
                'region_start' : region[1],
                'region_end' : region[2],
                'mean_coverage' : cov_mean
                }
        

            if deviant[0]:

                self.deviant_reads = {
                    
                    'nr_deviant' : len(deviant[0]),
                    
                    # 'deviation' : statistics.median(
                    #     [deviant[0][k][0] for k in deviant[0]]),
                    
                    'deviations' : [deviant[0][k][0] for k in deviant[0]],
                    
                    'start' : max(
                        [deviant[0][k][2] for k in deviant[0]]),
                    
                    'end' : min(
                        [deviant[0][k][3] for k in deviant[0]]),
                    
                    'mapqs' : {k:deviant[0][k][1] for k in deviant[0]}
                    
                    }
                
            else:
                self.deviant_reads = {'mapqs' : {}}
                    
            self.ref_support = deviant[1] if deviant[1] else {}
            
            
    def output_vcf(self, parameters):
        
        """
        Heterozygous TAPs, as they might arise from somatic transposition, 
        are not considered. If both REF and ALT support is found, this is 
        scored as 0/0, that is, presence of the annotated TE.
        
        This is because heterozygous calls are not reliable for short
        DNA transposons, which do not clearly recognizable deviant insert sizes.     
        
        Implement additional tests: 
            - presence of split reads at breakpoints
            - prob of finding n among m reads with deviation x, given isize mean and stdev
        
        
        """
        
        vcf_lines = []
        
        stats = {
            "0/1" : 0,
            "1/1" : 0,
            "0/0" : 0,
            "./." : 0
            }

        for i, site in enumerate(self.candidates):
            
            if not site:
                continue
            
            te = parameters.annotation[i]

            CHROM = te['chromosome']
            POS = te['start']

            REF = parameters.ref_contigs[CHROM][POS - 1]
            ALT = '<DEL:ME>'

            MEINFO = '%s,%i,%i,%s' % (te['id'], 1, te['end']-te['start'], te['strand'])
            SVLEN = -(te['end']-te['start'])
            
            DPADJ = site.context['mean_coverage']
            
            
            # Probability that deviant insert sizes differ from mean isize: one-sample t test
            if site.deviant_reads['mapqs']:
                
                tStat, pVal =  ttest_1samp(
                    site.deviant_reads['deviations'], parameters.isize_mean, axis=0
                    )
                
                PDEV = str(round(pVal, 3))
                
            else:
                PDEV = '1'
            
            INFO = 'MEINFO=%s;SVTYPE=DEL;SVLEN=%i;DPADJ=%i;PDEV=%s' % (MEINFO, SVLEN, DPADJ, PDEV)            
            
            
            """ Calculate genotype likelihood and quality
            
            """
            

            ref_Q = [site.ref_support[k] for k in site.ref_support]
            alt_Q = [site.deviant_reads['mapqs'][k] for k in site.deviant_reads['mapqs']]
            

            # Subsample to avoid numpy issues with extremely low or high values
            if len(ref_Q) > 200:
                ref_Q = sample(ref_Q, 200)
                
            if len(alt_Q) > 200:
                alt_Q = sample(alt_Q, 200)
            
            GT, GQ, PL = get_genotype(ref_Q, alt_Q)
            PL = [str(x) for x in PL]
            
            # Code as missing if GQ is 0
            if GQ == 0:
                GT = './.'
                
            if not parameters.include_invariant and GT in ['0/0', './.']:
                continue
            
            FORMAT = 'GT:GQ:AD:DP:PL'
            
            AD = '%i,%i' % (len(site.ref_support), len(site.deviant_reads['mapqs']))
            
            DP = str(len(site.deviant_reads['mapqs']) + len(site.ref_support))
            
            gt_field = '%s:%i:%s:%s:%s' % (GT, GQ, AD, DP, ','.join(PL))

            outline = [CHROM, POS, '.', REF, ALT, GQ, 'PASS', INFO, FORMAT, gt_field]
            vcf_lines.append(outline)
            
            stats[GT] += 1

        return vcf_lines, stats
            


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
        

        """ Traverse whole bam even if region is set, else discordant 
        mates mapping to TEs elsewhere will be missed!!
        """
        for read in pybam.fetch():

            clipped = is_softclipped(read, parameters.aln_len_SR)

            if (read.is_proper_pair and not clipped) or read.is_secondary:
                continue

            name = read.query_name


            """
            Splitreads. Only reads are considered where exactly one end
            is clipped off.
            """
            if clipped:

                #if AS-XS < parameters.uniq:
                if read.mapping_quality < parameters.mapq:
                    continue

                clseq = read.query_sequence
                if clseq.count('N')/float(len(clseq)) > 0.1:
                    continue

                cigar = read.cigartuples

                chrmsm = pybam.getrname(read.reference_id)
                strt, end = read.reference_start, read.reference_end
                
                if parameters.region:
                            
                    if chrmsm != parameters.region[0]:
                        continue
                    
                    if (strt not in range(parameters.region[1], parameters.region[2]) and
                        end not in range(parameters.region[1], parameters.region[2])):
                        continue
            
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

                # Anchor or not
                if read.mapping_quality < parameters.mapq:
                    uniq = 0
                else:
                    uniq = 1
                    
                if read.is_read1:
                    read_nr = 1
                else:
                    read_nr = 2

                if name in discordant_mates:
                    
                    uniq_mate = discordant_mates[name][2]
                    
                    # Anchor + multimapper: write multimapper to fasta
                    if uniq + uniq_mate == 1:
                        
                        mate_read = discordant_mates[name][0]

                        if uniq == 1:
                            anchor = read
                            multimapper = mate_read
                            
                        else:
                            anchor = mate_read
                            multimapper = read
                            
                        if parameters.region:
                            
                            chrmsm = pybam.getrname(anchor.reference_id)
                            
                            if chrmsm != parameters.region[0]:
                                continue
                            
                            if (anchor.reference_start not in range(parameters.region[1], parameters.region[2]) and
                                anchor.reference_end not in range(parameters.region[1], parameters.region[2])):
                                continue
                            
                        seq = multimapper.query_sequence
                        discordant_anchors[name] = anchor
                            
                        # More than 10% of Ns
                        if seq.count('N')/float(len(seq)) > 0.1:
                            continue

                        header = '>' + name
                        fasta.write(header + '\n' + seq + '\n')

                    del discordant_mates[name]

                else:
                    discordant_mates[name] = [read, read_nr, uniq]

        fasta.close()
        pybam.close()

        return [discordant_anchors, splitreads, split_positions]


def extract_deviant_reads(te,  region, parameters):

        """
        Returns a dictionary d[read name] = [isize, AlignedSegment objects]
        
        IMPROVE: 
            - only look up reads until start position of TE
            - alternative strategy for SE reads: count hard clipped reads as 
            evidence for a TAP?

        """
        
        deviant_read_pairs = {}
        ref_support = {}
    
        # Insert size threshold
        #q_upper = int(norm.ppf(0.95, parameters.isize_mean, parameters.isize_stdev))

        pybam = pysam.AlignmentFile(parameters.bamfile, "rb")

        for read in pybam.fetch(
                te['chromosome'], 
                region[1], 
                region[2],
                multiple_iterators=True):

            name = read.query_name
            
            # Filtering
            # bad read pair
            if (read.mapping_quality < parameters.mapq or 
                not read.is_paired or 
                read.is_secondary or 
                read.is_supplementary):
                continue

            # discard non-properly oriented reads
            if ( (read.mate_is_reverse and read.is_reverse) or 
                (not read.is_reverse and not read.mate_is_reverse) ):
                continue

            isize = abs(read.template_length) - (2*parameters.readlength)            
            
            # isize larger than any TE insertion
            if isize > parameters.longest + 1000:
                continue
            

            # ALT support
            
            # Mate already identified as deviant
            if name in deviant_read_pairs: 
                
                deviant_read_pairs[name].append(read.reference_start)
                deviant_read_pairs[name][1] += read.mapping_quality
            
            # New deviant pair
            elif ( read.reference_end <= te['start'] and read.next_reference_start >= te['end'] ):

                deviant_read_pairs[name] = [isize, read.mapping_quality, read.reference_end]
                
                
            # REF support
            else:
                
                if name in ref_support:
                    ref_support[name] += read.mapping_quality
                    
                else:
                    
                    # Allow few base pairs as error margin. If not, a single position
                    # overlap will be considered as supporting the reference allele ...
                    insert_interval = range(read.reference_end-2, read.next_reference_start+2)

                    if (  (te['start'] in insert_interval and not te['end'] in insert_interval) or
                          (te['end'] in insert_interval and not te['start'] in insert_interval) ):
                        
                        ref_support[read.query_name] = read.mapping_quality
                        
        pybam.close()
        
        # Remove if mate of deviant read pair has not been recovered
        remove = []

        for name in deviant_read_pairs:
            if len(deviant_read_pairs[name]) != 4:
                remove.append(name)

        for name in remove:
            deviant_read_pairs.pop(name)
            
        return deviant_read_pairs, ref_support


def remove_outliers(lista):
    # outliers defined as in R's boxplots
    q_1 = scipy.percentile(lista, 25)
    q_3 = scipy.percentile(lista, 75)

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


def minimap(queries, parameters, min_aln_len, k, w):

    """ Map discordant reads and split reads against TE consensus library. Mapq is the crucial
    output here, as it is later used to calculate genotype likelihoods.
    
    k: Minimizer k-mer length
    w: Minimizer window size [2/3 of k-mer length]. A minimizer is the smallest k-mer in a window of w consecutive k-mers.
    
    
    If TEs contain internal repeats, e.g LTR retrotransposons, this can result in multimappers and
    a mapq of 0. To correct for this, mapq is recalculated, setting s2 (f2 in Li 2018) to 0.

    mapq = 40 * (1-s2/s1) * min(1,m/10) * log(s1)

    Li 2018: "Minimap2: Pairwise alignment for nucleotide sequences"
    """

    cmd = [
        'minimap2', 
        '-xsr', 
        '--secondary=yes', 
        '-k', str(k), 
        '-w', str(w),
        parameters.targets, 
        queries
        ]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    output_raw = proc.stdout.read().splitlines()

    outfile = open('%s.k%i.paf' % (parameters.outname, k), 'w')

    minimapd = {}

    for line in output_raw:

        line = line.decode('utf-8')
        outfile.write(line + '\n')
          
        fields = line.strip().split()
        
        if int(fields[10]) < min_aln_len:
            continue
        
        readname = fields[0]

        d = {

            'strand' : fields[4],

            'target_name' : fields[5],
            'target_start' : int(fields[7]),
            'target_end' : int(fields[8]),
            'target_range' : set(range(int(fields[7]), int(fields[8]))),

            'aln_block_length' : int(fields[10]), # total nr of matches, mismatches, indels
            'mapping_quality' : int(fields[11]),

            'tp' : re.search(r'tp:A:(.*?)\t', line).group(1),
            'cm' : int(re.search(r'cm:i:(.*?)\t', line).group(1)),
            's1' : int(re.search(r's1:i:(.*?)\t', line).group(1)),
            's2' : int(re.search(r's1:i:(.*?)\t', line).group(1))

            }

        # Secondary alignments: if they align to the same element, recalibrate mapping quality.
        if readname in minimapd:

            if d['target_name'] == minimapd[readname]['target_name']:

                d_primary = minimapd[readname]

                # the min(1, m/10) term in the formula,
                # xmin = d_primary['mapping_quality'] / \
                #     (40 * (1 - (d_primary['s2'] / d_primary['s1'])) * np.log(d_primary['s1']))
                
                xmin = 1 if (minimapd[readname]['cm'] / 10) >= 1 else minimapd[readname]['cm']
                
                mapq_recal = 40 * xmin * np.log(d_primary['s1'])
                
                if mapq_recal > 60:
                    mapq_recal = 60
                
                minimapd[readname]['mapping_quality'] = mapq_recal

                # if s1 and s2 are equal, add positions form secondary alignment
                if d['s1'] == d_primary['s1']:

                    minimapd[readname]['target_range'].update(d['target_range'])

        else:
            minimapd[readname] = d

    return minimapd


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

    te_hits = copy.deepcopy(DR_TE_hits)

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
    """ If reads map to multiple TEs, get the TE with the highest cumulative mapping quality

    PROBLEM: mapq of 0 for reads mapping into repeats within elements, e.g. LTRs!!

    """
    score = -1
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

    return GT, GQ, Q_norm


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


