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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from joblib import Parallel, delayed
from scipy.stats import norm
from scripts import strumenti


class annotation_entry:
    def __init__(self, line, file_format):

        fields = line.strip().split('\t')

        if file_format == 'bed':
            self.chromosome = fields[0]
            self.start = int(fields[1])
            self.end = int(fields[2])
            self.id = fields[3]
            self.strand = fields[5]
            
        elif file_format == 'gff':
            self.chromosome = fields[0]
            self.start = int(fields[3])
            self.end = int(fields[4])
            self.id = re.search(r'Name=(.*?)[;\n]', line).group(1)
            self.feature = fields[2]
            self.strand = fields[6]

        else:
            print 'Check annotation format.'


class feature_summary:
    def __init__(self, fields):

        self.indx = fields[0]
        self.length = fields[1]
        
        # coverage in the surrounding region
        self.region_start = fields[2]
        self.region_end = fields[3]
        self.cov_mean = fields[4]
        self.cov_stdev = fields[5]

        # slots to be filled if there are deviant read pairs 
        self.nr_deviant = int()
        self.deviation = int()
        self.start = int()
        self.end = int()
        
        # slots to be filled if the feature should be reconstructed
        self.known = float()
        self.seq = str()
        
        
    def add_deviant_info(self, deviant_summary):
        
        self.nr_deviant = len(deviant_summary['name'])
        self.isize = statistics.median(deviant_summary['isizes']) 
        try:
            self.start = max(deviant_summary['f_strand_ends'])
            self.end = min(deviant_summary['r_strand_starts'])
        except ValueError:
            pass
        
        
    def reconstruct(self, region, bamfile, filters):
        
        seq = strumenti.consensus_from_bam(region, bamfile, filters)
        counts = Counter(seq)
        unknown = counts['N']
        known_proportion = 1 - (unknown/float(len(seq)))
        known = known_proportion
        return seq, known


    def write(self, annotation):
        
        i = self.indx

        outline = [annotation[i].chromosome,
                   annotation[i].start,
                   annotation[i].end,
                   annotation[i].id,
                   annotation[i].strand,
                   self.length,
                   self.isize,
                   self.nr_deviant,
                   self.start,
                   self.end,
                   self.cov_mean,
                   self.cov_stdev]
        
        return map(str, outline)
        

def extract_deviant_reads(region, isize_stats, bamfile, mapq):

    pybam = pysam.AlignmentFile(bamfile, "rb")
    
    isize_mean, isize_stdev, readlength = isize_stats
    q_upper = int(norm.ppf(0.99, isize_mean,isize_stdev))

    chromosome, start, end = region

    out = {'name':[], 'f_strand_ends':[], 'r_strand_starts':[], 'isizes':[]}
    
    skip = set()

    for read in pybam.fetch(chromosome, start, end):
        
        name = read.query_name
        # read pair already considered
        if name in skip:
            continue
        # bad read
        if read.is_secondary or not read.is_paired or read.mapq < mapq:
            skip.add(name)
            continue
        # insert size not deviant
        isize = abs(read.template_length) - (2*readlength)
        if not (isize > q_upper and isize < 30000):
            skip.add(name)
            continue
        
        """ Get 'neat deviant reads with 5' on forward and 3' on reverse
        """
        if read.mate_is_reverse and not read.is_reverse:
       
            out['name'].append(name)
            out['isizes'].append(isize)   
            out['f_strand_ends'].append(read.reference_end)
            out['r_strand_starts'].append(read.next_reference_start)                
            skip.add(name)
            
    pybam.close()       
    return out


def process_taps(annot_entry, isize_stats, bamfile, 
                 contigs, uniq, reconstruct):
    
    """ Check if there are deviant read-pairs around the annotated feature
    """
        
    i, feature = annot_entry
    isize_mean, isize_stdev, readlength = isize_stats

    region_start = feature.start - (isize_mean + isize_stdev) 
    region_end = feature.end + (isize_mean + isize_stdev)
    
    # Reset coordinates if they are 'outside' chromosomes
    if region_start < 0:
        region_start = 0
    if region_end > contigs[feature.chromosome]:
        region_end = contigs[feature.chromosome] - 1

    region = (feature.chromosome, region_start, region_end)

    feature_length = feature.end - feature.start
    deviant = extract_deviant_reads(region, isize_stats, bamfile, uniq)    
    
    if not deviant['name'] and not reconstruct:
        return None
    
    cov_mean, cov_stdev = strumenti.local_cov(region, bamfile)
    summary = feature_summary([i, 
                               feature_length, 
                               region_start, 
                               region_end,
                               cov_mean, 
                               cov_stdev])
   
    if deviant['name']:
        summary.add_deviant_info(deviant)
        
    return summary


def run_module(bamfile, 
               readinfo, 
               annotation_file, 
               reference,
               thresholds, 
               reconstruct,
               cpus):

    """ Read in annotation and read information
    """
    annot = []
    file_format = annotation_file.split('.')[-1]
    with open(annotation_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            annot_row = annotation_entry(line, file_format)
            annot.append(annot_row)

    with open(readinfo) as f:
        readlength = int(next(f).split('\t')[1])
        next(f)
        next(f)
        isize_mean = int(next(f).split('\t')[1])
        isize_stdev = int(next(f).split('\t')[1])

    isize_stats = (isize_mean, isize_stdev, readlength)
    uniq = thresholds[0]
    
    contigs = {seq_record.id: len(seq_record) for seq_record in SeqIO.parse(reference, "fasta")}

    baseq, mapq = 20, 30
    filters = [baseq, mapq]
    recs = []


    print('Extracting read pairs with deviant insert sizes around annotated features\n...')
    inputs = [(i, x) for i, x  in enumerate(annot)]
    out = Parallel(n_jobs=cpus)(delayed(process_taps)\
                   (i,
                    isize_stats,
                    bamfile,
                    contigs,
                    uniq,
                    reconstruct) for i in inputs)   
    
    header = ['chromosome',
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
              'cov_stdev']

    
    outfile = open('TAPs_candidates.tsv', 'w')
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
                start_overlap = strumenti.overlapping_reads(chromosome,
                                                            element_start,
                                                            bamfile, 20)
            
                end_overlap = strumenti.overlapping_reads(chromosome,
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
            five_bridge = strumenti.overlapping_reads(annot_row.chromosome, 
                                                      annot_row.start,
                                                      bamfile, 20)
            
            three_bridge = strumenti.overlapping_reads(annot_row.chromosome, 
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
                        
    