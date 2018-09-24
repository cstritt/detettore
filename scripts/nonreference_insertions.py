#!/usr/bin/env python

""" dettetore -a program to detect and characterize transposable element polymorphisms

Detect TE insertion polymorphisms (TIPs)

Copyright (C) 2018 C. Stritt
License: GNU General Public License v3. See LICENSE.txt for details.
"""

import pysam
from scripts import strumenti


def get_split_and_discordant_reads(bamfile, 
                                   min_aln_len,
                                   readlength,                                   
                                   uniq_threshold):
    
    """ Extract discordant read pairs in which at least one read maps uniquely 
    to the reference, i.e. AS - XS > uniqueness_threshold. Write uniquely 
    mapping reads to anchors dicitonary, their mate to fasta file.
    
    """
    pybam = pysam.AlignmentFile(bamfile, "rb")

    discordant_anchors = dict()
    discordant_mates = dict()   
    fasta = open('discordant.fasta', 'w')
    
    splitreads = dict()
    split_positions = dict()

    for read in pybam.fetch(): 
        
        clipped = is_softclipped(read, min_aln_len)
        
        if ((read.is_proper_pair and not clipped) or 
            read.is_secondary):
            continue
        
        name = read.query_name
        
        
        """ While BWA always outputs the XS tag, Bowtie2 only does it when
        there is a secondary alignment
        """
        AS = read.get_tag('AS')
        try:  
            XS = read.get_tag('XS')
        except KeyError:
            
            XS = 0


        """ Splitreads. Only reads are considered where exactly one end
        is clipped off.
        """
        if clipped:
            
            if AS-XS < uniq_threshold:
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
                interval = (end-readlength, end)

            elif cigar[-1][0] == 4:
                interval = (strt, strt+readlength)
                
            if chrmsm not in split_positions:
                split_positions[chrmsm] = []

            split_positions[chrmsm].append([name, interval[0], interval[1]])
 
        """ Discordant read pairs
        """
        if not read.is_proper_pair:
            
            if AS-XS < uniq_threshold:
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

    
def is_softclipped(read, min_clip_length):
    """ Return FALSE if read is not clipped.
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


def get_anchor_positions(hits, anchors, bamfile, isize):
    
    """ Go through the blast hits and get the corresponding anchor positions
    """

    pybam = pysam.AlignmentFile(bamfile, "rb")
    positions = dict() 
    
    for read in hits.keys():
        
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
        positions[k].sort(key=lambda x: int(x[1])) 
    
    pybam.close()
    return positions


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
                seq = strumenti.clip_seq(readseq, cigar)
                
                if len(seq) < min_clip_length:
                    continue
                outline = header + '\n' + seq + '\n'
                f.write(outline)
    f.close()
    return fasta


def run_module(bamfile,
               targets,
               thresholds,
               readinfo,
               cpus):


    with open(readinfo) as f:
            line = next(f)
            readlength = int(line.split("\t")[1])    
            
            
    with open(readinfo) as f:
        lines = f.readlines()
        readlength = int(lines[0].split("\t")[1].strip())
        isize_mean = int(lines[3].split("\t")[1].strip())        
    
            
    uniq, min_aln_len, min_perc_id, word_size = thresholds
               
    
    print('Getting candidate split reads and discordant read pairs from bam file ...')
    reads = get_split_and_discordant_reads(bamfile, 
                                           min_aln_len, 
                                           readlength, 
                                           uniq)
    
    discordant_anchors, splitreads, split_positions = reads
    
    print('Analysis begins with ' + str(len(discordant_anchors)) + 
          ' discordant pairs and ' + str(len(splitreads)) + ' splitreads') 
    
    
    # DISCORDANT READS
    print('Aligning discordant reads to target sequences\n...')
    strumenti.create_blastdb(targets)
    blast_out = strumenti.blastn('discordant.fasta', 
                                 min_perc_id,
                                 11,
                                 cpus)
    
    discordant_hits = strumenti.hit_dictionary(blast_out, min_aln_len, 'one_hit_only')
    print(' '.join([str(len(discordant_hits.keys())), 'anchor mates map to a TE.']))
    
    
    
    print('Clustering discordant reads\n...')   
    anchor_positions = get_anchor_positions(discordant_hits, discordant_anchors, 
                                            bamfile, isize_mean)  
    
    discordant_clusters = strumenti.cluster_reads(anchor_positions, 0.05)
    
    discordant_clusters = strumenti.prefilter_clusters(discordant_clusters, 4)
    
    for k in sorted(discordant_clusters.keys()):
        print(k + ': ' + str(len(discordant_clusters[k])) + ' clusters')
        
        
    # summarize discordant clusters
    summaries = strumenti.combine_hits_and_anchors(
            discordant_clusters,
            discordant_anchors,
            discordant_hits)
    
    cluster_summaries, discordant_cluster_positions = summaries

    
    # SPLITREADS
    print('Clustering split reads\n...')
    split_clusters = strumenti.cluster_reads(split_positions, 0.05)
    split_clusters = strumenti.prefilter_clusters(split_clusters, 4)
    
    for k in sorted(split_clusters.keys()):
        print(k + ': ' + str(len(split_clusters[k])) + ' clusters')
        
        
    print('Aligning split parts to target sequences\n...')
    fasta = write_clipped_to_fasta(split_clusters, splitreads, min_aln_len)

    blast_out_split = strumenti.blastn(fasta, min_perc_id, word_size, cpus)
    
    split_hits = strumenti.hit_dictionary(blast_out_split, min_aln_len, 'one_hit_only')    
    print(' '.join([str(len(split_hits.keys())), 'soft-clipped reads map to a TE.']))
    
    
    # COMBINE DISCORDANT READS AND SPLITREADS
    combined = strumenti.add_splitread_evidence(cluster_summaries,
                                                         split_clusters,
                                                         splitreads,
                                                         split_hits,
                                                         discordant_cluster_positions)
    
    strumenti.create_table(combined, bamfile)
    
