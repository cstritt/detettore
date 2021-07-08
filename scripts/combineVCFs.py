#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""" Combine single sample vcf files

Adds new INFO fields:
    
    #QD: quality by depth
    AC: allele count
    AF: allele frequency
    NHET: number of heterozygous genotypes
    
    
class args:
    def __init__(self):
        
        self.vcf_path = '/home/cristobal/TB/analyses/detettore/tb_1227/risultati/vcf'
        self.filter =  ['depth']
        self.tbl = False
        
class args:
    def __init__(self):
        
        self.vcf_path = '/home/cristobal/Bdis/analyses/detettore2/detetorre_Hannes_All_july2021'
        self.filter =  ['depth']
        self.tbl = False
        
args = args()

   
Created on Thu Jun 17 11:46:10 2021
@author: cristobal

"""

import argparse
import gzip
import os
import re
import sys
import time


def get_args():

    parser = argparse.ArgumentParser(
        description='Detect TE polymorphisms from paired-end read data')

    parser.add_argument(
        "-p",
        dest="vcf_path",
        required=True,
        help='Path to folder containing vcf.gz output of detettore.')
    
    parser.add_argument(
        "-f",
        dest="filter",
        nargs='+', 
        choices=['depth', 'min_supp2', 'require_split', 'nohet'],
        default=['depth', 'min_supp2'],
        help='Apply filters:\n\
            depth : discard variants in regions with a P(depth) < 0.05 \n\
            min_supp2: \n\
            require_split: \n\
            nohet : discard 0/1 genotypes')

    parser.add_argument(
        '--keep_invariant',
        action='store_true',
        help='Keep invariant sites.')

    args=parser.parse_args()

    return args


def parse_vcf_lines(handler, sample_order, depth_d, vard, imprecise):
    
    for line in handler:
        
        try:
            line = line.decode('utf-8')
        except AttributeError:
            pass
        
        if line.startswith("##"):
            continue

        fields = line.strip().split('\t')
        
        if line.startswith('#CHROM'):
            sample = fields[-1]
            sample_order.append(sample)
            depth_d[sample] = []
            continue
        
        # DPADJ
        try:
            dpadj = re.search(r'DPADJ=(.*?);', fields[7]).group(1)
        except AttributeError:
            dpadj = re.search(r'DPADJ=(.*?)$', fields[7]).group(1)
        depth_d[sample].append(int(dpadj))
         
        
        chromosome, pos = fields[0], int(fields[1])
        var_type = fields[4]
        
        meinfo = re.search(r'MEINFO=(.*?);', fields[7]).group(1)
        te = meinfo.split(',')[0]  
        
        imprec = True if 'IMPRECISE' in fields[7] else False
        
        
        uniq_id = (chromosome, pos, te)
        
        # TAPs and precise TIPs
        if (var_type == '<DEL:ME>') or (var_type == '<INS:ME>' and not imprec):
            
            if uniq_id not in vard[var_type]:
                
                vard[var_type][uniq_id] = {}
            
            vard[var_type][uniq_id][sample] = fields
            
        # Imprecise into separate list
        elif imprec:
            
            imprecise.append((sample, fields))
    
    
#%% MAIN
def main():

    args = get_args()
    
    
    #%% Load variants
    
    vard = { '<INS:ME>' : {}, '<DEL:ME>' : {} }
    
    depth_d = {}
    
    sample_order = []
    imprecise = []
    
    for subdir, dirs, files in os.walk(args.vcf_path):
        
        for f in files:
                   
            if f.endswith('.gz'):
                handler = gzip.open(subdir + '/' + f)
                
            else:
                handler = open(subdir + '/' + f)
                
            parse_vcf_lines(handler, sample_order, depth_d,vard, imprecise)
            
            handler.close()
                
            

    #%% Add imprecise variants
    
    """
    Join if the confidence interval of their position (CIPOS) overlaps with 
    precise or other imprecise variants    
    """        
                  
    for sample, fields in imprecise:
            
        chromosome, pos = fields[0], int(fields[1])                    
        meinfo = re.search(r'MEINFO=(.*?);', fields[7]).group(1)
        te = meinfo.split(',')[0]
        
        cipos = re.search(r'CIPOS=(.*?)$', fields[7]).group(1).split(',')
    
        closest = [int(), float('inf')]
           
        for i in range(int(cipos[0]), int(cipos[1])):
            
            if (chromosome, i, te) in vard['<INS:ME>']:
                
                if abs(pos-i) < closest[1]:
                    closest[0] = i
                    closest[1] = abs(pos-i)
                    
        if closest[0]:
            vard['<INS:ME>'][(chromosome, closest[0], te)][sample] = fields
     
        else:
            vard['<INS:ME>'][(chromosome, pos, te)] = {sample : fields}
                

    #%% Write combined vcf
    """
    New INFO fields:
        
        QD: quality by depth
        AC: allele count
        AF: allele frequency
    
    """
    
    metainfo = [
        '##fileFormat=VCFv4.2',
        '##fileDate=%s' % time.strftime("%d/%m/%Y"),
        '##source==detettore v2.0 combineVCF',

        '##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">',
        '##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count">',
        '##INFO=<ID=AF,Number=1,Type=Integer,Description="Allele frequency">',
        '##INFO=<ID=NHET,Number=1,Type=Integer,Description="Number of heterozygous genotypes">',
     
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">',

        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % ('\t'.join(sample_order))
        
        ]
    
    
    combined = []
    
    for var_type in ['<DEL:ME>', '<INS:ME>']:
        
        for site in vard[var_type]:
            
            variants = vard[var_type][site]
            samples = list(variants.keys())
    
            # Fields 1-5 do not change
            outline = variants[samples[0]][:5]
            
            # allele count and frequency
            genotypes = []
            
            for sample in sample_order:
                if sample in variants:
                    #GT = variants[sample][9]
                    GT = variants[sample][9].split(':')[0]
                    
                else:
                    #GT = '0/0:.:.,.:.'
                    GT = '0/0'
                    
                genotypes.append(GT)
            
            genotypes_str = [x.split(':')[0] for x in genotypes]
            NHET = genotypes_str.count('0/1')
        
            ref = ''.join(genotypes_str).count('0')
            AC = ''.join(genotypes_str).count('1')
            
            if AC + ref == 0:
                continue

            AF = AC / (AC + ref)
            
            INFO = 'MEINFO=%s;AC=%i;AF=%f;NHET=%s' % \
                (site[2], AC, AF, NHET)
            
        
            #outline = outline + [str(QUAL), 'PASS', INFO, 'GT:GQ:AD:DP'] + genotypes
            outline = outline + ['.', 'PASS', INFO, 'GT'] + genotypes
            
            combined.append(outline)
        
    combined = sorted(combined, key=lambda x: (x[0], int(x[1])))  
    

    #%% Write to stdout
    
    for line in metainfo:
        sys.stdout.write(line + '\n')
        
    for line in combined:
        sys.stdout.write('\t'.join(line) + '\n')
        
        
if __name__ == '__main__':
    main()
