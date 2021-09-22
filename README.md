[![GitHub Downloads](https://img.shields.io/github/downloads/cstritt/detettore/total.svg?style=social&logo=github&label=Download)](https://github.com/cstritt/detettore/releases)


*detettore* – a program to detect transposable element polymorphisms
====================================================================
June 2021

<img src="detettore_ad.png" alt="drawing" width="600"/>  

\
*detettore* uses reference-aligned paired-end reads to search for:

  - **TE insertion polymorphisms (TIPs)**, i.e. TEs absent in the
    reference genome but present in a sequenced individual
  - **TE absence polymorphisms (TAPs)**, i.e. TEs present in the
    reference but absent in a sequenced individual


## Table of Contents
1. [Installation](#install)
2. [Usage](#usage)
  1. [Single sample](#single)
  2. [Multiple samples](#multiple)
  3. [Combine VCF files](#combineVCFs)
3. [TE toolbox](#tools)
4. [New in version 2](#version)
5. [License](#license)


## <a name="install"></a>Installation
*detettore* is written in Python 3 and available on PyPI. To avoid conflicts
with dependencies, it is best to install it in a virtual environment:

``` bash
# Create environment called detettore
conda create -n detettore python=3.7

# Activate environment
source activate detettore

# Install detettore
pip install detettore

```
Or when not using Anaconda:
``` bash
# Create environment, where <location> is the path to the environment
virtualenv -p python3.7 <location>

# Activate environment
source <location>/bin/activate

# Install detettore
pip install detettore

```

## <a name="usage"></a>Usage


### <a name="single"></a>Single sample
Basic usage illustrated with the data in the example folder.

``` bash
detettore.py \
  -b example/reads.bam \
  -r example/reference.fasta \
  -a example/TE_annotation.gff \
  -t example/TE_consensusLib.fasta \
  -o example \
  -m tips taps \
  --require_split \
  --include_invariant

```

#### Explanation of command line parameters

| Parameter               | Explanation
|-                        |-
|**Input/Output**         |
|\-b                      | bam file with reference-aligned paired-end reads.
|\-r                      | Reference genome in fasta format.
|\-t                      | TE consensus sequences in fasta format.
|\-a                      | TE annotation in bed or gff format.
|\-o                      | Sample name, used as a prefix for output files.
|**Program settings**     |
|\-m                      | The module to run (tips, taps, or both, as above).
|\-c                      | Number of CPUs.
|\--region                | Restrict search to region chromosome:start:end.
|\--include_invariant     | Include conserved TEs in vcf output.
|\--require_split         | Discard variant candidates if no splitread evidence is present.
|\--keep                  | Keep intermediate files.
|**Thresholds**           |         
|\-q                      | Minimum mapping quality of reference-aligned anchor reads.
|\-lDR                    | Minimum alignment length for discordant read-pair mtarget hits. [50]
|\-lSR                    | Minimum alignment length for splitread target hits. [20]


Points to consider:
  - **Chromosome names** must be consistent in the different files, which can be a problem when files are downloaded from different sources.
  - The bam file should be as **unfiltered** as possible: files containing only properly paired or uniquely mapping reads, while useful for SNP calling, lack the information tapped by *detettore*.


#### Output
The command above will produce a file called example.vcf.gz containing TIPs and TAPs,
and a log file example.log. If --keep was set, a folder sample_tmp will also be present,
containing ...


### <a name="multiple"></a>Multiple samples
Here is an example workflow to call TE polymorphisms on multiple samples.
GNU parallel is used to run 10 samples in parallel, with one CPU per sample.
A different approach might be required on computing clusters with job queing.


``` bash
# Set input files which do not change for different samples
ref=/path/to/reference.fasta
annot=/path/to/TE_annotation.gff
telib=/path/to/TE_consensus.fasta

# Create a list of detettore commands, looping through a list containg paths to bam files.
# The files are assumed to be named sample.bam
while read bampath;
do

  sample=$(basename $bampath | cut -d'.' -f1)

  echo "\
  detettore.py \
    -b $bampath \
    -r $ref \
    -a $annot \
    -t $telib \
    -o $sample \
    -m tips taps \
    -c 1 \
    --require_split \
    --include_invariant"\
  >> run_detettore.cmds

done < paths_to_bam.txt

# Use GNU parallel to run 10 commands simultaneously. Save stderr and stdout to log files.
parallel -j10 < run_detettore.cmds 2> err.log > stdout.log

```

### <a name="combineVCFs"></a>Combine VCFs



## <a name="version"></a>**New in version 2**
- all output in vcf format
- output invariant sites
- genotype calling based on genotype likelihoods
- compatibility with minimap2
- tools for downstream analysis

## <a name="licence"></a>Licence
GNU General Public License v3. See LICENSE.txt for details.
