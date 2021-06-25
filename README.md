*detettore* – a program to detect transposable element polymorphisms
================
June 2021

<img src="detettore_ad.png" alt="drawing" width="500"/>


*detettore* is a program to detect and characterize transposable element
(TE) polymorphisms using reference-aligned paired-end reads.

The program searches for:

  - **TE insertion polymorphisms (TIPs)**, i.e. TEs absent in the
    reference genome but present in a sequenced individual
  - **TE absence polymorphisms (TAPs)**, i.e. TEs present in the
    reference but absent in a sequenced individual


**New in version 2**:

  - all output in vcf format
  - output invariant sites
  - genotype calling based on genotype likelihoods
  - compatibility with minimap2

  - tools for downstream analysis
  - compatibility with minimap2


## Contents
1. [Installation](#installation)
2. [Usage](#usage)
  1. [Single sample](#single)
  2. [Multiple samples](#multiple)
3. [License](#license)


## Installation <a name="installation"></a>

Download the *detettore* directory from Github.

``` bash
git clone https://github.com/cstritt/detettore
```


## Usage <a name="usage"></a>
### Single sample <a name="single"></a>
Basic usage illustrated with the data in the example folder.

``` bash

detettore.py \
  -b  ~/github/detettore/example/SRR10828834.bam \
  -t /home/cristobal/TB/data/ISfinder/IS.fna \
  -r ~/github/detettore/example/reference.fasta \
  -a ~/github/detettore/example/TE_annotation.gff \
  -q 20 \
  -c 4 \
  -o example \
  -m tips taps \
  --include_invariant

```

### Multiple samples <a name="multiple"></a>


``` bash


```



## Licence <a name="licence"></a>
GNU General Public License v3. See LICENSE.txt for details.
