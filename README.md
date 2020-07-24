*detettore* – a program to detect transposable element polymorphisms
================
April 2020

![](detettore_ad.png)

*Detettore* is a program to detect and characterize transposable element
(TE) polymorphisms using reference-aligned paired-end reads.

The program searches for:

  - **TE insertion polymorphisms (TIPs)**, i.e. TEs absent in the
    reference genome but present in a sequenced individual
  - **TE absence polymorphisms (TAPs)**, i.e. TEs present in the
    reference but absent in a sequenced individual

## Requirements

BLAST+ is the only external dependency
(<https://blast.ncbi.nlm.nih.gov>).

**Minimal data requirements**

  - Illumina paired-end reads aligned to a reference genome with BWA-MEM
    (Li 2013)  
  - TE consensus sequences  
  - reference genome  
  - TE annotation in bed or gff (only required for TAPs)

**Note on the paired-end reads**  
It is assumed that the paired-end reads were aligned to a reference
genome with BWA-MEM (Li 2013). The reason for this restriction is that
different aligners use different ways to indicate whether a read maps
uniquely to the reference or at multiple places.

## Download and install *detettore*

Download the *detettore* directory from Github.

``` bash
git clone https://github.com/cstritt/detettore
```

It is most convenient to install *detettore* in a Python virtual
environment. Here are two ways to set up a virtual environment.

``` bash
# When using the Anaconda distribution
conda create -n detettore python=3.7

# Without Anaconda, replacing <location> with the path to the environment
virtualenv -p python3.7 <location>
```

Finally, activate the environment and install *detettore*.

``` bash
# With Anaconda
conda activate detettore

# Without
source <location>/bin/activate    

# Install detettore
cd ~/programs/detettore
python setup.py install
```

## Usage

#### Command line parameters of detettore.py

| Parameter     | Explanation                                                                                                                                                                 |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| \-o           | Output folder name, e.g. the sample ID                                                                                                                                      |
| \-b           | Unfiltered (\!) bam file, i.e. including the not uniquely-mapping reads                                                                                                     |
| \-m           | The module to run. Can be both at the same time or just one \[tips, taps\]                                                                                                  |
| \-t           | TE consensus sequences in fasta format                                                                                                                                      |
| \-r           | Reference genome in fasta format                                                                                                                                            |
| \-a           | TE annotation in bed or gff format                                                                                                                                          |
| \-u           | Minimum difference between primary and secondary alignment score. Reads above the threshold are considered as mapping uniquely to the reference \[30\]                      |
| \-lSR           | Minimum length of soft-clipped read parts \[15\]                                                                                                                            |
| \-lDR           | Minimum alignment length of discordant read-pair target hits \[50\]      
| \-id          | Minimum sequence identity between reads and TE consensus sequences \[80\]                                                                                                   |
| \-ws          | Word size (minimum length of best perfect match) for blasting splitreads against TEs \[11\]                                                                                 |
| \-c           | Number of CPUs. The blast search can be run with multiple CPUs, as well as a time-consuming loop in the TAP module \[1\]                                                    |
| \-keep        | Keep fasta files with discordant read and clipped sequences, as well as the output of the blast search of theses seqences against the TE consensus sequences                |
| \-reconstruct | If the TAP module is used, TE sequences present in the reference and the sequenced individual can be reconstructed and written to a fasta file. Slow for large annotations. |

**Naming conventions**  
If a gff file is used in the TAP module, the TE name in the last column
of the gff is assumed to be formatted as “Name=TE\_name”.

## Structure of output files

*Detettore* provides a wealth of information for TIP candidates. Besides
the number of supporting discordant and split reads, the number of
bridge reads are reported, that is, the number of reads spanning the
insertion breakpoint and thus indicating heterozygosity or a false
positive. Furthermore, for each candidate the regional coverage is
shown, which allows to filter out repetitive and low-confidence
regions.

#### tips.tsv

| Header                 | Explanation                                                                                                                       |
| ---------------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| chromosome             |                                                                                                                                   |
| position               |                                                                                                                                   |
| TE                     |                                                                                                                                   |
| strand                 |                                                                                                                                   |
| nr\_supporting\_reads  | Total number of supporting reads. Evidence from the 5’ and the 3’ site is separated with a dash.                                  |
| nr\_discordant\_reads  |                                                                                                                                   |
| nr\_splitreads         |                                                                                                                                   |
| nr\_aligned\_positions | Nr of positions in the TE consensus sequence that are covered by the discordant reads and split reads reaching into the insertion |
| TSD                    | Target site duplication sequence                                                                                                  |
| nr\_bridge\_reads      | Number of reads which bridge the insertion breakpoint and point to heterozygous insertions or false positives                     |
| region\_start          | Coordinates of the region surrounding the TE insertion                                                                            |
| region\_end            |                                                                                                                                   |
| region\_cov\_mean      | Mean coverage in this region. Useful indicator of messiness and a good filtering criterion.                                       |
| region\_cov\_stdev     |                                                                                                                                   |

#### taps.tsv

| Header             | Explanation                                              |
| ------------------ | -------------------------------------------------------- |
| chromosome         |                                                          |
| TE\_start          | Start of the reference TE insertion                      |
| TE\_end            | End of the reference TE insertion                        |
| TE\_id             | ID of the reference TE                                   |
| TE\_length         | Length of the reference TE                               |
| isize              | Mean insert size of the deviant read-pairs around the TE |
| nr\_discordant     | Nr of deviant read-pairs                                 |
| absence\_start     | Beginning of the gap                                     |
| absence\_end       | End of the gap                                           |
| region\_start      | Same as for TIPs                                         |
| region\_end        |                                                          |
| region\_cov\_mean  |                                                          |
| region\_cov\_stdev |                                                          |


## Making sense of the output

Here two possible applications of *detettore* are shown.

1)  In the first, we have a have a candidate gene and want to know
    whether it contains TE insertion polymorphisms.  
2)  In the second, individual output files are pooled into a vcf file to
    study the variation a particular family generated in a natural
    population.

#### 1\. Search TIPs in and around a candidate region

``` r
# Define some candidate genes
candidate_genes <- list(
  'Bradi1g77130' = c('Bd1', 73734763,   73736724),
  'Bradi4g09872' = c('Bd4', 9394780,    9398745),
  'Bradi5g01366' = c('Bd5', 1292944,    1297431)
  )

# Load TIPs
tip_files <- list.files('~/github/detettore/example/tips', full.names = T, recursive = T, pattern="tips", include.dirs = T)

tips <- data.frame()

for (f in tip_files){

  t <- read.table(f, header=T)

  # Add accession name to table
  acc <- strsplit(f, split = "/")[[1]][8]
  t$acc <- as.factor(rep(acc, nrow(t)))

  tips <- rbind(tips, t)
}

# TIPs in and around (+/- 1000 bp) candiate genes
candidate_tips <- list()

for (gene in names(candidate_genes)){

  chr = candidate_genes[[gene]][1]

  # Add 1000 base pairs to gene coordinates to detect close TIPs
  start = as.numeric(candidate_genes[[gene]][2]) - 1000
  end = as.numeric(candidate_genes[[gene]][3]) + 1000
  candidate_tips[[gene]] = subset(tips, chromosome == chr & position %in% c(start:end))
}

# Show results
candidate_tips
candidate_tips$Bradi1g77130
```

#### 2\. Combining variants into a vcf file

*The variantCaller.py program is currently being updated.*

# License

GNU General Public License v3. See LICENSE.txt for details.
