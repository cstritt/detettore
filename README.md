#### detettore ####
===================

A program to detect and characterize transposable element (TE) polymorphisms using
reference-aligned paired-end reads. 

detettore searches for:
	a) TE insertion polymorphisms (TIPs), i.e. TEs absent in the reference genome but present in a sequenced individual
	b) TE absence polymorphisms (TAPs), i.e. TEs present in the reference but absent in a sequenced individual

With sequencing data for multiple individuals, detettore can summarize variants in a vcf file, which is a standard input format
for many population genetic tools. 



REQUIREMENTS
============

BLAST+ is the only external dependency (https://blast.ncbi.nlm.nih.gov).

Minimal data requirements:
	TIPs:
		Illumina paired-end reads aligned to a reference genome
		TE (consensus) sequences, available online (https://mobilednajournal.biomedcentral.com/databases)
		reference genome
	TAPs:
		TE annotation in bed or gff
		reference genome

Note on the paired-end reads:
It is assumed that the paired-end reads were aligned to a reference genome with BWA-MEM (Li 2013). 
The reason for this restriction is that detettore uses the AS (alignment score) and XS (secondary alignment score) tags in the bam file to determine whether a read maps uniquely to the reference genome. Other aligners might not produce these optional tags. It should be easy, however, to adapt the source code to different tags.



INSTALLATION
============

Detettore is written in Python 2 

Run the setup script by typing 
    	$ python setup.py install

This installs the necessary Python packages. If Python 3 is the system default ($ python --version), detettore can be run on a virtual environment. Create one using these commands:

    	$ virtualenv -p python2.7 <location>
    	$ source <location>/bin/activate

Or when using Anaconda:
	$ conda create -n myenv python=2.7
	$ source activate myenv

Then run setup.py on the activated environment.



USAGE
=====

The workflow of detettore is as follows:

	- detect candidate TIPs and TAPs in single individuals by by running detettore.py 
	- explore results in order to find reasonable filtering criteria (important!)
	- use variantcaller.py to apply filters and summarize TE polymorphisms for multiple individuals in a vcf file
	- alternativelly, if only single individuals are considered, use filter.py to get filtered results for individual samples


Example
-------

detettore.py \
  
# required arguments
  -o xy				# Output folder name. Ideally the sample ID
  -b xy.bam \			# Unfiltered (!) bam file, i.e. including the not nicely mapping reads
  -m tips taps \		# The module to run. Can be both at the same time or just one.
  -t TEconsensus.fasta \	
  -r ref.fa \	
  -a ref_TEannot_.gff \		# Bed format should work too

# filtering thresholds with defaults
  -u 30 \			# Minimum difference between primary and secondary alignment score. Reads above the threshold are considered as mapping uniquely to the reference
  -l 30 \			# Minimum length of soft-clipped read parts 
  -id 80 \			# Minimum sequence identity between reads and TE consensus sequences
					
# other optional parameters
  -c 4 \			# Number of CPUs. The blast search can be run with multiple CPUs, as well as a time-consuming loop in the TAP module
  -keep \			# Keep fasta files with discordant read and clipped sequences, as well as the output of the blast search of theses seqences against the TE consensus sequences
  -reconstruct			# If the TAP module is used, TE sequences present in the reference and the sequenced individual can be reconstructed and written to a fasta file. Very slow for large annotations!


This will run the the TIP and the TAP analysis and create the output folder xy with 4 files:
	xy_bamstats.txt			# bam read statistics	
	TIPs_candidates.tsv		# candidate non-reference TE insertions
	TIPs_supporting_reads.txt	# names of reads supporting the TIP candidates
	TAPs_candidates.tsv		# candidate TAPS

	
	
STRUCTURE OF OUTPUT FILES
=========================

An advantage of detettore compared to other TE detection tools is the amount of information in its raw output. Detecting TE polymorphisms relies on a series of ad-hoc criteria because, in contrast to SNP calling, there is neither a mutational model nor any way to model errors. To avoid hard-coding decisions which might work well in one study system but not others, detettore returns a spectrum of TIP candidates from representing very messy signatures to neat ones. In our experience, looking at this spectrum (using IGV, for example) is crucial to settle upon reasonable filteria criteria.



TIPs_candidates.tsv
-------------------

  chromosome
  position
  TE
  strand		
  
  nr_supporting_reads	# Total number of supporting reads. Evidence from the 5' and the 3' site is separated with a dash.
  nr_discordant_reads
  nr_splitreads
  
  nr_aligned_positions	# Nr of positions in the TE consensus sequence that are covered by the discordant reads and split reads reaching into the insertion
  TSD			# The target site duplication
  nr_bridge_reads	# Number of reads which bridge the insertion breakpoint and point to heterozygous insertions or false positives
  
  region_start		# Coordinates of the region surrounding the TE insertion
  region_end	
  region_cov_mean	# Mean coverage in this region. Useful indicator of messiness and a good filtering criterion.
  region_cov_stdev	



TAPs_candidates.tsv
-------------------
  chromosome
  TE_start		# Information about the reference TE insertion
  TE_end
  TE_id
  TE_length

  isize			# Mean insert size of the deviant read-pairs around the TE
  nr_discordant		# Nr of deviant read-pairs
  absence_start		# Beginning of the gap
  absence_end		# End of the gap
  
  region_start		# Same as for TIPs
  region_end
  region_cov_mean
  region_cov_stdev



FILTERING AND VARIANT CALLING
=============================

Single TIP candidate files can be filtered with filter.py. The same filters are implemented in variantcaller.py, which takes a list of TIP candidate
files as input. Such a list can be generated by typing 

	$ find $PWD -type f -name "TIPs_candidates.tsv"

in the project directory. Importantly, the output folder names will be the sample IDs in the vcf file.



LICENSE
=======
GNU General Public License v3. See LICENSE.txt for details.


