% [NGS data analysis course](http://ngscourse.github.io/)
% __Quality control & Adapter trimming__
% _(updated 02-03-2014)_


<!-- Common URLs: Tools -->

[fastqc]:http://www.bioinformatics.babraham.ac.uk/projects/fastqc "FastQC home page"
[cutadapt]:http://code.google.com/p/cutadapt "cutadapt home page"

<!-- Common URLs: File Formats -->

[fastq-format-wikipedia]:http://en.wikipedia.org/wiki/FASTQ_format  "FastQC in Wikipedia"
[fastq-format-nar]:http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217 "FastQC original NAR publication"

<!-- External URLs -->

[mirbase]:http://www.mirbase.org/ "miRBase: a searchable database of published miRNA"


Preliminaries
================================================================================

Software used in this practical:
--------------------------------

- [FastQC] : A quality control tool for high-throughput sequence data.
- [cutadapt] : A tool to remove adapter sequences from high-throughput sequencing data.

File formats explored:
----------------------

__FastQ__. See: 

- [Wikipedia][fastq-format-wikipedia].
- [NAR 2010][fastq-format-nar].

Data used in this practical
-------------------------------

- __f000_raw_mirna.fastq__: RNA-Seq of a microRNA sample.


Overview
================================================================================

1. Use [FastQC] to explore the raw data.
1. Use [cutadapt] to remove adapters.
1. Use [cutadapt] to filter reads based on quality.
1. Use [FastQC] to explore the filtered data.


Exercise
================================================================================

Create a new directory for the exercise and copy the raw data to it: 

<!-- clean directory
    rm -r data
-->

    mkdir data
	cd data
    cp Path_To_Your/COURSE_EXAMPLE_DATA/f000_raw_mirna.fastq .  # The final dot means "here", in the working directory.

<!--
    cp ../../../../COURSE_EXAMPLE_DATA/f000_raw_mirna.fastq .
-->


Explore the raw data using some Linux shell commands
--------------------------------------------------------------------------------

The file __f000_raw_mirna.fastq__ contains reads form a microRNA sequencing experiment. 
Use the command `head` to have a view of the first lines of the file:

	head f000_raw_mirna.fastq

Use `wc` to count how many reads there are in the file (remember you have to divide by 4)
	
    wc -l f000_raw_mirna.fastq


Explore the raw data quality using FastQC
--------------------------------------------------------------------------------

First create a directory to store the results of the fastqc analysis:

	mkdir f010_res_fastqc

Then execute fastqc storing the results in the created directory (option `-o`):

	fastqc -o f010_res_fastqc f000_raw_mirna.fastq

Find the results in the __fastqc_report.html__ file and discus them.

There are many _Overrepresented sequences_. 
Explore whether some of them correspond to miRNAs using [miRBase] search utilities (Search -> By sequence).



Handling adapters
--------------------------------------------------------------------------------

There are 2 known adapters used in this experiment: 

    CTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA
    CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC

Use the command grep to see whether they are still present in your data:

    grep "CTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA" f000_raw_mirna.fastq 
	grep "CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC" f000_raw_mirna.fastq 

Do the sequences appear systematically at the beginning or at the end of the reads?

\ 

But the adapters could also appear in the _reverse_, _complementary_ or _reverse complementary_ mode.

Compute the _reverse_, _complementary_ and the _reverse complementary_ sequences of the two adapters,
and find out which of them appear in your data.

To compute those sequences you can use some online resources as the one in:  
<http://www.bioinformatics.org/sms/rev_comp.html>


<!--
Or to use R-Bioconductor to compute their reverse, complementary and reverse complementary.

library (Biostrings)
myseq <- DNAString ("CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC")
reverse (myseq)
complement (myseq)
reverseComplement (myseq)

-->

Use grep form Linux shell to find out which of the versions of the adapter is in your data.


### Adapter 1

    grep CTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA f000_raw_mirna.fastq | wc -l  ## present in the sample (at the beginning of the reads)
	
    grep GACCCTTTAGTGGTATTTGCACTTTACAGAAACCTAAACCCTTAGAATATTCAAGACATACTCTGGTGAGATTTTT f000_raw_mirna.fastq | wc -l 
	
    grep TTTTTAGAGTGGTCTCATACAGAACTTATAAGATTCCCAAATCCAAAGACATTTCACGTTTATGGTGATTTCCCAG f000_raw_mirna.fastq | wc -l  ## present in the sample (at the end of the read) ... but not so numerous
	
    grep AAAAATCTCACCAGAGTATGTCTTGAATATTCTAAGGGTTTAGGTTTCTGTAAAGTGCAAATACCACTAAAGGGTC f000_raw_mirna.fastq | wc -l 

But sometimes the adapter does not appear complete. 
It may be there just the first part:

    grep CTGGGAAATCACCATAAACGTGAAATGTCTTTGGA f000_raw_mirna.fastq | wc -l 
	
    grep GACCCTTTAGTGGTATTTGCACTTTACAGAAACCT f000_raw_mirna.fastq | wc -l 
	
    grep TTTTTAGAGTGGTCTCATACAGAACTTATAAGATT f000_raw_mirna.fastq | wc -l 
	
    grep AAAAATCTCACCAGAGTATGTCTTGAATATTCTAA f000_raw_mirna.fastq | wc -l 

or the end part of it:

    grep GAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA f000_raw_mirna.fastq | wc -l 
	
    grep CTTAGAATATTCAAGACATACTCTGGTGAGATTTTT f000_raw_mirna.fastq | wc -l 
	
    grep ATCCAAAGACATTTCACGTTTATGGTGATTTCCCAG f000_raw_mirna.fastq | wc -l 
	
    grep TAGGTTTCTGTAAAGTGCAAATACCACTAAAGGGTC f000_raw_mirna.fastq | wc -l 


### Adapter 2

    grep CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC f000_raw_mirna.fastq | wc -l   ## present in the sample (at the end of the read) ... but not so numerous
	
    grep GAAAAAAAGCAGGAAAGGTGTTCTATATATTTCGGTTCTTTAGCTTTATGAAAGTTCAATGCCATTCG f000_raw_mirna.fastq | wc -l 
	
    grep GCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAAAAAAG f000_raw_mirna.fastq | wc -l   ## present in the sample (at the beginning of the reads)
	
    grep CGAATGGCATTGAACTTTCATAAAGCTAAAGAACCGAAATATATAGAACACCTTTCCTGCTTTTTTTC f000_raw_mirna.fastq | wc -l 

As before, sometimes the adapter does not appear complete. 
It may be there just the first part:

    grep CTTTTTTTCGTCCTTTCCACAAGATATATA f000_raw_mirna.fastq | wc -l 
	
    grep GAAAAAAAGCAGGAAAGGTGTTCTATATAT f000_raw_mirna.fastq | wc -l 
	
    grep GCTTACCGTAACTTGAAAGTATTTCGATTT f000_raw_mirna.fastq | wc -l 
	
    grep CGAATGGCATTGAACTTTCATAAAGCTAAA f000_raw_mirna.fastq | wc -l 

or the end part of it:

    grep AAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC f000_raw_mirna.fastq | wc -l 
	
    grep TTCGGTTCTTTAGCTTTATGAAAGTTCAATGCCATTCG f000_raw_mirna.fastq | wc -l 
	
    grep CTTGGCTTTATATATCTTGTGGAAAGGACGAAAAAAAG f000_raw_mirna.fastq | wc -l 
	
    grep GAACCGAAATATATAGAACACCTTTCCTGCTTTTTTTC f000_raw_mirna.fastq | wc -l 



Use cutadapt to make an adapter trimming of the reads.
--------------------------------------------------------------------------------

Check the options:

- __-a__ for adapter to the 3' end.
- __-g__ for adapter to the 5' end.


To get read of the the adapters found in our data we run [cutadapt] several times:

    cutadapt -g CTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA -o f020_mirna_trim1.fastq f000_raw_mirna.fastq
	
    cutadapt -a TTTTTAGAGTGGTCTCATACAGAACTTATAAGATTCCCAAATCCAAAGACATTTCACGTTTATGGTGATTTCCCAG -o f020_mirna_trim2.fastq f020_mirna_trim1.fastq
	
    cutadapt -g GCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAAAAAAG         -o f020_mirna_trim3.fastq f020_mirna_trim2.fastq
	
    cutadapt -a CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC         -o f020_mirna_trim4.fastq f020_mirna_trim3.fastq


Now you can `grep` again searching for the adapters

### Adapter 1

    grep CTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA f020_mirna_trim4.fastq | wc -l
	
    #grep GACCCTTTAGTGGTATTTGCACTTTACAGAAACCTAAACCCTTAGAATATTCAAGACATACTCTGGTGAGATTTTT f020_mirna_trim4.fastq | wc -l 
	
    grep TTTTTAGAGTGGTCTCATACAGAACTTATAAGATTCCCAAATCCAAAGACATTTCACGTTTATGGTGATTTCCCAG f020_mirna_trim4.fastq | wc -l
	
    #grep AAAAATCTCACCAGAGTATGTCTTGAATATTCTAAGGGTTTAGGTTTCTGTAAAGTGCAAATACCACTAAAGGGTC f020_mirna_trim4.fastq | wc -l 

### Adapter 2

    grep CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC f020_mirna_trim4.fastq | wc -l
	
    #grep GAAAAAAAGCAGGAAAGGTGTTCTATATATTTCGGTTCTTTAGCTTTATGAAAGTTCAATGCCATTCG f020_mirna_trim4.fastq | wc -l 
	
    grep GCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAAAAAAG f020_mirna_trim4.fastq | wc -l
	
    #grep CGAATGGCATTGAACTTTCATAAAGCTAAAGAACCGAAATATATAGAACACCTTTCCTGCTTTTTTTC f020_mirna_trim4.fastq | wc -l 



Explore the quality of the trimmed file using FastQC
--------------------------------------------------------------------------------

Check the data quality again using fastqc:

	mkdir f030_res_fastqc_trimmed
	fastqc -o f030_res_fastqc_trimmed f020_mirna_trim4.fastq


Some of the reads seems to be too short and some others may not have enough quality. 


Use cutadapt to filter reads by quality and length. 
--------------------------------------------------------------------------------

Check the options:

- __-q__ quality cutoff.
- __-m__ minimum length.
- __-M__ maximum length.


Run cutadapt for length and quality purge of the reads.

    cutadapt -m 17 -M 30 -q 10 -o f030_mirna_cut.fastq f020_mirna_trim4.fastq


Check the data quality again using fastqc:

	mkdir f040_res_fastqc_trimmed_purged
	
	fastqc -o f040_res_fastqc_trimmed_purged f030_mirna_cut.fastq


Explore again the _Overrepresented sequences_ in [mirbase] (Search -> By sequence).

Count how many reads are left for the analysis (divide by 4)

    wc -l f000_raw_mirna.fastq
	
    wc -l f030_mirna_cut.fastq
