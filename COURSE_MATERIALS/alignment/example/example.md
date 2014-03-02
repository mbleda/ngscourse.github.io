% [NGS data analysis course](http://ngscourse.github.io/)
% __DNA and RNA-seq NGS alignment__
% _(updated 28-02-2014)_

<!-- COMMON LINKS HERE -->

[BWA]: http://bio-bwa.sourceforge.net/ "BWA"
[Bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml "Bowtie2"
[TopHat]: http://tophat.cbcb.umd.edu/ "TopHat"
[STAR]: https://code.google.com/p/rna-star/ "STAR"
[SAMTools]: http://samtools.sourceforge.net/ "SAMtools"
[dwgsim]: http://sourceforge.net/apps/mediawiki/dnaa/index.php?title=Whole_Genome_Simulation "dwgsim"
[BEERS]: http://www.cbil.upenn.edu/BEERS/ "BEERS"
[Ensembl]: http://www.ensembl.org/index.html "Ensembl"

# Preliminaries

In this hands-on will learn how to align DNA and RNA-seq data with most widely used software today. Building a whole genome index requires a lot of RAM memory and almost one hour in a typical workstation, for this reason in this tutorial we will work with chromosome 21 to speed up the exercises. The same steps will be done for a whole genome alignment.

### NGS aligners used:

- [BWA]: BWA is a software package for mapping **DNA** low-divergent sequences against a large reference genome, such as the human genome.
- [Bowtie2]: *Bowtie 2* is an ultrafast and memory-efficient tool for aligning **DNA** sequencing reads to long reference sequences.
- [TopHat]: *TopHat* is a fast splice junction mapper for RNA-Seq reads. It aligns **RNA-Seq** reads to mammalian-sized genomes using the ultra high-throughput short read aligner Bowtie, and then analyzes the mapping results to identify splice junctions between exons.
- [STAR]: *STAR* aligns **RNA-seq** reads to a reference genome using uncompressed suffix arrays.

### Other software used in this hands-on:
- [SAMTools]: SAM Tools **provide various utilities** for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.
- [dwgsim]: dwgsim can perform whole **genome simulation**.
- [BEERS]: BEERS is a **simulation engine** for generating **RNA-Seq** data.

### File formats explored:

- [SAM](http://samtools.sourceforge.net/SAMv1.pdf): Sequence alignment format, plain text.
- [BAM](http://www.broadinstitute.org/igv/bam): Binary and compressed version of SAM


### Data used in this practical

Create a ```data``` folder in your *working directory* to store both the *reference genome* to be used (human chromosome 21) and *simulated datasets*:

    mkdir data

##### Download reference genome

Working with NGS data requires a high-end workstations and time for building the reference genome indexes and alignment. During this tutorial we will work only with chromosome 21 to speed up the runtimes. Go to the *Download* link at the top of [Ensembl] website and then to *Download data via FTP*, you go in only one step by going to:

    http://www.ensembl.org/info/data/ftp/index.html

You should see a species table with a Human (*Homo sapiens*) row and a *DNA (FASTA)* column, just go there and download the chromosome 21 (*Homo_sapiens.GRCh37.75.dna.chromosome.21.fa.gz*) and move it to your ```data``` folder:

    mv Homo_sapiens.GRCh37.75.dna.chromosome.21.fa.gz path_to_local_data

**NOTE:** For working with the whole genome the file to be downloaded is **Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz**	


##### Copy simulated datasets

For this hands-on we are going to use small DNA and RNA-seq datasets simulated from chromosome 21. Data has been already simulated using _dwgsim_ software from SAMtools for DNA and _BEERS_ for RNA-seq. You can copy from the shared resources into your ``data`` directory for the practical. Preparing the data directory:

    cp path_to_shared_data/* your_local_data/

The name of the folders and files describe the dataset, ie. ```dna_chr21_100_high``` stands for: _DNA_ type of data from _chromosome 21_ with _100_nt read lengths of _high_ quality. Where _high_ quality means 0.1% mutations and _low_ quality 1% mutations.

**NOTE:** If you want to learn how to simulate DNA and RNA-seq for other conditions go down to the end of this tutorial.


# Exercise 1: NGS Genomic DNA aligment

In this exercise we'll learn how to download, install, build the reference genome index and align in single-end and paired-end mode with the two most widely DNA aligners.


### BWA
[BWA] is probably the most used aligner for DNA. AS the documentation states it consists of three different algorithms: *BWA*, *BWA-SW* and *BWA-MEM*. The first algorithm, which is the oldest, is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA for 70-100bp Illumina reads.

All these three algorithms come in the same binary so only one download and installation is needed.

##### Download and install
You can click on ```SF download page``` link in the [BWA] page or click directly to:
    http://sourceforge.net/projects/bio-bwa/files/

Click in the last version of BWA and wait for a few seconds, the download will start. 


##### Build the index

    ./bwa index ../../data/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa

##### Aligning in SE and PE modes



    ./bwa aln index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 ../../data/dna_chr21_100_high/dna_chr21_100_high.bwa.read1.fastq -f ../../alignments/bwa/dna_chr21_100_high_se.sai
    ./bwa samse index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa ../../alignments/bwa/dna_chr21_100_high_se.sai ../../data/dna_chr21_100_high/dna_chr21_100_high.bwa.read1.fastq -f ../../alignments/bwa/dna_chr21_100_high_se.sam

    ./bwa aln index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 ../../data/dna_chr21_100_high/dna_chr21_100_high.bwa.read1.fastq -f ../../alignments/bwa/dna_chr21_100_high_pe1.sai
    ./bwa aln index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 ../../data/dna_chr21_100_high/dna_chr21_100_high.bwa.read2.fastq -f ../../alignments/bwa/dna_chr21_100_high_pe2.sai
    ./bwa sampe index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa ../../alignments/bwa/dna_chr21_100_high_pe1.sai ../../alignments/bwa/dna_chr21_100_high_pe2.sai ../../data/dna_chr21_100_high/dna_chr21_100_high.bwa.read1.fastq ../../data/dna_chr21_100_high/dna_chr21_100_high.bwa.read2.fastq -f ../../alignments/bwa/dna_chr21_100_high_pe.sam


 1036  ./bwa aln index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 ../../data/dna_chr21_100_low/dna_chr21_100_low.bwa.read1.fastq -f ../../alignments/bwa/dna_chr21_100_low_se.sai
 1037  ./bwa aln index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 ../../data/dna_chr21_100_low/dna_chr21_100_low.bwa.read1.fastq -f ../../alignments/bwa/dna_chr21_100_low_pe1.sai
 1038  ./bwa aln index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa -t 4 ../../data/dna_chr21_100_low/dna_chr21_100_low.bwa.read2.fastq -f ../../alignments/bwa/dna_chr21_100_low_pe2.sai
 1039  ./bwa samse index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa ../../alignments/bwa/dna_chr21_100_low_se.sai ../../data/dna_chr21_100_low/dna_chr21_100_low.bwa.read1.fastq -f ../../alignments/bwa/dna_chr21_100_low_se.sam
 1040  ./bwa sampe index/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa ../../alignments/bwa/dna_chr21_100_low_pe1.sai ../../alignments/bwa/dna_chr21_100_low_pe2.sai ../../data/dna_chr21_100_low/dna_chr21_100_low.bwa.read1.fastq ../../data/dna_chr21_100_low/dna_chr21_100_low.bwa.read2.fastq -f ../../alignments/bwa/dna_chr21_100_low_pe.sam

samtools view -S -b dna_chr21_100_high_pe.sam -o dna_chr21_100_high_pe.bam

##### 1. Prepare reference genome: generate the BWA index


Use ``BWA`` to index the the reference genome:

    bwa index -a bwtsw f000-reference.fa

where ``-a bwtsw`` specifies that we want to use the indexing algorithm that is capable of handling the whole human genome.


Use ``SAMTools`` to generate the fasta file index:

    samtools faidx f000-reference.fa

This creates a file called f000-reference.fa.fai, with one record per line for each of the contigs in the FASTA reference file.

Generate the sequence dictionary using ``Picard``:

    java -jar CreateSequenceDictionary.jar REFERENCE=f000-reference.fa OUTPUT=f000-reference.dict


2. Mark duplicates (using Picard)
--------------------------------------------------------------------------------

Run the following **Picard** command to mark duplicates:

    java -jar MarkDuplicates.jar INPUT=f000-paired_end.bam OUTPUT=f010-paired_end_noDup.bam METRICS_FILE=metrics.txt

This creates a sorted BAM file called ``f010-paired_end_noDup.bam`` with the same content as the input file, except that any duplicate reads are marked as such. It also produces a metrics file called ``metrics.txt`` containing (can you guess?) metrics.

Run the following **Picard** command to index the new BAM file:

    java -jar BuildBamIndex.jar INPUT=f010-paired_end_noDup.bam

Q1. How many reads are removed as duplicates from the files (hint view the on-screen output from the two commands)?
    
# Exercise 2: NGS RNA-seq aligment


4. Variant calling (using GATK - UnifiedGenotyper)
--------------------------------------------------------------------------------

SNPs and INDELS are called using separate instructions.

1. SNP calling

    java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R f000_reference.fa -I 7_realigned_aligned.bam -glm SNP -o 8_snp_variants.vcf

2. INDEL calling

    java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R 0_reference.fa -I 7_realigned_aligned.bam -glm INDEL -o 8_indel_variants.vcf


5. Introduce filters in the VCF file
--------------------------------------------------------------------------------

Example: filter SNPs with low confidence calling (QD < 12.0) and flag them as "LowConf".

    java -jar GenomeAnalysisTK.jar -T VariantFiltration -R 0_reference.fa -V 8_snp_variants.vcf --filterExpression "QD < 12.0" --filterName "LowConf" -o 9_snp_filtered.vcf


# Simulating NGS datasets

### DNA

    ./dwgsim-0.1.11/dwgsim -1 100 -2 100 -y 0 -N 1000000 ../data/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa ../data/dna_chr21_100_high/dna_chr21_100_high
    ./dwgsim-0.1.11/dwgsim -1 100 -2 100 -y 0 -N 1000000 -r 0.01 ../data/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa ../data/dna_chr21_100_low/dna_chr21_100_low


### RNA-seq


