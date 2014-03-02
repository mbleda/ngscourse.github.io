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

# Preliminaries

In this hands-on will learn how to align DNA and RNA-seq data with most widely used software. 

### NGS aligners used:

- [BWA]: BWA is a software package for mapping DNA low-divergent sequences against a large reference genome, such as the human genome.
- [Bowtie2]: Bowtie 2 is an ultrafast and memory-efficient tool for aligning DNA sequencing reads to long reference sequences.
- [TopHat]: TopHat is a fast splice junction mapper for RNA-Seq reads. It aligns RNA-Seq reads to mammalian-sized genomes using the ultra high-throughput short read aligner Bowtie, and then analyzes the mapping results to identify splice junctions between exons.
- [STAR]: STAR aligns RNA-seq reads to a reference genome using uncompressed suffix arrays.

### Other software used in this hands-on:
- [SAMTools]: SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.
- [dwgsim]: dwgsim can perform whole genome simulation.
- [BEERS]: BEERS is a simulation engine for generating RNA-Seq data.

### File formats explored:

- [SAM](http://samtools.sourceforge.net/SAMv1.pdf): Sequence alignment format, plain text.
- [BAM](http://www.broadinstitute.org/igv/bam): Binary and compressed version of SAM


### Data used in this practical

For this hands-on we are going to use DNA and RNA-seq simulated data from chromosome 21. Data has been already simulated using _dwgsim_ for DNA and _BEERS_ for RNA-seq. You can copy from the shared resources into your ``data`` directory for the practical. Preparing the data directory

    mkdir data
    cd data
    cp path-to_data/* .

The name of the folders describe the type of data, ie. ```dna_chr21_100_high``` stands for: _DNA_ data from _chromosome 21_ with _100_nt read lengths of _high_ quality. Where _high_ quality means 0.1% mutations and _low_ quality 1% mutations.

If you want to learn how to simulate DNA and RNA-seq go down to the end of this tutorial.


# Exercise 1: DNA aligment

<!-- Go to the directory where you have downoaded your data: 

    cd my_visual_data_dir  

In the following **folder** you wil find mapped sequencing data from a CEU trio (father, mother and child) from the 1000 Genomes Project:

    cd ~/ngscourse.github.io/COURSE_EXAMPLE_DATA/visualization/example_1
    
    ll

These datasets contain reads only for the [GABBR1](http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000204681;r=6:29523406-29601753) gene.
-->

### BWA
./bwa index ../../data/Homo_sapiens.GRCh37.75.dna.chromosome.21.fa

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
    

3. Local realignment around INDELS (using GATK)
--------------------------------------------------------------------------------

There are 2 steps to the realignment process:

1. Create a target list of intervals which need to be realigned

    #java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R f000-reference.fa -I f010-paired_end_noDup.bam -known gold_indels.vcf -o forIndelRealigner.intervals
    java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R f000-reference.fa -I f010-paired_end_noDup.bam -o f020-indelRealigner.intervals

2. Perform realignment of the target intervals

    java -jar GenomeAnalysisTK.jar -T IndelRealigner -R f000-reference.fa -I f010-paired_end_noDup.bam -targetIntervals f020-indelRealigner.intervals -o f030-paired_and_noDup_realigned.bam 

This creates a file called ``f030-paired_and_noDup_realigned.bam`` containing all the original reads, but with better local alignments in the regions that were realigned.


4. Base quality score recalibration (using GATK)
--------------------------------------------------------------------------------

1. Analyze patterns of covariation in the sequence dataset

    java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R reference.fa -I realigned_reads.bam -L 20 -knownSites dbsnp.vcf -knownSites gold_indels.vcf -o recal_data.table 

This creates a GATKReport file called recal_data.grp containing several tables. These tables contain the covariation data that will be used in a later step to recalibrate the base qualities of your sequence data.

It is imperative that you provide the program with a set of known sites, otherwise it will refuse to run. The known sites are used to build the covariation model and estimate empirical base qualities. For details on what to do if there are no known sites available for your organism of study, please see the online GATK documentation.

2. Do a second pass to analyze covariation remaining after recalibration

    java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R reference.fa -I realigned_reads.bam -L 20 -knownSites dbsnp.vcf -knownSites gold_indels.vcf -BQSR recal_data.table -o post_recal_data.table 

This creates another GATKReport file, which we will use in the next step to generate plots. Note the use of the -BQSR flag, which tells the GATK engine to perform on-the-fly recalibration based on the first recalibration data table.

3. Apply the recalibration to your sequence data

    java -jar GenomeAnalysisTK.jar -T PrintReads -R reference.fa -I realigned_reads.bam -L 20 -BQSR recal_data.table -o recal_reads.bam

This creates a file called recal_reads.bam containing all the original reads, but now with exquisitely accurate base substitution, insertion and deletion quality scores. By default, the original quality scores are discarded in order to keep the file size down. However, you have the option to retain them by adding the flag –emit_original_quals to the PrintReads command, in which case the original qualities will also be written in the file, tagged OQ.

Notice how this step uses a very simple tool, PrintReads, to apply the recalibration. What’s happening here is that we are loading in the original sequence data, having the GATK engine recalibrate the base qualities on-the-fly thanks to the -BQSR flag (as explained earlier), and just using PrintReads to write out the resulting data to the new file.


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


