% [NGS data analysis course](http://ngscourse.github.io/)
% __Variant calling__
% _(updated 28-02-2014)_

<!-- COMMON LINKS HERE -->

[BWA]: http://bio-bwa.sourceforge.net/ "BWA"
[SAMTools]: http://samtools.sourceforge.net/ "samtools"
[Picard]: http://picard.sourceforge.net/ "Picard"
[gatk-site]: http://www.broadinstitute.org/gatk/ "GATK"

Preliminaries
================================================================================

Software used in this practical:
--------------------------------

- [BWA] : BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome.
- [SAMTools] : SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.
- [Picard] : Picard comprises Java-based command-line utilities that manipulate SAM files, and a Java API (SAM-JDK) for creating new programs that read and write SAM files.
- [GATK] : (Genome Analysis Toolkit): A package to analyze next-generation re-sequencing data, primary focused on variant discovery and genotyping.


File formats explored:
----------------------

- [SAM](http://samtools.sourceforge.net/SAMv1.pdf)
- [BAM](http://www.broadinstitute.org/igv/bam)
- VCF Variant Call Format: see [1000 Genomes][vcf-format-1000ge] and [Wikipedia][vcf-format-wikipedia] specifications.


<!-- Data used in this practical
-------------------------------

- [mirbase_mature.fa](../../../COURSE_EXAMPLE_DATA/f010_mirbase_mature.fa): mature micro RNAs downloaded form mirbase

You can download them or copy them to your ``data`` directory for the practical

clean directory

    rm -r data
    mkdir data
    cd data
    cp ../../../../COURSE_EXAMPLE_DATA/f010_mirbase_mature.fa .


\ 

Find all data files for the course here: [COURSE_EXAMPLE_DATA](../../../COURSE_EXAMPLE_DATA)
-->


Exercise 2: Variant calling with single-end data
================================================================================

<!-- Go to the directory where you have downoaded your data: 

    cd my_visual_data_dir  

In the following **folder** you wil find mapped sequencing data from a CEU trio (father, mother and child) from the 1000 Genomes Project:

    cd ~/ngscourse.github.io/COURSE_EXAMPLE_DATA/variant_calling/
    
    ll

These datasets contain reads only for the [GABBR1](http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000204681;r=6:29523406-29601753) gene.
-->

1. Prepare reference genome: generate the BWA index
--------------------------------------------------------------------------------

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


