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


Exercise 2: Variant calling with single-end data
================================================================================

Go to the variant_calling folder in your course directory: 

    cd <my_course_directory>/variant_calling


1. Prepare reference genome: generate the fasta file index
--------------------------------------------------------------------------------

This step is no longer needed since we have already done it in [example1](http://ngscourse.github.io/COURSE_MATERIALS/variant_calling/example/010_example.html)

2. Prepare BAM file
--------------------------------------------------------------------------------

Go to the example1 folder:

    cd example2

Add the read group tag to our BAM file using ``Picard``:

    AddOrReplaceReadGroups.jar I=f000-dna_100_high_se.bam O=f010-dna_100_high_se_fixRG.bam RGID=group2 RGLB=lib2 RGPL=illumina RGSM=sample2 RGPU=unit2

We must sort the BAM file using ``samtools``:

    samtools sort f010-dna_100_high_se_fixRG.bam f020-dna_100_high_se_fixRG_sorted

Index the BAM file:

    samtools index f020-dna_100_high_se_fixRG_sorted.bam


3. Mark duplicates (using Picard)
--------------------------------------------------------------------------------

Mark and remove duplicates:

    MarkDuplicates.jar INPUT=f020-dna_100_high_se_fixRG_sorted.bam OUTPUT=f030-dna_100_high_se_fixRG_sorted_noDup.bam METRICS_FILE=metrics.txt

Index the new BAM file:

    BuildBamIndex.jar INPUT=f030-dna_100_high_se_fixRG_sorted_noDup.bam


4. Local realignment around INDELS (using GATK)
--------------------------------------------------------------------------------

There are 2 steps to the realignment process:

Create a target list of intervals which need to be realigned
  
    GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../f000-reference.fa -I f030-dna_100_high_se_fixRG_sorted_noDup.bam -o f040-indelRealigner.intervals

Perform realignment of the target intervals

    GenomeAnalysisTK.jar -T IndelRealigner -R ../f000-reference.fa -I f030-dna_100_high_se_fixRG_sorted_noDup.bam -targetIntervals f040-indelRealigner.intervals \
    -o f040-dna_100_high_se_fixRG_sorted_noDup_realigned.bam 


5. Base quality score recalibration (using GATK)
--------------------------------------------------------------------------------

Two steps:

Analyze patterns of covariation in the sequence dataset

    GenomeAnalysisTK.jar -T BaseRecalibrator -R ../f000-reference.fa -I f040-dna_100_high_se_fixRG_sorted_noDup_realigned.bam -knownSites ../f000-dbSNP_chr21.vcf \  
    -o f050-recalibration_data.table

Apply the recalibration to your sequence data

    GenomeAnalysisTK.jar -T PrintReads -R ../f000-reference.fa -I f040-dna_100_high_se_fixRG_sorted_noDup_realigned.bam -BQSR f050-recalibration_data.table \  
    -o f050-dna_100_high_se_fixRG_sorted_noDup_realigned_recalibrated.bam


6. Variant calling (using GATK - **UnifiedGenotyper**)
--------------------------------------------------------------------------------

**SNP calling**

    GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../f000-reference.fa -I f050-dna_100_high_se_fixRG_sorted_noDup_realigned_recalibrated.bam -glm SNP \  
    -o f060-dna_100_high_se_snps.vcf

**INDEL calling**

    GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../f000-reference.fa -I f050-dna_100_high_se_fixRG_sorted_noDup_realigned_recalibrated.bam -glm INDEL \  
    -o f060-dna_100_high_se_indels.vcf


7. Compare paired-end VCF against single-end VCF
--------------------------------------------------------------------------------

Open IGV and load a the single-end VCF we have generated in this tutorial (``f060-dna_100_high_se_snps.vcf``) and the the paired-end VCF generated in example1 (``f060-dna_100_high_pe_snps.vcf``).




