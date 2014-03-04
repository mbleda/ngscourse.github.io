% [NGS data analysis course](http://ngscourse.github.io/)
% __Variant calling__
% _(updated 28-02-2014)_

<!-- COMMON LINKS HERE -->

[BWA]: http://bio-bwa.sourceforge.net/ "BWA"
[SAMTools]: http://samtools.sourceforge.net/ "samtools"
[Picard]: http://picard.sourceforge.net/ "Picard"
[MuTect]: http://www.broadinstitute.org/cancer/cga/mutect_download "MuTect"

Preliminaries
================================================================================

Software used in this practical:
--------------------------------

- [BWA] : BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome.
- [SAMTools] : SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.
- [Picard] : Picard comprises Java-based command-line utilities that manipulate SAM files, and a Java API (SAM-JDK) for creating new programs that read and write SAM files.
- [MuTect] : method developed at the Broad Institute for the reliable and accurate identification of somatic point mutations in next generation sequencing data of cancer genomes.


File formats explored:
----------------------

- [SAM](http://samtools.sourceforge.net/SAMv1.pdf)
- [BAM](http://www.broadinstitute.org/igv/bam)
- VCF Variant Call Format: see [1000 Genomes][vcf-format-1000ge] and [Wikipedia][vcf-format-wikipedia] specifications.

Some previous configurations
================================================================================

Download MuTect software from http://www.broadinstitute.org/cancer/cga/mutect_download

You need to register before downloading. Extract the file in your desktop.


Just to run MuTect easier:

    echo "alias MuTect.jar='java -jar /home/Desktop/muTect-1.1.4-bin/muTect-1.1.4.jar'" >> ~/.bashrc
    source ~/.bashrc


Exercise 3: Somatic calling
================================================================================

1. Prepare BAM file
--------------------------------------------------------------------------------

Go to the example1 folder:

    cd example3

Add read groups:

    AddOrReplaceReadGroups.jar I=f000-normal.bam O=f010-normal_fixRG.bam RGID=group1 RGLB=lib1 RGPL=illumina RGSM=sample1 RGPU=unit1
    AddOrReplaceReadGroups.jar I=f000-tumor.bam O=f010-tumor_fixRG.bam RGID=group2 RGLB=lib1 RGPL=illumina RGSM=sample1 RGPU=unit2

Sort:

    samtools sort f010-normal_fixRG.bam f020-normal_fixRG_sorted
    samtools sort f010-tumor_fixRG.bam f020-tumor_fixRG_sorted

Index the BAM file:

    samtools index f020-normal_fixRG_sorted.bam
    samtools index f020-tumor_fixRG_sorted.bam


2. Somatic calling
--------------------------------------------------------------------------------

For brevity, we are not including BAM preprocessing steps. However, in real analysis it is recommended to include them.

    MuTect.jar --analysis_type MuTect --reference_sequence ../f000-reference.fa --dbsnp ../f000-dbSNP_chr21.vcf --input_file:normal f020-normal_fixRG_sorted.bam \  
    --input_file:tumor f020-tumor_fixRG_sorted.bam --out f030-call_stats.out --coverage_file f030-coverage.wig.txt

