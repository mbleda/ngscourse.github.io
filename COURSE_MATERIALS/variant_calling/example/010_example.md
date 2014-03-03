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


Some previous configurations
================================================================================

Just to run GATK easier:

    echo "alias GenomeAnalysisTK.jar='java -jar /home/Desktop/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar'" >> ~/.bashrc

To run Picard easier:

    echo "alias CreateSequenceDictionary.jar='java -jar /home/Desktop/picard-tools-1.108/picard-tools-1.108/CreateSequenceDictionary.jar'" >> ~/.bashrc
    echo "alias AddOrReplaceReadGroups.jar='java -jar /home/Desktop/picard-tools-1.108/picard-tools-1.108/AddOrReplaceReadGroups.jar'" >> ~/.bashrc
    echo "alias MarkDuplicates.jar='java -jar /home/Desktop/picard-tools-1.108/picard-tools-1.108/MarkDuplicates.jar'" >> ~/.bashrc
    echo "alias BuildBamIndex.jar='java -jar /home/Desktop/picard-tools-1.108/picard-tools-1.108/BuildBamIndex.jar'" >> ~/.bashrc


Exercise 1: Variant calling with paired-end data
================================================================================

Go to your course directory: 

    cd <my_course_directory>

Copy the data for the tutorials:

    cp -r /mounts/course_shares/Open_Share/variant_calling .
    ll

These datasets contain reads only for the chromosome 21.


1. Prepare reference genome: generate the fasta file index
--------------------------------------------------------------------------------

Use ``BWA`` to index the the reference genome:

    bwa index -a bwtsw f000-reference.fa

where ``-a bwtsw`` specifies that we want to use the indexing algorithm that is capable of handling the whole human genome.


Use ``SAMTools`` to generate the fasta file index:

    samtools faidx f000-reference.fa

This creates a file called f000-reference.fa.fai, with one record per line for each of the contigs in the FASTA reference file.


Generate the sequence dictionary using ``Picard``:

    CreateSequenceDictionary.jar REFERENCE=f000-reference.fa OUTPUT=f000-reference.dict


2. Prepare BAM file
--------------------------------------------------------------------------------

Go to the example1 folder:

    cd example1

The **read group** information is key for downstream GATK functionality. The GATK will not work without a read group tag. Make sure to enter as much metadata as you know about your data in the read group fields provided. For more information about all the possible fields in the @RG tag, take a look at the SAM specification.

    AddOrReplaceReadGroups.jar I=f000-dna_100_high_pe.bam O=f010-dna_100_high_pe_fixRG.bam RGID=group1 RGLB=lib1 RGPL=illumina RGSM=sample1 RGPU=unit1

We must sort and index the BAM file before processing it with Picard and GATK. To sort the bam file we use ``samtools``

    samtools sort f010-dna_100_high_pe_fixRG.bam f020-dna_100_high_pe_fixRG_sorted

Index the BAM file:

    samtools index f020-dna_100_high_pe_fixRG_sorted.bam


3. Mark duplicates (using Picard)
--------------------------------------------------------------------------------

Run the following **Picard** command to mark duplicates:

    MarkDuplicates.jar INPUT=f020-dna_100_high_pe_fixRG_sorted.bam OUTPUT=f030-dna_100_high_pe_fixRG_sorted_noDup.bam METRICS_FILE=metrics.txt

This creates a sorted BAM file called ``f030-dna_100_high_pe_fixRG_sorted_noDup.bam`` with the same content as the input file, except that any duplicate reads are marked as such. It also produces a metrics file called ``metrics.txt`` containing (can you guess?) metrics.

**QUESTION:** How many reads are removed as duplicates from the files (hint view the on-screen output from the two commands)?

Run the following **Picard** command to index the new BAM file:

    BuildBamIndex.jar INPUT=f030-dna_100_high_pe_fixRG_sorted_noDup.bam


4. Local realignment around INDELS (using GATK)
--------------------------------------------------------------------------------

There are 2 steps to the realignment process:

**First**, create a target list of intervals which need to be realigned
  
    GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../f000-reference.fa -I f030-dna_100_high_pe_fixRG_sorted_noDup.bam -o f040-indelRealigner.intervals

**Second**, perform realignment of the target intervals

    GenomeAnalysisTK.jar -T IndelRealigner -R ../f000-reference.fa -I f030-dna_100_high_pe_fixRG_sorted_noDup.bam -targetIntervals f040-indelRealigner.intervals \
    -o f040-dna_100_high_pe_fixRG_sorted_noDup_realigned.bam 

This creates a file called ``f040-dna_100_high_pe_fixRG_sorted_noDup_realigned.bam`` containing all the original reads, but with better local alignments in the regions that were realigned.


5. Base quality score recalibration (using GATK)
--------------------------------------------------------------------------------

Two steps:

**First**, analyze patterns of covariation in the sequence dataset

    GenomeAnalysisTK.jar -T BaseRecalibrator -R ../f000-reference.fa -I f040-dna_100_high_pe_fixRG_sorted_noDup_realigned.bam -knownSites ../f000-dbSNP_chr21.vcf \  
    -o f050-recalibration_data.table

This creates a GATKReport file called ``f050-recalibration_data.table`` containing several tables. These tables contain the covariation data that will be used in a later step to recalibrate the base qualities of your sequence data.

It is imperative that you provide the program with a set of **known sites**, otherwise it will refuse to run. The known sites are used to build the covariation model and estimate empirical base qualities. For details on what to do if there are no known sites available for your organism of study, please see the online GATK documentation.

**Second**, apply the recalibration to your sequence data

    GenomeAnalysisTK.jar -T PrintReads -R ../f000-reference.fa -I f040-dna_100_high_pe_fixRG_sorted_noDup_realigned.bam -BQSR f050-recalibration_data.table \  
    -o f050-dna_100_high_pe_fixRG_sorted_noDup_realigned_recalibrated.bam

This creates a file called ``f050-dna_100_high_pe_fixRG_sorted_noDup_realigned_recalibrated.bam`` containing all the original reads, but now with exquisitely accurate base substitution, insertion and deletion quality scores. By default, the original quality scores are discarded in order to keep the file size down. However, you have the option to retain them by adding the flag ``–emit_original_quals`` to the ``PrintReads`` command, in which case the original qualities will also be written in the file, tagged OQ.


6. Variant calling (using GATK - **UnifiedGenotyper**)
--------------------------------------------------------------------------------

SNPs and INDELS are called using separate instructions.

**SNP calling**

    GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../f000-reference.fa -I f050-dna_100_high_pe_fixRG_sorted_noDup_realigned_recalibrated.bam -glm SNP \  
    -o f060-dna_100_high_pe_snps.vcf

**INDEL calling**

    GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../f000-reference.fa -I f050-dna_100_high_pe_fixRG_sorted_noDup_realigned_recalibrated.bam -glm INDEL \  
    -o f060-dna_100_high_pe_indels.vcf


7. Introduce filters in the VCF file
--------------------------------------------------------------------------------

Example: filter SNPs with low confidence calling (QD < 12.0) and flag them as "LowConf".

    GenomeAnalysisTK.jar -T VariantFiltration -R ../f000-reference.fa -V f060-dna_100_high_pe_snps.vcf --filterExpression "QD < 12.0" --filterName "LowConf" \  
    -o f070-dna_100_high_pe_snps_filtered.vcf

The command ``--filterExpression`` will read the INFO field and check whether variants satisfy the requirement. If a variant does not satisfy your filter expression, the field FILTER will be filled with the indicated ``--filterName``. These commands can be called several times indicating different filtering expression (i.e: --filterName One --filterExpression "X < 1" --filterName Two --filterExpression "X > 2").

**QUESTION:** How many "LowConf" variants have we obtained?

    grep LowConf f070-dna_100_high_pe_snps_filtered.vcf | wc -l

