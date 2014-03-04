% [NGS data analysis course](http://ngscourse.github.io/)
% __Title of the practical__
% _(updated 23-02-2014)_

<!-- COMMON LINKS HERE -->

http://samstat.sourceforge.net/#install

Exercise
================================================================================

Create a new directory for the exercise and copy the raw data to it: 

<!-- clean directory
    rm -r data
-->

    mkdir data
	cd data
    cp Path_To_Your/COURSE_EXAMPLE_DATA/reference_genome/f000_genome_sequence.fa .     #copy the reference genome
	cp Path_To_Your/COURSE_EXAMPLE_DATA/reference_genome/f005_genome_annotation.gtf .  #copy the annotation for the genome
	cp Path_To_Your/COURSE_EXAMPLE_DATA/rna_seq/* .  # copy the reads

<!--
	cp ../../../../COURSE_EXAMPLE_DATA/reference_genome/f000_genome_sequence.fa .
    cp ../../../../COURSE_EXAMPLE_DATA/reference_genome/f005_genome_annotation.gtf .
    cp ../../../../COURSE_EXAMPLE_DATA/rna_seq/* .
-->

The final __dot__ in the lines above means "here", in the working directory. It is neccesary to write it.




<!--

NOTE: they are the same 

In dwgsim:
    -d INT: __inner distance__ between the two ends [500]

In tophat2:
   -r/--mate-inner-dist           <int>       [ default: 50 ]

tophat2 -r 300 -o cdna_top_spa_T1 Homo_sapiens.GRCh37.74.dna.chromosome.21 simulacion/f010_T1_read1.fastq simulacion/f010_T1_read2.fastq 
tophat2 -r 300 -o cdna_top_spa_T2 Homo_sapiens.GRCh37.74.dna.chromosome.21 simulacion/f010_T2_read1.fastq simulacion/f010_T2_read2.fastq 
tophat2 -r 300 -o cdna_top_spa_T3 Homo_sapiens.GRCh37.74.dna.chromosome.21 simulacion/f010_T3_read1.fastq simulacion/f010_T3_read2.fastq 

tophat2 -r 300 -o cdna_top_spa_C1 Homo_sapiens.GRCh37.74.dna.chromosome.21 simulacion/f010_C1_read1.fastq simulacion/f010_C1_read2.fastq 
tophat2 -r 300 -o cdna_top_spa_C2 Homo_sapiens.GRCh37.74.dna.chromosome.21 simulacion/f010_C2_read1.fastq simulacion/f010_C2_read2.fastq 
tophat2 -r 300 -o cdna_top_spa_C3 Homo_sapiens.GRCh37.74.dna.chromosome.21 simulacion/f010_C3_read1.fastq simulacion/f010_C3_read2.fastq 

samtools view cdna_top_spa_T1/accepted_hits.bam > cdnaT1.sam
samtools view cdna_top_spa_T2/accepted_hits.bam > cdnaT2.sam
samtools view cdna_top_spa_T3/accepted_hits.bam > cdnaT3.sam

samtools view cdna_top_spa_C1/accepted_hits.bam > cdnaC1.sam
samtools view cdna_top_spa_C2/accepted_hits.bam > cdnaC2.sam
samtools view cdna_top_spa_C3/accepted_hits.bam > cdnaC3.sam

cp cdna_top_spa_T1/accepted_hits.bam > cdnaT1.bam
cp cdna_top_spa_T2/accepted_hits.bam > cdnaT2.bam
cp cdna_top_spa_T3/accepted_hits.bam > cdnaT3.bam

cp cdna_top_spa_C1/accepted_hits.bam > cdnaC1.bam
cp cdna_top_spa_C2/accepted_hits.bam > cdnaC2.bam
cp cdna_top_spa_C3/accepted_hits.bam > cdnaC3.bam

wc -l *.sam

################################################################################

## cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> 

cuffdiff -o cuffdiff_simple_sam chr21.gtf cdnaT1.sam cdnaC2.sam
cuffdiff -o cuffdiff_simple_bam chr21.gtf cdnaT1.bam cdnaC2.bam   ##DOES NOT WORK WITH BAM !!!

cuffdiff -o cuffdiff_all_sam chr21.gtf   cdnaT1.sam,cdnaT2.sam,cdnaT3.sam   cdnaC1.sam,cdnaC2.sam,cdnaC3.sam
cuffdiff -o cuffdiff_all_bam chr21.gtf   cdnaT1.bam,cdnaT2.bam,cdnaT3.bam   cdnaC1.bam,cdnaC2.bam,cdnaC3.bam  ##DOES NOT WORK WITH BAM !!!

-->


Map against the reference genome using bowtie2
--------------------------------------------------------------------------------

Fist we need to __build an index__ for bowtie:

    bowtie2-build f000_genome_sequence.fa f001_bowtie_index


And now we can run the __alignments__ for the __paired end__ files:

	tophat2 -r 300 -o f020_C1_tophat_out   f001_bowtie_index   f010_C1_read1.fastq f010_C1_read2.fastq 
	tophat2 -r 300 -o f020_C2_tophat_out   f001_bowtie_index   f010_C2_read1.fastq f010_C2_read2.fastq 
	tophat2 -r 300 -o f020_C3_tophat_out   f001_bowtie_index   f010_C3_read1.fastq f010_C3_read2.fastq 

	tophat2 -r 300 -o f020_T1_tophat_out   f001_bowtie_index   f010_T1_read1.fastq f010_T1_read2.fastq 
	tophat2 -r 300 -o f020_T2_tophat_out   f001_bowtie_index   f010_T2_read1.fastq f010_T2_read2.fastq 
	tophat2 -r 300 -o f020_T3_tophat_out   f001_bowtie_index   f010_T3_read1.fastq f010_T3_read2.fastq 


<!-- Use [samtools] Convert to SAM (text file) in order to explore them: -->

<!--     samtools view f020_tophat_out_case/accepted_hits.bam > f025_accepted_case.sam -->
<!--     samtools view f020_tophat_out_cont/accepted_hits.bam > f025_accepted_cont.sam -->

    samtools view f020_C1_tophat_out/accepted_hits.bam > f030_C1.sam
    samtools view f020_C2_tophat_out/accepted_hits.bam > f030_C2.sam
    samtools view f020_C3_tophat_out/accepted_hits.bam > f030_C3.sam

    samtools view f020_T1_tophat_out/accepted_hits.bam > f030_T1.sam
    samtools view f020_T2_tophat_out/accepted_hits.bam > f030_T2.sam
    samtools view f020_T3_tophat_out/accepted_hits.bam > f030_T3.sam

    wc -l *.sam



<!-- Explore the first 2 lines of the SAM file -->

<!--     head -n 2 f025_accepted_case.sam -->
	

<!-- Find out the original reads in the fastq files: -->

<!--     grep "_436_38_1_0_0_0_1:0:0_1:0:0_0" f010* -->

<!-- See the sequence of the first read: -->

<!--     grep -A 2 "_436_38_1_0_0_0_1:0:0_1:0:0_0" f010* -->


<!-- and compare them with the corresponding alignments in the SAM file: -->

<!--     grep "_436_38_1_0_0_0_1:0:0_1:0:0_0" f025_accepted_case.sam -->


Differential expression using Cuffdiff
--------------------------------------------------------------------------------

__Using the information in the GTF file:__

Compare 2 samples:

    cuffdiff -o f040_dif_exp_two f005_genome_annotation.gtf   f030_T1.sam   f030_C1.sam


Compare several samples:

    cuffdiff -o f040_dif_exp_several f005_genome_annotation.gtf   f030_T1.sam,f030_T2.sam,f030_T3.sam   f030_C1.sam,f030_C2.sam,f030_C2.sam

Explore the results.

