rm -r data
mkdir data
cd data
cp Path_To_Your/COURSE_EXAMPLE_DATA/reference_genome/f000_genome_sequence.fa .     #copy the reference genome
cp Path_To_Your/COURSE_EXAMPLE_DATA/reference_genome/f005_genome_annotation.gtf .  #copy the annotation for the genome
cp Path_To_Your/COURSE_EXAMPLE_DATA/rna_seq/* .  # copy the reads
cp ../../../../COURSE_EXAMPLE_DATA/reference_genome/f000_genome_sequence.fa .
cp ../../../../COURSE_EXAMPLE_DATA/reference_genome/f005_genome_annotation.gtf .
cp ../../../../COURSE_EXAMPLE_DATA/rna_seq/* .
bowtie2-build f000_genome_sequence.fa f001_bowtie_index

tophat2 -o f020_C1_tophat_out   f001_bowtie_index   f010_C1_read1.fastq,f010_C1_read2.fastq 

tophat2 --read-mismatches 2 --read-gap-length 2 -r 300 -o f020_C1_tophat_out   f001_bowtie_index   f010_C1_read1.fastq,f010_C1_read2.fastq 
tophat2 --read-mismatches 2 --read-gap-length 2 -r 300 -o f020_C1_tophat_out   f001_bowtie_index   f010_C1_read1.fastq f010_C1_read2.fastq 
tophat2 --read-mismatches 2 --read-gap-length 2 -r 300 -o f020_C2_tophat_out   f001_bowtie_index   f010_C2_read1.fastq  f010_C2_read2.fastq 
tophat2 --read-mismatches 2 --read-gap-length 2 -r 300 -o f020_C3_tophat_out   f001_bowtie_index   f010_C3_read1.fastq  f010_C3_read2.fastq 
tophat2 --read-mismatches 2 --read-gap-length 2 -r 300 -o f020_T1_tophat_out   f001_bowtie_index   f010_T1_read1.fastq  f010_T1_read2.fastq 
tophat2 --read-mismatches 2 --read-gap-length 2 -r 300 -o f020_T2_tophat_out   f001_bowtie_index   f010_T2_read1.fastq  f010_T2_read2.fastq 
tophat2 --read-mismatches 2 --read-gap-length 2 -r 300 -o f020_T3_tophat_out   f001_bowtie_index   f010_T3_read1.fastq  f010_T3_read2.fastq 
samtools view f020_C1_tophat_out/accepted_hits.bam > f025_accepted_C1.sam
samtools view f020_C2_tophat_out/accepted_hits.bam > f025_accepted_C2.sam
samtools view f020_C3_tophat_out/accepted_hits.bam > f025_accepted_C3.sam
samtools view f020_T1_tophat_out/accepted_hits.bam > f025_accepted_T1.sam
samtools view f020_T2_tophat_out/accepted_hits.bam > f025_accepted_T2.sam
samtools view f020_T3_tophat_out/accepted_hits.bam > f025_accepted_T3.sam
wc -l *.sam

cuffdiff -o f030_cuffdiff_out f005_genome_annotation.gtf f020_T3_tophat_out/accepted_hits.bam f020_C2_tophat_out/accepted_hits.bam
	
cuffdiff -o f030_cuffdiff_out_sev f005_genome_annotation.gtf f020_T1_tophat_out/accepted_hits.bam f020_T2_tophat_out/accepted_hits.bam f020_T3_tophat_out/accepted_hits.bam f020_C1_tophat_out/accepted_hits.bam f020_C2_tophat_out/accepted_hits.bam f020_C3_tophat_out/accepted_hits.bam
