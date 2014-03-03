rm -r data
mkdir data
cd data
cp ../../../../COURSE_EXAMPLE_DATA/f010_mirbase_mature.fa .
cd my_visual_data_dir  
cd ~/ngscourse.github.io/COURSE_EXAMPLE_DATA/visualization/example_1

ll
bwa index -a bwtsw f000-reference.fa
samtools faidx f000-reference.fa
java -jar CreateSequenceDictionary.jar REFERENCE=f000-reference.fa OUTPUT=f000-reference.dict
java -jar MarkDuplicates.jar INPUT=f000-paired_end.bam OUTPUT=f010-paired_end_noDup.bam METRICS_FILE=metrics.txt
java -jar BuildBamIndex.jar INPUT=f010-paired_end_noDup.bam

#java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R f000-reference.fa -I f010-paired_end_noDup.bam -known gold_indels.vcf -o forIndelRealigner.intervals
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R f000-reference.fa -I f010-paired_end_noDup.bam -o f020-indelRealigner.intervals
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R f000-reference.fa -I f010-paired_end_noDup.bam -targetIntervals f020-indelRealigner.intervals -o f030-paired_and_noDup_realigned.bam 
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R reference.fa -I realigned_reads.bam -L 20 -knownSites dbsnp.vcf -knownSites gold_indels.vcf -o recal_data.table 
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R reference.fa -I realigned_reads.bam -L 20 -knownSites dbsnp.vcf -knownSites gold_indels.vcf -BQSR recal_data.table -o post_recal_data.table 
java -jar GenomeAnalysisTK.jar -T PrintReads -R reference.fa -I realigned_reads.bam -L 20 -BQSR recal_data.table -o recal_reads.bam
java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R f000_reference.fa -I 7_realigned_aligned.bam -glm SNP -o 8_snp_variants.vcf
java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R 0_reference.fa -I 7_realigned_aligned.bam -glm INDEL -o 8_indel_variants.vcf
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R 0_reference.fa -V 8_snp_variants.vcf --filterExpression "QD < 12.0" --filterName "LowConf" -o 9_snp_filtered.vcf
