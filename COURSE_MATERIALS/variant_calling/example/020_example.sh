cd <my_course_directory>/variant_calling
cd example2
AddOrReplaceReadGroups.jar I=f000-dna_100_high_se.bam O=f010-dna_100_high_se_fixRG.bam RGID=group2 RGLB=lib2 RGPL=illumina RGSM=sample2 RGPU=unit2
samtools sort f010-dna_100_high_se_fixRG.bam f020-dna_100_high_se_fixRG_sorted
samtools index f020-dna_100_high_se_fixRG_sorted.bam
MarkDuplicates.jar INPUT=f020-dna_100_high_se_fixRG_sorted.bam OUTPUT=f030-dna_100_high_se_fixRG_sorted_noDup.bam METRICS_FILE=metrics.txt
BuildBamIndex.jar INPUT=f030-dna_100_high_se_fixRG_sorted_noDup.bam
GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../f000-reference.fa -I f030-dna_100_high_se_fixRG_sorted_noDup.bam -o f040-indelRealigner.intervals
GenomeAnalysisTK.jar -T IndelRealigner -R ../f000-reference.fa -I f030-dna_100_high_se_fixRG_sorted_noDup.bam -targetIntervals f040-indelRealigner.intervals \
-o f040-dna_100_high_se_fixRG_sorted_noDup_realigned.bam 
GenomeAnalysisTK.jar -T BaseRecalibrator -R ../f000-reference.fa -I f040-dna_100_high_se_fixRG_sorted_noDup_realigned.bam -knownSites ../f000-dbSNP_chr21.vcf \  
-o f050-recalibration_data.table
GenomeAnalysisTK.jar -T PrintReads -R ../f000-reference.fa -I f040-dna_100_high_se_fixRG_sorted_noDup_realigned.bam -BQSR f050-recalibration_data.table \  
-o f050-dna_100_high_se_fixRG_sorted_noDup_realigned_recalibrated.bam
GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../f000-reference.fa -I f050-dna_100_high_se_fixRG_sorted_noDup_realigned_recalibrated.bam -glm SNP \  
-o f060-dna_100_high_se_snps.vcf
GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../f000-reference.fa -I f050-dna_100_high_se_fixRG_sorted_noDup_realigned_recalibrated.bam -glm INDEL \  
-o f060-dna_100_high_se_indels.vcf
