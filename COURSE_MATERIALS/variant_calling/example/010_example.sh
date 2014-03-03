echo "alias GenomeAnalysisTK.jar='java -jar /home/Desktop/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar'" >> ~/.bashrc
echo "alias CreateSequenceDictionary.jar='java -jar /home/Desktop/picard-tools-1.108/picard-tools-1.108/CreateSequenceDictionary.jar'" >> ~/.bashrc
echo "alias AddOrReplaceReadGroups.jar='java -jar /home/Desktop/picard-tools-1.108/picard-tools-1.108/AddOrReplaceReadGroups.jar'" >> ~/.bashrc
echo "alias MarkDuplicates.jar='java -jar /home/Desktop/picard-tools-1.108/picard-tools-1.108/MarkDuplicates.jar'" >> ~/.bashrc
echo "alias BuildBamIndex.jar='java -jar /home/Desktop/picard-tools-1.108/picard-tools-1.108/BuildBamIndex.jar'" >> ~/.bashrc
cd <my_course_directory>
cp -r /mounts/course_shares/Open_Share/variant_calling .
ll
bwa index -a bwtsw f000-reference.fa
samtools faidx f000-reference.fa
CreateSequenceDictionary.jar REFERENCE=f000-reference.fa OUTPUT=f000-reference.dict
cd example1
AddOrReplaceReadGroups.jar I=f000-dna_100_high_pe.bam O=f010-dna_100_high_pe_fixRG.bam RGID=group1 RGLB=lib1 RGPL=illumina RGSM=sample1 RGPU=unit1
samtools sort f010-dna_100_high_pe_fixRG.bam f020-dna_100_high_pe_fixRG_sorted
samtools index f020-dna_100_high_pe_fixRG_sorted.bam
MarkDuplicates.jar INPUT=f020-dna_100_high_pe_fixRG_sorted.bam OUTPUT=f030-dna_100_high_pe_fixRG_sorted_noDup.bam METRICS_FILE=metrics.txt
BuildBamIndex.jar INPUT=f030-dna_100_high_pe_fixRG_sorted_noDup.bam
GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../f000-reference.fa -I f030-dna_100_high_pe_fixRG_sorted_noDup.bam -o f040-indelRealigner.intervals
GenomeAnalysisTK.jar -T IndelRealigner -R ../f000-reference.fa -I f030-dna_100_high_pe_fixRG_sorted_noDup.bam -targetIntervals f040-indelRealigner.intervals \
-o f040-dna_100_high_pe_fixRG_sorted_noDup_realigned.bam 
GenomeAnalysisTK.jar -T BaseRecalibrator -R ../f000-reference.fa -I f040-dna_100_high_pe_fixRG_sorted_noDup_realigned.bam -knownSites ../f000-dbSNP_chr21.vcf \  
-o f050-recalibration_data.table
GenomeAnalysisTK.jar -T PrintReads -R ../f000-reference.fa -I f040-dna_100_high_pe_fixRG_sorted_noDup_realigned.bam -BQSR f050-recalibration_data.table \  
-o f050-dna_100_high_pe_fixRG_sorted_noDup_realigned_recalibrated.bam
GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../f000-reference.fa -I f050-dna_100_high_pe_fixRG_sorted_noDup_realigned_recalibrated.bam -glm SNP \  
-o f060-dna_100_high_pe_snps.vcf
GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../f000-reference.fa -I f050-dna_100_high_pe_fixRG_sorted_noDup_realigned_recalibrated.bam -glm INDEL \  
-o f060-dna_100_high_pe_indels.vcf
GenomeAnalysisTK.jar -T VariantFiltration -R ../f000-reference.fa -V f060-dna_100_high_pe_snps.vcf --filterExpression "QD < 12.0" --filterName "LowConf" \  
-o f070-dna_100_high_pe_snps_filtered.vcf
grep LowConf f070-dna_100_high_pe_snps_filtered.vcf | wc -l
