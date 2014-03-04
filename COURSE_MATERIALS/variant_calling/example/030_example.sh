echo "alias MuTect.jar='java -jar /home/Desktop/muTect-1.1.4-bin/muTect-1.1.4.jar'" >> ~/.bashrc
source ~/.bashrc
cd example3
AddOrReplaceReadGroups.jar I=f000-normal.bam O=f010-normal_fixRG.bam RGID=group1 RGLB=lib1 RGPL=illumina RGSM=sample1 RGPU=unit1
AddOrReplaceReadGroups.jar I=f000-tumor.bam O=f010-tumor_fixRG.bam RGID=group2 RGLB=lib1 RGPL=illumina RGSM=sample1 RGPU=unit2
samtools sort f010-normal_fixRG.bam f020-normal_fixRG_sorted
samtools sort f010-tumor_fixRG.bam f020-tumor_fixRG_sorted
samtools index f020-normal_fixRG_sorted.bam
samtools index f020-tumor_fixRG_sorted.bam
MuTect.jar --analysis_type MuTect --reference_sequence ../f000-reference.fa --dbsnp ../f000-dbSNP_chr21.vcf --input_file:normal f020-normal_fixRG_sorted.bam \  
--input_file:tumor f020-tumor_fixRG_sorted.bam --out f030-call_stats.out --coverage_file f030-coverage.wig.txt
