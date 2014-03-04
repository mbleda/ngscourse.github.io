echo "export PATH=$PATH:/home/Desktop/hpg-variant-1.0" >> ~/.bashrc
source ~/.bashrc
hpg-var-effect -v CHB.exon.2010_03.sites.vcf --outdir hpg-variant_results --config /home/parallels/course/variant_annotation/hpg-variant-1.0/
hpg-var-vcf split -v CEU.exon.2010_03.genotypes.vcf --criterion chromosome --config /home/parallels/course/variant_annotation/hpg-variant-1.0/
hpg-var-vcf filter -v CEU.exon.2010_03.genotypes.vcf --coverage 9000 --save-rejected --config /home/parallels/course/variant_annotation/hpg-variant-1.0/
hpg-var-vcf stats -v CEU.exon.2010_03.genotypes.vcf -config /home/parallels/course/variant_annotation/hpg-variant-1.0/
hpg-var-vcf stats -v CEU.exon.2010_03.genotypes.vcf --samples -config /home/parallels/course/variant_annotation/hpg-variant-1.0/
hpg-var-gwas tdt -v 4Kvariants_147samples.vcf -p 4Kvariants_147samples.ped -config /home/parallels/course/variant_annotation/hpg-variant-1.0/
hpg-var-gwas assoc --chisq -v 4Kvariants_147samples.vcf -p 4Kvariants_147samples.ped -config /home/parallels/course/variant_annotation/hpg-variant-1.0/
hpg-var-gwas assoc --fisher -v 4Kvariants_147samples.vcf -p 4Kvariants_147samples.ped -config /home/parallels/course/variant_annotation/hpg-variant-1.0/
