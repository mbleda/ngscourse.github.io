% RNA-Seq Data Analysis
% [NGS Data Analysis Course](http://ngscourse.github.io/)
% (updated 04-03-2014)


Links
================================================================================

http://cufflinks.cbcb.umd.edu/manual.html#cuffdiff_output

RPKM
================================================================================

FPKM

input
================================================================================

has to be SAM not BAM

Cuffdiff Input




Cuffdiff Output
================================================================================

<Cuffdiff Output>Cuffdiff Output

- isoforms
- genes
- cds
- 
- 

- isoforms: Transcripts
- genes: Gene
- cds: Coding sequence
- tss: Primary transcript

Cuffdiff Output II
================================================================================


FPKM tracking files
-------------------
FPKM of each transcript, primary transcript, and gene in each sample. 

Primary transcript and gene FPKMs are computed by summing the FPKMs of transcripts in each primary transcript group or gene group. ????

The results are output in FPKM tracking files in the format described here. There are four FPKM tracking files:


Count tracking files
-----------------------
estimate of the number of __fragments__ that originated from each transcript

primary transcript, and gene in each sample. Primary transcript and gene counts are computed by summing the counts of transcripts in each primary transcript group or gene group. The results are output in count tracking files in the format described here. There are four Count tracking files:



Differential expression tests
-----------------------------
This tab delimited file lists the results of differential expression testing between samples for spliced transcripts, primary transcripts, genes, and coding sequences. For each pair of samples x and y, four files are created
