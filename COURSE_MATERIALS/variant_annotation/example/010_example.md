e% [NGS data analysis course](http://ngscourse.github.io/)
% __Variant annotation__
% _(updated 28-02-2014)_

<!-- COMMON LINKS HERE -->

[HPG-Variant]: http://wiki.opencb.org/projects/hpg/doku.php?id=variant:overview "HPG Variant"


Preliminaries
================================================================================


Software used in this practical:
--------------------------------

- [HPG-Variant][HPG-Variant] : a complete suite of tools to work with genomic variation data, from VCF tools to variant profiling or genomic statistics.


File formats explored:
----------------------

- [VCF]()
- more


SAM/BAM

- links


Data used in this practical
-------------------------------

- [mirbase_mature.fa](../../../COURSE_EXAMPLE_DATA/f010_mirbase_mature.fa): mature micro RNAs downloaded form mirbase

You can download them or copy them to your ``data`` directory for the practical

<!-- clean directory
    rm -r data
-->

    mkdir data
	cd data
	cp ../../../../COURSE_EXAMPLE_DATA/f010_mirbase_mature.fa .


\ 

Find all data files for the course here: [COURSE_EXAMPLE_DATA](../../../COURSE_EXAMPLE_DATA)



Exercise
================================================================================

The file **f000_raw_mirna.fastq** in the data folder contains reads form a microRNA sequencing experiment.

Use `wc` to count how many reads there are in the file (remember you have to divide by 4)

    wc -l f000_raw_mirna.fastq


Raw Data Preprocessing
--------------------------------------------------------------------------------

Use fastqc to explore the quality of the data:

    fastqc f000_raw_mirna.fastq


Find the _Overrepresented sequences_ in [mirbase][mirbase-search].


There are 2 known adapters used in this experiment: 

    CTGGGAAATCACCATAAACGTGAAATGTCTTTGGATTTGGGAATCTTATAAGTTCTGTATGAGACCACTCTAAAAA
    CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC

Use R-Bioconductor to compute their reverse, complementary and reverse complementary.

    library (Biostrings)
    myseq <- DNAString ("CTTTTTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGC")
    reverse (myseq)
    complement (myseq)
    reverseComplement (myseq)



Further work
================================================================================

Some ideas

