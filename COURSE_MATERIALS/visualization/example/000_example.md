% [NGS data analysis course](http://ngscourse.github.io/)
% __Visualization of mapped reads__
% _(updated 04-06-2014)_

<!-- COMMON LINKS HERE -->

[IGV]: http://www.broadinstitute.org/igv/home "IGV"
[Samtools]: http://samtools.sourceforge.net/ "samtools"


Preliminaries
================================================================================

Software used in this practical:
--------------------------------

- [IGV] : The Integrative Genomics Viewer is a program for reading several types of indexed database information, including mapped reads and variant calls, and displaying them on a reference genome. It is invaluable as a tool for viewing and interpreting the "raw data" of many NGS data analysis pipelines.
- [samtools] : SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.


File formats explored:
----------------------

- [SAM](http://samtools.sourceforge.net/SAMv1.pdf)
- [BAM](http://www.broadinstitute.org/igv/bam)


Guided exercise: IGV introduction
================================================================================

Go to your course directory:

    cd <my_course_directory>

Copy the data for the tutorials:

    cp -r /mounts/course_shares/Open_Share/visualization .
    cd example_0
    ll

This BAM file contains reads only for chromosome 11.

Creating indexed files
--------------------------------------------------------------------------------

Use ``samtools`` to index the bam files:

    samtools index igv1.bam


Run IGV
--------------------------------------------------------------------------------

    igv

Downolad a referece genome
--------------------------------------------------------------------------------

By default, IGV loads Human hg18 assembly. However, we must work with the **same assembly used to mapped our reads**, in this case Human hg19.  
Genome assemblies for several species can be dowloaded using IGV:

- Go to ``Genomes`` --> ``Load Genome From Server...``  
Select **Human hg19**

    
Loading and browsing files
--------------------------------------------------------------------------------

- Go to ``File`` --> ``Load from file...``  
Select ``igv1.bam``

**Steps:**

1. Go to the location box and insert **chr11:996096-1036047** in the search box and hit ``Go``.
2. Zoom in until you find some SNPs - they might be in exons or introns.
3. Identify at least one example of a short insertion variant and deletion arround exon 4.
4. Load and look at the SNP track: ``File`` --> ``Load from server`` --> ``Annotations`` --> ``Variants and Repeats`` --> ``dbSNP``


