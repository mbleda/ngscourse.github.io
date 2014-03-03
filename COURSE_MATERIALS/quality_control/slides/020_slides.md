% Main Title of the Presentation
% [NGS Data Analysis Course](http://ngscourse.github.io/)
% (updated 22-02-2014)


FastQ Format
================================================================================

- Standard Format for NGS data
- Conversion can be done from _sff_, _fasta + qual_, ... 
- Extension of the Fasta format
- Text-based formats (easy to use!)
- If not compressed, it can be huge

\  

\centering

<http://en.wikipedia.org/wiki/FASTQ_format>


Quality measurements
================================================================================

Base-calling __error probabilities__ are reported by sequencers.

Usually in __Phred__ (quality) score.

Usually coded by ASCII characters


Phred score
------------

$$Q = -10 log_{10} P$$

$$P = 10^{\frac{-Q}{10}}$$

\ 

<http://en.wikipedia.org/wiki/Phred\_quality\_score#Definition>




NGS Data Preprocessing Steps
================================================================================

- File parsing: convert to __fastq__ format form __sff__, __fasta__ + __qual__ ...
- Split [multiplex](http://www.illumina.com/technology/multiplexing_sequencing_assay.ilmn "Multiplex Sequencing Assay") samples.

- Quality Control of the raw data.

- Filtering and trimming reads by quality.
- Adapter trimming

- Quality Control of the trimmed and filtered reads



Software
================================================================================


- __FastQC__:
    - quality control
    - some filtering ...

\ 

- __Cutadapt__: 
    - adapter trimming 
	- filter reads by length (short, long)
	- filter reads by quality




<!-- FastQC Images -->

Per Base Sequence Quality
================================================================================
![](images/per_base_quality)

Per Sequence Quality
================================================================================
![](images/per_sequence_quality)

Per Base Sequence Content
================================================================================
![](images/per_base_sequence_content)

Per Base GC Content
================================================================================
![](images/per_base_gc_content)

Per Sequence Nucleotide Content
================================================================================
![](images/per_sequence_gc_content)

Per Base N Content
===============================================================================
![](images/per_base_n_content)

Sequence Length Distribution
===============================================================================
![](images/sequence_length_distribution)

Duplicate Sequences Distribution
================================================================================
![](images/duplication_levels)

Overrepresented Kmers
================================================================================
![](images/kmer_profiles)
