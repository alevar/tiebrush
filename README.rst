TieBrush and TieCov - efficient methods for aggregating and summarizing aligned sequences across large datasets
===============================================================================================================

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: https://opensource.org/licenses/MIT
    :alt: MIT License

.. contents::
   :local:
   :depth: 2

Introduction
^^^^^^^^^^^^

TieBrush is a simple yet efficient method for merging redundant information from multiple alignment files, 
designed to enable rapid manipulation of extremely large sequencing datasets. The method is specifically 
designed to optimize investigations of RNA, whole-genome, exome and other types of sequencing experiments. 
TieBrush preserves much of the original information in a greatly condensed representation as a BAM file, 
which allows manipulation and extraction of dataset and subset-specific statistics using tools within 
the package, and which is also compatible with other common utilities.

This utility aims to merge/collapse "duplicate" read alignments (same location with the same CIGAR string),
across multiple sequencing samples (multiple input BAM files), adding custom SAM tags in order to keep
track of the "alignment multiplicity" count (how many times the same alignment is seen across all
input data) and "sample count" (how many samples show that same alignment).
The initial goal is to generate this composite BAM file which multiplexes read alignments
from many sequencing samples, painting a comprehensive "background" picture of read alignments
with their counts across many samples.

    Ales Varabyou, Geo Pertea, Christopher Pockrandt, Mihaela Pertea, TieBrush: an efficient method for aggregating and summarizing mapped reads across large datasets, Bioinformatics, 2021;, btab342, https://doi.org/10.1093/bioinformatics/btab342

Sahimi plot is largely based on the implementation from the MISO package. please cite both the TieBrush publication as well as the original MISO paper:

    Katz, Yarden, Eric T. Wang, Jacob Silterra, Schraga Schwartz, Bang Wong, Helga Thorvaldsdóttir, James T. Robinson, Jill P. Mesirov, Edoardo M. Airoldi, and Christopher B. Burge. "Quantitative visualization of alternative exon expression from RNA-seq data." Bioinformatics 31, no. 14 (2015): 2400-2402.

Installation
^^^^^^^^^^^^

Building from source
""""""""""""""""""""

If you want to build it from source, we recommend cloning the git repository as shown below to ensure
fixed releases of any dependencies are fetched and compiled with the software.

::

    $ git clone --recursive https://github.com/alevar/tiebrush.git
    $ cd tiebrush/
    $ mkdir build && cd build
    $ cmake -DCMAKE_BUILD_TYPE=Release ..
    $ make -j4
    $ make test #optional 
    $ make install

If you are using a very old version of Git (< 1.6.5) the flag ``--recursive`` does not exist.
In this case you need to clone the submodule separately  (``git submodule update --init --recursive``).

For a fully static build `-DTIEBRUSH_STATIC_BUILD=1` needs to be added to the list of arguments
in the cmake command.

By default `make install` will require administrative privileges. To specify a custom
installation path `-DCMAKE_INSTALL_PREFIX=<custom/installation/path>` needs to be added to
the list of arguments in the `cmake` command.

**Requirements**

Operating System
  GNU/Linux, Mac

Compiler
  GCC ≥ 4.8, LLVM/Clang ≥ 3.8

Build system
  CMake ≥ 2.8

Language support
  C++11

Methods
^^^^^^^

TieBrush
""""""""

Summarize and filter read alignments from multiple sequencing samples (taken as sorted BAM files).
This utility aims to merge/collapse "duplicate" read alignments (same location with the same
CIGAR string), across multiple sequencing samples (multiple input BAM files), adding custom SAM tags
in order to keep track of the "alignment multiplicity" count (how many times the same alignment is
seen across all input data) and "sample count" (how many samples show that same alignment).

The goal is to generate this composite BAM file which multiplexes read alignments from many sequencing
samples, painting a comprehensive "background" picture of read alignments with their counts across
many samples.

  tiebrush  [-h] -o OUTPUT [-L|-P|-E] [-S] [-M] [-N max_NH_value] [-Q min_mapping_quality] [-F FLAGS] ...

  Input arguments:

  ...        Input can be provided as a space-delimited list of filenames or as a text file containing a list of filenames, one per line

  Required arguments:

  -o        File for BAM output

  Optional arguments:

  -h, --help        Show this help message and exit
  --version         Show the program version end exit
  -L, --full        If enabled, only reads with the same CIGAR and MD strings will be grouped and collapsed. By default, TieBrush will consider the CIGAR string only when grouping reads
  -P, --clip        If enabled, reads will be grouped by clipped CIGAR string. In this mode 5S10M5S and 3S10M3S CIGAR strings will be grouped if the coordinates of the matching substring (10M) are the same between reads
  -E, --exon        If enabled, reads will be grouped if their exon boundaries are the same. This option discards any structural variants contained in mapped substrings of the read and only considers start and end coordinates of each non-splicing segment of the CIGAR string
  -S, --keep-supp   If enabled, supplementary alignments will be included in the collapsed groups of reads. By default, TieBrush removes any mappings not listed as primary (0x100). Note, that if enabled, each supplementary mapping will count as a separate read
  -M, --keep-unmap  If enabled, unmapped reads will be retained (uncollapsed) in the output. By default, TieBrush removes any unmapped reads
  -N                Maximum NH score (if available) to include.
  -Q                Minimum mapping quality to include.
  -F                Bits in SAM flag to use in read comparison. Only reads that have specified flags will be merged together (default: 0)

Note that options -L, -P and -E are mutually exclusive. 


SAM tags implemented
--------------------
1. __YC__:i:N stores the number of alignments that were merged into this alignment record (multiplicity count)
2. __YX__:i:N stores the number of samples that have this alignment (sample count)
3. __YD__:i:N keeps track of the maximum number of contiguous bases preceding the start of the read alignment in the samples(s) that it belongs to. In other words, if the current alignment is part of an exon-overlapping bundle (strand specific!), this value holds the maximum distance from the beginning of the bundle to the start of this alignment, across all samples having this alignment. If the alignment is not in a bundle (i.e. it is preceded by a uncovered region as it is not overlapped by any another alignment with a lower start position), in all the individual samples where that alignment is present, then the YD value is 0 and the tag is omitted from the output file produced by TieBrush. That means that all the alignments lacking a YD tag in the TieBrush output start at the very beginning of an exon-overlapping bundle (i.e. are not overlapped by a preceding alignment with a lower start coordinate).

If either YC or YX tags are missing (i.e. GBamRecord::__tag_int__() call returns 0) then the alignment is unique (when YC is 0) or only one sample has it (if YX is 0). The actual count in these cases is obviously 1.

TieCov
""""""

The TieCov utility can take the output file produced by TieBrush and can generate the following auxiliary files:

1. a BedGraph file with the coverage data (see http://genome.ucsc.edu/goldenPath/help/bedgraph.html); this file can be converted to BigWig (using bedGraphToBigWig) or to TDF format (using igvtools) in order to be loaded in IGV as an additional coverage track
2. a junction BED file which can be loaded directly in IGV as an additional junction track (http://software.broadinstitute.org/software/igv/splice_junctions)
3. a heatmap BED that uses color intensity to represent the number of samples that contain each position.

  tiecov [-s out.sample.bed] [-c out.coverage.bedgraph] [-j out.junctions.bed] [-W] input
  
  Input arguments (required):
  
  input  alignment file in SAM/BAM/CRAM format
  
  Optional arguments (at least one of -s/-c/-j must be specified):
  
  -s    output BED file with an estimate of the number of samples which contain alignments for each interval.
  -c    output BedGraph (or BedWig with '-W') file with coverage for all mapped bases.
  -j    output BED file with coverage of all splice-junctions in the input file.
  -W    save coverage in BigWig format. Default output is in Bed format.

TieWrap
"""""""

TieWrap is a small utility script provided to make running TieBrush on large datasets a bit easier.
Unlike TieBrush, TieWrap can be launched with as many input files as needed and will automatically
divide them into batches processing and combining batches to produce a single representation at the end.
All standard TieBrush arguments can be passed over to TieWrap. Additionally size of individual batches
as well as the concurrency parameters can be set explicitely.

  tiewrap.py [-h] -o OUTPUT [-L|-P|-E] [-S] [-M] [-N MAX_NH] [-Q MIN_MAP_QUAL] [-F FLAGS] [-t THREADS] [-b BATCH_SIZE] ...

  Input arguments:

  ...       Input can be provided as a space-delimited list of filenames or as a textfile containing a list of filenames one per each line.

  Required arguments:

  -o, --output          File for BAM output.

  Optional arguments:

  -h, --help            show this help message and exit
  -L, --full            If enabled, only reads with the same CIGAR and MD strings will be grouped and collapsed. By default, TieBrush will consider the CIGAR string only when grouping reads.
  -P, --clip            If enabled, reads will be grouped by clipped CIGAR string. In this mode 5S10M5S and 3S10M3S cigar strings will be grouped if the coordinates of the matching substring (10M) are the same between reads.
  -E, --exon            If enabled, reads will be grouped if their exon boundaries are the same. This option discards any structural variants contained in mapped substrings of the read and only considers start and end coordinates of each non-splicing segment of the CIGAR string.
  -S, --keep-supp       If enabled, supplementary alignments will be included in the collapsed groups of reads. By default, TieBrush removes any mappings not listed as primary (0x100). Note, that if enabled, each supplementary mapping will count as a separate read.
  -M, --keep-unmap      If enabled, unmapped reads will be retained (uncollapsed) in the output. By default, TieBrush removes any unmapped reads.
  -N, --max-nh          Maximum NH score of the reads to retain.
  -Q, --min-map-qual    Minimum mapping quality of the reads to retain.
  -F, --flags           Bits in SAM flag to use in read comparison. Only reads that have specified flags will be merged together (default: 0)
  -t, --threads         Number of threads to use.
  -b, --batch-size      Number of input files to process in a batch on each thread.
