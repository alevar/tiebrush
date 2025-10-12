# Copilot instructions for this repo

## What this project is
- C++11 CLI tools for large-scale read aggregation/summary on BAM/SAM/CRAM:
  - tiebrush: collapses duplicate alignments across many samples and writes a composite BAM with custom tags.
  - tiecov: derives coverage, junction BED, and optional BigWig from tiebrush (or plain) BAMs.
  - tiewrap.py: convenience batch runner; sashimi.py: plotting helper (Python, optional).

## Big picture and data flow
- Input: coordinate-sorted BAMs (all must share the same references/order).
- tiebrush streams records position-sorted across inputs, groups “same alignment” by strategy (CIGAR/full/clip/exons), and emits one record with:
  - YC: double multiplicity (can be fractional if --store-frac used on secondarys)
  - YX: int number of samples carrying the alignment
  - YD: int max upstream bundle distance per sample/strand (omitted if 0)
  - Sample provenance is stored as @CO SAMPLE:<abs path> lines in the BAM header.
- tiecov reads alignments, accumulates:
  - coverage (BedGraph or BigWig with -W),
  - splice junctions (BED),
  - sample-count heatmap (BedGraph; normalizes counts across samples).

## Build and dependencies
- Build system: CMake (C++11), external projects vendored at build time:
  - htslib (static libhts.a), libBigWig (for BigWig output).
- Options in CMakeLists.txt:
  - -DTIEBRUSH_BUILD_LIBBIGWIG=ON|OFF, -DTIEBRUSH_NATIVE_BUILD=ON, -DTIEBRUSH_STATIC_BUILD=ON.
- Required system libs: zlib, liblzma, libbz2, pthreads. Tests require samtools in PATH.
- Typical build: cmake -DCMAKE_BUILD_TYPE=Release . && make -j; executables: tiebrush, tiecov; install also tiewrap.py.

## Tests and workflows
- run_tests.sh uses example/ data and validates output with samtools view | diff. Fails if tiebrush/tiecov not built or samtools missing.
- VS Code tasks: Build, Clean & Build, Build testsam only, Run Tests.
- Example CLI:
  - tiebrush -o out.bam samples/s1.bam samples/s2.bam
  - tiecov -c out.cov.bedgraph -j out.juncs.bed -s out.sample.bedgraph out.bam

## Key internal APIs and patterns
- Core headers: src/GSam.h, src/GSam.cpp provide thin wrappers over htslib:
  - GSamReader/Writer stream records; GSamRecord exposes fields (start/end/exons are 1-based) and tag helpers:
    - add_int_tag("YX", int64), add_double_tag("YC", double), remove_tag("YD"), tag_int/tag_float/tag_str.
  - spliceStrand() resolves XS or minimap2 ts, adjusted by reverse flag.
- Grouping/merging logic (src/tiebrush.cpp):
  - Strategies: tMrgStratCIGAR, tMrgStratFull (CIGAR+MD), tMrgStratClip (trimmed soft clips), tMrgStratExon.
  - SPData aggregates duplicates, tracks per-sample bitset; rspacing computes YD via per-sample, per-strand segment lists.
  - TInputFiles (src/tmerge.*) validates headers, merges @CO SAMPLE lines, adds @PG TieBrush, and streams next() across files.
- Coverage/junctions (src/tiecov.cpp):
  - addCov/addJunction walk CIGAR; flushCoverage writes run-length segments; BigWig uses libBigWig bwAdd/Append.
  - Heatmap uses YX, then discretize+normalize to a fixed visual range.

## Repo conventions that matter
- Coordinates: GSamRecord.start/end/exons are 1-based; htslib core.pos is 0-based internally—be consistent when writing outputs.
- YC is double (fractional support), YX/YD integers. Remove YD if zero.
- Headers carry sample provenance as @CO SAMPLE:<abs path>. Use commons.h helpers to parse/load.
- CLI parsing via gclib/GArgs; mirror existing short/long forms and verbose flag patterns.
- Prefer gclib containers (GStr, GVec, GList) and keep streaming algorithms O(1) memory w.r.t. genome span.

## When extending
- Adding a new output tag: use GSamRecord::add_int_tag/add_double_tag; read via tag_int/tag_float. Keep types stable across tools.
- New CLI options: wire with GArgs in processOptions(), update USAGE strings, and keep mutually-exclusive strategy checks.
- Cross-tool data: ensure tiecov remains robust with plain BAMs (missing YC/YX/YD) by using defaults as in current code.
- Tests: add small BAMs under example/ and assert with diff_check in run_tests.sh style.

Files to study first
- CMakeLists.txt, src/tiebrush.cpp, src/tiecov.cpp, src/tmerge.{h,cpp}, src/GSam.{h,cpp}, src/commons.h, example/README.md.

Feedback
- If any of the above is unclear (e.g., header SAMPLE lines, YD semantics, or BigWig expectations), tell me what to refine and I’ll update this guide.
