# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Added
- New `kevlar gentrio` command for a more realistic similation of trios for testing and evaluation (#171).
- New filter for `kevlar alac` for discarding partitions with a small number of interesting k-mers (#189).

### Changed
- Replaced `pep8` with `pycodestyle` for enforcing code style in development (see #167).
- The `--refr` argument of the `kevlar dump` command is now optional, and when no reference is explicitly specified `kevlar dump` acts primarily as a BAM to Fastq converter (see #170).
- Split the functionality of the `count` subcommand: simple single-sample k-mer counting was kept in `count` with a much simplified interface, while the memory efficient multi-sample "masked counting" strategy was split out to a new subcommand `effcount` (see #185).
- Replaced `kevlar reaugment` with a more generalizable `kevlar augment` subcommand (see #188).

### Fixed
- Incorrect file names in the quick start documentation page (see 9f6bec06d4).
- The `kevlar alac` procedure now accepts a stream of read partitions (instead of a stream of reads) at the Python API level, and correctly handles a single partition-labeled sequence file at the CLI level (see #165).
- CIGARs that begin with I blocks (alternate allele contig is longer than reference locus) are now handled properly (see #191).
- Bug with how `kevlar alac` handles "no reference match" scenarios resolved (see #192).

## [0.3.0] - 2017-11-03

This release includes many new features, some refactoring of the core codebase, and the first end-to-end analysis workflow implemented in a single command.
Details are included below.

### Fixed
- Abundances reported by `kevlar filter` now correctly show re-computed proband k-mer abundances, not pre-filtering abundances (see #111).
- The `kevlar localize` and `kevlar call` procedures now handle multiple assembled contigs, calling variants from the best reference match for each contig (see #124, #126, and #147).

### Added
- New abundance screen now a part of `kevlar novel`. If any k-mer in a read is below some abundance threshold, the entire read is discarded (see #106).
- Better error reporting and handling of various issues with assembly, localization, and alignment (see #113, #114).
- Support for VCF output (see #130 and #144), including "windows" with all k-mers containing the reference allele (RW) and alternate allele (VW) to facilitate distinguishing inherited mutations from novel mutations (see #144 and #152).
- New subcommands
    - `alac`: assembles, localizes, aligns, and calls variants on a single partition basis
    - `simplex`: invokes the entire simplex analysis workflow

### Changed
- The `kevlar filter` procedure now handles both contamination and reference matches under a single "mask" interface (see #103).
- Explicitly dropped support for Python 2.7. Now supports only Python >=3.5 (see #125).
- Main methods for each core subcommand are now implemented as minimal wrappers around generator functions, to facilitate composing different steps of the workflow or invoking them from third-party Python code (see #95, #126, #133, #148, #149, #150, #159, #161).
- The home-grown greedy assembly implementation has been replaced by calls to the `fermi-lite` library, which is now bundled with kevlar (see #156).
- The default behavior of `kevlar partition` is now to output a single stream of reads.
  Writing each partition to a distinct file is still supported with the ``--split`` option.

### Removed
- The `kevlar collect` command and associated tests. Its functionality has now been fully distributed to other subcommands.
    - Read filtering to `kevlar filter`
    - Junction count contig assembly to `kevlar filter` as an optional mode

## [0.2.0] - 2017-07-21
### Added
- New subcommands
    - `partition`: group reads by shared interesting k-mers
    - `localize`: determine an assembled contig's location in the reference genome
    - `call`: align assembled contigs to reference and call variant
- Documentation suite in `docs/`, hosted at https://kevlar.readthedocs.io
- New third-party dependency `ksw2` for computing alignments. Wrapped with Cython, which is a new development-time dependency (but not install or run time).
- The `pandas` package is now a dependency, and `pysam` and `networkx` are now hard dependencies (rather than conditional).

### Fixed
- Bug with assembly when the order of a read pair was swapped and they had the opposite orientation (see #85).

## [0.1.0] - 2017-05-13
### Added
- Command-line interface with 8 subcommands
    - `dump`: discard reads that match reference completely
    - `count`: compute k-mer abundances for all samples
    - `novel`: identify "interesting" (potentially novel) k-mers
    - `filter`: re-compute k-mer abundances, discard false positives and contamination
    - `assemble`: assemble reads for a single variant
    - `collect`: collect and filter (legacy)
    - `mutate`: simulate variants on a genome
    - `reaugment`: re-attach interesting k-mer annotations to reads
- Extensive test suite
- Continuous integration configuration
