# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]
### Fixed
- Abundance list reported by `kevlar filter` now correctly show re-computed proband k-mer abundances, not pre-filtering abundances as before (see #111).
- The `kevlar localize` and `kevlar call` procedures now handle multiple assembled contigs and multiple reference matches (see #124 and #126).

### Added
- New abundance screen now a part of `kevlar novel`. If any k-mer in a read is below some abundance threshold, the entire read is discarded (see #106).
- Better error reporting and handling of various issues with assembly, localization, and alignment (see #113, #114).
- Preliminary support for VCF output.

### Changed
- The `kevlar filter` procedure now handles both contamination and reference matches under a single mask interface (see #103).
- Explicitly dropped support for Python 2.7. Now supports only Python >=3.5 (see #125).

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
