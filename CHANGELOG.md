# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

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
