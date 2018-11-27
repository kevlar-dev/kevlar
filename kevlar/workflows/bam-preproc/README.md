# Pre-processing procedure for mapped reads

Predicting *de novo* variants with kevlar requires one or more files per individual in Fasta or Fastq format.
Using error-corrected reads can [drastically reduce the amount of memory](https://standage.github.io/information-content-versus-data-volume-and-k-mer-counting-accuracy.html) required to count *k*-mers and detect variants.
However, data from human studies are often distributed as pre-aligned reads in BAM format.
The following procedure will convert the BAM files back into Fastq format, apply basic quality control, and correct sequencing errors in the reads.


## Software and hardware configuration

This procedure requires the following software environment.

- [Snakemake](https://bitbucket.org/snakemake/snakemake/), version 5.0 or newer
- [SAMtools and bgzip](https://github.com/samtools/samtools), version 1.9 or newer
- [fastp](https://github.com/OpenGene/fastp), tested with version 0.19.5
- [lighter](https://github.com/mourisl/Lighter), version 1.1.2 or newer

The file `config.json` contains the full paths of the BAM files.
Edit this as needed.
Note that only a single case sample is supported, but an arbitrary number of controls is supported.

**NOTE**: The number of required threads is hard-coded for some rules, and because these rules are usually grouped (run simultaneously through a "pipe" file) Snakemake's normal behavior of scaling threads down doesn't work as usual.
Some steps run 3 commands at a time with 24 threads per command, while others run 2 commands at a time with 36 threads per command (72 total threads).
I'd like to handle this more gracefully in the future (see [this thread](https://bitbucket.org/snakemake/snakemake/issues/1039/split-number-of-threads-for-a-rule)), but in the mean time it may be necessary to adjust the `threads` directive for the rules in the Snakefile to make better use of the resources available on your machine(s).


## Executing the procedure

Do a dry run to make sure there are no warnings or errors.

```
snakemake --configfile config.json --cores 72 --directory /path/to/desired/work/dir/ --printshellcmds --dryrun reads
```

If everything looks good, go ahead and execute the procedure.

```
snakemake --configfile config.json --cores 72 --directory /path/to/desired/work/dir/ --printshellcmds reads
```

**NOTE**: This procedure can be invoked from any directory on your system, assuming the following conditions are met.

- the file paths in the `config.json` are correct absolute paths
- the correct relative or absolute path to the `config.json` file is provided via the `--configfile` flag
- the `--snakefile` flag is appended to the command to specify the correct relative or absolute path to the `Snakefile` file
