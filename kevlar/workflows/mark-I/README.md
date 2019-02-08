# *De novo* variant discovery workflow: Mark I

The complete simplex *de novo* variant discovery workflow can be executed manually step-by-step, or as a single coordinated Snakemake workflow.
To invoke the Snakemake workflow:

- make sure you have Snakemake 5.0 or greater installed
- copy `kevlar/workflows/mark-I/Snakemake` and `kevlar/workflows/mark-I/config.json` to your working directory
- modify the copy of `config.json` and replace the existing filenames with the full paths of your input files
- execute the workflow with the `snakemake` command, using the examples below as a guide

```
# Do a dry run to make sure there are no warnings or errors.
snakemake --configfile myconfig.json --snakefile kevlar/workflows/mark-I/Snakefile --cores 64 --directory /path/to/desired/work/dir/ --printshellcmds --dryrun calls

# If everything looks good, go ahead an execute the procedure.
snakemake --configfile myconfig.json --snakefile kevlar/workflows/mark-I/Snakefile --cores 64 --directory /path/to/desired/work/dir/ --printshellcmds calls
```

Some notes on workflow configuration:

- An arbitrary number of Fasta/Fastq files is supported for each sample.
- Paired-end reads are not required, and any pairing information will be ignored.
- Parameters are tuned for error-corrected whole genome shotgun sequencing of human samples at ≈30x coverage. Without error correction much more memory will be required, or an alternative workflow will be required.
- kevlar will be accurate even when the *k*-mer counting false-positive rate (FPR) is fairly high in the case sample, as this is corrected in subsequent steps. Sensitivity will be lost if the *k*-mer counting FPR is too high in the control samples, however.
- The mask should include your reference genome and any potential sources of contamination. For example, UniVec includes many sources of technical contamination (such as adapters) that often cause problems.
- Accuracy can often be improved by filtering preliminary calls. Use the `varfilter` setting to provide a BED file of intervals from which to filter out variant predictions (by default, `varfilter: null` will disable this step). It's common to filter out variant calls in segmental duplications or SSRs, for example, as these are usually problematic. Also, filtering out common variants (using dbSNP, for example) will successfully remove inherited variants that are erroneously classified as *de novo* due to low coverage in the donor parent.
