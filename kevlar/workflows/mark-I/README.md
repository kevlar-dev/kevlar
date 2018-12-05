# *De novo* variant discovery workflow: Mark I

Make a copy of `config.json` and replace the existing filenames with the full paths of your input filenames.
An arbitrary number of Fasta/Fastq files is supported for each sample.
Paired-end reads are not required, but pairing information will be retained if provided for the case/proband sample (pairing info is always ignored for control samples).

The other parameters are reasonable for error-corrected human samples sequenced at ≈30x coverage.
Without error correction much more memory will be required, or an alternative workflow will be required.
kevlar will be accurate even when the *k*-mer counting false-positive rate (FPR) is fairly high in the case sample, as this is corrected in subsequent steps.
Sensitivity will be lost if the *k*-mer counting FPR is too high in the control samples, however.

```
# Do a dry run to make sure there are no warnings or errors.
snakemake --configfile myconfig.json --snakefile kevlar/workflows/mark-I/Snakefile --cores 64 --directory /path/to/desired/work/dir/ --printshellcmds --dryrun calls

# If everything looks good, go ahead an execute the procedure.
snakemake --configfile myconfig.json --snakefile kevlar/workflows/mark-I/Snakefile --cores 64 --directory /path/to/desired/work/dir/ --printshellcmds calls
```
