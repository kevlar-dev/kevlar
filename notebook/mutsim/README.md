# Download and preprocess reference genome

```bash
curl ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
    | gunzip -c \
    | aux/prep-genome.py --outfile GRCh38.noambig.fa - 
aux/split.py --numfiles 16 --out GRCh38.noambig.fa GRCh38.noambig.fa
```

# Count k-mers in the reference genome

Create a counttable that occupies 24G of memory.

```python
import khmer
ct = khmer.Counttable(31, 6e9, 4)
ct.consume_seqfile('GRCh38.noambig.fa')
ct.save('GRCh38.ct')
```

# Simulate mutations and collect novel k-mer distributions

On a machine with 128G, you can run 5 jobs simultaneously that require 24G of memory each.
Adjust accordingly depending on available memory and processors.

```bash
parallel --gnu --jobs 5 \
    bin/mut-hist -i 1000000 -k {2} -m 16 -r 0.1 -s 2018 -t {3} -z 5 GRCh38.noambig.fa.{1} GRCh38.ct \
        '>' k{2}-{3}-{1}.txt '2>' k{2}-{3}-{1}.log \
    ::: {1..16} \
    ::: 23 27 31 35 39 43 47 51 55 59 \
    ::: snv del

for muttype in snv del
do
    for k in 23 27 31 35 39 43 47 51 55 59
    do
        aux/combine.py k${k}-${muttype}-{1..16}.txt > k${k}-${muttype}.txt
    done
done
```
