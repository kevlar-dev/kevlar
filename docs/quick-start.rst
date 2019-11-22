Quick start
===========

If you have not already done so, install kevlar using :doc:`the following instructions <install>`, and Snakemake version 5.0 or greater.

The simplest way to execute kevlar's entire *de novo* variant discovery workflow is using the provided Snakemake workflow configuration.
Processing the example data set below should be able to run on a laptop in less than 5 minutes while consuming less than 200 Mb of RAM.
The results (``workdir/calls.scored.sorted.vcf.gz``) should include 5 variant calls: a 300 bp insertion and 4 single-nucleotide variants.

A :doc:`more detailed tutorial is available <tutorial>`, and a complete listing of all available configuration options for each kevlar command can be found in :doc:`the CLI documentation <cli>`, or by executing ``kevlar <subcommand> -h`` in the terminal.

----------

.. highlight:: none

.. code::

     # Download data
     curl -L https://osf.io/db82p/download -o mother.fq.gz
     curl -L https://osf.io/6vrnz/download -o father.fq.gz
     curl -L https://osf.io/wt5h8/download -o proband.fq.gz
     curl -L https://osf.io/35wgn/download -o refr.fa.gz
     bwa index refr.fa.gz

     # Download and format configuration file
     curl -L https://osf.io/86adm/download | sed "s:/home/user/Desktop:$(pwd):g" > helium-config.json

     # Invoke the workflow
     snakemake \
         --snakefile kevlar/workflows/mark-I/Snakefile \
         --configfile helium-config.json --cores 4 --directory workdir -p calls
