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

     curl -L https://s3-us-west-1.amazonaws.com/noble-trios/helium-mother-reads.fq.gz -o mother.fq.gz
     curl -L https://s3-us-west-1.amazonaws.com/noble-trios/helium-father-reads.fq.gz -o father.fq.gz
     curl -L https://s3-us-west-1.amazonaws.com/noble-trios/helium-proband-reads.fq.gz -o proband.fq.gz
     curl -L https://s3-us-west-1.amazonaws.com/noble-trios/helium-refr.fa.gz -o refr.fa.gz
     bwa index refr.fa.gz

     curl -L https://s3-us-west-1.amazonaws.com/noble-trios/helium-config.json | sed "s:/home/user/Desktop:$(pwd):g" > helium-config.json
     snakemake --snakefile kevlar/workflows/mark-I/Snakefile --configfile helium-config.json --cores 4 --directory workdir -p calls
