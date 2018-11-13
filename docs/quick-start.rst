Quick start
===========

If you have not already done so, install kevlar using :doc:`the following instructions <install>`.

This gives a crash course on running kevlar's simplex analysis workflow.
The ``kevlar simplex`` command should be able to run on a laptop in about 5 minutes while consuming less than 200 Mb of RAM for this demo data set.
The results should include 5 variant calls: a 300 bp insertion and 4 single-nucleotide variants.

A :doc:`more detailed tutorial is available <tutorial>`, and a complete listing of all available configuration options for each script can be found in :doc:`the CLI documentation <cli>`, or by executing ``kevlar <subcommand> -h`` in the terminal.

----------

.. highlight:: none

.. code::

     curl -L https://s3-us-west-1.amazonaws.com/noble-trios/helium-mother-reads.fq.gz -o mother.fq.gz
     curl -L https://s3-us-west-1.amazonaws.com/noble-trios/helium-father-reads.fq.gz -o father.fq.gz
     curl -L https://s3-us-west-1.amazonaws.com/noble-trios/helium-proband-reads.fq.gz -o proband.fq.gz
     curl -L https://s3-us-west-1.amazonaws.com/noble-trios/helium-refr.fa.gz -o refr.fa.gz
     bwa index refr.fa.gz

     kevlar simplex \
         --case proband.fq.gz --control mother.fq.gz --control father.fq.gz \
         --novel-memory 50M --filter-memory 1M --filter-fpr 0.005 --mask-memory 5M  \
         --mask-files refr.fa.gz \
         --threads 4 --ksize 31 \
         --out variant-calls.vcf \
         refr.fa.gz
