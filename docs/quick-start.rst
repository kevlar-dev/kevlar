Quick start
===========

If you have not already done so, install kevlar using :doc:`the following instructions <install>`.

This gives a crash course on running kevlar's simplex analysis workflow.
A :doc:`more detailed tutorial is available <tutorial>`, and a complete listing of all available configuration options for each script can be found in :doc:`the CLI documentation <cli>`, or by executing ``kevlar <subcommand> -h`` in the terminal.

----------

   .. code:: bash

        curl -L https://osf.io/e8jb3/download?version=1 -o mother.fq.gz
        curl -L https://osf.io/fuaty/download?version=1 -o father.fq.gz
        curl -L https://osf.io/f5trh/download?version=1 -o proband.fq.gz
        curl -L https://osf.io/58rwa/download?version=1 -o refr.fa.gz
        bwa index refr.fq.gz

        kevlar simplex \
            --case proband.fq.gz --case-min 6 \
            --control mother.fq.gz --control father.fq.gz --ctrl-max 0 \
            --novel-memory 5M --novel-fpr 0.6 --threads 4 \
            --filter-memory 1M --filter-fpr 0.005 \
            --mask-files refr.fa.gz --mask-memory 5M \
            --ksize 25 --out variant-calls.vcf
            refr.fa.gz
