Quick start
===========

If you have not already done so, install **kevlar** using :doc:`the following instructions <install>`.

A complete listing of all available configuration options for each script can be shown by executing ``kevlar <subcommand> -h`` in the terminal.

----------

#. Compute k-mer abundances

   .. code:: bash

       kevlar count \
           --case kevlar/tests/data/trio1/case1.fq \
           --controls kevlar/tests/data/trio1/ctrl[1,2].fq \
           --ksize 21 \
           --ctrl_max 0

#. Find "interesting" (potentially novel) k-mers

   .. code:: bash

       kevlar novel \
           --cases kevlar/tests/data/trio1/case1.counttable \
           --case_min 8 \
           --controls kevlar/tests/data/trio1/ctrl[1,2].counttable \
           --ctrl_max 0 \
           --ksize 21 \
           --out case1.novel.unfiltered.augfastq.gz
           kevlar/tests/data/trio1/case1.fq

#. Recompute k-mer abundances to discard false positives, partition reads by shared novel k-mers

   .. code:: bash

       kevlar filter \
           --refr kevlar/tests/data/bogus-genome/refr.fa \
           --contam kevlar/tests/data/bogus-genome/contam1.fa \
           --min-abund 8 \
           --ksize 21 \
           --aug-out case1.novel.filtered.augfastq.gz \
           --out case1.novel.filtered.fq.gz \
           --cc-prefix case1 \
           case1.novel.unfiltered.augfastq.gz

   Currently partitioning is done by ``kevlar filter``, but this will soon be handled by a dedicated ``kevlar partition`` command.

#. Assemble partitioned reads

   .. code:: bash

       kevlar assemble --out case1.cc0.augfasta case1.cc0.augfastq.gz
