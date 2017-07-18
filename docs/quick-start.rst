Quick start
===========

Currently, **kevlar** development is focused heavily on trio and quad experimental designs.
This document gives a bare-bones walkthrough of the **kevlar** workflow, from raw data to contig assembly.
Final variant calling comming soon.

Eventually, this should hopefully be cleaned up into a smaller sequence of commands, but for now things are changing frequently enough that thorough documentation isn't yet feasible.

----------

If you have not already done so, install **kevlar** using :doc:`the following instructions <install>`.

A complete listing of all available configuration options for each script can be found in :doc:`the CLI documentation <cli>`, or by executing ``kevlar <subcommand> -h`` in the terminal.

----------

#. Compute k-mer abundances

   .. code:: bash

       kevlar count \
           --case case.counttable case-1.fq.gz case-2.fq.gz \
           --control control1.counttable control-a-1.fq.gz control-a-2.fq.gz \
           --control control2.counttable control-b-1.fq.gz control-b-2.fq.gz \
           --ksize 21 \
           --ctrl-max 0

#. Find "interesting" (potentially novel) k-mers

   .. code:: bash

       kevlar novel \
           --case-counts case-1.counttable \
           --case-min 8 \
           --controls control-1.counttable control-2.counttable \
           --ctrl-max 0 \
           --ksize 21 \
           --out case-1.novel.unfiltered.augfastq.gz
           case-1.fq.gz

#. Recompute k-mer abundances to discard false positives

   .. code:: bash

       kevlar filter \
           --refr refr.fa.gz \
           --contam contaminants.fa \
           --min-abund 8 \
           --ksize 21 \
           --aug-out case-1.novel.filtered.augfastq.gz \
           --out case-1.novel.filtered.fq.gz \
           case-1.novel.unfiltered.augfastq.gz

#. Partition reads by shared novel k-mers

   .. code:: bash

       kevlar partition case-1-partition case-1.novel.filtered.augfastq.gz

#. Assemble partitioned reads

   .. code:: bash

       kevlar assemble \
           --out case-1-partition.cc0.augfasta \
           case-1-partition.cc0.augfastq.gz
