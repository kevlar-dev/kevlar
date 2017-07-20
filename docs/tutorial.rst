Tutorial
========

Sadly, this "tutorial" offers little more than Unix commands by way of exposition.
Once the software is a bit more stable we will invest in a more comprehensive tutorial.

.. code:: bash

    # 1. Dump reads that match reference genome perfecly
    kevlar dump GRCh38_full_analysis_set_plus_decoy_hla.fa NA91238.bam | gzip -c > NA19238.dump.fq.gz
    kevlar dump GRCh38_full_analysis_set_plus_decoy_hla.fa NA91239.bam | gzip -c > NA19239.dump.fq.gz
    kevlar dump GRCh38_full_analysis_set_plus_decoy_hla.fa NA91240.bam | gzip -c > NA19240.dump.fq.gz


    # 2. Count k-mers, find interesting k-mers, print out corresponding reads

    ## Option a: the old way
    kevlar novel \
        --case NA19240.dump.fq.gz \
        --control NA19238.dump.fq.gz \
        --control NA19239.dump.fq.gz \
        --case-min 11 --ctrl-max 0 \
        --memory 24G --max-fpr 0.5 \
        --ksize 31 \
        --out NA1940.novel.augfastq.gz

    ## Option b: the new way (not well tested)
    kevlar count \
        --case NA19240.counttable NA19240.dump.fq.gz \
        --control NA19238.counttable NA19238.dump.fq.gz
        --control NA19239.counttable NA19239.dump.fq.gz \
        --ctrl-max 0 \
        --memory 24G --max-fpr 0.5 \
        --ksize 31

    kevlar novel \
        --case-counts NA19240.counttable \
        --control-counts NA19238.counttable NA19239.counttable \
        --case-min 11 --ctrl-max 0 \
        --ksize 31
        --out NA1940.novel.augfastq.gz


    # 3. Filter out false positives, contamination
    kevlar filter \
        --refr GRCh38_full_analysis_set_plus_decoy_hla.fa --refr-memory 3G --refr-max-fpr 0.01 \
        --contam viral-contam.fa.gz --contam-memory 100M --contam-max-fpr 0.01 \
        --abund-memory 1G --abund-max-fpr 0.01 --min-abund 11 \
        --ksize 31 \
        --out NA19240.filtered.fq.gz --augout NA19240.filtered.augfastq.gz \
        NA1940.novel.augfastq.gz


    # 4. Get rid of super high abundance stuff
    load-into-counting.py --ksize 31 --max-memory-usage 1G NA19240.filtered.counttable NA19240.filtered.fq.gz

    slice-reads-by-coverage.py \
        --min-coverage 2 --max-coverage 250 \
        NA19240.filtered.counttable NA19240.filtered.fq.gz NA19240.sliced.fq

    kevlar reaugment --out NA19240.sliced.augfastq.gz NA19240.filtered.augfastq.gz NA19240.sliced.fq


    # 5. Group reads by shared interesting k-mers
    mkdir -p NA19240-ccs/
    kevlar partition NA19240-ccs/NA19240 NA19240.sliced.augfastq.gz


    # 6. Abundance trimming, assembly, and variant calling
    for i in {0..1234}  # replace 1234 with num(partitions)-1
    do
        gunzip -c NA19240-ccs/NA19240${i}.augfastq.gz | grep -v '#$' > NA19240-ccs/NA19240${i}.fq
        trim-low-abund.py -M 1M -k 31 -o NA19240-ccs/NA19240${i}-trim.fq NA19240-ccs/NA19240${i}.fq
        kevlar reaugment --out NA19240-ccs/NA19240${i}-trim.augfastq.gz NA19240-ccs/NA19240${i}.augfastq.gz NA19240-ccs/NA19240${i}.fq
        kevlar assemble --out NA19240-ccs/NA19240${i}-assembled.augfasta NA19240-ccs/NA19240${i}-trim.augfastq.gz
        kevlar localize --out NA19240-ccs/NA19240${i}-refr.fa NA19240-ccs/NA19240${i}-assembled.augfasta GRCh38_full_analysis_set_plus_decoy_hla.fa
        kevlar call --out NA19240-ccs/NA19240${i}-call.txt NA19240-ccs/NA19240${i}-assembled.augfasta NA19240-ccs/NA19240${i}-refr.fa
    done
