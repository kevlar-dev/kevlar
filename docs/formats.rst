File formats in **kevlar**
==========================

Although kevlar performs many operations on *k*-mers, read sequences are the primary currency of exchange between different stages of the analysis workflow.
kevlar supports reading from and writing to Fasta and Fastq files, and treats these identically since it does not use any base call quality information.
In most cases, kevlar should also be able to automatically detect whether an input file is gzip-compressed or not and handle it accordingly (no bzip2 support).

Broken-paired Fasta / Fastq files
---------------------------------

While kevlar does not require pairing information for variant discovery, it can be helpful in the final stages of variant calling.
The bioinformatics community uses two common conventions for encoding pairing information: paired files and interleaved files.
In paired files, the first record in file1 is paired with the first record in file2, the second record in file1 is paired with the second record in file2, and so on.
In an interleaved file, the first record is paired with the second record, the third record is paired with the fourth record, and so on.

If you want to retain and use pairing information, kevlar only supports reading pairing information from interleaved files.
kevlar also supports "broken paired" files, where single-end/orphaned reads are occasionally scattered in between paired reads in an interleaved file.

Augmented sequences
-------------------

"Interesing *k*-mers" are putatively novel *k*-mers that are high abundance in the proband/case sample(s) and effectively absent from control samples.
To facilitate reading and writing these "interesting *k*-mers" along with the reads to which they belong, kevlar uses an *augmented* version of the Fasta and Fastq formats.
Here is an example of an augmented Fastq file.

.. highlight:: none

.. code::

   @read1
   TTTTACCCGATGGGCGAGGTGAAATACTATGCCGATTTATTCTTACACAATTAAATTGCTAGTCCGGTTAGGGTTAGTTTGCGGCCTTCGTTCCAGCGCCGTGTT
   +
   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
        CCCGATGGGCGAGGTGAAA          18 1 0#
                                                                        AGGGTTAGTTTGCGGCCTT          11 0 0#
   @read2
   AAGAGATTGTCGCTTGCCCCGTAAAGGAATTAGACCGGGCGACCAGAGCCTATTAGTAGCCCGCGCCTGTAGCACAAACGACTTTCGTACTATTATTAGACGTCG
   +
   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
    AGAGATTGTCGCTTGCCCC          14 0 1#
     GAGATTGTCGCTTGCCCCG          12 0 0#
      AGATTGTCGCTTGCCCCGT          14 0 0#
   @read3
   GAGACCATAAACCAGCTCTTGGTACCGAAAGAACACCTATGAATAACCGTGAGTGCATGATTCCTGTGAAGAGATTGTCGCTTGCCCCGTAAAGGAATTAGACCG
   +
   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
                  CTCTTGGTACCGAAAGAAC          19 1 0#
                                                                        AGAGATTGTCGCTTGCCCC          14 0 1#
                                                                         GAGATTGTCGCTTGCCCCG          12 0 0#
                                                                          AGATTGTCGCTTGCCCCGT          14 0 0#
   @read4
   TCCGGTTAGGGTTAGTTTGCGGCCTTCGTTCCAGCGCCGTGTTGTTGCAATTTAATCCCGAGAAACCTCATGTAGCGGCTACTGGACCGCTGGGTAAGCTCAGAC
   +
   BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
          AGGGTTAGTTTGCGGCCTT          11 0 0#

As with a normal Fastq file, each record contains 4 lines to declare the read sequence and qualities.
However, these 4 lines are followed by one or more lines indicating the "interesing *k*-mers", showing their sequence followed by their abundance in each sample (case first, then controls), with a ``#`` as the final character.
Augmented Fastq files are easily converted to normal Fastq files by invoking a command like ``grep -v '#$' reads.augfastq > reads.fastq`` (same for augmented Fasta files).

The functions ``kevlar.parse_augmented_fastx`` and ``kevlar.print_augmented_fastx`` are used internally to read and write augmented Fastq/Fasta files.
However, these functions can easily be imported and called from third-party Python scripts as well.

Mate sequences
--------------

Although kevlar does not require pairing information, it can be used to improve calling when it's available.
The augmented Fastq/Fasta format also allows mate sequences to be associated with each record.
If a contig assembled from novel reads maps to multiple regions of the reference genome with the same score, this pairing information can be used to predict the most likely true variant.

The ``mateseq`` annotation should be placed after the first 4 lines of the record, and shown in the two records below.

.. code::

    @DraconisOccidentalisRead12/1
    CAGGTAGTGTGATGCCTCCAGCTTTGTTCTTTTGGCTTAGGATTGACTTGGCAATTCGGGCTCTTTTTTGGTTCCATATGAACTTTAAAGTAGTTTTTTC
    +
    8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
                              TTCTTTTGGCTTAGGATTGACTTGGCAATTC          7 0 1#
                               TCTTTTGGCTTAGGATTGACTTGGCAATTCG          5 0 1#
                                CTTTTGGCTTAGGATTGACTTGGCAATTCGG          5 0 1#
                                 TTTTGGCTTAGGATTGACTTGGCAATTCGGG          5 1 1#
                                  TTTGGCTTAGGATTGACTTGGCAATTCGGGC          5 0 1#
    #mateseq=CTGATAAGCAACTTCAGCAAAGTCTCAGGATACAAAATCAATGTACGAAAATCACAAGCGTTCTTATACACCAACAACAGACAAACAGAGAGCCAAATCA#
    @DraconisOccidentalisRead56/1
    TCTTGAATTCCCATGTGTTGTGGGAGGGACCCATTGGGAGGTAATTGAATCATGGGGGCACGTCTTTCCCATGCTGTTCTCATGATAGAGACTAAGTCTC
    +
    8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
         AATTCCCATGTGTTGTGGGAGGGACCCATTG          5 0 1#
          ATTCCCATGTGTTGTGGGAGGGACCCATTGG          5 0 1#
    #mateseq=ATTAGAAAAAAAAAGTGCATTCGTAAATGTCATAACAATAAAATTATACTCCAAGACTTTGTACAAGATGAAAGTAATATGAAGAAGGGGCTACAGGAAA#
           TTCCCATGTGTTGTGGGAGGGACCCATTGGG          5 0 1#
            TCCCATGTGTTGTGGGAGGGACCCATTGGGA          5 0 1#
