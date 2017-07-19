File formats in **kevlar**
==========================

Although **kevlar** performs many operations on *k*-mers, read sequences are the primary currency of exchange between different stages of the analysis workflow.
**kevlar** supports reading from and writing to Fasta and Fastq files, and treats these identically since it does not use any base call quality information.
In most cases, **kevlar** should also be able to automatically detect whether an input file is gzip-compressed or not and handle it accordingly (no bzip2 support).

Augmented sequences
-------------------

"Interesing *k*-mers" are putatively novel *k*-mers that are high abundance in the proband/case sample(s) and effectively absent from control samples.
To facilitate reading and writing these "interesting *k*-mers" along with the reads to which they belong, **kevlar** uses an *augmented* version of the Fasta and Fastq formats.
Here is an example of an augmented Fastq file.

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
