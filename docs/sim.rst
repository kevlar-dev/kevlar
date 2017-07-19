Simulating variants with **kevlar**
===================================

To facilitate testing, **kevlar** implements a simple command to apply simulated "mutations" to a reference sequence.
We have used this internally to simulate data sets for testing, to verify that **kevlar** can recover the simulated "mutation" or variant.

The command-line interface for ``kevlar mutate`` is very simple (for full details see `the CLI documentation <cli.html#kevlar-mutate>`_).

The "mutation file" format is described here by way of example.

.. code::

   seq1	2345915	del	141
   seq1	1022305	snv	2
   seq1	2062327	inv	429
   seq1	1234310	del	32
   seq1	388954	ins	TGTTTCCTTTCATACCCCACCAC
   seq1	2460047	snv	2

The mutations file is a plain text tabular data file with four fields separated by spaces or tabs.

- sequence ID
- variant starting position (0-based)
- variant type (currently supported types: ``snv``, ``ins``, ``del``, and ``inv`` for single-nucleotide variants, insertions, deletions, and inversions)
- value; represents lexicographic offset for SNVs, variant length for deletions and inversions, and inserted sequence for insertions.
