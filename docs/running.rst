Running **kevlar**
==================

The kevlar library can be invoked via a command-line interface or a Python API.

Command line interface
----------------------

Once installed, the kevlar software can be invoked from the shell using the ``kevlar`` command.
The kevlar command line interface (CLI) uses the *subcommand* pattern, in which a single master command supports several different operations by defining multiple subcommands (such as ``kevlar novel`` and ``kevlar partition``).
Comprehensive documentation of the kevlar CLI is available :doc:`here <cli>`.

Starting with version 1.0, the CLI will be under `semantic versioning <http://semver.org/>`_.


Python interface
----------------

The main procedure of each kevlar subcommand is implemented as a generator function and can be executed programmatically.
The following example shows how a standalone script or notebook would invoke the ``kevlar partition`` procedure.

.. code:: python

    fh = kevlar.open('novel-reads.augfastq.gz', 'r')
    readstream = kevlar.parse_augmented_fastx(fh)
    numreads_per_partition = list()
    for part in kevlar.partition.partition(readstream, maxabund=200):
        numreads = len(part)
        numreads_per_partition.append(numreads)
        # do some other uber cool calculations
    plt.hist(numreads_per_partition, bins=25)

It is also possible to mimic the behavior of the CLI with Python code.
The following example shows how to execute ``kevlar reaugment`` from a standalone script or notebook.

.. code:: python

   arglist = ['reaugment', '--out', 'new.augfastq', 'old.augfastq', 'new.fastq']
   args = kevlar.cli.parser().parse_args(arglist)
   kevlar.reaugment.main(args)

Other units of code in the kevlar package may also be useful for standalone Python programs.
However, the Python API is not yet slated for semantic versioning and is not as stable or well documented as the CLI.
Have fun and knock yourself out, but be prepared for changes in internal behavior in subsequent releases!
