Running **kevlar**
==================

The **kevlar** software implements a Python library for genetic sequence and variant analysis.
**kevlar**'s primary interface is invoked via the command line, but is also designed so that it can be seamlessly integrated into third-party Python programs.


Command line interface
----------------------

Once installed, the **kevlar** software can be invoked from the shell using the ``kevlar`` command.
The **kevlar** command line interface (CLI) uses the *subcommand* pattern, in which a single master command supports several different operations by defining multiple subcommands (such as ``kevlar novel`` and ``kevlar partition``).
Comprehensive documentation of the **kevlar** CLI is available :doc:`here <cli>`.

Starting with version 1.0, the CLI will be under `semantic versioning <http://semver.org/>`_.


Python interface
----------------

As a result of **kevlar**'s design to facilitate internal testing, the "main method" of each **kevlar** subcommand can easily be executed programmatically.
The following example shows how to execute ``kevlar reaugment`` from a standalone Python program.

.. code:: python

   import kevlar

   # Declare arguments just like you would on the command line
   arglist = ['reaugment', '-o', 'new.augfastq', 'old.augfastq', 'new.fastq']

   args = kevlar.cli.parser().parse_args(arglist)
   kevlar.reaugment.main(args)

Other units of code in the **kevlar** package may also be amenable to importing and executing programmatically.
However, the code internals are not under semantic versioning and by necessity will be less stable and have poorer documentation.
Have fun and knock yourself out, but be prepared for changes in internal behavior in subsequent releases!
