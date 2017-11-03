Comprehensive command-line interface reference
==============================================

The kevlar command-line interface is designed around a single command :code:`kevlar`.
From this one command, a variety of tasks and procedures can be invoked using several *subcommands*.

Once kevlar is installed, available subcommands can be listed by executing :code:`kevlar -h`.
To see instructions for running a specific subcommand, execute :code:`kevlar <subcommand> -h` (of course replacing :code:`subcommand` with the actual name of the subcommand).

kevlar simplex
--------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: simplex

kevlar alac
-----------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: alac

kevlar novel
------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: novel

kevlar filter
-------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: filter

kevlar assemble
---------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: assemble

kevlar localize
---------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: localize

kevlar call
-----------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: call

kevlar count
------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: count

kevlar dump
-----------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: dump

kevlar reaugment
----------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: reaugment

kevlar mutate
-------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: mutate
