Comprehensive command-line interface reference
==============================================

The kevlar command-line interface is designed around a single command :code:`kevlar`.
From this one command, a variety of tasks and procedures can be invoked using several *subcommands*.

Once kevlar is installed, available subcommands can be listed by executing :code:`kevlar -h`.
To see instructions for running a specific subcommand, execute :code:`kevlar <subcommand> -h` (of course replacing :code:`subcommand` with the actual name of the subcommand).

kevlar count
------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: count

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

kevlar cutout
-------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: cutout

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

kevlar simlike
--------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: simlike

kevlar alac
-----------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: alac

kevlar simplex
--------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: simplex

kevlar dump
-----------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: dump

kevlar augment
----------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: augment

kevlar mutate
-------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: mutate
