Comprehensive command-line interface reference
==============================================

The kevlar command-line interface is designed around a single command :code:`kevlar`.
From this one command, a variety of tasks and procedures can be invoked using several *subcommands*.

Once kevlar is installed, available subcommands can be listed by executing :code:`kevlar -h`.
To see instructions for running a specific subcommand, execute :code:`kevlar <subcommand> -h` (of course replacing :code:`subcommand` with the actual name of the subcommand).


.. _kevlar_count_api:

kevlar count
------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: count

.. _kevlar_novel_api:

kevlar novel
------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: novel

.. _kevlar_filter_api:

kevlar filter
-------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: filter

.. _kevlar_partition_api:

kevlar partition
----------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: partition

.. _kevlar_assemble_api:

kevlar assemble
---------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: assemble

.. _kevlar_localize_api:

kevlar localize
---------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: localize

.. _kevlar_call_api:

kevlar call
-----------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: call

.. _kevlar_simlike_api:

kevlar simlike
--------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: simlike

.. _kevlar_alac_api:

kevlar alac
-----------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: alac

.. _kevlar_unband_api:

kevlar unband
-------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: unband

.. _kevlar_augment_api:

kevlar augment
----------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: augment

.. _kevlar_mutate_api:

kevlar mutate
-------------

.. argparse::
   :module: kevlar.cli.__init__
   :func: parser
   :nodefault:
   :prog: kevlar
   :path: mutate
