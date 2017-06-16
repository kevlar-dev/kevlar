Installing **kevlar**
=====================

For the impatient
-----------------

.. code::

    virtualenv kevlar-env
    source kevlar-env/bin/activate  # Execute this with every new terminal session
    pip install git+https://github.com/dib-lab/khmer.git
    pip install biokevlar

Virtual environment
-------------------

We recommend installing **kevlar** and its dependencies in a dedicated `virtual environment <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_.
The command :code:`virtualenv kevlar-env` will create a new virtual environment, and only needs to be executed once.
The command :code:`source kevlar-env/bin/activate` will need to be re-executed any time you open a new session in your terminal.

.. note:: The label :code:`kevlar-env` can be replaced with an alternative label if desired.

Dependencies
------------

The **kevlar** software has three non-standard dependencies: the `networkx package <https://networkx.github.io/>`_, the `pysam <http://pysam.readthedocs.io/>`_ package, and the `khmer package <http://khmer.readthedocs.io/>`_.
Currently, the **kevlar** installation procedure handles the first two dependencies automatically, but since it relies on an unreleased version of **khmer** this last dependency must be install manually.

.. code::

    pip install git+https://github.com/dib-lab/khmer.git

Installation
------------

Once **khmer** is installed, **kevlar** can be installed with the :code:`pip` command.

.. code::

    pip install biokevlar

This installs the most recent stable release.
If you want to install the latest (possibly unstable) version, pip can install **kevlar** directly from GitHub.

.. code::

    pip install git+https://github.com/dib-lab/kevlar.git

If you want to test whether kevlar is installed and running correctly, use pytest.

.. code::

    pip install pytest
    pytest --pyargs kevlar

Development environment
-----------------------

If you'd like to contribute to **kevlar**'s development or simply poke around, the source code can be cloned from Github.
In addition to the dependencies listed above, a few additional dependencies are required for a complete development environment.
These can be installed with `make` for your convenience.

.. code::

    git clone https://github.com/dib-lab/kevlar.git
    cd kevlar
    make devenv

Hack away!
Feel free to ask questions or submit bug reports to the **kevlar** `issue tracker <https://github.com/dib-lab/kevlar/issues>`_.
