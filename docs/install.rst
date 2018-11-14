Installing **kevlar**
=====================

For the impatient
-----------------

If this isn't your first time in the wing, the following commands should be sufficient to install kevlar in the majority of cases.
Otherwise, we suggest reading through the entire installation instructions before beginning.

.. code::

    pip3 install pysam networkx pandas scipy git+https://github.com/dib-lab/khmer.git
    pip3 install biokevlar

Virtual environment
-------------------

We recommend installing kevlar and its dependencies in a dedicated virtual environment using `venv <https://docs.python.org/3/library/venv.html>`_ or `conda <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_.

- If you use ``venv``, the command ``python3 -m venv kevlar-env`` will create a new virtual environment, and only needs to be executed once.
  The command ``source kevlar-env/bin/activate`` will need to be re-executed any time you open a new session in your terminal.
- If you use ``conda``, the command ``conda create --name kevlar-env python=3.6`` will create a new virtual environment, and only needs to be executed once.
  The command ``source activate kevlar-env`` will need to be re-executed any time you open a new session in your terminal.

.. note:: The label ``kevlar-env`` can be replaced with an alternative label if desired.

Dependencies
------------

The kevlar package **requires Python 3** and has several dependencies that are not in the standard Python libraries.

- the `networkx package <https://networkx.github.io/>`_
- the `pysam module <http://pysam.readthedocs.io/>`_
- the `pandas library <http://pandas.pydata.org/>`_
- the `scipy library <https://www.scipy.org/>`_
- the `khmer package <http://khmer.readthedocs.io/>`_

Also, kevlar requires the `bwa <https://github.com/lh3/bwa>`_ and `samtools <https://github.com/samtools/samtools>`_ commands to be callable from your ``$PATH`` environmental variable.

When kevlar is installed from PyPI, most Python dependencies *should* handled automatically.
But since kevlar currently relies on an unreleased version of khmer this last dependency must be installed manually.

.. code::

    pip3 install git+https://github.com/dib-lab/khmer.git

.. note::

    According to `PEP 394 <https://www.python.org/dev/peps/pep-0394/>`_ a Python 3 distribution should include a ``pip3`` command for package management, but in some configurations this may not be true.
    If you've confirmed that Python 3 is installed correctly, you're probably safe using the ``pip`` command if ``pip3`` is unavailable.


.. note::

   In some cases pip cannot install all dependencies automatically, and so manual installation is required.

   .. code::

      pip3 install pysam>=0.11.2 networkx>=2.0 pandas scipy git+https://github.com/dib-lab/khmer.git

Installation
------------

Once the prerequisites are installed, kevlar can be installed with the pip.

.. code::

    pip3 install biokevlar

This installs the most recent stable release.
If you want to install the latest (possibly unstable) version, pip can install kevlar directly from GitHub.

.. code::

    pip3 install git+https://github.com/dib-lab/kevlar.git

To test whether kevlar is installed and running correctly, use `pytest <https://docs.pytest.org/>`_.

.. code::

    pip3 install pytest
    pytest --pyargs kevlar.tests

Development environment
-----------------------

If you'd like to contribute to kevlar's development or simply poke around, the source code can be cloned from Github.
In addition to the dependencies listed above, a few additional dependencies are required for a complete development environment.
These can be installed with ``make`` for your convenience.

.. code::

    git clone https://github.com/dib-lab/kevlar.git
    cd kevlar
    make devenv
    pip3 install -e .

Hack away!
Feel free to ask questions or submit bug reports to the kevlar `issue tracker <https://github.com/dib-lab/kevlar/issues>`_.
