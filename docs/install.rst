Installing **kevlar**
=====================

For the impatient
-----------------

.. code::

    virtualenv kevlar-env
    source kevlar-env/bin/activate  # Execute this with every new terminal session
    # pip install pysam  # You probably don't need this
    pip install networkx
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

Upcoming releases of **kevlar** will install dependencies automatically, but for now these must be installed manually.
Currently the only hard non-standard dependency is the **khmer** library—to be precise, its :code:`feature/consume_bitsplit` development branch.
The :code:`kevlar dump` command also depends on **pysam**, while :code:`kevlar filter` and :code:`kevlar assemble` depend on **networkx**.

.. code::

    pip install pysam
    pip install networkx
    pip install git+https://github.com/dib-lab/khmer.git

The **pysam** dependency will eventually be dropped.

Installation
------------

Once the dependencies are installed, **kevlar** can be installed with the :code:`pip` command.

.. code::

    pip install biokevlar

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
