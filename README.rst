|kevlar build status| |PyPI version| |Test coverage| |kevlar documentation| |Docker build status| |MIT licensed|

.. figure:: docs/_static/morpheus-kevlar.jpg
   :alt: What if I told you we don't need alignments to find variants?
   :scale: 75%

kevlar
======

Daniel Standage, 2016
https://kevlar.readthedocs.io

Welcome to **kevlar**, software for predicting *de novo* genetic variants without mapping reads to a reference genome!
kevlar's k-mer abundance based method calls single nucleotide variants (SNVs) as well as short, medium and long insertion/deletion variants (indels) simultaneously.
This software is free for use under the MIT license.

.. raw:: html

    <details>
        <summary>Where can I find kevlar online?</summary>
        <ul>
            <li>Source repository: https://github.com/dib-lab/kevlar</li>
            <li>Documentation: https://kevlar.readthedocs.io</li>
            <li>Stable releases: https://github.com/dib-lab/kevlar/releases</li>
            <li>Issue tracker: https://github.com/dib-lab/kevlar/issues</li>
        <ul>
        If you have questions or need help with kevlar, the <a href="https://github.com/dib-lab/kevlar">GitHub issue tracker</a> should be your first point of contact.
    </details>



How do I use kevlar?
--------------------

- Installation instructions: http://kevlar.readthedocs.io/en/latest/install.html
- Quick start guide: http://kevlar.readthedocs.io/en/latest/quick-start.html
- Tutorial: http://kevlar.readthedocs.io/en/latest/tutorial.html

**Note**: kevlar is currently focused almost entirely on finding novel germline variants in related individuals.
We hope to benchmark kevlar a wider range of experimental designs soon.

Contributing
------------

We welcome contributions to kevlar from the community!
If you're interested in modifying kevlar or contributing to its ongoing development, feel free to send us a message or submit a pull request!.

The kevlar software is a project of the `Lab for Data Intensive Biology <http://ivory.idyll.org/lab/>`__ at UC Davis.

.. |kevlar build status| image:: https://img.shields.io/travis/dib-lab/kevlar.svg
   :target: https://travis-ci.org/dib-lab/kevlar
   :alt: Status of the continuous integration build
.. |PyPI version| image:: https://img.shields.io/pypi/v/biokevlar.svg
   :target: https://pypi.python.org/pypi/genhub
   :alt: PyPI version
.. |Test coverage| image:: https://img.shields.io/codecov/c/github/dib-lab/kevlar.svg
   :target: https://codecov.io/github/dib-lab/kevlar
   :alt: Current code coverage from automated tests
.. |kevlar documentation| image:: https://readthedocs.org/projects/kevlar/badge/?version=latest&maxAge=900
   :target: http://kevlar.readthedocs.io/en/latest/?badge=latest
   :alt: Project documentation
.. |Docker build status| image:: https://quay.io/repository/dib-lab/kevlar/status
   :target: https://quay.io/repository/dib-lab/kevlar
   :alt: Docker image for cloud deployment
.. |MIT licensed| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://github.com/dib-lab/kevlar/blob/master/LICENSE
   :alt: Software license
