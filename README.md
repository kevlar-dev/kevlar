[![kevlar build status][travisbadge]](https://travis-ci.org/dib-lab/kevlar)
[![PyPI version][pypibadge]](https://pypi.python.org/pypi/biokevlar)
[![Test coverage][codecovbadge]](https://codecov.io/github/dib-lab/kevlar)
[![kevlar documentation][rtdbadge]](http://kevlar.readthedocs.io/en/latest/?badge=latest)
[![Docker build status][dockerbadge]](https://quay.io/repository/dib-lab/kevlar)
[![MIT licensed][licensebadge]](https://github.com/dib-lab/kevlar/blob/master/LICENSE)

<img src="docs/_static/morpheus-kevlar.jpg" alt=" What if I told you we don't need alignments to find variants?" width="400px" />

# kevlar

Daniel Standage, 2016  
https://kevlar.readthedocs.io

Welcome to **kevlar**, software for predicting *de novo* genetic variants without mapping reads to a reference genome!
kevlar's k-mer abundance based method calls single nucleotide variants (SNVs) as well as short, medium and long insertion/deletion variants (indels) simultaneously.
This software is free for use under the MIT license.

<details>
  <summary>Where can I find kevlar online?</summary>
  <ul>
    <li>Source repository: https://github.com/dib-lab/kevlar</li>
    <li>Documentation: https://kevlar.readthedocs.io</li>
    <li>Stable releases: https://github.com/dib-lab/kevlar/releases</li>
    <li>Issue tracker: https://github.com/dib-lab/kevlar/issues</li>
  </ul>

  If you have questions or need help with kevlar, the [GitHub issue tracker](https://github.com/dib-lab/kevlar) should be your first point of contact.
</details>

<details>
  <summary>How do I install kevlar?</summary>

  See [the kevlar documentation](http://kevlar.readthedocs.io/en/latest/install.html) for complete instructions, but the impatient can try the following.

  ```
  pip3 install git+https://github.com/dib-lab/khmer.git
  pip3 install biokevlar
  ```
</details>

<details>
  <summary>How do I use kevlar?</summary>
  <ul>
    <li>Installation instructions: http://kevlar.readthedocs.io/en/latest/install.html</li>
    <li>Quick start guide: http://kevlar.readthedocs.io/en/latest/quick-start.html</li>
    <li>Tutorial: http://kevlar.readthedocs.io/en/latest/tutorial.html</li>
  </ul>
</details>

<details>
  <summary>How can I contribute?</summary>
  
  We welcome contributions to the kevlar project from the community!
  If you're interested in modifying kevlar or contributing to its ongoing development, feel free to send us a message or submit a pull request!

  The kevlar software is a project of the [Lab for Data Intensive Biology](http://ivory.idyll.org/lab/) and the [Computational Genomics Lab](http://www.hormozdiarilab.org/) at UC Davis.
</details>


[travisbadge]: https://img.shields.io/travis/dib-lab/kevlar.svg
[pypibadge]: https://img.shields.io/pypi/v/biokevlar.svg
[codecovbadge]: https://img.shields.io/codecov/c/github/dib-lab/kevlar.svg
[rtdbadge]: https://readthedocs.org/projects/kevlar/badge/?version=latest&maxAge=900
[dockerbadge]: https://quay.io/repository/dib-lab/kevlar/status
[licensebadge]: https://img.shields.io/badge/license-MIT-blue.svg
