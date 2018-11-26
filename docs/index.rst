.. raw:: html

    <img src="_static/kevlar-logo.png"alt="kevlar logo" style="height: 150px; display: block" />

The **kevlar** library is a testbed for developing reference-free variant discovery methods for genomics.
The initial focus of the project has been discovery of *de novo* germline variants in simplex pedigrees, but we are laying the groundwork to support a wider range of experimental designs such as case/control studies.

We advertise kevlar as a "reference-free" method in that it does not require a reference genome to identify variant-spanning reads or to assemble these reads into contigs.
However, the reference genome is currently still used for making the final variant call by aligning each assembled contig to a small cutout of the reference genome.
One of the aspirations of the project moving forward is to reduce and eventually eliminate dependence on a reference genome completely.

kevlar is currently under heavy development, and internal features are not yet stable.
However, the core features and workflows are reasonably well tested, and leverage software components from various third-party libraries that are very well tested and widely used.

Documentation for **kevlar**
============================

.. toctree::
   :maxdepth: 1

   install
   running
   quick-start
   tutorial
   terms
   formats
   banding
   sim
   cli
   conduct

Links
-----

- `Github repository <https://github.com/dib-lab/kevlar>`_
- `License <https://github.com/dib-lab/kevlar/blob/master/LICENSE>`_
