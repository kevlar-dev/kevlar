.. raw:: html

    <img src="_static/kevlar-logo.png"alt="kevlar logo" style="height: 150px; display: block" />

The **kevlar** library is a testbed for developing reference-free variant discovery methods for genomics.
The initial focus of the project has been novel germline variant discovery in simplex pedigrees, but we are laying the groundwork to support more generalized case/control studies of large cohorts.

kevlar does not use a reference genome to identify variant associated reads or to assemble these reads into contigs representing each variant.
However, the reference genome is still used for making the final variant call by aligning each assembled contig to a small cutout of the reference genome.
The reference genome can also be used at various preliminary stages for filtering and improving kevlar's performance.

kevlar is currently under heavy development, and some features are not yet stable.
However, the core features and workflows are reasonably well tested, and leverage software components from various third-party libraries that are very well tested and widely used.

Documentation for **kevlar**
============================

.. toctree::
   :maxdepth: 1

   install
   running
   quick-start
   tutorial
   formats
   sim
   cli
   conduct

Links
-----

- `Github repository <https://github.com/dib-lab/kevlar>`_
- `License <https://github.com/dib-lab/kevlar/blob/master/LICENSE>`_
