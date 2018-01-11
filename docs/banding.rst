*K*-mer banding
===============

If memory is a limiting factor for *k*-mer counting, kevlar supports a scatter/gather approach to based on a strategy we call "*k*-mer banding."
In brief, kevlar can achieve an N-fold reduction in memory usage in exchange for counting *k*-mers N batches.
For each batch, kevlar ignores all *k*-mers except those whose hash values fall within a specified numerical range (band), reducing the memory required to achieve accurate *k*-mer counts.

The ``kevlar count``, ``kevlar effcount``, and ``kevlar novel`` commands support *k*-mer banding.
The output of multiple ``kevlar novel`` invocations can be combined using ``kevlar filter``.
