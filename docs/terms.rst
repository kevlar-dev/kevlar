Terminology
===========

Here is a brief glossary of terms used throughout the code and project.

- **novel k-mer**: a *k*-mer linked (or putatively linked) to a novel germline variant
- **ikmer**: short for *interesting k-mer*, a synonym for **novel k-mer** with an emphasis on its unverified status
- **partition**: a set of reads that share novel *k*-mers and are thus (putatively) associated with the same variant; sometimes this is abbreviated in the code using cc or CC, referring to the fact that these partitions are reflected as *connected components* in the *shared novel *k*-mers read graph; other abbreviations in the code include PART and CALLCLASS
- **augfastx**: sequences in Fasta or Fastq format, augmented with annotations indicating the position and abundance of interesting *k*-mers and mate sequences; the ``kevlar.parse_augmented_fastx`` and ``kevlar.print_augmented_fastx`` commands can be used to read and write data in augmented Fasta or Fastq format
- **contig**: in the context of kevlar, a contig almost always refers to a sequence assembled from a set of reads sharing novel *k*-mers and thus (putatively) spanning the same novel variant
- **reference cutout**: the algorithm kevlar uses to align contigs and call variants is not designed to map a short contig to a long chromosome sequence; therefore kevlar computes a "reference target sequence" or a "cutout" of the genome to which each contig is aligned for variant calling
