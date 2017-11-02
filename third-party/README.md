The `ksw2` directory contains two files from Heng Li's [ksw2](https://github.com/lh3/ksw2) repository: a header file and a file implementing one of the library's main alignment strategies.
The `KSW2_VERSION` file was added by me for sake of record keeping.

The `fermi-lite` directory contains the entire source code distribution of Heng Li's [fermi-lite](https://github.com/lh3/fermi-lite) library.
Compiling this successfully within the Python ecosystem required that the line `const uint8_t rle_auxtab[8];` be dropped from `rle.h`.
The `FERMILITE_VERSION` file was added by me for sake of record keeping.

After further reflection, maybe I should rename this directory from `third-party/` to `heng-li-libs/`. :-)
