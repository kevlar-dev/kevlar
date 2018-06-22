#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from setuptools import setup, Extension
import glob
import versioneer


ksw2 = Extension(
    'kevlar.alignment',
    sources=[
        'kevlar/alignment.c', 'third-party/ksw2/ksw2_extz.c', 'src/align.c'
    ],
    include_dirs=['inc/', 'third-party/ksw2/'],
    language='c',
)

fermilite = Extension(
    'kevlar.assembly',
    sources=['kevlar/assembly.c'] + glob.glob('third-party/fermi-lite/*.c'),
    include_dirs=['third-party/fermi-lite/'],
    extra_link_args=['-lz'],
    language='c',
)

sequencemod = Extension(
    'kevlar.sequence',
    sources=['kevlar/sequence.c'],
    language='c'
)

dependencies = [
    'pysam>=0.14', 'networkx>=2.0', 'pandas>=0.23', 'scipy>=1.1',
    'matplotlib>=2.2'
]

setup(name='biokevlar',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description=('Reference-free variant discovery scalable to large '
                   'eukaryotic genomes'),
      url='https://github.com/dib-lab/kevlar',
      author='Daniel Standage',
      author_email='daniel.standage@gmail.com',
      license='MIT',
      packages=['kevlar', 'kevlar.cli', 'kevlar.tests'],
      package_data={
          'kevlar': ['kevlar/tests/data/*', 'kevlar/tests/data/*/*']
      },
      include_package_data=True,
      ext_modules=[ksw2, fermilite, sequencemod],
      setup_requires=dependencies,
      install_requires=dependencies,
      entry_points={
          'console_scripts': ['kevlar = kevlar.__main__:main']
      },
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      zip_safe=True)
