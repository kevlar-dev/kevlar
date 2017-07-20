#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from setuptools import setup, Extension
import versioneer


align = Extension(
    'kevlar.alignment',
    sources=[
        'kevlar/alignment.c', 'third-party/ksw2/ksw2_extz.c', 'src/align.c'
    ],
    include_dirs=['inc/', 'third-party/ksw2/'],
    language='c',
)


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
      ext_modules=[align],
      setup_requires=['pysam', 'networkx', 'pandas'],
      install_requires=['pysam', 'networkx', 'pandas'],
      entry_points={
          'console_scripts': ['kevlar = kevlar.__main__:main']
      },
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      zip_safe=True)
