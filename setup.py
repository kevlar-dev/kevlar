#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2016 The Regents of the University of California
#
# This file is part of kevlar (http://github.com/standage/kevlar) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import setuptools
import versioneer


setuptools.setup(name='kevlar',
                 version=versioneer.get_version(),
                 cmdclass=versioneer.get_cmdclass(),
                 description=('Reference-free variant discovery in large '
                              'eukaryotic genomes'),
                 url='http://github.com/standage/kevlar',
                 author='Daniel Standage',
                 author_email='daniel.standage@gmail.com',
                 license='MIT',
                 packages=['kevlar'],
                 scripts=['bin/kevlar'],
                 #install_requires=['khmer>=2.1'],
                 classifiers=[
                    'Development Status :: 4 - Beta',
                    'Environment :: Console',
                    'License :: OSI Approved :: MIT License',
                    'Programming Language :: Python :: 2.7',
                    'Programming Language :: Python :: 3.5',
                    'Topic :: Scientific/Engineering :: Bio-Informatics'
                 ],
                 zip_safe=True)
