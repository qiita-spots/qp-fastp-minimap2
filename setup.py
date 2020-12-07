#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2020, Qiita development team.
#
# Distributed under the terms of the BSD 3-clause License License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from setuptools import setup
from glob import glob

__version__ = "0.1.0-dev"

classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: Libraries :: Python Modules
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""

with open('README.rst') as f:
    long_description = f.read()

classifiers = [s.strip() for s in classes.split('\n') if s]

setup(name='qp-fastp-minimap2 Qiita Plugin',
      version=__version__,
      long_description=long_description,
      license="BSD",
      description='Qiita Type Plugin: qp-fastp-minimap2',
      author="Qiita development team",
      author_email="qiita.help@gmail.com",
      url='https://github.com/qiita-spots/qp-fastp-minimap2',
      test_suite='nose.collector',
      packages=['qp_fastp_minimap2'],
      package_data={'qp_fastp_minimap2': ['support_files/config_file.cfg']},
      scripts=glob('scripts/*'),
      extras_require={'test': ["nose >= 0.10.1", "pep8"]},
      install_requires=['click >= 3.3', 'qiita_client'],
      classifiers=classifiers
      )
