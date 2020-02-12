#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import, print_function

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup

setup(
  name = 'CSTB_core',
  version = '1.0',
  license='BSD',
  description = 'Utility function for io and word searching and indexing',
  author = 'Guillaume Launay, Cecile Hilpert',
  author_email = 'pitooon@gmail.com, cecile.hilpert@gmail.com',
  url = 'https://github.com/MMSB-MOBI/CSTB_core', # use the URL to the github repo
  packages=find_packages('src'),
  package_dir={'': 'src'},
  include_package_data=True,
  zip_safe=False,
  py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
  download_url = 'https://github.com/MMSB-MOBI/CSTB_core/tarball/0.1', # I'll explain this in a second
  keywords = ['crispr', 'webserver'], # arbitrary keywords
  classifiers = [],
  install_requires=[
          'docopt', 'biopython'
      ]
  #data_files=[
  #          ('external', ['external/pathos.tar.bz']),
  #          ('bin', ['bin/module1.py']),
  #          ('conf',['conf/confModule1.json'])
  #    ]
  #  dependency_links = [
  #      "http://dev.danse.us/trac/pathos"
  #  ]
)
