#!/usr/bin/env python
# Setup file to compile tanimoto module into shared library for python

import os
from distutils.core import setup, Extension

FASTTANI_PATH = "/raid5/people/mysinger/code/tanimoto/trunk"
FASTTANI_C = os.path.join(FASTTANI_PATH, "fast_tanimoto.c")

local_path = os.path.dirname(os.path.abspath(__file__))
TANIMOTO_C = os.path.join(local_path, 'tanimoto.c')

setup (name='tanimoto',
       version='September-2011',
       author='Michael Mysinger',
       ext_modules = [Extension('tanimoto', [TANIMOTO_C, FASTTANI_C],
                     include_dirs = [FASTTANI_PATH],
                     library_dirs = [],
                     libraries = [])],
       )
