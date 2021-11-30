# Setup file to compile tanimoto module into shared library for python
import os
from distutils.core import setup, Extension

def dt_subdir(dirname):
    return os.path.join(os.environ["DY_ROOT"], dirname)

setup (name='fastdl',
       version='August-2006',
       author='Michael Mysinger',
       ext_modules = [Extension('fastdl', ['fast_daylight.c'],
                     include_dirs = [dt_subdir("include")],
                     library_dirs = [dt_subdir("lib")],
                     libraries = ["dt_smiles", "dt_finger"])],
       )
