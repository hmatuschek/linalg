#!/usr/bin/env python

from distutils.core import setup, Extension

from numpy.distutils.system_info import numpy_info
numpy_info = numpy_info();
include_dirs = numpy_info.get_include_dirs() + ['../src/'];

setup(name="Linalg",
      version="0.1",
      description="Python/NumPy interface to the Linalg C++ BLAS and LAPACK interface.",
      author="Hannes Matuschek",
      author_email="hmatuschek at gmail dot com",
      url="https://hmatuschek.github.com/linalg",
      ext_modules=[Extension('linalg', ['src/module.cc'],
                             include_dirs=include_dirs,
                             libraries=['blas', 'lapack', 'm'])])
