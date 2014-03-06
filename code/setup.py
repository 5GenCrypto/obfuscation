#!/usr/bin/env python2

from setuptools import setup, Extension

fastutils = Extension('fastutils',
                      libraries = ['gmp', 'gomp'],
                      extra_compile_args=['-fopenmp', '-g', '-Wall'],
                      sources = ['src/fastutils.c',
                                 'src/mpn_pylong.c',
                                 'src/mpz_pylong.c'])

setup(name = 'FastUtils',
      version = '1.0',
      description = 'Fast utilities',
      author = 'Alex J. Malozemoff',
      test_suite = "test.fastutils_unittest",
      ext_modules = [fastutils])
