#!/usr/bin/env python2

from setuptools import setup, Extension

fastutils = Extension('fastutils',
                      libraries = ['gmp'],
                      extra_compile_args=['-fopenmp'],
                      extra_link_args=['-lgomp'],
                      sources = ['src/fastutils.c'])

setup(name = 'FastUtils',
      version = '1.0',
      description = 'Fast utilities',
      author = 'Alex J. Malozemoff',
      test_suite = "test.fastutils_unittest",
      ext_modules = [fastutils])
