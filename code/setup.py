#!/usr/bin/env python2

from setuptools import setup, Extension

__version__ = '1.0'

fastutils = Extension(
    'fastutils',
    libraries = ['gmp', 'gomp'],
    extra_compile_args = ['-fopenmp', '-ggdb', '-Wall'],
    sources = ['src/fastutils.c',
               'src/mpn_pylong.c',
               'src/mpz_pylong.c']
)

setup(name = 'IndObfuscation',
      version = __version__,
      description = 'Indistinguishability obfuscation',
      author = 'Alex J. Malozemoff',
      packages = ['indobf', 'tests'],
      test_suite = "tests.fastutils_unittest",
      ext_modules = [fastutils])
