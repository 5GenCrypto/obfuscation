#!/usr/bin/env python2

from setuptools import setup, Extension, find_packages

__name__ = 'ind_obfuscation'
__version__ = '1.0'

fastutils = Extension(
    'fastutils',
    libraries = ['gmp', 'gomp'],
    extra_compile_args = ['-fopenmp', '-ggdb', '-Wall'],
    sources = ['src/fastutils.c',
               'src/mpn_pylong.c',
               'src/mpz_pylong.c']
)

setup(name = __name__,
      version = __version__,
      author = 'Alex J. Malozemoff',
      description = 'Indistinguishability obfuscation',
      url = 'https://github.com/amaloz/ind-obfuscation',
      packages = find_packages(),
      ext_modules = [fastutils],
      test_suite = 'unittests',
      classifiers = [
          'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'Frakework :: Sage',
          'Topic :: Security :: Cryptography',
          'License :: OSI Approved :: Free For Educational Use',
      ],
)
