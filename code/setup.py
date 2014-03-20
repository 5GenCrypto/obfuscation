#!/usr/bin/env python2

from setuptools import setup, Extension, find_packages

__name__ = 'ind_obfuscation'
__author__ = 'Alex J. Malozemoff'
__version__ = '1.0a1.dev'

fastutils = Extension(
    'indobf.fastutils',
    libraries = ['gmp', 'gomp'],
    extra_compile_args = ['-fopenmp', '-ggdb', '-Wall'],
    sources = ['src/fastutils.c',
               'src/mpn_pylong.c',
               'src/mpz_pylong.c']
)

setup(name = __name__,
      version = __version__,
      author = __author__,
      description = 'Indistinguishability obfuscation',
      url = 'https://github.com/amaloz/ind-obfuscation',
      package_data = {'circuits': ['*.circ']},
      packages = ['indobf', 'circuits'],
      ext_modules = [fastutils],
      test_suite = 'unittests',
      # entry_points = {
      #     'console_scripts': ['indobf = indobf.run:main']
      # },
      classifiers = [
          'Topic :: Security :: Cryptography',
          'Environment :: Console',
          'Development Status :: 2 - Pre-Alpha',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: Free For Educational Use',
          'Programming Language :: C',
          'Programming Language :: Sage',
      ],
)
