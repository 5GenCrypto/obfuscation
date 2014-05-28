#!/usr/bin/env python2

from setuptools import setup, Extension, find_packages

__name__ = 'ind_obfuscation'
__author__ = 'Alex J. Malozemoff'
__version__ = '0.9'

obfuscator = Extension(
    'indobf._obfuscator',
    libraries = [
        'gmp',
        'gomp',
    ],
    extra_compile_args = [
        '-fopenmp',
        '-O3',
        '-Wall',
    ],
    sources = [
        'src/_obfuscator.cpp',
        'src/mpn_pylong.cpp',
        'src/mpz_pylong.cpp',
        'src/utils.cpp',
    ]
)

setup(name = __name__,
      author = __author__,
      version = __version__,
      description = 'Indistinguishability obfuscation implementation',
      package_data = {'circuits': ['*.circ']},
      packages = ['indobf', 'circuits'],
      ext_modules = [obfuscator],
      test_suite = 't',
      classifiers = [
          'Topic :: Security :: Cryptography',
          'Environment :: Console',
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Operating System :: Unix',
          'License :: OSI Approved :: Free For Educational Use',
          'Programming Language :: C',
          'Programming Language :: Sage',
      ],
)
