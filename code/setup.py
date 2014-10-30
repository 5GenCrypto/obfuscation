#!/usr/bin/env python2

from setuptools import setup, Extension, find_packages

ATTACK = 0

__name__ = 'ind_obfuscation'
__author__ = 'Alex J. Malozemoff'
__version__ = '0.9'

libraries = [
    'gmp',
    'gomp',
]
compile_args = [
    '-fopenmp',
    '-O3',
    '-Wall',
]

if ATTACK:
    libraries.append('fplll')
    compile_args.append('-DATTACK')


obfuscator = Extension(
    'indobf._obfuscator',
    libraries=libraries,
    extra_compile_args=compile_args,
    sources=[
        'src/_obfuscator.cpp',
        'src/_zobfuscator.cpp',
        'src/mpn_pylong.cpp',
        'src/mpz_pylong.cpp',
        'src/pyutils.cpp',
        'src/utils.cpp',
    ]
)

setup(name=__name__,
      author=__author__,
      version=__version__,
      description='Cryptographic obfuscation implementation',
      package_data={'circuits': ['*.circ']},
      packages=['indobf', 'circuits'],
      ext_modules=[obfuscator],
      test_suite='t',
      classifiers=[
          'Topic :: Security :: Cryptography',
          'Environment :: Console',
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Operating System :: Unix',
          'License :: OSI Approved :: Free For Educational Use',
          'Programming Language :: C',
          'Programming Language :: Sage',
      ])
