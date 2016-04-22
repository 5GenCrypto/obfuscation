#!/usr/bin/env python2

from setuptools import setup, Extension

library_dirs = [
    '../src/.libs'
]

libraries = [
    'obf',
]
compile_args = [
    # XXX: Setting this to anything higher than 0 causes the program to crash
    # when reading in the directory string in obf_evaluate_wrapper.
    # Is this a python bug?
    '-O0',
    '-Wall',
    '-I../src/',
    '-L../src/'
]

zobfuscator = Extension(
    'obf._zobfuscator',
    library_dirs=library_dirs,
    libraries=libraries,
    extra_compile_args=compile_args,
    sources=[
        'src/zobfuscator_wrapper.cpp',
        'src/mpn_pylong.cpp',
        'src/mpz_pylong.cpp',
        'src/pyutils.cpp',
    ]
)

obfuscator = Extension(
    'obf._obfuscator',
    library_dirs=library_dirs,
    libraries=libraries,
    extra_compile_args=compile_args,
    sources=[
        'src/obfuscator_wrapper.cpp',
        'src/mpn_pylong.cpp',
        'src/mpz_pylong.cpp',
        'src/pyutils.cpp',
    ]
)

setup(name='obfuscator',
      author='Alex J. Malozemoff',
      author_email='amaloz@cs.umd.edu',
      version='0.2a0',
      description='Implementation of cryptographic program obfuscation',
      license='GPLv2',
      url='https://github.com/amaloz/obfuscation',
      packages=['obf'],
      ext_modules=[obfuscator, zobfuscator],
      scripts=['obfuscator'],
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
