#!/usr/bin/env sage

from setuptools import setup, Extension

libraries = [
    'gmp',
    'gomp',
]
compile_args = [
    '-fopenmp',
    '-O3',
    '-Wall',
    '-pthread',
    # '-DTHPOOL_DEBUG',
]

zobfuscator = Extension(
    'obf._zobfuscator',
    libraries=libraries,
    extra_compile_args=compile_args,
    sources=[
        'src/circuit.cpp',
        'src/clt_mlm.cpp',
        'src/_zobfuscator.cpp',
        'src/mpn_pylong.cpp',
        'src/mpz_pylong.cpp',
        'src/pyutils.cpp',
        'src/thpool.cpp',
        'src/thpool_fns.cpp',
        'src/utils.cpp',
    ]
)

obfuscator = Extension(
    'obf._obfuscator',
    libraries=libraries,
    extra_compile_args=compile_args,
    sources=[
        'src/clt_mlm.cpp',
        'src/_obfuscator.cpp',
        'src/mpn_pylong.cpp',
        'src/mpz_pylong.cpp',
        'src/pyutils.cpp',
        'src/thpool.cpp',
        'src/thpool_fns.cpp',
        'src/utils.cpp',
    ]
)

setup(name='obfuscator',
      author='Alex J. Malozemoff',
      author_email='amaloz@cs.umd.edu',
      version='1.0a0',
      description='Implementation of cryptographic program obfuscation',
      license='GPLv2',
      url='https://amaloz.github.io/obfuscation',
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
