#!/usr/bin/env python2

from setuptools import setup, Extension

library_dirs = [
    'src/.libs'
]

libraries = [
    'obf',
]
compile_args = [
    '-O3',
    '-Wall',
    '-Wextra',
    '-Isrc/',
    '-Lsrc/'
]

obfuscator = Extension(
    'pyobf._obfuscator',
    library_dirs=library_dirs,
    libraries=libraries,
    extra_compile_args=compile_args,
    sources=[
        'pywrapper/obfuscator_wrapper.cpp',
        'pywrapper/pyutils.cpp',
    ]
)

zobfuscator = Extension(
    'pyobf._zobfuscator',
    library_dirs=library_dirs,
    libraries=libraries,
    extra_compile_args=compile_args,
    sources=[
        'pywrapper/zobfuscator_wrapper.cpp',
        'pywrapper/pyutils.cpp',
    ]
)

setup(name='obfuscator',
      author='Alex J. Malozemoff',
      author_email='amaloz@cs.umd.edu',
      version='0.2a0',
      description='Implementation of cryptographic program obfuscation',
      license='GPLv2',
      url='https://github.com/amaloz/obfuscation',
      packages=['pyobf'],
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
