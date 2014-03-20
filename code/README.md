Implementation of Indistinguishability Obfuscation
==================================================

To run, you'll need sage, openmp, and libgmp installed.

To build, run (only works for python2.7)

> python setup.py test

If this succeeds, you should now be able to cd into `build/lib` and run

> sage indobf/run.py obf --test-circuit circuits/and.circ --verbose
