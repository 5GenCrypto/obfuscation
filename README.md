Implementation of Cryptographic Obfuscation
==================================================

This code presents an implementation of cryptographic obfuscation, as appearing
in the work:

"Implementing Cryptographic Program Obfuscation." Daniel Apon, Yan Huang,
Jonathan Katz, Alex J. Malozemoff. Cryptology ePrint Archive 2014/779.
https://eprint.iacr.org/2014/779.

We follow the general outlines of existing indistinguishability obfuscation
literature, mainly the work of Ananth, Gupta, Ishai, and Sahai (CCS 2014:
https://eprint.iacr.org/2014/222), using a graded encoding scheme based on the
code of Coron et al. (CRYPTO, 2013).

The implementation is in a mix of Python and C, using Sage
(http://sagemath.org), GNU GMP (https://gmplib.org) and OpenMP
(http://openmp.org).  The current code is still in development, but the master
branch should work most of the time.

<b>NEW</b>: We now have an (alpha) implementation of the Sahai-Zhandry scheme
(ePrint Cryptology Archive, 2014/773: https://eprint.iacr.org/2014/773)!

Instructions for building and running the code are in code/README.md, and
scripts for running experiments and processing the results are in the scripts/
directory.  Full traces of all the experiments found in the above paper can be
found in the scripts/results/ directory.

For any questions/comments, please e-mail amaloz at cs dot umd dot edu.
