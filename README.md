# Implementation of Cryptographic Program Obfuscation

This code presents several implementations of cryptographic program obfuscation,
as listed in the following table:

Scheme | Authors | ePrint Reference
------ | ------- | ----------------
SZ     | Sahai, Zhandry | [2014/773](https://eprint.iacr.org/2014/773)
Z      | Zimmerman | [2014/776](https://eprint.iacr.org/2014/776)

A discussion of the AGIS implementation (since removed due to SZ being better
across the board) appears in the work:

"Implementing Cryptographic Program Obfuscation." Daniel Apon, Yan Huang,
Jonathan Katz, Alex J. Malozemoff. Cryptology ePrint Archive 2014/779.
https://eprint.iacr.org/2014/779.

All schemes use the graded encoding scheme of Coron et al. (CRYPTO,
2013. https://eprint.iacr.org/2013/183).  The implementation is in a mix of
Python and C, using [Sage](http://sagemath.org), [GNU GMP](https://gmplib.org)
and [OpenMP](http://openmp.org).

Instructions for building and running the code are in code/README.md, and
scripts for running experiments and processing the results are in the scripts/
directory.

For any questions/comments, please e-mail amaloz at cs dot umd dot edu.
