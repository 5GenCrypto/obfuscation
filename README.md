# Implementation of Cryptographic Program Obfuscation

This code presents several implementations of cryptographic program obfuscation,
as listed in the following table:

Scheme | ePrint Reference | Status
------ | ---------------- | ------
Ananth, Gupta, Ishai, Sahai (AGIS) | [2014/222](https://eprint.iacr.org/2014/222) |
Sahai, Zhandry (SZ) | [2014/773](https://eprint.iacr.org/2014/773) |
Zimmerman (Z) | [2014/776](https://eprint.iacr.org/2014/776) | not working

For some background on cryptographic program obfuscation, see
https://eprint.iacr.org/2013/451.

A discussion of the AGIS implementation appears in the work:

"Implementing Cryptographic Program Obfuscation." Daniel Apon, Yan Huang,
Jonathan Katz, Alex J. Malozemoff. Cryptology ePrint Archive 2014/779.
https://eprint.iacr.org/2014/779.

All schemes use the graded encoding scheme of Coron et al. (CRYPTO,
2013. https://eprint.iacr.org/2013/183).  The implementation is in a mix of
Python and C, using [Sage](http://sagemath.org), [GNU GMP](https://gmplib.org)
and [OpenMP](http://openmp.org).  <b>The code is under active development, and
may not be stable!</b> I will try to keep the master branch in a working state
(most of the time), though.

Instructions for building and running the code are in code/README.md, and
scripts for running experiments and processing the results are in the scripts/
directory.  Full traces of all the experiments found in the above paper can be
found in the scripts/results/ directory.

For any questions/comments, please e-mail amaloz at cs dot umd dot edu.

## Challenges

In order to aid in the cryptanalysis of cryptographic program obfuscation, we
plan to release a series of obfuscated point functions (that is, functions that
output 0 on all inputs except one), with the goal of the attacker being to learn
the hidden point.  The following table lists the available/broken challenges.

Challenge | Link | Approach | Date Broken | Broken By
--------- | ---- | -------- | ----------- | ---------
14-bit point function; 60-bit security parameter | [Link](https://www.dropbox.com/s/85d03o0ny3b1c0c/point-14.circ.obf.60.zip) (23.96) | AGIS | 18 Oct 2014 | Daniel J. Bernstein, Andreas Huelsing, Tanja Lange, Ruben Niederhagen

