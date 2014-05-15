Implementation of Indistinguishability Obfuscation
==================================================

This code presents an implementation of indistinguishability obfuscation, as
appearing in the work:

"Implementing Provably Secure Program Obfuscation." Unpublished. 2014.

The code contains two implementations, one based on Barrington's theorem and the
other based on an approach by Sauerhoff et al. (STACS, 1999).  We follow the
general outlines of existing indistinguishability obfuscation literature (e.g.,
Ananth et al, ePrint 2014; Barak et al, EUROCRYPT 2014; Brakerski and Rothblum,
TCC 2014; Garg et al., FOCS 2013; etc.) using a graded encoding scheme based on
the code of Coron et al. (CRYPTO, 2013).

The implementation is in a mix of Python and C, using Sage
(http://sagemath.org), GNU GMP (https://gmplib.org) and OpenMP
(http://openmp.org).  The current code is still in prototype form, although it
is close to fully operational.

Instructions for building the code are in the code/ directory, and scripts for
running experiments and processing the results are in the scripts/ directory.
Full traces of all the experiments found in the above paper can be found in the
scripts/results/ directory.
