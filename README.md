# Cryptographic Program Obfuscation (v0.2a)

This project provides the following implementations of cryptographic program
obfuscation:

Scheme | Authors | ePrint Reference
------ | ------- | ----------------
SZ     | Sahai, Zhandry | [2014/773](https://eprint.iacr.org/2014/773)
Z      | Zimmerman | [2014/776](https://eprint.iacr.org/2014/776)

A discussion of the AGIS implementation (since removed due to SZ being better
across the board) appears in the work:

"Implementing Cryptographic Program Obfuscation." Daniel Apon, Yan Huang,
Jonathan Katz, Alex J. Malozemoff. Cryptology ePrint Archive 2014/779.
https://eprint.iacr.org/2014/779.

The implementation uses [libmmap](https://github.com/amaloz/gghlite-flint) to
support both the CLT13 and GGHlite graded encoding schemes for SZ; Z only works
with CLT13.

## Building

Run the following:

```
autoreconf -i
./configure
make
sudo make install
```

This installs the underlying obfuscation library `libobf` to your system.  To
install the python front-end, proceed as follows:

```
cd python
python2 setup.py test
```

This runs a bunch of test, all of which should hopefully pass.

You can then run the obfuscator by running
```
./obfuscator obf --test-circuit circuits/and.circ --secparam 16 -v
```

## Contact

For any questions/comments, please e-mail amaloz at cs dot umd dot edu.
