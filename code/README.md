# Implementation of Cryptographic Obfuscation

To run, you'll need sage (use 6.1.1 or later), openmp (should come with g++ by
default), and libgmp (use 6.0.0 or later).  The code has only been tested on
Linux, but should run fine on OS X.  Good luck to you Windows users.

To test, run
```
sage setup.py test
```

If this succeeds, you should now be able to test it as follows:
```
./obfuscator obf --test-circuit circuits/and.circ --secparam 8 --verbose
```

This should print out a bunch of stuff, and take no more than 5 seconds to run.
The important part is at the bottom, where with any luck you should see the word
"Pass" in green color.  If not, something is wrong; please e-mail the maintainer
to resolve this.

To use the Sahai-Zhandry scheme, use the `--sahai-zhandry` flag.  To use the
Zimmerman scheme, use the `--zimmerman` flag.  Using neither flag defaults to
the Ananth et al. scheme.  See `--help` for other options.
