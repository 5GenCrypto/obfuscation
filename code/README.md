# Implementation of Cryptographic Obfuscation

To run, you'll need python (use 2.7.X), sage (use 6.1.1 or later), openmp
(should come with g++ by default), libgmp (use 6.0.0 or later), and networkx
(use 1.8.1 or later) installed.  The code has only been tested on Linux, but
should run fine on OS X.  Good luck to you Windows users.

To build, run (only works for python2.7)

```
python setup.py test
```

If this succeeds, you should now be able to test it as follows:

```
sage obf/run.py obf --test-circuit circuits/and.circ --secparam 8 --verbose
```

This should print out a bunch of stuff, and take no more than 5 seconds to run.
The important part is at the bottom, where with any luck you should see the word
"Pass" in green color.  If not, something is wrong; please e-mail the maintainer
to resolve this.

There's code in the src/ directory which can be used to evaluate an existing
obfuscation, and is written entirely in C/C++ (so you do not need sage or python
installed).  <b>Note: this code is not up-to-date and may not work!</b> To use,
run "make" within the src/ directory and then run

```
./evaluate <obfuscation-directory> <input-as-a-bitstring>
```

So, for example, if you have an obfuscation in the directory obfuscation/, and
it takes 5-bit inputs, you can evaluate it on input 10101 as follows:

```
./evaluate obfuscation 10101
```

To enable the attack functionality, edit setup.py and set ATTACK = 1.  You will
need the fpLLL library (http://perso.ens-lyon.fr/damien.stehle/fplll/) for the
code to compile.  You can then run

```
sage obf/run.py obf --load-circuit circuits/and.circ --secparam 16 --attack --nslots 1
```

This obfuscates an AND circuit using only 1 slot of the plaintext space, and
then runs the attack specified in Section 3.1.2 of the paper.  It prints out two
values for g_1: one as generated during the obfuscation, and the other as
extracted from the attack.  Thus, they should both be equal if the attack
succeeds.

To use the Sahai-Zhandry scheme, use the `--sahai-zhandry` flag.  To use the
Zimmerman scheme, use the `--zimmerman` flag.  Using neither flag defaults to
the Ananth et al. scheme.
