# Cryptographic Program Obfuscation

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
./obfuscator obf --test circuits/and.circ --secparam 16 -v
```

## Contact

For any questions/comments, please e-mail amaloz at galois dot com.
