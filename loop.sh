#!/usr/bin/env bash

RESULT=0
while [ "$RESULT" -eq "0" ]; do
    ./obfuscator obf --test circuits/point-4.circ -v --nthreads 1 --ncores 1 --mmap CLT --secparam 8
    RESULT=$?
done
