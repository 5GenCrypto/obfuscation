#!/usr/bin/env bash

#
# Change these as needed
#
SAGE='sage'
PYTHON='python2'
CODE_DIR='../code'
CIRCUIT_DIR='../code/circuits'
LOG_DIR='runs'

# POINT=8

mkdir -p $LOG_DIR

# pushd $CIRCUIT_DIR
# ./point.py $POINT
# popd

secparam=24

circuit="fourandsnot.circ"
dir="$LOG_DIR/obliviate.$circuit.$secparam"
mkdir -p $dir


echo "* Running $circuit with security parameter $secparam (oblivious)"
$SAGE $CODE_DIR/indobf/run.py obf \
    --load-circuit $CIRCUIT_DIR/$circuit \
    --secparam $secparam \
    --obliviate \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-obv-obf-time.log
obf=$circuit.obf.$secparam
du --bytes $CIRCUIT_DIR/$obf/* \
    | tee $dir/$circuit-$secparam-obv-obf-size.log
eval='0000'
$SAGE $CODE_DIR/indobf/run.py obf \
    --load-obf $CIRCUIT_DIR/$obf \
    --eval $eval \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-obv-eval-time.log
rm -rf $CIRCUIT_DIR/$circuit.obf.$secparam

echo "* Running $circuit with security parameter $secparam (std)"
$SAGE $CODE_DIR/indobf/run.py obf \
    --load-circuit $CIRCUIT_DIR/$circuit \
    --secparam $secparam \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-obf-time.log
obf=$circuit.obf.$secparam
du --bytes $CIRCUIT_DIR/$obf/* \
    | tee $dir/$circuit-$secparam-obf-size.log
eval='0000'
$SAGE $CODE_DIR/indobf/run.py obf \
    --load-obf $CIRCUIT_DIR/$obf \
    --eval $eval \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-eval-time.log
rm -rf $CIRCUIT_DIR/$circuit.obf.$secparam
