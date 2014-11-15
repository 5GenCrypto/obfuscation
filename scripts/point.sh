#!/usr/bin/env bash

#
# Change these as needed
#
SAGE='sage'
PYTHON='python2'
CODE_DIR='../code'
CIRCUIT_DIR='../code/circuits'
LOG_DIR='runs'
SCHEME=''

mkdir -p $LOG_DIR

MIN=8
MAX=16
INC=4

echo "* Running point functions ($MIN -> $MAX)"

pushd $CIRCUIT_DIR
for point in `seq $MIN $INC $MAX`
do
    ./point.py $point
done
popd

secparam=52

dir="$LOG_DIR/point.$secparam"
mkdir -p $dir

for point in `seq $MIN $INC $MAX`
do
    circuit="point-$point.circ"
    echo "* Running $circuit with security parameter $secparam"
    $SAGE $CODE_DIR/obfuscator obf \
          --load-circuit $CIRCUIT_DIR/$circuit \
          --secparam $secparam \
          $SCHEME \
          --verbose 2>&1 | tee $dir/$circuit-$secparam-obf-time.log
    obf=$circuit.obf.$secparam
    du --bytes $CIRCUIT_DIR/$obf/* \
        | tee $dir/$circuit-$secparam-obf-size.log
    eval=`sed -n 1p $CIRCUIT_DIR/$circuit | awk '{ print $3 }'`
    echo $eval
    $SAGE $CODE_DIR/obfuscator obf \
          --load-obf $CIRCUIT_DIR/$obf \
          --eval $eval \
          $SCHEME \
          --verbose 2>&1 | tee $dir/$circuit-$secparam-eval-time.log
    rm -rf $CIRCUIT_DIR/$circuit.obf.$secparam
done
