#!/usr/bin/env bash

SAGE='/home/amaloz/Desktop/sage-6.1.1-x86_64-Linux/sage --python'
CODE_DIR='code'
CIRCUIT_DIR='code/circuits'
LOG_DIR='runs'

mkdir -p $LOG_DIR

################################################################################

MIN=24
MAX=24
echo "* Varying the security parameter ($MIN -> $MAX)"

circuit='id.circ'

dir="$LOG_DIR/secparam.$circuit"
mkdir -p $dir

for secparam in `seq $MIN 8 $MAX`
do
    echo "* Running $circuit with security parameter $secparam"
    $SAGE $CODE_DIR/indobf/run.py obf \
        --test-circuit $CIRCUIT_DIR/$circuit \
        --secparam $secparam \
        --verbose 2>&1 | tee $dir/$circuit-$secparam-time.log
    du --bytes $CIRCUIT_DIR/$circuit.obf.$secparam/* \
        | tee $dir/$circuit-$secparam-size.log
    rm -rf $CIRCUIT_DIR/$circuit.obf.$secparam
done

################################################################################

MIN=8
MAX=40
echo "* Running point functions ($MIN -> MAX)"

pushd $CIRCUIT_DIR
for point in `seq $MIN 8 $MAX`
do
    ./point.py $point
done
popd

secparam=24

dir="$LOG_DIR/point.$secparam"
mkdir -p $dir

for point in `seq $MIN 8 $MAX`
do
    circuit="point-$point.circ"
    echo "* Running $circuit with security parameter $secparam"
    $SAGE $CODE_DIR/indobf/run.py obf \
        --test-circuit $CIRCUIT_DIR/$circuit \
        --secparam $secparam \
        --verbose 2>&1 | tee $dir/$circuit-$secparam-time.log
    du --bytes $CIRCUIT_DIR/$circuit.obf.$secparam/* \
        | tee $dir/$circuit-$secparam-size.log
    rm -rf $CIRCUIT_DIR/$circuit.obf.$secparam
done

################################################################################
