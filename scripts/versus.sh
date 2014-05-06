#!/usr/bin/env bash

#
# Change these as needed
#
SAGE='sage'
PYTHON='python2'
CODE_DIR='../code'
CIRCUIT_DIR='../code/circuits'
LOG_DIR='runs'

POINT=8

mkdir -p $LOG_DIR

pushd $CIRCUIT_DIR
./point.py $POINT
popd

secparam=48

dir="$LOG_DIR/versus.$POINT.$secparam"
mkdir -p $dir

circuit="point-$point.circ"
echo "* Running $circuit with security parameter $secparam"
$SAGE $CODE_DIR/indobf/run.py obf \
    --load-circuit $CIRCUIT_DIR/$circuit \
    --secparam $secparam \
    --old \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-obf-time.log
obf=$circuit.obf.$secparam
du --bytes $CIRCUIT_DIR/$obf/* \
    | tee $dir/$circuit-$secparam-obf-size.log
eval=`$PYTHON -c "print('0' * $point)"`
$SAGE $CODE_DIR/indobf/run.py obf \
    --load-obf $CIRCUIT_DIR/$obf \
    --eval $eval \
    --verbose 2>&1 | tee $dir/$circuit-$secparam-eval-time.log
rm -rf $CIRCUIT_DIR/$circuit.obf.$secparam
