#!/usr/bin/env bash

#
# Change these as needed
#
SAGE='sage'
PYTHON='python2'
CODE_DIR='../code'
CIRCUIT_DIR='../code/circuits'
LOG_DIR='runs'

mkdir -p $LOG_DIR

if [[ "$1" = "" || "$2" = "" ]]
then
    echo "Usage: secparam.sh <circuit-name> <input-length>"
    exit 1
fi
circuit=$1
eval=`$PYTHON -c "print('0' * $2)"`

MIN=24
INC=16
MAX=160
echo "* Running circuit $circuit"
echo "* Varying the security parameter ($MIN -> $MAX)"

dir="$LOG_DIR/secparam.$circuit"
mkdir -p $dir

for secparam in `seq $MIN $INC $MAX`
do
    echo "* Running $circuit with security parameter $secparam"
    $SAGE $CODE_DIR/indobf/run.py obf \
        --load-circuit $CIRCUIT_DIR/$circuit \
        --secparam $secparam \
        --verbose 2>&1 | tee $dir/$circuit-$secparam-obf-time.log
    obf=$circuit.obf.$secparam
    du --bytes $CIRCUIT_DIR/$obf/* \
        | tee $dir/$circuit-$secparam-obf-size.log
    $SAGE $CODE_DIR/indobf/run.py obf \
        --load-obf $CIRCUIT_DIR/$obf \
        --eval $eval \
        --verbose 2>&1 | tee $dir/$circuit-$secparam-eval-time.log
    rm -rf $CIRCUIT_DIR/$circuit.obf.$secparam
done
