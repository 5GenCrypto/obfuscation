#!/usr/bin/env bash

#
# Change these as needed
#
SAGE='sage'
PYTHON='python2'
CODE_DIR='code'
CIRCUIT_DIR='code/circuits'
LOG_DIR='runs'

mkdir -p $LOG_DIR

################################################################################

MIN=24
MAX=160
echo "* Varying the security parameter ($MIN -> $MAX)"

circuit='id.circ'

dir="$LOG_DIR/secparam.$circuit"
mkdir -p $dir

for secparam in `seq $MIN 8 $MAX`
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
        --eval 1 \
        --verbose 2>&1 | tee $dir/$circuit-$secparam-eval-time.log
    rm -rf $CIRCUIT_DIR/$circuit.obf.$secparam
done

################################################################################

MIN=8
MAX=40
echo "* Running point functions ($MIN -> $MAX)"

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
        --load-circuit $CIRCUIT_DIR/$circuit \
        --secparam $secparam \
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
done

################################################################################
