#!/usr/bin/env python3
#
# Converts a circuit with AND and NOT gates to branching program.
#

import numpy as np
import sys


I = np.eye(5)
# The commutator (01234)
C = np.matrix('0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 1 0 0 0 0')
# R conjugates the commutator to alpha = (02143)
R = np.matrix('1 0 0 0 0; 0 0 1 0 0; 0 1 0 0 0; 0 0 0 0 1; 0 0 0 1 0')
Ri = np.linalg.inv(R)
# S conjugates the commutator to beta = (01342)
S = np.matrix('1 0 0 0 0; 0 1 0 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 1 0 0')
Si = np.linalg.inv(S)
# T conjugates the commutator to its inverse (43210)
T = np.matrix('0 0 0 0 1; 0 0 0 1 0; 0 0 1 0 0; 0 1 0 0 0; 1 0 0 0 0')
Ti = np.linalg.inv(T)


def ints(*args):
    return (int(arg) for arg in args)


def circ2bp(fname):
    ms = []
    with open(fname) as f:
        for line in f:
            num, rest = line.split(maxsplit=1)
            num = int(num)
            if rest.startswith('input'):
                ms.append("C:%d" % num)
            elif rest.startswith('gate') or rest.startswith('output'):
                if rest.startswith('gate'):
                    _, _, arity, _, rest = rest.split(maxsplit=4)
                else:
                    _, _, _, arity, _, rest = rest.split(maxsplit=5)
                arity = int(arity)
                if arity == 1:
                    _, a, b, _, _, _, in1, rest = rest.split(maxsplit=7)
                    a, b, in1 = ints(a, b, in1)
                    if a == 1 and b == 0:
                        # NOT gate
                        str = "Ti %s T C" % ms[in1]
                        ms.append(str)
                    else:
                        print("error: only support NOT so far")
                        exit(-1)
                elif arity == 2:
                    _, a, b, c, d, _, _, _, in1, in2, \
                        rest = rest.split(maxsplit=10)
                    a, b, c, d, in1, in2 = ints(a, b, c, d, in1, in2)
                    if a == 0 and b == 0 and c == 0 and d == 1:
                        # AND gate
                        str = "Ri %s R Si %s S inv(Ri %s R) inv(Si %s S)" \
                              % (ms[in1], ms[in2], ms[in1], ms[in2])
                        ms.append(str)
                    else:
                        print("error: only support AND so far")
                        exit(-1)
                else:
                    print("error: arity %d unsupported" % arity)
            else:
                print("error: unknown type")
                exit(-1)
    print(ms[-1])
    return ms[-1]


def splitme(str):
    parencount = 0
    newstr = ""
    for i in range(len(str)):
        if str[i] == ')':
            parencount -= 1
            newstr += str[i]
        elif str[i] == ' ' and parencount == 0:
            newstr += "%"
        elif str[i] == '(':
            parencount += 1
            newstr += str[i]
        else:
            newstr += str[i]
    return newstr.split('%')


def _eval_bp(bp, inp):
    ms = splitme(bp)
    comp = I
    for m in ms:
        if m == 'C':
            comp = comp * C
        elif m == 'R':
            comp = comp * R
        elif m == 'Ri':
            comp = comp * Ri
        elif m == 'S':
            comp = comp * S
        elif m == 'Si':
            comp = comp * Si
        elif m == 'T':
            comp = comp * T
        elif m == 'Ti':
            comp = comp * Ti
        elif m.startswith('inv'):
            str = m[m.find('(')+1:m.rfind(')')]
            comp = comp * np.linalg.inv(_eval_bp(str, inp))
        elif m.startswith('C:'):
            _, idx = m.split(':')
            idx = int(idx)
            comp = comp * (I if inp[idx] == '0' else C)
        else:
            print("error: unknown m %s" % m)
            exit(-1)
    return comp


def eval_bp(bp, inp):
    out = _eval_bp(bp, inp)
    if (out == np.eye(5)).all():
        return 0
    elif (out == C).all():
        return 1
    else:
        print("error: invalid return matrix %s" % out)
        exit(-1)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("error: invalid arguments")
        exit(-1)

    circuit = sys.argv[1]
    inp = sys.argv[2]
    bp = circ2bp(circuit)
    out = eval_bp(bp, inp)
    print("OUTPUT = %d" % out)
