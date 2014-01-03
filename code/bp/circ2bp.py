#!/usr/bin/env python3
#
# Converts a circuit with AND and NOT gates to branching program.
#

import numpy as np
import sys


I = np.eye(5)
# alpha = (02143)
A = np.matrix('0 0 1 0 0; 0 0 0 0 1; 0 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0')
# beta = (01342)
B = np.matrix('0 1 0 0 0; 0 0 0 1 0; 1 0 0 0 0; 0 0 0 0 1; 0 0 1 0 0')
# commutator = (01234)
C = np.matrix('0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 1 0 0 0 0')
# R conjugates the commutator to A
R = np.matrix('1 0 0 0 0; 0 0 1 0 0; 0 1 0 0 0; 0 0 0 0 1; 0 0 0 1 0')
Ri = np.linalg.inv(R)
# S conjugates the commutator to B
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
                        print("error: only support NOT so far:", line.strip())
                        exit(-1)
                elif arity == 2:
                    _, a, b, c, d, _, _, _, in1, in2, \
                        rest = rest.split(maxsplit=10)
                    a, b, c, d, in1, in2 = ints(a, b, c, d, in1, in2)
                    if a == 0 and b == 0 and c == 0 and d == 1:
                        # AND gate
                        if ms[in1].startswith('C:'):
                            _, idx = ms[in1].split(':')
                            s1 = "A:%s" % idx
                        else:
                            s1 = "Ri %s R" % ms[in1]
                        if ms[in2].startswith('C:'):
                            _, idx = ms[in2].split(':')
                            s2 = "B:%s" % idx
                        else:
                            s2 = "Si %s S" % ms[in2]
                        str = "%s %s inv(%s) inv(%s)" % (s1, s2, s1, s2)
                        ms.append(str)
                    else:
                        print("error: only support AND so far:", line.strip())
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
        if m.startswith('inv'):
            str = m[m.find('(')+1:m.rfind(')')]
            comp = comp * np.linalg.inv(_eval_bp(str, inp))
        elif m.find(':') != -1:
            matrix, idx = m.split(':')
            comp = comp * (I if inp[int(idx)] == '0' else eval(matrix))
        else:
            comp = comp * eval(m)
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
