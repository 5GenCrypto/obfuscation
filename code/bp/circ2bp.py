
import numpy as np
import sys


C = np.matrix('0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 1 0 0 0 0')
Cinv = np.linalg.inv(C)
R = np.matrix('1 0 0 0 0; 0 0 1 0 0; 0 1 0 0 0; 0 0 0 0 1; 0 0 0 1 0')
Rinv = np.linalg.inv(R)
S = np.matrix('1 0 0 0 0; 0 1 0 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 1 0 0')
Sinv = np.linalg.inv(S)
T = np.matrix('0 0 0 0 1; 0 0 0 1 0; 0 0 1 0 0; 0 1 0 0 0; 1 0 0 0 0')
Tinv = np.linalg.inv(T)


def toint(*args):
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
                    a, b, in1 = toint(a, b, in1)
                    if a == 1 and b == 0:
                        # NOT gate
                        str = "Tinv %s T C" % ms[in1]
                        ms.append(str)
                    else:
                        print("error: only support NOT so far")
                        exit(-1)
                elif arity == 2:
                    _, a, b, c, d, _, _, _, in1, in2, \
                        rest = rest.split(maxsplit=10)
                    a, b, c, d, in1, in2 = toint(a, b, c, d, in1, in2)
                    if a == 0 and b == 0 and c == 0 and d == 1:
                        # AND gate
                        str = "Rinv %s R Sinv %s S inv(Rinv %s R) inv(Sinv %s S)" \
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
    I = np.eye(5)
    comp = np.eye(5)
    for m in ms:
        if m == 'C':
            comp = comp * C
        elif m == 'R':
            comp = comp * R
        elif m == 'Rinv':
            comp = comp * Rinv
        elif m == 'S':
            comp = comp * S
        elif m == 'Sinv':
            comp = comp * Sinv
        elif m == 'T':
            comp = comp * T
        elif m == 'Tinv':
            comp = comp * Tinv
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
