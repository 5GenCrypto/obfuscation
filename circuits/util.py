import math, os, subprocess

def digit_to_char(digit):
    if digit < 10:
        return str(digit)
    return chr(ord('a') + digit - 10)

def str_base(number,base):
    if number < 0:
        return '-' + str_base(-number, base)
    (d, m) = divmod(number, base)
    if d > 0:
        return str_base(d, base) + digit_to_char(m)
    return digit_to_char(m)

def dary_repr(number, d, n):
    repr = list(str(str_base(number, d)))
    repr = (['0'] * (n - len(repr))) + repr
    return "".join(repr)

# turns a "10" in any base to a 001 000
def digit_dary_repr(num_str, d):
    L = []
    for x in num_str:
        L.append(format(int(x), 'b').zfill(int(math.ceil(math.log(d, 2)))))
    return "".join(L)

def run(lst):
    print('%s' % ' '.join(lst))
    with open(os.devnull, 'w') as fnull:
        return subprocess.call(lst, stdout=fnull, stderr=fnull)

