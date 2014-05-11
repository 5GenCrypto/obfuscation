import os

import pylab

def init():
    fig_width_pt = 240.0
    inches_per_pt = 1.0 / 72.27
    golden_mean = (pylab.sqrt(5) - 1.0) / 2.0
    fig_width = fig_width_pt * inches_per_pt
    fig_height = fig_width * golden_mean
    fig_size = [fig_width, fig_height]
    params = {'backend': 'ps',
              'axes.labelsize': 8,
              'text.fontsize': 8,
              'legend.fontsize': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'text.usetex': True,
              'figure.figsize': fig_size
    }
    pylab.rcParams.update(params)


def obftime(fname):
    secparam = [int(s) for s in fname.split('-') if s.isdigit()][0]
    with open(fname, 'r') as f:
        while True:
            try:
                line = f.next()
            except StopIteration:
                break
            if line.startswith('Obfuscation took'):
                a, b, c = line.rstrip().rsplit(' ', 2)
                if c == 'seconds':
                    return secparam, float(b)

def obfN(fname):
    secparam = [int(s) for s in fname.split('-') if s.isdigit()][0]
    with open(fname, 'r') as f:
        for line in f:
            if line.strip().startswith('N:'):
                _, b = line.rstrip().rsplit(' ', 1)
                return secparam, float(b)

def obfsize(fname):
    secparam = [int(s) for s in fname.split('-') if s.isdigit()][0]
    with open(fname, 'r') as f:
        total = 0
        for line in f:
            size, _ = line.strip().split('\t', 1)
            total += int(size)
        return secparam, int(total)

def evaltime(fname):
    secparam = [int(s) for s in fname.split('-') if s.isdigit()][0]
    with open(fname, 'r') as f:
        while True:
            try:
                line = f.next()
            except StopIteration:
                raise Exception('Could not find time info')
            if line.startswith('Took:'):
                a, b = line.rstrip().rsplit(' ', 1)
                return secparam, float(b)

def dir(directory, name, func):
    results = []
    for file in os.listdir(directory):
        if name in file:
            path = os.path.join(directory, file)
            result = func(path)
            results.append(result)
    results = sorted(results)
    xs = []
    ys = []
    for x, y in results:
        xs.append(x)
        ys.append(y)
    return xs, ys

def dir_obftime(directory):
    return dir(directory, 'obf-time', obftime)
def dir_obfsize(directory):
    return dir(directory, 'obf-size', obfsize)
def dir_evaltime(directory):
    return dir(directory, 'eval-time', evaltime)

    
