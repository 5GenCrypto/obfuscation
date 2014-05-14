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

def init_wide():
    fig_width_pt = 550.0
    inches_per_pt = 1.0 / 72.27
    golden_mean = (pylab.sqrt(5) - 1.0) / 5.0
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
    
def extract(line):
    return float(line.rstrip().rsplit(' ', 1)[1])

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

def extract_obf_time(fname):
    obf_count = 0
    obf_num = 0
    mlm = False
    lst = []
    with open(fname, 'r') as f:
        while True:
            try:
                line = f.next().strip()
            except StopIteration:
                break
            if line.startswith('Generating MLM'):
                mlm = True
            elif line.startswith('Generating p_i'):
                r = extract(line)
                print('p_i/g_i: %f' % round(r, 2))
                lst.append(r)
            elif line.startswith('Generating CRT'):
                r = extract(line)
                print('CRT: %f' % round(r, 2))
                lst.append(r)
            elif line.startswith('Generating z_i'):
                r = extract(line)
                print('z_i: %f' % round(r, 2))
                lst.append(r)
            elif line.startswith('Generating pzt'):
                r = extract(line)
                print('pzt: %f' % round(r, 2))
                lst.append(r)
            elif line.startswith('Took') and mlm:
                r = extract(line)
                print('MLM Total: %f' % round(r, 2))
                lst.append(r)
                mlm = False
            elif line.startswith('Randomizing BPs'):
                line = f.next().strip()
                r = extract(line)
                print('BP rand: %f' % round(r, 2))
                lst.append(r)
            elif line.startswith('Constructing bookend'):
                _ = f.next()
                a = f.next()
                _ = f.next()
                b = f.next()
                _ = f.next()
                total = extract(f.next())
                a = extract(a)
                b = extract(b)
                avg = (a + b) / 2
                print('Bookend enc (avg): %f' % round(avg, 2))
                print('Bookend enc (total): %f' % round(total, 2))
                lst.append(total)
            elif line.startswith('Obfuscating layer'):
                _ = f.next()
                _ = f.next()
                a = f.next().strip()
                obf_count += extract(a)
                obf_num += 1
            elif line.startswith('Obfuscation took'):
                print('Layer enc (avg): %f' % round(obf_count / obf_num, 1))
                print('Layer enc (total): %f' % round(obf_count, 1))
                lst.append(obf_count)
                r = extract(line)
                print('Obf Total: %f' % round(r, 2))
                lst.append(r)
            elif line.startswith('Max memory'):
                print('RAM (GB): %f' % round(extract(line) / 1024 / 1024, 1))
                break
    return lst


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

    
