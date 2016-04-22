__all__ = ['parse', 'ParseException']

class ParseException(Exception):
    pass

def _parse_param(line):
    try:
        _, param, value = line.split()
    except ValueError:
        raise ParseException("Invalid line '%s'" % line)
    param = param.lower()
    if param in ('nins', 'depth'):
        try:
            value = int(value)
        except ValueError:
            raise ParseException("Invalid value '%s'" % value)
        return {param: value}
    else:
        raise ParseException("Invalid parameter '%s'" % param)

def parse(fname, bp, f_inp_gate, f_gate, keyed=False):
    info = {'nlayers': 0}
    output = False
    with open(fname) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if line == '' or line.startswith('#'):
                continue
            elif line.startswith(':'):
                info.update(_parse_param(line))
                continue
            num, rest = line.split(None, 1)
            try:
                num = int(num)
            except ValueError:
                raise ParseException(
                    'Line %d: gate index not a number' % lineno)
            if rest.startswith('input'):
                if keyed:
                    _, inp = rest.split(None, 1)
                    f_inp_gate(bp, num, inp)
                else:
                    f_inp_gate(bp, num)
                info['nlayers'] += 1
            elif rest.startswith('gate') or rest.startswith('output'):
                if rest.startswith('output'):
                    if output:
                        raise ParseException(
                            'Line %d: only one output gate supported' % lineno)
                    else:
                        output = True
                _, gate, rest = rest.split(None, 2)
                inputs = [int(i) for i in rest.split()]
                try:
                    f_gate(bp, num, lineno, gate.upper(), inputs)
                except KeyError:
                    raise ParseException(
                        'Line %d: unsupported gate %s' % (lineno, gate))
                except TypeError:
                    raise ParseException(
                        'Line %d: incorrect number of arguments given' % lineno)
            else:
                raise ParseException('Line %d: unknown gate type' % lineno)
    if not output:
        raise ParseException('no output gate found')
    return bp[-1], info
