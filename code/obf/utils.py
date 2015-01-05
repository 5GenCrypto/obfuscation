#!/usr/bin/env sage -python

from __future__ import print_function
import functools, sys

def clr_error(s):
    return '\x1b[31m%s\x1b[0m' % s
def clr_warn(s):
    return '\x1b[33m%s\x1b[0m' % s
def clr_ok(s):
    return '\x1b[32m%s\x1b[0m' % s

def logger(s, end='\n', verbose=False):
    if verbose:
        print(s, end=end, file=sys.stderr)
        sys.stderr.flush()

def make_logger(verbose):
    return functools.partial(logger, verbose=verbose)
