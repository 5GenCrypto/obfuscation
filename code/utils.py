#!/usr/bin/env sage -python

from __future__ import print_function
import functools, sys

def logger(s, end='\n', verbose=False):
    if verbose:
        print(s, end=end)
        sys.stdout.flush()

def make_logger(verbose):
    return functools.partial(logger, verbose=verbose)
