#!/usr/bin/env sage -python

from gradedencoding import GradedEncoding

class Obfuscator(object):
    def __init__(self, bp, secparam=16):
        self.ge = GradedEncoding(secparam=secparam, kappa=len(bp))
        


