#!/usr/bin/env sage -python

from sage.all import *

class Group(object):
    def __init__(self):
        self.length = None
        self.G = None
        self.I = None
        self.C = None
        self.Cc = None
        self.conjugates = None

class S5(Group):
    def __init__(self):
        self.length = 5
        self.G = MatrixSpace(GF(2), self.length)
        self.I = self.G.one()
        self.C = self.G([[0, 1, 0, 0, 0],
                         [0, 0, 1, 0, 0],
                         [0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 1],
                         [1, 0, 0, 0, 0]])
        A = self.G([[0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1],
                    [0, 1, 0, 0, 0],
                    [1, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0]])
        B = self.G([[0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0],
                    [1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 1],
                    [0, 0, 1, 0, 0]])
        Ac = self.G([[1, 0, 0, 0, 0],
                     [0, 0, 1, 0, 0],
                     [0, 1, 0, 0, 0],
                     [0, 0, 0, 0, 1],
                     [0, 0, 0, 1, 0]])
        Aci = Ac.inverse()
        Aic = self.G([[1, 0, 0, 0, 0],
                      [0, 0, 0, 1, 0],
                      [0, 0, 0, 0, 1],
                      [0, 1, 0, 0, 0],
                      [0, 0, 1, 0, 0]])
        Aici = Aic.inverse()
        Bc = self.G([[1, 0, 0, 0, 0],
                     [0, 1, 0, 0, 0],
                     [0, 0, 0, 1, 0],
                     [0, 0, 0, 0, 1],
                     [0, 0, 1, 0, 0]])
        Bci = Bc.inverse()
        Bic = self.G([[1, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0],
                      [0, 0, 0, 0, 1],
                      [0, 0, 0, 1, 0],
                      [0, 1, 0, 0, 0]])
        Bici = Bic.inverse()
        self.Cc = self.G([[0, 0, 0, 0, 1],
                          [0, 0, 0, 1, 0],
                          [0, 0, 1, 0, 0],
                          [0, 1, 0, 0, 0],
                          [1, 0, 0, 0, 0]])
        self.conjugates = {
            'A': (Ac, Aci),
            'B': (Bc, Bci),
            'Ai': (Aic, Aici),
            'Bi': (Bic, Bici)
        }
    def __repr__(self):
        return 'S5'

class S6(Group):
    def __init__(self):
        self.length = 6
        self.G = MatrixSpace(GF(2), self.length)
        self.I = self.G.one()
        self.C = self.G([[0, 0, 0, 1, 0, 0],
                         [0, 0, 1, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0],
                         [1, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 0, 1]])
        A = self.G([[0, 1, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1]])
        B = self.G([[0, 0, 0, 1, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 1, 0]])
        Ac = self.G([[1, 0, 0, 0, 0, 0],
                     [0, 0, 1, 0, 0, 0],
                     [0, 0, 0, 1, 0, 0],
                     [0, 1, 0, 0, 0, 0],
                     [0, 0, 0, 0, 1, 0],
                     [0, 0, 0, 0, 0, 1]])
        Aci = Ac.inverse()
        Bc = self.G([[1, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 1, 0],
                     [0, 0, 0, 0, 0, 1],
                     [0, 0, 0, 1, 0, 0],
                     [0, 1, 0, 0, 0, 0],
                     [0, 0, 1, 0, 0, 0]])
        Bci = Bc.inverse()
        self.Cc = self.G.one()
        self.conjugates = {
            'A': (Ac, Aci),
            'B': (Bc, Bci),
            'Ai': (Ac, Aci),
            'Bi': (Bc, Bci)
        }
    def __repr__(self):
        return 'S6'
        
class SL3(Group):
    '''
    SL3 group.
    NOTE: Doesn't work with randomization!
    '''
    def __init__(self):
        self.length = 3
        self.G = MatrixSpace(GF(3), self.length)
        self.I = self.G.one()
        self.C = self.G([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
        self.Cc = self.G([[-1, 0, -1], [0, -1, 0], [0, 0, 1]])
        A = self.G([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
        B = self.G([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
        Ac = self.G([[-1, 0, 0], [0, 0, -1], [0, -1, 0]])
        Bc = self.G([[0, 1, 0], [1, 0, 1], [-1, 0, 1]])
        Bci = Bc.inverse()
        self.conjugates = {
            'A': (Ac, Ac),
            'B': (Bc, Bci),
            'Ai': (Ac, Ac),
            'Bi': (Bc, Bci)
        }
    def __repr__(self):
        return 'SL3'

groupmap = {
    'S5': S5,
    'S6': S6,
    'SL3': SL3,
}
