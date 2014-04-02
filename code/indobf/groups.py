import numpy as np

class Mmod(object):
    def __init__(self, m, p):
        if isinstance(m, Mmod):
            assert m.p <= p
            m = m.m
        self.m = np.mod(np.array(m, dtype=np.object), p)
        self.p = p
    def __str__(self):
        return str(self.m)
    def __mul__(self, other):
        assert self.p == other.p
        return Mmod(np.dot(self.m, other.m), self.p)
    def __rmul__(self, other):
        assert other.p == self.p
        return Mmod(np.dot(other.m, self.m), self.p)
    def __eq__(self, other):
        assert self.p == other.p
        return np.array_equal(self.m, other.m)

class Group(object):
    def __init__(self):
        self.length = None
        self.I = None
        self.C = None
        self.Cc = None
        self.conjugates = None

class S5(Group):
    def __init__(self):
        self.length = 5
        self.I = Mmod([[1L, 0L, 0L, 0L, 0L],
                       [0L, 1L, 0L, 0L, 0L],
                       [0L, 0L, 1L, 0L, 0L],
                       [0L, 0L, 0L, 1L, 0L],
                       [0L, 0L, 0L, 0L, 1L]], 2)
        self.C = Mmod([[0L, 1L, 0L, 0L, 0L],
                       [0L, 0L, 1L, 0L, 0L],
                       [0L, 0L, 0L, 1L, 0L],
                       [0L, 0L, 0L, 0L, 1L],
                       [1L, 0L, 0L, 0L, 0L]], 2)
        A = Mmod([[0L, 0L, 1L, 0L, 0L],
                  [0L, 0L, 0L, 0L, 1L],
                  [0L, 1L, 0L, 0L, 0L],
                  [1L, 0L, 0L, 0L, 0L],
                  [0L, 0L, 0L, 1L, 0L]], 2)
        B = Mmod([[0L, 1L, 0L, 0L, 0L],
                  [0L, 0L, 0L, 1L, 0L],
                  [1L, 0L, 0L, 0L, 0L],
                  [0L, 0L, 0L, 0L, 1L],
                  [0L, 0L, 1L, 0L, 0L]], 2)
        Ac = Mmod([[1L, 0L, 0L, 0L, 0L],
                   [0L, 0L, 1L, 0L, 0L],
                   [0L, 1L, 0L, 0L, 0L],
                   [0L, 0L, 0L, 0L, 1L],
                   [0L, 0L, 0L, 1L, 0L]], 2)
        Aic = Mmod([[1L, 0L, 0L, 0L, 0L],
                    [0L, 0L, 0L, 1L, 0L],
                    [0L, 0L, 0L, 0L, 1L],
                    [0L, 1L, 0L, 0L, 0L],
                    [0L, 0L, 1L, 0L, 0L]], 2)
        Bc = Mmod([[1L, 0L, 0L, 0L, 0L],
                   [0L, 1L, 0L, 0L, 0L],
                   [0L, 0L, 0L, 1L, 0L],
                   [0L, 0L, 0L, 0L, 1L],
                   [0L, 0L, 1L, 0L, 0L]], 2)
        Bci = Mmod([[1L, 0L, 0L, 0L, 0L],
                    [0L, 1L, 0L, 0L, 0L],
                    [0L, 0L, 0L, 0L, 1L],
                    [0L, 0L, 1L, 0L, 0L],
                    [0L, 0L, 0L, 1L, 0L]], 2)
        Bic = Mmod([[1L, 0L, 0L, 0L, 0L],
                    [0L, 0L, 1L, 0L, 0L],
                    [0L, 0L, 0L, 0L, 1L],
                    [0L, 0L, 0L, 1L, 0L],
                    [0L, 1L, 0L, 0L, 0L]], 2)
        Bici = Mmod([[1L, 0L, 0L, 0L, 0L],
                     [0L, 0L, 0L, 0L, 1L],
                     [0L, 1L, 0L, 0L, 0L],
                     [0L, 0L, 0L, 1L, 0L],
                     [0L, 0L, 1L, 0L, 0L]], 2)
        self.Cc = Mmod([[0L, 0L, 0L, 0L, 1L],
                        [0L, 0L, 0L, 1L, 0L],
                        [0L, 0L, 1L, 0L, 0L],
                        [0L, 1L, 0L, 0L, 0L],
                        [1L, 0L, 0L, 0L, 0L]], 2)
        self.conjugates = {
            'A': (Ac, Ac),
            'B': (Bc, Bci),
            'Ai': (Aic, Aic),
            'Bi': (Bic, Bici)
        }
    def __repr__(self):
        return 'S5'

class S6(Group):
    def __init__(self):
        self.length = 6
        self.I = Mmod([[1L, 0L, 0L, 0L, 0L, 0L],
                       [0L, 1L, 0L, 0L, 0L, 0L],
                       [0L, 0L, 1L, 0L, 0L, 0L],
                       [0L, 0L, 0L, 1L, 0L, 0L],
                       [0L, 0L, 0L, 0L, 1L, 0L],
                       [0L, 0L, 0L, 0L, 0L, 1L]], 2)
        self.C = Mmod([[0L, 0L, 0L, 1L, 0L, 0L],
                       [0L, 0L, 1L, 0L, 0L, 0L],
                       [0L, 1L, 0L, 0L, 0L, 0L],
                       [1L, 0L, 0L, 0L, 0L, 0L],
                       [0L, 0L, 0L, 0L, 1L, 0L],
                       [0L, 0L, 0L, 0L, 0L, 1L]], 2)
        A = Mmod([[0L, 1L, 0L, 0L, 0L, 0L],
                  [1L, 0L, 0L, 0L, 0L, 0L],
                  [0L, 0L, 0L, 1L, 0L, 0L],
                  [0L, 0L, 1L, 0L, 0L, 0L],
                  [0L, 0L, 0L, 0L, 1L, 0L],
                  [0L, 0L, 0L, 0L, 0L, 1L]], 2)
        B = Mmod([[0L, 0L, 0L, 1L, 0L, 0L],
                  [0L, 1L, 0L, 0L, 0L, 0L],
                  [0L, 0L, 1L, 0L, 0L, 0L],
                  [1L, 0L, 0L, 0L, 0L, 0L],
                  [0L, 0L, 0L, 0L, 0L, 1L],
                  [0L, 0L, 0L, 0L, 1L, 0L]], 2)
        Ac = Mmod([[1L, 0L, 0L, 0L, 0L, 0L],
                   [0L, 0L, 1L, 0L, 0L, 0L],
                   [0L, 0L, 0L, 1L, 0L, 0L],
                   [0L, 1L, 0L, 0L, 0L, 0L],
                   [0L, 0L, 0L, 0L, 1L, 0L],
                   [0L, 0L, 0L, 0L, 0L, 1L]], 2)
        Aci = Mmod([[1L, 0L, 0L, 0L, 0L, 0L],
                    [0L, 0L, 0L, 1L, 0L, 0L],
                    [0L, 1L, 0L, 0L, 0L, 0L],
                    [0L, 0L, 1L, 0L, 0L, 0L],
                    [0L, 0L, 0L, 0L, 1L, 0L],
                    [0L, 0L, 0L, 0L, 0L, 1L]], 2)
        Bc = Mmod([[1L, 0L, 0L, 0L, 0L, 0L],
                   [0L, 0L, 0L, 0L, 1L, 0L],
                   [0L, 0L, 0L, 0L, 0L, 1L],
                   [0L, 0L, 0L, 1L, 0L, 0L],
                   [0L, 1L, 0L, 0L, 0L, 0L],
                   [0L, 0L, 1L, 0L, 0L, 0L]], 2)
        self.Cc = Mmod([[0L, 0L, 0L, 1L, 0L, 0L],
                        [0L, 0L, 1L, 0L, 0L, 0L],
                        [0L, 1L, 0L, 0L, 0L, 0L],
                        [1L, 0L, 0L, 0L, 0L, 0L],
                        [0L, 0L, 0L, 0L, 1L, 0L],
                        [0L, 0L, 0L, 0L, 0L, 1L]], 2)
        self.conjugates = {
            'A': (Ac, Aci),
            'B': (Bc, Bc),
            'Ai': (Ac, Aci),
            'Bi': (Bc, Bc)
        }
    def __repr__(self):
        return 'S6'

# class SL3(Group):
#     '''
#     SL3 group.
#     NOTE: Doesn't work with randomization!
#     '''
#     def __init__(self):
#         self.length = 3
#         self.G = MatrixSpace(GF(3), self.length)
#         self.I = self.G.one()
#         self.C = self.G([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
#         self.Cc = self.G([[-1, 0, -1], [0, -1, 0], [0, 0, 1]])
#         A = self.G([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
#         B = self.G([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
#         Ac = self.G([[-1, 0, 0], [0, 0, -1], [0, -1, 0]])
#         Bc = self.G([[0, 1, 0], [1, 0, 1], [-1, 0, 1]])
#         Bci = Bc.inverse()
#         self.conjugates = {
#             'A': (Ac, Ac),
#             'B': (Bc, Bci),
#             'Ai': (Ac, Ac),
#             'Bi': (Bc, Bci)
#         }
#     def __repr__(self):
#         return 'SL3'

groupmap = {
    'S5': S5,
    'S6': S6,
    # 'SL3': SL3,
}
