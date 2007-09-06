r"""
PARI Groups
"""

import group
from sage.libs.all import pari, pari_gen
from sage.rings.all import Integer

class PariGroup(group.FiniteGroup):
    def __init__(self, x, degree=None):
        if not isinstance(x, pari_gen):
            raise TypeError, "x (=%s) must be a PARI gen"%x
        self.__x = x
        self.__degree = degree

    def __repr__(self):
        return "PARI group %s of degree %s"%(self.__x, self.__degree)

    def _pari_(self):
        return self.__x

    def degree(self):
        return self.__degree

    def order(self):
        return Integer(self.__x[0])

    def permutation_group(self):
        if self.__degree is None:
            raise NotImplementedError
        import perm_gps.permgroup_named
        return perm_gps.permgroup_named.TransitiveGroup(self.__degree, self.__x[2])

    _permgroup_ = permutation_group
