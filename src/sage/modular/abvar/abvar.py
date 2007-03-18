"""
Base class for modular abelian varieties
"""

###########################################################################
#       Copyright (C) 2004,2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from sage.structure.sage_object import SageObject

from hecke_operator             import HeckeOperator
from torsion_subgroup           import TorsionSubgroup

class ModularAbelianVariety(SageObject):
    def __init__(self, level, base_ring):
        self._level = level
        self._base_ring = base_ring

    def _repr_(self):
        return "Modular abelian variety of level %s over %s"%(self._level, self._base_ring)

    def dimension(self):
        raise NotImplementedError

    def level(self):
        return self._level

    def base_ring(self):
        return self._base_ring

    def homology(self):
        raise NotImplementedError

    def hecke_operator(self, n):
        try:
            return self._hecke_operator[n]
        except AttributeError:
            self._hecke_operator = {}
        except KeyError:
            pass
        Tn = HeckeOperator(self, n)
        self._hecke_operator[n] = Tn
        return Tn

    def torsion_subgroup(self):
        try:
            return self._torsion_subgroup
        except AttributeError:
            T = TorsionSubgroup(self)
            self._torsion_subgroup = T
            return T

    def cuspidal_subgroup(self):
        try:
            return self._cuspidal_subgroup
        except AttributeError:
            T = CuspidalSubgroup(self)
            self._cuspidal_subgroup = T
            return T

    def rational_cuspidal_subgroup(self):
        try:
            return self._rational_cuspidal_subgroup
        except AttributeError:
            T = RationalCuspidalSubgroup(self)
            self._rational_cuspidal_subgroup = T
            return T



