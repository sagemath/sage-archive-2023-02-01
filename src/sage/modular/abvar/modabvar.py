"""
Modular abelian varieties
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import weakref

import sage.modular.hecke.all as hecke
import sage.modular.modsym as modsym
import sage.rings.all as rings

def J0(N):
    J = ModularAbelianVariety(N)
    J._AbelianVariety__repr = "J0(%s)"%N
    return J

_objsAbelianVariety = {}
class _uniqAbelianVariety(object):
    def __new__(cls, level=1, weight=2, character=None, base_ring=rings.RationalField()):
        key = (level, weight, character, base_ring)
        if _objsAbelianVariety.has_key(key):
            x = _objsAbelianVariety[key]()
            if x != None: return x
        R = object.__new__(cls)
        _objsAbelianVariety[key] = weakref.ref(R)
        return R

class ModularAbelianVariety(_uniqAbelianVariety, hecke.HeckeModule_generic):
    def __init__(self, level=1, weight=2, character=None, base_ring=rings.RationalField()):
        hecke.HeckeModule.__init__(self)
        self.__repr = None
        if isinstance(level, modsym.ModularSymbols):
            M = level
            self.__level = M.level()
            self.__weight = M.weight()
            self.__character = M.character()
            self.__base_ring = M.base_ring()
            self.__homology = Homology(M)
            return
        self.__level = level
        self.__weight = weight
        self.__character = character
        self.__base_ring = base_ring
        self.__homology = Homology(modsym.ModularSymbols(\
            level=level,weight=weight,character=character,base_ring=base_ring, sign=0))

    def __repr__(self):
        if self.__repr == None:
            self.__repr="Modular abelian variety (N=%s,k=%s,eps=%s)"%\
                         (self.__level, self.__weight, self.__character)
        return self.__repr
    def dimension(self):
        raise NotImplementedError
    def level(self):
        return self.__level
    def weight(self):
        return self.__weight
    def base_ring(self):
        return self.__base_ring
    def homology(self):
        return self.__homology
    def hecke_operator(self, n):
        """Hecke operator on modular abelian variety."""
        raise NotImplementedError
    def torsion_multiple(self, stop=50):
        raise NotImplementedError
    def cuspidal_subgroup(self):
        raise NotImplementedError
    def rational_cuspidal_subgroup(self):
        raise NotImplementedError


class TorsionPoint:
    def __init__(self, parent, x):
        self.__parent = parent
        self.__x = x
    def __repr__(self):
        return "Torsion point '%s' on %s"%\
               (self.__parent, self.__x)

class HomSpace:
    def __init__(self, domain, codomain, x=None):
        self.__domain = domain
        self.__codomain = codomain
        self.__x = x
    def __repr__(self):
        return "Space of homomorphisms from %s to %s defined by %s"%\
               (self.__domain, self.__codomain, self.__x)

class LSeries:
    def __init__(self, parent):
        self.__parent = parent
    def __repr__(self):
        return "LSeries L(%s,s)"%self.__parent
    def __call__(self, x):
        raise NotImplementedError
    def rational_part(self, x):
        raise NotImplementedError

class Homology:
    def __init__(self, modsym):
        self.__modsym = modsym
    def __repr__(self):
        return "Homology associated to %s"%self.__modsym
