"""
Points on schemes
"""

#*******************************************************************************
#  Copyright (C) 2006 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.structure.element import Element

########################################################
# Base class for points on a scheme, either topological
# or defined by a morphism.
########################################################

class SchemePoint(Element):
    def __init__(self, S):
        """
        INPUT:


        -  ``S`` - a scheme
        """
        Element.__init__(self, S.point_set())
        self.__S = S

    def scheme(self):
        return self.__S

########################################################
# Topological points on a scheme
########################################################

def is_SchemeTopologicalPoint(x):
    return isinstance(x, SchemeTopologicalPoint)

class SchemeTopologicalPoint(SchemePoint):
    pass

class SchemeTopologicalPoint_affine_open(SchemeTopologicalPoint):
    def __init__(self, u, x):
        """
        INPUT:


        -  ``u`` - morphism with domain U an affine scheme

        -  ``x`` - point on U
        """
        SchemePoint.__init__(self, u.codomain())
        self.__u = u
        self.__x = x

    def _repr_(self):
        return "Point on %s defined by x in U, where:\n  U: %s\n  x: %s"%(\
                   self.scheme(), self.embedding_of_affine_open().domain(),
                                  self.point_on_affine())

    def point_on_affine(self):
        """
        Return the scheme point on the affine open U.
        """
        return self.__x

    def affine_open(self):
        """
        Return the affine open subset U.
        """
        return self.__u.domain()

    def embedding_of_affine_open(self):
        """
        Return the embedding from the affine open subset U into this
        scheme.
        """
        return self.__u


class SchemeTopologicalPoint_prime_ideal(SchemeTopologicalPoint):
    def __init__(self, S, P):
        """
        INPUT:


        -  ``S`` - an affine scheme

        -  ``P`` - a prime ideal of the coordinate ring of S
        """
        SchemeTopologicalPoint.__init__(self, S)
        self.__P = P

    def _repr_(self):
        return "Point on %s defined by the prime ideal"%(self.scheme(),
                                                         self.prime_ideal())
    def prime_ideal(self):
        return self.__P

########################################################
# Points on a scheme defined by a morphism
########################################################

def is_SchemeRationalPoint(x):
    return isinstance(x, SchemeRationalPoint)

class SchemeRationalPoint(SchemePoint):
    def __init__(self, f):
        """
        INPUT:


        -  ``f`` - a morphism of schemes
        """
        SchemePoint.__init__(self, f.codomain())
        self.__f = f

    def _repr_(self):
        return "Point on %s defined by the morphism %s"%(self.scheme(),
                                                         self.morphism())

    def morphism(self):
        return self.__f
