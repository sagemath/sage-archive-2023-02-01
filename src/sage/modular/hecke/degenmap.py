"""
Degeneracy maps
"""

#*****************************************************************************
#       SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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


import morphism
import homspace

class DegeneracyMap(morphism.HeckeModuleMorphism_matrix):
    """
    A degneracy map between Hecke modules of different levels.

    EXAMPLES:
    We construct a number of degeneracy maps.

        sage: M = ModularSymbols(33)
        sage: d = M.degeneracy_map(11)
        sage: d
        Hecke module morphism degeneracy map corresponding to f(q) |--> f(q) defined by the matrix
        (not printing 9 x 3 matrix)
        Domain: Full Modular Symbols space for Gamma_0(33) of weight 2 with sign ...
        Codomain: Full Modular Symbols space for Gamma_0(11) of weight 2 with sign ...
        sage: d.t()
        1
        sage: d = M.degeneracy_map(11,3)
        sage: d.t()
        3

    The parameter d must be a divisor of the quotient of the two levels.

        sage: d = M.degeneracy_map(11,2)
        Traceback (most recent call last):
        ...
        ValueError: The level of self (=33) must be a divisor or multiple of level (=11), and t (=2) must be a divisor of the quotient.

    Degeneracy maps can also go from lower level to higher level.

        sage: M.degeneracy_map(66,2)
        Hecke module morphism degeneracy map corresponding to f(q) |--> f(q^2) defined by the matrix
        (not printing 9 x 25 matrix)
        Domain: Full Modular Symbols space for Gamma_0(33) of weight 2 with sign ...
        Codomain: Full Modular Symbols space for Gamma_0(66) of weight 2 with sign ...
    """
    def __init__(self, matrix, domain, codomain, t):
        self.__t = t
        H = homspace.HeckeModuleHomspace(domain, codomain)
        if t == 1:
            pow = ""
        else:
            pow = "^%s"%t
        name = "degeneracy map corresponding to f(q) |--> f(q%s)"%(pow)
        morphism.HeckeModuleMorphism_matrix.__init__(self, H, matrix, name)

    def t(self):
        """
        Return the divisor of the quotient of the two levels
        associated to the degeneracy map.

        EXAMPLES:
            sage: M = ModularSymbols(33)
            sage: d = M.degeneracy_map(11,3)
            sage: d.t()
            3
            sage: d = M.degeneracy_map(11,1)
            sage: d.t()
            1
        """
        return self.__t

    def image(self):
        V = morphism.HeckeModuleMorphism_matrix.image(self)
        if self.domain().is_anemic_hecke_module():
            V._is_anemic_hecke_module = True
        return V

    def kernel(self):
        V = morphism.HeckeModuleMorphism_matrix.kernel(self)
        if self.domain().is_anemic_hecke_module():
            V._is_anemic_hecke_module = True
        return V


