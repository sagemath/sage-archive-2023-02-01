"""
Morphism of Hecke modules

AUTHORS:

- William Stein
"""

#*****************************************************************************
#       SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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


import sage.misc.misc as misc
import sage.categories.all as cat
from sage.modules.matrix_morphism import MatrixMorphism

# We also define other types of Hecke-module morphisms that aren't
# specified by a matrix.  E.g., Hecke operators, or maybe morphisms on
# modular abelian varieties (which are specified by matrices, but on
# integral homology). All morphisms derive from heckeModuleMorphism.

def is_HeckeModuleMorphism(x):
    return isinstance(x, HeckeModuleMorphism)

def is_HeckeModuleMorphism_matrix(x):
    return isinstance(x, HeckeModuleMorphism_matrix)

class HeckeModuleMorphism(cat.Morphism):
    pass

class HeckeModuleMorphism_matrix(MatrixMorphism, HeckeModuleMorphism):
    """
    Morphisms of Hecke modules when the morphism is given by a matrix.
    """
    def __init__(self, parent, A, name=''):
        """
        INPUT:


        -  ``parent`` - ModularSymbolsHomspace

        -  ``A`` - Matrix

        -  ``name`` - str (defaults to ") name of the morphism
           (used for printing)
        """
        if not isinstance(name, str):
            raise TypeError, "name must be a string"
        self.__name = name
        MatrixMorphism.__init__(self, parent, A)

    def name(self, new=None):
        if new is None:
            return self.__name
        self.__name = new

    def _repr_(self):
        if max(self.matrix().nrows(),self.matrix().ncols()) > 5:
            mat = "(not printing %s x %s matrix)"%(self.matrix().nrows(), self.matrix().ncols())
        else:
            mat = str(self.matrix())
        name = self.__name
        if name != '':
            name += ' '
        return "Hecke module morphism %sdefined by the matrix\n%s\nDomain: %s\nCodomain: %s"%(\
                name, mat, misc.strunc(self.domain()), misc.strunc(self.codomain()))

    def __mul__(self, right):
        if not isinstance(right, HeckeModuleMorphism_matrix):
            try:
                right = right.hecke_module_morphism()
            except AttributeError:
                R = self.base_ring()
                return self.parent()(self.matrix() * R(right))
        H = right.domain().Hom(self.codomain())
        return HeckeModuleMorphism_matrix(H, right.matrix() * self.matrix())

