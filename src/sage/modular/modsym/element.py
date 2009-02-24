"""
A single element of an ambient space of modular symbols.
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

import operator

import sage.modules.free_module_element
import sage.modules.module_element as module_element
import sage.misc.misc as misc
import sage.structure.element as element
import sage.structure.formal_sum as formal_sum
import ambient
import sage.modular.hecke.all as hecke
import sage.misc.latex as latex

_print_mode = "manin"

def is_ModularSymbolsElement(x):
    return isinstance(x, ModularSymbolsElement)

def set_modsym_print_mode(mode="manin"):
    """
    Set the mode for printing of elements of modular symbols spaces.

    INPUT:


    -  ``mode`` - a string. The possibilities are as
       follows:

    -  ``'manin'`` - (the default) formal sums of Manin
       symbols [P(X,Y),(u,v)]

    -  ``'modular'`` - formal sums of Modular symbols
       P(X,Y)\*alpha,beta, where alpha and beta are cusps

    -  ``'vector'`` - as vectors on the basis for the
       ambient space
    """
    mode = str(mode).lower()
    if not (mode in ['manin', 'modular', 'vector']):
        raise ValueError, "mode must be one of 'manin', 'modular', or 'vector'"
    global _print_mode
    _print_mode = mode

class ModularSymbolsElement(hecke.HeckeModuleElement):
    """
    An element of a space of modular symbol.
    """
    def __init__(self, parent, x, check=True):
        """
        INPUT:


        -  ``parent`` - a space of modular symbols

        -  ``x`` - a free module element that represents the
           modular symbol in terms of a basis for the ambient space (not in
           terms of a basis for parent!)
        """
        if check:
            if not isinstance(parent, ambient.ModularSymbolsAmbient):
                raise TypeError, "parent must be an ambient space of modular symbols."
            if not isinstance(x, sage.modules.free_module_element.FreeModuleElement):
                raise TypeError, "x must be a free module element."
            if x.degree() != parent.degree():
                raise TypeError, "x (of degree %s) must be of degree the same as the degree of the parent (of degree %s)."%(x.degree(), parent.degree())
        hecke.HeckeModuleElement.__init__(self, parent, x)

    def __is_compatible(self, other):
        return isinstance(other, ModularSymbolsElement) and self.parent() == other.parent()

    def _add_(self, right):
        return ModularSymbolsElement(self.parent(), self.element() + right.element(), check=False)

    def __cmp__(self, other):
        return self.element().__cmp__(other.element())

    def _repr_(self):
        if _print_mode == "vector":
            return str(self.element())
        elif _print_mode == "manin":
            m = self.manin_symbol_rep()
        elif _print_mode == "modular":
            m = self.modular_symbol_rep()
        c = [x[0] for x in m]
        v = [x[1] for x in m]
        return misc.repr_lincomb(v, c)

    def _latex_(self):
        if _print_mode == "vector":
            return latex(self.element())
        elif _print_mode == "manin":
            m = self.manin_symbol_rep()
        elif _print_mode == "modular":
            m = self.modular_symbol_rep()
        c = [x[0] for x in m]
        v = [x[1] for x in m]
        return latex.repr_lincomb(v, c)


    # TODO -- use module machinery
    def __mul__(self, right):
        return ModularSymbolsElement(self.parent(), self.element()*right)

    def __neg__(self):
        return ModularSymbolsElement(self.parent(), -self.element())

    def _sub_(self, right):
        return ModularSymbolsElement(self.parent(), self.element() - right.element())

    def coordinate_vector(self):
        if self.parent().is_ambient():
            return self.element()
        return self.parent().embedded_vector_space().coordinate_vector(self.element())

    def list(self):
        return self.element().list()

    def manin_symbol_rep(self):
        """
        Returns a representation of self as a formal sum of Manin symbols.

        (The result is cached for future use.)
        """
        try:
            return self.__manin_symbols
        except AttributeError:
            A = self.parent()
            v = self.element()
            manin_symbols = A.ambient_hecke_module().manin_symbols_basis()
            F = formal_sum.FormalSums(A.base_ring())
            ms = F([(v[i], manin_symbols[i]) for i in \
                  range(v.degree()) if v[i] != 0], check=False, reduce=False)
            self.__manin_symbols = ms
        return self.__manin_symbols

    def modular_symbol_rep(self):
        """
        Returns a representation of self as a formal sum of modular
        symbols.

        (The result is cached for future use.)
        """
        try:
            return self.__modular_symbols
        except AttributeError:
            A = self.parent()
            v = self.manin_symbol_rep()
            if v == 0:
                return v
            w = [c * x.modular_symbol_rep() for c, x in v]
            return sum(w)
        return self.__modular_symbols


