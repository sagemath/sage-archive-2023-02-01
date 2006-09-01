r"""
Modular Forms with Character

EXAMPLES:
    sage: eps = DirichletGroup(13).0
    sage: M = ModularForms(eps^2, 2); M
    Modular Forms space of dimension 3, character [zeta6] and weight 2 over Cyclotomic Field of order 6 and degree 2

    sage: S = M.cuspidal_submodule(); S
    Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3, character [zeta6] and weight 2 over Cyclotomic Field of order 6 and degree 2
    sage: S.modular_symbols()
    Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2

We create a spaces associated to Dirichlet characters of modulus 225:
    sage: e = DirichletGroup(225).0
    sage: e.order()
    6
    sage: e.base_ring()
    Cyclotomic Field of order 60 and degree 16
    sage: M = ModularForms(e,3)


Notice that the base ring is ``minimized'':
    sage: M
    Modular Forms space of dimension 66, character [zeta6, 1] and weight 3
    over Cyclotomic Field of order 6 and degree 2

If we don't want the base ring to change, we can explicitly specify it:
    sage: ModularForms(e, 3, e.base_ring())
    Modular Forms space of dimension 66, character [zeta6, 1] and weight 3
    over Cyclotomic Field of order 60 and degree 16

Next we create a space associated to a Dirichlet character of order 20:
    sage: e = DirichletGroup(225).1
    sage: e.order()
    20
    sage: e.base_ring()
    Cyclotomic Field of order 60 and degree 16
    sage: M = ModularForms(e,17); M
    Modular Forms space of dimension 484, character [1, zeta20] and
    weight 17 over Cyclotomic Field of order 20 and degree 8

We compute the Eisenstein subspace, which is fast even though
the dimension of the space is large (since an explicit basis
of $q$-expansions has not been computed yet).
    sage: M.eisenstein_submodule()
    Eisenstein subspace of dimension 8 of Modular Forms space of
    dimension 484, character [1, zeta20] and weight 17 over Cyclotomic Field of order 20 and degree 8

    sage: M.cuspidal_submodule()
    Cuspidal subspace of dimension 476 of Modular Forms space of dimension 484, character [1, zeta20] and weight 17 over Cyclotomic Field of order 20 and degree 8

"""

#########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import sage.rings.all as rings

import sage.modular.congroup as congroup
import sage.modular.dims as dims
import sage.modular.dirichlet as dirichlet
import sage.modular.modsym.modsym as modsym

import ambient
import cuspidal_submodule
import eisenstein_submodule
import space
import submodule

class ModularFormsAmbient_eps(ambient.ModularFormsAmbient):
    """
    A space of modular forms with character.
    """
    def __init__(self, character, weight=2, base_ring=None):
        """
        Create an ambient modular forms space with character.

        INPUT:
            weight -- int
            character -- dirichlet.DirichletCharacter
            base_ring -- base field
        """
        if not dirichlet.is_DirichletCharacter(character):
            raise TypeError, "character (=%s) must be a Dirichlet character"%character
        if base_ring==None: base_ring=character.base_ring()
        if character.base_ring() != base_ring:
            character = character.change_ring(base_ring)
        group = congroup.Gamma1(character.modulus())
        base_ring = character.base_ring()
        ambient.ModularFormsAmbient.__init__(self, group, weight, base_ring, character)

    def _repr_(self):
        return "Modular Forms space of dimension %s, character %s and weight %s over %s"%(
            self.dimension(), self.character(), self.weight(), self.base_ring())

    def cuspidal_submodule(self):
        """
        Return the cuspidal submodule of this ambient space of modular forms.

        EXAMPLES:
            sage: eps = DirichletGroup(4).0
            sage: M = ModularForms(eps, 5); M
            Modular Forms space of dimension 3, character [-1] and weight 5 over Rational Field
            sage: M.cuspidal_submodule()
            Cuspidal subspace of dimension 1 of Modular Forms space of dimension 3, character [-1] and weight 5 over Rational Field
        """
        try:
            return self.__cuspidal_submodule
        except AttributeError:
            self.__cuspidal_submodule = cuspidal_submodule.CuspidalSubmodule_eps(self)
        return self.__cuspidal_submodule

    def change_ring(self, base_ring):
        return ModularFormsAmbient_eps(character = self.character(),
                                       weight    = self.weight(),
                                       base_ring = base_ring)

    def modular_symbols(self, sign=0):
        """
        Return corresponding space of modular symbols with given sign.

        EXAMPLES:
            sage: eps = DirichletGroup(13).0
            sage: M = ModularForms(eps^2, 2)
            sage: M.modular_symbols()
            Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
            sage: M.modular_symbols(1)
            Modular Symbols space of dimension 3 and level 13, weight 2, character [zeta6], sign 1, over Cyclotomic Field of order 6 and degree 2
            sage: M.modular_symbols(-1)
            Modular Symbols space of dimension 1 and level 13, weight 2, character [zeta6], sign -1, over Cyclotomic Field of order 6 and degree 2
            sage: M.modular_symbols(2)
            Traceback (most recent call last):
            ...
            ValueError: sign must be -1, 0, or 1
        """
        sign = rings.Integer(sign)
        try:
            return self.__modsym[sign]
        except AttributeError:
            self.__modsym = {}
        except KeyError:
            pass
        self.__modsym[sign] = modsym.ModularSymbols(\
            self.character(),
            weight = self.weight(),
            sign = sign,
            base_ring = self.base_ring())
        return self.__modsym[sign]

    def eisenstein_submodule(self):
        try:
            return self.__eisenstein_submodule
        except AttributeError:
            self.__eisenstein_submodule = eisenstein_submodule.EisensteinSubmodule_eps(self)
        return self.__eisenstein_submodule

    ####################################################################
    # Computations of Dimensions
    ####################################################################
    def _dim_cuspidal(self):
        try:
            return self.__the_dim_cuspidal
        except AttributeError:
            self.__the_dim_cuspidal = dims.dimension_cusp_forms_eps(
                self.character(), self.weight())
        return self.__the_dim_cuspidal

    def _dim_eisenstein(self):
        try:
            return self.__the_dim_eisenstein
        except AttributeError:
            self.__the_dim_eisenstein = dims.dimension_eis_eps(self.character(), self.weight())
        return self.__the_dim_eisenstein

    def _dim_new_cuspidal(self):
        try:
            return self.__the_dim_new_cuspidal
        except AttributeError:
            self.__the_dim_new_cuspidal = dims.dimension_new_cusp_forms(
                self.character(), self.weight())
        return self.__the_dim_new_cuspidal

    def _dim_new_eisenstein(self):
        try:
            return self.__the_dim_new_eisenstein
        except AttributeError:
            if self.character().is_trivial() and self.weight() == 2:
                if rings.is_prime(self.level()):
                    d = 1
                else:
                    d = 0
            else:
                E = self.eisenstein_series()
                d = len([g for g in E if g.new_level() == self.level()])
            self.__the_dim_new_eisenstein = d
        return self.__the_dim_new_cuspidal


