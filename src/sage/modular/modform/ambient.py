"""
Modular Forms
"""

#########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

# system packages
import math
import weakref

# SAGE packages
import sage.rings.all as rings
import sage.modular.congroup as congroup
import sage.misc.db as db
import sage.modular.dims as dims
import sage.modular.dirichlet as dirichlet
import sage.modular.hecke.all as hecke
import sage.misc.misc as misc
import sage.modular.modsym.all as modsym
import sage.modules.free_module as free_module
import sage.modules.free_module_element as free_module_element
import sage.rings.all as rings

from sage.misc.all import latex

import cuspidal_submodule
import defaults
import eisenstein_submodule
import eis_series
import space
import submodule

class ModularFormsAmbient(space.ModularFormsSpace,
                          hecke.AmbientHeckeModule):
    """
    An ambient space of modular forms.
    """
    def __init__(self, group, weight, base_ring, character=None):
        if not isinstance(group, congroup.CongruenceSubgroup):
            raise TypeError, 'group (=%s) is a group'%group
        weight = rings.Integer(weight)

        if character is None and isinstance(group, congroup.Gamma0):
            character = dirichlet.TrivialCharacter(group.level(), base_ring)

        space.ModularFormsSpace.__init__(self, group, weight, character, base_ring)
        hecke.AmbientHeckeModule.__init__(self, base_ring, self.dimension(), group.level(), weight)

    def _repr_(self):
        return "Modular Forms space of dimension %s for %s of weight %s over %s"%(
                self.dimension(), self.group(), self.weight(), self.base_ring())

    def change_ring(self, base_ring):
        """
        Change the base ring of this space of modular forms.

        INPUT:
            R -- ring

        EXAMPLES:
            sage: M = ModularForms(Gamma0(37),2)
            sage: M.basis()
            [
            q + q^3 - 2*q^4 + O(q^6),
            q^2 + 2*q^3 - 2*q^4 + q^5 + O(q^6),
            1 + 2/3*q + 2*q^2 + 8/3*q^3 + 14/3*q^4 + 4*q^5 + O(q^6)
            ]

        The basis after changing the base ring is the reduction modulo
        $3$ of an integral basis.
            sage: M3 = M.change_ring(GF(3))
            sage: M3.basis()
            [
            ?????
            ]
        """
        import constructor
        M = constructor.ModularForms(self.group(), self.weight(), base_ring, prec=self.prec())
        return M

    def dimension(self):
        try:
            return self.__dimension
        except AttributeError:
            self.__dimension = self._dim_eisenstein() + self._dim_cuspidal()
            return self.__dimension

    def rank(self):
        return self.dimension()

    def ambient_space(self):
        return self

    def is_ambient(self):
        """
        Return True if this an ambient space of modular forms.

        This is an ambient space, so this function always returns True.

        EXAMPLES:
            sage: ModularForms(11).is_ambient()
            True
            sage: CuspForms(11).is_ambient()
            False
        """
        return True

    def modular_symbols(self, sign=0):
        """
        Return the corresponding space of modular symbols with the given sign.

        EXAMPLES:
            sage: S = ModularForms(11,2)
            sage: S.modular_symbols()
            Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: S.modular_symbols(sign=1)
            Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: S.modular_symbols(sign=-1)
            Modular Symbols space of dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Rational Field

            sage: ModularForms(1,12).modular_symbols()
            Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field
        """
        sign = rings.Integer(sign)
        try:
            return self.__modular_symbols[sign]
        except AttributeError:
            self.__modular_symbols = {}
        except KeyError:
            pass
        M = modsym.ModularSymbols(group = self.group(),
                                  weight = self.weight(),
                                  sign = sign,
                                  base_ring = self.base_ring())
        self.__modular_symbols[sign] = M
        return M

    def module(self):
        if hasattr(self, "__module"): return self.__module
        self.__module = free_module.VectorSpace(self.base_ring(),
                                                 self.dimension())
        return self.__module

    def prec(self, new_prec=None):
        """
        Set or get default initial precision for printing modular forms.

        INPUT:
            new_prec -- positive integer (default: None)

        OUTPUT:
            if new_prec is None, returns the current precision.

        EXAMPLES:
            sage: M = ModularForms(1,12, prec=3)
            sage: M.prec()
            3

            sage: M.basis()
            [
            q - 24*q^2 + O(q^3),
            1 + 65520/691*q + 134250480/691*q^2 + O(q^3)
            ]

            sage: M.prec(5)
            sage: M.basis()
            [
            q - 24*q^2 + 252*q^3 - 1472*q^4 + O(q^5),
            1 + 65520/691*q + 134250480/691*q^2 + 11606736960/691*q^3 + 274945048560/691*q^4 + O(q^5)
            ]

        """
        if new_prec:
            self.__prec = new_prec
        try:
            return self.__prec
        except AttributeError:
            self.__prec = defaults.DEFAULT_PRECISION
        return self.__prec

    def set_precision(self, n):
        if n < 0:
            raise ValueError, "n (=%s) must be >= 0"%n
        self.__prec = rings.Integer(n)

    ####################################################################
    # Computation of Special Submodules
    ####################################################################
    def cuspidal_submodule(self):
        """
        EXAMPLES:
            sage: ModularForms(Gamma1(13)).cuspidal_submodule()
            Cuspidal subspace of dimension 2 of Modular Forms space of dimension 13 for
            Congruence Subgroup Gamma1(13) of weight 2 over Rational Field
        """
        try:
            return self.__cuspidal_submodule
        except AttributeError:
            self.__cuspidal_submodule = cuspidal_submodule.CuspidalSubmodule(self)
        return self.__cuspidal_submodule

    def eisenstein_submodule(self):
        try:
            return self.__eisenstein_submodule
        except AttributeError:
            self.__eisenstein_submodule = eisenstein_submodule.EisensteinSubmodule(self)
        return self.__eisenstein_submodule

    def new_submodule(self, p=None):
        try:
            return self.__new_submodule[p]
        except AttributeError:
            self.__new_submodule = {}
        except KeyError:
            pass
        if not p is None:
            p = rings.Integer(p)
            if not p.is_prime():
                raise ValueError, "p (=%s) must be a prime or None."%p

        if p is None:
            M = self._full_new_submodule()
            self.__new_submodule[None] = M
            return M
        else:
            M = self._new_submodule(p)
            self.__new_submodule[p] = M
            return M



    def _full_new_submodule(self):
        s = self._dim_new_cuspidal()
        e = self._dim_new_eisenstein()
        d = self._dim_cuspidal()
        B = range(s) + range(d, d+e)
        V = self.module()
        W = V.submodule([V.gen(i) for i in B])
        return submodule.ModularFormsSubmodule(self, W)

    def _new_submodule(self, p):
        raise NotImplementedError

    def _q_expansion(self, element, prec):
        B = self.q_expansion_basis(prec)
        f = self._q_expansion_zero()
        for i in range(len(element)):
            if element[i] != 0:
                f += element[i] * B[i]
        return f


    ####################################################################
    # Computations of Dimensions
    ####################################################################
    def _dim_cuspidal(self):
        try:
            return self.__the_dim_cuspidal
        except AttributeError:
            self.__the_dim_cuspidal = dims.dimension_cusp_forms(self.group(), self.weight())
        return self.__the_dim_cuspidal

    def _dim_eisenstein(self):
        try:
            return self.__the_dim_eisenstein
        except AttributeError:
            if self.weight() == 1:
                self.__the_dim_eisenstein = len(self.eisenstein_params())
            else:
                self.__the_dim_eisenstein = dims.dimension_eis(self.group(), self.weight())
        return self.__the_dim_eisenstein

    def _dim_new_cuspidal(self):
        try:
            return self.__the_dim_new_cuspidal
        except AttributeError:
            self.__the_dim_new_cuspidal = dims.dimension_new_cusp_forms_group(
                self.group(), self.weight())
        return self.__the_dim_new_cuspidal

    def _dim_new_eisenstein(self):
        try:
            return self.__the_dim_new_eisenstein
        except AttributeError:
            if isinstance(self.group(), congroup.Gamma0) and self.weight() == 2:
                if rings.is_prime(self.level()):
                    d = 1
                else:
                    d = 0
            else:
                E = self.eisenstein_series()
                d = len([g for g in E if g.new_level() == self.level()])
            self.__the_dim_new_eisenstein = d
        return self.__the_dim_new_cuspidal


    ####################################################################
    # Computations of all Eisenstein series in self
    ####################################################################

    def eisenstein_params(self):
        try:
            return self.__eisenstein_params
        except AttributeError:
            eps = self.character()
            if eps == None:
                if isinstance(self.group(), congroup.Gamma1):
                    eps = self.level()
                else:
                    raise NotImplementedError
            params = eis_series.compute_eisenstein_params(eps, self.weight())
            self.__eisenstein_params = params
        return self.__eisenstein_params

    def eisenstein_series(self):
        """
        Return all Eisenstein series associated to this space.

            sage: ModularForms(27,2).eisenstein_series()
            [
            q^3 + O(q^6),
            q - 3*q^2 + 7*q^4 - 6*q^5 + O(q^6),
            1/12 + q + 3*q^2 + q^3 + 7*q^4 + 6*q^5 + O(q^6),
            1/3 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6),
            13/12 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6)
            ]

            sage: ModularForms(Gamma1(5),3).eisenstein_series()
            [
            -1/5*zeta4 - 2/5 + q + (4*zeta4 + 1)*q^2 + (-9*zeta4 + 1)*q^3 + (4*zeta4 - 15)*q^4 + O(q^6),
            q + (zeta4 + 4)*q^2 + (-zeta4 + 9)*q^3 + (4*zeta4 + 15)*q^4 + O(q^6),
            1/5*zeta4 - 2/5 + q + (-4*zeta4 + 1)*q^2 + (9*zeta4 + 1)*q^3 + (-4*zeta4 - 15)*q^4 + O(q^6),
            q + (-zeta4 + 4)*q^2 + (zeta4 + 9)*q^3 + (-4*zeta4 + 15)*q^4 + O(q^6)
            ]

            sage: eps = DirichletGroup(13).0^2
            sage: ModularForms(eps,2).eisenstein_series()
            [
            7/13*zeta6 - 18/13 + q + (-2*zeta6 + 3)*q^2 + (3*zeta6 - 2)*q^3 + (-6*zeta6 + 3)*q^4 + -4*q^5 + O(q^6),
            q + (zeta6 + 2)*q^2 + (-zeta6 + 3)*q^3 + (3*zeta6 + 3)*q^4 + 4*q^5 + O(q^6)
            ]
        """
        return self.eisenstein_submodule().eisenstein_series()

    def _compute_q_expansion_basis(self, prec):
        """
        EXAMPLES:


        """
        S = self.cuspidal_submodule()
        E = self.eisenstein_submodule()
        B_S = S._compute_q_expansion_basis(prec)
        B_E = E._compute_q_expansion_basis(prec)
        return B_S + B_E

