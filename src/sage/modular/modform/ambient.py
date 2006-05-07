"""
Modular Forms

[[overview of what is here.]]

EXAMPLES:

"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

# system packages
import math
import weakref

# SAGE packages
import sage.rings.arith as arith
import sage.modular.congroup as congroup
import sage.misc.db as db
import sage.modular.dims as dims
import sage.modular.dirichlet as dirichlet
import sage.modular.hecke.all as hecke
import sage.misc.misc as misc
import sage.modular.modsym as modsym
import sage.modules.free_module as free_module
import sage.modules.free_module_element as free_module_element
import sage.rings.all as rings

from sage.misc.all import latex

import defaults

import space

class ModularFormsAmbient(space.ModularFormsSpace):
    """
    An ambient space of modular forms.
    """
    def __init__(self, group, weight, ring):
        if not isinstance(group, congroup.CongruenceSubgroup):
            group = congroup.Gamma0(group)

        if isinstance(group, congroup.Gamma0):
            character = dirichlet.TrivialCharacter(group.level(), base_field)
        else:
            character = None

        space.ModularFormsSpace.__init__(self, group, weight, character, base_field)

    def _repr_(self):
        return "Space of modular forms on %s of weight %s and dimension %s over %s"%(
                self.group(), self.weight(), self.dimension(), self.base_field())

    def change_ring(self, ring):
        import constructor
        return constructor.ModularForms(self.group(), self.weight(), ring)

    def dimension(self):
        if hasattr(self, "__dimension"): return self.__dimension
        self.__dimension = self.dim_eisenstein() + self.dim_cuspidal()
        return self.__dimension

    def ambient_space(self):
        return self

    def is_ambient(self):
        return True

    def modular_symbols(self):
        if hasattr(self, "__modsym"): return self.__modsym
        self.__modsym = modsym.ModularSymbols(self.weight(), self.group(),
                                                  self.base_field())
        return self.__modsym

    def vector_space(self):
        if hasattr(self, "__vector_space"): return self.__vector_space
        self.__vector_space = free_module.VectorSpace(self.base_field(),
                                                 self.dimension())
        return self.__vector_space

    def prec(self, set=None):
        """
        Set or get default initial precision for printing modular forms.
        """
        if set == None:
            if hasattr(self, "__prec"): return self.__prec
            self.__prec = defaults.DEFAULT_PRECISION
            return self.__prec
        self.__prec = set

    ####################################################################
    # Computation of q-expansions
    ####################################################################
    def qexp(self, vector, prec):
        """
        Compute the $q$-expansion to precision prec of the linear
        combination of the basis for this space given by the vector.
        """
        R = rings.PowerSeriesRing(self.base_field(), name='q')
        return rings.q
        raise NotImplementedError


    ####################################################################
    # Computation of Special Submodules
    ####################################################################

    def cuspidal_submodule(self):
        try:
            return self.__cuspidal_submodule
        except AttributeError:
            n = self.dim_cuspidal()
            V = self.vector_space()
            W = V.submodule([V.gen(i) for i in range(n)])
            self.__cuspidal_submodule = ModularFormsSubmodule(self, W)
        return self.__cuspidal_submodule

    def eisenstein_submodule(self):
        try:
            return self.__eisenstein_submodule
        except AttributeError:
            n = self.dim_cuspidal()
            d = self.dimension()
            V = self.vector_space()
            W = V.submodule([V.gen(i) for i in range(n, d)])
            self.__eisenstein_submodule = ModularFormsSubmodule(self, W)
        return self.__eisenstein_submodule

    def new_submodule(self):
        try:
            return self.__new_submodule
        except AttributeError:
            s = self.dim_new_cuspidal()
            e = self.dim_new_eisenstein()
            d = self.dim_cuspidal()
            B = range(s) + range(d, d+e)
            V = self.vector_space()
            W = V.submodule([V.gen(i) for i in B])
            self.__new_submodule = ModularFormsSubmodule(self, W)
        return self.__new_submodule


    ####################################################################
    # Computations of Dimensions
    ####################################################################
    def dim_cuspidal(self):
        try:
            return self.__the_dim_cuspidal
        except AttributeError:
            self.__the_dim_cuspidal = dims.dimension_cusp_forms(self.group(), self.weight())
        return self.__the_dim_cuspidal

    def dim_eisenstein(self):
        try:
            return self.__the_dim_eisenstein
        except AttributeError:
            if self.weight() == 1:
                self.__the_dim_eisenstein = len(self.eisenstein_params())
            else:
                self.__the_dim_eisenstein = dims.dimension_eis(self.group(), self.weight())
        return self.__the_dim_eisenstein

    def dim_new_cuspidal(self):
        try:
            return self.__the_dim_new_cuspidal
        except AttributeError:
            self.__the_dim_new_cuspidal = dims.dimension_new_cusp_forms_group(
                self.group(), self.weight())
        return self.__the_dim_new_cuspidal

    def dim_new_eisenstein(self):
        try:
            return self.__the_dim_new_eisenstein
        except AttributeError:
            if isinstance(self.group(), congroup.Gamma0) and self.weight() == 2:
                if arith.is_prime(self.level()):
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
            params = compute_eisenstein_params(eps, self.weight())
            self.__eisenstein_params = params
        return self.__eisenstein_params

    def eisenstein_series(self):
        try:
            return self.__eisenstein_series
        except AttributeError:
            pass
        time = misc.verbose("Finding eisenstein series.")

        misc.verbose("Finished finding eisenstein series.", time)

        self.__eisenstein_series = []
        c = self.dim_cuspidal()

        V = self.vector_space()
        params = self.eisenstein_params()
        assert V.dimension() - c == len(params)

        for i in range(len(params)):
            chi, psi, t = params[i]
            if chi.base_ring() != self.base_field():
                F = chi.base_ring()
                M = self.change_ring(F)
                v = V.change_ring(F).gen(c+i)
            else:
                M = self
                v = V.gen(c+i)
            E = EisensteinSeries(M, v, t, chi, psi)
            self.__eisenstein_series.append(E)
        return self.__eisenstein_series

class ModularFormsAmbient_g0_Q(ModularFormsAmbient):
    pass

class ModularFormsAmbient_g1_Q(ModularFormsAmbient):
    pass

class ModularFormsAmbient_eps(ModularFormsAmbient):
    """
    A space of modular forms with character.
    """
    def __init__(self, character, weight=2, base_field=None):
        """
        weight -- int
        character -- dirichlet.DirichletCharacter
        base_field -- base field
        """
        if base_field==None: base_field=character.base_ring()
        if character.base_ring() != base_field:
            character = character.change_ring(base_field)
        self.__weight, self.__character = weight, character
        group = congroup.Gamma1(character.modulus())
        base_field = character.base_ring()
        ModularFormsSpace.__init__(self, group, weight, character, base_field)

    def __repr__(self):
        return "Space of Modular Forms with character %s, weight %s, and dimension %s over %s"%(
                    self.character(), self.weight(), self.dimension(), self.base_field())

    def change_ring(self, F):
        return ModularFormsWithCharacter(self.weight(), self.character(), F)

    def modular_symbols(self):
        try:
            return self.__modsym
        except AttributeError:
            self.__modsym = modsym.ModularSymbolsWithCharacter(\
                self.character(), self.weight(), self.base_field())
        return self.__modsym

    ####################################################################
    # Computations of Dimensions
    ####################################################################
    def dim_cuspidal(self):
        try:
            return self.__the_dim_cuspidal
        except AttributeError:
            self.__the_dim_cuspidal = dims.dimension_cusp_forms_eps(
                self.character(), self.weight())
        return self.__the_dim_cuspidal

    def dim_eisenstein(self):
        try:
            return self.__the_dim_eisenstein
        except AttributeError:
            if self.character().is_trivial():
                edim = dims.dimension_eis(self.group(), self.weight())
            else:
                # calling self.eisenstein_series() is guaranteed to have the
                # side affect of setting self.__the_dim_eisenstein...
                self.eisenstein_series()  # todo: more direct way?
                edim = len(self.eisenstein_params())
            self.__the_dim_eisenstein = edim
        return self.__the_dim_eisenstein

    def dim_new_cuspidal(self):
        try:
            return self.__the_dim_new_cuspidal
        except AttributeError:
            self.__the_dim_new_cuspidal = dims.dimension_new_cusp_forms(
                self.character(), self.weight())
        return self.__the_dim_new_cuspidal

    def dim_new_eisenstein(self):
        try:
            return self.__the_dim_new_eisenstein
        except AttributeError:
            if self.character().is_trivial() and self.weight() == 2:
                if arith.is_prime(self.level()):
                    d = 1
                else:
                    d = 0
            else:
                E = self.eisenstein_series()
                d = len([g for g in E if g.new_level() == self.level()])
            self.__the_dim_new_eisenstein = d
        return self.__the_dim_new_cuspidal


