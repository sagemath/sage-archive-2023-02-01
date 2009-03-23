"""
Ambient Hecke modules
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

import degenmap
import module
import submodule

import sage.modules.all

import sage.rings.all

import sage.misc.misc as misc

import sage.rings.arith as arith

import sage.matrix.matrix_space as matrix_space
from   sage.matrix.constructor import matrix

from sage.modular.arithgroup.all import Gamma0 # for Sturm bound

def is_AmbientHeckeModule(x):
    return isinstance(x, AmbientHeckeModule)

class AmbientHeckeModule(module.HeckeModule_free_module):
    """
    Ambient Hecke module.
    """
    def __init__(self, base_ring, rank, level, weight):
        rank = sage.rings.all.Integer(rank)
        if rank < 0:
            raise ValueError, "rank (=%s) must be nonnegative"%rank
        self.__rank = rank
        module.HeckeModule_free_module.__init__(self, base_ring, level, weight)

    def rank(self):
        """
        Return the rank of this ambient Hecke module.

        OUTPUT:
            Integer

        EXAMPLES::

            sage: M = sage.modular.hecke.ambient_module.AmbientHeckeModule(QQ, 3, 11, 2); M
            Ambient Hecke module of rank Rational Field over 3
            sage: M.rank()
            3
        """
        return self.__rank

    def __add__(self, other):
        if not isinstance(other, module.HeckeModule_free_module):
            raise TypeError, "other (=%s) must be a Hecke module."%other
        if other.ambient_hecke_module() == self:
            return self
        raise ArithmeticError, "Sum only defined for subspaces of a common ambient Hecke module."

    def _repr_(self):
        return "Ambient Hecke module of rank %s over %s"%(self.base_ring(), self.rank())


    def _degeneracy_raising_matrix(self, N):
        """
        Matrix of the degeneracy map (with t = 1) to level N, where N is a
        multiple of the level.
        """
        raise NotImplementedError

    def _degeneracy_lowering_matrix(self, N, t):
        """
        Matrix of the degeneracy map of index t to level N, where N is a
        divisor of the level.
        """
        raise NotImplementedError

    def _hecke_image_of_ith_basis_element(self, n, i):
        """
        Return the image under the Hecke operator T_n of the i-th basis
        element.
        """
        return self.hecke_operator(n)(self.gen(i))


    def _set_dual_free_module(self, V):
        pass  # setting dual free module of ambient space is not necessary


    def ambient_hecke_module(self):
        return self

    def complement(self):
        """
        Return the largest Hecke-stable complement of this space.
        EXAMPLES::

            sage: M=ModularSymbols(11,2,1)
            sage: M
            Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: M.complement()
            Modular Symbols subspace of dimension 0 of Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: C=M.cuspidal_subspace()
            sage: C
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: C.complement()
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
        """
        return self.zero_submodule()

    def decomposition_matrix(self):
        r"""
        Returns the matrix whose columns form a basis for the canonical
        sorted decomposition of self coming from the Hecke operators.

        If the simple factors are `D_0, \ldots, D_n`, then the
        first few columns are an echelonized basis for `D_0`, the
        next an echelonized basis for `D_1`, the next for
        `D_2`, etc.
        """
        try:
            return self.__decomposition_matrix_cache
        except AttributeError:
            rows = []
            for A in self.decomposition():
                for x in A.basis():
                    rows.append(x.list())
            A = matrix_space.MatrixSpace(self.base_ring(),self.rank())(rows)
            self.__decomposition_matrix_cache = A
            return self.__decomposition_matrix_cache

    def decomposition_matrix_inverse(self):
        """
        Returns the inverse of the matrix returned by
        decomposition_matrix().
        """
        try:
            return self.__decomposition_matrix_inverse_cache
        except AttributeError:
            self.__decomposition_matrix_inverse_cache = ~self.decomposition_matrix()
            return self.__decomposition_matrix_inverse_cache

    def degeneracy_map(self, level, t=1):
        """
        The t-th degeneracy map from self to the corresponding Hecke module
        of the given level. The level of self must be a divisor or multiple
        of level, and t must be a divisor of the quotient.

        INPUT:


        -  ``level`` - int, the level of the codomain of the
           map (positive int).

        -  ``t`` - int, the parameter of the degeneracy map,
           i.e., the map is related to `f(q)` - `f(q^t)`.


        OUTPUT: A morphism from self to corresponding the Hecke module of
        given level.

        EXAMPLES::

            sage: M = ModularSymbols(11,sign=1)
            sage: d1 = M.degeneracy_map(33); d1
            Hecke module morphism degeneracy map corresponding to f(q) |--> f(q) defined by the matrix
            (not printing 2 x 6 matrix)
            Domain: Modular Symbols space of dimension 2 for Gamma_0(11) of weight ...
            Codomain: Modular Symbols space of dimension 6 for Gamma_0(33) of weight ...
            sage: M.degeneracy_map(33,3).matrix()
            [ 3  2  2  0 -2  1]
            [ 0  2  0 -2  0  0]
            sage: M = ModularSymbols(33,sign=1)
            sage: d2 = M.degeneracy_map(11); d2.matrix()
            [  1   0]
            [  0 1/2]
            [  0  -1]
            [  0   1]
            [ -1   0]
            [ -1   0]
            sage: (d2*d1).matrix()
            [4 0]
            [0 4]

        ::

            sage: M = ModularSymbols(3,12,sign=1)
            sage: M.degeneracy_map(1)
            Hecke module morphism degeneracy map corresponding to f(q) |--> f(q) defined by the matrix
            [1 0]
            [0 0]
            [0 1]
            [0 1]
            [0 1]
            Domain: Modular Symbols space of dimension 5 for Gamma_0(3) of weight ...
            Codomain: Modular Symbols space of dimension 2 for Gamma_0(1) of weight ...

        ::

            sage: S = M.cuspidal_submodule()
            sage: S.degeneracy_map(1)
            Hecke module morphism defined by the matrix
            [1 0]
            [0 0]
            [0 0]
            Domain: Modular Symbols subspace of dimension 3 of Modular Symbols space ...
            Codomain: Modular Symbols space of dimension 2 for Gamma_0(1) of weight ...

        ::

            sage: D = ModularSymbols(10,4).cuspidal_submodule().decomposition()
            sage: D
            [
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 10 for Gamma_0(10) of weight 4 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 10 for Gamma_0(10) of weight 4 with sign 0 over Rational Field
            ]
            sage: D[1].degeneracy_map(5)
            Hecke module morphism defined by the matrix
            [   0    0   -1    1]
            [   0  1/2  3/2   -2]
            [   0   -1    1    0]
            [   0 -3/4 -1/4    1]
            Domain: Modular Symbols subspace of dimension 4 of Modular Symbols space ...
            Codomain: Modular Symbols space of dimension 4 for Gamma_0(5) of weight ...
        """
        level = int(level); t = int(t)

        err = False
        if self.level() % level == 0:
            quo = self.level() // level
            if quo % t != 0:
                err = True
        elif level % self.level() == 0:
            quo = level // self.level()
            if quo % t != 0:
                err = True
        else:
            err = True
        if err:
            raise ValueError, ("The level of self (=%s) must be a divisor or multiple of " + \
                               "level (=%s), and t (=%s) must be a divisor of the quotient.")%\
                               (self.level(), level, t)

        key = (level,t)
        try:
            self._degeneracy_maps
        except AttributeError:
            self._degeneracy_maps = {}

        if self._degeneracy_maps.has_key(key):
            return self._degeneracy_maps[key]

        eps = self.character()
        if not (eps is None) and level % eps.conductor() != 0:
            raise ArithmeticError, "The conductor of the character of this space " + \
                  "(=%s) must be divisible by the level (=%s)."%\
                  (eps.conductor(), level)

        M = self.hecke_module_of_level(level)

        if M.rank() == 0:

            A = matrix_space.MatrixSpace(self.base_ring(), self.rank(),0)(0)

        elif self.level() % level == 0:  # lower the level

            A = self._degeneracy_lowering_matrix(level, t)

        elif level % self.level() == 0:  # raise the level
            if t != 1:

                # use Hecke operator and t=1 case.
                d1 = self.degeneracy_map(level, 1).matrix()
                T = M.hecke_matrix(t)
                A = (~self.base_ring()(t)) * d1 * T

            else:

                A = self._degeneracy_raising_matrix(level)

        d = degenmap.DegeneracyMap(A, self, M, t)
        self._degeneracy_maps[key] = d
        return d

    def dual_free_module(self):
        return self.free_module()

    def fcp(self, n, var='x'):
        """
        Returns the factorization of the characteristic polynomial of the
        Hecke operator `T_n` of index `n`.

        INPUT:


        -  ``ModularSymbols self`` - space of modular symbols
           invariant under the Hecke operator of index n.

        -  ``int n`` - a positive integer.

        -  ``var`` - variable of polynomiall


        OUTPUT:


        -  ``list`` - list of the pairs (g,e), where g is an
           irreducible factor of the characteristic polynomial of T_n, and e
           is its multiplicity.


        EXAMPLES::

            sage: m = ModularSymbols(23, 2, sign=1)
            sage: m.fcp(2)
            (x - 3) * (x^2 + x - 1)
            sage: m.hecke_operator(2).charpoly('x').factor()
            (x - 3) * (x^2 + x - 1)
        """
        n = int(n)
        if n <= 0:
            raise ArithmeticError, "n (=%s) must be positive"%n
        return self.hecke_operator(n).fcp(var)

    def free_module(self):
        """
        Return the free module underlying this ambient Hecke module.
        """
        try:
            return self.__free_module
        except AttributeError:
            M = sage.modules.all.FreeModule(self.base_ring(), self.rank())
            self.__free_module = M
            return M

    def hecke_bound(self):
        """
        Return an integer B such that the Hecke operators `T_n`,
        for `n\leq B`, generate the full Hecke algebra as a module
        over the base ring. Note that we include the `n` with
        `n` not coprime to the level.
        """
        misc.verbose("WARNING: ambient.py -- hecke_bound; returning unproven guess.")
        return Gamma0(self.level()).sturm_bound(self.weight()) + 2*self.group().dimension_eis(self.weight()) + 5

    def hecke_module_of_level(self, level):
        raise NotImplementedError


    def hecke_images(self, i, v):
        """
        Return images of the `i`-th standard basis vector under the
        Hecke operators `T_p` for all integers in `v`.

        INPUT:


        -  ``i`` - nonnegative integer

        -  ``v`` - a list of positive integer


        OUTPUT:


        -  ``matrix`` - whose rows are the Hecke images


        EXAMPLES::

            sage: M = ModularSymbols(DirichletGroup(13).0, 3)
            sage: M.T(2)(M.0).element()
            (zeta12 + 4, 0, -1, 1)
            sage: M.hecke_images(0, [1,2])
            [         1          0          0          0]
            [zeta12 + 4          0         -1          1]
        """
        try:
            return self._hecke_images(i, v)
        except (AttributeError, NotImplementedError):
            pass
        # Use slow generic algorithm
        x = self.gen(i)
        X = [self.hecke_operator(n).apply_sparse(x).element() for n in v]
        return matrix(self.base_ring(), X)

    def intersection(self, other):
        """
        Returns the intersection of self and other, which must both lie in
        a common ambient space of modular symbols.

        EXAMPLES::

            sage: M = ModularSymbols(43, sign=1)
            sage: A = M[0] + M[1]
            sage: B = M[1] + M[2]
            sage: A.rank(), B.rank()
            (2, 3)
            sage: C = A.intersection(B); C.rank()  # TODO
            1
        """
        if not isinstance(other, module.HeckeModule_free_module):
            raise TypeError, "other (=%s) must be a Hecke module."%other
        if self.ambient_hecke_module() != other.ambient_hecke_module():
            raise ArithmeticError, "Intersection only defined for subspaces of a common ambient Hecke module."
        return other  # since self is ambient, so the intersection must equal other.

    def is_ambient(self):
        r"""
        Returns True if and only if self is an ambient Hecke module.

        .. warning::

           self can only be ambient by being of type
           AmbientHeckeModule.

           For example, decomposing a simple ambient space yields a
           single factor, and that factor is *not* considered an
           ambient space.

        EXAMPLES::

            sage: m = ModularSymbols(10)
            sage: m.is_ambient()
            True

        ::

            sage: a = m[0]  # the unique simple factor
            sage: a == m
            True
            sage: a.is_ambient()
            False
        """
        return True

    def is_full_hecke_module(self, compute=True):
        """
        Returns True if this space is invariant under the action of all
        Hecke operators, even those that divide the level.
        """
        return True

    def is_new(self, p=None):
        try:
            if self.__is_new.has_key(p):
                return self.__is_new[p]
        except AttributeError:
            pass
        self.new_submodule(p)
        return self.__is_new[p]

    def is_old(self, p=None):
        try:
            if self.__is_old.has_key(p):
                return self.__is_old[p]
        except AttributeError:
            pass
        self.old_submodule(p)
        return self.__is_old[p]

    def is_submodule(self, V):
        """
        Returns True if and only if self is a submodule of V.
        """
        if not isinstance(V, module.HeckeModule_free_module):
            raise TypeError, "V must be a Hecke module"
        if not V.is_ambient():
            return False
        return V.ambient_space() == self

    def linear_combination_of_basis(self, v):
        return self(v)

    def new_submodule(self, p=None):
        """
        Returns the new or p-new submodule of self.

        INPUT:


        -  ``p`` - (default: None); if not None, return only
           the p-new submodule.


        OUTPUT: the new or p-new submodule of self

        EXAMPLES::

            sage: m = ModularSymbols(33); m.rank()
            9
            sage: m.new_submodule().rank()
            3
            sage: m.new_submodule(3).rank()
            4
            sage: m.new_submodule(11).rank()
            8
        """
        try:
            if self.__is_new[p]:
                return self
        except AttributeError:
            self.__is_new = {}
        except KeyError:
            pass

        if self.rank() == 0:
            self.__is_new[p] = True
            return self
        try:
            return self.__new_submodule[p]
        except AttributeError:
            self.__new_submodule = {}
        except KeyError:
            pass

        # Construct the degeneracy map d.
        N = self.level()
        d = None
        eps = self.character()
        if eps == None:
            f = 1
        else:
            f = eps.conductor()
        if p == None:
            D = arith.prime_divisors(N)
        else:
            if N % p != 0:
                raise ValueError, "p must divide the level."
            D = [p]
        for q in D:
            if ((N//q) % f) == 0:
                NN = N//q
                d1 = self.degeneracy_map(NN,1).matrix()
                if d is None:
                    d = d1
                else:
                    d = d.augment(d1)
                d = d.augment(self.degeneracy_map(NN,q).matrix())
            #end if
        #end for
        if d is None or d == 0:
            self.__is_new[p] = True
            return self
        else:
            self.__is_new[p] = False
        ns = self.submodule(d.kernel(), check=False)
        ns.__is_new = {p:True}
        ns._is_full_hecke_module = True
        self.__new_submodule[p] = ns
        return ns

    def nonembedded_free_module(self):
        return self.free_module()

    def old_submodule(self, p=None):
        """
        Returns the old or p-old submodule of self.

        INPUT:


        -  ``p`` - (default: None); if not None, return only
           the p-old submodule.


        OUTPUT: the old or p-old submodule of self

        EXAMPLES::

            sage: m = ModularSymbols(33); m.rank()
            9
            sage: m.old_submodule().rank()
            7
            sage: m.old_submodule(3).rank()
            6
            sage: m.new_submodule(11).rank()
            8

        ::

            sage: e = DirichletGroup(16)([-1, 1])
            sage: M = ModularSymbols(e, 3, sign=1); M
            Modular Symbols space of dimension 4 and level 16, weight 3, character [-1, 1], sign 1, over Rational Field
            sage: M.old_submodule()
            Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 4 and level 16, weight 3, character [-1, 1], sign 1, over Rational Field
        """
        try:
            if self.__is_old[p]:
                return self
        except AttributeError:
            self.__is_old = {}
        except KeyError:
            pass

        if self.rank() == 0:
            self.__is_old[p] = True
            return self
        try:
            return self.__old_submodule[p]
        except AttributeError:
            self.__old_submodule = {}
        except KeyError:
            pass

        # Construct the degeneracy map d.
        N = self.level()
        d = None

        eps = self.character()
        if eps is None:
            f = 1
        else:
            f = eps.conductor()

        if p is None:
            D = arith.prime_divisors(N)
        else:
            if N % p != 0:
                raise ValueError, "p must divide the level."
            D = [p]

        for q in D:
            NN = N//q
            if NN % f == 0:
                M = self.hecke_module_of_level(NN)
                d1 = M.degeneracy_map(N, 1).matrix()
                if d is None:
                    d = d1
                else:
                    d = d.stack(d1)
                d = d.stack(M.degeneracy_map(N, q).matrix())
            #end if
        #end for
        os = self.submodule(d.image(), check=False)
        self.__is_old[p] = (os == self)

        os.__is_old = {p:True}
        os._is_full_hecke_module = True
        self.__old_submodule[p] = os
        return os


    def submodule(self, M, Mdual=None, check=True):
        """
        Return the Hecke submodule of self defined by the free module M.
        """
        if check:
            if not sage.modules.all.is_FreeModule(M):
                V = self.free_module()
                if isinstance(M, (list,tuple)):
                    M = V.span([V(x.element()) for x in M])
                else:
                    M = V.span(M)
            if not M.is_submodule(self.free_module()):
                raise TypeError, "M must be a submodule of the free module associated to this module."
            if M == self.free_module():
                return self
        return self._submodule_class()(self, M, Mdual, check=check)

    def _submodule_class(self):
        return submodule.HeckeSubmodule

    def submodule_from_nonembedded_module(self, V, Vdual=None, check=True):
        """
        INPUT:


        -  ``V`` - submodule of ambient free module of the same
           rank as the rank of self.

        -  ``Vdual`` - used to pass in dual submodule

        -  ``check`` - whether to check that submodule is Hecke
           equivariant


        OUTPUT: Hecke submodule of self
        """
        return self.submodule(V, Vdual, check=check)

    def submodule_generated_by_images(self, M):
        """
        Return the submodule of this ambient modular symbols space
        generated by the images under all degeneracy maps of M. The space M
        must have the same weight, sign, and group or character as this
        ambient space.
        """
        S = self.zero_submodule()
        if self.level() % M.level() == 0:
            D = arith.divisors(self.level() // M.level())
        elif M.level() % self.level() == 0:
            D = arith.divisors(M.level() // self.level())
        else:
            D = []
        for t in D:
            d = M.degeneracy_map(self.level(), t)
            if d.codomain() != self:
                raise ArithmeticError, "incompatible spaces of modular symbols"
            S += d.image()

        if self.is_full_hecke_module(compute=False):
            S._is_full_hecke_module = True

        return S


