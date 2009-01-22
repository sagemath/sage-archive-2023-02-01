"""
Submodule of a Hecke module.
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

import sage.structure.factorization
import sage.rings.arith as arith
import sage.misc.misc as misc

import sage.modules.all

import module
import ambient_module

from sage.rings.polynomial.polynomial_ring import polygen

def is_HeckeSubmodule(x):
    return isinstance(x, HeckeSubmodule)

class HeckeSubmodule(module.HeckeModule_free_module):
    """
    Submodule of a Hecke module.
    """
    def __init__(self, ambient, submodule, dual_free_module=None, check=True):
        if not isinstance(ambient, ambient_module.AmbientHeckeModule):
            raise TypeError, "ambient must be an ambient Hecke module"
        if not sage.modules.all.is_FreeModule(submodule):
            raise TypeError, "submodule must be a free module"
        if not submodule.is_submodule(ambient.free_module()):
            raise ValueError, "submodule must be a submodule of the ambient free module"

        if check:
            if not ambient._is_hecke_equivariant_free_module(submodule):
                raise ValueError, "The submodule must be invariant under all Hecke operators."

        self.__ambient = ambient
        self.__submodule = submodule
        module.HeckeModule_free_module.__init__(self,
                                  ambient.base_ring(), ambient.level(), ambient.weight())
        if not (dual_free_module is None):
            if not sage.modules.all.is_FreeModule(dual_free_module):
                raise TypeError, "dual_free_module must be a free module"
            if dual_free_module.rank () != submodule.rank():
                raise ArithmeticError, "dual_free_module must have the same rank as submodule"
            self.__dual_free_module = dual_free_module


    def _repr_(self):
        return "Rank %s submodule of a Hecke module of level %s"%(
                      self.rank(), self.level())

    def __add__(self, other):
        if not isinstance(other, module.HeckeModule_free_module):
            raise TypeError, "other (=%s) must be a Hecke module."%other
        if self.ambient() != other.ambient():
            raise ArithmeticError, "Sum only defined for submodules of a common ambient space."
        if other.is_ambient():
            return other
        # Neither is ambient
        M = self.free_module() + other.free_module()
        return self.ambient().submodule(M, check=False)

    def __call__(self, x, check=True):
        """
        Coerce x into the ambient module and checks that x is in this submodule.

        EXAMPLES:
            sage: M = ModularSymbols(37)
            sage: S = M.cuspidal_submodule()
            sage: M([0,oo])
            -(1,0)
            sage: S([0,oo])
            Traceback (most recent call last):
            ...
            TypeError: x does not coerce to an element of this Hecke module
            sage: S([-1/23,0])
            (1,23)
        """
        z = self.ambient_hecke_module()(x)
        if check:
            if not z.element() in self.__submodule:
                raise TypeError, "x does not coerce to an element of this Hecke module"
        return z

    def __cmp__(self, other):
        if not isinstance(other, module.HeckeModule_free_module) or self.ambient() != other.ambient():
            return -1
        return cmp(self.free_module(), other.free_module())

    ################################
    # Semi-Private functions
    ################################
    def _compute_dual_hecke_matrix(self, n):
        A = self.ambient_hecke_module().dual_hecke_matrix(n)
        check =  arith.gcd(self.level(), n) != 1
        return A.restrict(self.dual_free_module(), check=check)

    def _compute_hecke_matrix(self, n):
        A = self.ambient_hecke_module().hecke_matrix(n)
        check = arith.gcd(self.level(), n) != 1
        return A.restrict(self.free_module(), check=check)

    def _compute_atkin_lehner_matrix(self, d):
        A = self.ambient_hecke_module()._compute_atkin_lehner_matrix(d)
        return A.restrict(self.free_module(), check=True)

    def _set_dual_free_module(self, V):
        if V.degree() != self.ambient_hecke_module().rank():
            raise ArithmeticError, "The degree of V must equal the rank of the ambient space."
        if V.rank() != self.rank():
            raise ArithmeticError, "The rank of V must equal the rank of self."
        self.__dual_free_module = V

    def _set_dual_free_module_from_nonembedded_module(self, V):
        """
        INPUT:
            V -- submodule of ambient free module of the same rank as the
                 rank of self.
        OUTPUT:
            Hecke submodule of self
        """
        M_V = V.matrix()
        E   = self.dual_free_module()
        M_E = E.matrix()
        A   = M_Vdual * M_E
        self.__dual_free_module = A.row_space()



    ################################
    # Public functions
    ################################

    def ambient_hecke_module(self):
        return self.__ambient

    def ambient(self):
        return self.__ambient

    def complement(self, bound=None):
        """
        Return the largest Hecke-stable complement of this space.

        EXAMPLES:
            sage: M = ModularSymbols(15, 6).cuspidal_subspace()
            sage: M.complement()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 20 for Gamma_0(15) of weight 6 with sign 0 over Rational Field
            sage: E = EllipticCurve("128a")
            sage: ME = E.modular_symbol_space()
            sage: ME.complement()
            Modular Symbols subspace of dimension 17 of Modular Symbols space of dimension 18 for Gamma_0(128) of weight 2 with sign 1 over Rational Field
        """
        try:
            return self.__complement
        except AttributeError:
            pass

        if self.is_ambient():
            return self.ambient_hecke_module().zero_submodule()

        if self.is_zero():
            return self.ambient_hecke_module()

        if self.is_full_hecke_module():
            anemic = False
        else:
            anemic = True

        # TODO: optimize in some cases by computing image of complementary factor instead of kernel...?
        misc.verbose("computing")
        N = self.level()
        A = self.ambient_hecke_module()
        V = A.free_module()
        p = 2
        if bound is None:
            bound = A.hecke_bound()
        while True:
            if anemic:
                while N % p == 0: p = arith.next_prime(p)
            misc.verbose("using T_%s"%p)
            f = self.hecke_polynomial(p)
            T = A.hecke_matrix(p)
            g = T.charpoly('x')
            V = T.kernel_on(V, poly=g//f, check=False)
            if V.rank() + self.rank() <= A.rank():
                break
            p = arith.next_prime(p)
            if p > bound:  # to avoid computing hecke bound unless necessary
                break

        if V.rank() + self.rank() == A.rank():
            C = A.submodule(V, check=False)
            self.__complement = C
            return C

        # first attempt to compute the complement failed, we now try
        # the following naive approach: decompose the ambient space,
        # decompose self, and sum the pieces of ambient that are not
        # subspaces of self
        misc.verbose("falling back on naive algorithm")
        D = A.decomposition()
        C = A.zero_submodule()
        for X in D:
            if self.intersection(X).dimension() == 0:
                C = C + X
        if C.rank() + self.rank() == A.rank():
            self.__complement = C
            return C

        # failed miserably
        raise RuntimeError, "Computation of complementary space failed (cut down to rank %s, but should have cut down to rank %s)."%(V.rank(), A.rank()-self.rank())


    def degeneracy_map(self, level, t=1):
        """
        The t-th degeneracy map from self to the space of ambient
        modular symbols of the given level.  The level of self must be
        a divisor or multiple of level, and t must be a divisor of the
        quotient.

        INPUT:
            level -- int, the level of the codomain of the map (positive int).
            t  -- int, the parameter of the degeneracy map, i.e., the map is
                  related to $f(q)$ |--> $f(q^t)$.

        OUTPUT:
            A linear function from self to the space of modular symbols
            of given level with the same weight, character, sign,
            etc., as this space.

        EXAMPLES:
            sage: D = ModularSymbols(10,4).cuspidal_submodule().decomposition(); D
            [
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 10 for Gamma_0(10) of weight 4 with sign 0 over Rational Field,
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 10 for Gamma_0(10) of weight 4 with sign 0 over Rational Field
            ]
            sage: d = D[1].degeneracy_map(5); d
            Hecke module morphism defined by the matrix
            [   0    0   -1    1]
            [   0  1/2  3/2   -2]
            [   0   -1    1    0]
            [   0 -3/4 -1/4    1]
            Domain: Modular Symbols subspace of dimension 4 of Modular Symbols space ...
            Codomain: Modular Symbols space of dimension 4 for Gamma_0(5) of weight ...

            sage: d.rank()
            2
            sage: d.kernel()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 10 for Gamma_0(10) of weight 4 with sign 0 over Rational Field
            sage: d.image()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(5) of weight 4 with sign 0 over Rational Field
        """
        d = self.ambient_hecke_module().degeneracy_map(level, t)
        return d.restrict_domain(self)


    def dual_free_module(self, bound=None, anemic=True, use_star=True):
        r"""
        Compute embedded dual free module if possible.  In general
        this won't be possible, e.g., if this space is not Hecke
        equivariant, possibly if it is not cuspidal, or if the
        characteristic is not 0.  In all these cases we raise a
        RuntimeError exception.

        If use_star is True (which is the default), we also use the
        +/- eigenspaces for the star operator to find the dual free
        module of self. If the self does not have a star involution,
        use_star will automatically be set to True.

        EXAMPLES:
            sage: M = ModularSymbols(11, 2)
            sage: M.dual_free_module()
            Vector space of dimension 3 over Rational Field
            sage: Mpc = M.plus_submodule().cuspidal_submodule()
            sage: Mcp = M.cuspidal_submodule().plus_submodule()
            sage: Mcp.dual_free_module() == Mpc.dual_free_module()
            True
            sage: Mpc.dual_free_module()
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [  1 5/2   5]

            sage: M = ModularSymbols(35,2).cuspidal_submodule()
            sage: M.dual_free_module(use_star=False)
            Vector space of degree 9 and dimension 6 over Rational Field
            Basis matrix:
            [   1    0    0    0   -1    0    0    4   -2]
            [   0    1    0    0    0    0    0 -1/2  1/2]
            [   0    0    1    0    0    0    0 -1/2  1/2]
            [   0    0    0    1   -1    0    0    1    0]
            [   0    0    0    0    0    1    0   -2    1]
            [   0    0    0    0    0    0    1   -2    1]

            sage: M = ModularSymbols(40,2)
            sage: Mmc = M.minus_submodule().cuspidal_submodule()
            sage: Mcm = M.cuspidal_submodule().minus_submodule()
            sage: Mcm.dual_free_module() == Mmc.dual_free_module()
            True
            sage: Mcm.dual_free_module()
            Vector space of degree 13 and dimension 3 over Rational Field
            Basis matrix:
            [ 0  1  0  0  0  0  1  0 -1 -1  1 -1  0]
            [ 0  0  1  0 -1  0 -1  0  1  0  0  0  0]
            [ 0  0  0  0  0  1  1  0 -1  0  0  0  0]

            sage: M = ModularSymbols(43).cuspidal_submodule()
            sage: S = M[0].plus_submodule() + M[1].minus_submodule()
            sage: S.dual_free_module(use_star=False)
            Traceback (most recent call last):
            ...
            RuntimeError: Computation of embedded dual vector space failed (cut down to rank 6, but should have cut down to rank 3).
            sage: S.dual_free_module().dimension() == S.dimension()
            True
        """
        try:
            return self.__dual_free_module
        except AttributeError:
            pass

        misc.verbose("computing")

        A = self.ambient_hecke_module()

        if self.dimension() == 0:
            self.__dual_free_module = A.zero_submodule()
            return self.__dual_free_module

        if A.dimension() == self.dimension():
            self.__dual_free_module = A.free_module()
            return self.__dual_free_module

        # ALGORITHM: Compute the char poly of each Hecke operator on
        # the submodule, then use it to cut out a submodule of the
        # dual.  If the dimension cuts down to the dimension of self
        # terminate with success.  If it stays larger beyond the Sturm
        # bound, raise a RuntimeError exception.

        # In the case that the sign of self is not 1, we need to use
        # the star involution as well as the Hecke operators in order
        # to find the dual of self.
        #
        # Note that one needs to comment out the line caching the
        # result of this computation below in order to get meaningful
        # timings.

        # If the star involution doesn't make sense for self, then we
        # can't use it.
        if not hasattr(self, 'star_eigenvalues'):
            use_star = False

        if use_star:
            # If the star involution has both + and - eigenspaces on self,
            # then we compute the dual on each eigenspace, then put them
            # together.
            if len(self.star_eigenvalues()) == 2:
                V = self.plus_submodule(compute_dual = False).dual_free_module() + \
                    self.minus_submodule(compute_dual = False).dual_free_module()
                return V

            # At this point, we know that self is an eigenspace for star.
            V = A.sign_submodule(self.sign()).dual_free_module()
        else:
            V = A.free_module()

        N = self.level()
        p = 2
        if bound is None:
            bound = A.hecke_bound()
        while True:
            if anemic:
                while N % p == 0: p = arith.next_prime(p)
            misc.verbose("using T_%s"%p)
            f = self.hecke_polynomial(p)
            T = A.dual_hecke_matrix(p)
            V = T.kernel_on(V, poly=f, check=False)
            if V.dimension() <= self.dimension():
                break
            p = arith.next_prime(p)
            if p > bound:
                break

        if V.rank() == self.rank():
            self.__dual_free_module = V
            return V
        else:
            # Failed to reduce V to the appropriate dimension
            raise RuntimeError, "Computation of embedded dual vector space failed " + \
                  "(cut down to rank %s, but should have cut down to rank %s)."%(V.rank(), self.rank())


    def free_module(self):
        """
        Return the free module corresponding to self.

        EXAMPLES:
            sage: M = ModularSymbols(33,2).cuspidal_subspace() ; M
            Modular Symbols subspace of dimension 6 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
            sage: M.free_module()
            Vector space of degree 9 and dimension 6 over Rational Field
            Basis matrix:
            [ 0  1  0  0  0  0  0 -1  1]
            [ 0  0  1  0  0  0  0 -1  1]
            [ 0  0  0  1  0  0  0 -1  1]
            [ 0  0  0  0  1  0  0 -1  1]
            [ 0  0  0  0  0  1  0 -1  1]
            [ 0  0  0  0  0  0  1 -1  0]
        """
        return self.__submodule

    def module(self):
        r"""
        Alias for \code{self.free_module()}.

        EXAMPLES:
            sage: M = ModularSymbols(17,4).cuspidal_subspace()
            sage: M.free_module() is M.module()
            True
        """
        return self.free_module()

    def intersection(self, other):
        """
        Returns the intersection of self and other, which must both
        lie in a common ambient space of modular symbols.

        EXAMPLES:
            sage: M = ModularSymbols(43, sign=1)
            sage: A = M[0] + M[1]
            sage: B = M[1] + M[2]
            sage: A.dimension(), B.dimension()
            (2, 3)
            sage: C = A.intersection(B); C.dimension()
            1

        TESTS:
            sage: M = ModularSymbols(1,80)
            sage: M.plus_submodule().cuspidal_submodule().sign() # indirect doctest
            1
        """
        if self.ambient_hecke_module() != other.ambient_hecke_module():
            raise ArithmeticError, "Intersection only defined for subspaces of"\
                  + " a common ambient modular symbols space."
        if other.is_ambient():
            return self
        if self.is_ambient():
            return other

        # Neither is ambient
        V = self.free_module().intersection(other.free_module())
        M = self.ambient_hecke_module().submodule(V,check=False)

        ## if sign is nonzero, the intersection will be, too
        if self.sign():
            M._set_sign(self.sign())
        elif other.sign():
            M._set_sign(other.sign())

        return M

    def is_ambient(self):
        r"""
        Return \code{True} if self is an ambient space of modular
        symbols.

        EXAMPLES:
            sage: M = ModularSymbols(17,4)
            sage: M.cuspidal_subspace().is_ambient()
            False
            sage: A = M.ambient_hecke_module()
            sage: S = A.submodule(A.basis())
            sage: sage.modular.hecke.submodule.HeckeSubmodule.is_ambient(S)
            True
        """
        return self.free_module() == self.ambient_hecke_module().free_module()

    def is_new(self, p=None):
        """
        Returns True if this Hecke module is p-new.  If p is None,
        returns True if it is new.
        """
        try:
            return self.__is_new[p]
        except AttributeError:
            self.__is_new = {}
        except KeyError:
            pass
        N = self.ambient_hecke_module().new_submodule(p)
        self.__is_new[p] = self.is_submodule(N)
        return self.__is_new[p]

    def is_old(self, p=None):
        """
        Returns True if this Hecke module is p-old.  If p is None,
        returns True if it is old.
        """
        try:
            return self.__is_old[p]
        except AttributeError:
            self.__is_old = {}
        except KeyError:
            pass
        O = self.ambient_hecke_module().old_submodule(p)
        self.__is_old[p] = self.is_submodule(O)
        return self.__is_old[p]

    def is_submodule(self, V):
        """
        Returns True if and only if self is a submodule of V.
        """
        if not isinstance(V, module.HeckeModule_free_module):
            return False
        return self.ambient_hecke_module() == V.ambient_hecke_module() and \
               self.free_module().is_subspace(V.free_module())

    def linear_combination_of_basis(self, v):
        """
        Return the linear combination of the basis of self given by
        the entries of v.
        """
        x = self.free_module().linear_combination_of_basis(v)
        return self.__ambient(x)

    def new_submodule(self, p=None):
        """
        Return the new or p-new submodule of this space of modular symbols.
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

        S = self.ambient_hecke_module().new_submodule(p)
        ns = S.intersection(self)
        if ns.rank() == self.rank():
            self.__is_new[p] = True
        ns.__is_new = {p:True}
        self.__new_submodule[p] = ns
        return ns

    def nonembedded_free_module(self):
        return self.free_module().nonembedded_free_module()

    def old_submodule(self, p=None):
        """
        Return the old or p-old submodule of this space of modular symbols.

        EXAMPLES:
        We compute the old and new submodules of $\sS_2(\Gamma_0(33))$.

            sage: M = ModularSymbols(33); S = M.cuspidal_submodule(); S
            Modular Symbols subspace of dimension 6 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
            sage: S.old_submodule()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
            sage: S.new_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
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

        S = self.ambient_hecke_module().old_submodule(p)
        os = S.intersection(self)
        if os.rank() == self.rank():
            self.__is_old[p] = True
        os.__is_old = {p:True}
        self.__old_submodule[p] = os
        return os

    def rank(self):
        return self.__submodule.rank()

    def dimension(self):
        return self.rank()

    def submodule(self, M, Mdual=None, check=True):
        """
        Construct a submodule of self from the embedded free module M.
        """
        if check:
            if not sage.modules.all.is_FreeModule(M):
                V = self.ambient_module().free_module()
                if isinstance(M, (list,tuple)):
                    M = V.span([V(x.element()) for x in M])
                else:
                    M = V.span(M)
            if not M.is_submodule(self.free_module()):
                raise TypeError, "M (=%s) must be a submodule of the free module (=%s) associated to this module."%(M, self.free_module())

        return self.ambient().submodule(M, Mdual, check=check)

    def submodule_from_nonembedded_module(self, V, Vdual=None, check=True):
        """
        INPUT:
            V -- submodule of ambient free module of the same rank as the
                 rank of self.
            check -- whether to check that V is Hecke equivariant.

        OUTPUT:
            Hecke submodule of self
        """
        E = self.free_module()
        M_V = V.matrix()
        M_E = E.matrix()
        # We encode the operation of taking the linear combinations of
        # the basis of E given by the basis of V as a single matrix
        # multiplication, since matrix multiplication is (presumed to be)
        # so fast, and their are asymptotically fast algorithms.
        A = M_V * M_E
        V = A.row_space()
        if not (Vdual is None):
            E = self.dual_free_module()
            M_Vdual = Vdual.matrix()
            M_E = E.matrix()
            A = M_Vdual * M_E
            Vdual = A.row_space()
        return self.ambient_hecke_module().submodule(V, Vdual, check=check)

