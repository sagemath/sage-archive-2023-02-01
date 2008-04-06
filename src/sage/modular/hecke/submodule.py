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
            misc.verbose("using T_%s"%p)
            if anemic:
                while N % p == 0: p = arith.next_prime(p)
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
        else:
            # failed
            raise RuntimeError, "Computation of complementary space failed (cut down to rank %s, but should have cut down to rank %s)."%(V.rank(), self.rank())


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


    def dual_free_module(self, bound=None, anemic=True):
        """
        Compute embedded dual free module if possible.  In general
        this won't be possible, e.g., if this space is not Hecke
        equivariant, possibly if it is not cuspidal, or if the
        characteristic is not 0.  In all these cases we raise a
        RuntimeError exception.
        """
        try:
            return self.__dual_free_module
        except AttributeError:
            if self.dimension() == 0:
                self.__dual_free_module = self.ambient_hecke_module().dual_free_module().zero_submodule()
                return self.__dual_free_module

            # ALGORITHM: Compute the char poly of each Hecke operator on the
            # submodule, then use it to cut out a submodule of the dual.  If
            # the dimension cuts down to the dimension of self terminate
            # with success.  If it stays bigger beyond the bound (Sturm)
            # bound, raise a RuntimeError exception.
            misc.verbose("computing")
            N = self.level()
            A = self.ambient_hecke_module()
            if A.dimension() == self.dimension():
                self.__dual_free_module = A.free_module()
                return self.__dual_free_module
            V = A.free_module()
            p = 2
            if bound is None:
                bound = A.hecke_bound()
            while True:
                misc.verbose("using T_%s"%p)
                if anemic:
                    while N % p == 0: p = arith.next_prime(p)
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
                # failed
                raise RuntimeError, "Computation of embedded dual vector space failed " + \
                      "(cut down to rank %s, but should have cut down to rank %s)."%(V.rank(), self.rank())


    def free_module(self):
        return self.__submodule

    def module(self):
        return self.__submodule

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

