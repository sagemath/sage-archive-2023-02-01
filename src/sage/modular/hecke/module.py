"""
Hecke modules
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@ucsd.edu>
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


import sage.rings.arith as arith
import sage.misc.misc as misc
import sage.modules.module
import sage.categories.all
import sage.structure.factorization
from sage.structure.all import Sequence
import sage.matrix.matrix_space as matrix_space

import random

import algebra
import element
import homspace

from sage.modules.all import FreeModule

def is_HeckeModule(x):
    return isinstance(x, HeckeModule_generic)

class HeckeModule_generic(sage.modules.module.Module):
    """
    A very general Hecke module.

    All Hecke module classes derive from this class---spaces of
    modular symbols (free modules), modular forms (finite-rank free
    modules), modular abelian varieties (infinitely divisible groups),
    torsion submodules of abelian varieties (finite groups), etc.
    """
    def __init__(self, base_ring, level):
        if not sage.rings.all.is_CommutativeRing(base_ring):
            raise TypeError, "base_ring must be commutative ring"
        self.__base_ring = base_ring

        level = int(level)
        if level <= 0:
            raise ValueError, "level (=%s) must be positive"%level
        self.__level = level

    def __hash__(self):
        return hash((self.__base_ring, self.__level))

    def __cmp__(self, other):
        if not isinstance(other, HeckeModule_generic):
            return -1
        return cmp((self.__level, self.__base_ring), (other.__level, other.__base_ring))

    def _compute_hecke_matrix_prime_power(self, n, p, r):
        # convert input arguments to int's.
        (n,p,r) = (int(n), int(p), int(r))
        if not arith.is_prime(p):
            raise ArithmeticError, "p must be a prime"
        # T_{p^r} := T_p * T_{p^{r-1}} - eps(p)p^{k-1} T_{p^{r-2}}.
        pow = p**(r-1)
        if not self._hecke_matrices.has_key(pow):
            # The following will force computation of T_{p^s}
            # for all s<=r-1, except possibly s=0.
            self._hecke_matrices[pow] = self._compute_hecke_matrix(pow)
        if not self._hecke_matrices.has_key(1):
            self._hecke_matrices[1] = self._compute_hecke_matrix(1)
        Tp = self._hecke_matrices[p]
        Tpr1 = self._hecke_matrices[pow]
        eps = self.character()
        if eps is None:
            raise NotImplementedError, "_compute_hecke_matrix_prime_power must be overloaded in a derived class"
        k = self.weight()
        Tpr2 = self._hecke_matrices[pow/p]
        return Tp*Tpr1 - eps(p)*(p**(k-1)) * Tpr2

    def _compute_hecke_matrix_general_product(self, n, F):
        n = int(n)
        prod = None
        for p, r in F:
            pow = int(p**r)
            if not self._hecke_matrices.has_key(pow):
                self._hecke_matrices[pow] = self._compute_hecke_matrix(pow)
            if prod == None:
                prod = self._hecke_matrices[pow]
            else:
                prod *= self._hecke_matrices[pow]
        return prod

    def _compute_dual_hecke_matrix(self, n):
        return self.hecke_matrix(n).transpose()

    def _compute_hecke_matrix(self, n):
        n = int(n)
        if n<1:
            raise ValueError, "Hecke operator T_%s is not defined."%n
        if n==1:
            Mat = matrix_space.MatrixSpace(self.base_ring(),self.rank())
            return Mat(1)

        if arith.is_prime(n):
            return self._compute_hecke_matrix_prime(n)

        F = arith.factor(n)
        if len(F) == 1:  # nontrivial prime power case
            return self._compute_hecke_matrix_prime_power(n, F[0][0],F[0][1])

        else:
            return self._compute_hecke_matrix_general_product(n, F)

    def _compute_hecke_matrix_prime(self, p):
        """
        Compute and return the matrix of the p-th Hecke operator.
        """
        raise NotImplementedError

    def anemic_hecke_algebra(self):
        """
        Return the Hecke algebra associated to this Hecke module.

        EXAMPLES:
            sage: T = ModularSymbols(1,12).hecke_algebra()
            sage: A = ModularSymbols(1,12).anemic_hecke_algebra()
            sage: T == A
            False
            sage: A
            Anemic Hecke algebra acting on Modular Symbols space of dimension 3 for Gamma_0(1) of weight 12 with sign 0 over Rational Field

            sage: A.is_anemic()
            True
        """
        try:
            return self.__anemic_hecke_algebra
        except AttributeError:
            self.__anemic_hecke_algebra = algebra.AnemicHeckeAlgebra(self)
            return self.__anemic_hecke_algebra

    def base_ring(self):
        return self.__base_ring

    def basis_matrix(self):
        return self.free_module().basis_matrix()

    def category(self):
        return sage.categories.all.HeckeModules(self.base_ring())

    def character(self):
        return None

    def dimension(self):
        return self.rank()

    def hecke_algebra(self):
        """
        Return the Hecke algebra associated to this Hecke module.

        EXAMPLES:
            sage: T = ModularSymbols(Gamma1(5),3).hecke_algebra()
            sage: T
            Full Hecke algebra acting on Modular Symbols space of dimension 4 for Gamma_1(5) of weight 3 with sign 0 and over Rational Field
            sage: T.is_anemic()
            False

            sage: M = ModularSymbols(37,sign=1)
            sage: E, A, B = M.decomposition()
            sage: A.hecke_algebra() == B.hecke_algebra()
            False
        """
        try:
            return self.__hecke_algebra
        except AttributeError:
            self.__hecke_algebra = algebra.HeckeAlgebra(self)
            return self.__hecke_algebra

    def is_zero(self):
        """
        Return True if this modular symbols space has dimension 0.
        """
        return self.dimension() == 0

    def is_full_hecke_module(self):
        """
        Return True if this space is invariant under all Hecke
        operators.
        """
        try:
            return self._is_full_hecke_module
        except AttributeError:
            pass

        # now compute whether invariant under Hecke operators of index
        # dividing the level
        misc.verbose("Determining if Hecke module is full.")
        N = self.level()
        for p in arith.prime_divisors(N):
            if not self.is_hecke_invariant(p):
                self._is_full_hecke_module = False
                return False
        self._is_full_hecke_module = True
        return True

    def is_hecke_invariant(self, n):
        """
        Return True if self is invariant under the Hecke operator $T_n$.
        """
        if arith.gcd(n, self.level()) == 1:
            return True
        try:
            self.hecke_operator(n)
        except ArithmeticError:
            return False
        return True

    def level(self):
        """
        Returns the level of this modular symbols space.

        INPUT:
           ModularSymbols self -- an arbitrary space of modular symbols

        OUTPUT:
           int -- the level

        EXAMPLES:
            sage: m = ModularSymbols(20)
            sage: m.level()
            20
        """
        return self.__level

    def submodule(self, X):
        raise NotImplementedError


class HeckeModule_free_module(HeckeModule_generic):
    """
    A Hecke module modeled on a free module over a ring.
    """
    def __init__(self, base_ring, level, weight):
        HeckeModule_generic.__init__(self, base_ring, level)
        self.__weight = weight

    def __cmp__(self, other):
        if not isinstance(other, HeckeModule_free_module):
            return -1
        c = HeckeModule_generic.__cmp__(self, other)
        if c: return c
        return cmp(self.__weight, other.__weight)

    def __contains__(self, x):
        if not element.is_HeckeModuleElement(x):
            return False
        if x.parent() == self:  # easy
            return True
        return x.element() in self.free_module()

    def __getitem__(self, n):
        n = int(n)
        D = self.decomposition()
        if n < 0 or n >= len(D):
            raise IndexError, "index (=%s) must be between 0 and %s"%(n, len(D)-1)
        return D[n]

    def __hash__(self):
        return hash((self.__weight, self.level(), self.base_ring(), str(self)))

    def __len__(self):
        return len(self.decomposition())

    def _eigen_nonzero(self):
        try:
            return self.__eigen_nonzero
        except AttributeError:
            pass
        # Find the smallest integer i such that the i-th entries
        # of the entries of a basis for the dual vector space
        # are not all 0.  Then return the i-th basis vector for
        # the ambient space of modular symbols.
        A = self.ambient_hecke_module()
        V = self.dual_free_module()
        B = V.basis()
        for i in range(A.rank()):
            for b in B:
                if b[i] != 0:
                    self.__eigen_nonzero = i
                    return i
        assert False, 'bug in _eigen_nonzero'

    def _eigen_nonzero_element(self, n=1):
        r"""
        Return $T_n(x)$ where $x$ is a sparse modular symbol such that
        the image of $x$ is nonzero under the dual projection map
        associated to this space, and $T_n$ is the $n$-th Hecke
        operator.
        """
        if self.rank() == 0:
            raise ArithmeticError, "the rank of self must be positive"
        A = self.ambient_hecke_module()
        i = self._eigen_nonzero()
        return A._hecke_image_of_ith_basis_vector(n, i)

    def _hecke_image_of_ith_basis_vector(self, n, i):
        r"""
        Return $T_n(e_i)$, where $e_i$ is the $i$th basis vector of
        the ambient space.
        """
        T = self.hecke_operator(n)
        return T.apply_sparse(self.gen(i))

    def _element_eigenvalue(self, x):
        if not element.is_HeckeModuleElement(x):
            raise TypeError, "x must be a Hecke module element."
        if not x in self.ambient_hecke_module():
            raise ArithmeticError, "x must be in the ambient Hecke module."
        v = self.dual_eigenvector()
        return v.dot_product(x.element())

    def _set_factor_number(self, i):
        self.__factor_number = i

    def ambient(self):
        return self.ambient_hecke_module()

    def ambient_module(self):
        return self.ambient_hecke_module()

    def atkin_lehner_operator(self, d=None):
        """
        Return the Atkin-Lehner operator $W_d$ on this space, if
        defined, where $d$ is a divisor of the level $N$ such that
        $N/d$ and $d$ are coprime.

        EXAMPLES:
            sage: M = ModularSymbols(11)
            sage: w = M.atkin_lehner_operator()
            sage: w
            Hecke module morphism Atkin-Lehner operator W_11 defined by the matrix
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0 -1]
            Domain: Modular Symbols space of dimension 3 for Gamma_0(11) of weight ...
            Codomain: Modular Symbols space of dimension 3 for Gamma_0(11) of weight ...
            sage: M = ModularSymbols(Gamma1(13))
            sage: w = M.atkin_lehner_operator()
            sage: w.fcp()
            (x - 1)^7 * (x + 1)^8

            sage: M = ModularSymbols(33)
            sage: S = M.cuspidal_submodule()
            sage: S.atkin_lehner_operator()
            Hecke module morphism Atkin-Lehner operator W_33 defined by the matrix
            (not printing 6 x 6 matrix)
            Domain: Modular Symbols subspace of dimension 6 of Modular Symbols space ...
            Codomain: Modular Symbols subspace of dimension 6 of Modular Symbols space ...

            sage: S.atkin_lehner_operator(3)
            Hecke module morphism Atkin-Lehner operator W_3 defined by the matrix
            (not printing 6 x 6 matrix)
            Domain: Modular Symbols subspace of dimension 6 of Modular Symbols space ...
            Codomain: Modular Symbols subspace of dimension 6 of Modular Symbols space ...

            sage: N = M.new_submodule()
            sage: N.atkin_lehner_operator()
            Hecke module morphism Atkin-Lehner operator W_33 defined by the matrix
            [  1 2/5 4/5]
            [  0  -1   0]
            [  0   0  -1]
            Domain: Modular Symbols subspace of dimension 3 of Modular Symbols space ...
            Codomain: Modular Symbols subspace of dimension 3 of Modular Symbols space ...
        """
        if d is None:
            d = self.level()
        d = int(d)
        if self.level() % d != 0:
            raise ArithmeticError, "d (=%s) must be a divisor of the level (=%s)"%(d,self.level())

        N = self.level()
        for p, e in arith.factor(d):
            v = arith.valuation(N, p)
            if e < v:
                d *= p**(v-e)
        d = int(d)
        try:
            return self.__atkin_lehner_operator[d]
        except AttributeError:
            self.__atkin_lehner_operator = {}
        except KeyError:
            pass
        Wmat = self._compute_atkin_lehner_matrix(d)
        H = self.endomorphism_ring()
        W = H(Wmat, "Atkin-Lehner operator W_%s"%d)
        self.__atkin_lehner_operator[d] = W
        return W

    def basis(self):
        """
        Returns a basis for self.
        """
        try:
            return self.__basis
        except AttributeError:
            self.__basis = self.gens()
        return self.__basis

    def decomposition(self, bound=None, anemic=True, compute_dual=False):
        """
        Returns the maximal decomposition of this Hecke module under
        the action of Hecke operators of index coprime to the level.
        This is the finest decomposition of self that we can obtain
        using factors obtained by taking kernels of Hecke operators.

        Each factor in the decomposition is a Hecke submodule obtained
        as the kernel of $f(T_n)^r$ acting on self, where n is coprime
        to the level and $r=1$.  If anemic if False, instead choose
        $r$ so that $f(X)^r$ exactly divides the characteristic
        polynomial.

        INPUT:
            anemic -- bool (default: True),  if True, use only Hecke operators
                      of index coprime to the level.

            compute_dual -- bool (default: False) also compute dual subspaces
                       along the way.  These are useful for many algorithms.
                       This is only allowed for ambient Hecke modules.

            bound -- int or None, (default: None).  If None, use all Hecke operators
                     up to the Sturm bound, and hence obtain the same result as
                     one would obtain by using every element of the Hecke ring.
                     If a fixed integer, decompose using only Hecke operators
                     T_p, with p prime, up to bound.

        OUTPUT:
            list -- a list of subspaces of self.

        """
        if compute_dual and not self.is_ambient():
            raise ValueError, "compute_dual can only be True for ambient spaces."
        if not isinstance(anemic, bool):
            raise TypeError, "anemic must be of type bool."

        key = (bound, anemic)

        try:
            if self.__decomposition[key] != None:
                return self.__decomposition[key]
        except AttributeError:
            self.__decomposition = {}
        except KeyError:
            pass
        if self.rank() == 0:
            self.__decomposition[key] = Sequence([], immutable=True, cr=True)
            return self.__decomposition[key]

        time = misc.verbose("Decomposing %s"%self)
        T = self.ambient_hecke_module().hecke_algebra()
        if bound is None:
            bound = self.hecke_bound()
        D = Sequence([], cr=True)
        U = [self.free_module()]
        if compute_dual:
            Udual = [self.free_module()]
        p = 2
        while len(U) > 0 and p <= bound:
            misc.verbose(mesg="p=%s"%p,t=time)
            if anemic:
                while arith.GCD(p, self.level()) != 1:
                    p = arith.next_prime(p)
            misc.verbose("Decomposition using p=%s"%p)
            t = T.hecke_operator(p).matrix()
            Uprime = []
            Udualprime = []
            for i in range(len(U)):
                if self.base_ring().characteristic() == 0 and self.level()%p != 0:
                    is_diagonalizable = True
                else:
                    is_diagonalizable = False
                X = t.decomposition_of_subspace(U[i], is_diagonalizable=is_diagonalizable)
                if compute_dual:
                    Xdual = t.transpose().decomposition_of_subspace(Udual[i], is_diagonalizable=is_diagonalizable)
                    if len(X) != len(Xdual):
                        raise RuntimeError, "Unable to compute compatible dual decomposition."
                for i in range(len(X)):
                    W, is_irred = X[i]
                    if is_irred:
                        A = self.submodule(W)
                        if compute_dual:
                            if X[i][1] != is_irred:
                                raise RuntimeError, "Unable to compute compatible dual decomposition."
                            A._set_dual_free_module(Xdual[i][0])
                        D.append(A)
                    else:
                        Uprime.append(W)
                        if compute_dual:
                            Udualprime.append(Xdual[i][0])
            # end for
            p = arith.next_prime(p)
            U = Uprime
            Udual = Udualprime
        #end while
        for i in range(len(U)):
            if compute_dual:
                A = self.submodule(U[i], Udual[i])
            else:
                A = self.submodule(U[i])
            D.append(A)
        for A in D:
            if anemic:
                A.__is_splittable_anemic = False
                A.__is_splittable = False
            else:
                A.__is_splittable = False
        self.__is_splittable = len(D) > 1
        if anemic:
            self.__is_splittable_anemic = len(D) > 1

        D.sort()
        D.set_immutable()
        self.__decomposition[key] = D
        for i in range(len(D)):
            self.__decomposition[key][i]._set_factor_number(i)
        return self.__decomposition[key]

    def degree(self):
        return self.free_module().degree()

    def dual_eigenvector(self):
        """
        Return an eigenvector for the Hecke operators acting on the
        linear dual of this space.  This eigenvector will have entries
        in an extension of the base ring of degree equal to the
        dimension of this space.

        INPUT:
            The input space must be simple.

        OUTPUT:
            A vector with entries possibly in an extension of the base
            ring.  This vector is an eigenvector for all Hecke operators
            acting via their transpose.

        NOTES:
            (1) The answer is cached so subsequent calls always return
            the same vector.   However, the algorithm is randomized,
            so calls during another session may yield a different
            eigenvector.  This function is used mainly for computing
            systems of Hecke eigenvalues.

            (2) One can also view a dual eigenvector as defining (via
            dot product) a functional phi from the ambient space of
            modular symbols to a field.  This functional phi is an
            eigenvector for the dual action of Hecke operators on
            functionals.
        """
        try:
            return self.__dual_eigenvector
        except AttributeError:
            pass

        if not self.is_simple():
            raise ArithmeticError, "self must be simple"

        # Find a Hecke operator that acts irreducibly on this space:
        p = 2
        t = self.dual_hecke_matrix(p)
        while True:
            f = t.charpoly()
            if f.is_irreducible():
                break
            p = arith.next_prime(p)
            t += random.choice([-2,-1,1,2]) * self.dual_hecke_matrix(p)

        # Write down the eignvector.
        # Write f(x) = (x-alpha)*g(x), where alpha is a root
        # of f(x).
        n = f.degree()
        if n > 1:
            R = f.parent()
            K = R.quotient(f, 'alpha')    # Let K be the quotient R/(f),
                                          # with generator printed "alpha".
            alpha = K.gen()
            beta = ~alpha   # multiplicative inverse of alpha
            c = [-f[0]*beta]
            for i in range(1,n-1):
                c.append((c[i-1] - f[i])*beta)
            c.append( K(1) )
        else:
            K = self.base_ring()
            c = [1]

        # The entries of c are the coefficients of g (as stated in
        # William Stein's Ph.D. thesis, Section 3.5.3).  We compute
        # g(t)v for a some vector v, and get an eigenvector.
        V = FreeModule(K, n)
        t = t.change_ring(K)      # coerce t to be over K.
        for j in range(n):
            v = V.gen(j)
            I = t.iterates(v, n)  # iterates v, v*t, v*t^2, ...
            w = V(0)
            for i in range(n):
                w += c[i]*V(I[i].list())
            if w != 0:
                break

        # Now w is an eigenvector for the action of the Hecke
        # operators on the subspace.  We need an eigenvector
        # in the original space, so we take the linear combination
        # of the basis for the embedded dual vector space given
        # by the entries of w.
        Vdual = self.dual_free_module().change_ring(K)
        w_lift = Vdual.linear_combination_of_basis(w)

        # Finally rescale so the dot product of this vector and
        # the _eigen_nonzero_element is 1.
        alpha = w_lift.dot_product(self._eigen_nonzero_element().element())
        w_lift = w_lift * (~alpha)

        self.__dual_eigenvector = w_lift
        return self.__dual_eigenvector

    def dual_hecke_matrix(self, n):
        """
        The matrix of the $n$-th Hecke operator acting on the dual
        embedded representation of self.
        """
        n = int(n)
        try:
            self._dual_hecke_matrices
        except AttributeError:
            self._dual_hecke_matrices = {}
        if not self._dual_hecke_matrices.has_key(n):
            T = self._compute_dual_hecke_matrix(n)
            self._dual_hecke_matrices[n] = T
        return self._dual_hecke_matrices[n]

    def eigenvalue(self, n):
        """
        Assuming that self is a simple space, return the eigenvalue of
        the $n$th Hecke operator on self.

        NOTES:
        (1) In fact there are $d$ systems of eigenvalues associated to
        self, where $d$ is the rank of self.  Each of the systems of
        eigenvalues is conjugate over the base field.  This function
        chooses one of the systems and consistently returns
        eigenvalues from that system.  Thus these are the coefficients
        $a_n$ for $n\geq 1$ of a modular eigenform attached to self.

        (2) This function works even for Eisenstein subspaces, though
        it will not give the constant coefficient of one of the
        corresponding Eisenstein series (i.e., the generalized
        Bernoulli number).
        """
        if not self.is_simple():
            raise ArithmeticError, "self must be simple"
        n = int(n)
        try:
            return self.__eigenvalues[n]
        except AttributeError:
            self.__eigenvalues = {}
        except KeyError:
            pass
        if n <= 0:
            raise IndexError, "n must be a positive integer"

        if n == 1:
            a1 = self.base_ring()(1)
            self.__eigenvalues[1] = a1
            return a1

        if arith.is_prime(n):
            Tn_e = self._eigen_nonzero_element(n)
            an = self._element_eigenvalue(Tn_e)
            self.__eigenvalues[n] = an
            return an

        # Now use the Hecke eigenvalue recurrence, since arithmetic in
        # a field is faster than computing Heilbronn matrices for
        # non-prime n and doing some big sum (i.e., computing T_n(e)).
        # Also by computing using the recurrence on eigenvalues
        # we use information about divisors.
        F = arith.factor(n)
        prod = None
        for p, r in F:
            (p, r) = (int(p), int(r))
            pow = p**r
            if not self.__eigenvalues.has_key(pow):
                # TODO: Optimization -- do something much more intelligent in case character is not defined.
                # For example, compute it using diamond operators <d>
                eps = self.character()
                if eps is None:
                    Tn_e = self._eigen_nonzero_element(pow)
                    self.__eigenvalues[pow] = self._element_eigenvalue(Tn_e)
                else:
                    # a_{p^r} := a_p * a_{p^{r-1}} - eps(p)p^{k-1} a_{p^{r-2}}
                    apr1 = self.eigenvalue(pow//p)
                    ap = self.eigenvalue(p)
                    k = self.weight()
                    apr2 = self.eigenvalue(pow//(p*p))
                    apow = ap*apr1 - eps(p)*(p**(k-1)) * apr2
                    self.__eigenvalues[pow] = apow
            if prod == None:
                prod = self.__eigenvalues[pow]
            else:
                prod *= self.__eigenvalues[pow]
        self.__eigenvalues[n] = prod
        return prod

    def factor_number(self):
        """
        If this Hecke module was computed via a decomposition of
        another Hecke module, this is the corresponding number.
        Otherwise return -1.
        """
        try:
            return self.__factor_number
        except AttributeError:
            return -1

    def gen(self, n):
        return self(self.free_module().gen(n))

    def hecke_matrix(self, n):
        """
        The matrix of the $n$-th Hecke operator acting on given basis.
        """
        n = int(n)
        if n <= 0:
            raise IndexError, "n must be positive."
        try:
            self._hecke_matrices
        except AttributeError:
            self._hecke_matrices = {}
        if not self._hecke_matrices.has_key(n):
            T = self._compute_hecke_matrix(n)
            T.set_immutable()
            self._hecke_matrices[n] = T
        return self._hecke_matrices[n]

    def hecke_operator(self, n):
        """
        Returns the n-th Hecke operator $T_n$.

        INPUT:
           ModularSymbols self -- Hecke equivariant space of
                                  modular symbols
           int n -- an integer at least 1.
        """
        return self.hecke_algebra().hecke_operator(n)

    def T(self, n):
        r"""
        Returns the $n$-th Hecke operator $T_n$.  This function is a
        synonym for \code{hecke_operator}.
        """
        return self.hecke_operator(n)

    def hecke_polynomial(self, n):
        """
        Return the characteristic polynomial of the n-th Hecke operator
        acting on this space.

        INPUT:
            n -- integer
        OUTPUT:
            a polynomial
        """
        return self.hecke_operator(n).charpoly()

    def is_simple(self):
        raise NotImplementedError

    def is_splittable(self):
        """
        Returns True if and only if only it is possible to split
        off a nontrivial generalized eigenspace of self as the
        kernel of some Hecke operator.
        """
        if not hasattr(self, "__is_splittable"):
            self.decomposition(anemic=False)
        return self.__is_splittable

    def is_submodule(self, other):
        if not isinstance(other, HeckeModule_free_module):
            return False
        return self.ambient_free_module() == other.ambient_free_module() and \
               self.free_module().is_submodule(other.free_module())

    def is_splittable_anemic(self):
        """
        Returns true if and only if only it is possible to split
        off a nontrivial generalized eigenspace of self as the
        kernel of some Hecke operator of index coprime to the level.
        """
        if not hasattr(self,"__is_splittable_anemic"):
            self.decomposition(anemic=True)
        return self.__is_splittable_anemic

    def ngens(self):
        return self.rank()

    def projection(self):
        r"""
        Return the projection map from the ambient space to self.


        ALGORITHM:
            Let $B$ be the matrix whose columns are got by
            concatenating together a basis for the factors of the
            ambient space.  Then the projection matrix onto self is
            the submatrix of $B^{-1}$ got from the rows corresponding
            to self, i.e., if the basis vectors for self appear as
            columns $n$ through $m$ of $B$, then the projection matrix
            is got from rows $n$ through $m$ of $B^{-1}$.  This is
            because projection with respect to the B basis is just
            given by an $m-n+1$ row slice $P$ of a diagonal matrix D
            with 1's in the $n$ through $m$ positions, so projection
            with respect to the standard basis is given by $P\cdot
            B^{-1}$, which is just rows $n$ through $m$ of $B^{-1}$.
        """

        # Compute the Hecke-stable projection map pi from the ambient
        # space of M to M.  Computing the projection map is the same
        # as writing the ambient space as a direct sum of M and its
        # Hecke-stable complement, which is the old subspace plus the
        # other new factors, then *inverting*.  With the projection
        # map in hand, we can compute Hecke operators directly on M
        # fairly quickly without having to compute them on the whole
        # ambient space.  Of course, computing this inverse is way too
        # much work to be useful in general (!).  (I sort of learned
        # this trick from Joe Wetherell, or at least he was aware of
        # it when I mentioned it to him in an airport once.  It's also
        # sort of like a trick Cremona uses in his book for elliptic
        # curves.)  It's not a very good trick though.

        try:
            return self.__projection
        except AttributeError:
            i = self.factor_number()
            if i == -1:
                raise NotImplementedError,\
                      "Computation of projection only implemented "+\
                      "for decomposition factors."
            A = self.ambient_hecke_module()
            B = A.decomposition_matrix_inverse()
            n = sum([A[j].rank() for j in range(i)])
            C = B.matrix_from_columns(range(n,n+self.rank()))
            H = homspace.HeckeModuleHomspace(A, self)
            pi = H(C, "Projection"%self)
            self.__projection = pi
            return self.__projection

    def system_of_eigenvalues(self, n):
        r"""
        Assuming that self is a simple space of modular symbols, return
        the eigenvalues $[a_1, \ldots, a_nmax]$ of the Hecke operators
        on self.  See \code{self.eigenvalue(n)} for more details.

        EXAMPLES:
        We compute eigenvalues for newforms of level 62.
            sage: M = ModularSymbols(62,2,sign=-1)
            sage: S = M.cuspidal_submodule().new_submodule()
            sage: [A.system_of_eigenvalues(3) for A in S.decomposition()]  # random output

            [[1, 1, 0], [1, -1, -alpha - 1]]


        Next we define a function that does the above:
            sage: def b(N,k=2):
            ...    t=cputime()
            ...    S = ModularSymbols(N,k,sign=-1).cuspidal_submodule().new_submodule()
            ...    for A in S.decomposition():
            ...        print N, A.system_of_eigenvalues(5)

            sage: b(63)
            63 [1, 1, 0, -1, 2]
            63 [1, alpha, 0, 1, -2*alpha]
        """
        return [self.eigenvalue(m) for m in range(1,n+1)]

    def weight(self):
        """
        Returns the weight of this modular symbols space.

        INPUT:
           ModularSymbols self -- an arbitrary space of modular symbols

        OUTPUT:
           int -- the weight

        EXAMPLES:
            sage: m = ModularSymbols(20, weight=2)
            sage: m.weight()
            2
        """
        return self.__weight



    def zero_submodule(self):
        """
        Return the zero submodule of self.
        """
        return self.submodule(self.free_module().zero_submodule())


