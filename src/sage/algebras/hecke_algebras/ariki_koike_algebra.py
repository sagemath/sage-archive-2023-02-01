"""
Ariki-Koike Algebras

AUTHORS:

- Travis Scrimshaw (2016-04): initial version
"""

#*****************************************************************************
#  Copyright (C) 2016 Travis Scrimshaw <tscrimsh at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.all import ZZ
from sage.categories.algebras import Algebras
from sage.categories.rings import Rings
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutations
from sage.sets.family import Family

class ArikiKoikeAlgebra(CombinatorialFreeModule):
    r"""
    The Ariki-Koike algebra `H_{r,n}(q)`.

    Has basis `\{L_1^c_i \cdots L_n^c_n T_w \mid w \in S_n, 0 \leq c_i < r\}`.

    INPUT:

    - ``r`` -- the maximum power of `L_i`
    - ``n`` -- the rank `S_n`
    - ``q`` -- (optional) an invertible element in a commutative ring;
      the default is `q \in R[q,q^{-1}]`, where `R` is the ring containing
      the variables ``u``
    - ``u`` -- (optional) the variables `u_1, \dotsc, u_r`; the
      default is the generators of `\ZZ[u_1, \dotsc, u_r]`
    - ``R`` -- (optional) a commutative ring containing ``q`` and ``u``; the
      default is the parent of `q` and `u_1, \dotsc, u_r`

    EXAMPLES:

    REFERENCES:

    .. [AK94] Susumu Ariki and Kazuhiko Koike.
       *A Hecke algebra of* `(\ZZ / r\ZZ) \wr \mathfrak{S}_n`
       *and construction of its irreducible representations*.
       Advances in Mathematics **106**, (1994) pp. 216-243.

    .. [HM16] Jun Hu and Andrew Mathas.
       *Seminormal forms and cyclotomic quiver Hecke algebras of type* `A`.
       Mathematische Annalen **364**, (2016) pp. 1189-1254.
       :arxiv:`1304.0906`.
    """
    @staticmethod
    def __classcall_private__(cls, r, n, q=None, u=None, R=None):
        """
        Standardize input to ensure a unqiue representation.

        TESTS::

            sage: H1 = algebras.ArikiKoike(4, 3)
            sage: S = PolynomialRing(ZZ, 'u', 4)
            sage: R.<q> = LaurentPolynomialRing(S)
            sage: H2 = algebras.ArikiKoike(4, 3, q=q)
            sage: H3 = algebras.ArikiKoike(4, 3, q, S.gens(), R)
            sage: H1 is H2
            True
            sage: H2 is H3
            True
        """
        if u is None:
            if q is not None:
                R = q.parent()
            if R is None:
                R = PolynomialRing(ZZ, 'u', r)
                u = R.gens()
                if q is None:
                    R = LaurentPolynomialRing(R, 'q')
                    q = R.gen()
            else:
                u = PolynomialRing(ZZ, 'u', r).gens()
                if q is None:
                    q = 'q'
        else:
            if not isinstance(u, (list,tuple)):
                u = [u]*r
            if R is None:
                from sage.structure.element import get_coercion_model
                cm = get_coercion_model()
                if q is None:
                    R = cm.common_parent(*[val.parent() for val in u])
                    R = LaurentPolynomialRing(R, 'q')
                    q = R.gen()
                else:
                    R = cm.common_parent(q.parent(), *[val.parent() for val in u])
            elif q is None:
                q = 'q'
            u = [R(val) for val in u]
        if R not in Rings().Commutative():
            raise TypeError("base ring must be a commutative ring")
        q = R(q)
        u = tuple(u)
        return super(ArikiKoikeAlgebra, cls).__classcall__(cls, r, n, q, u, R)

    def __init__(self, r, n, q, u, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: TestSuite(H).run()
        """
        self._r = r
        self._n = n
        self._q = q
        self._u = u
        self._Pn = Permutations(n)
        import itertools
        C = itertools.product(*([range(r)]*n))
        indices = list(itertools.product(C, self._Pn))
        cat = Algebras(R).WithBasis()
        CombinatorialFreeModule.__init__(self, R, indices, prefix='T',
                                         category=cat)
        self._assign_names(self.algebra_generators().keys())

    def _repr_(self):
        """ 
        Return a string representation of ``self``.

        EXAMPLES::

            sage: algebras.ArikiKoike(5, 2)
            Ariki-Koike algebra of rank 5 and order 2
             with q=q and u=(u0, u1, u2, u3, u4)
             over Univariate Laurent Polynomial Ring in q
             over Multivariate Polynomial Ring in u0, u1, u2, u3, u4
             over Integer Ring
        """
        return "Ariki-Koike algebra of rank {} and order {} with q={} and u={} over {}".format(
            self._r, self._n, self._q, self._u, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 2)
            sage: latex(H)
            \mathcal{H}_{5,2}(q)
        """
        return "\\mathcal{H}_{%s,%s}(%s)"%(self._r, self._n, self._q)

    def _repr_term(self, m):
        """
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(4, 3)
            sage: H._repr_term( ((1, 0, 2), Permutation([3,2,1])) )
            'L1*L3^2*T[2,1,2]'
        """
        gen_str = lambda e: '' if e == 1 else '^%s'%e
        lhs = '*'.join('L%s'%(j+1) + gen_str(i) for j,i in enumerate(m[0]) if i > 0)
        redword = m[1].reduced_word()
        if not redword:
            if not lhs:
                return '1'
            return lhs
        rhs = 'T[{}]'.format(','.join(str(i) for i in redword))
        if not lhs:
            return rhs
        return lhs + '*' + rhs

    def _latex_term(self, m):
        r"""
        Return a latex representation for the basis element indexed by ``m``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(4, 3)
            sage: H._latex_term( ((1, 0, 2), Permutation([3,2,1])) )
            'L_{1} L_{3}^{2} T_{2} T_{1} T_{2}'
        """
        gen_str = lambda e: '' if e == 1 else '^{%s}'%e
        lhs = ' '.join('L_{%s}'%(j+1) + gen_str(i) for j,i in enumerate(m[0]) if i > 0)
        redword = m[1].reduced_word()
        if not redword:
            if not lhs:
                return '1'
            return lhs
        return lhs + ' ' + ' '.join("T_{%d}"%i for i in redword)

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: dict(H.algebra_generators())
            {'L1': L1, 'L2': L2, 'L3': L3, 'T1': T[1], 'T2': T[2]}
        """
        one = self._Pn.one()
        zero = [0]*self._n
        d = {}
        for i in range(self._n):
            r = list(zero) # Make a copy
            r[i] = 1
            d['L%s'%(i+1)] = self.monomial( (tuple(r), one) )
        G = self._Pn.group_generators()
        for i in range(1, self._n):
            d['T%s'%i] = self.monomial( (tuple(zero), G[i]) )
        return Family(sorted(d), lambda i: d[i])

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: H.gens()
            (L1, L2, L3, T[1], T[2])
        """
        return tuple(self.algebra_generators())

    @cached_method
    def one_basis(self):
        """
        Return the index of the basis element of `1`.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: H.one_basis()
            ((0, 0, 0), [1, 2, 3])
        """
        one = self._Pn.one()
        zero = [0]*self._n
        return (tuple(zero), one)

    def q(self):
        """
        Return the variable `q`.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: H.q()
            q
        """
        return self._q

    def u(self):
        """
        Return the variables `u`.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 3)
            sage: H.u()
            (u0, u1, u2, u3, u4)
        """
        return self._u

    def T(self, i=None):
        """
        Return the generator(s) `T_i`.

        INPUT:

        - ``i`` -- (default: ``None``) the generator `T_i` or if ``None``,
          then the list of all generators `T_i`

        EXAMPLES::

            sage: H = algebras.ArikiKoike(8, 3)
            sage: H.T(1)
            T[1]
            sage: H.T()
            [T[1], T[2]]
        """
        G = self.algebra_generators()
        if i is None:
            return [G['T%s'%j] for j in range(1, self._n)]
        return G['T%s'%i]

    def L(self, i=None):
        """
        Return the generator(s) `L_i`.

        INPUT:

        - ``i`` -- (default: ``None``) the generator `L_i` or if ``None``,
          then the list of all generators `L_i`

        EXAMPLES::

            sage: H = algebras.ArikiKoike(8, 3)
            sage: H.L(2)
            L2
            sage: H.L()
            [L1, L2, L3]
        """
        G = self.algebra_generators()
        if i is None:
            return [G['L%s'%j] for j in range(1, self._n+1)]
        return G['L%s'%i]

    @cached_method
    def product_on_basis(self, m1, m2):
        """
        Return the product of the basis elements indexed by ``m1`` and ``m2``.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(6, 3)
            sage: m = ((1, 0, 2), Permutations(3)([2,1,3]))
            sage: H.product_on_basis(m, m)
            L1*L2*L3^4

            sage: H = algebras.ArikiKoike(4, 3)
            sage: L1,L2,L3,T1,T2 = H.algebra_generators()
            sage: L1 * T1 * L1^2 * T1
            (q^-1)*L1*L2^2 + (q^-1-1)*L1^2*L2*T[1]
            sage: L1^2 * T1 * L1^2 * T1
            (q^-1)*L1^2*L2^2 + (q^-1-1)*L1^3*L2*T[1]
            sage: L1^3 * T1 * L1^2 * T1
            (u0*u1*u2*u3*q^-1-u0*u1*u2*u3)*L2*T[1]
             + ((-u0*u1*u2-u0*u1*u3-u0*u2*u3-u1*u2*u3)*q^-1+(u0*u1*u2+u0*u1*u3+u0*u2*u3+u1*u2*u3))*L1*L2*T[1]
             + ((u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3)*q^-1+(-u0*u1-u0*u2-u1*u2-u0*u3-u1*u3-u2*u3))*L1^2*L2*T[1]
             + ((-u0-u1-u2-u3)*q^-1+(u0+u1+u2+u3))*L1^3*L2*T[1]
             + (q^-1)*L1^3*L2^2

            sage: H = algebras.ArikiKoike(4, 3)
            sage: L1^2 * T1 * L1^3 * T1
            (u0*u1*u2*u3*q^-1-u0*u1*u2*u3)*L2*T[1]
             + ((-u0*u1*u2-u0*u1*u3-u0*u2*u3-u1*u2*u3)*q^-1+(u0*u1*u2+u0*u1*u3+u0*u2*u3+u1*u2*u3))*L1*L2*T[1]
             + ((u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3)*q^-1+(-u0*u1-u0*u2-u1*u2-u0*u3-u1*u3-u2*u3))*L1^2*L2*T[1]
             + (q^-2)*L1^2*L2^3
             + ((-u0-u1-u2-u3)*q^-1+(u0+u1+u2+u3))*L1^3*L2*T[1]
             + (q^-2-q^-1)*L1^3*L2^2*T[1]

            sage: L1^2 * T1*T2*T1 * L2 * L3 * T2
            (q-2*q^2+q^3)*L1^2*L2*L3 - (1-2*q+2*q^2-q^3)*L1^2*L2*L3*T[2]
             - (q^2-q^3)*L1^3*L3*T[1] + (q-q^2+q^3)*L1^3*L3*T[1,2]
        """
        L1,T1 = m1
        L2,T2 = m2
        # Product is of the form L1*T1*L2*T2
        id_perm = self._Pn.one()
        # T1 is trivial, so commuting we just have L1*L2*T2
        if T1 == id_perm:
            L = [(L1[i] + L2[i]) % self._r for i in range(self._n)]
            # The l variables all commute, so we don't have to worry about
            #   the order in which we multiply
            ret = self.monomial((tuple(L), T2))
            return prod(self._Li_r_power(i+1) for i in range(self._n)
                        if L1[i] + L2[i] >= self._r) * ret
        # Let T1 = T_{i_1} ... T_{i_k}. To compute the product, we do:
        # 1 - commute T_{i_k} L2 = x.
        # 2 - Multiply L1 * (T_{i_1} ... T_{i_{k-1}}) * x * T2
        if sum(L2) == 0:
            if T2 == id_perm:
                return self.monomial((L1, T1))
            wd = T2.reduced_word()
            return (self._product_basis_Ti((L1, T1), wd[0])
                    * self.monomial( (tuple([0]*len(L1)), self._Pn.from_reduced_word(wd[1:])) ))
        wd = T1.reduced_word()
        return (self.monomial((L1, self._Pn.from_reduced_word(wd[:-1])))
                * self._product_Ti_L(wd[-1], L2)
                * self.monomial((tuple([0]*len(L1)), T2)))

    def _product_Ti_L(self, i, L):
        r"""
        Return the product `T_i L`.

        We compute `T_i (L_1^{k_1} \cdots L_n^{k_n})`
        by using Lemma 3.3 and Proposition 3.4 of [AK94]_.

        INPUT:

        - ``i`` -- the index `i`
        - ``l`` -- the tuple `(k_1, \dotsc, k_n)`

        EXAMPLES::

            sage: H = algebras.ArikiKoike(5, 4)
            sage: H._product_Ti_L(2, (0, 3, 2, 4))
            (q^-1-1)*L2^2*L3^3*L4^4 + (q^-1)*L2^2*L3^3*L4^4*T[2]
            sage: H._product_Ti_L(3, (0, 3, 2, 4))
            -(1-q)*L2^3*L3^2*L4^4 - (q-q^2)*L2^3*L3^3*L4^3
             + q^2*L2^3*L3^4*L4^2*T[3]
            sage: H._product_Ti_L(1, (0, 3, 2, 4))
            -(1-q)*L2^3*L3^2*L4^4 - (q-q^2)*L1*L2^2*L3^2*L4^4
             - (q^2-q^3)*L1^2*L2*L3^2*L4^4 + q^3*L1^3*L3^2*L4^4*T[1]
        """
        ki = L[i-1]
        ki1 = L[i]
        Lp = list(L) # Make a mutable copy
        d = abs(ki - ki1)
        Lp[i-1] = Lp[i] = min(ki, ki1)
        if ki > ki1:
            ret = self._product_Ti_Lj_k(i, i, ki - ki1)
        elif ki1 > ki:
            ret = self._product_Ti_Lj_k(i, i+1, ki1 - ki)
        else:
            ret = self.T(i)
        return self.monomial((tuple(Lp), self._Pn.one())) * ret

    def _product_Ti_Lj_k(self, i, j, k):
        r"""
        Return the product `T_i L_j^k`.

        We use Lemma 3.3 of [AK94]_ to commute `T_i` and `L_j^k`:

        .. MATH::

            \begin{aligned}
            T_{i-1} L_i^k & = q^k L_{i-1}^k T_{i-1} + (q - 1)
            \left( \sum_{m=1}^k q^{m-1} L_{i-1}^{m-1} L_i^{k-m+1} \right),
            \\ T_i L_i^k & = q^{-k} L_{i+1}^k T_i + (q^{-1} - 1)
            \left( \sum_{m=1}^k q^{m-1} L_i^{k-m} L_{i+1}^m \right).
            \end{aligned}

        INPUT:

        - ``i`` -- the index `i`
        - ``j`` -- the index `j`
        - ``k`` -- the exponent `k`

        EXAMPLES::

            sage: H = algebras.ArikiKoike(4, 3)
            sage: H._product_Ti_Lj_k(2, 2, 1)
            (q^-1-1)*L3 + (q^-1)*L3*T[2]
            sage: H._product_Ti_Lj_k(2, 1, 2)
            L1^2*T[2]
            sage: H._product_Ti_Lj_k(2, 2, 2)
            (q^-2-q^-1)*L3^2 + (q^-2)*L3^2*T[2] + (q^-1-1)*L2*L3
            sage: H._product_Ti_Lj_k(2, 3, 2)
            -(1-q)*L3^2 - (q-q^2)*L2*L3 + q^2*L2^2*T[2]
        """
        G = self._Pn.group_generators()
        if i != j - 1 and i != j:
            L = [0]*self._n
            L[j-1] = k
            return self.monomial((tuple(L), G[i]))
        q = self._q
        def L(i, j, pow_i, pow_j):
            ret = [0]*self._n
            ret[i-1] = pow_i # -1 for indexing
            ret[j-1] = pow_j # -1 for indexing
            return tuple(ret)
        one = self.base_ring().one()
        id_perm = self._Pn.one()
        if i == j - 1:
            c = q - one
            d = {(L(i, j, m-1, k-m+1), id_perm): q**(m-1) * c for m in range(1, k+1)}
            d[L(i, i, k, k), G[i]] = q**k
            return self._from_dict(d, remove_zeros=False)
        # Else j == i
        c = ~q - one
        d = {(L(i, i+1, k-m, m), id_perm): q**(1-m) * c for m in range(1, k+1)}
        d[L(i+1, i+1, k, k), G[i]] = q**-k
        return self._from_dict(d, remove_zeros=False)

    def _product_basis_Ti(self, m, i):
        r"""
        Return the product `l T_w T_i`.

        If the quadratic relation is `T_i^2 = (q - 1) T_i + q`,
        then we have

        .. MATH::

            T_w T_i = \begin{cases}
            T_{ws_i} & \text{if } \ell(ws_i) = \ell(w) + 1, \\
            q T_{ws_i} + (q - 1) T_w & \text{if }
            \ell(w s_i) = \ell(w) - 1.
            \end{cases}

        INPUT:

        - ``m`` -- a pair ``[t, w]``, where ``t`` encodes the monomial
          and ``w``  is an element of the permutation group
        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: H = algebras.ArikiKoike(4, 3)
            sage: m = ((1, 0, 2), Permutations(3)([2,1,3]))
            sage: H._product_basis_Ti(m, 1)
            q*L1*L3^2 - (1-q)*L1*L3^2*T[1]
            sage: L1,L2,L3,T1,T2 = H.algebra_generators()
            sage: (L1 * L3^2) * (T1 * T1)
            q*L1*L3^2 - (1-q)*L1*L3^2*T[1]
        """
        L, w = m
        # We have to flip the side due to Sage's multiplication
        #   convention for permutations
        wi = w.apply_simple_reflection(i, side="left")
        if not w.has_descent(i, side="left"):
            return self.monomial((L, wi))

        d = {(L, wi): self._q, (L, w): self._q - self.base_ring().one()}
        return self._from_dict(d, remove_zeros=False)

    @cached_method
    def _Li_r_power(self, i):
        """
        Return `L_i^r`.

        EXAMPLES::

            sage: H = algebras.ArikiKoike(4, 3)
            sage: for i in range(1,4): H._Li_r_power(i)
            u0*u1*u2*u3 + ((-u0*u1*u2-u0*u1*u3-u0*u2*u3-u1*u2*u3))*L1
             + ((u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3))*L1^2 + ((-u0-u1-u2-u3))*L1^3
            u0*u1*u2*u3*q^4 + (-u0*u1*u2*u3*q^3+u0*u1*u2*u3*q^4)*T[1]
             + ((-u0*u1*u2-u0*u1*u3-u0*u2*u3-u1*u2*u3)*q^3)*L2
             + ((u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3)*q^2)*L2^2
             + ((-u0-u1-u2-u3)*q)*L2^3
             + ((u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3)*q^2+(-u0*u1-u0*u2-u1*u2-u0*u3-u1*u3-u2*u3)*q^3)*L1*L2*T[1]
             + ((-u0-u1-u2-u3)*q+(u0+u1+u2+u3)*q^2)*L1*L2^2*T[1] - (1-q)*L1*L2^3*T[1]
             + ((-u0-u1-u2-u3)*q^2+(u0+u1+u2+u3)*q^3)*L1^2*L2*T[1]
             - (q-q^2)*L1^2*L2^2*T[1] - (q^2-q^3)*L1^3*L2*T[1]
            u0*u1*u2*u3*q^8 + (-u0*u1*u2*u3*q^7+u0*u1*u2*u3*q^8)*T[2]
             + (-u0*u1*u2*u3*q^6+u0*u1*u2*u3*q^7)*T[2,1,2]
             + ((-u0*u1*u2-u0*u1*u3-u0*u2*u3-u1*u2*u3)*q^6)*L3
             + ((u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3)*q^4)*L3^2
             + ((-u0-u1-u2-u3)*q^2)*L3^3
             + ((u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3)*q^4+(-u0*u1-u0*u2-u1*u2-u0*u3-u1*u3-u2*u3)*q^5)*L2*L3*T[2]
             + ((-u0-u1-u2-u3)*q^2+(u0+u1+u2+u3)*q^3)*L2*L3^2*T[2] - (1-q)*L2*L3^3*T[2]
             + ((-u0-u1-u2-u3)*q^3+(u0+u1+u2+u3)*q^4)*L2^2*L3*T[2]
             - (q-q^2)*L2^2*L3^2*T[2] - (q^2-q^3)*L2^3*L3*T[2]
             + ((u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3)*q^4+(-2*u0*u1-2*u0*u2-2*u1*u2-2*u0*u3-2*u1*u3-2*u2*u3)*q^5+(u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3)*q^6)*L1*L3*T[1,2]
             + ((u0*u1+u0*u2+u1*u2+u0*u3+u1*u3+u2*u3)*q^4+(-u0*u1-u0*u2-u1*u2-u0*u3-u1*u3-u2*u3)*q^5)*L1*L3*T[2,1,2]
             + ((-u0-u1-u2-u3)*q^2+(2*u0+2*u1+2*u2+2*u3)*q^3+(-u0-u1-u2-u3)*q^4)*L1*L3^2*T[1,2]
             + ((-u0-u1-u2-u3)*q^2+(u0+u1+u2+u3)*q^3)*L1*L3^2*T[2,1,2]
             - (1-2*q+q^2)*L1*L3^3*T[1,2] - (1-q)*L1*L3^3*T[2,1,2]
             + ((-u0-u1-u2-u3)*q^3+(2*u0+2*u1+2*u2+2*u3)*q^4+(-u0-u1-u2-u3)*q^5)*L1*L2*L3*T[1,2]
             - (q-2*q^2+q^3)*L1*L2*L3^2*T[1,2] - (q^2-2*q^3+q^4)*L1*L2^2*L3*T[1,2]
             + ((-u0-u1-u2-u3)*q^4+(2*u0+2*u1+2*u2+2*u3)*q^5+(-u0-u1-u2-u3)*q^6)*L1^2*L3*T[1,2]
             + ((-u0-u1-u2-u3)*q^4+(u0+u1+u2+u3)*q^5)*L1^2*L3*T[2,1,2]
             - (q^2-2*q^3+q^4)*L1^2*L3^2*T[1,2] - (q^2-q^3)*L1^2*L3^2*T[2,1,2]
             - (q^3-2*q^4+q^5)*L1^2*L2*L3*T[1,2] - (q^4-2*q^5+q^6)*L1^3*L3*T[1,2]
             - (q^4-q^5)*L1^3*L3*T[2,1,2]
        """
        one = self._Pn.one()
        L = lambda exp: tuple([exp if j == i else 0 for j in range(1, self._n+1)])
        if i == 1:
            z = PolynomialRing(self.base_ring(), 'z').gen()
            p = list(prod(z - val for val in self._u))[:-1]
            zero = self.base_ring().zero()
            return self._from_dict({(L(exp), one): coeff for exp,coeff in enumerate(p)
                                    if coeff != zero},
                                   remove_zeros=False)
        # We need to multiply things in the correct order in order to not
        #   trigger an infinite recursion
        m = (L(self._r-1), one)
        return self.T(i-1) * (self.L(i-1) * (self.T(i-1) * self.monomial(m)))

    @cached_method
    def inverse_T(self, i):
        r"""
        Return the inverse of the generator `T_i`.

        From the quadratic relation, we have

        .. MATH::

            T_i^{-1} = q^{-1} T_i + (q^{-1} - 1).

        EXAMPLES::

            sage: H = algebras.ArikiKoike(3, 4)
            sage: [H.inverse_T(i) for i in range(1, 4)]
            [(q^-1-1) + (q^-1)*T[1],
             (q^-1-1) + (q^-1)*T[2],
             (q^-1-1) + (q^-1)*T[3]]

        TESTS::

            sage: H = algebras.ArikiKoike(4, 4)
            sage: all(H.inverse_T(i) * H.T(i) == H.one() for i in range(1, 4))
            True
            sage: all(H.T(i) * H.inverse_T(i) == H.one() for i in range(1, 4))
            True
        """
        c = ~self._q - self.base_ring().one()
        m = self.T(i).leading_support()
        return self._from_dict({m: ~self._q, self.one_basis(): c})

    class Element(CombinatorialFreeModule.Element):
        def inverse(self):
            r"""
            Return the inverse if ``self`` is a basis element.

            EXAMPLES::

                sage: H = algebras.ArikiKoike(3, 4)
                sage: t = prod(H.T()); t
                T[1,2,3]
                sage: t.inverse()
                (q^-3-3*q^-2+3*q^-1-1) + (q^-3-2*q^-2+q^-1)*T[3]
                 + (q^-3-2*q^-2+q^-1)*T[2] + (q^-3-q^-2)*T[3,2]
                 + (q^-3-2*q^-2+q^-1)*T[1] + (q^-3-q^-2)*T[1,3]
                 + (q^-3-q^-2)*T[2,1] + (q^-3)*T[3,2,1]
            """
            if len(self) != 1:
                raise NotImplementedError("inverse only implemented for monomials")
            l,w = self.support_of_term()
            if sum(l) != 0:
                raise NotImplementedError("inverse only implemented for monomials in T variables")
            H = self.parent()
            return ~self[l,w] * H.prod(H.inverse_T(i) for i in reversed(w.reduced_word()))

        __invert__ = inverse

