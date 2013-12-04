r"""
Yangians

AUTHORS:

- Travis Scrimshaw (2013-10-08): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod

#from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from sage.rings.all import ZZ
from sage.rings.infinity import infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.family import Family
from sage.sets.positive_integers import PositiveIntegers
from sage.monoids.indexed_monoid import IndexedFreeAbelianMonoid
from sage.combinat.cartesian_product import CartesianProduct
from sage.combinat.free_module import CombinatorialFreeModule

class Yangian(CombinatorialFreeModule):
    r"""
    The Yangian `Y(\mathfrak{gl}_n)`.

    .. NOTE::

        We using the grading defined by `\deg t_{ij}^{(r)} = r - 1` as opposed
        to the natural filtering of `\deg t_{ij}^{(r)} = r`.

    INPUT:

    - ``base_ring`` -- the base ring
    - ``n`` -- the size `n`
    - ``level`` -- (optional) the level of the Yangian
    - ``variable_name`` -- (default: ``'t'``) the name of the variable

    EXAMPLES::

        sage: Y = Yangian(QQ, 4)
    """
    @staticmethod
    def __classcall_private__(cls, base_ring, n, level=None, variable_name='t'):
        """
        Return the correct parent based upon input.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y2 = Yangian(QQ, 4)
            sage: Y is Y2
            True
            sage: YL = Yangian(QQ, 4, 3)
            sage: YL2 = Yangian(QQ, 4, 3)
            sage: YL is YL2
            True
        """
        if level is not None:
            return YangianLevel(base_ring, n, level, variable_name)
        # We need to specify the parameter name for pickling, so it doesn't pass
        #   ``variable_name`` as ``level``
        return super(Yangian, cls).__classcall__(cls, base_ring, n, variable_name=variable_name)

    def __init__(self, base_ring, n, variable_name):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: TestSuite(Y).run(skip="_test_antipode")
        """
        self._n = n
        category = HopfAlgebrasWithBasis(base_ring).Graded()
        self._index_set = tuple(range(1,n+1))
        # The keys for the basis are tuples of these indices
        # TODO: A parent for all sequences of a given base set for the bases
        indices = CartesianProduct(PositiveIntegers(), self._index_set, self._index_set)
        # We note that the generators are non-commutative, but we always sort
        #   them, so they are, in effect, indexed by the free abelian monoid
        basis_keys = IndexedFreeAbelianMonoid(indices, bracket=False, prefix=variable_name)
        CombinatorialFreeModule.__init__(self, base_ring, basis_keys,
                                         prefix=variable_name, category=category)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Yangian(QQ, 4)
            The Yangian of gl(4) over Rational Field
        """
        return "The Yangian of gl({}) over {}".format(self._n, self.base_ring())

    def _repr_term(self, m):
        """
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: I = Y._indices
            sage: Y._repr_term(I.gen((3,1,2)) * I.gen((4,3,1)))
            't(3)[1,2]*t(4)[3,1]'
        """
        if len(m) == 0:
            return '1'
        prefix = self.prefix()
        return '*'.join('*'.join(prefix + '({})[{},{}]'.format(r,i,j)
                                 for x in range(exp))
                        for (r,i,j), exp in m._sorted_items())

    def _latex_term(self, m):
        r"""
        Return a `\LaTeX` representation of the basis element indexed
        by ``m``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: I = Y._indices
            sage: Y._latex_term(I.gen((3,1,2)) * I.gen((4,3,1)))
            't^{(3)}_{1,2} t^{(4)}_{3,1}'
        """
        if len(m) == 0:
            return '1'
        prefix = self.prefix()
        return ' '.join(' '.join(prefix + '^{{({})}}_{{{},{}}}'.format(r,i,j) for x in range(exp))
                        for (r,i,j), exp in m._sorted_items())

    def gen(self, r, i, j):
        """
        Return the generator `t^{(r)}_{ij}` of ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.gen(2, 1, 3)
            t(2)[1,3]
            sage: Y.gen(12, 2, 1)
            t(12)[2,1]
            sage: Y.gen(0, 1, 1)
            1
            sage: Y.gen(0, 1, 3)
            0
        """
        if r == 0:
            if i == j:
                return self.one()
            return self.zero()
        m = self._indices.gen((r,i,j))
        return self._from_dict({m: self.base_ring().one()}, remove_zeros=False)

    def ngens(self):
        """
        Return the number of generators of ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.ngens()
            +Infinity
        """
        return infinity

    @cached_method
    def one_basis(self):
        """
        Return the basis index of the element `1`.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.one_basis()
            1
        """
        return self._indices.one()

    def is_commutative(self):
        """
        Check if ``self`` is a commutative algebra.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.is_commutative()
            True
        """
        return True

    def degree_on_basis(self, m):
        """
        Return the degree of the monomial index by ``m``.

        The degree of `t_{ij}^{(r)}` is equal to `r - 1`.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.gen(2,1,1).degree()
            1
            sage: x = Y.gen(5,2,3)^4
            sage: x.degree()
            16
            sage: elt = Y.gen(10,3,1) * Y.gen(2,1,1) * Y.gen(1,2,4); elt
            -t(1)[3,1]*t(10)[1,1]*t(10)[1,1] + t(1)[2,4]*t(11)[3,1]
             + t(1)[1,1]*t(1)[2,4]*t(10)[3,1] + t(1)[2,4]*t(10)[3,1]
             + t(2)[1,1]*t(10)[3,1]*t(10)[3,1]
            sage: for s in elt.support(): s, Y.degree_on_basis(s)
            (t[1, 1, 1]*t[1, 2, 4]*t[10, 3, 1], 9)
            (t[1, 2, 4]*t[10, 3, 1], 9)
            (t[1, 2, 4]*t[11, 3, 1], 10)
            (t[1, 3, 1]*t[10, 1, 1]^2, 18)
            (t[2, 1, 1]*t[10, 3, 1]^2, 19)
        """
        return sum(max(0, r[0][0] - 1) * r[1] for r in m._sorted_items())

    def dimension(self):
        """
        Return the dimension of ``self``, which is `\infty`.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.dimension()
            +Infinity
        """
        return infinity

    @cached_method
    def product_on_basis(self, x, y):
        """
        Return the product of two monomials given by ``x`` and ``y``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.gen(12, 2, 1) * Y.gen(2, 1, 1) # indirect doctest
            -t(1)[2,1]*t(12)[1,1] + t(13)[2,1] + t(2)[1,1]*t(12)[2,1]
             + t(1)[1,1]*t(12)[2,1] + t(12)[2,1]
        """
        # If x or y indexed by (), it is 1, so return the other
        if len(x) == 0:
            return self.monomial(y)
        if len(y) == 0:
            return self.monomial(x)
        # If it's smaller, just add it to the front
        if x.trailing_support() <= y.leading_support():
            return self.monomial(x * y)

        # The computation is done on generators, so apply generators one at
        #   a time until all have been applied
        if len(x) != 1:
            cur = self.monomial(y)
            for gen,exp in reversed(list(x)):
                for i in range(exp):
                    cur = self.monomial(gen) * cur
            return cur

        # If we are both generators, then apply the basic computation
        if len(y) == 1:
            return self.product_on_gens(tuple(x.support()[0]), tuple(y.support()[0]))

        # Otherwise we need to commute it along
        rhs = y._sorted_items()
        if rhs[0][1] == 1:
            rhs.pop(0)
        else:
            rhs[0] = (rhs[0][0], rhs[0][1]-1)
        rem_y = self._indices.element_class(self._indices, dict(rhs))
        return self.product_on_gens(tuple(x.support()[0]), tuple(rhs[0][0])) * self.monomial(rem_y)

    @cached_method
    def product_on_gens(self, a, b):
        r"""
        Return the product on two generators indexed by ``a`` and ``b``.

        We assume `(r, i, j) \geq (s, k, l)`, and we start with the basic
        relation:

        .. MATH::

            [t_{ij}^{(r)}, t_{kl}^{(s)}] - [t_{ij}^{(r-1)}, t_{kl}^{(s+1)}]
            = t_{kj}^{(r-1)} t_{il}^{(s)} - t_{kj}^{(s)} t_{il}^{(r-1)}.

        Solving for the first term and using induction we get:

        .. MATH::

            [t_{ij}^{(r)}, t_{kl}^{(s)}] = \sum_{a=1}^s \left(
            t_{kj}^{(a-1)} t_{il}^{(r+s-a)} - t_{kj}^{(r+s-a)}
            t_{il}^{(a-1)} \right).

        Next applying induction on this we get

        .. MATH::

            t_{ij}^{(r)} t_{kl}^{(s)} = t_{kl}^{(s)} t_{ij}^{(r)} +
            \sum t_{ab}^{(m)} t_{cd}^{(l)}

        where `m + l < r + s` and `t_{ab}^{(m)} < t_{cd}^{(l)}`.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.product_on_gens((2,1,1), (12,2,1))
            t(2)[1,1]*t(12)[2,1]
            sage: Y.gen(2, 1, 1) * Y.gen(12, 2, 1)
            t(2)[1,1]*t(12)[2,1]
            sage: Y.product_on_gens((12,2,1), (2,1,1))
            -t(1)[2,1]*t(12)[1,1] + t(13)[2,1] + t(2)[1,1]*t(12)[2,1]
             + t(1)[1,1]*t(12)[2,1] + t(12)[2,1]
            sage: Y.gen(12, 2, 1) * Y.gen(2, 1, 1)
            -t(1)[2,1]*t(12)[1,1] + t(13)[2,1] + t(2)[1,1]*t(12)[2,1]
             + t(1)[1,1]*t(12)[2,1] + t(12)[2,1]
        """
        I = self._indices
        if a <= b:
            return self.monomial(I.gen(a) * I.gen(b))
        mid = self.zero() # This is the special term for x = 1
        if b[1] == a[2]:
            mid += self.monomial( I.gen([a[0]+b[0]-1, a[1], b[2]]) )
        if a[1] == b[2]:
            mid -= self.monomial( I.gen([a[0]+b[0]-1, b[1], a[2]]) )
        return self.monomial(I.gen(b) * I.gen(a)) + mid + self.sum(
                self.monomial( I.gen([x-1, b[1], a[2]]) * I.gen([a[0]+b[0]-x, a[1], b[2]]) )
                - self.product_on_gens((a[0]+b[0]-x, b[1], a[2]), (x-1, a[1], b[2]))
                for x in range(2, b[0]+1))

    def coproduct_on_basis(self, m):
        """
        Return the coproduct on the basis element indexed by ``m``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.gen(2,1,1).coproduct() # indirect doctest
            1 # t(2)[1,1] + t(1)[1,1] # t(1)[1,1] + t(1)[1,2] # t(1)[2,1]
             + t(1)[1,3] # t(1)[3,1] + t(1)[1,4] # t(1)[4,1] + t(2)[1,1] # 1
            sage: Y.gen(2,3,1).coproduct()
            1 # t(2)[3,1] + t(1)[3,1] # t(1)[1,1] + t(1)[3,2] # t(1)[2,1]
             + t(1)[3,3] # t(1)[3,1] + t(1)[3,4] # t(1)[4,1] + t(2)[3,1] # 1
            sage: Y.gen(2,2,3).coproduct()
            1 # t(2)[2,3] + t(1)[2,1] # t(1)[1,3] + t(1)[2,2] # t(1)[2,3]
             + t(1)[2,3] # t(1)[3,3] + t(1)[2,4] # t(1)[4,3] + t(2)[2,3] # 1
        """
        T = self.tensor_square()
        I = self._indices
        return T.prod(T.monomial( (I.one(), I.gen([a[0],a[1],a[2]])) )
                      + T.monomial( (I.gen([a[0],a[1],a[2]]), I.one()) )
                      + T.sum_of_terms([( (I.gen([s,a[1],k]), I.gen([a[0]-s,k,a[2]])), 1 )
                                        for k in range(1, self._n+1)
                                        for s in range(1, a[0])])
                      for a,exp in m._sorted_items() for p in range(exp))

    def counit_on_basis(self, m):
        """
        Return the antipode on the basis element indexed by ``m``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4)
            sage: Y.gen(2,3,1).counit() # indirect doctest
            0
            sage: Y.gen(0,0,0).counit()
            1
        """
        if len(m) == 0:
            return self.base_ring().one()
        return self.base_ring().zero()

class YangianLevel(Yangian):
    r"""
    The Yangian `Y_{\ell}(\mathfrak{gl_n})` of level `\ell`.
    """
    def __init__(self, base_ring, n, level, variable_name):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4, 3)
            sage: TestSuite(Y).run(skip="_test_antipode")
        """
        self._level = level
        Yangian.__init__(self, base_ring, n, variable_name)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Yangian(QQ, 4, 3)
            The Yangian of level 3 of gl(4) over Rational Field
        """
        return "The Yangian of level {} of gl({}) over {}".format(self._level, self._n, self.base_ring())

    def level(self):
        """
        Return the level of ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 3, 5)
            sage: Y.level()
            5
        """
        return self._level

    def defining_polynomial(self, i, j, u=None):
        """
        Return the defining polynomial of ``i`` and ``j``.

        The defining polynomial is given by:

        .. MATH::

            T_{ij}(u) = \delta_{ij} u^{\ell} + \sum_{k=1}^{\ell} t_{ij}^{(k)}
            u^{\ell-k}.

        EXAMPLES::

            sage: Y = Yangian(QQ, 3, 5)
            sage: Y.defining_polynomial(3, 2)
            t(1)[3,2]*u^4 + t(2)[3,2]*u^3 + t(3)[3,2]*u^2 + t(4)[3,2]*u + t(5)[3,2]
            sage: Y.defining_polynomial(1, 1)
            u^5 + t(1)[1,1]*u^4 + t(2)[1,1]*u^3 + t(3)[1,1]*u^2 + t(4)[1,1]*u + t(5)[1,1]
        """
        if u is None:
            u = PolynomialRing(self.base_ring(), 'u').gen(0)
        ell = self._level
        return sum(self.gen(k, i, j) * u**(ell-k) for k in range(ell+1))

    def quantum_determinant(self, u=None):
        """
        Return the quantum determinant of ``self``.

        The quantum determinant is defined by:

        .. MATH::

            qdet(u) = \sum_{\sigma \in S_n} (-1)^{\sigma}
            \prod_{k=1}^n T_{\sigma(k),k}(u - k + 1).

        EXAMPLES::

            sage: Y = Yangian(QQ, 2, 2)
            sage: Y.quantum_determinant()
            u^4 + (-2 + t(1)[1,1] + t(1)[2,2])*u^3
             + (t(2)[1,1] + 1 - t(1)[1,1] + t(2)[2,2] + t(1)[1,1]*t(1)[2,2]
                - t(1)[1,2]*t(1)[2,1] - 2*t(1)[2,2])*u^2
             + (-t(1)[1,2]*t(2)[2,1] + t(1)[1,1]*t(2)[2,2]
                + t(1)[1,2]*t(1)[2,1] - t(2)[1,1] + t(1)[2,2]*t(2)[1,1] - t(2)[2,2]
                - t(1)[1,1]*t(1)[2,2] - t(1)[2,1]*t(2)[1,2] + t(1)[2,2])*u
             + t(1)[1,2]*t(2)[2,1] - t(2)[1,2]*t(2)[2,1] - t(3)[2,2]
                - t(1)[1,1]*t(2)[2,2] + t(3)[1,1] + t(2)[1,1]*t(2)[2,2] + t(2)[2,2]
        """
        if u is None:
            u = PolynomialRing(self.base_ring(), 'u').gen(0)
        from sage.combinat.permutation import Permutations
        n = self._n
        return sum(p.sign() * prod(self.defining_polynomial(p[k], k+1, u - k)
                                   for k in range(n))
                   for p in Permutations(n))

    def gen(self, r, i, j):
        """
        Return the generator `t^{(r)}_{ij}` of ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4, 3)
            sage: Y.gen(2, 1, 3)
            t(2)[1,3]
            sage: Y.gen(12, 2, 1)
            0
            sage: Y.gen(0, 1, 1)
            1
            sage: Y.gen(0, 1, 3)
            0
        """
        if r > self._level:
            return self.zero()
        return Yangian.gen(self, r, i, j)

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 2, 2)
            sage: Y.gens()
            (t(1)[1,1], t(2)[1,1], t(1)[1,2], t(2)[1,2], t(1)[2,1],
             t(2)[2,1], t(1)[2,2], t(2)[2,2])
        """
        return tuple(self.gen(r, i, j)
                     for i in range(1, self._n+1)
                     for j in range(1, self._n+1)
                     for r in range(1, self._level+1))

    def ngens(self):
        """
        Return the number of generators of ``self``.

        EXAMPLES::

            sage: Y = Yangian(QQ, 4, 3)
            sage: Y.ngens()
            48
        """
        return self._level * self._n**2

    @cached_method
    def product_on_basis(self, x, y):
        """
        Return the product of two monomials given by ``x`` and ``y``.

        .. SEEALSO::

            :meth:`Yangian.product_on_gens()`

        EXAMPLES::

            sage: Y = Yangian(QQ, 4, 3)
            sage: elt = Y.gen(3,2,1) * Y.gen(1,1,3)
            sage: elt * Y.gen(1, 1, 2) # indirect doctest
            t(1)[1,3]*t(3)[2,1]*t(3)[2,1] - t(3)[1,3] + t(1)[1,3]*t(3)[2,2]
             + t(1)[1,2]*t(3)[2,3] - t(1)[1,3]*t(3)[1,1]
        """
        ret = super(YangianLevel, self).product_on_basis(x, y)
        return self._from_dict({m:c for m,c in ret
                               if len(m) == 0 or m.trailing_support()[0] <= self._level})

    @cached_method
    def product_on_gens(self, a, b):
        r"""
        Return the product on two generators indexed by ``a`` and ``b``.

        .. SEEALSO::

            :meth:`Yangian.product_on_gens()`

        EXAMPLES::

            sage: Y = Yangian(QQ, 4, 3)
            sage: Y.gen(1,2,2) * Y.gen(2,1,3) # indirect doctest
            t(1)[2,2]*t(2)[1,3]
            sage: Y.gen(1,2,1) * Y.gen(2,1,3) # indirect doctest
            t(1)[2,1]*t(2)[1,3]
            sage: Y.gen(3,2,1) * Y.gen(1,1,3) # indirect doctest
            t(3)[2,3] + t(1)[1,3]*t(3)[2,1]
        """
        ret = super(YangianLevel, self).product_on_gens(a, b)
        return self._from_dict({m:c for m,c in ret
                               if len(m) == 0 or m.trailing_support()[0] <= self._level})

