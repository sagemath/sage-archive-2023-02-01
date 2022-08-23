r"""
`q`-Commuting Polynomials

AUTHORS:

- Travis Scrimshaw (2022-08-23): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 Travis Scrimshaw  <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.rings.infinity import infinity
from sage.categories.algebras import Algebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.monoids.free_abelian_monoid import FreeAbelianMonoid, FreeAbelianMonoid_class

class qCommutingPolynomials(CombinatorialFreeModule):
    r"""
    The algebra of `q`-commuting polynomials.

    Let `R` be a commutative ring, and fix an element `q \in R`. We say two
    distinct variables `x` and `y` `q`-*commute* if they satisfy the relation

    .. MATH::

        x y = q \cdot y x.

    These form a graded `R`-algebra with a natural basis given by monomials
    written in increasing order. These then satisfy a `q`-analog of the
    classical binomial coefficient theorem:

    .. MATH::

        (x + y)^n = \sum_{k=0}^n \binom{n}{k}_q x^k y^{n-k}.

    EXAMPLES::

        sage: q = ZZ['q'].fraction_field().gen()
        sage: R.<x,y> = algebras.qCommutingPolynomials(q)

    We verify the `q`-binomial theorem::

        sage: f = (x + y)^10
        sage: all(f[b] == q_binomial(10, b.list()[0]) for b in f.support())
        True
    """
    @staticmethod
    def __classcall_private__(cls, q, n=None, base_ring=None, names=None):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R1.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R2 = algebras.qCommutingPolynomials(q, base_ring=q.parent(), names='x,y,z')
            sage: R3 = algebras.qCommutingPolynomials(q, names=['x', 'y', 'z'])
            sage: R1 is R2 is R3
            True
        """
        if base_ring is not None:
            q = base_ring(q)

        if isinstance(n, FreeAbelianMonoid_class):
            indices = n
        else:
            indices = FreeAbelianMonoid(n, names)
        return super().__classcall__(cls, q, indices)

    def __init__(self, q, indices):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: TestSuite(R).run()
        """
        self._q = q
        base_ring = q.parent()
        category = Algebras(base_ring).WithBasis().Graded()
        CombinatorialFreeModule.__init__(self, base_ring, indices,
                                         bracket=False, prefix='',
                                         sorting_key=qCommutingPolynomials._term_key,
                                         names=indices.variable_names(), category=category)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R
            q-commuting polynomial ring in x, y, z over Fraction Field of
             Univariate Polynomial Ring in q over Integer Ring
        """
        names = ", ".join(self.variable_names())
        return "{}-commuting polynomial ring in {} over {}".format(self._q, names, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: latex(R)
            \mathrm{Frac}(\Bold{Z}[q])[x, y, z]_{q}
        """
        from sage.misc.latex import latex
        names = ", ".join(self.variable_names())
        return "{}[{}]_{{{}}}".format(latex(self.base_ring()), names, self._q)

    @staticmethod
    def _term_key(x):
        r"""
        Compute a key for ``x`` for comparisons.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: elt = (x*y^3*z^2).leading_support()
            sage: R._term_key(elt)
            (6, [2, 3, 1])
        """
        L = x.list()
        L.reverse()
        return (sum(L), L)

    def gen(self, i):
        """
        Return the ``i``-generator of ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.gen(0)
            x
            sage: R.gen(2)
            z
        """
        return self.monomial(self._indices.gen(i))

    @cached_method
    def gens(self):
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.gens()
            (x, y, z)
        """
        return tuple([self.monomial(g) for g in self._indices.gens()])

    @cached_method
    def algebra_generators(self):
        r"""
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.algebra_generators()
            Finite family {'x': x, 'y': y, 'z': z}
        """
        d = {v: self.gen(i) for i,v in enumerate(self.variable_names())}
        return Family(self.variable_names(), d.__getitem__, name="generator")

    @cached_method
    def one_basis(self):
        r"""
        Return the basis index of the element `1`.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.one_basis()
            1
        """
        return self._indices.one()

    def degree_on_basis(self, m):
        r"""
        Return the degree of the monomial index by ``m``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.degree_on_basis(R.one_basis())
            0
            sage: f = (x + y)^3 + z^3
            sage: f.degree()
            3
        """
        return sum(m.list())

    def dimension(self):
        r"""
        Return the dimension of ``self``, which is `\infty`.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y,z> = algebras.qCommutingPolynomials(q)
            sage: R.dimension()
            +Infinity
        """
        return infinity

    @cached_method
    def product_on_basis(self, x, y):
        r"""
        Return the product of two monomials given by ``x`` and ``y``.

        EXAMPLES::

            sage: q = ZZ['q'].fraction_field().gen()
            sage: R.<x,y> = algebras.qCommutingPolynomials(q)
            sage: R.product_on_basis(x.leading_support(), y.leading_support())
            x*y
            sage: R.product_on_basis(y.leading_support(), x.leading_support())
            q*x*y

            sage: x * y
            x*y
            sage: y * x
            q*x*y
            sage: y^2 * x
            q^2*x*y^2
            sage: y * x^2
            q^2*x^2*y
            sage: x * y * x
            q*x^2*y
            sage: y^2 * x^2
            q^4*x^2*y^2
            sage: (x + y)^2
            x^2 + (q+1)*x*y + y^2
            sage: (x + y)^3
            x^3 + (q^2+q+1)*x^2*y + (q^2+q+1)*x*y^2 + y^3
            sage: (x + y)^4
            x^4 + (q^3+q^2+q+1)*x^3*y + (q^4+q^3+2*q^2+q+1)*x^2*y^2 + (q^3+q^2+q+1)*x*y^3 + y^4
        """
        # Special case for multiplying by 1
        if x == self.one_basis():
            return self.monomial(y)
        if y == self.one_basis():
            return self.monomial(x)

        Lx = x.list()
        Ly = y.list()

        # This could be made more efficient
        qpow = sum(exp * sum(Ly[:i]) for i,exp in enumerate(Lx))
        return self.term(x * y, self._q ** qpow)

