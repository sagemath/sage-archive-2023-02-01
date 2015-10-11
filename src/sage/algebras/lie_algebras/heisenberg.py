"""
Heisenberg Algebras

AUTHORS:

- Travis Scrimshaw (2013-08-13): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.misc import repr_lincomb
from sage.structure.indexed_generators import IndexedGenerators

from sage.algebras.algebra import Algebra
from sage.algebras.lie_algebras.lie_algebra import InfinitelyGeneratedLieAlgebra, \
    LieAlgebraFromAssociative, FinitelyGeneratedLieAlgebra
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.categories.lie_algebras import LieAlgebras
from sage.combinat.cartesian_product import CartesianProduct
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.infinity import infinity
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.positive_integers import PositiveIntegers
from sage.sets.set import Set

class HeisenbergAlgebra_abstract(IndexedGenerators):
    """
    The common methods for the (non-matrix) Heisenberg algebras.
    """
    def __init__(self, I):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo) # indirect doctest
        """
        IndexedGenerators.__init__(self, I, prefix='', bracket=False,
                                   latex_bracket=False, string_quotes=False)
    def p(self, i):
        """
        The generator `p_i`.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.p(2)
            p2
        """
        return self.element_class(self, {('p', i): self.base_ring().one()})

    def q(self, i):
        """
        The generator `q_i`.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.q(2)
            q2
        """
        return self.element_class(self, {('q', i): self.base_ring().one()})

    def z(self):
        """
        The generator `z`.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.z()
            z
        """
        return self.element_class(self, {'z': self.base_ring().one()})

    def bracket_on_basis(self, x, y):
        """
        Return the bracket of basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 3)
            sage: p1 = ('p', 1)
            sage: q1 = ('q', 1)
            sage: H.bracket_on_basis(p1, q1)
            z
        """
        if x == 'z' or y == 'z':
            return self.zero()
        if x[0] == 'p' and y[0] == 'q' and x[1] == y[1]:
            return self.z()
        return self.zero()

    def _repr_term(self, m):
        r"""
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 3)
            sage: H._repr_term(('p', 1))
            'p1'
            sage: H._repr_term('z')
            'z'
        """
        if isinstance(m, str):
            return m
        return m[0] + str(m[1])

    def _latex_term(self, m):
        r"""
        Return a string representation of the term indexed by ``m``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 3)
            sage: H._latex_term(('p', 1))
            'p_{1}'
            sage: H._latex_term('z')
            'z'
        """
        if isinstance(m, str):
            return m
        return "%s_{%s}" % m # else it is a tuple of length 2

class HeisenbergAlgebra_fd:
    """
    Common methods for finite-dimensional Heisenberg algebras.
    """
    def __init__(self, n):
        """
        Initialize ``self``.

        INPUT:

        - ``n`` -- the rank

        TESTS::

            sage: H = lie_algebras.Heisenberg(QQ, 3) # indirect doctest
        """
        self._n = n

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 2)
            sage: H.gens()
            (p1, p2, q1, q2)
        """
        L  = [self.p(i) for i in range(1, self._n+1)]
        L += [self.q(i) for i in range(1, self._n+1)]
        #L += [self.z()]
        return tuple(L)

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 2)
            sage: H.gen(0)
            p1
            sage: H.gen(3)
            q2
        """
        return self.gens()[i]

    @cached_method
    def lie_algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 1)
            sage: H.lie_algebra_generators()
            Finite family {'q1': q1, 'p1': p1}
        """
        d = {}
        for i in range(1, self._n+1):
            d['p%s'%i] = self.p(i)
            d['q%s'%i] = self.q(i)
        return Family(d)

    @cached_method
    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 1)
            sage: H.basis()
            Finite family {'q1': q1, 'p1': p1, 'z': z}
        """
        d = {}
        for i in range(1, self._n+1):
            d['p%s'%i] = self.p(i)
            d['q%s'%i] = self.q(i)
        d['z'] = self.z()
        return Family(d)

class HeisenbergAlgebra(HeisenbergAlgebra_fd, HeisenbergAlgebra_abstract,
                        FinitelyGeneratedLieAlgebra):
    """
    A Heisenberg algebra defined using structure coefficients.

    The Heisenberg algebra is the Lie algebra generated by `p_i`, `q_i`, and
    `z` with the following relations:

    .. MATH::

        [p_i, q_j] = \delta_{ij} z, \quad [p_i, z] = [q_i, z] = [p_i, p_j]
        = [q_i, q_j] = 0.

    .. NOTE::

        The relations `[p_i, q_j] = \delta_{ij} z`, `[p_i, z] = 0`, and
        `[q_i, z] = 0` are known as canonical commutation relations. See
        :wikipedia:`Canonical_commutation_relations`.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the rank of the Heisenberg algebra

    REFERENCES:

    - :wikipedia:`Heisenberg_algebra`

    EXAMPLES::

        sage: L = lie_algebras.Heisenberg(QQ, 2)
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 2)
            sage: TestSuite(L).run()
        """
        HeisenbergAlgebra_fd.__init__(self, n)
        names = ['p%s'%i for i in range(1,n+1)] + ['q%s'%i for i in range(1,n+1)] #+ ['z']
        names = tuple(names)
        FinitelyGeneratedLieAlgebra.__init__(self, R, names=names, index_set=names,
            category=LieAlgebras(R).FiniteDimensional().WithBasis())
        HeisenbergAlgebra_abstract.__init__(self, names)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.Heisenberg(QQ, 3)
            Heisenberg algebra of rank 3 over Rational Field
        """
        return "Heisenberg algebra of rank {0} over {1}".format(self._n, self.base_ring())

class InfiniteHeisenbergAlgebra(HeisenbergAlgebra_abstract, InfinitelyGeneratedLieAlgebra):
    """
    The infinite Heisenberg algebra.

    This is the Heisenberg algebra on an infinite number of generators. See
    :class:`HeisenbergAlgebra` for more information.
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: TestSuite(L).run()
        """
        S = CartesianProduct(PositiveIntegers(), ['p','q'])
        cat = LieAlgebras(R).WithBasis()
        InfinitelyGeneratedLieAlgebra.__init__(self, R, index_set=S, category=cat)
        HeisenbergAlgebra_abstract.__init__(self, S)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.Heisenberg(QQ, oo)
            The infinite Heisenberg algebra over Rational Field
        """
        return "The infinite Heisenberg algebra over {}".format(self.base_ring())

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L._an_element_()
            z + p2 + q2 - 1/2*q3
        """
        c = self.base_ring().an_element()
        return self.p(2) + self.q(2) - c * self.q(3) + self.z()

    def lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.lie_algebra_generators()
            Lazy family (generator map(i))_{i in Cartesian product of
                                            Positive integers, ['p', 'q']}
        """
        return Family(self._indices, lambda x: self.monomial(tuple(x)),
                      name='generator map')

    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.basis()
            Lazy family (basis map(i))_{i in Disjoint union of
             Family ({'z'}, Cartesian product of Positive integers, ['p', 'q'])}
        """
        S = CartesianProduct(PositiveIntegers(), ['p','q'])
        I = DisjointUnionEnumeratedSets([Set(['z']), S])
        return Family(I, self.monomial, name="basis map")

#######################################################
## Finite rank Heisenberg algebra using matrices

class HeisenbergAlgebra_matrix(HeisenbergAlgebra_fd, LieAlgebraFromAssociative):
    r"""
    A Heisenberg algebra represented using matrices.

    The Heisenberg algebra over `R` is a Lie algebra which is defined as
    the Lie algebra of the matrices:

    .. MATH::

        \begin{bmatrix}
        0 & p^T & z \\
        0 & 0_n & q \\
        0 & 0 & 0
        \end{bmatrix}

    where `p, q \in R^n` and `0_n` in the `n \times n` zero matrix. It has a
    basis of

    .. MATH::

        \begin{aligned}
        p_i & = \begin{bmatrix}
        0 & e_i^T & 0 \\
        0 & 0_n & 0 \\
        0 & 0 & 0
        \end{bmatrix},
        \\ q_i & = \begin{bmatrix}
        0 & 0 & 0 \\
        0 & 0_n & e_i \\
        0 & 0 & 0
        \end{bmatrix},
        \\ z & = \begin{bmatrix}
        0 & 0 & z \\
        0 & 0_n & 0 \\
        0 & 0 & 0
        \end{bmatrix},
        \end{aligned}

    where `\{e_i\}` is the standard basis.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the size of the matrices

    EXAMPLES::

        sage: L = lie_algebras.Heisenberg(QQ, 1, representation="matrix")
        sage: p = L.p(1)
        sage: q = L.q(1)
        sage: z = L.bracket(p, q); z
        [0 0 1]
        [0 0 0]
        [0 0 0]
        sage: z == L.z()
        True
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 2, representation="matrix")
            sage: TestSuite(L).run()
        """
        HeisenbergAlgebra_fd.__init__(self, n)
        MS = MatrixSpace(R, n+2, sparse=True)
        one = R.one()
        p = tuple(MS({(0,i):one}) for i in range(1, n+1))
        q = tuple(MS({(i,n+1):one}) for i in range(1, n+1))
        names = tuple('p%s'%i for i in range(1,n+1))
        names = names + tuple('q%s'%i for i in range(1,n+1)) + ('z',)
        cat = LieAlgebras(R).FiniteDimensional().WithBasis()
        LieAlgebraFromAssociative.__init__(self, MS, p + q + (MS({(0,n+1):one}),),
                                           names=names, index_set=names, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.Heisenberg(QQ, 3, representation="matrix")
            Heisenberg algebra of rank 3 over Rational Field
        """
        return "Heisenberg algebra of rank {} over {}".format(self._n, self.base_ring())

    def p(self, i):
        r"""
        Return the generator `p_i`.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation="matrix")
            sage: L.p(1)
            [0 1 0]
            [0 0 0]
            [0 0 0]
        """
        return self._gens['p%s'%i]

    def q(self, i):
        r"""
        Return the generator `q_i`.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation="matrix")
            sage: L.q(1)
            [0 0 0]
            [0 0 1]
            [0 0 0]
        """
        return self._gens['q%s'%i]

    def z(self):
        """
        Return the generator `z`.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation="matrix")
            sage: L.z()
            [0 0 1]
            [0 0 0]
            [0 0 0]
        """
        return self._gens['z']

