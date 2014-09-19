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
from sage.algebras.algebra import Algebra
from sage.algebras.lie_algebras.lie_algebra import InfinitelyGeneratedLieAlgebra, \
    LieAlgebraFromAssociative, FinitelyGeneratedLieAlgebra
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement, LieGenerator
from sage.categories.lie_algebras import LieAlgebras
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.infinity import infinity
from sage.sets.family import Family

class HeisenbergGenerator(LieGenerator):
    """
    Generator for the Heisenberg algebra.
    """
    def __init__(self, name, i=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.heisenberg import HeisenbergGenerator
            sage: gen = HeisenbergGenerator('p', 2)
            sage: TestSuite(gen).run()
        """
        self._i = i
        LieGenerator.__init__(self, name)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.heisenberg import HeisenbergGenerator
            sage: HeisenbergGenerator('p', 2)
            p2
            sage: HeisenbergGenerator('d')
            d
        """
        if self._i is None:
            return self._name
        return self._name + repr(self._i)

    def _latex_(self):
        r"""
        Return a `LaTeX` representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.heisenberg import HeisenbergGenerator
            sage: latex(HeisenbergGenerator('p', 2))
            p_{2}
            sage: latex(HeisenbergGenerator('d'))
            d
        """
        if self._i is None:
            return self._name
        return self._name + "_{{{}}}".format(self._i)

    def __eq__(self, rhs):
        """
        Compare equals.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.heisenberg import HeisenbergGenerator
            sage: p2 = HeisenbergGenerator('p', 2)
            sage: p2 == HeisenbergGenerator('p', 2)
            True
            sage: p2 == HeisenbergGenerator('q', 2)
            False
            sage: p2 == HeisenbergGenerator('p', 4)
            False
        """
        return isinstance(rhs, LieGenerator) and self._name == rhs._name \
                and self._i == rhs._i

    def __ne__(self, rhs):
        """
        Compare equals.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.heisenberg import HeisenbergGenerator
            sage: p2 = HeisenbergGenerator('p', 2)
            sage: p2 != HeisenbergGenerator('p', 2)
            False
            sage: p2 != HeisenbergGenerator('q', 2)
            True
            sage: p2 != HeisenbergGenerator('p', 4)
            True
        """
        return not self.__eq__(rhs)

    def __lt__(self, rhs):
        """
        Compare less than.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.heisenberg import HeisenbergGenerator
            sage: p2 = HeisenbergGenerator('p', 2)
            sage: q2 = HeisenbergGenerator('q', 2)
            sage: p3 = HeisenbergGenerator('p', 3)
            sage: p2 < p3
            True
            sage: p2 < q2
            True
        """
        if isinstance(rhs, HeisenbergGenerator):
            return self._name < rhs._name or (self._name == rhs._name and self._i < rhs._i)
        return False

class HeisenbergAlgebra_abstract:
    """
    The common methods for the (non-matrix) Heisenberg algebras.
    """
    def p(self, i):
        """
        The generator `p_i`.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.p(2)
            p2
        """
        return self.element_class(self, {HeisenbergGenerator('p', i): self.base_ring().one()})

    def q(self, i):
        """
        The generator `q_i`.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.q(2)
            q2
        """
        return self.element_class(self, {HeisenbergGenerator('q', i): self.base_ring().one()})

    def z(self):
        """
        The generator `z`.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.z()
            z
        """
        return self.element_class(self, {HeisenbergGenerator('z'): self.base_ring().one()})

    def bracket_on_basis(self, x, y):
        """
        Return the bracket of basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.heisenberg import HeisenbergGenerator
            sage: H = lie_algebras.Heisenberg(QQ, 3)
            sage: p1 = HeisenbergGenerator('p', 1)
            sage: q1 = HeisenbergGenerator('q', 1)
            sage: H.bracket_on_basis(p1, q1)
            z
        """
        if x._name == 'p' and y._name == 'q' and x._i == y._i:
            return self.z()
        return self.zero()

    Element = LieAlgebraElement

class HeisenbergAlgebra(HeisenbergAlgebra_abstract, FinitelyGeneratedLieAlgebra):
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
        self._n = n
        names = ['p%s'%i for i in range(1,n+1)] + ['q%s'%i for i in range(1,n+1)] + ['z']
        FinitelyGeneratedLieAlgebra.__init__(self, R, names,
            category=LieAlgebras(R).FiniteDimensional().WithBasis())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.Heisenberg(QQ, 3)
            Heisenberg algebra of rank 3 over Rational Field
        """
        return "Heisenberg algebra of rank {0} over {1}".format(self._n, self.base_ring())

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 2)
            sage: H.gens()
            (p1, p2, q1, q2, z)
        """
        L  = [self.p(i) for i in range(1, self._n+1)]
        L += [self.q(i) for i in range(1, self._n+1)]
        L += [self.z()]
        return tuple(L)

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 2)
            sage: H.gen(2)
            q1
            sage: H.gen(4)
            z
        """
        return self.gens()[i]

    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 1)
            sage: H.algebra_generators()
            Finite family {'q1': q1, 'p1': p1}
        """
        return Family({self.variable_names()[i]: x for i,x in enumerate(self.gens()[:-1])})

    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: H = lie_algebras.Heisenberg(QQ, 1)
            sage: H.basis()
            Finite family {'q1': q1, 'p1': p1, 'z': z}
        """
        return Family({self.variable_names()[i]: x for i,x in enumerate(self.gens())})

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
        InfinitelyGeneratedLieAlgebra.__init__(self, R, category=LieAlgebras(R).WithBasis())

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
            p2 - 1/2*q2 + z
        """
        c = self.base_ring().an_element()
        return self.p(2) - c * self.q(2) + self.z()

    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, oo)
            sage: L.basis()
            Disjoint union of Family
             (Lazy family (<bound method InfiniteHeisenbergAlgebra_with_category.p of The infinite Heisenberg algebra over Rational Field>(i))_{i in Positive integers},
              Lazy family (<bound method InfiniteHeisenbergAlgebra_with_category.q of The infinite Heisenberg algebra over Rational Field>(i))_{i in Positive integers},
              Family (z,))
        """
        from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
        from sage.sets.positive_integers import PositiveIntegers
        from sage.sets.family import Family
        p = Family(PositiveIntegers(), self.p)
        q = Family(PositiveIntegers(), self.q)
        z = Family([self.z()])
        return DisjointUnionEnumeratedSets([p, q, z])

#######################################################
## Finite rank Heisenberg algebra using matrices

class HeisenbergAlgebra_matrix(LieAlgebraFromAssociative, HeisenbergAlgebra):
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
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 2, representation="matrix")
            sage: TestSuite(L).run()
        """
        self._n = n
        MS = MatrixSpace(R, n+2, sparse=True)
        one = R.one()
        p = [MS({(0,i):one}) for i in range(1, n+1)]
        q = [MS({(i,n+1):one}) for i in range(1, n+1)]
        names = ['p%s'%i for i in range(n)] + ['q%s'%i for i in range(n)] + ['z']
        LieAlgebraFromAssociative.__init__(self, R, MS, tuple(p + q + [MS({(0,n+1):one})]), tuple(names))
        self._p = self.gens()[:n]
        self._q = self.gens()[n:2*n]
        self._z = self.gen(2*n)

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
            sage: L.p(0)
            [0 1 0]
            [0 0 0]
            [0 0 0]
        """
        return self._p[i]

    def q(self, i):
        r"""
        Return the generator `q_i`.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation="matrix")
            sage: L.q(0)
            [0 0 0]
            [0 0 1]
            [0 0 0]
        """
        return self._q[i]

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
        return self._z

    def bracket_on_basis(self, x, y):
        """
        Return the bracket of basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation="matrix")
            sage: p = matrix(QQ, 3, 3, {(0,1):1})
            sage: q = matrix(QQ, 3, 3, {(1,2):1})
            sage: L.bracket_on_basis(p, q)
            [0 0 1]
            [0 0 0]
            [0 0 0]
        """
        return self(x*y - y*x)

