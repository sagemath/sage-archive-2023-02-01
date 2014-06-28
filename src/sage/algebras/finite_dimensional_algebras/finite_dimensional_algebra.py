"""
Finite-Dimensional Algebras
"""

#*****************************************************************************
#  Copyright (C) 2011 Johan Bosman <johan.g.bosman@gmail.com>
#  Copyright (C) 2011, 2013 Peter Bruin <peter.bruin@math.uzh.ch>
#  Copyright (C) 2011 Michiel Kosters <kosters@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from finite_dimensional_algebra_element import FiniteDimensionalAlgebraElement
from finite_dimensional_algebra_ideal import FiniteDimensionalAlgebraIdeal

from sage.categories.all import FiniteDimensionalAlgebrasWithBasis
from sage.matrix.constructor import Matrix
from sage.matrix.matrix import is_Matrix
from sage.modules.free_module_element import vector
from sage.rings.ring import Algebra

from sage.misc.cachefunc import cached_method
from functools import reduce


class FiniteDimensionalAlgebra(Algebra):
    """
    Create a finite-dimensional `k`-algebra from a multiplication table.

    INPUT:

    - ``k`` -- a field

    - ``table`` -- a list of matrices

    - ``names`` -- (default: ``'e'``) string; names for the basis
      elements

    - ``assume_associative`` -- (default: ``False``) boolean; if
      ``True``, then methods requiring associativity assume this
      without checking

    - ``category`` -- (default: ``FiniteDimensionalAlgebrasWithBasis(k)``)
      the category to which this algebra belongs

    The list ``table`` must have the following form: there exists a
    finite-dimensional `k`-algebra of degree `n` with basis
    `(e_1, \ldots, e_n)` such that the `i`-th element of ``table`` is the
    matrix of right multiplication by `e_i` with respect to the basis
    `(e_1, \ldots, e_n)`.

    EXAMPLES::

        sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
        sage: A
        Finite-dimensional algebra of degree 2 over Finite Field of size 3

        sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
        sage: B
        Finite-dimensional algebra of degree 3 over Rational Field
    """

    def __init__(self, k, table, names='e', assume_associative=False, category=None):
        """
        TESTS::

            sage: A = FiniteDimensionalAlgebra(QQ, [])
            sage: A
            Finite-dimensional algebra of degree 0 over Rational Field
            sage: type(A)
            <class 'sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra.FiniteDimensionalAlgebra_with_category'>
            sage: TestSuite(A).run()

            sage: B = FiniteDimensionalAlgebra(GF(7), [Matrix([1])])
            sage: B
            Finite-dimensional algebra of degree 1 over Finite Field of size 7
            sage: TestSuite(B).run()

            sage: C = FiniteDimensionalAlgebra(CC, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: C
            Finite-dimensional algebra of degree 2 over Complex Field with 53 bits of precision
            sage: TestSuite(C).run()

            sage: FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]])])
            Traceback (most recent call last):
            ...
            ValueError: input is not a multiplication table

            sage: D.<a,b> = FiniteDimensionalAlgebra(RR, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [-1, 0]])])
            sage: D.gens()
            (a, b)

            sage: E = FiniteDimensionalAlgebra(QQ, [Matrix([0])])
            sage: E.gens()
            (e,)
        """
        n = len(table)
        self._table = [b.base_extend(k) for b in table]
        if not all([is_Matrix(b) and b.dimensions() == (n, n) for b in table]):
            raise ValueError("input is not a multiplication table")
        self._assume_associative = assume_associative
        # No further validity checks necessary!
        if category is None:
            category = FiniteDimensionalAlgebrasWithBasis(k)
        Algebra.__init__(self, base_ring=k, names=names, category=category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TEST::

            sage: FiniteDimensionalAlgebra(RR, [Matrix([1])])._repr_()
            'Finite-dimensional algebra of degree 1 over Real Field with 53 bits of precision'
        """
        return "Finite-dimensional algebra of degree {} over {}".format(self.degree(), self.base_ring())

    def _coerce_map_from_(self, S):
        """
        TESTS::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A.has_coerce_map_from(ZZ)
            True
            sage: A.has_coerce_map_from(GF(3))
            True
            sage: A.has_coerce_map_from(GF(5))
            False
            sage: A.has_coerce_map_from(QQ)
            False
        """
        return S == self or (self.base_ring().has_coerce_map_from(S) and self.is_unitary())

    Element = FiniteDimensionalAlgebraElement

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([0])])
            sage: a = A(0)
            sage: a.parent()
            Finite-dimensional algebra of degree 1 over Rational Field
            sage: A(1)
            Traceback (most recent call last):
            ...
            TypeError: algebra is not unitary

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B(17)
            17*e0 + 17*e2
        """
        return self.element_class(self, x)

    # This is needed because the default implementation
    # assumes that the algebra is unitary.
    from_base_ring = _element_constructor_

    def _Hom_(self, B, category):
        """
        Construct a homset of ``self`` and ``B``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([1])])
            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A._Hom_(B, A.category())
            Set of Homomorphisms from Finite-dimensional algebra of degree 1 over Rational Field to Finite-dimensional algebra of degree 2 over Rational Field
        """
        if isinstance(B, FiniteDimensionalAlgebra):
            category = FiniteDimensionalAlgebrasWithBasis(self.base_ring()).or_subcategory(category)
            from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_morphism import FiniteDimensionalAlgebraHomset
            return FiniteDimensionalAlgebraHomset(self, B, category=category)
        return super(FiniteDimensionalAlgebra, self)._Hom_(B, category)

    def ngens(self):
        """
        Return the number of generators of ``self``, i.e., the degree
        of ``self`` over its base field.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A.ngens()
            2
        """
        return len(self._table)

    degree = ngens

    @cached_method
    def gen(self, i):
        """
        Return the `i`-th basis element of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A.gen(0)
            e0
        """
        return self.element_class(self, [j == i for j in xrange(self.ngens())])

    def basis(self):
        """
        Return a list of the basis elements of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A.basis()
            [e0, e1]
        """
        return list(self.gens())

    def __iter__(self):
        """
        Iterates over the elements of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: list(A)
            [0, e0, 2*e0, e1, e0 + e1, 2*e0 + e1, 2*e1, e0 + 2*e1, 2*e0 + 2*e1]

        This is used in the :class:`Testsuite`'s when ``self`` is
        finite.
        """
        if not self.is_finite():
            raise NotImplementedError("object does not support iteration")
        V = self.zero_element().vector().parent()
        for v in V:
            yield self(v)

    def _ideal_class_(self, n=0):
        """
        Return the ideal class of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A._ideal_class_()
            <class 'sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_ideal.FiniteDimensionalAlgebraIdeal'>
        """
        return FiniteDimensionalAlgebraIdeal

    def table(self):
        """
        Return the multiplication table of ``self``, as a list of
        matrices for right multiplication by the basis elements.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A.table()
            [
            [1 0]  [0 1]
            [0 1], [0 0]
            ]
        """
        return self._table

    @cached_method
    def left_table(self):
        """
        Return the list of matrices for left multiplication by the
        basis elements.

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]), Matrix([[0,1],[-1,0]])])
            sage: B.left_table()
            [
            [1 0]  [ 0  1]
            [0 1], [-1  0]
            ]
        """
        B = self.table()
        n = self.degree()
        return [Matrix([B[j][i] for j in xrange(n)]) for i in xrange(n)]

    def __cmp__(self, other):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: B = FiniteDimensionalAlgebra(GF(5), [Matrix([0])])
            sage: cmp(A, A)
            0
            sage: cmp(B, B)
            0
            sage: cmp(A, B)
            1
        """
        if not isinstance(other, FiniteDimensionalAlgebra):
            return cmp(type(self), type(other))
        if self.base_ring() == other.base_ring():
            return cmp(self.table(), other.table())
        else:
            return 1

    def base_extend(self, F):
        """
        Return ``self`` base changed to ``F``.

        EXAMPLES::

            sage: C = FiniteDimensionalAlgebra(GF(2), [Matrix([1])])
            sage: k.<y> = GF(4)
            sage: C.base_extend(k)
            Finite-dimensional algebra of degree 1 over Finite Field in y of size 2^2
        """
        # Base extension of the multiplication table is done by __init__.
        return FiniteDimensionalAlgebra(F, self.table())

    @cached_method
    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(7), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [2, 3]])])
            sage: A.cardinality()
            49

            sage: B = FiniteDimensionalAlgebra(RR, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [2, 3]])])
            sage: B.cardinality()
            +Infinity

            sage: C = FiniteDimensionalAlgebra(RR, [])
            sage: C.cardinality()
            1
        """
        n = self.degree()
        if n == 0:
            return 1
        return self.base_ring().cardinality() ** n

    def ideal(self, gens=None, given_by_matrix=False):
        """
        Return the right ideal of ``self`` generated by ``gens``.

        INPUT:

        - ``A`` -- a :class:`FiniteDimensionalAlgebra`

        - ``gens`` -- (default: None) - either an element of ``A`` or a
          list of elements of ``A``, given as vectors, matrices, or
          FiniteDimensionalAlgebraElements.  If ``given_by_matrix`` is ``True``, then
          ``gens`` should instead be a matrix whose rows form a basis
          of an ideal of ``A``.

        - ``given_by_matrix`` -- boolean (default: ``False``) - if
          ``True``, no checking is done

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A.ideal(A([1,1]))
            Ideal (e0 + e1) of Finite-dimensional algebra of degree 2 over Finite Field of size 3
        """
        return self._ideal_class_()(self, gens=gens,
                                    given_by_matrix=given_by_matrix)

    @cached_method
    def is_associative(self):
        """
        Return ``True`` if ``self`` is associative.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]), Matrix([[0,1],[-1,0]])])
            sage: A.is_associative()
            True

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,1]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,1], [0,0,0], [1,0,0]])])
            sage: B.is_associative()
            False

            sage: e = B.basis()
            sage: (e[1]*e[2])*e[2]==e[1]*(e[2]*e[2])
            False
        """
        B = self.table()
        n = self.degree()
        for i in xrange(n):
            for j in xrange(n):
                eiej = B[j][i]
                if B[i]*B[j] != sum(eiej[k] * B[k] for k in xrange(n)):
                    return False
        return True

    @cached_method
    def is_commutative(self):
        """
        Return ``True`` if ``self`` is commutative.

        EXAMPLES::

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B.is_commutative()
            True

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,0,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,1,0], [0,0,1]])])
            sage: C.is_commutative()
            False
        """
        # Equivalent to self.table() == self.left_table()
        B = self.table()
        for i in xrange(self.degree()):
            for j in xrange(i):
                if B[j][i] != B[i][j]:
                    return False
        return True

    def is_finite(self):
        """
        Return ``True`` if the cardinality of ``self`` is finite.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(7), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [2, 3]])])
            sage: A.is_finite()
            True

            sage: B = FiniteDimensionalAlgebra(RR, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [2, 3]])])
            sage: B.is_finite()
            False

            sage: C = FiniteDimensionalAlgebra(RR, [])
            sage: C.is_finite()
            True
        """
        return self.degree() == 0 or self.base_ring().is_finite()

    @cached_method
    def is_unitary(self):
        """
        Return ``True`` if ``self`` has a two-sided multiplicative
        identity element.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(QQ, [])
            sage: A.is_unitary()
            True

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]), Matrix([[0,1], [-1,0]])])
            sage: B.is_unitary()
            True

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[0,0], [0,0]]), Matrix([[0,0], [0,0]])])
            sage: C.is_unitary()
            False

            sage: D = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]), Matrix([[1,0], [0,1]])])
            sage: D.is_unitary()
            False

        .. NOTE::

            If a finite-dimensional algebra over a field admits a left identity,
            then this is the unique left identity, and it is also a
            right identity.
        """
        k = self.base_ring()
        n = self.degree()
        # B is obtained by concatenating the elements of
        # self.table(), and v by concatenating the rows of
        # the n times n identity matrix.
        B = reduce(lambda x, y: x.augment(y),
                   self.table(), Matrix(k, n, 0))
        v = vector(Matrix.identity(k, n).list())
        try:
            self._one = B.solve_left(v)
            return True
        except ValueError:
            return False

    def is_zero(self):
        """
        Return ``True`` if ``self`` is the zero ring.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(QQ, [])
            sage: A.is_zero()
            True

            sage: B = FiniteDimensionalAlgebra(GF(7), [Matrix([0])])
            sage: B.is_zero()
            False
        """
        return self.degree() == 0

    def one(self):
        """
        Return the multiplicative identity element of ``self``, if it
        exists.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(QQ, [])
            sage: A.one()
            0

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0], [0,1]]), Matrix([[0,1], [-1,0]])])
            sage: B.one()
            e0

            sage: C = FiniteDimensionalAlgebra(QQ, [Matrix([[0,0], [0,0]]), Matrix([[0,0], [0,0]])])
            sage: C.one()
            Traceback (most recent call last):
            ...
            TypeError: algebra is not unitary
        """
        if not self.is_unitary():
            raise TypeError("algebra is not unitary")
        else:
            return self(self._one)

    def random_element(self, *args, **kwargs):
        """
        Return a random element of ``self``.

        Optional input parameters are propagated to the ``random_element``
        method of the underlying :class:`VectorSpace`.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A.random_element()  # random
            e0 + 2*e1

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B.random_element(num_bound=1000)  # random
            215/981*e0 + 709/953*e1 + 931/264*e2
        """
        return self(self.zero_element().vector().parent().random_element(*args, **kwargs))

    def _is_valid_homomorphism_(self, other, im_gens):
        """
        TESTS::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([1])])
            sage: Hom(A, B)(Matrix([[1], [0]]))
            Morphism from Finite-dimensional algebra of degree 2 over Rational Field to Finite-dimensional algebra of degree 1 over Rational Field given by matrix
            [1]
            [0]
            sage: Hom(B, A)(Matrix([[1, 0]]))
            Morphism from Finite-dimensional algebra of degree 1 over Rational Field to Finite-dimensional algebra of degree 2 over Rational Field given by matrix
            [1 0]
            sage: H = Hom(A, A)
            sage: H(Matrix.identity(QQ, 2))
            Morphism from Finite-dimensional algebra of degree 2 over Rational Field to Finite-dimensional algebra of degree 2 over Rational Field given by matrix
            [1 0]
            [0 1]
            sage: H(Matrix([[1, 0], [0, 0]]))
            Morphism from Finite-dimensional algebra of degree 2 over Rational Field to Finite-dimensional algebra of degree 2 over Rational Field given by matrix
            [1 0]
            [0 0]
            sage: H(Matrix([[1, 0], [1, 1]]))
            Traceback (most recent call last):
            ...
            ValueError: relations do not all (canonically) map to 0 under map determined by images of generators.
            sage: Hom(B, B)(Matrix([[2]]))
            Traceback (most recent call last):
            ...
            ValueError: relations do not all (canonically) map to 0 under map determined by images of generators.
        """
        assert len(im_gens) == self.degree()

        B = self.table()
        for i,gi in enumerate(im_gens):
            for j,gj in enumerate(im_gens):
                eiej = B[j][i]
                if (sum([other(im_gens[k]) * v for k,v in enumerate(eiej)])
                        != other(gi) * other(gj)):
                    return False
        return True

    def quotient_map(self, ideal):
        """
        Return the quotient of ``self`` by ``ideal``.

        INPUT:

        - ``ideal`` -- a ``FiniteDimensionalAlgebraIdeal``

        OUTPUT:

        - :class:`~sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_morphism.FiniteDimensionalAlgebraMorphism`;
          the quotient homomorphism

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: q0 = A.quotient_map(A.zero_ideal())
            sage: q0
            Morphism from Finite-dimensional algebra of degree 2 over Finite Field of size 3 to Finite-dimensional algebra of degree 2 over Finite Field of size 3 given by matrix
            [1 0]
            [0 1]
            sage: q1 = A.quotient_map(A.ideal(A.gen(1)))
            sage: q1
            Morphism from Finite-dimensional algebra of degree 2 over Finite Field of size 3 to Finite-dimensional algebra of degree 1 over Finite Field of size 3 given by matrix
            [1]
            [0]
        """
        k = self.base_ring()
        f = ideal.basis_matrix().transpose().kernel().basis_matrix().echelon_form().transpose()
        pivots = f.pivot_rows()
        table = []
        for p in pivots:
            v = vector(k, self.degree())
            v[p] = 1
            v = self.element_class(self, v)
            table.append(f.solve_right(v.matrix() * f))
        B = FiniteDimensionalAlgebra(k, table)
        return self.hom(f, codomain=B, check=False)

    def maximal_ideal(self):
        """
        Compute the maximal ideal of the local algebra ``self``.

        .. NOTE::

            ``self`` must be unitary, commutative, associative and local
            (have a unique maximal ideal).

        OUTPUT:

        - :class:`~sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_ideal.FiniteDimensionalAlgebraIdeal`;
          the unique maximal ideal of ``self``.  If ``self`` is not a local
          algebra, a ``ValueError`` is raised.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A.maximal_ideal()
            Ideal (0, e1) of Finite-dimensional algebra of degree 2 over Finite Field of size 3

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B.maximal_ideal()
            Traceback (most recent call last):
            ...
            ValueError: algebra is not local
        """
        if self.degree() == 0:
            raise ValueError("the zero algebra is not local")
        if not(self.is_unitary() and self.is_commutative()
               and (self._assume_associative or self.is_associative())):
            raise TypeError("algebra must be unitary, commutative and associative")
        gens = []
        for x in self.gens():
            f = x.characteristic_polynomial().factor()
            if len(f) != 1:
                raise ValueError("algebra is not local")
            if f[0][1] > 1:
                gens.append(f[0][0](x))
        return FiniteDimensionalAlgebraIdeal(self, gens)

    def primary_decomposition(self):
        """
        Return the primary decomposition of ``self``.

        .. NOTE::

            ``self`` must be unitary, commutative and associative.

        OUTPUT:

        - a list consisting of the quotient maps ``self`` -> `A`,
          with `A` running through the primary factors of ``self``

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A.primary_decomposition()
            [Morphism from Finite-dimensional algebra of degree 2 over Finite Field of size 3 to Finite-dimensional algebra of degree 2 over Finite Field of size 3 given by matrix [1 0]
            [0 1]]

            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1,0,0], [0,1,0], [0,0,0]]), Matrix([[0,1,0], [0,0,0], [0,0,0]]), Matrix([[0,0,0], [0,0,0], [0,0,1]])])
            sage: B.primary_decomposition()
            [Morphism from Finite-dimensional algebra of degree 3 over Rational Field to Finite-dimensional algebra of degree 1 over Rational Field given by matrix [0]
            [0]
            [1], Morphism from Finite-dimensional algebra of degree 3 over Rational Field to Finite-dimensional algebra of degree 2 over Rational Field given by matrix [1 0]
            [0 1]
            [0 0]]
        """
        k = self.base_ring()
        n = self.degree()
        if n == 0:
            return []
        if not (self.is_unitary() and self.is_commutative()
                and (self._assume_associative or self.is_associative())):
            raise TypeError("algebra must be unitary, commutative and associative")
        # Start with the trivial decomposition of self.
        components = [Matrix.identity(k, n)]
        for b in self.table():
            # Use the action of the basis element b to refine our
            # decomposition of self.
            components_new = []
            for c in components:
                # Compute the matrix of b on the component c, find its
                # characteristic polynomial, and factor it.
                b_c = c.solve_left(c * b)
                fact = b_c.characteristic_polynomial().factor()
                if len(fact) == 1:
                    components_new.append(c)
                else:
                    for f in fact:
                        h, a = f
                        e = h(b_c) ** a
                        ker_e = e.kernel().basis_matrix()
                        components_new.append(ker_e * c)
            components = components_new
        quotients = []
        for i in range(len(components)):
            I = Matrix(k, 0, n)
            for j,c in enumerate(components):
                if j != i:
                    I = I.stack(c)
            quotients.append(self.quotient_map(self.ideal(I, given_by_matrix=True)))
        return quotients

    def maximal_ideals(self):
        """
        Return a list consisting of all maximal ideals of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A.maximal_ideals()
            [Ideal (e1) of Finite-dimensional algebra of degree 2 over Finite Field of size 3]

            sage: B = FiniteDimensionalAlgebra(QQ, [])
            sage: B.maximal_ideals()
            []
        """
        P = self.primary_decomposition()
        return [f.inverse_image(f.codomain().maximal_ideal()) for f in P]

