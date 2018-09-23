r"""
The Baker-Campbell-Hausdorff formula

AUTHORS:

- Eero Hakavuori (2018-09-23): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.algebras.lie_algebras.lie_algebra import LieAlgebra
from sage.arith.misc import bernoulli
from sage.categories.lie_algebras import LieAlgebras
from sage.combinat.integer_vector import IntegerVectorsConstraints
from sage.functions.other import factorial
from sage.rings.rational_field import QQ
from sage.structure.element import canonical_coercion


class BCH_iterator:
    r"""
    An iterator which returns successive terms of the
    Baker-Campbell-Hausdorff formula.

    INPUT:

    - ``X`` -- (optional) an element of a Lie algebra
    - ``Y`` -- (optional) an element of a Lie algebra

    The BCH formula is an expression for `\log(\exp(X)\exp(Y))` as a sum of Lie
    brackets of ``X`` and ``Y`` with rational coefficients. In arbitrary Lie
    algebras, the infinite sum is only guaranteed to converge for ``X`` and
    ``Y`` close to zero.

    If the elements ``X`` and ``Y`` are not given, then the iterator will
    return successive terms of the abstract BCH formula, i.e., the BCH formula
    for the generators of the free Lie algebra on 2 generators.

    If the Lie algebra containing ``X`` and ``Y`` is not nilpotent, the
    iterator will output infinitely many elements. If the Lie algebra is
    nilpotent, the number of elements outputted is equal to the nilpotency step.

    EXAMPLES:

    The terms of the abstract BCH formula up to fifth order brackets::

        sage: from sage.algebras.lie_algebras.bch import BCH_iterator
        sage: bch = BCH_iterator()
        sage: bch.next()
        X + Y
        sage: bch.next()
        1/2*[X, Y]
        sage: bch.next()
        1/12*[X, [X, Y]] + 1/12*[[X, Y], Y]
        sage: bch.next()
        1/24*[X, [[X, Y], Y]]
        sage: bch.next()
        -1/720*[X, [X, [X, [X, Y]]]] + 1/180*[X, [X, [[X, Y], Y]]]
        + 1/360*[[X, [X, Y]], [X, Y]] + 1/180*[X, [[[X, Y], Y], Y]]
        + 1/120*[[X, Y], [[X, Y], Y]] - 1/720*[[[[X, Y], Y], Y], Y]

    For nilpotent Lie algebras the BCH formula only has finitely many terms::

        sage: L = LieAlgebra(QQ, 2, step=3)
        sage: L.inject_variables()
        Defining X_1, X_2, X_12, X_112, X_122
        sage: [Z for Z in BCH_iterator(X_1, X_2)]
        [X_1 + X_2, 1/2*X_12, 1/12*X_112 + 1/12*X_122]
        sage: [Z for Z in BCH_iterator(X_1 + X_2, X_12)]
        [X_1 + X_2 + X_12, 1/2*X_112 - 1/2*X_122, 0]

    The elements ``X`` and ``Y`` don't need to be elements of the same Lie
    algebra if there is a coercion from one to the other::

        sage: L = LieAlgebra(QQ, 3, step=2)
        sage: L.inject_variables()
        Defining X_1, X_2, X_3, X_12, X_13, X_23
        sage: S = L.subalgebra(X_1, X_2)
        sage: bch1 = [Z for Z in BCH_iterator(S(X_1), S(X_2))]; bch1
        [X_1 + X_2, 1/2*X_12]
        sage: bch1[0].parent() == S
        True
        sage: bch2 = [Z for Z in BCH_iterator(S(X_1), X_3)]; bch2
        [X_1 + X_3, 1/2*X_13]
        sage: bch2[0].parent() == L
        True

    The BCH formula requires a coercion from the rationals::

        sage: L.<X,Y,Z> = LieAlgebra(ZZ, 2, step=2)
        sage: bch = BCH_iterator(X, Y)
        Traceback (most recent call last):
        ...
        TypeError: the BCH formula is not well defined since Integer Ring has no coercion from Rational Field

    TESTS:

    Compare to the BCH formula up to degree 5 given by wikipedia::

        sage: from sage.algebras.lie_algebras.bch import BCH_iterator
        sage: bch = BCH_iterator()
        sage: L = bch._lie_algebra
        sage: computed_BCH = L.sum(bch.next() for k in range(5))
        sage: X, Y = L.graded_basis(1)
        sage: wikiBCH = X + Y + 1/2*L[X,Y] + 1/12*(L[X,[X,Y]] + L[Y,[Y,X]])
        sage: wikiBCH += -1/24*L[Y,[X,[X,Y]]]
        sage: wikiBCH += -1/720*(L[Y,[Y,[Y,[Y,X]]]] + L[X,[X,[X,[X,Y]]]])
        sage: wikiBCH += 1/360*(L[X,[Y,[Y,[Y,X]]]] + L[Y,[X,[X,[X,Y]]]])
        sage: wikiBCH += 1/120*(L[Y,[X,[Y,[X,Y]]]] + L[X,[Y,[X,[Y,X]]]])
        sage: computed_BCH == wikiBCH
        True

    ALGORITHM:

    The BCH formula `\log(\exp(X)\exp(Y)) = \sum_k Z_k` is computed starting
    from `Z_1 = X + Y`, by the recursion

    .. MATH::

        (m+1)Z_{m+1} =  \frac{1}{2}[X - Y, Z_m] + \sum_{2\leq 2p \leq m}\frac{B_{2p}}{(2p)!}\sum_{k_1+\dots+k_{2p}=m}[Z_{k_1}, [\dots [Z_{k_{2p}}, X + Y]\dots],

    where `B_{2p}` are the Bernoulli numbers, see Lemma 2.15.3. in [Var1984]_.

    .. WARNING::

        The time needed to compute each successive term increases exponentially.
        For example on one machine iterating through `Z_{11},...,Z_{18}` for a
        free Lie algebra, computing each successive term took 4-5 times longer,
        going from 0.1s for `Z_{11}` to 21 minutes for `Z_{18}`.
    """

    def __init__(self, X=None, Y=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: from sage.algebras.lie_algebras.bch import BCH_iterator
            sage: BCH_iterator()
            Iterator for BCH(X, Y) in Free Lie algebra generated by (X, Y) over Rational Field in the Lyndon basis
        """
        if X is None or Y is None:
            self._lie_algebra = LieAlgebra(QQ, ['X', 'Y']).Lyndon()
            X, Y = self._lie_algebra.lie_algebra_generators()
        else:
            X, Y = canonical_coercion(X, Y)
            self._lie_algebra = X.parent()

        R = self._lie_algebra.base_ring()
        if not R.has_coerce_map_from(QQ):
            raise TypeError("the BCH formula is not well defined since %s "
                            "has no coercion from %s" % (R, QQ))

        self._xdif = X - Y
        self._Z = [0, X + Y]  # 1-based indexing for convenience
        self._m = 0

    def __iter__(self):
        r"""
        Return the iterator of terms of the BCH formula.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.bch import BCH_iterator
            sage: bch = BCH_iterator()
            sage: l = []
            sage: for k, Z in enumerate(bch):
            ....:     if k==3: break
            ....:     l.append(Z)
            sage: l
            [X + Y, 1/2*[X, Y], 1/12*[X, [X, Y]] + 1/12*[[X, Y], Y]]
        """
        return self

    def __repr__(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.bch import BCH_iterator
            sage: BCH_iterator()
            Iterator for BCH(X, Y) in Free Lie algebra generated by (X, Y) over Rational Field in the Lyndon basis
            sage: L = lie_algebras.Heisenberg(QQ, 1)
            sage: p, q, z = L.basis()
            sage: BCH_iterator(p + z, q)
            Iterator for BCH(p1 + z, q1) in Heisenberg algebra of rank 1 over Rational Field
        """
        X = (self._Z[1] + self._xdif) / 2
        Y = (self._Z[1] - self._xdif) / 2
        return "Iterator for BCH(%s, %s) in %s" % (X, Y, self._lie_algebra)

    def next(self):
        r"""
        Return the next terms of the BCH formula.

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.bch import BCH_iterator
            sage: bch = BCH_iterator()
            sage: bch.next()
            X + Y
            sage: bch.next()
            1/2*[X, Y]
        """
        self._m += 1
        if self._m == 1:
            return self._Z[1]

        L = self._lie_algebra

        if L in LieAlgebras(L.base_ring()).Nilpotent() and self._m > L.step():
            raise StopIteration

        # apply the recursion formula of [Var1984]
        Zm = QQ(1) / QQ(2 * self._m) * L.bracket(self._xdif, self._Z[-1])
        for p in range(1, (self._m - 1) / 2 + 1):
            partitions = IntegerVectorsConstraints(self._m - 1, 2 * p, min_part=1)
            coeff = QQ(bernoulli(2 * p)) / QQ(self._m * factorial(2 * p))
            for kvec in partitions:
                W = self._Z[1]
                for k in kvec:
                    W = L.bracket(self._Z[k], W)
                Zm += coeff * W

        self._Z.append(Zm)
        return Zm
