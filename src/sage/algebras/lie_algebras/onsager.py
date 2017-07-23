"""
Onsager Algebra

AUTHORS:

- Travis Scrimshaw (2017-07): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.lie_algebras import LieAlgebras
from sage.structure.indexed_generators import IndexedGenerators
from sage.sets.family import Family
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.algebras.lie_algebras.lie_algebra import LieAlgebraWithGenerators

class OnsagerAlgebra(LieAlgebraWithGenerators, IndexedGenerators):
    r"""
    The Onsager (Lie) algebra.

    The Onsager (Lie) algebra `\mathcal{O}` is a Lie algebra with
    generators `A_0, A_1` that satisfy

    .. MATH::

        [A_0, [A_0, [A_0, A_1]]] = -4 [A_0, A_1],
        \qquad
        [A_1, [A_1, [A_1, A_0]]] = -4 [A_1, A_0].

    There exist a basis `\{A_m, G_n \mid m \in \ZZ, n \in \ZZ_{>0}\}`
    for `\mathcal{O}` with structure coefficients

    .. MATH::

        [A_m, A_{m'}] = G_{m-m'},
        \qquad
        [G_n, G_{n'}] = 0,
        \qquad
        [G_n, A_m] = 2A_{m-n} - 2A_{m+n},

    where `m > m'`.

    The Onsager algebra is isomorphic to the subalgebra of the affine
    Lie algebra `\widehat{\mathfrak{sl}}_2[t]` that is invariant under
    the Chevalley involution. In particular, we have

    .. MATH::

        A_i \mapsto f_1 \otimes t^i - e_1 \otimes t^{-i},
        \qquad
        G_i \mapsto h_1 \otimes t^{-i} - h_1 \otimes t^i.

    EXAMPLES:

    We construct the Onsager algebra and do some basic computations::

        sage: O = lie_algebras.OnsagerAlgebra(QQ)
        sage: O.inject_variables()
        Defining A0, A1

    We verify the defining relations::

        sage: O([A0, [A0, [A0, A1]]]) == -4 * O([A0, A1])
        True
        sage: O([A1, [A1, [A1, A0]]]) == -4 * O([A1, A0])
        True

    We check the embedding into `\widehat{\mathfrak{sl}}_2`::

        sage: L = LieAlgebra(QQ, cartan_type=['A',1,1])
        sage: B = L.basis()
        sage: al = RootSystem(['A',1]).root_lattice().simple_root(1)
        sage: ac = al.associated_coroot()
        sage: def emb_A(i): return B[-al,i] - B[al,-i]
        sage: def emb_G(i): return B[ac,i] - B[ac,-i]
        sage: a0 = emb_A(0)
        sage: a1 = emb_A(1)
        sage: L([a0, [a0, [a0, a1]]]) == -4 * L([a0, a1])
        True
        sage: L([a1, [a1, [a1, a0]]]) == -4 * L([a1, a0])
        True

        sage: all(emb_G(n).bracket(emb_A(m)) == 2*emb_A(m-n) - 2*emb_A(m+n)
        ....:     for m in range(-10, 10) for n in range(1,10))
        True
        sage: all(emb_A(m).bracket(emb_A(mp)) == emb_G(m-mp)
        ....:     for m in range(-10,10) for mp in range(m-10, m))
        True

    REFERENCES:

    - [Onsager1944]_
    - [DG1982]_
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: TestSuite(O).run()
        """
        cat = LieAlgebras(R).WithBasis()
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        IndexedGenerators.__init__(self, FiniteEnumeratedSet([0,1]))
        LieAlgebraWithGenerators.__init__(self, R, index_set=self._indices,
                                          names=('A0', 'A1'), category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.OnsagerAlgebra(QQ)
            Onsager algebra over Rational Field
        """
        return "Onsager algebra over {}".format(self.base_ring())

    def _latex_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: latex(O)
            \mathcal{O}_{\Bold{Q}}
        """
        from sage.misc.latex import latex
        return "\\mathcal{{O}}_{{{}}}".format(latex(self.base_ring()))

    def _repr_generator(self, m):
        """
        Return a string representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O._repr_generator((0,-2))
            'A[-2]'
            sage: O._repr_generator((1,4))
            'G[4]'
        """
        if m[0] == 0:
            return 'A[{}]'.format(m[1])
        return 'G[{}]'.format(m[1])

    def _latex_generator(self, m):
        r"""
        Return a LaTeX representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O._latex_generator((0,-2))
            'A_{-2}'
            sage: O._latex_generator((1,4))
            'G_{4}'
        """
        if m[0] == 0:
            return 'A_{{{}}}'.format(m[1])
        return 'G_{{{}}}'.format(m[1])

    # For compatibility with CombinatorialFreeModuleElement
    _repr_term = _repr_generator
    _latex_term = _latex_generator

    @cached_method
    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O.basis()
            Lazy family (Onsager monomial(i))_{i in
             Disjoint union of Family (Integer Ring, Positive integers)}
        """
        from sage.rings.all import ZZ
        from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
        from sage.sets.positive_integers import PositiveIntegers
        I = DisjointUnionEnumeratedSets([ZZ, PositiveIntegers()],
                                        keepkey=True, facade=True)
        return Family(I, self.monomial, name='Onsager monomial')

    @cached_method
    def lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O.lie_algebra_generators()
            Finite family {'A1': A[1], 'A0': A[0]}
        """
        d = {"A0": self.basis()[0,0], "A1": self.basis()[0,1]}
        return Family(self._names, d.__getitem__)

    def bracket_on_basis(self, x, y):
        """
        Return the bracket of basis elements indexed by ``x`` and ``y``
        where ``x < y``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O.bracket_on_basis((1,3), (1,9))  # [G, G]
            0
            sage: O.bracket_on_basis((0,8), (1,13))  # [A, G]
            -2*A[-5] + 2*A[21]
            sage: O.bracket_on_basis((0,-9), (0, 7))  # [A, A]
            -G[16]
        """
        if x[0] == 1:
            # From < property, we have x[1] == 1
            # Therefore, we have [G_n, G_{n'}] = 0
            return self.zero()
        R = self.base_ring()
        if y[0] == 1: # [A_m, G_n] = -(2A_{m-n} - 2A_{m+n})
            d = {(0, x[1]-y[1]): R(-2), (0, x[1]+y[1]): R(2)}
            return self.element_class(self, d)
        # [A_m, A_{m'}] = -G_{m' - m}, where m < m'
        return self.element_class(self, {(1, y[1]-x[1]): -R.one()})

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O.an_element()
            -2*A[-3] + A[2] + 3*G[2]
        """
        B = self.basis()
        return B[0,2] - 2*B[0,-3] + 3*B[1,2]

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: O = lie_algebras.OnsagerAlgebra(QQ)
            sage: O.some_elements()
            [A[0], A[2], A[-1], G[4], -2*A[-3] + A[2] + 3*G[2]]
        """
        B = self.basis()
        return [B[0,0], B[0,2], B[0,-1], B[1,4], self.an_element()]

    Element = LieAlgebraElement

