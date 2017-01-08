r""" 
Cellular Algebra Basis
"""

#*****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.algebras import Algebras

class CellularBasis(CombinatorialFreeModule):
    r"""
    The cellular basis of an algebra.

    INPUT:

    - ``A`` -- the cellular algebra

    EXAMPLES:

    We compute a cellular basis and do some basic computations::

        sage: S = SymmetricGroupAlgebra(QQ, 3)
        sage: C = S.cellular_basis()
        sage: C
        Cellular basis of Symmetric group algebra of order 3
         over Rational Field
        sage: len(C.basis())
        6
        sage: len(S.basis())
        6
        sage: a,b,c,d,e,f = C.basis()
        sage: a
        C([3], [[1, 2, 3]], [[1, 2, 3]])
        sage: c
        C([2, 1], [[1, 3], [2]], [[1, 2], [3]])
        sage: d
        C([2, 1], [[1, 2], [3]], [[1, 3], [2]])
        sage: a * a
        C([3], [[1, 2, 3]], [[1, 2, 3]])
        sage: a * c
        0
        sage: d * c
        C([2, 1], [[1, 2], [3]], [[1, 2], [3]])
        sage: c * d
        C([2, 1], [[1, 3], [2]], [[1, 3], [2]])
        sage: S(a)
        1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1]
         + 1/6*[3, 1, 2] + 1/6*[3, 2, 1]
        sage: S(d)
        1/4*[1, 3, 2] - 1/4*[2, 3, 1] + 1/4*[3, 1, 2] - 1/4*[3, 2, 1]
        sage: B = list(S.basis())
        sage: B[2]
        [2, 1, 3]
        sage: C(B[2])
        -C([1, 1, 1], [[1], [2], [3]], [[1], [2], [3]])
         + C([2, 1], [[1, 2], [3]], [[1, 2], [3]])
         - C([2, 1], [[1, 3], [2]], [[1, 3], [2]])
         + C([3], [[1, 2, 3]], [[1, 2, 3]])
    """
    def __init__(self, A):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: TestSuite(C).run()
        """
        self._algebra = A
        I = [(la, s, t) for la in A.cell_poset()
             for s in A.cell_module_indices(la)
             for t in A.cell_module_indices(la)]

        # TODO: Use instead A.category().Realizations() so
        #   operations are defined by coercion?
        cat = Algebras(A.category().base_ring()).FiniteDimensional().WithBasis().Cellular()
        CombinatorialFreeModule.__init__(self, A.base_ring(), I,
                                         prefix='C', bracket=False,
                                         category=cat)

        # Register coercions
        if A._to_cellular_element is not NotImplemented:
            to_cellular = A.module_morphism(A._to_cellular_element, codomain=self,
                                            category=cat)
        if A._from_cellular_index is NotImplemented:
            from_cellular = ~to_cellular
        else:
            from_cellular = self.module_morphism(A._from_cellular_index, codomain=A,
                                                 category=cat)
            if A._to_cellular_element is NotImplemented:
                to_cellular = ~from_cellular
        to_cellular.register_as_coercion()
        from_cellular.register_as_coercion()

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: S.cellular_basis()
            Cellular basis of Symmetric group algebra of order 3 over Rational Field
        """
        return "Cellular basis of {}".format(self._algebra)

    def _latex_term(self, x):
        r"""
        Return a latex representation of the term indexed by ``x``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: s = Tableau([[1,2],[3]])
            sage: C._latex_term((Partition([2,1]), s, s))
            'C^{...}_{\\left(...\\right)}'
        """
        from sage.misc.latex import latex
        la = x[0]
        m = (x[1], x[2])
        # m contains "non-LaTeXed" strings, use string representation
        sla = latex(la)
        if sla.find('\\text{\\textt') != -1:
            sla = str(la)
        sm = latex(m)
        if sm.find('\\text{\\textt') != -1:
            sm = str(m)
        return "C^{%s}_{%s}"%(sla, sm)

    def cellular_basis_of(self):
        """
        Return the defining algebra of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: C.cellular_basis_of() is S
            True
        """
        return self._algebra

    def cell_poset(self):
        """
        Return the cell poset of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: C.cell_poset()
            Finite poset containing 3 elements
        """
        return self._algebra.cell_poset()

    def cell_module_indices(self, la):
        """
        Return the indices of the cell module of ``self``
        indexed by ``la`` .

        This is the finite set `M(\lambda)`.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: C.cell_module_indices([2,1])
            Standard tableaux of shape [2, 1]
        """
        return self._algebra.cell_module_indices(la)

    def cellular_basis(self):
        """
        Return the cellular basis of ``self``, which is ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: C.cellular_basis() is C
            True
        """
        return self

    @cached_method
    def one(self):
        """
        Return the element `1` in ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: C.one()
            C([1, 1, 1], [[1], [2], [3]], [[1], [2], [3]])
             + C([2, 1], [[1, 2], [3]], [[1, 2], [3]])
             + C([2, 1], [[1, 3], [2]], [[1, 3], [2]])
             + C([3], [[1, 2, 3]], [[1, 2, 3]])
        """
        return self(self._algebra.one())

    @cached_method
    def product_on_basis(self, x, y):
        """
        Return the product of basis indices by ``x`` and ``y``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: la = Partition([2,1])
            sage: s = StandardTableau([[1,2],[3]])
            sage: t = StandardTableau([[1,3],[2]])
            sage: C.product_on_basis((la, s, t), (la, s, t))
            0
        """
        A = self._algebra
        return self(A(self.monomial(x)) * A(self.monomial(y)))

