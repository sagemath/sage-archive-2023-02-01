r"""
Nilpotent Lie algebras

AUTHORS:

- Eero Hakavuori (2018-08-16): initial version
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

from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
from sage.categories.lie_algebras import LieAlgebras


class NilpotentLieAlgebra_dense(LieAlgebraWithStructureCoefficients):
    r"""
    A nilpotent Lie algebra `L` over a base ring.

    INPUT:

    - ``R`` -- the base ring
    - ``s_coeff`` -- a dictionary of structural coefficients
    - ``names`` -- (default:``None``) list of strings to use as names of basis
      elements; if ``None``, the names will be inferred from the structural
      coefficients
    - ``index_set`` -- (default:``None``) list of hashable and comparable
      elements to use for indexing
    - ``step`` -- (optional) an integer; the nilpotency step of the
      Lie algebra if known; otherwise it will be computed when needed
    - ``category`` -- (optional) a subcategory of finite dimensional
      nilpotent Lie algebras with basis

    EXAMPLES:

    The input to a :class:`NilpotentLieAlgebra_dense` should be of the
    same form as to a
    :class:`~sage.algebras.lie_algebras.structure_coefficients.LieAlgebraWithStructureCoefficients`::

        sage: L.<X,Y,Z,W> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
        sage: L
        Nilpotent Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
        sage: L[X, Y]
        Z
        sage: L[X, W]
        0

    If the parameter ``names`` is omitted, then the terms appearing in the
    structural coefficients are used as names::

        sage: L = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True); L
        Nilpotent Lie algebra on 3 generators (X, Y, Z) over Rational Field

    TESTS::

        sage: L = LieAlgebra(QQ, {('X','Y'): {'Z': 1},
        ....:                     ('X','Z'): {'W': 1},
        ....:                     ('Y','Z'): {'T': 1}}, nilpotent=True)
        sage: TestSuite(L).run()
    """
    @staticmethod
    def __classcall_private__(cls, R, s_coeff, names=None, index_set=None, **kwds):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES:

        If the variable order is specified, the order of structural
        coefficients does not matter::

            sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra_dense
            sage: L1.<x,y,z> = NilpotentLieAlgebra_dense(QQ, {('x','y'): {'z': 1}})
            sage: L2.<x,y,z> = NilpotentLieAlgebra_dense(QQ, {('y','x'): {'z': -1}})
            sage: L1 is L2
            True

        If the variables are implicitly defined by the structural coefficients,
        the ordering may be different and the Lie algebras will be considered
        different::

            sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra_dense
            sage: L1 = NilpotentLieAlgebra_dense(QQ, {('x','y'): {'z': 1}})
            sage: L2 = NilpotentLieAlgebra_dense(QQ, {('y','x'): {'z': -1}})
            sage: L1
            Nilpotent Lie algebra on 3 generators (x, y, z) over Rational Field
            sage: L2
            Nilpotent Lie algebra on 3 generators (y, x, z) over Rational Field
            sage: L1 is L2
            False
        """
        if not names:
            # extract names from structural coefficients
            names = []
            for (X, Y), d in s_coeff.items():
                if X not in names: names.append(X)
                if Y not in names: names.append(Y)
                for k in d:
                    if k not in names: names.append(k)

        from sage.structure.indexed_generators import standardize_names_index_set
        names, index_set = standardize_names_index_set(names, index_set)
        s_coeff = LieAlgebraWithStructureCoefficients._standardize_s_coeff(
            s_coeff, index_set)

        return super(NilpotentLieAlgebra_dense, cls).__classcall__(
            cls, R, s_coeff, names, index_set, **kwds)

    def __init__(self, R, s_coeff, names, index_set, step=None,
                 category=None, **kwds):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L.<X,Y,Z,W> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
            sage: TestSuite(L).run()
        """
        if step:
            self._step = step

        cat = LieAlgebras(R).FiniteDimensional().WithBasis().Nilpotent()
        cat = cat.or_subcategory(category)
        LieAlgebraWithStructureCoefficients.__init__(self, R, s_coeff,
                                                     names, index_set,
                                                     category=cat, **kwds)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L.<X,Y,Z,W> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
            sage: L
            Nilpotent Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
        """
        return "Nilpotent %s" % (super(NilpotentLieAlgebra_dense, self)._repr_())
