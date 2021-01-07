r"""
Symplectic vector spaces

AUTHORS:

- Tobias Diez (2020): initial version

TESTS::

    sage: import pytest
    sage: pytest.main(["symplectic_vector_space_test.py"])
    TODO: add output
"""

# *****************************************************************************
#       Copyright (C) 2020 Tobias Diez
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from typing import Optional, Tuple

from sage.manifolds.differentiable.symplectic_form import SymplecticForm, SymplecticFormParal
from sage.manifolds.differentiable.examples.euclidean import EuclideanSpace
from sage.rings.real_mpfr import RR
from sage.categories.manifolds import Manifolds


class SymplecticVectorSpace(EuclideanSpace):
    r"""
    A symplectic vector space is a vector space over `\RR` equipped with a constant symplectic form.
    """

    _symplectic_form: SymplecticForm

    def __init__(self, dimension: int, name: Optional[str] = None, latex_name: Optional[str] = None,
                 coordinates: str ='Cartesian', symbols: Optional[str] = None, symplectic_name: Optional[str] = 'omega',
                 symplectic_latex_name: Optional[str] = None, start_index: int = 1,
                 base_manifold: Optional['SymplecticVectorSpace'] = None, names: Optional[Tuple[str]] = None):
        r"""
        INPUT:

        - ``dimension`` -- dimension of the space over the real field (has to be even)
        - ``name`` -- name (symbol) given to the underlying vector space; if ``None``, the name will be set to ``'V'``
        - ``latex_name`` -- LaTeX symbol to denote the underlying vector space; if ``None``, it is set to ``name``
        - ``coordinates`` -- (default: ``'Cartesian'``) the
            type of coordinates to be initialized at the Euclidean space
            creation; allowed values are

            - ``'Cartesian'`` (canonical coordinates on `\RR^{2n}`)
            - ``'polar'`` for ``dimension=2`` only (see
                :meth:`~sage.manifolds.differentiable.examples.euclidean.EuclideanPlane.polar_coordinates`)

        - ``symbols`` -- the coordinate text symbols and LaTeX symbols, with the same conventions as the
            argument ``coordinates`` in :class:`~sage.manifolds.differentiable.chart.RealDiffChart`, namely
            ``symbols`` is a string of coordinate fields separated by a blank
            space, where each field contains the coordinate's text symbol and
            possibly the coordinate's LaTeX symbol (when the latter is different
            from the text symbol), both symbols being separated by a colon
            (``:``); if ``None``, the symbols will be automatically generated
            according to the value of ``coordinates``
        - ``symplectic_name`` -- name (symbol) given to the symplectic form
        - ``symplecic_latex_name`` -- LaTeX symbol to denote the symplectic form;
            if none is provided, it is set to ``symplectic_name``
        - ``start_index`` -- lower value of the range of
            indices used for "indexed objects" in the vector space, e.g.
            coordinates of a chart
        - ``base_manifold`` -- if not ``None``, the created object is then an open subset
            of ``base_manifold``
        - ``names`` -- (default: ``None``) unused argument, except if
            ``symbols`` is not provided; it must then be a tuple containing
            the coordinate symbols (this is guaranteed if the shortcut operator
            ``<,>`` is used)
            If ``names`` is specified, then ``dimension`` does not have to be specified.

        EXAMPLES:

        Standard symplectic form on `\RR^2`::

            sage: from sage.manifolds.differentiable.examples.symplectic_vector_space import SymplecticVectorSpace
            sage: M.<q, p> = SymplecticVectorSpace(2, symplectic_name='omega')
            sage: omega = M.symplectic_form()
            sage: omega.display()
            omega = -dq/\dp
        """
        dim_half = dimension // 2

        if names is not None and symbols is None:
            symbols = ' '.join(names)

        if symbols is None:
            if dim_half == 1:
                symbols = r"q:q p:p"
            else:
                symbols_list = [f"q{i}:q^{i} p{i}:p_{i}" for i in range(1, dim_half + 1)]
                symbols = ' '.join(symbols_list)

        if name is None:
            name = "V"

        category = Manifolds(RR).Smooth()

        EuclideanSpace.__init__(self, dimension, name, latex_name=latex_name,
                                coordinates=coordinates, symbols=symbols, start_index=start_index,
                                base_manifold=base_manifold, category=category, init_coord_methods=None)

        self._symplectic_form = SymplecticFormParal(self, symplectic_name, symplectic_latex_name)
        for i in range(0, dim_half):
            q_index = 2 * i + 1
            self._symplectic_form.set_comp()[q_index, q_index + 1] = -1

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES:

            sage: from sage.manifolds.differentiable.examples.symplectic_vector_space import SymplecticVectorSpace
            sage: V.<q, p> = SymplecticVectorSpace(2, symplectic_name='omega'); V
            2-dimensional symplectic vector space V
        """
        return f"{self._dim}-dimensional symplectic vector space {self._name}"

    def symplectic_form(self) -> SymplecticForm:
        r"""
        Return the symplectic form.

        EXAMPLES:

        Standard symplectic form on `\RR^2`::

            sage: from sage.manifolds.differentiable.examples.symplectic_vector_space import SymplecticVectorSpace
            sage: M.<q, p> = SymplecticVectorSpace(2, symplectic_name='omega')
            sage: omega = M.symplectic_form()
            sage: omega.display()
            omega = -dq/\dp
        """
        return self._symplectic_form
