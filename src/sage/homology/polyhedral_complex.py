# -*- coding: utf-8 -*-
r"""
Finite polyhedral complexes

This module implements the basic structure of finite polyhedral complexes.

A polyhedral complex `PC` is a collection of polyhedra in a certain ambient
space `\RR^n` such that

- If a poyhedron `P` is in `PC`, then all the faces of `P` are in `PC`.

- If polyhedra `P` and `Q` are in `PC`, then `P \cap Q` is either empty or a face of both `P` and `Q`.

In this context, a "polyhedron" means the geometric realization
of a polyhedron. This is in contrast to :mod:`simplicial complex
<sage.homology.simplicial_complex>`, whose cells are abstract simplices.
The concept of a polyhedral complex generalizes that of a **geometric**
simplicial complex.

.. note::

   This class derives from
   :class:`~sage.homology.cell_complex.GenericCellComplex`, and so
   inherits its methods.  Some of those methods are not listed here;
   see the :mod:`Generic Cell Complex <sage.homology.cell_complex>`
   page instead.

AUTHORS:

- Yuan Zhou (2021-04): initial implementation

List of PolyhedralComplex methods
---------------------------------

**Maximal cells and cells**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~PolyhedralComplex.maximal_cells` | Return the dictionary of the maximal cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.maximal_cell_iterator` | Return an iterator over maximal cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.maximal_cells_sorted` | Return the sorted list of all maximal cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.n_maximal_cells` | List the maximal cells of dimension `n` in this polyhedral complex.
    :meth:`~PolyhedralComplex._n_maximal_cells_sorted` | Return the sorted list of maximal cells of dim `n` in this complex.
    :meth:`~PolyhedralComplex.has_maximal_cell` | Return ``True`` if the given cell is a maximal cell in this complex.
    :meth:`~PolyhedralComplex.cells` | Return the dictionary of the cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.cell_iterator` | Return an iterator over cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.cells_sorted` | Return the sorted list of all cells in this polyhedral complex.
    :meth:`~sage.homology.cell_complex.GenericCellComplex.n_cells` | List the cells of dimension `n` in this polyhedral complex.
    :meth:`~PolyhedralComplex._n_cells_sorted` | Return the sorted list of `n`-cells in this polyhedral complex.
    :meth:`~PolyhedralComplex.has_cell` | Return ``True`` if the given cell is in this polyhedral complex.
    :meth:`~PolyhedralComplex.face_poset` | Return the poset of nonempty cells in the polyhedral complex.
    :meth:`~PolyhedralComplex.relative_boundary_cells` | List the maximal cells on the boundary of the polyhedral complex.

**Properties of the polyhedral complex**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~PolyhedralComplex.dimension` | Return the dimension of the polyhedral complex.
    :meth:`~PolyhedralComplex.ambient_dimension` | Return the ambient dimension of the polyhedral complex.
    :meth:`~PolyhedralComplex.is_pure` | Return ``True`` if the polyhedral complex is pure.
    :meth:`~PolyhedralComplex.is_full_dimensional` | Return ``True`` if the polyhedral complex is full dimensional.
    :meth:`~PolyhedralComplex.is_compact` | Return ``True`` if the polyhedral complex is bounded.
    :meth:`~PolyhedralComplex.is_connected` | Return ``True`` if the polyhedral complex is connected.
    :meth:`~PolyhedralComplex.is_subcomplex` | Return ``True`` if this complex is a subcomplex of the other.
    :meth:`~PolyhedralComplex.is_convex` | Return ``True`` if the polyhedral complex is convex.

**New polyhedral complexes from old ones**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~PolyhedralComplex.connected_component` | Return the connected component containing a cell as a subcomplex.
    :meth:`~PolyhedralComplex.connected_components` | Return the connected components of this polyhedral complex.
    :meth:`~PolyhedralComplex.n_skeleton` | Return the `n`-skeleton of this polyhedral complex.
    :meth:`~PolyhedralComplex.stratify` | Return the (pure) subcomplex formed by the maximal cells of dim `n` in this complex.
    :meth:`~PolyhedralComplex.boundary_subcomplex` | Return the boundary subcomplex of this polyhedral complex.
    :meth:`~PolyhedralComplex.product` | Return the (Cartesian) product of this polyhedral complex with another one.
    :meth:`~PolyhedralComplex.disjoint_union` | Return the disjoint union of this polyhedral complex with another one.
    :meth:`~PolyhedralComplex.join` | Return the join of this polyhedral complex with another one.

**Miscellaneous**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~PolyhedralComplex.plot` | Return a Graphic object showing the plot of polyhedral complex.
    :meth:`~PolyhedralComplex.graph` | Return a directed graph corresponding to the 1-skeleton of this polyhedral complex, given that it is bounded.
    :meth:`~PolyhedralComplex.union_as_polyhedron` | Return a ``Polyhedron`` which is the union of cells in this polyhedral complex, given that it is convex.

Classes and functions
---------------------
"""

# ****************************************************************************
#       Copyright (C) 2021 Yuan Zhou <yuan.zhou@uky.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from sage.homology.cell_complex import GenericCellComplex
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.geometry.polyhedron.base import is_Polyhedron
from sage.modules.free_module_element import vector
from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.graphs.graph import Graph


class PolyhedralComplex(GenericCellComplex):
    r"""
    Define a PolyhedralComplex.

    INPUT:

    - ``maximal_cells`` -- a list, a tuple, or a dictionary (indexed by
      dimension) of cells of the Complex. Each cell is of class
      :class:`Polyhedron` of the same ambient dimension. To set up a
      :class:PolyhedralComplex, it is sufficient to provide the maximal
      faces. Use keyword argument partial=``True`` to set up a partial
      polyhedral complex, which is a subset of the faces (viewed as
      relatively open) of a polyhedral complex that is not necessarily
      closed under taking intersection.

    - ``maximality_check`` -- boolean; default ``True``;
      if it is ``True``, then the constructor checks that each given
      maximal cell is indeed maximal, and ignores those that are not.


    - ``face_to_face_check`` -- boolean; default ``False``;
      if it is ``True``, then the constructor checks whether the cells
      are face-to-face, and it raises a ValueError if they are not.

    - ``is_mutable`` and ``is_immutable`` -- boolean; default ``True`` and
      ``False``, respectively; Set ``is_mutable=False`` or ``is_immutable=True``
      to make this polyhedral complex immutable.

    EXAMPLES::

        sage: pc = PolyhedralComplex([
        ....: Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1/7, 2/7)]),
        ....: Polyhedron(vertices=[(1/7, 2/7), (0, 0), (0, 1/4)])])
        sage: [p.Vrepresentation() for p in pc.cells_sorted()]
        [(A vertex at (0, 0), A vertex at (0, 1/4), A vertex at (1/7, 2/7)),
         (A vertex at (0, 0), A vertex at (1/3, 1/3), A vertex at (1/7, 2/7)),
         (A vertex at (0, 0), A vertex at (0, 1/4)),
         (A vertex at (0, 0), A vertex at (1/7, 2/7)),
         (A vertex at (0, 0), A vertex at (1/3, 1/3)),
         (A vertex at (0, 1/4), A vertex at (1/7, 2/7)),
         (A vertex at (1/3, 1/3), A vertex at (1/7, 2/7)),
         (A vertex at (0, 0),),
         (A vertex at (0, 1/4),),
         (A vertex at (1/7, 2/7),),
         (A vertex at (1/3, 1/3),)]
        sage: pc.plot()                # not tested
        sage: pc.is_pure()
        True
        sage: pc.is_full_dimensional()
        True
        sage: pc.is_compact()
        True
        sage: pc.boundary_subcomplex()
        Polyhedral complex with 4 maximal cells
        sage: pc.is_convex()
        True
        sage: pc.union_as_polyhedron().Hrepresentation()
        (An inequality (1, -4) x + 1 >= 0,
         An inequality (-1, 1) x + 0 >= 0,
         An inequality (1, 0) x + 0 >= 0)
        sage: pc.face_poset()
        Finite poset containing 11 elements
        sage: pc.is_connected()
        True
        sage: pc.connected_component() == pc
        True

    TESTS:

    Check that non-maximal cells are ignored if ``maximality_check=True``::

        sage: pc = PolyhedralComplex([
        ....:        Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
        ....:        Polyhedron(vertices=[(1, 2), (0, 0)]) ])
        sage: pc.maximal_cells()
        {2: {A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices}}

    Check that non face-to-face can be detected::

        sage: PolyhedralComplex([
        ....:        Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
        ....:        Polyhedron(vertices=[(2, 2), (0, 0)]) ],
        ....:        face_to_face_check=True)
        Traceback (most recent call last):
        ...
        ValueError: The given cells are not face-to-face

    Check that all the cells must have the same ambient dimension::

        sage: PolyhedralComplex([
        ....:        Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
        ....:        Polyhedron(vertices=[[2], [0]]) ])
        Traceback (most recent call last):
        ...
        ValueError: The given cells are not polyhedrain the same ambient space.
    """
    def __init__(self, maximal_cells=None, maximality_check=True,
                 face_to_face_check=False, is_mutable=True, is_immutable=False):
        r"""
        Define a PolyhedralComplex.

        See ``PolyhedralComplex`` for more information.

        EXAMPLES::

            sage: PolyhedralComplex([Polyhedron(vertices=[(1, 1), (0, 0)])])
            Polyhedral complex with 1 maximal cells
        """
        def cells_list_to_cells_dict(cells_list):
            cells_dict = {}
            for cell in cells_list:
                d = cell.dimension()
                if d in cells_dict:
                    cells_dict[d].add(cell)
                else:
                    cells_dict[d] = set([cell])
            return cells_dict

        if maximal_cells is None:
            cells_dict = {}
        elif isinstance(maximal_cells, (list, tuple)):
            cells_dict = cells_list_to_cells_dict(maximal_cells)
        elif isinstance(maximal_cells, dict):
            cells_dict = {k: set(l) for (k, l) in maximal_cells.items()}
        else:
            raise ValueError
        if not cells_dict:
            self._dim = -1
            self._ambient_dim = -1
        else:
            self._dim = max(cells_dict.keys())
            self._ambient_dim = next(iter(cells_dict[self._dim])).ambient_dim()
        self._maximal_cells = cells_dict
        if not all((is_Polyhedron(cell) and
                   cell.ambient_dim() == self._ambient_dim)
                   for cell in self.maximal_cell_iterator()):
            raise ValueError("The given cells are not polyhedra" +
                             "in the same ambient space.")
        # initialize the attributes
        self._is_convex = None
        self._polyhedron = None
        self._maximal_cells_sorted = None    # needed for hash
        self._cells = None
        self._face_poset = None

        if maximality_check:
            cells = self.cells()    # compute self._cells and self._face_poset
            self._maximal_cells = cells_list_to_cells_dict(
                                      self._face_poset.maximal_elements())
        if face_to_face_check:
            poset = self.face_poset()
            maximal_cells = poset.maximal_elements()    # a list
            for i in range(len(maximal_cells)):
                p = maximal_cells[i]
                for j in range(i, len(maximal_cells)):
                    q = maximal_cells[j]
                    r = p.intersection(q)
                    if not (r.is_empty() or (r in poset) and
                            poset.is_gequal(p, r) and poset.is_gequal(q, r)):
                        raise ValueError("The given cells are not face-to-face")
        # For now, a polyhedral complex is immutable
        # TODO: is_mutable and is_immutable parameters and set_immutable method.
        self._is_immutable = True

    def cells(self, subcomplex=None):
        """
        The cells of this polyhedral complex, in the form of a dictionary:
        the keys are integers, representing dimension, and the value
        associated to an integer `d` is the set of `d`-cells.  If the
        optional argument ``subcomplex`` is present, then return only
        the cells which are *not* in the subcomplex.

        :param subcomplex: a subcomplex of this polyhedral complex.  Return
           the cells which are not in this subcomplex.
        :type subcomplex: optional, default None

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: list(pc.cells().keys())
            [2, 1, 0]
        """
        if subcomplex is not None:
            raise NotImplementedError
        if self._cells is not None:
            return self._cells
        maximal_cells = self.maximal_cells()
        cells = {}
        covers = {}
        for k in range(self._dim, -1, -1):
            if k in maximal_cells:
                if k not in cells:
                    cells[k] = set([])
                cells[k].update(maximal_cells[k])
            if k in cells:
                for cell in cells[k]:
                    if cell not in covers:
                        covers[cell] = []
                    for facet in cell.facets():
                        p = facet.as_polyhedron()
                        if p not in covers:
                            covers[p] = []
                        covers[p].append(cell)
                        if (k-1) not in cells:
                            cells[k-1] = set([])
                        cells[k-1].add(p)
        from sage.combinat.posets.posets import Poset
        self._face_poset = Poset(covers)
        self._cells = cells
        return self._cells

    def cell_iterator(self, increasing=True):
        """
        An iterator for the cells in this polyhedral complex.

        INPUT:

        - ``increasing`` -- (optional, default ``True``) if ``True``, return
          cells in increasing order of dimension, thus starting with the
          zero-dimensional cells. Otherwise it returns cells in decreasing
          order of dimension.

        .. NOTE::

            Among the cells of a fixed dimension, there is no sorting.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: len(list(pc.cell_iterator()))
            11
        """
        cells = self.cells()
        dim_index = range(0, self.dimension() + 1)
        if not increasing:
            dim_index = reversed(dim_index)
        for d in dim_index:
            if d in cells:
                for c in cells[d]:
                    yield c

    def _n_cells_sorted(self, n, subcomplex=None):
        """
        Sorted list of cells of dimension ``n`` of this polyhedral complex.
        If the optional argument ``subcomplex`` is present, then
        return the ``n``-dimensional cells which are *not* in the
        subcomplex.

        :param n: the dimension
        :type n: non-negative integer
        :param subcomplex: a subcomplex of this cell complex. Return
           the cells which are not in this subcomplex.
        :type subcomplex: optional, default ``None``

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: [p.Vrepresentation() for p in pc._n_cells_sorted(1)]
            [(A vertex at (0, 0), A vertex at (0, 2)),
             (A vertex at (0, 0), A vertex at (1, 1)),
             (A vertex at (0, 0), A vertex at (1, 2)),
             (A vertex at (0, 2), A vertex at (1, 2)),
             (A vertex at (1, 1), A vertex at (1, 2))]
            sage: pc._n_cells_sorted(3)
            []
        """
        n_cells = self.n_cells(n, subcomplex)
        return sorted(n_cells,
                      key=lambda p: (p.vertices(), p.rays(), p.lines()))

    def cells_sorted(self, subcomplex=None):
        """
        The sorted list of the cells of this polyhedral complex
        in non-increasing dimensions.
        If the optional argument ``subcomplex`` is present, then return only
        the cells which are *not* in the subcomplex.

        :param subcomplex: a subcomplex of this polyhedral complex.  Return
           the cells which are not in this subcomplex.
        :type subcomplex: optional, default None

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: len(pc.cells_sorted())
            11
            sage: pc.cells_sorted()[0].Vrepresentation()
            (A vertex at (0, 0), A vertex at (0, 2), A vertex at (1, 2))
        """
        cells = []
        for n in range(self._dim, -1, -1):
            cells += self._n_cells_sorted(n, subcomplex)
        return cells

    def maximal_cells(self):
        """
        The maximal cells of this polyhedral complex, in the form of a
        dictionary: the keys are integers, representing dimension, and the
        value associated to an integer `d` is the set of `d`-maximal cells.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: len(pc.maximal_cells()[2])
            2
            sage: 1 in pc.maximal_cells()
            False

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:              maximality_check=False)
            sage: len(pc_invalid.maximal_cells()[1])
            1
        """
        return self._maximal_cells

    def maximal_cell_iterator(self, increasing=False):
        """
        An iterator for the maximal cells in this polyhedral complex.

        INPUT:

        - ``increasing`` -- (optional, default ``False``) if ``True``, return
          maximal cells in increasing order of dimension.
          Otherwise it returns cells in decreasing order of dimension.

        .. NOTE::

            Among the cells of a fixed dimension, there is no sorting.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: len(list(pc.maximal_cell_iterator()))
            2

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:              maximality_check=False)
            sage: len(list(pc_invalid.maximal_cell_iterator()))
            3
        """
        maximal_cells = self.maximal_cells()
        dim_index = range(-1, self.dimension() + 1)
        if not increasing:
            dim_index = reversed(dim_index)
        for d in dim_index:
            if d in maximal_cells:
                for c in maximal_cells[d]:
                    yield c

    def n_maximal_cells(self, n):
        """
        List of maximal cells of dimension ``n`` of this polyhedral complex.

        :param n: the dimension
        :type n: non-negative integer

        .. NOTE::

            The resulting list need not be sorted. If you want a sorted
            list of `n`-cells, use :meth:`_n_maximal_cells_sorted`.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: len(pc.n_maximal_cells(2))
            2
            sage: len(pc.n_maximal_cells(1))
            0

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:              maximality_check=False)
            sage: len(pc_invalid.n_maximal_cells(1))
            1
        """
        if n in self.maximal_cells():
            return list(self.maximal_cells()[n])
        else:
            return []

    def _n_maximal_cells_sorted(self, n):
        """
        Sorted list of maximal cells of dimension ``n`` of this polyhedral
        complex.

        :param n: the dimension
        :type n: non-negative integer

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: pc._n_maximal_cells_sorted(2)[0].vertices_list()
            [[0, 0], [0, 2], [1, 2]]
        """
        n_maximal_cells = self.n_maximal_cells(n)
        return sorted(n_maximal_cells,
                      key=lambda p: (p.vertices(), p.rays(), p.lines()))

    def maximal_cells_sorted(self):
        """
        Return the sorted list of the maximal cells of this polyhedral complex
        by non-increasing dimensions.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: [p.vertices_list() for p in pc.maximal_cells_sorted()]
            [[[0, 0], [0, 2], [1, 2]], [[0, 0], [1, 1], [1, 2]]]
        """
        if self._maximal_cells_sorted is None:
            maximal_cells = []
            for n in range(self._dim, -1, -1):
                maximal_cells += self._n_maximal_cells_sorted(n)
            self._maximal_cells_sorted = maximal_cells
        return self._maximal_cells_sorted

    def has_maximal_cell(self, c):
        """
        Return whether the given cell ``c`` is a maximal cell of ``self``.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: pc.has_maximal_cell(p1)
            True
            sage: pc.has_maximal_cell(p3)
            False

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:              maximality_check=False)
            sage: pc_invalid.has_maximal_cell(p3)
            True
        """
        d = c.dimension()
        return (c in self.n_maximal_cells(d))

    def has_cell(self, c):
        """
        Return whether the given cell ``c`` is a cell of ``self``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2])
            sage: pc.has_cell(p3)
            True
            sage: pc.has_cell(Polyhedron(vertices=[(0, 0)]))
            True
        """
        d = c.dimension()
        return (c in self.n_cells(d))

    def dimension(self):
        """
        The dimension of this cell complex: the maximum
        dimension of its cells.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....:        Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....:        Polyhedron(vertices=[(1, 2), (0, 2)]) ])
            sage: pc.dimension()
            2
            sage: empty_pc = PolyhedralComplex([])
            sage: empty_pc.dimension()
            -1
        """
        return self._dim

    def ambient_dimension(self):
        """
        The ambient dimension of this cell complex: the ambient
        dimension of each of its cells.

        EXAMPLES::

            sage: pc = PolyhedralComplex([Polyhedron(vertices=[(1, 2, 3)])])
            sage: pc.ambient_dimension()
            3
            sage: empty_pc = PolyhedralComplex([])
            sage: empty_pc.ambient_dimension()
            -1
        """
        return self._ambient_dim

    def plot(self, **kwds):
        """
        Return a plot of the polyhedral complex, if it is of dim at most 2.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2])
            sage: pc.plot()
            Graphics object consisting of 10 graphics primitives
        """
        if self.dimension() > 2:
            raise ValueError("Cannot plot in high dimension")
        return sum(cell.plot(**kwds) for cell in self.maximal_cell_iterator())

    def is_pure(self):
        """
        Test if this polyhedral complex is pure.

        A polyhedral complex is pure if and only if all of its maximal cells
        have the same dimension.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: pc.is_pure()
            True

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:              maximality_check=False)
            sage: pc_invalid.is_pure()
            False
        """
        return len(self._maximal_cells) == 1

    def is_full_dimensional(self):
        """
        Return whether this polyhedral complex is full-dimensional:
        its dimension is equal to its ambient dimension.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: pc.is_full_dimensional()
            True
            sage: PolyhedralComplex([p3]).is_full_dimensional()
            False
        """
        return self._dim == self._ambient_dim

    def __hash__(self):
        """
        Compute the hash value of ``self`` using its ``maximal_cells_sorted``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])
            sage: pc1 = PolyhedralComplex([p1, p2])
            sage: hash(pc1) == hash(pc1)
            True
            sage: pc2 = PolyhedralComplex([p2, p1])
            sage: hash(pc1) == hash(pc2)
            True
        """
        if not self._is_immutable:
            raise ValueError("This polyhedral complex must be immutable" +
                             "Call set_immutable().")
        return hash(tuple(self.maximal_cells_sorted()))

    def __eq__(self, right):
        """
        Two polyhedral complexes are equal iff their maximal cells are equal.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])
            sage: pc1 = PolyhedralComplex([p1, p2])
            sage: pc1 == pc1
            True
            sage: pc2 = PolyhedralComplex([p2, p1])
            sage: pc1 == pc2
            True
        """
        return isinstance(right, PolyhedralComplex) and (
               self.maximal_cells_sorted() == right.maximal_cells_sorted())

    def __ne__(self, right):
        """
        Return ``True`` if ``self`` and ``right`` are not equal.

        EXAMPLES::

            sage: pc1 = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)])])
            sage: pc2 = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])])
            sage: pc1 != pc2
            True
        """
        return not self.__eq__(right)

    def __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)])
            sage: p2 = copy(p1)
            sage: p1 == p2
            True
        """
        return PolyhedralComplex(self._maximal_cells, maximality_check=False)

    def _an_element_(self):
        """
        Return a (maximal) cell of this complex.

        EXAMPLES::

            sage: PolyhedralComplex()._an_element_()
            ()
            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])])
            sage: pc._an_element_().vertices_list()
            [[0, 0], [0, 1/2], [1, 2]]
        """
        try:
            return next(self.maximal_cell_iterator(increasing=False))
        except StopIteration:
            return ()

    def __contains__(self, x):
        """
        True if ``x`` is a polyhedron which is contained in this complex.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])
            sage: pc = PolyhedralComplex([p1, p2])
            sage: (p1 in pc) and (p2 in pc)
            True
            sage: Polyhedron(vertices=[(1, 2), (0, 0)]) in pc
            True
            sage: Polyhedron(vertices=[(1, 1), (0, 0)]) in pc
            False
            sage: Polyhedron(vertices=[(0, 0)]) in pc
            True
            sage: (0, 0) in pc  # not a polyhedron
            False
        """
        if not is_Polyhedron(x):
            return False
        dim = x.dimension()
        return dim in self.cells() and x in self.cells()[dim]

    def __call__(self, x):
        """
        If ``x`` is a polyhedron in this complex, return it.
        Otherwise, raise a ``ValueError``.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])])
            sage: pc(Polyhedron(vertices=[(1, 2), (0, 0)]))
            A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: pc(Polyhedron(vertices=[(1, 1)]))
            Traceback (most recent call last):
            ...
            ValueError: the polyhedron is not in this complex
        """
        if x not in self:
            raise ValueError('the polyhedron is not in this complex')
        return x

    def face_poset(self):
        r"""
        The face poset of this polyhedral complex, the poset of
        nonempty cells, ordered by inclusion.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])])
            sage: poset = pc.face_poset()
            sage: poset
            Finite poset containing 11 elements
            sage: d = {i:i.vertices_matrix() for i in poset}
            sage: poset.plot(element_labels=d)                    # not tested

        TESTS on nonbounded polyhedral complex::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)]),
            ....: Polyhedron(vertices=[(-1/2, -1/2)], lines=[(1, -1)]),
            ....: Polyhedron(rays=[(1, 0)])])
            sage: poset = pc.face_poset()
            sage: poset
            Finite poset containing 13 elements
            sage: d = {i:''.join([str(v)+'\n'
            ....:        for v in i.Vrepresentation()]) for i in poset}
            sage: poset.show(element_labels=d, figsize=15)        # not tested
            sage: pc = PolyhedralComplex([
            ....: Polyhedron(rays=[(1,0),(0,1)]),
            ....: Polyhedron(rays=[(-1,0),(0,1)]),
            ....: Polyhedron(rays=[(-1,0),(0,-1)]),
            ....: Polyhedron(rays=[(1,0),(0,-1)])])
            sage: pc.face_poset()
            Finite poset containing 9 elements
        """
        if self._face_poset is None:
            cells = self.cells()    # poset is obtained and cached in cells()
        return self._face_poset

    def is_subcomplex(self, other):
        r"""
        Return True if ``self`` is a subcomplex of ``other``.

        :param other: a polyhedral complex

        Each maximal cell of ``self`` must be a cell of ``other``
        for this to be True.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])
            sage: p3 = Polyhedron(vertices=[(0, 0), (1, 0)])
            sage: pc = PolyhedralComplex([p1, Polyhedron(vertices=[(1, 0)])])
            sage: pc.is_subcomplex(PolyhedralComplex([p1, p2, p3]))
            True
            sage: pc.is_subcomplex(PolyhedralComplex([p1, p2]))
            False
        """
        other_cells = other.cells()
        for (d, stratum) in self.maximal_cells().items():
            if not stratum.issubset(other_cells.get(d, set([]))):
                return False
        return True

    def is_compact(self):
        """
        Test for boundedness of the polyhedral complex

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)])
            sage: p2 = Polyhedron(rays=[(1, 0)])
            sage: PolyhedralComplex([p1]).is_compact()
            True
            sage: PolyhedralComplex([p1, p2]).is_compact()
            False
        """
        return all(p.is_compact() for p in self.maximal_cell_iterator())

    def graph(self):
        """
        The 1-skeleton of this polyhedral complex, as a graph.
        The vertices of the graph are of type ``vector``.
        Raise NotImplementedError if the polyhedral complex is unbounded.

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: g = pc.graph(); g
            Graph on 4 vertices
            sage: g.vertices()
            [(0, 0), (0, 2), (1, 1), (1, 2)]
            sage: g.edges(labels=False)
            [((0, 0), (0, 2)), ((0, 0), (1, 1)), ((0, 0), (1, 2)), ((0, 2), (1, 2)), ((1, 1), (1, 2))]
            sage: PolyhedralComplex([Polyhedron(rays=[(1,1)])]).graph()
            Traceback (most recent call last):
            ...
            NotImplementedError: The polyhedral complex is unbounded.

        Wrong answer due to ``maximality_check=False``::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: PolyhedralComplex([p1, p2]).is_pure()
            True
            sage: PolyhedralComplex([p2, p3], maximality_check=True).is_pure()
            True
            sage: PolyhedralComplex([p2, p3], maximality_check=False).is_pure()
            False
        """
        if not self.is_compact():
            raise NotImplementedError("The polyhedral complex is unbounded.")
        edges = self.n_cells(1)
        d = {}
        for e in edges:
            v, max_e = sorted(e.vertices_matrix().columns())
            if v in d:
                d[v].append(max_e)
            else:
                d[v] = [max_e]
        for v in self.n_maximal_cells(0):
            d[v] = []
        return Graph(d)

    def is_connected(self):
        """
        True if this cell complex is connected.

        EXAMPLES::

            sage: pc1 = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: pc1.is_connected()
            True
            sage: pc2 = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(0, 2)])])
            sage: pc2.is_connected()
            False
            sage: pc3 = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 1/2)]),
            ....: Polyhedron(vertices=[(-1/2, -1/2)], lines=[(1, -1)]),
            ....: Polyhedron(rays=[(1, 0)])])
            sage: pc3.is_connected()
            False
            sage: pc4 = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1/3, 1/3), (0, 0), (1, 2)]),
            ....: Polyhedron(rays=[(1, 0)])])
            sage: pc4.is_connected()
            True
        """
        if self.is_compact():
            return self.graph().is_connected()    # faster than using poset?
        else:
            return self.face_poset().is_connected()

    def connected_component(self, cell=None):
        """
        Return the connected component of this polyhedral complex
        containing ``cell``. If ``cell`` is omitted, then return
        the connected component containing the self._an_element_.
        (If the polyhedral complex is empty or if it does not contain the
        given cell, raise an error.)

        EXAMPLES::

            sage: t1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: t2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: v1 = Polyhedron(vertices=[(1, 1)])
            sage: v2 = Polyhedron(vertices=[(0, 2)])
            sage: v3 = Polyhedron(vertices=[(-1, 0)])
            sage: o =  Polyhedron(vertices=[(0, 0)])
            sage: r = Polyhedron(rays=[(1, 0)])
            sage: l = Polyhedron(vertices=[(-1, 0)], lines=[(1, -1)])
            sage: pc1 = PolyhedralComplex([t1, t2])
            sage: pc1.connected_component() == pc1
            True
            sage: pc1.connected_component(v1) == pc1
            True
            sage: pc2 = PolyhedralComplex([t1, v2])
            sage: pc2.connected_component(t1) == PolyhedralComplex([t1])
            True
            sage: pc2.connected_component(o) == PolyhedralComplex([t1])
            True
            sage: pc2.connected_component(v3)
            Traceback (most recent call last):
            ...
            ValueError: the polyhedral complex does not contain the given cell
            sage: pc2.connected_component(r)
            Traceback (most recent call last):
            ...
            ValueError: the polyhedral complex does not contain the given cell
            sage: pc3 = PolyhedralComplex([t1, t2, r])
            sage: pc3.connected_component(v2) == pc3
            True
            sage: pc4 = PolyhedralComplex([t1, t2, r, l])
            sage: pc4.connected_component(o) == pc3
            True
            sage: pc4.connected_component(v3)
            Traceback (most recent call last):
            ...
            ValueError: the polyhedral complex does not contain the given cell
            sage: pc5 = PolyhedralComplex([t1, t2, r, l, v3])
            sage: pc5.connected_component(v3) == PolyhedralComplex([v3])
            True
            sage: PolyhedralComplex([]).connected_component()
            Traceback (most recent call last):
            ...
            ValueError: the empty polyhedral complex has no connected components
        """
        if self.dimension() == -1:
            raise ValueError(
                "the empty polyhedral complex has no connected components")
        if cell is None:
            cell = self._an_element_()
        if self.is_compact():    # use graph (faster than poset?)
            if not cell.is_compact():
                raise ValueError(
                    "the polyhedral complex does not contain the given cell")
            v = cell.vertices_matrix().columns()[0]
            g = self.graph()
            if v not in g:
                raise ValueError(
                    "the polyhedral complex does not contain the given cell")
            vertices = g.connected_component_containing_vertex(v)
            facets = [f for f in self.maximal_cell_iterator()
                      if any(vf in f.vertices_matrix().columns()
                             for vf in vertices)]
        else:    # use face_poset
            g = self.face_poset().hasse_diagram()
            if cell not in g:
                raise ValueError(
                    "the polyhedral complex does not contain the given cell")
            faces = g.connected_component_containing_vertex(cell)
            facets = [f for f in self.maximal_cell_iterator()
                      if f in faces]
        return PolyhedralComplex(facets, maximality_check=False)

    def connected_components(self):
        """
        Return the connected components of this polyhedral complex,
        as list of (sub-)PolyhedralComplexes.

        EXAMPLES::

            sage: t1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: t2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: v1 = Polyhedron(vertices=[(1, 1)])
            sage: v2 = Polyhedron(vertices=[(0, 2)])
            sage: v3 = Polyhedron(vertices=[(-1, 0)])
            sage: o =  Polyhedron(vertices=[(0, 0)])
            sage: r = Polyhedron(rays=[(1, 0)])
            sage: l = Polyhedron(vertices=[(-1, 0)], lines=[(1, -1)])
            sage: pc1 = PolyhedralComplex([t1, t2])
            sage: len(pc1.connected_components())
            1
            sage: pc2 = PolyhedralComplex([t1, v2])
            sage: len(pc2.connected_components())
            2
            sage: pc3 = PolyhedralComplex([t1, t2, r])
            sage: len(pc3.connected_components())
            1
            sage: pc4 = PolyhedralComplex([t1, t2, r, l])
            sage: len(pc4.connected_components())
            2
            sage: pc5 = PolyhedralComplex([t1, t2, r, l, v3])
            sage: len(pc5.connected_components())
            3
            sage: PolyhedralComplex([]).connected_components()
            Traceback (most recent call last):
            ...
            ValueError: the empty polyhedral complex has no connected components
        """
        if self.dimension() == -1:
            raise ValueError(
                "the empty polyhedral complex has no connected components")
        if self.is_compact():    # use graph (faster than poset)?
            g = self.graph()
            lists_of_vertices = g.connected_components(sort=False)
            lists_of_facets = [[f for f in self.maximal_cell_iterator()
                                if any(vf in f.vertices_matrix().columns()
                                       for vf in vertices)]
                               for vertices in lists_of_vertices]
        else:    # use face_poset
            g = self.face_poset().hasse_diagram()
            lists_of_faces = g.connected_components(sort=False)
            lists_of_facets = [
                [f for f in self.maximal_cell_iterator() if f in faces]
                for faces in lists_of_faces]
        results = [PolyhedralComplex(facets, maximality_check=False)
                   for facets in lists_of_facets]
        return results

    def n_skeleton(self, n):
        """
        The `n`-skeleton of this polyhedral complex.

        The `n`-skeleton of a polyhedral complex is obtained by discarding
        all of the cells in dimensions larger than `n`.

        :param n: non-negative integer

        .. SEEALSO::

            :meth:`stratify`

        EXAMPLES::

            sage: pc = PolyhedralComplex([
            ....: Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)]),
            ....: Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])])
            sage: pc.n_skeleton(2)
            Polyhedral complex with 2 maximal cells
            sage: pc.n_skeleton(1)
            Polyhedral complex with 5 maximal cells
            sage: pc.n_skeleton(0)
            Polyhedral complex with 4 maximal cells
        """
        if n >= self.dimension():
            return copy(self)
        facets = [f for f in self.maximal_cell_iterator() if f.dimension() < n]
        facets.extend(self.n_cells(n))
        return PolyhedralComplex(facets, maximality_check=False)

    def stratify(self, n):
        """
        Return the pure sub-polyhedral complex which is constructed from the
        n-dimensional maximal cells of this polyhedral complex.

        .. SEEALSO::

            :meth:`n_skeleton`

        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: pc = PolyhedralComplex([p1, p2, p3])
            sage: pc.stratify(2) == pc
            True
            sage: pc.stratify(1)
            Polyhedral complex with 0 maximal cells

        Wrong answer due to ``maximality_check=False``::

            sage: pc_invalid = PolyhedralComplex([p1, p2, p3],
            ....:              maximality_check=False)
            sage: pc_invalid.stratify(1)
            Polyhedral complex with 1 maximal cells
        """
        n_faces = self.n_maximal_cells(n)
        return PolyhedralComplex(n_faces, maximality_check=False)

    def boundary_subcomplex(self):
        """
        Return the sub-polyhedral complex that is the boundary of ``self``.

        A point `P` is on the boundary of a set `S` if `P` is in the
        closure of `S` but not in the interoir of `S`.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: bd = PolyhedralComplex([p1, p2]).boundary_subcomplex()
            sage: len(bd.n_maximal_cells(2))
            0
            sage: len(bd.n_maximal_cells(1))
            4
            sage: pt = PolyhedralComplex([p3])
            sage: pt.boundary_subcomplex() == pt
            True

        Test on polyhedral complex which is not pure::

            sage: pc_non_pure = PolyhedralComplex([p1, p3])
            sage: pc_non_pure.boundary_subcomplex() == pc_non_pure.n_skeleton(1)
            True

        Test with ``maximality_check == False``::

            sage: pc_invalid = PolyhedralComplex([p2, p3],
            ....:              maximality_check=False)
            sage: pc_invalid.boundary_subcomplex() == pc_invalid.n_skeleton(1)
            True

        Test unbounded cases::

            sage: pc1 = PolyhedralComplex([
            ....: Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,0], [0,1]])])
            sage: pc1.boundary_subcomplex() == pc1.n_skeleton(1)
            True
            sage: pc1b = PolyhedralComplex([Polyhedron(
            ....: vertices=[[1,0,0], [0,1,0]], rays=[[1,0,0],[0,1,0]])])
            sage: pc1b.boundary_subcomplex() == pc1b
            True
            sage: pc2 = PolyhedralComplex([
            ....:        Polyhedron(vertices=[[-1,0], [1,0]], lines=[[0,1]])])
            sage: pc2.boundary_subcomplex() == pc2.n_skeleton(1)
            True
            sage: pc3 = PolyhedralComplex([
            ....: Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,0], [0,1]]),
            ....: Polyhedron(vertices=[[1,0], [0,-1]], rays=[[1,0], [0,-1]])])
            sage: pc3.boundary_subcomplex() == pc3.n_skeleton(1)
            False
        """
        if self.is_full_dimensional():
            return PolyhedralComplex(self.relative_boundary_cells())
        else:
            return copy(self)

    def relative_boundary_cells(self):
        r"""
        Return the maximal cells of the relative-boundary sub-complex.

        A point `P` is in the relative boundary of a set `S` if `P` is in the
        closure of `S` but not in the relative interoir of `S`.


        .. WARNING::

            This may give the wrong answer if the polyhedral complex
            was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(1, 2), (0, 2)])
            sage: p4 = Polyhedron(vertices=[(2, 2)])
            sage: pc = PolyhedralComplex([p1, p2])
            sage: rbd_cells = pc.relative_boundary_cells()
            sage: len(rbd_cells)
            4
            sage: all(p.dimension() == 1 for p in rbd_cells)
            True
            sage: pc_lower_dim = PolyhedralComplex([p3])
            sage: [p.vertices() for p in pc_lower_dim.relative_boundary_cells()]
            [(A vertex at (0, 2),), (A vertex at (1, 2),)]

        Test on polyhedral complex which is not pure::

            sage: pc_non_pure = PolyhedralComplex([p1, p3, p4])
            sage: set(pc_non_pure.relative_boundary_cells()) == set(
            ....: [f.as_polyhedron() for f in p1.faces(1)] + [p3, p4])
            True

        Test with ``maximality_check == False``::

            sage: pc_invalid = PolyhedralComplex([p2, p3],
            ....:              maximality_check=False)
            sage: set(pc_invalid.relative_boundary_cells()) == set(
            ....:                    [f.as_polyhedron() for f in p2.faces(1)])
            True

        Test unbounded case::

            sage: pc3 = PolyhedralComplex([
            ....: Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,0], [0,1]]),
            ....: Polyhedron(vertices=[[1,0], [0,-1]], rays=[[1,0], [0,-1]])])
            sage: len(pc3.relative_boundary_cells())
            4
        """
        d = self.dimension()
        poset = self.face_poset()
        faces = self.n_cells(d - 1)
        ans = [face for face in faces if len(poset.upper_covers(face)) == 1]
        if not self.is_pure():
            ans += [p for p in poset.maximal_elements() if p.dimension() < d]
        return ans

    def is_convex(self):
        """
        Return whether the set of points in ``self`` is a convex set.

        When ``self`` is convex, the union of its cells is a Polyhedron.

        .. SEEALSO::

            :meth:`union_as_polyhedron`

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(0, 0), (1, 1), (2, 0)])
            sage: p4 = Polyhedron(vertices=[(2, 2)])
            sage: PolyhedralComplex([p1, p2]).is_convex()
            True
            sage: PolyhedralComplex([p1, p3]).is_convex()
            False
            sage: PolyhedralComplex([p1, p4]).is_convex()
            False

        Test unbounded cases::

            sage: pc1 = PolyhedralComplex([
            ....: Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,0], [0,1]])])
            sage: pc1.is_convex()
            True
            sage: pc2 = PolyhedralComplex([
            ....:        Polyhedron(vertices=[[-1,0], [1,0]], lines=[[0,1]])])
            sage: pc2.is_convex()
            True
            sage: pc3 = PolyhedralComplex([
            ....: Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,0], [0,1]]),
            ....: Polyhedron(vertices=[[1,0], [0,-1]], rays=[[1,0], [0,-1]])])
            sage: pc3.is_convex()
            False
            sage: pc4 = PolyhedralComplex([Polyhedron(rays=[[1,0], [-1,1]]),
            ....:                         Polyhedron(rays=[[1,0], [-1,-1]])])
            sage: pc4.is_convex()
            False

        The whole 3d space minus the first orthant is not convex::

            sage: pc5 = PolyhedralComplex([
            ....:          Polyhedron(rays=[[1,0,0], [0,1,0], [0,0,-1]]),
            ....:          Polyhedron(rays=[[1,0,0], [0,-1,0], [0,0,-1]]),
            ....:          Polyhedron(rays=[[1,0,0], [0,-1,0], [0,0,1]]),
            ....:          Polyhedron(rays=[[-1,0,0], [0,-1,0], [0,0,-1]]),
            ....:          Polyhedron(rays=[[-1,0,0], [0,-1,0], [0,0,1]]),
            ....:          Polyhedron(rays=[[-1,0,0], [0,1,0], [0,0,-1]]),
            ....:          Polyhedron(rays=[[-1,0,0], [0,1,0], [0,0,1]])])
            sage: pc5.is_convex()
            False

        Test some non-full-dimensional examples::

            sage: l = PolyhedralComplex([Polyhedron(vertices=[(1, 2), (0, 2)])])
            sage: l.is_convex()
            True
            sage: pc1b = PolyhedralComplex([Polyhedron(
            ....: vertices=[[1,0,0], [0,1,0]], rays=[[1,0,0],[0,1,0]])])
            sage: pc1b.is_convex()
            True
            sage: pc4b = PolyhedralComplex([
            ....: Polyhedron(rays=[[1,0,0], [-1,1,0]]),
            ....: Polyhedron(rays=[[1,0,0], [-1,-1,0]])])
            sage: pc4b.is_convex()
            False
        """
        if self._is_convex is not None:
            return self._is_convex
        if not self.is_pure():
            self._is_convex = False
            return False
        d = self.dimension()
        if not self.is_full_dimensional():
            # if max cells must lie in different subspaces, can't be convex.
            from sage.modules.free_module import span
            f = self.n_maximal_cells(d)[0]
            affine_space = span(f.equations_list(), f.base_ring())
            for f in self.n_maximal_cells(d)[1::]:
                if span(f.equations_list(), f.base_ring()) != affine_space:
                    self._is_convex = False
                    return False
        # orient the (relative) boundary halfspaces toward a strict convex
        # combination of the vertices. Then check if all vertices are contained
        # After making sure that the affine hulls of the cells are the same,
        # it does not matter that is not full dimensional.
        boundaries = self.relative_boundary_cells()
        vertices = set([])
        rays = set([])
        lines = set([])
        for cell in boundaries:
            # it suffices to consider only vertices on the boundaries
            # Note that a line (as polyhedron) has vertex too
            for v in cell.vertices_list():
                vv = vector(v)
                vv.set_immutable()
                vertices.add(vv)
        for cell in self.n_maximal_cells(d):
            for r in cell.rays_list():
                rr = vector(r)
                rr.set_immutable()
                rays.add(rr)
            for li in cell.lines_list():
                ll = vector(li)
                ll.set_immutable()
                lines.add(ll)
        center = sum(vertices) / len(vertices)
        for cell in boundaries:
            for equation in cell.equations_list():
                coeff = vector(equation[1::])
                const = equation[0]
                if const + coeff * center == 0:
                    sign = 0
                elif const + coeff * center > 0:
                    sign = 1
                    for v in vertices:
                        if const + coeff * v < 0:
                            self._is_convex = False
                            return False
                elif const + coeff * center < 0:
                    sign = -1
                    for v in vertices:
                        if const + coeff * v > 0:
                            self._is_convex = False
                            return False
                for r in rays:
                    if sign == 0:
                        sign = coeff * r
                    else:
                        if sign * (coeff * r) < 0:
                            self._is_convex = False
                            return False
                # lines are in the affine space of each boundary cell already
        self._is_convex = True
        self._polyhedron = Polyhedron(vertices=vertices, rays=rays, lines=lines)
        return True

    def union_as_polyhedron(self):
        """
        Assuming the polyhedral complex is convex, return it as a Polyhedron.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(0, 0), (1, 1), (2, 0)])
            sage: P = PolyhedralComplex([p1, p2]).union_as_polyhedron()
            sage: P.vertices_list()
            [[0, 0], [0, 2], [1, 1], [1, 2]]
            sage: PolyhedralComplex([p1, p3]).union_as_polyhedron()
            Traceback (most recent call last):
            ...
            ValueError: The polyhedral complex is not convex.
        """
        if not self.is_convex():
            raise ValueError("The polyhedral complex is not convex.")
        return self._polyhedron

    def product(self, right):
        """
        The (Cartesian) product of this polyhedral complex with another one.

        :param right: the other polyhedral complex (the right-hand
           factor)

        :return: the product ``self x right``

        EXAMPLES::

            sage: pc = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])])
            sage: pc_square = pc.product(pc)
            sage: pc_square
            Polyhedral complex with 1 maximal cells
            sage: next(pc_square.maximal_cell_iterator()).vertices()
            (A vertex at (0, 0),
             A vertex at (0, 1),
             A vertex at (1, 0),
             A vertex at (1, 1))
        """
        maximal_cells = [f.product(g) for f in self.maximal_cell_iterator()
                         for g in right.maximal_cell_iterator()]
        return PolyhedralComplex(maximal_cells, maximality_check=False)

    def disjoint_union(self, right):
        """
        The disjoint union of this polyhedral complex with another one.

        :param right: the other polyhedral complex (the right-hand factor)

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: p2 = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: p3 = Polyhedron(vertices=[(0, 0), (1, 1), (2, 0)])
            sage: pc1 = PolyhedralComplex([p1, p3])
            sage: pc2 = PolyhedralComplex([p2])
            sage: pc = pc1.disjoint_union(pc2)
            sage: set(pc.maximal_cell_iterator()) == set([p1, p2, p3])
            True
        """
        maximal_cells = list(self.maximal_cell_iterator()) + list(
                        right.maximal_cell_iterator())
        return PolyhedralComplex(maximal_cells, maximality_check=True,
                                 face_to_face_check=True)

    def join(self, right):
        """
        The join of this polyhedral complex with another one.

        :param right: the other polyhedral complex (the right-hand factor)

        EXAMPLES::

            sage: pc = PolyhedralComplex([Polyhedron(vertices=[[0], [1]])])
            sage: pc_join = pc.join(pc)
            sage: pc_join
            Polyhedral complex with 1 maximal cells
            sage: next(pc_join.maximal_cell_iterator()).vertices()
            (A vertex at (0, 0, 0),
             A vertex at (0, 0, 1),
             A vertex at (0, 1, 1),
             A vertex at (1, 0, 0))
        """
        maximal_cells = [f.join(g) for f in self.maximal_cell_iterator()
                         for g in right.maximal_cell_iterator()]
        return PolyhedralComplex(maximal_cells, maximality_check=False)

    ############################################################
    # abstract methods not implemented in generic cell complexe
    ############################################################

    def wedge(self, right):
        raise NotImplementedError

    ############################################################
    # chain complexes, homology
    ############################################################
    def chain_complex(self, subcomplex=None, augmented=False,
                      verbose=False, check=True, dimensions=None,
                      base_ring=ZZ, cochain=False):
        raise NotImplementedError

    def alexander_whitney(self, cell, dim_left):
        raise NotImplementedError

    ############################################################
    # end of chain complexes, homology
    ############################################################

# TODO: mutable complex: add and update stuff incrementally
# TODO: replace one cell by its triangulation and adapt other cells
# TODO: graph of maximal cells by wall-crossing # use poset.meet instead
# TODO: SimplicialComplex to PolyhedralComplex: geometric realization
# TODO: learn about the boundary stuff of chain complex
# TODO: Polyhedral Arrangement to PolyhedralComplex using #25122
