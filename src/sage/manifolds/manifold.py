r"""
Topological Manifolds

Given a topological field `K` (in most applications, `K = \RR` or
`K = \CC`) and a non-negative integer `n`, a *topological manifold of
dimension* `n` *over K* is a topological space `M` such that

- `M` is a Hausdorff space,
- `M` is second countable,
- every point in `M` has a neighborhood homeomorphic to `K^n`.

Topological manifolds are implemented via the class
:class:`TopologicalManifold`. Open subsets of topological manifolds
are also implemented via :class:`TopologicalManifold`, since they are
topological manifolds by themselves.

In the current setting, topological manifolds are mostly described by
means of charts (see :class:`~sage.manifolds.chart.Chart`).

:class:`TopologicalManifold` serves as a base class for more specific
manifold classes.

The user interface is provided by the generic function
:func:`~sage.manifolds.manifold.Manifold`, with
with the argument ``structure`` set to ``'topological'``.

.. RUBRIC:: Example 1: the 2-sphere as a topological manifold of dimension
  2 over `\RR`

One starts by declaring `S^2` as a 2-dimensional topological manifold::

    sage: M = Manifold(2, 'S^2', structure='topological')
    sage: M
    2-dimensional topological manifold S^2

Since the base topological field has not been specified in the argument list
of ``Manifold``, `\RR` is assumed::

    sage: M.base_field()
    Real Field with 53 bits of precision
    sage: dim(M)
    2

Let us consider the complement of a point, the "North pole" say; this is an
open subset of `S^2`, which we call `U`::

    sage: U = M.open_subset('U'); U
    Open subset U of the 2-dimensional topological manifold S^2

A standard chart on `U` is provided by the stereographic projection from the
North pole to the equatorial plane::

    sage: stereoN.<x,y> = U.chart(); stereoN
    Chart (U, (x, y))

Thanks to the operator ``<x,y>`` on the left-hand side, the coordinates
declared in a chart (here `x` and `y`), are accessible by their names;
they are Sage's symbolic variables::

    sage: y
    y
    sage: type(y)
    <type 'sage.symbolic.expression.Expression'>

The South pole is the point of coordinates `(x, y) = (0, 0)` in the above
chart::

    sage: S = U.point((0,0), chart=stereoN, name='S'); S
    Point S on the 2-dimensional topological manifold S^2

Let us call `V` the open subset that is the complement of the South pole and
let us introduce on it the chart induced by the stereographic projection from
the South pole to the equatorial plane::

    sage: V = M.open_subset('V'); V
    Open subset V of the 2-dimensional topological manifold S^2
    sage: stereoS.<u,v> = V.chart(); stereoS
    Chart (V, (u, v))

The North pole is the point of coordinates `(u, v) = (0, 0)` in this chart::

    sage: N = V.point((0,0), chart=stereoS, name='N'); N
    Point N on the 2-dimensional topological manifold S^2

To fully construct the manifold, we declare that it is the union of `U`
and `V`::

    sage: M.declare_union(U,V)

and we provide the transition map between the charts ``stereoN`` =
`(U, (x, y))` and ``stereoS`` = `(V, (u, v))`, denoting by `W` the
intersection of `U` and `V` (`W` is the subset of `U` defined by
`x^2 + y^2 \neq 0`, as well as the subset of `V` defined by
`u^2 + v^2 \neq 0`)::

    sage: stereoN_to_S = stereoN.transition_map(stereoS, [x/(x^2+y^2), y/(x^2+y^2)],
    ....:                          intersection_name='W', restrictions1= x^2+y^2!=0,
    ....:                                                 restrictions2= u^2+v^2!=0)
    sage: stereoN_to_S
    Change of coordinates from Chart (W, (x, y)) to Chart (W, (u, v))
    sage: stereoN_to_S.display()
    u = x/(x^2 + y^2)
    v = y/(x^2 + y^2)

We give the name ``W`` to the Python variable representing `W = U \cap V`::

    sage: W = U.intersection(V)

The inverse of the transition map is computed by the method
:meth:`sage.manifolds.chart.CoordChange.inverse`::

    sage: stereoN_to_S.inverse()
    Change of coordinates from Chart (W, (u, v)) to Chart (W, (x, y))
    sage: stereoN_to_S.inverse().display()
    x = u/(u^2 + v^2)
    y = v/(u^2 + v^2)

At this stage, we have four open subsets on `S^2`::

    sage: M.list_of_subsets()
    [2-dimensional topological manifold S^2,
     Open subset U of the 2-dimensional topological manifold S^2,
     Open subset V of the 2-dimensional topological manifold S^2,
     Open subset W of the 2-dimensional topological manifold S^2]

`W` is the open subset that is the complement of the two poles::

    sage: N in W or S in W
    False

The North pole lies in `V` and the South pole in `U`::

    sage: N in V, N in U
    (True, False)
    sage: S in U, S in V
    (True, False)

The manifold's (user) atlas contains four charts, two of them
being restrictions of charts to a smaller domain::

    sage: M.atlas()
    [Chart (U, (x, y)), Chart (V, (u, v)),
     Chart (W, (x, y)), Chart (W, (u, v))]

Let us consider the point of coordinates `(1, 2)` in the chart ``stereoN``::

    sage: p = M.point((1,2), chart=stereoN, name='p'); p
    Point p on the 2-dimensional topological manifold S^2
    sage: p.parent()
    2-dimensional topological manifold S^2
    sage: p in W
    True

The coordinates of `p` in the chart ``stereoS`` are computed by letting
the chart act on the point::

    sage: stereoS(p)
    (1/5, 2/5)

Given the definition of `p`, we have of course::

    sage: stereoN(p)
    (1, 2)

Similarly::

    sage: stereoS(N)
    (0, 0)
    sage: stereoN(S)
    (0, 0)


.. RUBRIC:: Example 2: the Riemann sphere as a topological manifold of
  dimension 1 over `\CC`

We declare the Riemann sphere `\CC^*` as a 1-dimensional topological manifold
over `\CC`::

    sage: M = Manifold(1, 'C*', structure='topological', field='complex'); M
    Complex 1-dimensional topological manifold C*

We introduce a first open subset, which is actually
`\CC = \CC^*\setminus\{\infty\}` if we interpret `\CC^*` as the
Alexandroff one-point compactification of `\CC`::

    sage: U = M.open_subset('U')

A natural chart on `U` is then nothing but the identity map of `\CC`, hence
we denote the associated coordinate by `z`::

    sage: Z.<z> = U.chart()

The origin of the complex plane is the point of coordinate `z = 0`::

    sage: O = U.point((0,), chart=Z, name='O'); O
    Point O on the Complex 1-dimensional topological manifold C*

Another open subset of `\CC^*` is `V = \CC^*\setminus\{O\}`::

    sage: V = M.open_subset('V')

We define a chart on `V` such that the point at infinity is the point of
coordinate `0` in this chart::

    sage: W.<w> = V.chart(); W
    Chart (V, (w,))
    sage: inf = M.point((0,), chart=W, name='inf', latex_name=r'\infty')
    sage: inf
    Point inf on the Complex 1-dimensional topological manifold C*

To fully construct the Riemann sphere, we declare that it is the union
of `U` and `V`::

    sage: M.declare_union(U,V)

and we provide the transition map between the two charts as `w = 1 / z`
on `A = U \cap V`::

    sage: Z_to_W = Z.transition_map(W, 1/z, intersection_name='A',
    ....:                           restrictions1= z!=0, restrictions2= w!=0)
    sage: Z_to_W
    Change of coordinates from Chart (A, (z,)) to Chart (A, (w,))
    sage: Z_to_W.display()
    w = 1/z
    sage: Z_to_W.inverse()
    Change of coordinates from Chart (A, (w,)) to Chart (A, (z,))
    sage: Z_to_W.inverse().display()
    z = 1/w

Let consider the complex number `i` as a point of the Riemann sphere::

    sage: i = M((I,), chart=Z, name='i'); i
    Point i on the Complex 1-dimensional topological manifold C*

Its coordinates w.r.t. the charts ``Z`` and ``W`` are::

    sage: Z(i)
    (I,)
    sage: W(i)
    (-I,)

and we have::

    sage: i in U
    True
    sage: i in V
    True

The following subsets and charts have been defined::

    sage: M.list_of_subsets()
    [Open subset A of the Complex 1-dimensional topological manifold C*,
     Complex 1-dimensional topological manifold C*,
     Open subset U of the Complex 1-dimensional topological manifold C*,
     Open subset V of the Complex 1-dimensional topological manifold C*]
    sage: M.atlas()
    [Chart (U, (z,)), Chart (V, (w,)), Chart (A, (z,)), Chart (A, (w,))]


AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2015): structure described via
  :class:`~sage.manifolds.structure.TopologicalStructure` or
  :class:`~sage.manifolds.structure.RealTopologicalStructure`


REFERENCES:

.. [Lee11] J.M. Lee : *Introduction to Topological Manifolds*,
   2nd ed., Springer (New York) (2011).
.. [Lee13] J.M. Lee : *Introduction to Smooth Manifolds*,
   2nd ed., Springer (New York) (2013)
.. [KN63] S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*,
   vol. 1, Interscience Publishers (New York) (1963).
.. [Huybrechts05] D. Huybrechts : *Complex Geometry*,
   Springer (Berlin) (2005).
"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh@umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.fields import Fields
from sage.categories.manifolds import Manifolds
from sage.rings.all import CC
from sage.rings.real_mpfr import RR, RealField_class
from sage.rings.complex_field import ComplexField_class
from sage.misc.prandom import getrandbits
from sage.rings.integer import Integer
from sage.manifolds.subset import ManifoldSubset
from sage.manifolds.structure import TopologicalStructure, \
                                     RealTopologicalStructure

#############################################################################
## Class

class TopologicalManifold(ManifoldSubset):
    r"""
    Topological manifold over a topological field `K`.

    Given a topological field `K` (in most applications, `K = \RR` or
    `K = \CC`) and a non-negative integer `n`, a *topological manifold of
    dimension* `n` *over K* is a topological space `M` such that

    - `M` is a Hausdorff space,
    - `M` is second countable, and
    - every point in `M` has a neighborhood homeomorphic to `K^n`.

    This is a Sage *parent* class, the corresponding *element*
    class being :class:`~sage.manifolds.point.ManifoldPoint`.

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``field`` -- field `K` on which the manifold is
      defined; allowed values are

      - ``'real'`` or an object of type ``RealField`` (e.g., ``RR``) for
        a manifold over `\RR`
      - ``'complex'`` or an object of type ``ComplexField`` (e.g., ``CC``)
        for a manifold over `\CC`
      - an object in the category of topological fields (see
        :class:`~sage.categories.fields.Fields` and
        :class:`~sage.categories.topological_spaces.TopologicalSpaces`)
        for other types of manifolds

    - ``structure`` -- manifold structure (see
      :class:`~sage.manifolds.structure.TopologicalStructure` or
      :class:`~sage.manifolds.structure.RealTopologicalStructure`)
    - ``ambient`` -- (default: ``None``) if not ``None``, must be a
      topological manifold; the created object is then an open subset of
      ``ambient``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the manifold; if none are provided, it is set to ``name``
    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the manifold, e.g., coordinates
      in a chart
    - ``category`` -- (default: ``None``) to specify the category; if
      ``None``, ``Manifolds(field)`` is assumed (see the category
      :class:`~sage.categories.manifolds.Manifolds`)
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.subset.ManifoldSubset`
      would return the previously constructed object corresponding to these
      arguments)

    EXAMPLES:

    A 4-dimensional topological manifold (over `\RR`)::

        sage: M = Manifold(4, 'M', latex_name=r'\mathcal{M}', structure='topological')
        sage: M
        4-dimensional topological manifold M
        sage: latex(M)
        \mathcal{M}
        sage: type(M)
        <class 'sage.manifolds.manifold.TopologicalManifold_with_category'>
        sage: M.base_field()
        Real Field with 53 bits of precision
        sage: dim(M)
        4

    The input parameter ``start_index`` defines the range of indices
    on the manifold::

        sage: M = Manifold(4, 'M', structure='topological')
        sage: list(M.irange())
        [0, 1, 2, 3]
        sage: M = Manifold(4, 'M', structure='topological', start_index=1)
        sage: list(M.irange())
        [1, 2, 3, 4]
        sage: list(Manifold(4, 'M', structure='topological', start_index=-2).irange())
        [-2, -1, 0, 1]

    A complex manifold::

        sage: N = Manifold(3, 'N', structure='topological', field='complex'); N
        Complex 3-dimensional topological manifold N

    A manifold over `\QQ`::

        sage: N = Manifold(6, 'N', structure='topological', field=QQ); N
        6-dimensional topological manifold N over the Rational Field

    A manifold over `\QQ_5`, the field of 5-adic numbers::

        sage: N = Manifold(2, 'N', structure='topological', field=Qp(5)); N
        2-dimensional topological manifold N over the 5-adic Field with capped
         relative precision 20

    A manifold is a Sage *parent* object, in the category of topological
    manifolds over a given topological field (see
    :class:`~sage.categories.manifolds.Manifolds`)::

        sage: isinstance(M, Parent)
        True
        sage: M.category()
        Category of manifolds over Real Field with 53 bits of precision
        sage: from sage.categories.manifolds import Manifolds
        sage: M.category() is Manifolds(RR)
        True
        sage: M.category() is Manifolds(M.base_field())
        True
        sage: M in Manifolds(RR)
        True
        sage: N in Manifolds(Qp(5))
        True

    The corresponding Sage *elements* are points::

        sage: X.<t, x, y, z> = M.chart()
        sage: p = M.an_element(); p
        Point on the 4-dimensional topological manifold M
        sage: p.parent()
        4-dimensional topological manifold M
        sage: M.is_parent_of(p)
        True
        sage: p in M
        True

    The manifold's points are instances of class
    :class:`~sage.manifolds.point.ManifoldPoint`::

        sage: isinstance(p, sage.manifolds.point.ManifoldPoint)
        True

    Since an open subset of a topological manifold `M` is itself a
    topological manifold, open subsets of `M` are instances of the class
    :class:`TopologicalManifold`::

        sage: U = M.open_subset('U'); U
        Open subset U of the 4-dimensional topological manifold M
        sage: isinstance(U, sage.manifolds.manifold.TopologicalManifold)
        True
        sage: U.base_field() == M.base_field()
        True
        sage: dim(U) == dim(M)
        True
        sage: U.category()
        Join of Category of subobjects of sets and Category of manifolds over
         Real Field with 53 bits of precision

    The manifold passes all the tests of the test suite relative to its
    category::

        sage: TestSuite(M).run()

    .. SEEALSO::

        :mod:`sage.manifolds.manifold`
    """
    def __init__(self, n, name, field, structure, ambient=None,
                 latex_name=None, start_index=0, category=None,
                 unique_tag=None):
        r"""
        Construct a topological manifold.

        TESTS::

            sage: M = Manifold(3, 'M', latex_name=r'\mathbb{M}',
            ....:              structure='topological', start_index=1)
            sage: M
            3-dimensional topological manifold M
            sage: latex(M)
            \mathbb{M}
            sage: dim(M)
            3
            sage: X.<x,y,z> = M.chart()
            sage: TestSuite(M).run()

        Tests for open subsets::

            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: TestSuite(U).run()
            sage: U.category() is M.category().Subobjects()
            True

        """
        # Initialization of the attributes _dim, _field, _field_type:
        self._dim = n
        if field == 'real':
            self._field = RR
            self._field_type = 'real'
        elif field == 'complex':
            self._field = CC
            self._field_type = 'complex'
        else:
            if field not in Fields():
                raise TypeError("the argument 'field' must be a field")
            self._field = field
            if isinstance(field, RealField_class):
                self._field_type = 'real'
            elif isinstance(field, ComplexField_class):
                self._field_type = 'complex'
            else:
                self._field_type = 'neither_real_nor_complex'
        # Structure and category:
        self._structure = structure
        if ambient is None:
            ambient = self
            category = Manifolds(self._field).or_subcategory(category)
            category = self._structure.subcategory(category)
        else:
            category = ambient.category().Subobjects()
        # Initialization as a manifold set:
        ManifoldSubset.__init__(self, ambient, name, latex_name=latex_name,
                                category=category)
        self._is_open = True
        self._open_covers.append([self])  # list of open covers of self
        #
        if not isinstance(start_index, (int, Integer)):
            raise TypeError("the starting index must be an integer")
        self._sindex = start_index
        #
        self._atlas = []  # list of charts defined on subsets of self
        self._top_charts = []  # list of charts defined on subsets of self
                        # that are not subcharts of charts on larger subsets
        self._def_chart = None  # default chart
        self._charts_by_coord = {} # dictionary of charts whose domain is self
                                   # (key: string formed by the coordinate
                                   #  symbols separated by a white space)
        self._coord_changes = {} # dictionary of transition maps (key: pair of
                                 # of charts)
        # List of charts that individually cover self, i.e. whose
        # domains are self (if non-empty, self is a coordinate domain):
        self._covering_charts = []

    def _repr_(self):
        r"""
        Return a string representation of the manifold.

        TESTS::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: M._repr_()
            '3-dimensional topological manifold M'
            sage: repr(M)  # indirect doctest
            '3-dimensional topological manifold M'
            sage: M  # indirect doctest
            3-dimensional topological manifold M
            sage: M = Manifold(3, 'M', structure='topological', field='complex')
            sage: M._repr_()
            'Complex 3-dimensional topological manifold M'
            sage: M = Manifold(3, 'M', structure='topological', field=QQ)
            sage: M._repr_()
            '3-dimensional topological manifold M over the Rational Field'

        If the manifold is actually an open subset of a larger manifold, the
        string representation is different::

            sage: U = M.open_subset('U')
            sage: U._repr_()
            'Open subset U of the 3-dimensional topological manifold M
             over the Rational Field'
        """
        if self is self._manifold:
            if self._field_type == 'real':
                return "{}-dimensional {} manifold {}".format(self._dim,
                                                          self._structure.name,
                                                          self._name)
            elif self._field_type == 'complex':
                return "Complex {}-dimensional {} manifold {}".format(self._dim,
                                                           self._structure.name,
                                                           self._name)
            return "{}-dimensional {} manifold {} over the {}".format(self._dim,
                                                           self._structure.name,
                                                           self._name,
                                                           self._field)
        else:
            return "Open subset {} of the {}".format(self._name, self._manifold)

    def _an_element_(self):
        r"""
        Construct some point on the manifold.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M._an_element_(); p
            Point on the 2-dimensional topological manifold M
            sage: p.coord()
            (0, 0)
            sage: U = M.open_subset('U', coord_def={X: y>1}); U
            Open subset U of the 2-dimensional topological manifold M
            sage: p = U._an_element_(); p
            Point on the 2-dimensional topological manifold M
            sage: p in U
            True
            sage: p.coord()
            (0, 2)
            sage: V = U.open_subset('V', coord_def={X.restrict(U): x<-pi})
            sage: p = V._an_element_(); p
            Point on the 2-dimensional topological manifold M
            sage: p in V
            True
            sage: p.coord()
            (-pi - 1, 2)

        """
        from sage.rings.infinity import Infinity
        if self._def_chart is None:
            return self.element_class(self)
        # Attempt to construct a point in the domain of the default chart
        chart = self._def_chart
        if self._field_type == 'real':
            coords = []
            for coord_range in chart._bounds:
                xmin = coord_range[0][0]
                xmax = coord_range[1][0]
                if xmin == -Infinity:
                    if xmax == Infinity:
                        x = 0
                    else:
                        x = xmax - 1
                else:
                    if xmax == Infinity:
                        x = xmin + 1
                    else:
                        x = (xmin + xmax)/2
                coords.append(x)
        else:
            coords = self._dim*[0]
        if not chart.valid_coordinates(*coords):
            # Attempt to construct a point in the domain of other charts
            if self._field_type == 'real':
                for ch in self._atlas:
                    if ch is self._def_chart:
                        continue # since this case has already been attempted
                    coords = []
                    for coord_range in ch._bounds:
                        xmin = coord_range[0][0]
                        xmax = coord_range[1][0]
                        if xmin == -Infinity:
                            if xmax == Infinity:
                                x = 0
                            else:
                                x = xmax - 1
                        else:
                            if xmax == Infinity:
                                x = xmin + 1
                            else:
                                x = (xmin + xmax)/2
                        coords.append(x)
                    if ch.valid_coordinates(*coords):
                        chart = ch
                        break
                else:
                    # A generic element with specific coordinates could not be
                    # automatically generated, due to too complex cooordinate
                    # conditions. An element without any coordinate set is
                    # returned instead:
                    return self.element_class(self)
            else:
                # Case of manifolds over a field different from R
                for ch in self._atlas:
                    if ch is self._def_chart:
                        continue # since this case has already been attempted
                    if ch.valid_coordinates(*coords):
                        chart = ch
                        break
                else:
                    return self.element_class(self)
        # The point is constructed with check_coords=False since the check
        # has just been performed above:
        return self.element_class(self, coords=coords, chart=chart,
                                  check_coords=False)

    def __contains__(self, point):
        r"""
        Check whether a point is contained in the manifold.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,2), chart=X)
            sage: M.__contains__(p)
            True
            sage: p in M  # indirect doctest
            True
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: U.__contains__(p)
            True
            sage: p in U  # indirect doctest
            True
            sage: V = U.open_subset('V', coord_def={X.restrict(U): y<0})
            sage: V.__contains__(p)
            False
            sage: p in V  # indirect doctest
            False

        """
        # for efficiency, a quick test first:
        if point.parent() is self:
            return True
        if point.parent().is_subset(self):
            return True
        for chart in self._atlas:
            if chart in point._coordinates:
                if chart.valid_coordinates( *(point._coordinates[chart]) ):
                    return True
        for chart in point._coordinates:
            for schart in chart._subcharts:
                if schart in self._atlas and schart.valid_coordinates(
                                          *(point._coordinates[chart]) ):
                    return True
        return False

    def open_subset(self, name, latex_name=None, coord_def={}):
        r"""
        Create an open subset of the manifold.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology. It is a topological
        manifold by itself. Hence the returned object is an instance of
        :class:`TopologicalManifold`.

        INPUT:

        - ``name`` -- name given to the open subset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote
          the subset; if none are provided, it is set to ``name``
        - ``coord_def`` -- (default: {}) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts on the manifold and values the symbolic expressions formed
          by the coordinates to define the subset

        OUTPUT:

        - the open subset, as an instance of :class:`TopologicalManifold`

        EXAMPLES:

        Creating an open subset of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.open_subset('A'); A
            Open subset A of the 2-dimensional topological manifold M

        As an open subset of a topological manifold, ``A`` is itself a
        topological manifold, on the same topological field and of the same
        dimension as ``M``::

            sage: isinstance(A, sage.manifolds.manifold.TopologicalManifold)
            True
            sage: A.base_field() == M.base_field()
            True
            sage: dim(A) == dim(M)
            True
            sage: A.category() is M.category().Subobjects()
            True

        Creating an open subset of ``A``::

            sage: B = A.open_subset('B'); B
            Open subset B of the 2-dimensional topological manifold M

        We have then::

            sage: A.subsets()  # random (set output)
            {Open subset B of the 2-dimensional topological manifold M,
             Open subset A of the 2-dimensional topological manifold M}
            sage: B.is_subset(A)
            True
            sage: B.is_subset(M)
            True

        Defining an open subset by some coordinate restrictions: the open
        unit disk in `\RR^2`::

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1}); U
            Open subset U of the 2-dimensional topological manifold R^2

        Since the argument ``coord_def`` has been set, ``U`` is automatically
        provided with a chart, which is the restriction of the Cartesian one
        to ``U``::

            sage: U.atlas()
            [Chart (U, (x, y))]

        Therefore, one can immediately check whether a point belongs
        to ``U``::

            sage: M.point((0,0)) in U
            True
            sage: M.point((1/2,1/3)) in U
            True
            sage: M.point((1,2)) in U
            False

        """
        resu = TopologicalManifold(self._dim, name, self._field,
                                   self._structure, ambient=self._manifold,
                                   latex_name=latex_name,
                                   start_index=self._sindex)
        resu._supersets.update(self._supersets)
        for sd in self._supersets:
            sd._subsets.add(resu)
        self._top_subsets.add(resu)
        # Charts on the result from the coordinate definition:
        for chart, restrictions in coord_def.iteritems():
            if chart not in self._atlas:
                raise ValueError("the {} does not belong to ".format(chart) +
                                 "the atlas of {}".format(self))
            chart.restrict(resu, restrictions)
        # Transition maps on the result inferred from those of self:
        for chart1 in coord_def:
            for chart2 in coord_def:
                if chart2 != chart1 and (chart1, chart2) in self._coord_changes:
                    self._coord_changes[(chart1, chart2)].restrict(resu)
        return resu

    def get_chart(self, coordinates, domain=None):
        r"""
        Get a chart from its coordinates.

        The chart must have been previously created by the method
        :meth:`chart`.

        INPUT:

        - ``coordinates`` --  single string composed of the coordinate symbols
          separated by a space
        - ``domain`` -- (default: ``None``) string containing the name of the
          chart's domain, which must be a subset of the current manifold; if
          ``None``, the current manifold is assumed

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.chart.Chart` (or of the subclass
          :class:`~sage.manifolds.chart.RealChart`) representing the chart
          corresponding to the above specifications

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: M.get_chart('x y')
            Chart (M, (x, y))
            sage: M.get_chart('x y') is X
            True
            sage: U = M.open_subset('U', coord_def={X: (y!=0,x<0)})
            sage: Y.<r, ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: M.atlas()
            [Chart (M, (x, y)), Chart (U, (x, y)), Chart (U, (r, ph))]
            sage: M.get_chart('x y', domain='U')
            Chart (U, (x, y))
            sage: M.get_chart('x y', domain='U') is X.restrict(U)
            True
            sage: U.get_chart('r ph')
            Chart (U, (r, ph))
            sage: M.get_chart('r ph', domain='U')
            Chart (U, (r, ph))
            sage: M.get_chart('r ph', domain='U') is Y
            True

        """
        if domain is None:
            dom = self
        else:
            dom = self.get_subset(domain)
        try:
            return dom._charts_by_coord[coordinates]
        except KeyError:
            raise KeyError("the coordinates '{}' ".format(coordinates) +
                           "do not correspond to any chart with " +
                           "the {} as domain".format(dom))

    def dimension(self):
        r"""
        Return the dimension of the manifold over its base field.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: M.dimension()
            2

        A shortcut is ``dim()``::

            sage: M.dim()
            2

        The Sage global function ``dim`` can also be used::

            sage: dim(M)
            2

        """
        return self._dim

    dim = dimension

    def base_field(self):
        r"""
        Return the field on which the manifold is defined.

        OUTPUT:

        - a topological field

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: M.base_field()
            Real Field with 53 bits of precision
            sage: M = Manifold(3, 'M', structure='topological', field='complex')
            sage: M.base_field()
            Complex Field with 53 bits of precision
            sage: M = Manifold(3, 'M', structure='topological', field=QQ)
            sage: M.base_field()
            Rational Field

        """
        return self._field

    def base_field_type(self):
        r"""
        Return the type of topological field on which the manifold is defined.

        OUTPUT:

        - a string describing the field, with three possible values:

          - ``'real'`` for the real field `\RR`
          - ``'complex'`` for the complex field `\CC`
          - ``'neither_real_nor_complex'`` for a field different from `\RR`
            and `\CC`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: M.base_field_type()
            'real'
            sage: M = Manifold(3, 'M', structure='topological', field='complex')
            sage: M.base_field_type()
            'complex'
            sage: M = Manifold(3, 'M', structure='topological', field=QQ)
            sage: M.base_field_type()
            'neither_real_nor_complex'

        """
        return self._field_type

    def start_index(self):
        r"""
        Return the first value of the index range used on the manifold.

        This is the parameter ``start_index`` passed at the construction of
        the manifold.

        OUTPUT:

        - the integer `i_0` such that all indices of indexed objects on the
          manifold range from `i_0` to `i_0 + n - 1`, where `n` is the
          manifold's dimension

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: M.start_index()
            0
            sage: M = Manifold(3, 'M', structure='topological', start_index=1)
            sage: M.start_index()
            1

        """
        return self._sindex

    def irange(self, start=None):
        r"""
        Single index generator.

        INPUT:

        - ``start`` -- (default: ``None``) initial value `i_0` of the index;
          if none are provided, the value returned by :meth:`start_index()`
          is assumed

        OUTPUT:

        - an iterable index, starting from `i_0` and ending at
          `i_0 + n - 1`, where `n` is the manifold's dimension

        EXAMPLES:

        Index range on a 4-dimensional manifold::

            sage: M = Manifold(4, 'M', structure='topological')
            sage: for i in M.irange():
            ....:     print i,
            ....:
            0 1 2 3
            sage: for i in M.irange(2):
            ....:     print i,
            ....:
            2 3
            sage: list(M.irange())
            [0, 1, 2, 3]

        Index range on a 4-dimensional manifold with starting index=1::

            sage: M = Manifold(4, 'M', structure='topological', start_index=1)
            sage: for i in M.irange():
            ....:     print i,
            ....:
            1 2 3 4
            sage: for i in M.irange(2):
            ....:    print i,
            ....:
            2 3 4

        In general, one has always::

            sage: M.irange().next() == M.start_index()
            True

        """
        si = self._sindex
        imax = self._dim + si
        if start is None:
            i = si
        else:
            i = start
        while i < imax:
            yield i
            i += 1

    def index_generator(self, nb_indices):
        r"""
        Generator of index series.

        INPUT:

        - ``nb_indices`` -- number of indices in a series

        OUTPUT:

        - an iterable index series for a generic component with the specified
          number of indices

        EXAMPLES:

        Indices on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological', start_index=1)
            sage: for ind in M.index_generator(2):
            ....:     print ind
            ....:
            (1, 1)
            (1, 2)
            (2, 1)
            (2, 2)

        Loops can be nested::

            sage: for ind1 in M.index_generator(2):
            ....:     print ind1, " : ",
            ....:     for ind2 in M.index_generator(2):
            ....:         print ind2,
            ....:     print ""
            ....:
            (1, 1)  :  (1, 1) (1, 2) (2, 1) (2, 2)
            (1, 2)  :  (1, 1) (1, 2) (2, 1) (2, 2)
            (2, 1)  :  (1, 1) (1, 2) (2, 1) (2, 2)
            (2, 2)  :  (1, 1) (1, 2) (2, 1) (2, 2)

        """
        si = self._sindex
        imax = self._dim - 1 + si
        ind = [si for k in range(nb_indices)]
        ind_end = [si for k in range(nb_indices)]
        ind_end[0] = imax+1
        while ind != ind_end:
            yield tuple(ind)
            ret = 1
            for pos in range(nb_indices-1,-1,-1):
                if ind[pos] != imax:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = imax + 1 # end point reached
                    else:
                        ind[pos] = si
                        ret = 1

    def atlas(self):
        r"""
        Return the list of charts that have been defined on the manifold.

        EXAMPLES:

        Let us consider `\RR^2` as a 2-dimensional manifold::

            sage: M = Manifold(2, 'R^2', structure='topological')

        Immediately after the manifold creation, the atlas is empty, since no
        chart has been defined yet::

            sage: M.atlas()
            []

        Let us introduce the chart of Cartesian coordinates::

            sage: c_cart.<x,y> = M.chart()
            sage: M.atlas()
            [Chart (R^2, (x, y))]

        The complement of the half line `\{y = 0, x \geq 0\}`::

            sage: U = M.open_subset('U', coord_def={c_cart: (y!=0,x<0)})
            sage: U.atlas()
            [Chart (U, (x, y))]
            sage: M.atlas()
            [Chart (R^2, (x, y)), Chart (U, (x, y))]

        Spherical (polar) coordinates on ``U``::

            sage: c_spher.<r, ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: U.atlas()
            [Chart (U, (x, y)), Chart (U, (r, ph))]
            sage: M.atlas()
            [Chart (R^2, (x, y)), Chart (U, (x, y)), Chart (U, (r, ph))]

        .. SEEALSO::

            :meth:`top_charts`

        """
        return list(self._atlas) # Make a (shallow) copy

    def top_charts(self):
        r"""
        Return the list of charts defined on subsets of the current manifold
        that are not subcharts of charts on larger subsets.

        OUTPUT:

        - list of charts defined on open subsets of the manifold but not on
          larger subsets

        EXAMPLES:

        Charts on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: Y.<u,v> = U.chart()
            sage: M.top_charts()
            [Chart (M, (x, y)), Chart (U, (u, v))]

        Note that the (user) atlas contains one more chart: ``(U, (x,y))``,
        which is not a "top" chart::

            sage: M.atlas()
            [Chart (M, (x, y)), Chart (U, (x, y)), Chart (U, (u, v))]

        .. SEEALSO::

            :meth:`atlas` for the complete list of charts defined on the
            manifold.

        """
        return list(self._top_charts) # Make a (shallow) copy

    def default_chart(self):
        r"""
        Return the default chart defined on the manifold.

        Unless changed via :meth:`set_default_chart`, the *default chart*
        is the first one defined on a subset of the manifold (possibly itself).

        OUTPUT:

        - instance of :class:`~sage.manifolds.chart.Chart`
          representing the default chart

        EXAMPLES:

        Default chart on a 2-dimensional manifold and on some subsets::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: M.chart('x y')
            Chart (M, (x, y))
            sage: M.chart('u v')
            Chart (M, (u, v))
            sage: M.default_chart()
            Chart (M, (x, y))
            sage: A = M.open_subset('A')
            sage: A.chart('t z')
            Chart (A, (t, z))
            sage: A.default_chart()
            Chart (A, (t, z))

        """
        return self._def_chart

    def set_default_chart(self, chart):
        r"""
        Changing the default chart on ``self``.

        INPUT:

        - ``chart`` -- a chart (must be defined on some subset ``self``)

        EXAMPLES:

        Charts on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: M.default_chart()
            Chart (M, (x, y))
            sage: M.set_default_chart(c_uv)
            sage: M.default_chart()
            Chart (M, (u, v))

        """
        from chart import Chart
        if not isinstance(chart, Chart):
            raise TypeError("{} is not a chart".format(chart))
        if chart._domain is not self:
            if self.is_manifestly_coordinate_domain():
                raise TypeError("the chart domain must coincide with " +
                                "the {}".format(self))
            if chart not in self._atlas:
                raise ValueError("the chart must be defined on the " +
                                 "{}".format(self))
        self._def_chart = chart

    def coord_change(self, chart1, chart2):
        r"""
        Return the change of coordinates (transition map) between two charts
        defined on the manifold.

        The change of coordinates must have been defined previously, for
        instance by the method
        :meth:`~sage.manifolds.chart.Chart.transition_map`.

        INPUT:

        - ``chart1`` -- chart 1
        - ``chart2`` -- chart 2

        OUTPUT:

        - instance of :class:`~sage.manifolds.chart.CoordChange`
          representing the transition map from chart 1 to chart 2

        EXAMPLES:

        Change of coordinates on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: c_xy.transition_map(c_uv, (x+y, x-y)) # defines the coord. change
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
            sage: M.coord_change(c_xy, c_uv) # returns the coord. change defined above
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))

        """
        if (chart1, chart2) not in self._coord_changes:
            raise TypeError("the change of coordinates from " +
                            "{} to {}".format(chart1, chart2) + " has not " +
                            "been defined on the {}".format(self))
        return self._coord_changes[(chart1, chart2)]

    def coord_changes(self):
        r"""
        Return the changes of coordinates (transition maps) defined on
        subsets of the manifold.

        OUTPUT:

        - dictionary of changes of coordinates, with pairs of charts as keys

        EXAMPLES:

        Various changes of coordinates on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, [x+y, x-y])
            sage: M.coord_changes()
            {(Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))}
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: M.coord_changes()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))}
            sage: c_rs.<r,s> = M.chart()
            sage: uv_to_rs = c_uv.transition_map(c_rs, [-u+2*v, 3*u-v])
            sage: M.coord_changes()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (r, s)),
             (Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))}
            sage: xy_to_rs = uv_to_rs * xy_to_uv
            sage: M.coord_changes()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (r, s)),
             (Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v)),
             (Chart (M, (x, y)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (r, s))}

        """
        return self._coord_changes

    def is_manifestly_coordinate_domain(self):
        r"""
        Return ``True`` if the manifold is known to be the domain of some
        coordinate chart and ``False`` otherwise.

        If ``False`` is returned, either the manifold cannot be the domain of
        some coordinate chart or no such chart has been declared yet.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: X.<x,y> = U.chart()
            sage: U.is_manifestly_coordinate_domain()
            True
            sage: M.is_manifestly_coordinate_domain()
            False
            sage: Y.<u,v> = M.chart()
            sage: M.is_manifestly_coordinate_domain()
            True

        """
        return bool(self._covering_charts)

    def chart(self, coordinates='', names=None):
        r"""
        Define a chart, the domain of which is the manifold.

        A *chart* is a pair `(U, \varphi)`, where `U` is the current
        manifold and `\varphi: U \rightarrow V \subset K^n`
        is a homeomorphism from `U` to an open subset `V` of `K^n`, `K`
        being the field on which the manifold is defined.

        The components `(x^1, \ldots, x^n)` of `\varphi`, defined by
        `\varphi(p) = (x^1(p), \ldots, x^n(p)) \in K^n` for any point
        `p \in U`, are called the *coordinates* of the chart `(U, \varphi)`.

        See :class:`~sage.manifolds.chart.Chart` for a complete
        documentation.

        INPUT:

        - ``coordinates`` --  (default: ``''`` (empty string)) string
          defining the coordinate symbols and ranges, see below
        - ``names`` -- (default: ``None``) unused argument, except if
          ``coordinates`` is not provided; it must then be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        The coordinates declared in the string ``coordinates`` are
        separated by ``' '`` (whitespace) and each coordinate has at most three
        fields, separated by a colon (``':'``):

        1. The coordinate symbol (a letter or a few letters).
        2. (optional, only for manifolds over `\RR`) The interval `I`
           defining the coordinate range: if not provided, the coordinate
           is assumed to span all `\RR`; otherwise `I` must be provided
           in the form ``(a,b)`` (or equivalently ``]a,b[``)
           The bounds ``a`` and ``b`` can be ``+/-Infinity``, ``Inf``,
           ``infinity``, ``inf`` or ``oo``. For *singular* coordinates,
           non-open intervals such as ``[a,b]`` and
           ``(a,b]`` (or equivalently ``]a,b]``) are allowed. Note that
           the interval declaration must not contain any space character.
        3. (optional) The LaTeX spelling of the coordinate; if not provided
           the coordinate symbol given in the first field will be used.

        The order of the fields 2 and 3 does not matter and each of them can
        be omitted. If it contains any LaTeX expression, the string
        ``coordinates`` must be declared with the prefix 'r' (for "raw") to
        allow for a proper treatment of the backslash character (see
        examples below). If no interval range and no LaTeX spelling is to
        be provided for any coordinate, the argument ``coordinates`` can be
        omitted when the shortcut operator ``<,>`` is used via Sage
        preparser (see examples below).

        OUTPUT:

        - the created chart, as an instance of
          :class:`~sage.manifolds.chart.Chart` or of the subclass
          :class:`~sage.manifolds.chart.RealChart` for manifolds over `\RR`.

        EXAMPLES:

        Chart on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: X = U.chart('x y'); X
            Chart (U, (x, y))
            sage: X[0]
            x
            sage: X[1]
            y
            sage: X[:]
            (x, y)

        The declared coordinates are not known at the global level::

            sage: y
            Traceback (most recent call last):
            ...
            NameError: name 'y' is not defined

        They can be recovered by the operator ``[:]`` applied to the chart::

            sage: (x, y) = X[:]
            sage: y
            y
            sage: type(y)
            <type 'sage.symbolic.expression.Expression'>

        But a shorter way to proceed is to use the operator ``<,>`` in the
        left-hand side of the chart declaration (there is then no need to
        pass the string 'x y' to chart())::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: X.<x,y> = U.chart(); X
            Chart (U, (x, y))

        Indeed, the declared coordinates are then known at the global level::

            sage: y
            y
            sage: (x,y) == X[:]
            True

        Actually the instruction ``X.<x,y> = U.chart()`` is
        equivalent to the combination of the two instructions
        ``X = U.chart('x y')`` and ``(x,y) = X[:]``.

        See the documentation of class
        :class:`~sage.manifolds.chart.Chart` for more examples,
        especially regarding the coordinates ranges and restrictions.

        """
        return self._structure.chart(self, coordinates=coordinates, names=names)

    def is_open(self):
        """
        Return if ``self`` is an open set.

        In the present case (manifold or open subset of it), always
        return ``True``.

        TEST::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: M.is_open()
            True

        """
        return True

##############################################################################
## Constructor function

def Manifold(dim, name, latex_name=None, field='real', structure='smooth',
             start_index=0, **extra_kwds):
    r"""
    Construct a manifold of a given type over a topological field `K`.

    Given a topological field `K` (in most applications, `K = \RR` or
    `K = \CC`) and a non-negative integer `n`, a *topological manifold of
    dimension* `n` *over K* is a topological space `M` such that

    - `M` is a Hausdorff space,
    - `M` is second countable, and
    - every point in `M` has a neighborhood homeomorphic to `K^n`.

    A *real manifold* is a manifold over `\RR`. A *differentiable* (resp.
    *smooth*, resp. *analytic*) is a real manifold such that all transition
    maps are *differentiable* (resp. *smooth*, resp. *analytic*).

    INPUT:

    - ``dim`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the manifold; if none are provided, it is set to ``name``
    - ``field`` -- (default: ``'real'``) field `K` on which the
      manifold is defined; allowed values are

      - ``'real'`` or an object of type ``RealField`` (e.g. ``RR``) for a
        manifold over `\RR`
      - ``'complex'`` or an object of type ``ComplexField`` (e.g. ``CC``)
        for a manifold over `\CC`
      - an object in the category of topological fields (see
        :class:`~sage.categories.fields.Fields` and
        :class:`~sage.categories.topological_spaces.TopologicalSpaces`)
        for other types of manifolds

    - ``structure`` -- (default: ``'smooth'``) to specify the structure or
      type of manifold; allowed values are

      - ``'topological'`` or ``'top'`` for a topological manifold
      - ``'differentiable'`` or ``'diff'`` for a differentiable manifold
      - ``'smooth'`` for a smooth manifold
      - ``'analytic'`` for an analytic manifold

    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the manifold, e.g. coordinates
      in a chart
    - ``extra_kwds`` -- keywords meaningful only for some specific types
      of manifolds

    OUTPUT:

    - a manifold of the specified type

    EXAMPLES:

    A 3-dimensional real topological manifold::

        sage: M = Manifold(3, 'M', structure='topological'); M
        3-dimensional topological manifold M

    Given the default value of the parameter ``field``, the above is
    equivalent to::

        sage: M = Manifold(3, 'M', structure='topological', field='real'); M
        3-dimensional topological manifold M

    A complex topological manifold::

        sage: M = Manifold(3, 'M', structure='topological', field='complex'); M
        Complex 3-dimensional topological manifold M

    A topological manifold over `\QQ`::

        sage: M = Manifold(3, 'M', structure='topological', field=QQ); M
        3-dimensional topological manifold M over the Rational Field

    See the documentation of class
    :class:`~sage.manifolds.manifold.TopologicalManifold` for more
    detailed examples.

    .. RUBRIC:: Uniqueness of manifold objects

    Suppose we construct a manifold named `M`::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: X.<x,y> = M.chart()

    At some point, we change our mind and would like to restart with a new
    manifold, using the same name `M` and keeping the previous manifold for
    reference::

        sage: M_old = M  # for reference
        sage: M = Manifold(2, 'M', structure='topological')

    This results in a brand new object::

        sage: M.atlas()
        []

    The object ``M_old`` is intact::

        sage: M_old.atlas()
        [Chart (M, (x, y))]

    Both objects have the same display::

        sage: M
        2-dimensional topological manifold M
        sage: M_old
        2-dimensional topological manifold M

    but they are different::

        sage: M != M_old
        True

    Let us introduce a chart on ``M``, using the same coordinate symbols
    as for ``M_old``::

        sage: X.<x,y> = M.chart()

    The charts are displayed in the same way::

        sage: M.atlas()
        [Chart (M, (x, y))]
        sage: M_old.atlas()
        [Chart (M, (x, y))]

    but they are actually different::

        sage: M.atlas()[0] != M_old.atlas()[0]
        True

    Moreover, the two manifolds ``M`` and ``M_old`` are still considered
    distinct::

        sage: M != M_old
        True

    This reflects the fact that the equality of manifold objects holds only
    for identical objects, i.e. one has ``M1 == M2`` if, and only if,
    ``M1 is M2``. Actually, the manifold classes inherit from
    :class:`~sage.misc.fast_methods.WithEqualityById`::

        sage: isinstance(M, sage.misc.fast_methods.WithEqualityById)
        True
    """
    from time import time
    # Some sanity checks
    if not isinstance(dim, (int, Integer)):
        raise TypeError("the manifold dimension must be an integer")
    if dim < 1:
        raise ValueError("the manifold dimension must be strictly positive")

    if structure in ['topological', 'top']:
        if field == 'real' or isinstance(field, RealField_class):
            structure = RealTopologicalStructure()
        else:
            structure = TopologicalStructure()
    else:
        raise NotImplementedError("manifolds of type {} are not ".format(structure) +
                                  "implemented")
    return TopologicalManifold(dim, name, field, structure,
                               latex_name=latex_name, start_index=start_index,
                               unique_tag=getrandbits(128)*time())

