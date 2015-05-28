r"""
Topological manifolds

Given a topological field `K` (in most applications, `K = \RR` or
`K = \CC`) and a non-negative integer `n`, a *topological manifold of
dimension* `n` *over K* is a topological space `M` such that

- `M` is a Hausdorff space,
- `M` is second countable,
- every point in `M` has a neighborhood homeomorphic to `K^n`

Topological manifolds are implemented via the class :class:`TopManifold`.
In the current setting, topological manifolds are mostly described by means of
charts (see :class:`~sage.manifolds.chart.Chart`).

:class:`TopManifold` serves as a base class for more specific manifold classes.

.. RUBRIC:: Example 1: the 2-sphere as a topological manifold of dimension
  2 over `\RR`

One starts by declaring `S^2` as a 2-dimensional topological manifold::

    sage: M = TopManifold(2, 'S^2')
    sage: M
    2-dimensional topological manifold S^2

Since the base topological field has not been specified in the argument list
of ``TopManifold``, `\RR` is assumed::

    sage: M.base_field()
    'real'
    sage: dim(M)
    2

Let us consider the complement of a point, the "North pole" say; this is an
open subset of `S^2`, which we call U::

    sage: U = M.open_subset('U'); U
    Open subset U of the 2-dimensional topological manifold S^2

A standard chart on U is provided by the stereographic projection from the
North pole to the equatorial plane::

    sage: stereoN.<x,y> = U.chart(); stereoN
    Chart (U, (x, y))

Thanks to the operator ``<x,y>`` on the left-hand side, the coordinates
declared in a chart (here x and y), are accessible by their names; they are
Sage's symbolic variables::

    sage: y
    y
    sage: type(y)
    <type 'sage.symbolic.expression.Expression'>

The South pole is the point of coordinates `(x,y)=(0,0)` in the above
chart::

    sage: S = U.point((0,0), chart=stereoN, name='S'); S
    Point S on the 2-dimensional topological manifold S^2

Let us call V the open subset that is the complement of the South pole and let
us introduce on it the chart induced by the stereographic projection from
the South pole to the equatorial plane::

    sage: V = M.open_subset('V'); V
    Open subset V of the 2-dimensional topological manifold S^2
    sage: stereoS.<u,v> = V.chart(); stereoS
    Chart (V, (u, v))

The North pole is the point of coordinates `(u,v)=(0,0)` in this chart::

    sage: N = V.point((0,0), chart=stereoS, name='N'); N
    Point N on the 2-dimensional topological manifold S^2

To fully construct the manifold, we declare that it is the union of U
and V::

    sage: M.declare_union(U,V)

and we provide the transition map between the charts ``stereoN`` = `(U, (x, y))`
and ``stereoS`` = `(V, (u, v))`, denoting by W the intersection of U and V
(W is the subset of U defined by `x^2+y^2\not=0`, as well as the subset of V
defined by`u^2+v^2\not=0`)::

    sage: stereoN_to_S = stereoN.transition_map(stereoS, [x/(x^2+y^2), y/(x^2+y^2)],
    ....:                intersection_name='W', restrictions1= x^2+y^2!=0, restrictions2= u^2+v^2!=0)
    sage: stereoN_to_S
    Change of coordinates from Chart (W, (x, y)) to Chart (W, (u, v))
    sage: stereoN_to_S.display()
    u = x/(x^2 + y^2)
    v = y/(x^2 + y^2)

We give the name W for the Python variable representing `W=U\cap V`::

    sage: W = U.intersection(V)

The inverse of the transition map is computed by the method inverse()::

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
    [Chart (U, (x, y)), Chart (V, (u, v)), Chart (W, (x, y)), Chart (W, (u, v))]

Let us consider the point of coordinates (1,2) in the chart ``stereoN``::

    sage: p = M.point((1,2), chart=stereoN, name='p'); p
    Point p on the 2-dimensional topological manifold S^2
    sage: p.parent()
    2-dimensional topological manifold S^2
    sage: p in W
    True

The coordinates of `p` in the chart ``stereoS`` are::

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

    sage: M = TopManifold(1, 'C*', field='complex'); M
    Complex 1-dimensional topological manifold C*
    sage: U = M.open_subset('U')
    sage: Z.<z> = U.chart()
    sage: O = U.point((0,), chart=Z, name='O'); O
    Point O on the Complex 1-dimensional topological manifold C*
    sage: V = M.open_subset('V')
    sage: W.<w> = V.chart(); W
    Chart (V, (w,))
    sage: inf = M.point((0,), chart=W, name='inf', latex_name=r'\infty')
    sage: inf
    Point inf on the Complex 1-dimensional topological manifold C*
    sage: Z_to_W = Z.transition_map(W, 1/z, intersection_name='A', restrictions1= z!=0, restrictions2= w!=0)
    sage: Z_to_W
    Change of coordinates from Chart (A, (z,)) to Chart (A, (w,))
    sage: Z_to_W.display()
    w = 1/z
    sage: Z_to_W.inverse()
    Change of coordinates from Chart (A, (w,)) to Chart (A, (z,))
    sage: Z_to_W.inverse().display()
    z = 1/w
    sage: M.list_of_subsets()
    [Open subset A of the Complex 1-dimensional topological manifold C*,
     Complex 1-dimensional topological manifold C*,
     Open subset U of the Complex 1-dimensional topological manifold C*,
     Open subset V of the Complex 1-dimensional topological manifold C*]
    sage: M.atlas()
    [Chart (U, (z,)), Chart (V, (w,)), Chart (A, (z,)), Chart (A, (w,))]


AUTHORS:

- Eric Gourgoulhon (2015): initial version

REFERENCES:

- J.M. Lee : *Introduction to Topological Manifolds*, 2nd ed., Springer (New
  York) (2011)

"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.fields import Fields
from sage.manifolds.subset import TopManifoldOpenSubset

class TopManifold(TopManifoldOpenSubset):
    r"""
    Topological manifold over a topological field `K`.

    Given a topological field `K` (in most applications, `K = \RR` or
    `K = \CC`) and a non-negative integer `n`, a *topological manifold of
    dimension* `n` *over K* is a topological space `M` such that

    - `M` is a Hausdorff space,
    - `M` is second countable,
    - every point in `M` has a neighborhood homeomorphic to `K^n`

    This is a Sage *parent* class, the corresponding *element*
    class being :class:`~sage.manifolds.point.TopManifoldPoint`.

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      manifold; if none is provided, it is set to ``name``
    - ``field`` -- (default: 'real') field `K` on which the manifold is
      defined; allowed values are

        - 'real' for a manifold over `\RR`
        - 'complex' for a manifold over `\CC`
        - any object in the category of fields (see
          :class:`~sage.categories.fields.Fields`) for more general manifolds

    - ``start_index`` -- (default: 0) integer; lower bound of the range of
      indices used for "indexed objects" on the manifold, e.g. coordinates
      in a chart.

    EXAMPLES:

    A 4-dimensional topological manifold (over `\RR`)::

        sage: M = TopManifold(4, 'M', latex_name=r'\mathcal{M}')
        sage: M
        4-dimensional topological manifold M
        sage: latex(M)
        \mathcal{M}
        sage: M.base_field()
        'real'

    The input parameter ``start_index`` defines the range of indices on the
    manifold::

        sage: M = TopManifold(4, 'M')
        sage: list(M.irange())
        [0, 1, 2, 3]
        sage: M = TopManifold(4, 'M', start_index=1)
        sage: list(M.irange())
        [1, 2, 3, 4]
        sage: list(TopManifold(4, 'M', start_index=-2).irange())
        [-2, -1, 0, 1]

    A complex manifold::

        sage: N = TopManifold(3, 'N', field='complex'); N
        Complex 3-dimensional topological manifold N

    A manifold over `\QQ`::

        sage: N = TopManifold(6, 'N', field=QQ); N
        6-dimensional topological manifold N over the Rational Field

    A manifold is a Sage *parent* object, in the category of sets::

        sage: isinstance(M, Parent)
        True
        sage: M.category()
        Category of sets
        sage: M in Sets()
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
    :class:`~sage.manifolds.point.TopManifoldPoint`::

        sage: isinstance(p, sage.manifolds.point.TopManifoldPoint)
        True

    Manifolds are unique, as long as they are created with the same arguments::

        sage: M is TopManifold(4, 'M', start_index=1)
        True
        sage: M is TopManifold(4, 'M')
        False
        sage: M is TopManifold(4, 'M', latex_name='M', start_index=1)
        False

    The manifold passes all the tests of the test suite relative to the
    category of Sets::

        sage: TestSuite(M).run(verbose=True)
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass


    """
    def __init__(self, n, name, latex_name=None, field='real', start_index=0):
        r"""
        Construct a topological manifold.

        TESTS::

            sage: M = TopManifold(3, 'M', latex_name=r'\mathbb{M}', start_index=1)
            sage: M
            3-dimensional topological manifold M
            sage: latex(M)
            \mathbb{M}
            sage: dim(M)
            3
            sage: X.<x,y,z> = M.chart()
            sage: TestSuite(M).run()

        """
        # Initialization of the attributes _dim, _field and _start_index:
        self._init_dim(n, field, start_index)
        # Initialization as an open set of itself:
        TopManifoldOpenSubset.__init__(self, self, name, latex_name=latex_name)

    def _init_dim(self, n, field, start_index):
        r"""
        Check and initialization of the manifold attributes ``_dim``,
        ``_field`` and ``_start_index``.

        This method shall be called by the method ``__init__`` of any derived
        class of :class:`TopManifold` (if :meth:`TopManifold.__init__` is not
        called).

        """
        from sage.rings.integer import Integer
        if not isinstance(n, (int, Integer)):
            raise TypeError("the manifold dimension must be an integer")
        if n<1:
            raise ValueError("the manifold dimension must be strictly " +
                             "positive")
        self._dim = n
        if field not in ['real', 'complex']:
            if field not in Fields():
                raise TypeError("the argument 'field' must be a field")
        self._field = field
        if not isinstance(start_index, (int, Integer)):
            raise TypeError("the start index must be an integer")
        self._sindex = start_index


    def _repr_(self):
        r"""
        String representation of the manifold.

        TESTS::

            sage: M = TopManifold(3, 'M')
            sage: M._repr_()
            '3-dimensional topological manifold M'
            sage: repr(M)
            '3-dimensional topological manifold M'
            sage: M
            3-dimensional topological manifold M
            sage: M = TopManifold(3, 'M', field='complex')
            sage: M._repr_()
            'Complex 3-dimensional topological manifold M'
            sage: M = TopManifold(3, 'M', field=QQ)
            sage: M._repr_()
            '3-dimensional topological manifold M over the Rational Field'

        """
        if self._field == 'real':
            return "{}-dimensional topological manifold {}".format(self._dim,
                                                                   self._name)
        elif self._field == 'complex':
            return "Complex {}-dimensional topological manifold {}".format(
                                                         self._dim, self._name)
        return "{}-dimensional topological manifold {} over the {}".format(
                                            self._dim, self._name, self._field)

    def _latex_(self):
        r"""
        LaTeX representation of the manifold.

        TESTS::

            sage: M = TopManifold(3, 'M')
            sage: M._latex_()
            'M'
            sage: latex(M)
            M
            sage: M = TopManifold(3, 'M', latex_name=r'\mathcal{M}')
            sage: M._latex_()
            '\\mathcal{M}'
            sage: latex(M)
            \mathcal{M}

        """
        return self._latex_name

    def dimension(self):
        r"""
        Return the dimension of the manifold over its base field.

        EXAMPLE::

            sage: M = TopManifold(2, 'M')
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
        Return the field on which the manifolds is defined.

        OUTPUT:

        - a field or one of the strings 'real' and 'complex', since there
          is no exact representation of the fields `\RR` and `\CC` in Sage

        EXAMPLES::

            sage: M = TopManifold(3, 'M')
            sage: M.base_field()
            'real'
            sage: M = TopManifold(3, 'M', field='complex')
            sage: M.base_field()
            'complex'
            sage: M = TopManifold(3, 'M', field=QQ)
            sage: M.base_field()
            Rational Field

        """
        return self._field

    def irange(self, start=None):
        r"""
        Single index generator.

        INPUT:

        - ``start`` -- (default: ``None``) initial value of the index; if none is
          provided, ``self._sindex`` is assumed

        OUTPUT:

        - an iterable index, starting from ``start`` and ending at
          ``self._sindex + self._dim -1``

        EXAMPLES:

        Index range on a 4-dimensional manifold::

            sage: M = TopManifold(4, 'M')
            sage: for i in M.irange():
            ...       print i,
            ...
            0 1 2 3
            sage: for i in M.irange(2):
            ...       print i,
            ...
            2 3
            sage: list(M.irange())
            [0, 1, 2, 3]

        Index range on a 4-dimensional manifold with starting index=1::

            sage: M = TopManifold(4, 'M', start_index=1)
            sage: for i in M.irange():
            ...       print i,
            ...
            1 2 3 4
            sage: for i in M.irange(2):
            ...      print i,
            ...
            2 3 4

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

            sage: M = TopManifold(2, 'M', start_index=1)
            sage: for ind in M.index_generator(2):
            ...       print ind
            ...
            (1, 1)
            (1, 2)
            (2, 1)
            (2, 2)

        Loops can be nested::

            sage: for ind1 in M.index_generator(2):
            ...       print ind1, " : ",
            ...       for ind2 in M.index_generator(2):
            ...           print ind2,
            ...       print ""
            ...
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
