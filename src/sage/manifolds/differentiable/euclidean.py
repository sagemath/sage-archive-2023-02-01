r"""
Euclidean spaces

An *Euclidean space of dimension* `n` is a Riemannian manifold diffeomorphic to
`\RR^n` equipped with a flat metric.

Euclidean spaces are implemented via the following classes

- :class:`EuclideanSpaceGeneric`

  - :class:`EuclideanPlane` (case `n=2`)
  - :class:`EuclideanSpace3dim` (case `n=3`)

The user interface is provided by the generic function :func:`EuclideanSpace`.

.. RUBRIC:: Example 1: the Euclidean plane

We start by declaring the Euclidean plane ``E``, with ``(x,y)`` as
Cartesian coordinates::

    sage: E.<x,y> = EuclideanSpace(2)
    sage: E
    Euclidean plane E^2
    sage: dim(E)
    2

``E`` is automatically endowed with the chart of Cartesian coordinates::

    sage: E.atlas()
    [Chart (E^2, (x, y))]
    sage: X_cartesian = E.default_chart(); X_cartesian
    Chart (E^2, (x, y))

Thanks to the use of ``<x,y>`` when declaring ``E``, the coordinates `(x,y)`
have been injected in the global name space, i.e. the Python variables ``x``
and ``y`` have been created and are available to form symbolic expressions::

    sage: y
    y
    sage: type(y)
    <type 'sage.symbolic.expression.Expression'>
    sage: assumptions()
    [x is real, y is real]

The metric tensor of ``E`` is predefined::

    sage: g = E.metric(); g
    Riemannian metric g on the Euclidean plane E^2
    sage: g.display()
    g = dx*dx + dy*dy
    sage: g[:]
    [1 0]
    [0 1]

It is a *flat* metric, i.e. it has a vanishing Riemann tensor::

    sage: g.riemann()
    Tensor field Riem(g) of type (1,3) on the Euclidean plane E^2
    sage: g.riemann().display()
    Riem(g) = 0

Polar coordinates are introduced by::

    sage: X_polar.<r,ph> = E.polar_coordinates()
    sage: X_polar
    Chart (E^2, (r, ph))

``E`` is now endowed with two coordinate charts::

    sage: E.atlas()
    [Chart (E^2, (x, y)), Chart (E^2, (r, ph))]

The ranges of the coordinates introduced so far are::

    sage: X_cartesian.coord_range()
    x: (-oo, +oo); y: (-oo, +oo)
    sage: X_polar.coord_range()
    r: (0, +oo); ph: (0, 2*pi)

The transition map from polar coordinates to Cartesian ones is::

    sage: E.coord_change(X_polar, X_cartesian).display()
    x = r*cos(ph)
    y = r*sin(ph)

while the reverse is::

    sage: E.coord_change(X_cartesian, X_polar).display()
    r = sqrt(x^2 + y^2)
    ph = arctan2(y, x)

At this stage, ``E`` is endowed with three vector frames::

    sage: E.frames()
    [Coordinate frame (E^2, (e_x,e_y)),
     Coordinate frame (E^2, (d/dr,d/dph)),
     Vector frame (E^2, (e_r,e_ph))]

The third one is the standard orthonormal frame associated with polar
coordinates, as we can check from the metric coefficients in it::

    sage: F_polar = E.polar_frame(); F_polar
    Vector frame (E^2, (e_r,e_ph))
    sage: g[F_polar,:]
    [1 0]
    [0 1]

The expression of the metric tensor in terms of the polar coordinates is::

    sage: g.display(X_polar.frame(), X_polar)
    g = dr*dr + r^2 dph*dph

A vector field on ``E``::

    sage: v = E.vector_field(-y, x, name='v'); v
    Vector field v on the Euclidean plane E^2
    sage: v.display()
    v = -y e_x + x e_y
    sage: v[:]
    [-y, x]

By default, the components of ``v``, as returned by ``display`` or the bracket
operator, refer to the Cartesian frame on ``E``; to get the components with
respect to the orthonormal polar frame, one has to specify it explicitly,
generally along with the polar chart for the coordinate expression of the
components::

    sage: v.display(F_polar, X_polar)
    v = r e_ph
    sage: v[F_polar,:,X_polar]
    [0, r]

A scalar field on ``E``::

    sage: f = E.scalar_field(x*y, name='f'); f
    Scalar field f on the Euclidean plane E^2
    sage: f.display()
    f: E^2 --> R
       (x, y) |--> x*y
       (r, ph) |--> r^2*cos(ph)*sin(ph)

The gradient of ``f``::

    sage: w = grad(f); w
    Vector field grad(f) on the Euclidean plane E^2
    sage: w.display()
    grad(f) = y e_x + x e_y
    sage: w.display(F_polar, X_polar)
    grad(f) = 2*r*cos(ph)*sin(ph) e_r + (2*cos(ph)^2 - 1)*r e_ph

The dot product of two vector fields::

    sage: s = v.dot(w); s
    Scalar field v.grad(f) on the Euclidean plane E^2
    sage: s.display()
    v.grad(f): E^2 --> R
       (x, y) |--> x^2 - y^2
       (r, ph) |--> (2*cos(ph)^2 - 1)*r^2
    sage: s.expr()
    x^2 - y^2

The norm is related to the dot product by the standard formula::

    sage: norm(v)^2 == v.dot(v)
    True

The divergence of the vector field ``v``::

    sage: s = div(v); s
    Scalar field div(v) on the Euclidean plane E^2
    sage: s.display()
    div(v): E^2 --> R
       (x, y) |--> 0
       (r, ph) |--> 0


AUTHORS:

- Eric Gourgoulhon (2018): initial version

"""

#*****************************************************************************
#       Copyright (C) 2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.functions.trig import cos, sin, atan2
from sage.functions.other import sqrt
from sage.misc.latex import latex
from sage.manifolds.differentiable.pseudo_riemannian import \
                                                       PseudoRiemannianManifold

###############################################################################

class EuclideanSpaceGeneric(PseudoRiemannianManifold):
    r"""
    Euclidean space.

    INPUT:

    - ``n`` -- positive integer; dimension of the space over the real field
    - ``name`` -- (default: ``None``) string; name (symbol) given to the
      Euclidean space; if ``None``, the name will be set to ``'E^n'``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the space; if ``None``, it is set to ``'\mathbb{E}^{n}'`` if
      ``name`` is  ``None`` and to ``name`` otherwise
    - ``coordinates`` -- (default: ``'Cartesian'``) type of coordinates to
      be initialized at the Euclidean space creation
    - ``symbols`` -- (default: ``None``) string defining the coordinate text
      symbols and LaTeX symbols, with the same conventions as the argument
      ``coordinates`` in
      :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; if ``None``,
      the symbols will be automatically generated
    - ``metric_name`` -- (default: ``'g'``) string; name (symbol) given to the
      Euclidean metric tensor
    - ``metric_latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the Euclidean metric tensor; if none is provided, it is set to
      ``metric_name``
    - ``start_index`` -- (default: 1) integer; lower value of the range of
      indices used for "indexed objects" in the Euclidean space, e.g.
      coordinates of a chart
    - ``ambient`` -- (default: ``None``) if not ``None``, must be an
      Euclidean space; the created object is then an open subset of ``ambient``
    - ``category`` -- (default: ``None``) to specify the category; if ``None``,
      ``Manifolds(RR).Differentiable()`` (or ``Manifolds(RR).Smooth()``
      if ``diff_degree`` = ``infinity``) is assumed (see the category
      :class:`~sage.categories.manifolds.Manifolds`)
    - ``names`` -- (default: ``None``) unused argument, except if
      ``symbols`` is not provided; it must then be a tuple containing
      the coordinate symbols (this is guaranteed if the shortcut operator
      ``<,>`` is used)
    - ``init_coord_methods`` -- (default: ``None``) dictionary of methods
      to initialize the various type of coordinates, with each key being a
      string describing the type of coordinates; to be used by derived classes
      only
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`
      would return the previously constructed object corresponding to these
      arguments)

    """
    def __init__(self, n, name=None, latex_name=None,
                 coordinates='Cartesian', symbols=None, metric_name='g',
                 metric_latex_name=None, start_index=1, ambient=None,
                 category=None, names=None, init_coord_methods=None,
                 unique_tag=None):
        r"""
        Construct an Euclidean space.

        TESTS::

            sage: E = EuclideanSpace(4); E
            4-dimensional Euclidean space E^4
            sage: E.metric()
            Riemannian metric g on the 4-dimensional Euclidean space E^4
            sage: TestSuite(E).run()

        """
        if name is None:
            name = 'E^{}'.format(n)
            if latex_name is None:
                latex_name = r'\mathbb{E}^{' + str(n) + '}'
        PseudoRiemannianManifold.__init__(self, n, name, metric_name=metric_name,
                                          signature=n, ambient=ambient,
                                          latex_name=latex_name,
                                          metric_latex_name=metric_latex_name,
                                          start_index=start_index,
                                          category=category)
        if names is not None and symbols is None:
            symbols = ''
            for x in names:
                symbols += x + ' '
            symbols = symbols[:-1]
        if symbols is None:
            if n == 1:
                if coordinates == 'Cartesian':
                    symbols = 'x'
                else:
                    raise TypeError("unkown coordinate type")
            elif n > 3:
                if coordinates == 'Cartesian':
                    symbols = ''
                    for i in self.irange():
                        symbols += "x{}".format(i) + r":x_{" + str(i) + r"} "
                    symbols = symbols[:-1]
                else:
                    raise TypeError("unkown coordinate type")
            else:
                raise NotImplementedError("dimension not implemented yet")
        self._cartesian_chart = None  # to be constructed later if necessary
        if init_coord_methods is None:
            self._init_coordinates = {'Cartesian':
                                      self._init_coordinates_cartesian}
        else:
            self._init_coordinates = init_coord_methods
        self._init_coordinates[coordinates](symbols)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: E = EuclideanSpace(4)
            sage: E._repr_()
            '4-dimensional Euclidean space E^4'
            sage: E  # indirect doctest
            4-dimensional Euclidean space E^4

        """
        return "{}-dimensional Euclidean space {}".format(self._dim,
                                                          self._name)

    def _first_ngens(self, n):
        r"""
        Return the list of coordinates of the default chart.

        This is useful only for the use of Sage preparser::

            sage: preparse("E.<x,y,z> = EuclideanSpace(3)")
            "E = EuclideanSpace(Integer(3), names=('x', 'y', 'z',)); (x, y, z,) = E._first_ngens(3)"

        TESTS::

            sage: E = EuclideanSpace(2)
            sage: E._first_ngens(2)
            (x, y)
            sage: E.<u,v> = EuclideanSpace(2)
            sage: E._first_ngens(2)
            (u, v)

        """
        return self._def_chart[:]

    def _init_coordinates_cartesian(self, symbols):
        r"""
        Construct the chart of Cartesian coordinates and initialize the
        components of the metric tensor in it.

        TESTS::

            sage: E = EuclideanSpace(2)
            sage: E._init_coordinates_cartesian('x y')

        """
        chart = self.chart(coordinates=symbols)
        self._cartesian_chart = chart
        frame = chart.frame()
        # Renaming (d/dx, d/dy, ...) to (e_x, e_y, ...):
        coords = chart[:]
        frame.set_name('e',
                       indices=tuple(str(x) for x in coords),
                       latex_indices=tuple(latex(x) for x in coords))
        g = self.metric()
        gc = g.add_comp(frame)
        for i in self.irange():
            gc[i, i, chart] = 1
        nabla = g.connection(init_coef=False)  # False to avoid any computation
        nabla.add_coef(frame)  # initialize a zero set of coefficients

    def cartesian_coordinates(self, symbols=None, names=None):
        r"""
        Return the chart of Cartesian coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as
          the argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the Cartesian chart has not been already defined; if
          ``None`` the symbols are generated as `(x_1,\ldots,x_n)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of Cartesian coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(4)
            sage: X = E.cartesian_coordinates(); X
            Chart (E^4, (x1, x2, x3, x4))
            sage: X[2]
            x2
            sage: X[:]
            (x1, x2, x3, x4)

        An example where the Cartesian coordinates have not been previously
        created::

            sage: E.<r,phi> = EuclideanSpace(2, coordinates='polar')
            sage: E.atlas()  # only polar coordinates exist
            [Chart (E^2, (r, phi))]
            sage: E.cartesian_coordinates(symbols='u v')
            Chart (E^2, (u, v))
            sage: E.atlas()  # the Cartesian chart has been added to the atlas
            [Chart (E^2, (r, phi)), Chart (E^2, (u, v))]

        Note that if the Cartesian coordinates have been already initialized,
        the argument ``symbols`` has not effect::

            sage: E.cartesian_coordinates(symbols='X Y')
            Chart (E^2, (u, v))

        """
        if self._cartesian_chart is None:
            if symbols is None:
                symbols = ''
                if names is None:
                    for i in self.irange():
                        symbols += "x{}".format(i) + r":x_{" + str(i) + r"} "
                else:
                    for x in names:
                        symbols += x + ' '
                symbols = symbols[:-1]
            self._init_coordinates_cartesian(symbols)
        return self._cartesian_chart

    def cartesian_frame(self):
        r"""
        Return the orthonormal vector frame associated with Cartesian
        coordinates.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.vectorframe.CoordFrame`

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.cartesian_frame()
            Coordinate frame (E^2, (e_x,e_y))
            sage: E.cartesian_frame()[1]
            Vector field e_x on the Euclidean plane E^2
            sage: E.cartesian_frame()[:]
            (Vector field e_x on the Euclidean plane E^2,
             Vector field e_y on the Euclidean plane E^2)

        """
        if self._cartesian_chart is None:
            self.cartesian_coordinates()  # creates the Cartesian chart
        return self._cartesian_chart.frame()

    def vector_field(self, *args, **kwargs):
        r"""
        Define a vector field on ``self``.

        INPUT:

        - ``args`` -- components of the vector field with respect to the
          vector frame specified by the argument ``frame`` or a dictionary
          of components (see examples below)
        - ``frame`` -- (default: ``None``) vector frame in which the components
          are given; if ``None``, the default vector frame on ``self`` is
          assumed
        - ``chart`` -- (default: ``None``) coordinate chart in which the
          components are expressed; if ``None``, the default chart on ``self``
          is assumed
        - ``name`` -- (default: ``None``) name given to the vector field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          vector field; if none is provided, the LaTeX symbol is set to
          ``name``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that `N = M`
          and that `\Phi` is the identity map (case of a vector field
          *on* `M`), otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.vectorfield.VectorFieldParal`)
          representing the defined vector field

        EXAMPLES:

        A vector field in the Euclidean plane::

            sage: E.<x,y> = EuclideanSpace(2)
            sage: v = E.vector_field(x*y, x+y)
            sage: v.display()
            x*y e_x + (x + y) e_y
            sage: v[:]
            [x*y, x + y]

        A name can be provided; it is then used for displaying the vector
        field::

            sage: v = E.vector_field(x*y, x+y, name='v')
            sage: v.display()
            v = x*y e_x + (x + y) e_y

        It is also possible to initialize a vector field from a vector of
        symbolic expressions::

            sage: v = E.vector_field(vector([x*y, x+y]))
            sage: v.display()
            x*y e_x + (x + y) e_y

        If the components are relative to a vector frame different from the
        default one (here the Cartesian frame `(e_x,e_y)`), the vector
        frame has to be specified explicitly::

            sage: F_polar = E.polar_frame(); F_polar
            Vector frame (E^2, (e_r,e_ph))
            sage: v = E.vector_field(1, 0, frame=F_polar)
            sage: v.display(F_polar)
            e_r
            sage: v.display()
            x/sqrt(x^2 + y^2) e_x + y/sqrt(x^2 + y^2) e_y

        The argument ``chart`` must be used to specify in which coordinate
        chart the components are expressed::

            sage: X_polar.<r, ph> = E.polar_coordinates()
            sage: v = E.vector_field(0, r, frame=F_polar, chart=X_polar)
            sage: v.display(F_polar, X_polar)
            r e_ph
            sage: v.display()
            -y e_x + x e_y

        It is also possible to pass the components as a dictionary, with
        a pair (vector frame, chart) as a key::

            sage: v = E.vector_field({(F_polar, X_polar): (0, r)})
            sage: v.display(F_polar, X_polar)
            r e_ph

        The key can be reduced to the vector frame if the chart is the default
        one::

            sage: v = E.vector_field({F_polar: (0, 1)})
            sage: v.display(F_polar)
            e_ph

        Finally, it is possible to construct the vector field without
        initializing any component::

            sage: v = E.vector_field(); v
            Vector field on the Euclidean plane E^2

        The components can then by set in a second stage, via the square
        bracket operator, the unset components being assumed to be zero::

            sage: v[1] = x*y
            sage: v.display()  # v[2] is zero
            x*y e_x
            sage: v[2] = x+y
            sage: v.display()
            x*y e_x + (x + y) e_y

        The above is equivalent to::

            sage: v[:] = x*y, x+y
            sage: v.display()
            x*y e_x + (x + y) e_y

        The square bracket operator can also be used to set components in a
        vector frame that is not the default one::

            sage: v = E.vector_field(name='v')
            sage: v[F_polar, 1, X_polar] = r
            sage: v.display(F_polar, X_polar)
            v = r e_r
            sage: v.display()
            v = x e_x + y e_y

        """
        name = kwargs.get('name')
        latex_name = kwargs.get('latex_name')
        dest_map = kwargs.get('dest_map')
        resu = super(EuclideanSpaceGeneric, self).vector_field(name=name,
                                      latex_name=latex_name, dest_map=dest_map)
        if args:
            # Some components are to be initialized
            args0 = args[0]
            if isinstance(args0, dict):
                for frame, components in args0.items():
                    chart = None
                    if isinstance(frame, tuple):
                        # a pair (frame, chart) has been provided:
                        frame, chart = frame
                    resu.add_comp(frame)[:, chart] = components
            else:
                if hasattr(args0, '__getitem__') and len(args0) == self._dim:
                    # args0 is a list/vector of components
                    args = args0
                elif len(args) != self._dim:
                    raise ValueError("{} components must ".format(self._dim) +
                                     "be provided")
                    # args is the tuple of components in a specific frame
                frame = kwargs.get('frame')
                chart = kwargs.get('chart')
                resu.add_comp(frame)[:, chart] = args
        return resu


###############################################################################

class EuclideanPlane(EuclideanSpaceGeneric):
    r"""
    Euclidean plane.

    INPUT:

    - ``name`` -- (default: ``None``) string; name (symbol) given to the
      Euclidean plane; if ``None``, the name will be set to ``'E^2'``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      Euclidean plane; if ``None``, it is set to ``'\mathbb{E}^{2}'`` if
      ``name`` is  ``None`` and to ``name`` otherwise
    - ``coordinates`` -- (default: ``'Cartesian'``) type of coordinates to
      be initialized at the Euclidean plane creation
    - ``symbols`` -- (default: ``None``) string defining the coordinate text
      symbols and LaTeX symbols, with the same conventions as the argument
      ``coordinates`` in
      :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; if ``None``,
      the symbols will be automatically generated
    - ``metric_name`` -- (default: ``'g'``) string; name (symbol) given to the
      Euclidean metric tensor
    - ``metric_latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the Euclidean metric tensor; if none is provided, it is set to
      ``metric_name``
    - ``start_index`` -- (default: 1) integer; lower value of the range of
      indices used for "indexed objects" in the Euclidean plane, e.g.
      coordinates of a chart
    - ``ambient`` -- (default: ``None``) if not ``None``, must be an
      Euclidean plane; the created object is then an open subset of ``ambient``
    - ``category`` -- (default: ``None``) to specify the category; if ``None``,
      ``Manifolds(RR).Differentiable()`` (or ``Manifolds(RR).Smooth()``
      if ``diff_degree`` = ``infinity``) is assumed (see the category
      :class:`~sage.categories.manifolds.Manifolds`)
    - ``names`` -- (default: ``None``) unused argument, except if
      ``symbols`` is not provided; it must then be a tuple containing
      the coordinate symbols (this is guaranteed if the shortcut operator
      ``<,>`` is used)
    - ``init_coord_methods`` -- (default: ``None``) dictionary of methods
      to initialize the various type of coordinates, with each key being a
      string describing the type of coordinates; to be used by derived classes
      only
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`
      would return the previously constructed object corresponding to these
      arguments)

    """
    def __init__(self, name=None, latex_name=None, coordinates='Cartesian',
                 symbols=None, metric_name='g', metric_latex_name=None,
                 start_index=1, ambient=None, category=None, names=None,
                 unique_tag=None):
        r"""
        Construct an Euclidean plane.

        TESTS::

            sage: E = EuclideanSpace(2); E
            Euclidean plane E^2
            sage: E.metric()
            Riemannian metric g on the Euclidean plane E^2
            sage: TestSuite(E).run()

        """
        if names is not None and symbols is None:
            symbols = ''
            for x in names:
                symbols += x + ' '
            symbols = symbols[:-1]
            if coordinates == 'polar':
                if names[1] in ['p', 'ph', 'phi']:
                    symbols += ':\\phi'
                elif names[1] in ['t', 'th', 'theta']:
                    symbols += ':\\theta'
        if symbols is None:
            if coordinates == 'Cartesian':
                symbols = 'x y'
            elif coordinates == 'polar':
                symbols = 'r ph:\\phi'
            else:
                raise TypeError("unkown coordinate type")
        self._polar_chart = None   # to be constructed later if necessary
        self._polar_frame = None   #
        init_coord_methods = {'Cartesian': self._init_coordinates_cartesian,
                              'polar': self._init_coordinates_polar}
        EuclideanSpaceGeneric.__init__(self, 2, name=name,
                                       latex_name=latex_name,
                                       coordinates=coordinates,
                                       symbols=symbols,
                                       metric_name=metric_name,
                                       metric_latex_name=metric_latex_name,
                                       start_index=start_index,
                                       ambient=ambient, category=category,
                                       init_coord_methods=init_coord_methods)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: E = EuclideanSpace(2)
            sage: E._repr_()
            'Euclidean plane E^2'
            sage: E  # indirect doctest
            Euclidean plane E^2
            sage: E = EuclideanSpace(2, name='E')
            sage: E._repr_()
            'Euclidean plane E'

        """
        return "Euclidean plane {}".format(self._name)

    def _init_coordinates_polar(self, symbols):
        r"""
        Construct the chart of polar coordinates and initialize the
        components of the metric tensor in it.

        TESTS::

            sage: E = EuclideanSpace(2)
            sage: E.atlas()
            [Chart (E^2, (x, y))]
            sage: E._init_coordinates_polar(symbols=r"R Phi:\Phi")
            sage: E.atlas()
            [Chart (E^2, (x, y)), Chart (E^2, (R, Phi))]

        """
        coords = symbols.split()  # list of strings, one per coordinate
        # Adding the coordinate ranges:
        coordinates = coords[0] + ':(0,+oo) ' + coords[1] + ':(0,2*pi)'
        chart = self.chart(coordinates=coordinates)
        self._polar_chart = chart
        frame = chart.frame()
        # Initialization of the metric components and the associated
        # Christoffel symbols
        g = self.metric()
        gc = g.add_comp(frame)
        i0 = self._sindex
        i1 = i0 + 1
        r, ph = chart[:]
        gc[i0, i0, chart] = 1
        gc[i1, i1, chart] = r**2
        nabla = g.connection()
        nabla.coef(frame)
        # Orthonormal frame associated with polar coordinates:
        to_orthonormal = self.automorphism_field()
        to_orthonormal[frame, i0, i0, chart] = 1
        to_orthonormal[frame, i1, i1, chart] = 1/r
        oframe = frame.new_frame(to_orthonormal, 'e',
                                 indices=(str(r), str(ph)),
                                 latex_indices=(latex(r), latex(ph)))
        self._polar_frame = oframe
        g.comp(oframe)
        nabla.coef(oframe)

    def _transition_polar_cartesian(self):
        r"""
        Transitions between polar and Cartesian coordinates.

        TESTS::

            sage: E = EuclideanSpace(2)
            sage: E._init_coordinates_polar(symbols=r"r ph:\phi")
            sage: E.atlas()
            [Chart (E^2, (x, y)), Chart (E^2, (r, ph))]
            sage: E.coord_changes() # no transition map has been defined yet
            {}
            sage: E._transition_polar_cartesian()
            sage: E.coord_changes()  # random (dictionary output)
            {(Chart (E^2, (r, ph)),
              Chart (E^2, (x, y))): Change of coordinates from
              Chart (E^2, (r, ph)) to Chart (E^2, (x, y)),
             (Chart (E^2, (x, y)),
              Chart (E^2, (r, ph))): Change of coordinates from
              Chart (E^2, (x, y)) to Chart (E^2, (r, ph))}

            """
        # Transition maps polar chart <-> Cartesian chart
        chart_cart = self._cartesian_chart
        chart_pol = self._polar_chart
        x, y = chart_cart[:]
        r, ph = chart_pol[:]
        pol_to_cart = chart_pol.transition_map(chart_cart,
                                               [r*cos(ph), r*sin(ph)])
        pol_to_cart.set_inverse(sqrt(x**2+y**2), atan2(y,x))
        # Automorphisms ortho polar frame <-> Cartesian frame
        oframe = self._polar_frame
        cframe = chart_cart.frame()
        sframe = chart_pol.frame()
        changes = self._frame_changes
        cframe_to_oframe = changes[(sframe, oframe)] * \
                           changes[(cframe, sframe)]
        oframe_to_cframe = changes[(sframe, cframe)] * \
                           changes[(oframe, sframe)]
        changes[(cframe, oframe)] = cframe_to_oframe
        changes[(oframe, cframe)] = oframe_to_cframe
        vmodule = self.vector_field_module()
        vmodule._basis_changes[(cframe, oframe)] = cframe_to_oframe
        vmodule._basis_changes[(oframe, cframe)] = oframe_to_cframe

    def cartesian_coordinates(self, symbols=None, names=None):
        r"""
        Return the chart of Cartesian coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as the
          argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the Cartesian chart has not been already defined; if
          ``None`` the symbols are generated as `(x,y)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of Cartesian coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.cartesian_coordinates()
            Chart (E^2, (x, y))

        An example where the Cartesian coordinates have not been previously
        created::

            sage: E = EuclideanSpace(2, coordinates='polar')
            sage: E.atlas()  # only polar coordinates exist
            [Chart (E^2, (r, ph))]
            sage: E.cartesian_coordinates(symbols='X Y')
            Chart (E^2, (X, Y))
            sage: E.atlas()  # the Cartesian chart has been added to the atlas
            [Chart (E^2, (r, ph)), Chart (E^2, (X, Y))]

        The coordinates themselves are obtained via the square bracket
        operator::

            sage: E.cartesian_coordinates()[1]
            X
            sage: E.cartesian_coordinates()[2]
            Y
            sage: E.cartesian_coordinates()[:]
            (X, Y)

        It is also possible to use the operator ``<,>`` to set symbolic
        variable containing the coordinates::

            sage: E = EuclideanSpace(2, coordinates='polar')
            sage: cartesian_coord.<u,v> = E.cartesian_coordinates()
            sage: cartesian_coord
            Chart (E^2, (u, v))
            sage: u,v
            (u, v)

        The command ``cartesian_coord.<u,v> = E.cartesian_coordinates()``
        is actually a shortcut for::

            sage: cartesian_coord = E.cartesian_coordinates(symbols='u v')
            sage: u, v = cartesian_coord[:]

        """
        if self._cartesian_chart is None:
            if symbols is None:
                if names is None:
                    symbols = 'x y'
                else:
                    symbols = names[0] + ' ' + names[1]
            self._init_coordinates_cartesian(symbols)
            if self._polar_chart:
                self._transition_polar_cartesian()
        return self._cartesian_chart

    def polar_coordinates(self, symbols=None, names=None):
        r"""
        Return the chart of polar coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as the
          argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the polar chart has not been already defined; if
          ``None`` the symbols are generated as `(r,\phi)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of polar coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.polar_coordinates()
            Chart (E^2, (r, ph))
            sage: latex(E.polar_coordinates())
            \left(\mathbb{E}^{2},(r, {\phi})\right)

        The coordinates can be obtained via the square bracket operator::

            sage: E.polar_coordinates()[1]
            r
            sage: E.polar_coordinates()[2]
            ph
            sage: E.polar_coordinates()[:]
            (r, ph)

        They can also be obtained via the operator ``<,>``::

            sage: E = EuclideanSpace(2)
            sage: polar.<r,ph> = E.polar_coordinates(); polar
            Chart (E^2, (r, ph))
            sage: r, ph
            (r, ph)

        Actually, ``polar.<r,ph> = E.polar_coordinates()`` is a shortcut for::

            sage: polar = E.polar_coordinates()
            sage: r, ph = polar[:]

        The coordinate symbols can be customized::

            sage: E = EuclideanSpace(2)
            sage: E.polar_coordinates(symbols=r"r th:\theta")
            Chart (E^2, (r, th))
            sage: latex(E.polar_coordinates())
            \left(\mathbb{E}^{2},(r, {\theta})\right)

        Note that if the polar coordinates have been already initialized, the
        argument ``symbols`` has no effect::

            sage: E.polar_coordinates(symbols=r"R Th:\Theta")
            Chart (E^2, (r, th))

        """
        if self._polar_chart is None:
            if symbols is None:
                if names is None:
                    symbols = 'r ph:\\phi'
                else:
                    symbols = names[0] + ' ' + names[1]
                    if names[1] in ['p', 'ph', 'phi']:
                        symbols += ':\\phi'
                    elif names[1] in ['t', 'th', 'theta']:
                        symbols += ':\\theta'
            self._init_coordinates_polar(symbols)
            if self._cartesian_chart:
                self._transition_polar_cartesian()
        return self._polar_chart

    def polar_frame(self):
        r"""
        Return the orthonormal vector frame associated with polar
        coordinates.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.polar_frame()
            Vector frame (E^2, (e_r,e_ph))
            sage: E.polar_frame()[1]
            Vector field e_r on the Euclidean plane E^2
            sage: E.polar_frame()[:]
            (Vector field e_r on the Euclidean plane E^2,
             Vector field e_ph on the Euclidean plane E^2)

        """
        if self._polar_frame is None:
            self.polar_coordinates()  # creates the polar chart
        return self._polar_frame

###############################################################################

def EuclideanSpace(n, name=None, latex_name=None, coordinates='Cartesian',
                   symbols=None, metric_name='g', metric_latex_name=None,
                   start_index=1, names=None):
    r"""
    Construct an Euclidean space.

    INPUT:

        - ``n`` -- positive integer; dimension of the space over the real field
        - ``name`` -- (default: ``None``) string; name (symbol) given to the
          Euclidean space; if ``None``, the name will be set to ``'E^n'``
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
          denote the Euclidean space; if ``None``, it is set to
          ``'\mathbb{E}^{n}'`` if ``name`` is  ``None`` and to ``name``
          otherwise
        - ``coordinates`` -- (default: ``'Cartesian'``) type of coordinates to
          be initialized at the Euclidean space creation
        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as the
          argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`;
          if ``None``, the symbols will be automatically generated
        - ``metric_name`` -- (default: ``'g'``) string; name (symbol) given to
          the Euclidean metric tensor
        - ``metric_latex_name`` -- (default: ``None``) string; LaTeX symbol to
          denote the Euclidean metric tensor; if none is provided, it is set to
          ``metric_name``
        - ``start_index`` -- (default: 1) integer; lower value of the range of
          indices used for "indexed objects" in the Euclidean space, e.g.
          coordinates of a chart
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must then be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

    """
    from sage.misc.prandom import getrandbits
    from time import time
    if n == 2:
        return EuclideanPlane(name=name, latex_name=latex_name,
                              coordinates=coordinates, symbols=symbols,
                              metric_name=metric_name,
                              metric_latex_name=metric_latex_name,
                              start_index=start_index, names=names,
                              unique_tag=getrandbits(128)*time())
    return EuclideanSpaceGeneric(n, name=name, latex_name=latex_name,
                              coordinates=coordinates, symbols=symbols,
                              metric_name=metric_name,
                              metric_latex_name=metric_latex_name,
                              start_index=start_index, names=names,
                              unique_tag=getrandbits(128)*time())

