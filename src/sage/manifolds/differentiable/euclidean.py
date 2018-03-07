r"""
Euclidean spaces

An *Euclidean space of dimension* `n` is a Riemannian manifold diffeomorphic to
`\RR^n` equipped with a flat metric.

Euclidean spaces are implemented via the following classes

- :class:`EuclideanSpaceGeneric`

  - :class:`EuclideanPlane` (case `n=2`)
  - :class:`Euclidean3dimSpace` (case `n=3`)

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
    sage: cartesian = E.default_chart(); cartesian
    Chart (E^2, (x, y))

Thanks to the use of ``<x,y>`` when declaring ``E``, the coordinates `(x,y)`
have been injected in the global namespace, i.e. the Python variables ``x``
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

Polar coordinates `(r,\phi)` are introduced by::

    sage: polar.<r,ph> = E.polar_coordinates()
    sage: polar
    Chart (E^2, (r, ph))

``E`` is now endowed with two coordinate charts::

    sage: E.atlas()
    [Chart (E^2, (x, y)), Chart (E^2, (r, ph))]

The ranges of the coordinates introduced so far are::

    sage: cartesian.coord_range()
    x: (-oo, +oo); y: (-oo, +oo)
    sage: polar.coord_range()
    r: (0, +oo); ph: (0, 2*pi)

The transition map from polar coordinates to Cartesian ones is::

    sage: E.coord_change(polar, cartesian).display()
    x = r*cos(ph)
    y = r*sin(ph)

while the reverse is::

    sage: E.coord_change(cartesian, polar).display()
    r = sqrt(x^2 + y^2)
    ph = arctan2(y, x)

At this stage, ``E`` is endowed with three vector frames::

    sage: E.frames()
    [Coordinate frame (E^2, (e_x,e_y)),
     Coordinate frame (E^2, (d/dr,d/dph)),
     Vector frame (E^2, (e_r,e_ph))]

The third one is the standard orthonormal frame associated with polar
coordinates, as we can check from the metric coefficients in it::

    sage: polar_frame = E.polar_frame(); polar_frame
    Vector frame (E^2, (e_r,e_ph))
    sage: g[polar_frame,:]
    [1 0]
    [0 1]

The expression of the metric tensor in terms of polar coordinates is::

    sage: g.display(polar.frame(), polar)
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

    sage: v.display(polar_frame, polar)
    v = r e_ph
    sage: v[polar_frame,:,polar]
    [0, r]

A scalar field on ``E``::

    sage: f = E.scalar_field(x*y, name='f'); f
    Scalar field f on the Euclidean plane E^2
    sage: f.display()
    f: E^2 --> R
       (x, y) |--> x*y
       (r, ph) |--> r^2*cos(ph)*sin(ph)

The gradient of ``f``::

    sage: from sage.manifolds.operators import * # to get grad, div, etc.
    sage: w = grad(f); w
    Vector field grad(f) on the Euclidean plane E^2
    sage: w.display()
    grad(f) = y e_x + x e_y
    sage: w.display(polar_frame, polar)
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

.. RUBRIC:: Example 2: Vector calculus in the Euclidean 3-space

We start by declaring the 3-dimensional Euclidean space ``E``, with
``(x,y,z)`` as Cartesian coordinates::

    sage: E.<x,y,z> = EuclideanSpace(3)
    sage: E
    Euclidean space E^3

A simple vector field on ``E``::

    sage: v = E.vector_field(-y, x, 0, name='v')
    sage: v.display()
    v = -y e_x + x e_y
    sage: v[:]
    [-y, x, 0]

The Euclidean norm of ``v``::

    sage: s = norm(v); s
    Scalar field |v| on the Euclidean space E^3
    sage: s.display()
    |v|: E^3 --> R
       (x, y, z) |--> sqrt(x^2 + y^2)
    sage: s.expr()
    sqrt(x^2 + y^2)

The divergence of ``v`` is zero::

    sage: from sage.manifolds.operators import *
    sage: div(v)
    Scalar field div(v) on the Euclidean space E^3
    sage: div(v).display()
    div(v): E^3 --> R
       (x, y, z) |--> 0

while its curl is a constant vector field along `e_z`::

    sage: w = curl(v); w
    Vector field curl(v) on the Euclidean space E^3
    sage: w.display()
    curl(v) = 2 e_z

The gradient of a scalar field::

    sage: f = E.scalar_field(sin(x*y*z), name='f')
    sage: u = grad(f); u
    Vector field grad(f) on the Euclidean space E^3
    sage: u.display()
    grad(f) = y*z*cos(x*y*z) e_x + x*z*cos(x*y*z) e_y + x*y*cos(x*y*z) e_z

The curl of a gradient is zero::

    sage: curl(u).display()
    curl(grad(f)) = 0

The dot product of two vector fields::

    sage: s = u.dot(v); s
    Scalar field grad(f).v on the Euclidean space E^3
    sage: s.expr()
    (x^2 - y^2)*z*cos(x*y*z)

The cross product of two vector fields::

    sage: a = u.cross(v); a
    Vector field grad(f) x v on the Euclidean space E^3
    sage: a.display()
    grad(f) x v = -x^2*y*cos(x*y*z) e_x - x*y^2*cos(x*y*z) e_y
     + 2*x*y*z*cos(x*y*z) e_z

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
      the symbols will be automatically generated, depending on the value
      of ``coordinates``
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

    EXAMPLES:

    A 4-dimensional Euclidean space::

        sage: E = EuclideanSpace(4); E
        4-dimensional Euclidean space E^4
        sage: latex(E)
        \mathbb{E}^{4}

    ``E`` belongs to the class :class:`EuclideanSpaceGeneric` (actually to
    a subclass of it via SageMath's category mechanism)::

        sage: type(E)
        <class 'sage.manifolds.differentiable.euclidean.EuclideanSpaceGeneric_with_category'>

    ``E`` is a real smooth manifold of dimension 4::

        sage: E.category()
        Category of smooth manifolds over Real Field with 53 bits of precision
        sage: dim(E)
        4

    It is endowed with a default coordinate chart, which is that of Cartesian
    coordinates `(x_1,x_2,x_3,x_3)`::

        sage: E.atlas()
        [Chart (E^4, (x1, x2, x3, x4))]
        sage: E.default_chart()
        Chart (E^4, (x1, x2, x3, x4))
        sage: E.default_chart() is E.cartesian_coordinates()
        True

    ``E`` is also endowed with a default metric tensor::

        sage: g = E.metric(); g
        Riemannian metric g on the 4-dimensional Euclidean space E^4
        sage: g.display()
        g = dx1*dx1 + dx2*dx2 + dx3*dx3 + dx4*dx4

    This makes ``E`` be a pseudo-Riemannian manifold. Actually the class
    :class:`EuclideanSpaceGeneric` implementing ``E`` inherits from the class
    :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`::

        sage: isinstance(E, sage.manifolds.differentiable.pseudo_riemannian.
        ....:            PseudoRiemannianManifold)
        True

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
            gc[[i, i]] = self.one_scalar_field()
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
            sage: latex(X[:])
            \left({x_{1}}, {x_{2}}, {x_{3}}, {x_{4}}\right)

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

            sage: polar_frame = E.polar_frame(); polar_frame
            Vector frame (E^2, (e_r,e_ph))
            sage: v = E.vector_field(1, 0, frame=polar_frame)
            sage: v.display(polar_frame)
            e_r
            sage: v.display()
            x/sqrt(x^2 + y^2) e_x + y/sqrt(x^2 + y^2) e_y

        The argument ``chart`` must be used to specify in which coordinate
        chart the components are expressed::

            sage: polar.<r, ph> = E.polar_coordinates()
            sage: v = E.vector_field(0, r, frame=polar_frame, chart=polar)
            sage: v.display(polar_frame, polar)
            r e_ph
            sage: v.display()
            -y e_x + x e_y

        It is also possible to pass the components as a dictionary, with
        a pair (vector frame, chart) as a key::

            sage: v = E.vector_field({(polar_frame, polar): (0, r)})
            sage: v.display(polar_frame, polar)
            r e_ph

        The key can be reduced to the vector frame if the chart is the default
        one::

            sage: v = E.vector_field({polar_frame: (0, 1)})
            sage: v.display(polar_frame)
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
            sage: v[polar_frame, 1, polar] = r
            sage: v.display(polar_frame, polar)
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
        if coordinates not in ['Cartesian', 'polar']:
            raise TypeError("unkown coordinate type")
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
        self._polar_chart = None  # to be constructed later if necessary
        self._polar_frame = None  # orthonormal frame associated to polar coord
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
        if coordinates == 'polar':
            # The default frame is the polar coordinate frame; we change it
            # to the orthornomal polar frame
            self.set_default_frame(self.polar_frame())

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
            sage: E._init_coordinates_polar(r"R Phi:\Phi")
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
        i1 = self._sindex
        i2 = i1 + 1
        r, ph = chart[:]
        gc[i1, i1, chart] = 1
        gc[i2, i2, chart] = r**2
        # Orthonormal frame associated with polar coordinates:
        to_orthonormal = self.automorphism_field()
        to_orthonormal[frame, i1, i1, chart] = 1
        to_orthonormal[frame, i2, i2, chart] = 1/r
        oframe = frame.new_frame(to_orthonormal, 'e',
                                 indices=(str(r), str(ph)),
                                 latex_indices=(latex(r), latex(ph)))
        self._polar_frame = oframe
        g.comp(oframe)

    def _transition_polar_cartesian(self):
        r"""
        Transitions between polar and Cartesian coordinates.

        TESTS::

            sage: E = EuclideanSpace(2)
            sage: E._init_coordinates_polar(r"r ph:\phi")
            sage: E.atlas()
            [Chart (E^2, (x, y)), Chart (E^2, (r, ph))]
            sage: E.coord_changes() # no transition map has been defined yet
            {}
            sage: E._transition_polar_cartesian()
            sage: len(E.coord_changes())  # 2 transition maps have been created
            2

        Tests of the change-of-frame formulas::

            sage: polar = E.polar_coordinates()
            sage: polar_f = E.polar_frame()
            sage: cart_f = E.cartesian_frame()
            sage: E.change_of_frame(cart_f, polar_f)[polar_f,:,polar]
            [ cos(ph) -sin(ph)]
            [ sin(ph)  cos(ph)]
            sage: E.change_of_frame(polar_f, cart_f)[polar_f,:,polar]
            [ cos(ph)  sin(ph)]
            [-sin(ph)  cos(ph)]

        """
        # Transition maps polar chart <-> Cartesian chart
        chart_cart = self._cartesian_chart
        chart_pol = self._polar_chart
        x, y = chart_cart[:]
        r, ph = chart_pol[:]
        pol_to_cart = chart_pol.transition_map(chart_cart,
                                               [r*cos(ph), r*sin(ph)])
        pol_to_cart.set_inverse(sqrt(x**2+y**2), atan2(y,x))
        # Automorphisms orthonormal polar frame <-> Cartesian frame
        oframe = self._polar_frame
        cframe = chart_cart.frame()
        sframe = chart_pol.frame()
        chg = self._frame_changes
        cframe_to_oframe = chg[(sframe, oframe)] * chg[(cframe, sframe)]
        oframe_to_cframe = chg[(sframe, cframe)] * chg[(oframe, sframe)]
        chg[(cframe, oframe)] = cframe_to_oframe
        chg[(oframe, cframe)] = oframe_to_cframe
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

        Note that if the Cartesian coordinates have been already initialized,
        the argument ``symbols`` has no effect::

            sage: E.cartesian_coordinates(symbols='x y')
            Chart (E^2, (X, Y))

        The coordinate variables are returned by the square bracket operator::

            sage: E.cartesian_coordinates()[1]
            X
            sage: E.cartesian_coordinates()[2]
            Y
            sage: E.cartesian_coordinates()[:]
            (X, Y)

        It is also possible to use the operator ``<,>`` to set symbolic
        variable containing the coordinates::

            sage: E = EuclideanSpace(2, coordinates='polar')
            sage: cartesian.<u,v> = E.cartesian_coordinates()
            sage: cartesian
            Chart (E^2, (u, v))
            sage: u,v
            (u, v)

        The command ``cartesian.<u,v> = E.cartesian_coordinates()`` is actually
        a shortcut for::

            sage: cartesian = E.cartesian_coordinates(symbols='u v')
            sage: u, v = cartesian[:]

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

        The relation to Cartesian coordinates is::

            sage: E.coord_change(E.polar_coordinates(),
            ....:                E.cartesian_coordinates()).display()
            x = r*cos(ph)
            y = r*sin(ph)

        The coordinate variables are returned by the square bracket operator::

            sage: E.polar_coordinates()[1]
            r
            sage: E.polar_coordinates()[2]
            ph
            sage: E.polar_coordinates()[:]
            (r, ph)

        They can also be obtained via the operator ``<,>``::

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
            self.polar_coordinates()  # creates the polar chart and the
                                      # associated orthonormal frame
        return self._polar_frame


###############################################################################

class Euclidean3dimSpace(EuclideanSpaceGeneric):
    r"""
    3-dimensional Euclidean space.

    INPUT:

    - ``name`` -- (default: ``None``) string; name (symbol) given to the
      Euclidean 3-space; if ``None``, the name will be set to ``'E^3'``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      Euclidean 3-space; if ``None``, it is set to ``'\mathbb{E}^{3}'`` if
      ``name`` is  ``None`` and to ``name`` otherwise
    - ``coordinates`` -- (default: ``'Cartesian'``) type of coordinates to
      be initialized at the Euclidean 3-space creation
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
      indices used for "indexed objects" in the Euclidean 3-space, e.g.
      coordinates of a chart
    - ``ambient`` -- (default: ``None``) if not ``None``, must be an
      Euclidean 3-space; the created object is then an open subset of
      ``ambient``
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
        Construct an Euclidean 3-space.

        TESTS::

            sage: E = EuclideanSpace(3); E
            Euclidean space E^3
            sage: E.metric()
            Riemannian metric g on the Euclidean space E^3
            sage: TestSuite(E).run()

        """
        if coordinates not in ['Cartesian', 'spherical', 'cylindrical']:
            raise TypeError("unkown coordinate type")
        if names is not None and symbols is None:
            # We add the LaTeX symbols when relevant:
            if coordinates == 'spherical':
                names = list(names)
                if names[1] in ['t', 'th', 'theta']:
                    names[1] = names[1] + ':\\theta'
                if names[2] in ['p', 'ph', 'phi']:
                    names[2] = names[2] + ':\\phi'
            elif coordinates == 'cylindrical':
                names = list(names)
                if names[0] in ['rh', 'rho']:
                    names[0] = names[0] + ':\\rho'
                if names[1] in ['p', 'ph', 'phi']:
                    names[1] = names[1] + ':\\phi'
            symbols = ''
            for x in names:
                symbols += x + ' '
            symbols = symbols[:-1]
        if symbols is None:
            if coordinates == 'Cartesian':
                symbols = 'x y z'
            elif coordinates == 'spherical':
                symbols = 'r th:\\theta ph:\\phi'
            elif coordinates == 'cylindrical':
                symbols = 'r ph:\\phi z'
        self._spherical_chart = None    # to be constructed later if necessary
        self._spherical_frame = None    # orthonormal frame
        self._cylindrical_chart = None  #
        self._cylindrical_frame = None  # orthonormal frame
        init_coord_methods = {'Cartesian': self._init_coordinates_cartesian,
                              'spherical': self._init_coordinates_spherical,
                              'cylindrical': self._init_coordinates_cylindrical}
        EuclideanSpaceGeneric.__init__(self, 3, name=name,
                                       latex_name=latex_name,
                                       coordinates=coordinates,
                                       symbols=symbols,
                                       metric_name=metric_name,
                                       metric_latex_name=metric_latex_name,
                                       start_index=start_index,
                                       ambient=ambient, category=category,
                                       init_coord_methods=init_coord_methods)
        if coordinates == 'spherical':
            # The default frame is the spherical coordinate frame; we change it
            # to the orthornomal spherical frame
            self.set_default_frame(self.spherical_frame())
        if coordinates == 'cylindrical':
            # The default frame is the cylindrical coordinate frame; we change
            # it to the orthornomal cylindrical frame
            self.set_default_frame(self.cylindrical_frame())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: E = EuclideanSpace(3)
            sage: E._repr_()
            'Euclidean space E^3'
            sage: E  # indirect doctest
            Euclidean space E^3
            sage: E = EuclideanSpace(3, name='E')
            sage: E._repr_()
            'Euclidean space E'

        """
        return "Euclidean space {}".format(self._name)

    def _init_coordinates_spherical(self, symbols):
        r"""
        Construct the chart of spherical coordinates and initialize the
        components of the metric tensor in it.

        TESTS::

            sage: E = EuclideanSpace(3)
            sage: E.atlas()
            [Chart (E^3, (x, y, z))]
            sage: E._init_coordinates_spherical(r"R Th:\Theta Ph:\Phi")
            sage: E.atlas()
            [Chart (E^3, (x, y, z)), Chart (E^3, (R, Th, Ph))]

        """
        coords = symbols.split()  # list of strings, one per coordinate
        # Adding the coordinate ranges:
        coordinates = coords[0] + ':(0,+oo) ' + coords[1] + ':(0,pi) ' + \
                      coords[2] + ':(0,2*pi)'
        chart = self.chart(coordinates=coordinates)
        self._spherical_chart = chart
        frame = chart.frame()
        # Initialization of the metric components and the associated
        # Christoffel symbols
        g = self.metric()
        gc = g.add_comp(frame)
        i1 = self._sindex
        i2 = i1 + 1
        i3 = i1 + 2
        r, th, ph = chart[:]
        gc[i1, i1, chart] = 1
        gc[i2, i2, chart] = r**2
        gc[i3, i3, chart] = (r*sin(th))**2
        # Orthonormal frame associated with spherical coordinates:
        to_orthonormal = self.automorphism_field()
        to_orthonormal[frame, i1, i1, chart] = 1
        to_orthonormal[frame, i2, i2, chart] = 1/r
        to_orthonormal[frame, i3, i3, chart] = 1/(r*sin(th))
        oframe = frame.new_frame(to_orthonormal, 'e',
                                 indices=(str(r), str(th), str(ph)),
                                 latex_indices=(latex(r), latex(th), latex(ph)))
        self._spherical_frame = oframe
        g.comp(oframe)

    def _init_coordinates_cylindrical(self, symbols):
        r"""
        Construct the chart of cylindrical coordinates and initialize the
        components of the metric tensor in it.

        TESTS::

            sage: E = EuclideanSpace(3)
            sage: E.atlas()
            [Chart (E^3, (x, y, z))]
            sage: E._init_coordinates_cylindrical(r"r ph:\phi z")
            sage: E.atlas()
            [Chart (E^3, (x, y, z)), Chart (E^3, (r, ph, z))]

        """
        coords = symbols.split()  # list of strings, one per coordinate
        # Adding the coordinate ranges:
        coordinates = coords[0] + ':(0,+oo) ' + coords[1] + ':(0,2*pi) ' +  \
                      coords[2]
        chart = self.chart(coordinates=coordinates)
        self._cylindrical_chart = chart
        frame = chart.frame()
        # Initialization of the metric components and the associated
        # Christoffel symbols
        g = self.metric()
        gc = g.add_comp(frame)
        i1 = self._sindex
        i2 = i1 + 1
        i3 = i1 + 2
        rh, ph, z = chart[:]
        gc[i1, i1, chart] = 1
        gc[i2, i2, chart] = rh**2
        gc[i3, i3, chart] = 1
        # Orthonormal frame associated with cylindrical coordinates:
        to_orthonormal = self.automorphism_field()
        to_orthonormal[frame, i1, i1, chart] = 1
        to_orthonormal[frame, i2, i2, chart] = 1/rh
        to_orthonormal[frame, i3, i3, chart] = 1
        oframe = frame.new_frame(to_orthonormal, 'e',
                                 indices=(str(rh), str(ph), str(z)),
                                 latex_indices=(latex(rh), latex(ph), latex(z)))
        self._cylindrical_frame = oframe
        g.comp(oframe)

    def _transition_spherical_cartesian(self):
        r"""
        Transitions between spherical and Cartesian coordinates.

        TESTS::

            sage: E = EuclideanSpace(3)
            sage: E._init_coordinates_spherical(r"r th:\theta ph:\phi")
            sage: E.atlas()
            [Chart (E^3, (x, y, z)), Chart (E^3, (r, th, ph))]
            sage: E.coord_changes() # no transition map has been defined yet
            {}
            sage: E._transition_spherical_cartesian()
            sage: len(E.coord_changes())  # 2 transition maps have been created
            2

        Tests of the change-of-frame formulas::

            sage: spher = E.spherical_coordinates()
            sage: spher_f = E.spherical_frame()
            sage: cart_f = E.cartesian_frame()
            sage: E.change_of_frame(cart_f, spher_f)[spher_f,:,spher]
            [cos(ph)*sin(th) cos(ph)*cos(th)        -sin(ph)]
            [sin(ph)*sin(th) cos(th)*sin(ph)         cos(ph)]
            [        cos(th)        -sin(th)               0]
            sage: E.change_of_frame(spher_f, cart_f)[spher_f,:,spher]
            [cos(ph)*sin(th) sin(ph)*sin(th)         cos(th)]
            [cos(ph)*cos(th) cos(th)*sin(ph)        -sin(th)]
            [       -sin(ph)         cos(ph)               0]

        """
        # Transition maps spherical chart <-> Cartesian chart
        chart_cart = self._cartesian_chart
        chart_spher = self._spherical_chart
        x, y, z = chart_cart[:]
        r, th, ph = chart_spher[:]
        spher_to_cart = chart_spher.transition_map(chart_cart,
                             [r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th)])
        spher_to_cart.set_inverse(sqrt(x**2+y**2+z**2),
                                  atan2(sqrt(x**2+y**2),z), atan2(y, x))
        # Automorphisms orthonormal spherical frame <-> Cartesian frame
        oframe = self._spherical_frame
        cframe = chart_cart.frame()
        sframe = chart_spher.frame()
        chg = self._frame_changes
        cframe_to_oframe = chg[(sframe, oframe)] * chg[(cframe, sframe)]
        oframe_to_cframe = chg[(sframe, cframe)] * chg[(oframe, sframe)]
        chg[(cframe, oframe)] = cframe_to_oframe
        chg[(oframe, cframe)] = oframe_to_cframe
        vmodule = self.vector_field_module()
        vmodule._basis_changes[(cframe, oframe)] = cframe_to_oframe
        vmodule._basis_changes[(oframe, cframe)] = oframe_to_cframe

    def _transition_cylindrical_cartesian(self):
        r"""
        Transitions between cylindrical and Cartesian coordinates.

        TESTS::

            sage: E = EuclideanSpace(3)
            sage: E._init_coordinates_cylindrical(r"r ph:\phi z")
            sage: E.atlas()
            [Chart (E^3, (x, y, z)), Chart (E^3, (r, ph, z))]
            sage: E.coord_changes() # no transition map has been defined yet
            {}
            sage: E._transition_cylindrical_cartesian()
            sage: len(E.coord_changes())  # 2 transition maps have been created
            2

        Tests of the change-of-frame formulas::

            sage: cylind = E.cylindrical_coordinates()
            sage: cylind_f = E.cylindrical_frame()
            sage: cart_f= E.cartesian_frame()
            sage: E.change_of_frame(cart_f, cylind_f)[cylind_f,:,cylind]
            [ cos(ph) -sin(ph)        0]
            [ sin(ph)  cos(ph)        0]
            [       0        0        1]
            sage: E.change_of_frame(cylind_f, cart_f)[cylind_f,:,cylind]
            [ cos(ph)  sin(ph)        0]
            [-sin(ph)  cos(ph)        0]
            [       0        0        1]

        """
        # Transition maps cylindrical chart <-> Cartesian chart
        chart_cart = self._cartesian_chart
        chart_cylind = self._cylindrical_chart
        x, y, z = chart_cart[:]
        rh, ph, z = chart_cylind[:]
        cylind_to_cart = chart_cylind.transition_map(chart_cart,
                                                   [rh*cos(ph), rh*sin(ph), z])
        cylind_to_cart.set_inverse(sqrt(x**2+y**2), atan2(y, x), z)
        # Automorphisms orthonormal cylindrical frame <-> Cartesian frame
        oframe = self._cylindrical_frame
        cframe = chart_cart.frame()
        sframe = chart_cylind.frame()
        chg = self._frame_changes
        cframe_to_oframe = chg[(sframe, oframe)] * chg[(cframe, sframe)]
        oframe_to_cframe = chg[(sframe, cframe)] * chg[(oframe, sframe)]
        chg[(cframe, oframe)] = cframe_to_oframe
        chg[(oframe, cframe)] = oframe_to_cframe
        vmodule = self.vector_field_module()
        vmodule._basis_changes[(cframe, oframe)] = cframe_to_oframe
        vmodule._basis_changes[(oframe, cframe)] = oframe_to_cframe

    def _transition_spherical_cylindrical(self):
        r"""
        Transitions between spherical and cylindrical coordinates.

        TESTS::

            sage: E = EuclideanSpace(3, coordinates='cylindrical')
            sage: E._init_coordinates_spherical(r"r th:\theta ph:\phi")
            sage: E.atlas()
            [Chart (E^3, (r, ph, z)), Chart (E^3, (r, th, ph))]
            sage: E.coord_changes()  # no transition map has been defined yet
            {}
            sage: E._transition_spherical_cylindrical()
            sage: len(E.coord_changes())  # 2 transition maps have been created
            2

        Tests of the change-of-frame formulas::

            sage: spher = E.spherical_coordinates()
            sage: spher_frame = E.spherical_frame()
            sage: cylind_frame = E.cylindrical_frame()
            sage: E.change_of_frame(cylind_frame, spher_frame)[:, spher]
            [ sin(th)  cos(th)        0]
            [       0        0        1]
            [ cos(th) -sin(th)        0]
            sage: E.change_of_frame(spher_frame, cylind_frame)[:, spher]
            [ sin(th)        0  cos(th)]
            [ cos(th)        0 -sin(th)]
            [       0        1        0]

        """
        # Transition maps spherical chart <-> cylindrical chart
        cylind = self._cylindrical_chart
        spher = self._spherical_chart
        rh, ph, z = cylind[:]
        r, th, ph = spher[:]
        spher_to_cylind = spher.transition_map(cylind,
                                               [r*sin(th), ph, r*cos(th)])
        spher_to_cylind.set_inverse(sqrt(rh**2 + z**2), atan2(rh,z), ph)
        # Automorphism orthon. cylindrical frame -> orthon. spherical frame
        cf = cylind.frame() # coordinate cylindrical frame
        sf = spher.frame()  # coordinate spherical frame
        ocf = self._cylindrical_frame # orthonormal cylindrical frame
        osf = self._spherical_frame   # orthonormal spherical frame
        chg = self._frame_changes
        oc_to_os = chg[(sf, osf)] * chg[(cf, sf)] * chg[(ocf, cf)]
        # A priori oc_to_os has been computed only in sf frame; its components
        # in osf frame are computed by
        cmp_osf = oc_to_os.comp(osf)
        # The matrix of oc_to_os in ocf frame is equal to that in osf frame:
        cmp_ocf = oc_to_os.add_comp(ocf)
        for i in self.irange():
            for j in self.irange():
                cmp_ocf[[i,j]] = cmp_osf[[i,j]]
        # Automorphism orthon. spherical frame -> orthon. cylindrical frame
        os_to_oc = chg[(cf, ocf)] * chg[(sf, cf)] * chg[(osf, sf)]
        # A priori oc_to_os has been computed only in cf frame; its components
        # in ocf frame are computed by
        cmp_ocf = os_to_oc.comp(ocf)
        # The matrix of os_to_oc in osf frame is equal to that in ocf frame:
        cmp_osf = os_to_oc.add_comp(osf)
        for i in self.irange():
            for j in self.irange():
                cmp_osf[[i,j]] = cmp_ocf[[i,j]]
        # Storage of the results:
        chg[(ocf, osf)] = oc_to_os
        chg[(osf, ocf)] = os_to_oc
        vmodule = self.vector_field_module()
        vmodule._basis_changes[(ocf, osf)] = oc_to_os
        vmodule._basis_changes[(osf, ocf)] = os_to_oc

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
          ``None`` the symbols are generated as `(x,y,z)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of Cartesian coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(3)
            sage: E.cartesian_coordinates()
            Chart (E^3, (x, y, z))

        An example where the Cartesian coordinates have not been previously
        created::

            sage: E = EuclideanSpace(3, coordinates='spherical')
            sage: E.atlas()  # only spherical coordinates exist
            [Chart (E^3, (r, th, ph))]
            sage: E.cartesian_coordinates(symbols='X Y Z')
            Chart (E^3, (X, Y, Z))
            sage: E.atlas()  # the Cartesian chart has been added to the atlas
            [Chart (E^3, (r, th, ph)), Chart (E^3, (X, Y, Z))]

        The coordinate variables are returned by the square bracket operator::

            sage: E.cartesian_coordinates()[1]
            X
            sage: E.cartesian_coordinates()[3]
            Z
            sage: E.cartesian_coordinates()[:]
            (X, Y, Z)

        It is also possible to use the operator ``<,>`` to set symbolic
        variable containing the coordinates::

            sage: E = EuclideanSpace(3, coordinates='spherical')
            sage: cartesian.<u,v,w> = E.cartesian_coordinates()
            sage: cartesian
            Chart (E^3, (u, v, w))
            sage: u, v, w
            (u, v, w)

        The command ``cartesian.<u,v,w> = E.cartesian_coordinates()``
        is actually a shortcut for::

            sage: cartesian = E.cartesian_coordinates(symbols='u v w')
            sage: u, v, w = cartesian[:]

        """
        if self._cartesian_chart is None:
            if symbols is None:
                if names is None:
                    symbols = 'x y z'
                else:
                    symbols = names[0] + ' ' + names[1] + ' ' + names[2]
            self._init_coordinates_cartesian(symbols)
            if self._spherical_chart:
                self._transition_spherical_cartesian()
            if self._cylindrical_chart:
                self._transition_cylindrical_cartesian()
        return self._cartesian_chart

    def spherical_coordinates(self, symbols=None, names=None):
        r"""
        Return the chart of spherical coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as the
          argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the spherical chart has not been already defined; if
          ``None`` the symbols are generated as `(r,\theta,\phi)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of spherical coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(3)
            sage: E.spherical_coordinates()
            Chart (E^3, (r, th, ph))
            sage: latex(E.spherical_coordinates())
            \left(\mathbb{E}^{3},(r, {\theta}, {\phi})\right)

        The relation to Cartesian coordinates is::

            sage: E.coord_change(E.spherical_coordinates(),
            ....:                E.cartesian_coordinates()).display()
            x = r*cos(ph)*sin(th)
            y = r*sin(ph)*sin(th)
            z = r*cos(th)

        The coordinate variables are returned by the square bracket operator::

            sage: E.spherical_coordinates()[1]
            r
            sage: E.spherical_coordinates()[3]
            ph
            sage: E.spherical_coordinates()[:]
            (r, th, ph)

        They can also be obtained via the operator ``<,>``::

            sage: spherical.<r,th,ph> = E.spherical_coordinates()
            sage: spherical
            Chart (E^3, (r, th, ph))
            sage: r, th, ph
            (r, th, ph)

        Actually, ``spherical.<r,th,ph> = E.spherical_coordinates()`` is a
        shortcut for::

            sage: spherical = E.spherical_coordinates()
            sage: r, th, ph = spherical[:]

        The coordinate symbols can be customized::

            sage: E = EuclideanSpace(3)
            sage: E.spherical_coordinates(symbols=r"R T:\Theta F:\Phi")
            Chart (E^3, (R, T, F))
            sage: latex(E.spherical_coordinates())
            \left(\mathbb{E}^{3},(R, {\Theta}, {\Phi})\right)

        Note that if the spherical coordinates have been already initialized,
        the argument ``symbols`` has no effect::

            sage: E.spherical_coordinates(symbols=r"r th:\theta ph:\phi")
            Chart (E^3, (R, T, F))

        """
        if self._spherical_chart is None:
            if symbols is None:
                if names is None:
                    symbols = 'r th:\\theta ph:\\phi'
                else:
                    names = list(names)
                    if names[1] in ['t', 'th', 'theta']:
                        names[1] = names[1] + ':\\theta'
                    if names[2] in ['p', 'ph', 'phi']:
                        names[2] = names[2] + ':\\phi'
                    symbols = names[0] + ' ' + names[1] + ' ' + names[2]
            self._init_coordinates_spherical(symbols)
            if self._cartesian_chart:
                self._transition_spherical_cartesian()
            if self._cylindrical_chart:
                self._transition_spherical_cylindrical()
        return self._spherical_chart

    def spherical_frame(self):
        r"""
        Return the orthonormal vector frame associated with spherical
        coordinates.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`

        EXAMPLES::

            sage: E = EuclideanSpace(3)
            sage: E.spherical_frame()
            Vector frame (E^3, (e_r,e_th,e_ph))
            sage: E.spherical_frame()[1]
            Vector field e_r on the Euclidean space E^3
            sage: E.spherical_frame()[:]
            (Vector field e_r on the Euclidean space E^3,
             Vector field e_th on the Euclidean space E^3,
             Vector field e_ph on the Euclidean space E^3)

        The spherical frame expressed in terms of the Cartesian one::

            sage: for e in E.spherical_frame():
            ....:     e.display(E.cartesian_frame(), E.spherical_coordinates())
            ....:
            e_r = cos(ph)*sin(th) e_x + sin(ph)*sin(th) e_y + cos(th) e_z
            e_th = cos(ph)*cos(th) e_x + cos(th)*sin(ph) e_y - sin(th) e_z
            e_ph = -sin(ph) e_x + cos(ph) e_y

        """
        if self._spherical_frame is None:
            self.spherical_coordinates()  # creates the spherical chart and the
                                          # associated orthonormal frame
        return self._spherical_frame

    def cylindrical_coordinates(self, symbols=None, names=None):
        r"""
        Return the chart of cylindrical coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as the
          argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the cylindrical chart has not been already defined; if
          ``None`` the symbols are generated as `(\rho,\phi,z)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of cylindrical coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(3)
            sage: E.cylindrical_coordinates()
            Chart (E^3, (rh, ph, z))
            sage: latex(E.cylindrical_coordinates())
            \left(\mathbb{E}^{3},({\rho}, {\phi}, z)\right)

        The relation to Cartesian coordinates is::

            sage: E.coord_change(E.cylindrical_coordinates(),
            ....:                E.cartesian_coordinates()).display()
            x = rh*cos(ph)
            y = rh*sin(ph)
            z = z

        The coordinate variables are returned by the square bracket operator::

            sage: E.cylindrical_coordinates()[1]
            rh
            sage: E.cylindrical_coordinates()[3]
            z
            sage: E.cylindrical_coordinates()[:]
            (rh, ph, z)

        They can also be obtained via the operator ``<,>``::

            sage: cylindrical.<rh,ph,z> = E.cylindrical_coordinates()
            sage: cylindrical
            Chart (E^3, (rh, ph, z))
            sage: rh, ph, z
            (rh, ph, z)

        Actually, ``cylindrical.<rh,ph,z> = E.cylindrical_coordinates()`` is a
        shortcut for::

            sage: cylindrical = E.cylindrical_coordinates()
            sage: rh, ph, z = cylindrical[:]

        The coordinate symbols can be customized::

            sage: E = EuclideanSpace(3)
            sage: E.cylindrical_coordinates(symbols=r"R Phi:\Phi Z")
            Chart (E^3, (R, Phi, Z))
            sage: latex(E.cylindrical_coordinates())
            \left(\mathbb{E}^{3},(R, {\Phi}, Z)\right)

        Note that if the cylindrical coordinates have been already initialized,
        the argument ``symbols`` has no effect::

            sage: E.cylindrical_coordinates(symbols=r"rh:\rho ph:\phi z")
            Chart (E^3, (R, Phi, Z))

        """
        if self._cylindrical_chart is None:
            if symbols is None:
                if names is None:
                    symbols = 'rh:\\rho ph:\\phi z'
                else:
                    names = list(names)
                    if names[0] in ['rh', 'rho']:
                        names[0] = names[0] + ':\\rho'
                    if names[1] in ['p', 'ph', 'phi']:
                        names[1] = names[1] + ':\\phi'
                    symbols = names[0] + ' ' + names[1] + ' ' + names[2]
            self._init_coordinates_cylindrical(symbols)
            if self._cartesian_chart:
                self._transition_cylindrical_cartesian()
            if self._spherical_chart:
                self._transition_spherical_cylindrical()
        return self._cylindrical_chart

    def cylindrical_frame(self):
        r"""
        Return the orthonormal vector frame associated with cylindrical
        coordinates.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`

        EXAMPLES::

            sage: E = EuclideanSpace(3)
            sage: E.cylindrical_frame()
            Vector frame (E^3, (e_rh,e_ph,e_z))
            sage: E.cylindrical_frame()[1]
            Vector field e_rh on the Euclidean space E^3
            sage: E.cylindrical_frame()[:]
            (Vector field e_rh on the Euclidean space E^3,
             Vector field e_ph on the Euclidean space E^3,
             Vector field e_z on the Euclidean space E^3)

        The cylindrical frame expressed in terms of the Cartesian one::

            sage: for e in E.cylindrical_frame():
            ....:     e.display(E.cartesian_frame(), E.cylindrical_coordinates())
            ....:
            e_rh = cos(ph) e_x + sin(ph) e_y
            e_ph = -sin(ph) e_x + cos(ph) e_y
            e_z = e_z

        """
        if self._cylindrical_frame is None:
            self.cylindrical_coordinates()  # creates the cylindrical chart and
                                            # the associated orthonormal frame
        return self._cylindrical_frame


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
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`, namely
          ``symbols`` is a string of coordinate fields separated by a blank
          space, where each field contains the coordinate's text symbol and
          possibly the coordinate's LaTeX symbol (when the latter is different
          from the text symbol), both symbols being separated by a colon
          (``:``); if ``None``, the symbols will be automatically generated
          according to the value of ``coordinates``
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

    EXAMPLES:

    Constructing a 2-dimensional Euclidean space::

        sage: E = EuclideanSpace(2); E
        Euclidean plane E^2

    Each call to :func:`EuclideanSpace` creates a different object::

        sage: E1 = EuclideanSpace(2)
        sage: E1 is E
        False
        sage: E1 == E
        False

    The LaTeX symbol of the Euclidean space is by default `\mathbb{E}^n`,
    where `n` is the dimension::

        sage: latex(E)
        \mathbb{E}^{2}

    But both the name and LaTeX names of the Euclidean space can be
    customized::

        sage: F = EuclideanSpace(2, name='F', latex_name=r'\mathcal{F}'); F
        Euclidean plane F
        sage: latex(F)
        \mathcal{F}

    By default, an Euclidean space is created with a single coordinate chart:
    that of Cartesian coordinates::

        sage: E.atlas()
        [Chart (E^2, (x, y))]
        sage: E.cartesian_coordinates()
        Chart (E^2, (x, y))
        sage: E.default_chart() is E.cartesian_coordinates()
        True

    The coordinate variables can be initialized, as the Python variables ``x``
    and ``y``, by::

        sage: x, y = E.cartesian_coordinates()[:]

    However, it is possible to both construct the Euclidean space and
    initialize the coordinate variables in a single stage, thanks to
    SageMath operator ``<,>``::

        sage: E.<x,y> = EuclideanSpace(2)

    The coordinate symbols can be customized::

        sage: E.<X,Y> = EuclideanSpace(2)
        sage: E.cartesian_coordinates()
        Chart (E^2, (X, Y))

    By default, the LaTeX symbols of the coordinates coincide with the text
    ones::

        sage: latex(X+Y)
        X + Y

    However, it is possible to customize them, via the argument ``symbols``,
    which must be a string, usually prefixed by ``r`` (for *raw* string, in
    order to allow for the backslash character of LaTeX expressions). This
    string contains the coordinate fields separated by a blank space; each
    field contains the coordinate's text symbol and possibly the coordinate's
    LaTeX symbol (when the latter is different from the text symbol), both
    symbols being separated by a colon (``:``)::

        sage: E.<xi,ze> = EuclideanSpace(2, symbols=r"xi:\xi ze:\zeta")
        sage: E.cartesian_coordinates()
        Chart (E^2, (xi, ze))
        sage: latex(xi+ze)
        {\xi} + {\zeta}

    Thanks to the argument ``coordinates``, an Euclidean space can be
    constructed with curvilinear coordinates initialized instead of the
    Cartesian ones::

        sage: E.<r,ph> = EuclideanSpace(2, coordinates='polar')
        sage: E.atlas()   # no Cartesian coordinates have been constructed
        [Chart (E^2, (r, ph))]
        sage: polar = E.polar_coordinates(); polar
        Chart (E^2, (r, ph))
        sage: E.default_chart() is polar
        True
        sage: latex(r+ph)
        {\phi} + r

    The Cartesian coordinates, along with the transition maps to and from
    the curvilinear coordinates, can be constructed at any time by::

        sage: cartesian.<x,y> = E.cartesian_coordinates()
        sage: E.atlas()  # both polar and Cartesian coordinates now exist
        [Chart (E^2, (r, ph)), Chart (E^2, (x, y))]

    The transition maps have been initialized by the command
    ``E.cartesian_coordinates()``::

        sage: E.coord_change(polar, cartesian).display()
        x = r*cos(ph)
        y = r*sin(ph)
        sage: E.coord_change(cartesian, polar).display()
        r = sqrt(x^2 + y^2)
        ph = arctan2(y, x)

    The default name of the Euclidean metric tensor is `g`::

        sage: E.metric()
        Riemannian metric g on the Euclidean plane E^2
        sage: latex(_)
        g

    But this can be customized::

        sage: E = EuclideanSpace(2, metric_name='h')
        sage: E.metric()
        Riemannian metric h on the Euclidean plane E^2
        sage: latex(_)
        h
        sage: E = EuclideanSpace(2, metric_latex_name=r'\mathbf{g}')
        sage: E.metric()
        Riemannian metric g on the Euclidean plane E^2
        sage: latex(_)
        \mathbf{g}

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
    if n == 3:
        return Euclidean3dimSpace(name=name, latex_name=latex_name,
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
