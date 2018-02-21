r"""
Euclidean spaces

An *Euclidean space of dimension `n`* is a Riemannian manifold diffeomorphic to
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
    sage: E.default_chart()
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

    sage: r, ph = E.polar_coordinates()

They belong to the polar chart::

    sage: E.polar_chart()
    Chart (E^2, (r, ph))

which is part of the atlas of ``E``::

    sage: E.atlas()
    [Chart (E^2, (x, y)), Chart (E^2, (r, ph))]

The ranges of the coordinates introduced so far are::

    sage: E.polar_chart().coord_range()
    r: (0, +oo); ph: (0, 2*pi)
    sage: E.cartesian_chart().coord_range()
    x: (-oo, +oo); y: (-oo, +oo)

The transition map from polar coordinates to Cartesian one is::

    sage: E.coord_change(E.polar_chart(), E.cartesian_chart()).display()
    x = r*cos(ph)
    y = r*sin(ph)

while the reverse is::

    sage: E.coord_change(E.cartesian_chart(), E.polar_chart()).display()
    r = sqrt(x^2 + y^2)
    ph = arctan2(y, x)

At this stage, ``E`` is endowed with three vector frames::

    sage: E.frames()
    [Coordinate frame (E^2, (d/dx,d/dy)),
     Coordinate frame (E^2, (d/dr,d/dph)),
     Vector frame (E^2, (e_r,e_ph))]

The third one is the standard orthonormal frame associated with polar
coordinates.

The expression of the metric tensor in terms of the polar coordinates is::

    sage: g.display(E.polar_chart().frame(), E.polar_chart())
    g = dr*dr + r^2 dph*dph

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
        self._named_charts = {}
        self._named_frames = {}
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
        self._named_charts['Cartesian'] = chart
        frame = chart.frame()
        self._named_frames['Cartesian'] = frame
        g = self.metric()
        gc = g.add_comp(frame)
        for i in self.irange():
            gc[i, i, chart] = 1
        nabla = g.connection(init_coef=False)  # False to avoid any computation
        nabla.add_coef(frame)  # initialize a zero set of coefficients

    def cartesian_chart(self, symbols=None):
        r"""
        Return the chart of Cartesian coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as
          the argument ``coordinates`` in  :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the Cartesian chart has not been already defined; if
          ``None`` the symbols are generated as `(x_1,\ldots,x_n)`.

        OUTPUT:

        - the chart of Cartesian coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(4)
            sage: X = E.cartesian_chart(); X
            Chart (E^4, (x1, x2, x3, x4))
            sage: X in E.atlas()
            True

        """
        if 'Cartesian' not in self._named_charts:
            if symbols is None:
                symbols = ''
                for i in self.irange():
                    symbols += "x{}".format(i) + r":x_{" + str(i) + r"} "
                symbols = symbols[:-1]
            self._init_coordinates_cartesian(symbols)
        return self._named_charts['Cartesian']

    def cartesian_coordinates(self, symbols=None):
        r"""
        Return Cartesian coordinates, possibly creating them if they do not
        already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as
          the argument ``coordinates`` in  :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the Cartesian chart has not been already defined; if
          ``None`` the symbols are generated as `(x_1,\ldots,x_n)`.

        OUTPUT:

        - tuple of symbolic variables representing the Cartesian coordinates

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.cartesian_coordinates()
            (x, y)

        An example where the Cartesian coordinates have not been previously
        created::

            sage: F.<r,phi> = EuclideanSpace(2, coordinates='polar')
            sage: F.atlas()  # only polar coordinates exist
            [Chart (E^2, (r, phi))]
            sage: F.cartesian_coordinates(symbols='u v')
            (u, v)
            sage: F.atlas()  # the Cartesian chart has been added to the atlas
            [Chart (E^2, (r, phi)), Chart (E^2, (u, v))]

        Note that if the Cartesian coordinates have been previously created,
        the argument ``symbols`` has not effect::

            sage: E.cartesian_coordinates(symbols='u v')
            (x, y)
            sage: E.atlas()
            [Chart (E^2, (x, y))]

        """
        return self.cartesian_chart(symbols=symbols)[:]

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
        self._named_charts['polar'] = chart
        frame = chart.frame()
        self._named_frames['polar_coord'] = frame
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
        self._named_frames['polar_ortho'] = oframe
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
        chart_cart = self._named_charts['Cartesian']
        chart_pol = self._named_charts['polar']
        x, y = chart_cart[:]
        r, ph = chart_pol[:]
        pol_to_cart = chart_pol.transition_map(chart_cart,
                                               [r*cos(ph), r*sin(ph)])
        pol_to_cart.set_inverse(sqrt(x**2+y**2), atan2(y,x))
        # Automorphisms ortho polar frame <-> Cartesian frame
        oframe = self._named_frames['polar_ortho']
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

    def polar_chart(self, symbols=None):
        r"""
        Return the chart of polar coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as
          the argument ``coordinates`` in  :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the polar chart has not been already defined; if
          ``None`` the symbols are generated as `(r,\phi)`.

        OUTPUT:

        - the chart of polar coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.polar_chart()
            Chart (E^2, (r, ph))
            sage: latex(E.polar_chart())
            \left(\mathbb{E}^{2},(r, {\phi})\right)

        Customizing the coordinate symbols::

            sage: E = EuclideanSpace(2)
            sage: E.polar_chart(symbols=r"r th:\theta")
            Chart (E^2, (r, th))
            sage: latex(E.polar_chart())
            \left(\mathbb{E}^{2},(r, {\theta})\right)

        Note that if the polar chart has been already defined, the argument
        ``symbols`` has no effect::

            sage: E.polar_chart(symbols=r"R Th:\Theta")
            Chart (E^2, (r, th))

        """
        if 'polar' not in self._named_charts:
            if symbols is None:
                symbols = 'r ph:\\phi'
            self._init_coordinates_polar(symbols)
            if 'Cartesian' in self._named_charts:
                self._transition_polar_cartesian()
        return self._named_charts['polar']

    def polar_coordinates(self, symbols=None):
        r"""
        Return polar coordinates, possibly creating them if they do not already
        exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as
          the argument ``coordinates`` in  :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the polar chart has not been already defined; if
          ``None`` the symbols are generated as `(r,\phi)`.

        OUTPUT:

        - tuple of symbolic variables representing the Cartesian coordinates

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.polar_coordinates()
            (r, ph)

        The coordinate symbols can be customized::

            sage: E = EuclideanSpace(2)
            sage: E.polar_coordinates(symbols=r"r th:\theta")
            (r, th)
            sage: latex(_)
            \left(r, {\theta}\right)

        Note that if the polar coordinates have been already initialized, the
        argument ``symbols`` has no effect::

            sage: E.polar_coordinates(symbols=r"R Ph:\Phi")
            (r, th)

        """
        return self.polar_chart(symbols=symbols)[:]

    def cartesian_chart(self, symbols=None):
        r"""
        Return the chart of Cartesian coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as
          the argument ``coordinates`` in  :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the Cartesian chart has not been already defined; if
          ``None`` the symbols are generated as `(x,y)`.

        OUTPUT:

        - the chart of Cartesian coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.cartesian_chart()
            Chart (E^2, (x, y))

        An example where the Cartesian coordinates have not been previously
        created::

            sage: E = EuclideanSpace(2, coordinates='polar')
            sage: E.atlas()  # only polar coordinates exist
            [Chart (E^2, (r, ph))]
            sage: E.cartesian_chart(symbols='X Y')
            Chart (E^2, (X, Y))
            sage: E.atlas()  # the Cartesian chart has been added to the atlas
            [Chart (E^2, (r, ph)), Chart (E^2, (X, Y))]

        """
        if 'Cartesian' not in self._named_charts:
            if symbols is None:
                symbols = 'x y'
            self._init_coordinates_cartesian(symbols)
            if 'polar' in self._named_charts:
                self._transition_polar_cartesian()
        return self._named_charts['Cartesian']


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

