r"""
Hyperbolic Points

This module implements points in hyperbolic space of arbitrary dimension.
It also contains the implementations for specific models of
hyperbolic geometry.

This module also implements ideal points in hyperbolic space of arbitrary
dimension. It also contains the implementations for specific models
of hyperbolic geometry.

Note that not all models of hyperbolic space are bounded, meaning that
the ideal boundary is not the topological boundary of the set underlying
tho model.  For example, the unit disk model is bounded with boundary
given by the unit sphere.  The hyperboloid model is not bounded.

AUTHORS:

- Greg Laun (2013): initial version

EXAMPLES:

We can construct points in the upper half plane model, abbreviated
UHP for convenience::

    sage: HyperbolicPlane.UHP.point(2 + I)
    Point in UHP I + 2
    sage: g = HyperbolicPlane.UHP.point(3 + I)
    sage: g.dist(HyperbolicPlane.UHP.point(I))
    arccosh(11/2)

We can also construct boundary points in the upper half plane model::

    sage: HyperbolicPlane.UHP.point(3)
    Boundary point in UHP 3

Points on the boundary are infinitely far from interior points::

    sage: HyperbolicPlane.UHP.point(3).dist(HyperbolicPlane.UHP.point(I))
    +Infinity
"""

#***********************************************************************
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

from sage.structure.element import Element
from sage.symbolic.pynac import I
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.latex import latex
from sage.geometry.hyperbolic_space.hyperbolic_isometry import HyperbolicIsometry
from sage.rings.infinity import infinity
from sage.rings.all import RR, CC

from sage.misc.lazy_import import lazy_import
lazy_import('sage.functions.other', 'real')
lazy_import('sage.modules.free_module_element', 'vector')

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_methods',
            ['HyperbolicAbstractMethods', 'HyperbolicMethodsUHP'])
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_geodesic', 'HyperbolicGeodesic')


class HyperbolicPoint(Element):
    r"""
    Abstract base class for hyperbolic points.  This class should never
    be instantiated.

    INPUT:

    - ``model`` -- the model of the hyperbolic space
    - ``coordinates`` -- the coordinates of a hyperbolic point in the
      appropriate model
    - ``is_boundary`` -- whether the point is a boundary point
    - ``check`` -- (default: ``True``) if ``True``, then check to make sure
      the coordinates give a valid point the the model

    EXAMPLES:

    Note that the coordinate representation does not differentiate the
    different models::

        sage: p = HyperbolicPlane().UHP().get_point(.2 + .3*I); p
        Point in UHP 0.200000000000000 + 0.300000000000000*I

        sage: q = HyperbolicPlane().PD().get_point(0.2 + 0.3*I); q
        Point in PD 0.200000000000000 + 0.300000000000000*I

        sage: p == q
        False

        sage: bool(p.coordinates() == q.coordinates())
        True

    Similarly for boundary points::

        sage: p = HyperbolicPlane().UHP().get_point(1); p
        Boundary point in UHP 1

        sage: q = HyperbolicPlane().PD().get_point(1); q
        Boundary point in PD 1

        sage: p == q
        False

        sage: bool(p.coordinates() == q.coordinates())
        True

    It is an error to specify a point that does not lie in the
    appropriate model::

        sage: HyperbolicPlane().UHP().get_point(0.2 - 0.3*I)
        Traceback (most recent call last):
        ...
        ValueError: 0.200000000000000 - 0.300000000000000*I is not a valid point in the UHP model

        sage: HyperbolicPlane().PD().get_point(1.2)
        Traceback (most recent call last):
        ...
        ValueError: 1.20000000000000 is not a valid point in the PD model

        sage: HyperbolicPlane().KM().get_point((1,1))
        Traceback (most recent call last):
        ...
        ValueError: (1, 1) is not a valid point in the KM model

        sage: HyperbolicPlane().HM().get_point((1, 1, 1))
        Traceback (most recent call last):
        ...
        ValueError: (1, 1, 1) is not a valid point in the HM model

    It is an error to specify an interior point of hyperbolic space as a
    boundary point::

        sage: HyperbolicPlane().UHP().get_point(0.2 + 0.3*I, boundary=True)
        Traceback (most recent call last):
        ...
        ValueError: 0.200000000000000 + 0.300000000000000*I is not a valid boundary point in the UHP model
    """
    def __init__(self, model, coordinates, is_boundary, check=True, **graphics_options):
        r"""
        See ``HyperbolicPoint`` for full documentation.

        EXAMPLES::

            sage: p = HyperbolicPlane().UHP().get_point(I)
            sage: TestSuite(p).run()
        """
        if is_boundary:
            if not model.is_bounded():
                raise NotImplementedError("boundary points are not implemented in the {0} model".format(model.short_name()))
            if check and not model.boundary_point_in_model(coordinates):
                raise ValueError(
                    "{0} is not a valid".format(coordinates) +
                    " boundary point in the {0} model".format(model.short_name()))
        elif check and not model.point_in_model(coordinates):
            raise ValueError(
                "{0} is not a valid".format(coordinates) +
                " point in the {0} model".format(model.short_name()))

        if isinstance(coordinates, tuple):
            coordinates = vector(coordinates)
        self._coordinates = coordinates
        self._bdry = is_boundary
        self._graphics_options = graphics_options

        Element.__init__(self, model)

    #####################
    # "Private" Methods #
    #####################

    @lazy_attribute
    def _cached_coordinates(self):
        r"""
        The representation of the current point used for calculations.

        EXAMPLES::

            sage: A = HyperbolicPlane().HM().get_point((0, 0, 1))
            sage: A._cached_coordinates
            I
        """
        R = self.parent().realization_of().a_realization()
        return R(self).coordinates()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_point(3 + 4*I)
            Point in UHP 4*I + 3

            sage: HyperbolicPlane().PD().get_point(1/2 + I/2)
            Point in PD 1/2*I + 1/2

            sage: HyperbolicPlane().KM().get_point((1/2, 1/2))
            Point in KM (1/2, 1/2)

            sage: HyperbolicPlane().HM().get_point((0,0,1))
            Point in HM (0, 0, 1)

            sage: HyperbolicPlane().UHP().get_point(infinity)
            Boundary point in UHP +Infinity

            sage: HyperbolicPlane().PD().get_point(-1)
            Boundary point in PD -1

            sage: HyperbolicPlane().KM().get_point((0, -1))
            Boundary point in KM (0, -1)
        """
        if self._bdry:
            base = "Boundary point"
        else:
            base = "Point"
        return base + " in {0} {1}".format(self.parent().short_name(), self._coordinates)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: p = HyperbolicPlane().UHP().get_point(0)
            sage: latex(p)
            0
            sage: q = HyperbolicPlane().HM().get_point((0,0,1))
            sage: latex(q)
            \left(0,\,0,\,1\right)
        """
        return latex(self._coordinates)

    def __eq__(self, other):
        r"""
        Return ``True`` if ``self`` is equal to ``other``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: p1 = HyperbolicPlane().UHP().get_point(1 + I)
            sage: p2 = HyperbolicPlane().UHP().get_point(2 + I)
            sage: p1 == p2
            False
            sage: p1 == p1
            True

        ::

            sage: p1 = HyperbolicPlane().PD().get_point(0)
            sage: p2 = HyperbolicPlane().PD().get_point(1/2 + 2*I/3)
            sage: p1 == p2
            False
            sage: p1 == p1
            True

        ::

            sage: p1 = HyperbolicPlane().KM().get_point((0,0))
            sage: p2 = HyperbolicPlane().KM().get_point((0, 1/2))
            sage: p1 == p2
            False

        ::

            sage: p1 = HyperbolicPlane().HM().get_point((0,0,1))
            sage: p2 = HyperbolicPlane().HM().get_point((0,0,1/1))
            sage: p1 == p2
            True
        """
        return (self.parent() is other.parent()
                and bool(self._coordinates == other._coordinates))

    # TODO: Add a test that this works with isometries
    def __rmul__(self, other):
        r"""
        Implement the action of matrices on points of hyperbolic space.

        EXAMPLES::

            sage: A = matrix(2, [0, 1, 1, 0])
            sage: A * HyperbolicPlane.UHP.point(2 + I)
            Point in UHP 1/5*I + 2/5
            sage: B = diagonal_matrix([-1, -1, 1])
            sage: B * HyperbolicPlane.HM.point((0, 1, sqrt(2)))
            Point in HM (0, -1, sqrt(2))
        """
        from sage.matrix.matrix import is_Matrix
        if is_Matrix(other):
            A = self.parent().get_isometry(other)
            return A * self
        elif isinstance(other, HyperbolicIsometry):
            return other * self
        else:
            raise TypeError("unsupported operand type(s) for *:"
                            "{0} and {1}".format(self, other))

    #######################
    # Setters and Getters #
    #######################

    def coordinates(self):
        r"""
        Return the coordinates of the point.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_point(2 + I).coordinates()
            I + 2

            sage: HyperbolicPlane().PD().get_point(1/2 + 1/2*I).coordinates()
            1/2*I + 1/2

            sage: HyperbolicPlane().KM().get_point((1/3, 1/4)).coordinates()
            (1/3, 1/4)

            sage: HyperbolicPlane().HM().get_point((0,0,1)).coordinates()
            (0, 0, 1)
        """
        return self._coordinates

    def model(self):
        r"""
        Return the model to which the :class:`HyperbolicPoint` belongs.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: HyperbolicPlane().UHP().get_point(I).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelUHP'>

            sage: HyperbolicPlane().PD().get_point(0).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelPD'>

            sage: HyperbolicPlane().KM().get_point((0,0)).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelKM'>

            sage: HyperbolicPlane().HM().get_point((0,0,1)).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelHM'>
        """
        return self.parent()

    def is_boundary(self):
        """
        Return ``True`` if ``self`` is a boundary point.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: p = PD.get_point(0.5+.2*I)
            sage: p.is_boundary()
            False
            sage: p = PD.get_point(I)
            sage: p.is_boundary()
            True
        """
        return self._bdry

    def update_graphics(self, update=False, **options):
        r"""
        Update the graphics options of a :class:`HyperbolicPoint`.
        If ``update`` is ``True``, update rather than overwrite.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: p = HyperbolicPointUHP(I); p.graphics_options()
            {}

            sage: p.update_graphics(color = "red"); p.graphics_options()
            {'color': 'red'}

            sage: p.update_graphics(color = "blue"); p.graphics_options()
            {'color': 'blue'}

            sage: p.update_graphics(True, size = 20); p.graphics_options()
            {'color': 'blue', 'size': 20}
        """
        if not update:
            self._graphics_options = {}
        self._graphics_options.update(**options)

    def graphics_options(self):
        r"""
        Return the graphics options of the current point.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: p = HyperbolicPointUHP(2 + I, color="red")
            sage: p.graphics_options()
            {'color': 'red'}
        """
        return self._graphics_options

    ####################################
    # Methods implemented in _HMethods #
    ####################################

    def symmetry_in(self):
        r"""
        Return the involutary isometry fixing the given point.

        EXAMPLES::

            sage: z = HyperbolicPlane().UHP().get_point(3 + 2*I)
            sage: z.symmetry_in()
            Isometry in UHP
            [  3/2 -13/2]
            [  1/2  -3/2]

            sage: HyperbolicPlane().UHP().get_point(I).symmetry_in()
            Isometry in UHP
            [ 0 -1]
            [ 1  0]

            sage: HyperbolicPlane().PD().get_point(0).symmetry_in()
            Isometry in PD
            [-I  0]
            [ 0  I]

            sage: HyperbolicPlane().KM().get_point((0, 0)).symmetry_in()
            Isometry in KM
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]

            sage: HyperbolicPlane().HM().get_point((0,0,1)).symmetry_in()
            Isometry in HM
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]

            sage: p = HyperbolicPlane().UHP().random_element()
            sage: A = p.symmetry_in()
            sage: A*p == p
            True

            sage: A.orientation_preserving()
            True

            sage: A*A == HyperbolicPlane().UHP().isometry(identity_matrix(2))
            True
        """
        R = self.parent().realization_of().a_realization()
        A = self._HMethods.symmetry_in(self._cached_coordinates)
        A = self._HMethods.model().isometry_to_model(A, self.model_name())
        return self.parent().get_isometry(A)

    ###########
    # Display #
    ###########

    def show(self, boundary=True, **options):
        r"""
        EXAMPLES::

            sage: HyperbolicPlane().PD().get_point(0).show()
            sage: HyperbolicPlane().KM().get_point((0,0)).show()
            sage: HyperbolicPlane().HM().get_point((0,0,1)).show()
        """
        p = self.coordinates()
        if p == infinity:
            raise NotImplementedError("can't draw the point infinity")

        opts = dict([('axes', False), ('aspect_ratio',1)])
        opts.update(self.graphics_options())
        opts.update(options)

        if self._bdry: # It is a boundary point
            p = numerical_approx(p)
            pic = point((p, 0), **opts)
            if boundary:
                bd_pic = self._model.get_background_graphic(bd_min = p - 1,
                                                            bd_max = p + 1)
                pic = bd_pic + pic
            return pic

        # It is an interior point
        from sage.misc.functional import numerical_approx
        if p in RR:
            p = CC(p)
        elif hasattr(p, 'iteritems') or hasattr(p, '__iter__'):
            p = map(numerical_approx, p)
        else:
            p = numerical_approx(p)
        from sage.plot.point import point
        pic = point(p, **opts)
        if boundary:
            bd_pic = self.parent().get_background_graphic()
            pic = bd_pic + pic
        return pic

class HyperbolicPointUHP(HyperbolicPoint):
    r"""
    A point in the UHP model.

    INPUT:

    - the coordinates of a point in the unit disk in the complex plane `\CC`

    EXAMPLES::

        sage: HyperbolicPlane().UHP().get_point(2*I)
        Point in UHP 2*I

        sage: HyperbolicPlane().UHP().get_point(1)
        Boundary point in UHP 1

    """
    def symmetry_in(self):
        r"""
        Return the involutary isometry fixing the given point.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.symmetry_in(3 + 2*I)
            [  3/2 -13/2]
            [  1/2  -3/2]
        """
        p = self._coordinates
        x, y = real(p), imag(p)
        if y > 0:
            return matrix(2,[x/y,-(x**2/y) - y,1/y,-(x/y)])

    def show(self, boundary=True, **options):
        r"""
        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_point(I).show()

            sage: HyperbolicPlane().UHP().get_point(0).show()
            sage: HyperbolicPlane().UHP().get_point(infinity).show()
            Traceback (most recent call last):
            ...
            NotImplementedError: can't draw the point infinity
        """
        # FIXME: Something didn't get put into the UHP point's show() properly
        opts = dict([('axes', False),('aspect_ratio',1)])
        opts.update(self.graphics_options())
        opts.update(options)
        from sage.misc.functional import numerical_approx
        p = self.coordinates() + 0*I
        p = numerical_approx(p)
        from sage.plot.point import point
        pic = point(p, **opts)
        if boundary:
            cent = real(p)
            bd_pic = self.parent().get_background_graphic(bd_min = cent - 1,
                                                          bd_max = cent + 1)
            pic = bd_pic + pic
        return pic

class HyperbolicPointPD(HyperbolicPoint):
    r"""
    Create a point in the PD model.

    INPUT:

    - the coordinates of a point in the unit disk in the complex plane `\CC`

    EXAMPLES::

        sage: HyperbolicPlane().PD().get_point(0)
        Point in PD 0

        sage: HyperbolicPlane().PD().get_point(1)
        Boundary point in PD 1
    """

class HyperbolicPointKM(HyperbolicPoint):
    r"""
    Create a point in the KM model.

    INPUT:

    - the coordinates of a point in the unit disk in the real plane `\RR^2`

    EXAMPLES::

        sage: HyperbolicPlane().KM().get_point((0,0))
        Point in KM (0, 0)

        sage: HyperbolicPlane().KM().get_point((1,0))
        Boundary point in KM (1, 0)
    """

class HyperbolicPointHM(HyperbolicPoint):
    r"""
    Create a point in the HM model.

    INPUT:

    - the coordinates of a point in the hyperboloid given
      by `x^2 + y^2 - z^2 = -1`

    EXAMPLES::

        sage: HyperbolicPlane().HM().get_point((0,0,1))
        Point in HM (0, 0, 1)

        sage: HyperbolicPlane().HM().get_point((1,0,0), is_boundary=True)
        Traceback (most recent call last):
        ...
        NotImplementedError: boundary points are not implemented in the HM model
    """

