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
from sage.misc.lazy_import import lazy_import
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.latex import latex
from sage.geometry.hyperbolic_space.hyperbolic_isometry import HyperbolicIsometry

lazy_import('sage.rings.all', ['RR', 'CC'])
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

        sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
        sage: p = HyperbolicPointUHP(.2 + .3*I); p
        Point in UHP 0.200000000000000 + 0.300000000000000*I

        sage: q = HyperbolicPointPD(0.2 + 0.3*I); q
        Point in PD 0.200000000000000 + 0.300000000000000*I

        sage: p == q
        False

        sage: bool(p.coordinates() == q.coordinates())
        True

    Similarly for boundary points::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import *
        sage: p = HyperbolicBdryPointUHP(1); p
        Boundary point in UHP 1

        sage: q = HyperbolicBdryPointPD(1); q
        Boundary point in PD 1

        sage: p == q
        False

        sage: bool(p.coordinates() == q.coordinates())
        True

    It is an error to specify a point that does not lie in the
    appropriate model::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
        sage: HyperbolicPointUHP(0.2 - 0.3*I)
        Traceback (most recent call last):
        ...
        ValueError: 0.200000000000000 - 0.300000000000000*I is not a valid point in the UHP model

        sage: HyperbolicPointPD(1.2)
        Traceback (most recent call last):
        ...
        ValueError: 1.20000000000000 is not a valid point in the PD model

        sage: HyperbolicPointKM((1,1))
        Traceback (most recent call last):
        ...
        ValueError: (1, 1) is not a valid point in the KM model

        sage: HyperbolicPointHM((1, 1, 1))
        Traceback (most recent call last):
        ...
        ValueError: (1, 1, 1) is not a valid point in the HM model

    It is an error to specify an interior point of hyperbolic space as a
    boundary point::

        sage: HyperbolicBdryPointUHP(0.2 + 0.3*I)
        Traceback (most recent call last):
        ...
        ValueError: 0.200000000000000 + 0.300000000000000*I is not a valid boundary point in the UHP model
    """
    def __init__(self, model, coordinates, is_boundary, check=True, **graphics_options):
        r"""
        See ``HyperbolicPoint`` for full documentation.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: p = HyperbolicPointUHP(I)
            sage: p
            Point in UHP I
        """
        if is_boundary:
            if not model.is_bounded():
                raise NotImplementedError("boundary points are not implemented in the {0} model".format(model.short_name()))
            if check and not model.bdry_point_in_model(coordinates):
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
        For example, if the current model uses the HyperbolicMethodsUHP
        class, then ``_cached_coordinates`` will hold the upper half plane
        representation of ``self.coordinates()``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: A = HyperbolicPointHM((0, 0, 1))
            sage: A._cached_coordinates
            I
        """
        return self.parent().point_to_model(self.coordinates(),
                                            self._HMethods.model_name())
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

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: HyperbolicPointUHP(2 + I).coordinates()
            I + 2

            sage: HyperbolicPointPD(1/2 + 1/2*I).coordinates()
            1/2*I + 1/2

            sage: HyperbolicPointKM((1/3, 1/4)).coordinates()
            (1/3, 1/4)

            sage: HyperbolicPointHM((0,0,1)).coordinates()
            (0, 0, 1)
        """
        return self._coordinates

    def model(self):
        r"""
        Return the model to which the :class:`HyperbolicPoint` belongs.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: HyperbolicPointUHP(I).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelUHP'>

            sage: HyperbolicPointPD(0).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelPD'>

            sage: HyperbolicPointKM((0,0)).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelKM'>

            sage: HyperbolicPointHM((0,0,1)).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelHM'>
        """
        return self.parent()

    def model_name(self):
        r"""
        Return the short name of the hyperbolic model.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: HyperbolicPointUHP(I).model_name()
            'UHP'

            sage: HyperbolicPointPD(0).model_name()
            'PD'

            sage: HyperbolicPointKM((0,0)).model_name()
            'KM'

            sage: HyperbolicPointHM((0,0,1)).model_name()
            'HM'
        """
        return self.parent().short_name

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

    def dist(self, other):
        r"""
        Calculate the hyperbolic distance between two points in the same
        model.

        INPUT:

        - ``other`` -- a hyperbolic point in the same model as ``self``

        OUTPUT:

        - the hyperbolic distance

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: p1 = HyperbolicPointUHP(5 + 7*I)
            sage: p2 = HyperbolicPointUHP(1.0 + I)
            sage: p1.dist(p2)
            2.23230104635820

            sage: p1 = HyperbolicPointPD(0)
            sage: p2 = HyperbolicPointPD(I/2)
            sage: p1.dist(p2)
            arccosh(5/3)

            sage: p1.to_model('UHP').dist(p2.to_model('UHP'))
            arccosh(5/3)

            sage: p1 = HyperbolicPointKM((0, 0))
            sage: p2 = HyperbolicPointKM((1/2, 1/2))
            sage: numerical_approx(p1.dist(p2))
            0.881373587019543

            sage: p1 = HyperbolicPointHM((0,0,1))
            sage: p2 = HyperbolicPointHM((1,0,sqrt(2)))
            sage: numerical_approx(p1.dist(p2))
            0.881373587019543

        Distance between a point and itself is 0::

            sage: p = HyperbolicPointUHP(47 + I)
            sage: p.dist(p)
            0
        """
        tmp_other = other.to_model(self.model_name())
        if isinstance(other, HyperbolicPoint):
            return self._HMethods.point_dist(self._cached_coordinates, tmp_other._cached_coordinates)
        elif isinstance(other, HyperbolicGeodesic):
            return self._HMethods.geod_dist_from_point(
                *(other._cached_endpoints + self._cached_coordinates
              ))

    def symmetry_in(self):
        r"""
        Return the involutary isometry fixing the given point.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: z = HyperbolicPointUHP(3 + 2*I)
            sage: z.symmetry_in()
            Isometry in UHP
            [  3/2 -13/2]
            [  1/2  -3/2]

            sage: HyperbolicPointUHP(I).symmetry_in()
            Isometry in UHP
            [ 0 -1]
            [ 1  0]

            sage: HyperbolicPointPD(0).symmetry_in()
            Isometry in PD
            [-I  0]
            [ 0  I]

            sage: HyperbolicPointKM((0, 0)).symmetry_in()
            Isometry in KM
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]

            sage: HyperbolicPointHM((0,0,1)).symmetry_in()
            Isometry in HM
            [-1  0  0]
            [ 0 -1  0]
            [ 0  0  1]

            sage: p = HyperbolicPointUHP.random_element()
            sage: A = p.symmetry_in()
            sage: A*p == p
            True

            sage: A.orientation_preserving()
            True

            sage: A*A == HyperbolicPlane.UHP.isometry(identity_matrix(2))
            True
        """
        A = self._HMethods.symmetry_in(self._cached_coordinates)
        A = self._HMethods.model().isometry_to_model(A, self.model_name())
        return self.parent().get_isometry(A)

    ###########
    # Display #
    ###########

    def show(self, boundary=True, **options):
        r"""
        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: HyperbolicPointUHP(I).show()
            sage: HyperbolicPointPD(0).show()
            sage: HyperbolicPointKM((0,0)).show()
            sage: HyperbolicPointHM((0,0,1)).show()

            sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import *
            sage: HyperbolicBdryPointUHP(0).show()
            sage: HyperbolicBdryPointUHP(infinity).show()
            Traceback (most recent call last):
            ...
            NotImplementedError: can't draw the point infinity
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

        sage: from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPointUHP
        sage: HyperbolicPointUHP(2*I)
        Point in UHP 2*I

        sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import HyperbolicBdryPointUHP
        sage: q = HyperbolicBdryPointUHP(1); q
        Boundary point in UHP 1

    """
    _HMethods = HyperbolicMethodsUHP

    def show(self, boundary=True, **options):
        r"""
        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: HyperbolicPointUHP(I).show()
            sage: HyperbolicPointPD(0).show()
            sage: HyperbolicPointKM((0,0)).show()
            sage: HyperbolicPointHM((0,0,1)).show()
        """
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

        sage: from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPointPD
        sage: HyperbolicPointPD(0)
        Point in PD 0

        sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import HyperbolicBdryPointPD
        sage: q = HyperbolicBdryPointPD(1); q
        Boundary point in PD 1
    """
    _HMethods = HyperbolicMethodsUHP

class HyperbolicPointKM(HyperbolicPoint):
    r"""
    Create a point in the KM model.

    INPUT:

    - the coordinates of a point in the unit disk in the real plane `\RR^2`

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPointKM
        sage: HyperbolicPointKM((0,0))
        Point in KM (0, 0)

        sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import HyperbolicBdryPointKM
        sage: q = HyperbolicBdryPointKM((1,0)); q
        Boundary point in KM (1, 0)
    """
    _HMethods = HyperbolicMethodsUHP

class HyperbolicPointHM(HyperbolicPoint):
    r"""
    Create a point in the HM model.

    INPUT:

    - the coordinates of a point in the hyperboloid given
      by `x^2 + y^2 - z^2 = -1`

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPointHM
        sage: HyperbolicPointHM((0,0,1))
        Point in HM (0, 0, 1)

        sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import HyperbolicBdryPointHM
        sage: q = HyperbolicBdryPointHM((1,0,0)); q
        Traceback (most recent call last):
        ...
        NotImplementedError: boundary points are not implemented in the HM model
    """
    _HMethods = HyperbolicMethodsUHP

