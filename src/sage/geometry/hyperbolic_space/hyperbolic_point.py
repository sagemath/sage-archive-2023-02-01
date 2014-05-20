r"""
Hyperbolic Points

This module implements the abstract base class for points in 
hyperbolic space of arbitrary dimension.  It also contains the
implementations for specific models of hyperbolic geometry.

AUTHORS:

- Greg Laun (2013): initial version

EXAMPLES:

We can construct points in the upper half plane model, abbreviated
UHP for convenience::

    sage: UHP.point(2 + I)
    Point in UHP I + 2
    sage: g = UHP.point(3 + I)
    sage: g.dist(UHP.point(I))
    arccosh(11/2)
"""

#***********************************************************************
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

from sage.structure.sage_object import SageObject
from sage.symbolic.pynac import I
from sage.misc.lazy_import import lazy_import
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.latex import latex

lazy_import('sage.rings.all', ['RR', 'CC'])
lazy_import('sage.functions.other', 'real')
lazy_import('sage.modules.free_module_element', 'vector')

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_methods',
            ['HyperbolicAbstractMethods', 'HyperbolicMethodsUHP'])
lazy_import('sage.geometry.hyperbolic_space.model_factory', 'ModelFactory')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_factory',
            ['HyperbolicAbstractFactory', 'HyperbolicFactoryUHP',
             'HyperbolicFactoryPD', 'HyperbolicFactoryKM',
             'HyperbolicFactoryHM'])
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_geodesic', 'HyperbolicGeodesic')


class HyperbolicPoint(SageObject):
    r"""
    Abstract base class for hyperbolic points.  This class should never
    be instantiated.

    INPUT:

    - the coordinates of a hyperbolic point in the appropriate model

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
    """
    HFactory =  HyperbolicAbstractFactory
    HMethods = HyperbolicAbstractMethods

    def __init__(self, coordinates, **graphics_options):
        r"""
        See ``HyperbolicPoint`` for full documentation.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: p = HyperbolicPointUHP(I)
            sage: p
            Point in UHP I
        """
        if self.model().point_in_model(coordinates):
            if isinstance(coordinates, tuple):
                coordinates = vector(coordinates)
            self._coordinates = coordinates
        else:
            raise ValueError(
                "{0} is not a valid".format(coordinates) +
                " point in the {0} model".format(self.model().short_name))
        self._graphics_options = graphics_options

    #####################
    # "Private" Methods #
    #####################

    @lazy_attribute
    def _cached_coordinates(self):
        r"""
        The representation of the current point used for calculations.
        For example, if the current model uses the HyperbolicMethodsUHP
        class, then _cached_coordinates will hold the upper half plane
        representation of self.coordinates().

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: A = HyperbolicPointHM((0, 0, 1))
            sage: A._cached_coordinates
            I
        """
        return self.model().point_to_model(self.coordinates(),
                                           self.HMethods.model_name())
    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: HyperbolicPointUHP(3 + 4*I)
            Point in UHP 4*I + 3

            sage: HyperbolicPointPD(1/2 + I/2)
            Point in PD 1/2*I + 1/2

            sage: HyperbolicPointKM((1/2, 1/2))
            Point in KM (1/2, 1/2)

            sage: HyperbolicPointHM((0,0,1))
            Point in HM (0, 0, 1)
        """
        return "Point in {0} {1}".format(self.model_name(), self.coordinates())

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: p = UHP.point(0)
            sage: latex(p)
            0
            sage: q = HM.point((0,0,1))
            sage: latex(q)
            \left(0,\,0,\,1\right)
        """
        return latex(self.coordinates())

    def __eq__(self, other):
        r"""
        Return ``True`` if ``self`` is equal to ``other``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: p1 = HyperbolicPointUHP(1 + I)
            sage: p2 = HyperbolicPointUHP(2 + I)
            sage: p1 == p2
            False
            sage: p1 == p1
            True

        ::

            sage: p1 = HyperbolicPointPD(0)
            sage: p2 = HyperbolicPointPD(1/2 + 2*I/3)
            sage: p1 == p2
            False
            sage: p1 == p1
            True

        ::

            sage: p1 = HyperbolicPointKM((0,0))
            sage: p2 = HyperbolicPointKM((0, 1/2))
            sage: p1 == p2
            False

        ::

            sage: p1 = HyperbolicPointHM((0,0,1))
            sage: p2 = HyperbolicPointHM((0,0,1/1))
            sage: p1 == p2
            True
        """
        return (self.model_name() == other.model_name()
                and bool(self.coordinates() == other.coordinates()))

    def __rmul__(self, other):
        r"""
        Implement the action of matrices on points of hyperbolic space.

        EXAMPLES::

            sage: A = matrix(2, [0, 1, 1, 0])
            sage: A*UHP.point(2 + I)
            Point in UHP 1/5*I + 2/5
            sage: B = diagonal_matrix([-1, -1, 1])
            sage: B*HM.point((0, 1, sqrt(2)))
            Point in HM (0, -1, sqrt(2))
        """
        from sage.matrix.matrix import is_Matrix
        if is_Matrix(other):
            A = self.HFactory.get_isometry(other)
            return A*self
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

    @classmethod
    def model(cls):
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
        return cls.HFactory.get_model()

    @classmethod
    def model_name(cls):
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
        return cls.model().short_name

    def to_model(self, model_name):
        r"""
        Convert the current object to image in another model.

        INPUT:

        - ``model_name`` -- a string representing the image model

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: HyperbolicPointUHP(I).to_model('PD')
            Point in PD 0
        """
        factory = ModelFactory.find_factory(model_name)
        coordinates = self.model().point_to_model(self.coordinates(), model_name)
        return factory.get_point(coordinates)

    def update_graphics(self, update=False, **options):
        r"""
        Update the graphics options of a HyperbolicPoint.  If ``update`` is
        ``True``, update rather than overwrite.

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

    ###################################
    # Methods implemented in HMethods #
    ###################################

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
            return self.HMethods.point_dist(self._cached_coordinates, tmp_other._cached_coordinates)
        elif isinstance(other, HyperbolicGeodesic):
            return self.HMethods.geod_dist_from_point(
                *(other._cached_endpoints + self._cached_coordinates
              ))

    def symmetry_in (self):
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

            sage: A*A == UHP.isometry(identity_matrix(2))
            True
        """
        A = self.HMethods.symmetry_in(self._cached_coordinates)
        A = self.HMethods.model().isometry_to_model(A, self.model_name())
        return self.HFactory.get_isometry(A)

    @classmethod
    def random_element(cls, **kwargs):
        r"""
        Return a random point in the upper half
        plane.  The points are uniformly distributed over the rectangle
        `[-10, 10] \times [-10 i, 10 i]`.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_point import *
            sage: p = HyperbolicPointUHP.random_element()
            sage: bool((p.coordinates().imag()) > 0)
            True

            sage: p = HyperbolicPointPD.random_element()
            sage: PD.point_in_model(p.coordinates())
            True

            sage: p = HyperbolicPointKM.random_element()
            sage: KM.point_in_model(p.coordinates())
            True

            sage: p = HyperbolicPointHM.random_element().coordinates()
            sage: bool((p[0]**2 + p[1]**2 - p[2]**2 - 1) < 10**-8)
            True
        """
        p = cls.HMethods.random_point(**kwargs)
        p = cls.HMethods.model().point_to_model(p, cls.model().short_name)
        return cls.HFactory.get_point(p)

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
        """
        opts = dict([('axes', False),('aspect_ratio',1)])
        opts.update(self.graphics_options())
        opts.update(options)
        from sage.misc.functional import numerical_approx
        p = self.coordinates()
        if p in RR:
            p = CC(p)
        elif hasattr(p, 'iteritems') or hasattr(p, '__iter__'):
            p = map(numerical_approx, p)
        else:
            p = numerical_approx(p)
        from sage.plot.point import point
        pic = point(p, **opts)
        if boundary:
            bd_pic = self.HFactory.get_background_graphic()
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
    """
    HFactory = HyperbolicFactoryUHP
    HMethods = HyperbolicMethodsUHP

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
            bd_pic = self.HFactory.get_background_graphic(bd_min = cent - 1,
                                                          bd_max = cent + 1)
            pic = bd_pic + pic
        return pic

class HyperbolicPointPD (HyperbolicPoint):
    r"""
    Create a point in the PD model.

    INPUT:

    - the coordinates of a point in the unit disk in the complex plane `\CC`

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPointPD
        sage: HyperbolicPointPD(0)
        Point in PD 0
    """
    HFactory = HyperbolicFactoryPD
    HMethods = HyperbolicMethodsUHP

class HyperbolicPointKM(HyperbolicPoint):
    r"""
    Create a point in the KM model.

    INPUT:

    - the coordinates of a point in the unit disk in the real plane `\RR^2`

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPointKM
        sage: HyperbolicPointKM((0,0))
        Point in KM (0, 0)
    """
    HFactory = HyperbolicFactoryKM
    HMethods = HyperbolicMethodsUHP

class HyperbolicPointHM (HyperbolicPoint):
    r"""
    Create a point in the HM model.

    INPUT:

    - the coordinates of a point in the hyperboloid given
      by `x^2 + y^2 - z^2 = -1`

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPointHM
        sage: HyperbolicPointHM((0,0,1))
        Point in HM (0, 0, 1)
    """
    HFactory = HyperbolicFactoryHM
    HMethods = HyperbolicMethodsUHP

