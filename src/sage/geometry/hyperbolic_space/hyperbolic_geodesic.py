r"""
Hyperbolic Geodesics

This module implements the abstract base class for geodesics in 
hyperbolic space of arbitrary dimension.  It also contains the
implementations for specific models of hyperbolic geometry.

AUTHORS:

- Greg Laun (2013): initial version

EXAMPLES:

We can construct geodesics in the upper half plane model, abbreviated
UHP for convenience::

    sage: HyperbolicPlane().UHP().get_geodesic(2, 3)
    Geodesic in UHP from 2 to 3
    sage: g = HyperbolicPlane().UHP().get_geodesic(I, 3 + I)
    sage: g.length()
    arccosh(11/2)

Geodesics are oriented, which means that two geodesics with the same
graph will only be equal if their starting and ending points are
the same::

    sage: HyperbolicPlane().UHP().get_geodesic(1,2) == HyperbolicPlane().UHP().get_geodesic(2,1)
    False

.. TODO::

    Implement a parent for all geodesics of the hyperbolic plane?
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
from sage.geometry.hyperbolic_space.hyperbolic_constants import EPSILON
from sage.rings.infinity import infinity
from sage.rings.all import CC, RR
from sage.symbolic.constants import pi

lazy_import('sage.functions.other', 'real')
lazy_import('sage.functions.other', 'imag')

lazy_import('sage.functions.trig', 'cos')
lazy_import('sage.functions.trig', 'sin')
lazy_import('sage.plot.line', 'line')
lazy_import('sage.modules.free_module_element', 'vector')
lazy_import('sage.functions.other','sqrt')
lazy_import('sage.functions.hyperbolic', 'cosh')
lazy_import('sage.functions.hyperbolic', 'sinh')
lazy_import('sage.functions.hyperbolic', 'arcsinh')

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_methods',
            ['HyperbolicAbstractMethods', 'HyperbolicMethodsUHP'])

class HyperbolicGeodesic(SageObject):
    r"""
    Abstract base class for oriented geodesics that are not necessarily
    complete.

    INPUT:

    - ``start`` -- a HyperbolicPoint or coordinates of a point in
      hyperbolic space representing the start of the geodesic

    - ``end`` -- a HyperbolicPoint or coordinates of a point in
      hyperbolic space representing the end of the geodesic

    EXAMPLES::

        sage: UHP = HyperbolicPlane().UHP()
        sage: g = UHP.get_geodesic(UHP.get_point(I), UHP.get_point(2 + I))
        sage: g = UHP.get_geodesic(I, 2 + I)
    """

    #####################
    # "Private" Methods #
    #####################

    def __init__(self, model, start, end, **graphics_options):
        r"""
        See :class:`HyperbolicGeodesic` for full documentation.

        EXAMPLES ::

            sage: HyperbolicPlane().UHP().get_geodesic(I, 2 + I)
            Geodesic in UHP from I to I + 2
        """
        self._model = model
        self._start = start
        self._end = end
        self._graphics_options = graphics_options

    @lazy_attribute
    def _cached_start(self):
        r"""
        The representation of the start point used for calculations.

        EXAMPLES::

            sage: HyperbolicPlane().PD().get_geodesic(0, 1/2)._cached_start
            I
        """
        M = self._model.realization_of().a_realization()
        return M(self.start()).coordinates()

    @lazy_attribute
    def _cached_end(self):
        r"""
        The representation of the end point used for calculations.

        EXAMPLES::

            sage: HyperbolicPlane().PD().get_geodesic(0, 1/2)._cached_end
            3/5*I + 4/5
        """
        M = self._model.realization_of().a_realization()
        return M(self.end()).coordinates()

    @lazy_attribute
    def _cached_endpoints(self):
        r"""
        The representation of the endpoints used for calculations.

        EXAMPLES::

            sage: A = HyperbolicPlane().PD().get_geodesic(0, 1/2)
            sage: A._cached_endpoints
            [I, 3/5*I + 4/5]
        """
        return [self._cached_start, self._cached_end]

    @lazy_attribute
    def _complete(self):
        r"""
        Return whether the geodesic is complete.  This is used for
        geodesics in non-bounded models.  For thse models,
        ``self.complete()`` simply sets ``_complete`` to ``True``.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_geodesic(1, -12)._complete
            True
            sage: HyperbolicPlane().UHP().get_geodesic(I, 2 + I)._complete
            False
            sage: g = HyperbolicPlane().HM().get_geodesic((0,0,1), (0,1, sqrt(2)))
            sage: g._complete
            False
            sage: g.complete()._complete
            True
        """ 
        if self._model.is_bounded():
            return (self._start.is_boundary() and self._end.is_boundary())
        return False #All non-bounded geodesics start life incomplete.

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_geodesic(3 + 4*I, I)
            Geodesic in UHP from 4*I + 3 to I

            sage: HyperbolicPlane().PD().get_geodesic(1/2 + I/2, 0)
            Geodesic in PD from 1/2*I + 1/2 to 0

            sage: HyperbolicPlane().KM().get_geodesic((1/2, 1/2), (0, 0))
            Geodesic in KM from (1/2, 1/2) to (0, 0)

            sage: HyperbolicPlane().HM().get_geodesic((0,0,1), (0, 1, sqrt(Integer(2))))
            Geodesic in HM from (0, 0, 1) to (0, 1, sqrt(2))
        """
        return "Geodesic in {0} from {1} to {2}".format(self._model.short_name(),
                self._start.coordinates(), self._end.coordinates())

    def __eq__(self, other):
        r"""
        Return ``True`` if ``self`` is equal to other as an oriented geodesic.

        EXAMPLES::

            sage: g1 = HyperbolicPlane().UHP().get_geodesic(I, 2*I)
            sage: g2 = HyperbolicPlane().UHP().get_geodesic(2*I,I)
            sage: g1 == g2
            False
            sage: g1 == g1
            True
        """
        return (self._model is other._model
                and self._start == other._start
                and self._end == other._end)

    #######################
    # Setters and Getters #
    #######################

    def start(self):
        r"""
        Return the starting point of the geodesic.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(I, 3*I)
            sage: g.start()
            Point in UHP I
        """
        return self._start

    def end(self):
        r"""
        Return the starting point of the geodesic.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(I, 3*I)
            sage: g.end()
            Point in UHP 3*I
        """
        return self._end

    def endpoints(self):
        r"""
        Return a list containing the start and endpoints.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(I, 3*I)
            sage: g.endpoints()
            [Point in UHP I, Point in UHP 3*I]
        """
        return [self._start, self._end]

    def model(self):
        r"""
        Return the model to which the :class:`HyperbolicGeodesic` belongs.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_geodesic(I, 2*I).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelUHP'>

            sage: HyperbolicPlane().PD().get_geodesic(0, I/2).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelPD'>

            sage: HyperbolicPlane().KM().get_geodesic((0, 0), (0, 1/2)).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelKM'>

            sage: HyperbolicPlane().HM().get_geodesic((0, 0, 1), (0, 1, sqrt(2))).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelHM'>
        """
        return self._model

    def model_name(self):
        r"""
        Return the short name of the hyperbolic model.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_geodesic(I, 2*I).model_name()
            'UHP'

            sage: HyperbolicPlane().PD().get_geodesic(0, I/2).model_name()
            'PD'

            sage: HyperbolicPlane().KM().get_geodesic((0, 0), (0, 1/2)).model_name()
            'KM'

            sage: HyperbolicPlane().HM().get_geodesic((0, 0, 1), (0, 1, sqrt(2))).model_name()
            'HM'
        """
        return self.model().short_name()

    def to_model(self, model):
        r"""
        Convert the current object to image in another model.

        INPUT:

        - ``model`` -- the image model

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_geodesic(I, 2*I).to_model('PD')
            Geodesic in PD from 0 to 1/3*I
        """
        if not model.is_bounded() and self.is_complete():
            g = self.uncomplete()
            return g.to_model(model).complete()
        start = model(self._start)
        end = model(self._end)
        g = model.get_geodesic(start, end)
        if not self._model.is_bounded() and model.is_bounded() and self.is_complete():
            # Converting from a non-bounded model to a bounded model
            return g.complete()
        return g

    def graphics_options(self):
        r"""
        Return the graphics options of ``self``.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(I, 2*I, color="red")
            sage: g.graphics_options()
            {'color': 'red'}
        """
        return self._graphics_options

    def update_graphics(self, update=False, **options):
        r"""
        Update the graphics options of ``self``.

        INPUT:

        - ``update`` -- if ``True``, the original option are updated
          rather than overwritten

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(I, 2*I)
            sage: g.graphics_options()
            {}

            sage: g.update_graphics(color = "red"); g.graphics_options()
            {'color': 'red'}

            sage: g.update_graphics(color = "blue"); g.graphics_options()
            {'color': 'blue'}

            sage: g.update_graphics(True, size = 20); g.graphics_options()
            {'color': 'blue', 'size': 20}
        """
        if not update:
            self._graphics_options = {}
        self._graphics_options.update(**options)

    ###################
    # Boolean Methods #
    ###################

    def is_complete(self):
        r"""
        Return ``True`` if ``self`` is a complete geodesic (that is, both
        endpoints are on the ideal boundary) and ``False`` otherwise.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_geodesic(I, 2*I).is_complete()
            False

            sage: HyperbolicPlane().UHP().get_geodesic(0, I).is_complete()
            False

            sage: HyperbolicPlane().UHP().get_geodesic(0, infinity).is_complete()
            True
        """
        return self._complete

    def is_asymptotically_parallel(self, other):
        r"""
        Return ``True`` if ``self`` and ``other`` are asymptotically
        parallel and ``False`` otherwise.

        INPUT:

        - ``other`` -- a hyperbolic geodesic

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: h = HyperbolicPlane().UHP().get_geodesic(-2,4)
            sage: g.is_asymptotically_parallel(h)
            True

        Ultraparallel geodesics are not asymptotically parallel::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: h = HyperbolicPlane().UHP().get_geodesic(-1,4)
            sage: g.is_asymptotically_parallel(h)
            False

        No hyperbolic geodesic is asymptotically parallel to itself::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: g.is_asymptotically_parallel(g)
            False
        """
        p1, p2 = self.complete().endpoints()
        q1, q2 = other.complete().endpoints()
        return ((self != other) and ((p1 in [q1, q2]) or (p2 in [q1,q2]))
                and self.model() is other.model())

    def is_ultra_parallel(self,other):
        r"""
        Return ``True`` if ``self`` and ``other`` are ultra parallel
        and ``False`` otherwise.

        INPUT:

        - ``other`` -- a hyperbolic geodesic

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_geodesic import *
            sage: g = HyperbolicPlane().UHP().get_geodesic(0,1)
            sage: h = HyperbolicPlane().UHP().get_geodesic(-3,3)
            sage: g.is_ultra_parallel(h)
            True

        ::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: h = HyperbolicPlane().UHP().get_geodesic(2,6)
            sage: g.is_ultra_parallel(h)
            False

        ::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: g.is_ultra_parallel(g)
            False
        """
        [R_self, R_other] = [k.reflection_in() for k in [self,other]]
        return (R_self*R_other).classification() == 'hyperbolic'

    def is_parallel(self, other):
        r"""
        Return ``True`` if the two given hyperbolic geodesics are either
        ultra parallel or asymptotically parallel and``False`` otherwise.

        INPUT:

        - ``other`` -- a hyperbolic geodesic in any model

        OUTPUT:

        ``True`` if the given geodesics are either ultra parallel or
        asymptotically parallel, ``False`` if not.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: h = HyperbolicPlane().UHP().get_geodesic(5,12)
            sage: g.is_parallel(h)
            True

        ::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: h = HyperbolicPlane().UHP().get_geodesic(-2,4)
            sage: g.is_parallel(h)
            True

        No hyperbolic geodesic is either ultra parallel or
        asymptotically parallel to itself::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: g.is_parallel(g)
            False
        """
        [R_self, R_other] = [k.reflection_in() for k in [self,other]]
        return (R_self*R_other).classification() in ['parabolic', 'hyperbolic']

    ###################################
    # Methods implemented in _HMethods #
    ###################################

    def ideal_endpoints(self):
        r"""
        Return the ideal endpoints in bounded models.  Raise a
        ``NotImplementedError`` in models that are not bounded.

        EXAMPLES::

            sage: H = HyperbolicPlane()
            sage: H.UHP().get_geodesic(1 + I, 1 + 3*I).ideal_endpoints()
            [Boundary point in UHP 1, Boundary point in UHP +Infinity]

            sage: H.PD().get_geodesic(0, I/2).ideal_endpoints()
            [Boundary point in PD -I, Boundary point in PD I]

            sage: H.KM().get_geodesic((0,0), (0, 1/2)).ideal_endpoints()
            [Boundary point in KM (0, -1), Boundary point in KM (0, 1)]

            sage: H.HM().get_geodesic((0,0,1), (1, 0, sqrt(2))).ideal_endpoints()
            Traceback (most recent call last):
            ...
            NotImplementedError: boundary points are not implemented in the HM model
        """
        if not self.model().is_bounded():
            raise NotImplementedError("boundary points are not implemented in the "
                                          + "{0} model".format(self.model_name()))
        if self.is_complete():
            return self.endpoints()
        ends = self._HMethods.boundary_points(*self._cached_endpoints)
        ends = [self._HMethods.model().point_to_model(k, self.model_name()) for
                k in ends]
        return [self._model.get_bdry_point(k) for k in ends]

    def complete(self):
        r"""
        Return the ideal endpoints in bounded models.  Raise a
        ``NotImplementedError`` in models that are not bounded.

        EXAMPLES::

            sage: H = HyperbolicPlane()
            sage: H.UHP().get_geodesic(1 + I, 1 + 3*I).complete()
            Geodesic in UHP from 1 to +Infinity

            sage: H.PD().get_geodesic(0, I/2).complete()
            Geodesic in PD from -I to I

            sage: H.KM().get_geodesic((0,0), (0, 1/2)).complete()
            Geodesic in KM from (0, -1) to (0, 1)

            sage: H.HM().get_geodesic((0,0,1), (1, 0, sqrt(2))).complete()
            Geodesic in HM from (0, 0, 1) to (1, 0, sqrt(2))

            sage: H.HM().get_geodesic((0,0,1), (1, 0, sqrt(2))).complete().is_complete()
            True
        """
        if self._model.is_bounded():
            return self._model.get_geodesic(*self.ideal_endpoints())

        from copy import copy
        g = copy(self)
        g._complete = True
        return g

    def uncomplete(self):
        r"""
        Return a geodesic segment whose completion is the same as that
        of ``self``.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(I, 2 + 3*I)
            sage: g.uncomplete()
            Geodesic in UHP from I to 3*I + 2

            sage: g.uncomplete().complete() == g.complete()
            True

            sage: h = HyperbolicPlane().UHP().get_geodesic(2, 3)
            sage: h.uncomplete().complete()
            Geodesic in UHP from 2 to 3
        """
        if not self.is_complete():
            return self
        ends = self._HMethods.uncomplete(*self._cached_endpoints)
        ends = [self._HMethods.model().point_to_model(k, self.model_name())
                for k in ends]
        return self._model.get_geodesic(*ends)

    def reflection_in(self):
        r"""
        Return the involution fixing ``self``.

        EXAMPLES::

            sage: H = HyperbolicPlane()
            sage: H.UHP().get_geodesic(2,4).reflection_in()
            Isometry in UHP
            [ 3 -8]
            [ 1 -3]

            sage: H.PD().get_geodesic0, I).reflection_in()
            Isometry in PD
            [ 0 -1]
            [ 1  0]

            sage: H.HM().get_geodesic((0,0), (0,1)).reflection_in()
            Isometry in KM
            [-1  0  0]
            [ 0  1  0]
            [ 0  0  1]

            sage: A = H.HM().get_geodesic((0,0,1), (1,0, n(sqrt(2)))).reflection_in()
            sage: B = diagonal_matrix([1, -1, 1])
            sage: bool((B - A.matrix()).norm() < 10**-9)
            True
        """
        A = self._HMethods.reflection_in(*self._cached_endpoints)
        A = self._HMethods.model().isometry_to_model(A, self.model_name())
        return self._model.get_isometry(A)

    def common_perpendicular(self, other, **graphics_options):
        r"""
        Return the unique hyperbolic geodesic perpendicular to two given
        geodesics, if such a geodesic exists.  If none exists, raise a
        ``ValueError``.

        INPUT:

        - ``other`` -- a hyperbolic geodesic in the same model as ``self``

        OUTPUT:

        - a hyperbolic geodesic

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2,3)
            sage: h = HyperbolicPlane().UHP().get_geodesic(4,5)
            sage: g.common_perpendicular(h)
            Geodesic in UHP from 1/2*sqrt(3) + 7/2 to -1/2*sqrt(3) + 7/2

        It is an error to ask for the common perpendicular of two
        intersecting geodesics::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2,4)
            sage: h = HyperbolicPlane().UHP().get_geodesic(3, infinity)
            sage: g.common_perpendicular(h)
            Traceback (most recent call last):
            ...
            ValueError: geodesics intersect, no common perpendicular exists
        """
        if not self.is_parallel(other):
            raise ValueError('geodesics intersect, no common perpendicular exists')
        perp_ends = self._HMethods.common_perpendicular(
            *(self._cached_endpoints + other._cached_endpoints))
        M = self._HMethods.model()
        perp_ends = [M.point_to_model(k, self.model_name())
                     for k in perp_ends]
        return self._model.get_geodesic(*perp_ends, **graphics_options)

    def intersection(self, other, **graphics_options):
        r"""
        Return the point of intersection of two geodesics (if such a
        point exists).

        The option ``as_complete`` determines whether we test for the
        completed geodesics to intersect, or just the segments.

        INPUT:

        - ``other`` -- a hyperbolic geodesic in the same model as ``self``

        OUTPUT:

        - a hyperbolic point

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g = UHP.get_geodesic(3,5)
            sage: h = UHP.get_geodesic(4,7)
            sage: g.intersection(h)
            Point in UHP 2/3*sqrt(-2) + 13/3

        If the given geodesics do not intersect, raise an error::

            sage: g = UHP.get_geodesic(4,5)
            sage: h = UHP.get_geodesic(5,7)
            sage: g.intersection(h)
            Traceback (most recent call last):
            ...
            ValueError: geodesics don't intersect

        If the given geodesics are identical, return that geodesic::

            sage: g = UHP.get_geodesic(4+I,18*I)
            sage: h = UHP.get_geodesic4+I,18*I)
            sage: g.intersection(h)
            Geodesic in UHP from I + 4 to 18*I
        """
        if self == other:
            return self
        elif self.is_parallel(other):
            raise ValueError("geodesics don't intersect")
        inters = self._HMethods.intersection(*(self._cached_endpoints +
                                              other._cached_endpoints))
        if len(inters) == 2:
            return self
        elif len(inters) == 1:
            inters =  self._HMethods.model().point_to_model(inters[0],
                                                           self.model_name())
            return self._model.get_point(inters, **graphics_options)
        else:
            raise ValueError("can't calculate the intersection of"
                             "{1} and {2}".format(self, other))


    def perpendicular_bisector(self, **graphics_options):
        r"""
        Return the perpendicular bisector of ``self`` if ``self`` has
        finite length.  Here distance is hyperbolic distance.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().random_geodesic()
            sage: h = g.perpendicular_bisector()
            sage: bool(h.intersection(g).coordinates() - g.midpoint().coordinates() < 10**-9)
            True

        Complete geodesics cannot be bisected::

            sage: g = HyperbolicPlane().UHP().get_geodesic(0, 1)
            sage: g.perpendicular_bisector()
            Traceback (most recent call last):
            ...
            ValueError: perpendicular bisector is not defined for complete geodesics
        """
        UHP = self._model.realization_of().UHP
        P = self.to_model(UHP).perpendicular_bisector(**graphics_options)
        return P.to_model(self._model)

    def midpoint(self, **graphics_options):
        r"""
        Return the (hyperbolic) midpoint of a hyperbolic line segment.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().random_geodesic()
            sage: m = g.midpoint()
            sage: end1, end2 = g.endpoints()
            sage: bool(abs(m.dist(end1) - m.dist(end2)) < 10**-9)
            True

        Complete geodesics have no midpoint::

            sage: HyperbolicPlane().UHP().get_geodesic(0,2).midpoint()
            Traceback (most recent call last):
            ...
            ValueError: midpoint is not defined for complete geodesics
        """
        UHP = self._model.realization_of().UHP
        P = self.to_model(UHP).midpoint(**graphics_options)
        return self._model(P)

    def dist(self, other):
        r"""
        Return the hyperbolic distance from a given hyperbolic geodesic
        to another geodesic or point.

        INPUT:

        - ``other`` -- a hyperbolic geodesic or hyperbolic point in
          the same model

        OUTPUT:

        - the hyperbolic distance

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2, 4.0)
            sage: h = HyperbolicPlane().UHP().get_geodesic(5, 7.0)
            sage: bool(abs(g.dist(h).n() - 1.92484730023841) < 10**-9)
            True

        If the second object is a geodesic ultraparallel to the first,
        or if it is a point on the boundary that is not one of the
        first object's endpoints, then return +infinity::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2, 2+I)
            sage: p = HyperbolicPlane().UHP().get_point(5)
            sage: g.dist(p)
            +Infinity
        """
        return self._model.dist(self, other)

    def angle(self, other):
        r"""
        Return the angle  between any two given geodesics if they
        intersect.

        - ``other`` -- a hyperbolic geodesic in the same model as ``self``

        OUTPUT:

        - the angle in radians between the two given geodesics

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2, 4)
            sage: h = HyperbolicPlane().UHP().get_geodesic(3, 3+I)
            sage: g.angle(h)
            1/2*pi

        It is an error to ask for the angle of two geodesics that do not
        intersect::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2, 4)
            sage: h = HyperbolicPlane().UHP().get_geodesic(5, 7)
            sage: g.angle(h)
            Traceback (most recent call last):
            ...
            ValueError: geodesics do not intersect

        If the geodesics are identical, return angle 0::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2, 4)
            sage: g.angle(g)
            0
        """
        if self.is_parallel(other):
            raise ValueError("geodesics do not intersect")
        if not (self.is_complete() and other.is_complete()):
            try:
                # Make sure the segments intersect.
                self.intersection(other)
            except ValueError:
                print("Warning: Geodesic segments do not intersect. "
                      "The angle between them is not defined.\n"
                      "Returning the angle between their completions.")
                return self.complete().angle(other.complete())
        return self._HMethods.angle(*(self._cached_endpoints +
                                     other._cached_endpoints))

    def length(self):
        r"""
        Return the Hyperbolic length of the hyperbolic line segment.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2 + I, 3 + I/2)
            sage: g.length()
            arccosh(9/4)
        """
        return self._model._point_dist(self._start.coordinates(), self._end.coordinates())

#####################################################################
## UHP geodesics

class HyperbolicGeodesicUHP(HyperbolicGeodesic):
    r"""
    Create a geodesic in the upper half plane model.

    INPUT:

    - ``start`` -- a :class:`HyperbolicPoint` or coordinates of a point
      in hyperbolic space representing the start of the geodesic

    - ``end`` -- a :class:`HyperbolicPoint` or coordinates of a point
      in hyperbolic space representing the end of the geodesic

    EXAMPLES::

        sage: UHP = HyperbolicPlane().UHP()
        sage: g = UHP.get_geodesic(UHP.point(I), UHP.point(2 + I))
        sage: g = UHP.get_geodesic(I, 2 + I)
    """
    def reflection_in(self, start, end):
        r"""
        Return the matrix of the involution fixing the geodesic through
        ``start`` and ``end``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.reflection_in(0, 1)
            [ 1  0]
            [ 2 -1]
            sage: HyperbolicMethodsUHP.reflection_in(I, 2*I)
            [ 1  0]
            [ 0 -1]
        """
        x, y = [real(k) for k in self.boundary_points(start, end)]
        if x == infinity:
            M = matrix(2, [[1,-2*y],[0,-1]])
        elif y == infinity:
            M = matrix(2, [[1,-2*x],[0,-1]])
        else:
            M = matrix(2, [[(x+y)/(y-x),- 2*x*y/(y-x)], [2/(y-x), -(x+y)/(y-x)]])
        return self._model.get_isometry(M)

    def show(self, boundary=True, **options):
        r"""
        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_geodesic(0, 1).show()
        """
        opts = dict([('axes', False), ('aspect_ratio',1)])
        opts.update(self.graphics_options())
        opts.update(options)
        end_1, end_2 = [CC(k.coordinates()) for k in self.endpoints()]
        bd_1, bd_2 = [CC(k.coordinates()) for k in self.ideal_endpoints()]
        if (abs(real(end_1) - real(end_2)) < EPSILON) \
                or CC(infinity) in [end_1, end_2]: #on same vertical line
            # If one of the endpoints is infinity, we replace it with a
            # large finite  point
            if end_1 == CC(infinity):
                end_1 = (real(end_2) ,(imag(end_2) + 10))
                end_2 = (real(end_2), imag(end_2))
            elif end_2 == CC(infinity):
                end_2 = (real(end_1), (imag(end_1) + 10))
                end_1 = (real(end_1),imag(end_1))
            pic = line((end_1,end_2), **opts)
            if boundary:
                cent = min(bd_1,bd_2)
                bd_dict = {'bd_min': cent - 3, 'bd_max': cent + 3}
                bd_pic = self._model.get_background_graphic(**bd_dict)
                pic = bd_pic + pic
                return pic
        else:
            center = (bd_1 + bd_2)/2 # Circle center
            radius = abs(bd_1 - bd_2)/2
            theta1 = CC(end_1 - center).arg()
            theta2 = CC(end_2 - center).arg()
            if abs(theta1 - theta2) < EPSILON:
                theta2 += pi
            [theta1, theta2] = sorted([theta1,theta2])
            from sage.calculus.var import var
            from sage.plot.plot import parametric_plot
            x = var('x')
            pic= parametric_plot((radius*cos(x) + real(center),radius*sin(x) +
                                  imag(center)), (x, theta1, theta2), **opts)
            if boundary:
                # We want to draw a segment of the real line.  The
                # computations below compute the projection of the
                # geodesic to the real line, and then draw a little
                # to the left and right of the projection.
                shadow_1, shadow_2 = [real(k) for k in [end_1,end_2]]
                midpoint = (shadow_1 + shadow_2)/2
                length = abs(shadow_1 - shadow_2)
                bd_dict = {'bd_min': midpoint - length, 'bd_max': midpoint +
                           length}
                bd_pic = self._model.get_background_graphic(**bd_dict)
                pic = bd_pic + pic
            return pic

    def perpendicular_bisector(self, **graphics_options): #UHP
        r"""
        Return the perpendicular bisector of the hyperbolic geodesic with
        endpoints ``start`` and ``end`` if that geodesic has finite length.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: a, b = [HyperbolicMethodsUHP.random_point() for k in range(2)]
            sage: g = [a, b]
            sage: h = HyperbolicMethodsUHP.perpendicular_bisector(*g)
            sage: bool(HyperbolicMethodsUHP.intersection(*(g + h))[0] - HyperbolicMethodsUHP.midpoint(*g) < 10**-9)
            True

        Infinite geodesics cannot be bisected::

            sage: HyperbolicMethodsUHP.perpendicular_bisector(0, 1)
            Traceback (most recent call last):
            ...
            ValueError: the length must be finite
        """
        if self.length() == infinity:
            raise ValueError("the length must be finite")
        d = self._dist_points(self._start.coordinates(), self._end.coordinates()) / 2
        end_1, end_2 = self.boundary_points(self._start, self._end)
        S = self._to_std_geod(end_1, self._start, end_2)
        T1 = matrix(2,[exp(d/2),0,0,exp(-d/2)])*
        T2 = matrix(2,[cos(pi/4),-sin(pi/4),sin(pi/4),cos(pi/4)])
        H = S.inverse() * (T1 * T2) * S
        L = [self._model.get_point(mobius_transform(H, k)) for k in [end_1, end_2]]
        return self._model.get_geodesic(L[0], L[1], **graphics_options)

    def midpoint(self): #UHP
        r"""
        Return the (hyperbolic) midpoint of ``self`` if .

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: a, b = [HyperbolicMethodsUHP.random_point() for k in range(2)]
            sage: g = [a, b]
            sage: m = HyperbolicMethodsUHP.midpoint(*g)
            sage: d1 =HyperbolicMethodsUHP.point_dist(m, g[0])
            sage: d2 = HyperbolicMethodsUHP.point_dist(m, g[1])
            sage: bool(abs(d1 - d2) < 10**-9)
            True

        Infinite geodesics have no midpoint::

            sage: HyperbolicMethodsUHP.midpoint(0, 2)
            Traceback (most recent call last):
            ...
            ValueError: the length must be finite
        """
        if self.length() == infinity:
            raise ValueError("the length must be finite")

        start = self._start.coordinates()
        end = self._end.coordinates()
        d = self._dist_points(start, end) / 2
        end_1, end_2 = self.boundary_points(self._start, self._end)
        S = self._to_std_geod(end_1, self._start, end_2)
        T = matrix(2, [[exp(half_dist), 0], [0, 1]])
        M = S.inverse() * T * S
        if ((real(start - end) < EPSILON) or
                (abs(real(start - end)) < EPSILON and
                 imag(start - end) < EPSILON)):
            end_p = start
        else:
            end_p = end
        return mobius_transform(M, end_p)

    def angle(self, other): #UHP
        r"""
        Return the angle between any two given completed geodesics if
        they intersect.

        INPUT:

        -``other`` -- a hyperbolic geodesic in the UHP model

        OUTPUT:

        - the angle in radians between the two given geodesics

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: numerical_approx(HyperbolicMethodsUHP.angle(2, 4, 3, 3 + I))
            1.57079632679490

        If the geodesics are identical, return angle 0::

            sage: HyperbolicMethodsUHP.angle(2, 4, 2, 4)
            0
        """
        start2,end2 = other.endpoints()
        (p_1,p_2) = sorted(self._model.boundary_points(start_1, end_1))
        (q_1,q_2) = sorted(self._model.boundary_points(start_2, end_2))
        # if the geodesics are equal, the angle between them is 0
        if (abs(p_1 - q_1) < EPSILON \
                and abs(p_2 - q_2) < EPSILON):
            return 0
        elif p_2 != infinity: # geodesic not a straight line
            # So we send it to the geodesic with endpoints [0, oo]
            T = self._crossratio_matrix(p_1, (p_1+p_2)/2, p_2)
        else:
            # geodesic is a straight line, so we send it to the geodesic
            # with endpoints [0,oo]
            T = self._crossratio_matrix(p_1, p_1 +1, p_2)
        # b_1 and b_2 are the endpoints of the image of other
        b_1, b_2 = [mobius_transform(T, k) for k in [q_1, q_2]]
        # If other is now a straight line...
        if (b_1 == infinity or b_2 == infinity):
            # then since they intersect, they are equal
            return 0
        else:
            return real(arccos((b_1+b_2)/abs(b_2-b_1)))


class HyperbolicGeodesicPD(HyperbolicGeodesic):
    r"""
    Create a geodesic in the Poincare disk model.

    INPUT:

    - ``start`` -- a :class:`HyperbolicPoint` or coordinates of a
      point in hyperbolic space representing the start of the geodesic

    - ``end`` -- a :class:`HyperbolicPoint` or coordinates of a point
      in hyperbolic space representing the end of the geodesic

    EXAMPLES::

        sage: PD = HyperbolicPlane().PD()
        sage: g = PD.get_geodesic(PD.get_point(I), PD.get_point(I/2))
        sage: g = PD.get_geodesic(I, I/2)
    """
    _HMethods = HyperbolicMethodsUHP

    def show(self, boundary=True, **options):
        r"""
        EXAMPLES::

            sage: HyperbolicPlane().PD().get_geodesic(0, 1).show()
        """
        opts = dict([('axes', False), ('aspect_ratio',1)])
        opts.update(self.graphics_options())
        opts.update(options)
        end_1, end_2 = [CC(k.coordinates()) for k in self.endpoints()]
        bd_1, bd_2 = [CC(k.coordinates()) for k in self.ideal_endpoints()]
        # Check to see if it's a line
        if bool (real(bd_1)*imag(bd_2) - real(bd_2)*imag(bd_1))**2 < EPSILON:
            pic = line([(real(bd_1),imag(bd_1)),(real(bd_2),imag(bd_2))],
                       **opts)
        else:
            # If we are here, we know it's not a line
            # So we compute the center and radius of the circle
            center = (1/(real(bd_1)*imag(bd_2)-real(bd_2)*imag(bd_1))*
                ((imag(bd_2)-imag(bd_1)) + (real(bd_1)-real(bd_2))*I))
            radius = RR(abs(bd_1 - center)) # abs is Euclidean distance
            # Now we calculate the angles for the parametric plot
            theta1 = CC(end_1- center).arg()
            theta2 = CC(end_2 - center).arg()
            if theta2 < theta1:
                theta1, theta2 = theta2, theta1
            from sage.calculus.var import var
            from sage.plot.plot import parametric_plot
            x = var('x')
            mid = (theta1 + theta2)/2.0
            if (radius*cos(mid) + real(center))**2 + \
               (radius*sin(mid) + imag(center))**2 > 1.0:
                # Swap theta1 and theta2
                tmp = theta1 + 2*pi
                theta1 = theta2
                theta2 = tmp
                pic = parametric_plot((radius*cos(x) + real(center),
                                       radius*sin(x) + imag(center)),
                                      (x, theta1, theta2), **opts)

            else:
                pic = parametric_plot((radius*cos(x) + real(center),
                                   radius*sin(x) + imag(center)),
                                  (x, theta1, theta2), **opts)
        if boundary:
            bd_pic = self._model.get_background_graphic()
            pic = bd_pic + pic
        return pic


class HyperbolicGeodesicKM(HyperbolicGeodesic):
    r"""
    Create a geodesic in the Klein disk model.

    INPUT:

    - ``start`` -- a :class:`HyperbolicPoint` or coordinates of a
      point in hyperbolic space representing the start of the geodesic

    - ``end`` -- a :class:`HyperbolicPoint` or coordinates of a point
      in hyperbolic space representing the end of the geodesic

    EXAMPLES::

        sage: KM = HyperbolicPlane().KM()
        sage: g = KM.get_geodesic(KM.point((0,1)), KM.point((0,1/2)))
        sage: g = KM.get_geodesic((0,1), (0,1/2))
    """
    _HMethods = HyperbolicMethodsUHP

    def show(self, boundary=True, **options):
        r"""
        EXAMPLES::

            sage: HyperbolicPlane().KM().get_geodesic((0,0), (1,0)).show()
        """
        from sage.plot.line import line
        opts = dict ([('axes', False), ('aspect_ratio', 1)])
        opts.update(self.graphics_options())
        end_1,end_2 = [k.coordinates() for k in self.endpoints()]
        pic = line([end_1,end_2], **opts)
        if boundary:
            bd_pic = self._model.get_background_graphic()
            pic = bd_pic + pic
        return pic


class HyperbolicGeodesicHM(HyperbolicGeodesic):
    r"""
    Create a geodesic in the hyperboloid model.

    INPUT:

    - ``start`` -- a :class:`HyperbolicPoint` or coordinates of a
      point in hyperbolic space representing the start of the geodesic

    - ``end`` -- a :class:`HyperbolicPoint` or coordinates of a point
      in hyperbolic space representing the end of the geodesic

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_geodesic import *
        sage: HM = HyperbolicPlane().HM()
        sage: g = HM.get_geodesic(HM.point((0, 0, 1)), HM.point((0,1,sqrt(2))))
        sage: g = HM.get_geodesic((0, 0, 1), (0, 1, sqrt(2)))
    """
    _HMethods = HyperbolicMethodsUHP

    def show(self, show_hyperboloid=True, **graphics_options):
        r"""
        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_geodesic import *
            sage: g = HyperbolicPlane().HM().random_geodesic()
            sage: g.show()
        """
        from sage.calculus.var import var
        (x,y,z) = var('x,y,z')
        opts = self.graphics_options()
        opts.update(graphics_options)
        v1,u2 = [vector(k.coordinates()) for k in self.endpoints()]
        # Lorentzian Gram Shmidt.  The original vectors will be
        # u1, u2 and the orthogonal ones will be v1, v2.  Except
        # v1 = u1, and I don't want to declare another variable,
        # hence the odd naming convention above.
        # We need the Lorentz dot product of v1 and u2.
        v1_ldot_u2 = u2[0]*v1[0] + u2[1]*v1[1] - u2[2]*v1[2]
        v2 = u2 + v1_ldot_u2*v1
        v2_norm = sqrt(v2[0]**2 + v2[1]**2 - v2[2]**2)
        v2 = v2/v2_norm
        v2_ldot_u2 =  u2[0]*v2[0] + u2[1]*v2[1] - u2[2]*v2[2]
        # Now v1 and v2 are Lorentz orthogonal, and |v1| = -1, |v2|=1
        # That is, v1 is unit timelike and v2 is unit spacelike.
        # This means that cosh(x)*v1 + sinh(x)*v2 is unit timelike.
        hyperbola = cosh(x)*v1 + sinh(x)*v2
        endtime = arcsinh(v2_ldot_u2)
        from sage.plot.plot3d.all import parametric_plot3d
        pic = parametric_plot3d(hyperbola,(x,0, endtime),**graphics_options)
        if show_hyperboloid:
            bd_pic = self._model.get_background_graphic()
            pic = bd_pic + pic
        return pic

