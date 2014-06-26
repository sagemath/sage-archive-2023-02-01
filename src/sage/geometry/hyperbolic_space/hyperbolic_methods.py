r"""
Hyperbolic Methods

This module should not be used directly by users.  It is provided for
developers of sage.

This module implements computational methods for some models of
hyperbolic space.  The methods may operate on points, geodesics, or
isometries of hyperbolic space.  However, instead of taking
HyperbolicPoint, HyperbolicGeodesic, or HyperbolicIsometry objects as
input, they instead take the coordinates of points, the endpoints of
geodesics, or matrices.  Similarly, they output coordinates or matrices
rather than Hyperbolic objects.

The methods are factored out of the :class:`HyperbolicPoint`,
:class:`HyperbolicGeodesic`, and :class:`HyperbolicIsometry` classes
to allow the implementation of additional
models of hyperbolic space with minimal work.  For example, to implement
a model of 2-dimensional hyperbolic space, as long as one provides an
isometry of that model with the upper half plane, one can use the upper
half plane methods to do computations.  This prevents, for example,
having to work out an efficient algorithm for computing the midpoint of
a geodesic in every model.  Isometries are implemented in the
HyperbolicModel module, and it is primarily in that model that new code
must be added to implement a new model of hyperbolic space.

In practice, all of the current models of 2 dimensional hyperbolic space
use the upper half plane model for their computations.  This can lead to
some problems, such as long coordinate strings for symbolic points.  For
example, the vector ``(1, 0, sqrt(2))`` defines a point in the hyperboloid
model.  Performing mapping this point to the upper half plane and
performing computations there may return with vector whose components
are unsimplified strings have several ``sqrt(2)``'s.  Presently, this
drawback is outweighed by the rapidity with which new models can be
implemented.

AUTHORS:

- Greg Laun (2013): Refactoring, rewrites, all docstrings.
- Rania Amer (2011): some UHP and PD methods.
- Jean-Philippe Burelle (2011): some UHP and PD methods.
- Zach Groton (2011): some UHP and PD methods.
- Greg Laun (2011): some UHP and PD methods.
- Jeremy Lent (2011): some UHP and PD methods.
- Leila Vaden (2011): some UHP and PD methods.
- Derrick Wigglesworth (2011): some UHP and PD methods.
- Bill Goldman (2011): many UHP and PD methods, implemented in Mathematica.
"""

#***********************************************************************
#
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.lazy_import import lazy_import
from sage.symbolic.pynac import I
from sage.functions.all import exp, cos, sin, arccosh, arccos, sqrt, sign
from sage.functions.all import imag, real
from sage.matrix.all import matrix
from sage.rings.all import Integer, RR, RDF, infinity
from sage.geometry.hyperbolic_space.hyperbolic_constants import EPSILON
from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModel
from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
from sage.geometry.hyperbolic_space.hyperbolic_model import mobius_transform


class HyperbolicAbstractMethods(UniqueRepresentation):
    r"""
    The abstract base class for hyperbolic methods.  Primarily serving
    as a list of methods that must be implemented.
    """
    HModel = HyperbolicModel

    @classmethod
    def model(cls):
        r"""
        Return the class of the underlying hyperbolic model.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelUHP'>
        """
        return cls.HModel

    @classmethod
    def model_name(cls):
        r"""
        Return the short name of the underlying hyperbolic model.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.model_name()
            'UHP'
        """
        return cls.HModel.short_name


class HyperbolicMethodsUHP (HyperbolicAbstractMethods):
    r"""
    Hyperbolic methods for the UHP model of hyperbolic space.
    """
    HModel = HyperbolicModelUHP

#################
# Point Methods #
#################

    @classmethod
    def point_dist(cls, p1, p2):
        r"""
        Compute the distance between two points in the Upper Half Plane
        using the hyperbolic metric.

        EXAMPLES::

           sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
           sage: HyperbolicMethodsUHP.point_dist(4.0*I,I)
           1.38629436111989
        """
        cls.model().point_test(p1)
        cls.model().point_test(p2)
        num = (real(p2) - real(p1))**2 + (imag(p2) - imag(p1))**2
        denom = 2*imag(p1)*imag(p2)
        if denom == 0:
            return infinity
        else:
            return arccosh(1 + num/denom)

    @classmethod
    def symmetry_in(cls, p):
        r"""
        Return the involutary isometry fixing the given point.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.symmetry_in(3 + 2*I)
            [  3/2 -13/2]
            [  1/2  -3/2]
        """
        cls.model().point_test(p)
        x, y = real(p), imag(p)
        if y > 0:
            return matrix(2,[x/y,-(x**2/y) - y,1/y,-(x/y)])

    @classmethod
    def random_point(cls, **kwargs):
        r"""
        Return a random point in the upper half
        plane.  The points are uniformly distributed over the rectangle
        `[-10,10] \times [-10i,10i]`.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: p = HyperbolicMethodsUHP.random_point()
            sage: bool((p.imag()) > 0)
            True
        """
        # TODO use **kwargs to allow these to be set
        real_min = -10
        real_max = 10
        imag_min = 0
        imag_max = 10
        return RR.random_element(min = real_min ,max=real_max) +\
            I*RR.random_element(min = imag_min,max = imag_max)

####################
# Geodesic Methods #
####################

    @classmethod
    def boundary_points(cls, p1, p2):
        r"""
        Given two points ``p1`` and ``p2`` in the hyperbolic plane,
        determine the endpoints of the complete hyperbolic geodesic
        through ``p1`` and ``p2``.

        INPUT:

        - ``p1``, ``p2`` -- points in the hyperbolic plane.

        OUTPUT:

        - a list of boundary points.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.boundary_points(I, 2*I)
            [0, +Infinity]
            sage: HyperbolicMethodsUHP.boundary_points(1 + I, 2 + 4*I)
            [-sqrt(65) + 9, sqrt(65) + 9]
        """
        if p1 == p2:
            raise ValueError(str(p1)  + " and " + str(p2) +
                             " are not distinct.")
        [x1, x2] = [real(k) for k in [p1,p2]]
        [y1,y2] = [imag(k) for k in [p1,p2]]
        # infinity is the first endpoint, so the other ideal endpoint
        # is just the real part of the second coordinate
        if p1 == infinity:
            return [p1, x2]
        # Same idea as above
        elif p2 == infinity:
            return [x1, p2]
        # We could also have a vertical line with two interior points
        elif x1 == x2:
            return [x1, infinity]
        # Otherwise, we have a semicircular arc in the UHP
        else:
            c = ((x1+x2)*(x2-x1)+(y1+y2)*(y2-y1))/(2*(x2-x1))
            r = sqrt((c - x1)**2 + y1**2)
            return [c-r, c + r]

    @classmethod
    def reflection_in(cls, start, end):
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
        x,y = [real(k) for k in cls.boundary_points(start, end)]
        if x == infinity:
            M = matrix(2,[1,-2*y,0,-1])
        elif y == infinity:
            M = matrix(2,[1,-2*x,0,-1])
        else:
            M = matrix(2,[(x+y)/(y-x),-2*x*y/(y-x),2/(y-x),-(x+y)/(y-x)])
        return M

    @classmethod
    def common_perpendicular(cls, start_1, end_1, start_2, end_2):
        r"""
        Return the unique hyperbolic geodesic perpendicular to two given
        geodesics, if such a geodesic exists.  If none exists, raise a
        ValueError.

        INPUT:

        - ``other`` -- a hyperbolic geodesic in current model.

        OUTPUT:

        - a hyperbolic geodesic.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.common_perpendicular(2, 3, 4, 5)
             [1/2*sqrt(3) + 7/2, -1/2*sqrt(3) + 7/2]

        It is an error to ask for the common perpendicular of two
        intersecting geodesics::

            sage: HyperbolicMethodsUHP.common_perpendicular(2, 4, 3, infinity)
            Traceback (most recent call last):
            ...
            ValueError: Geodesics intersect. No common perpendicular exists.
        """
        A = cls.reflection_in(start_1, end_1)
        B = cls.reflection_in(start_2, end_2)
        C = A*B
        if cls.classification(C) != 'hyperbolic':
            raise ValueError("Geodesics intersect. " +
                                 "No common perpendicular exists.")
        return cls.fixed_point_set(C)

    @classmethod
    def intersection(cls, start_1, end_1, start_2, end_2):
        r"""
        Return the point of intersection of two complete geodesics
        (if such a point exists).

        INPUT:

        - ``other`` -- a hyperbolic geodesic in the current model.

        OUTPUT:

        - a hyperbolic point.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.intersection(3, 5, 4, 7)
            [2/3*sqrt(-2) + 13/3]

        If the given geodesics do not intersect, the function raises an
        error::

            sage: HyperbolicMethodsUHP.intersection(4, 5, 5, 7)
            Traceback (most recent call last):
            ...
            ValueError: Geodesics don't intersect.

        If the given geodesics are identical, return that
        geodesic::

            sage: HyperbolicMethodsUHP.intersection(4 + I, 18*I, 4 + I, 18*I)
            [-1/8*sqrt(114985) - 307/8, 1/8*sqrt(114985) - 307/8]
        """
        start_1, end_1 = sorted(cls.boundary_points(start_1, end_1))
        start_2, end_2 = sorted(cls.boundary_points(start_2, end_2))
        if start_1 == start_2 and end_1 == end_2: # Unoriented geods are same
            return [start_1, end_1]
        if start_1 == start_2:
            return start_1
        elif end_1 == end_2:
            return end_1
        A = cls.reflection_in(start_1, end_1)
        B = cls.reflection_in(start_2, end_2)
        C = A*B
        if cls.classification(C) in ['hyperbolic', 'parabolic']:
            raise ValueError("Geodesics don't intersect.")
        return cls.fixed_point_set(C)

    @classmethod
    def perpendicular_bisector(cls, start, end):
        r"""
        Return the perpendicular bisector of the hyperbolic geodesic
        with endpoints start and end if that geodesic has finite length.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: a, b = [HyperbolicMethodsUHP.random_point() for k in range(2)]
            sage: g = [a, b]
            sage: h = HyperbolicMethodsUHP.perpendicular_bisector(*g)
            sage: bool(HyperbolicMethodsUHP.intersection(*(g + h))[0] - HyperbolicMethodsUHP.midpoint(*g) < 10**-9)
            True

        Complete geodesics cannot be bisected::

            sage: HyperbolicMethodsUHP.perpendicular_bisector(0, 1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular
        """
        from sage.symbolic.constants import pi
        from sage.functions.log import exp
        start, end = sorted((start, end))
        d = cls.point_dist(start, end)/2
        end_1, end_2 = cls.boundary_points(start, end)
        T = (matrix(2,[exp(d/2),0,0,exp(-d/2)])*
             matrix(2,[cos(pi/4),-sin(pi/4),sin(pi/4),cos(pi/4)]))
        S= cls._to_std_geod(end_1, start, end_2)
        H = S.inverse()*T*S
        return [mobius_transform(H ,k) for k in [end_1, end_2]]

    @classmethod
    def midpoint(cls, start, end):
        r"""
        Return the (hyperbolic) midpoint of a hyperbolic line segment.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: a, b = [HyperbolicMethodsUHP.random_point() for k in range(2)]
            sage: g = [a, b]
            sage: m = HyperbolicMethodsUHP.midpoint(*g)
            sage: d1 =HyperbolicMethodsUHP.point_dist(m, g[0])
            sage: d2 = HyperbolicMethodsUHP.point_dist(m, g[1])
            sage: bool(abs(d1 - d2) < 10**-9)
            True

        Complete geodesics have no midpoint::

            sage: HyperbolicMethodsUHP.midpoint(0, 2)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular
        """
        half_dist = cls.point_dist(start, end)/2
        end_1,end_2 = cls.boundary_points(start, end)
        S = cls._to_std_geod(end_1, start , end_2)
        T = matrix(2,[exp(half_dist), 0, 0, 1])
        M = S.inverse()*T*S
        if ((real(start - end) < EPSILON) or
                (abs(real(start - end)) < EPSILON and
                 imag(start - end) < EPSILON)):
            end_p = start
        else:
            end_p = end
        end_p = mobius_transform (M, end_p)
        return end_p

    @classmethod
    def geod_dist_from_point(cls, start, end, p):
        r"""
        Return the hyperbolic distance from a given hyperbolic geodesic
        and a hyperbolic point.

        INPUT:

        - ``other`` -- a hyperbolic point in the same model.

        OUTPUT:

        - the hyperbolic distance.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.geod_dist_from_point(2, 2 + I, I)
            arccosh(1/10*sqrt(5)*((sqrt(5) - 1)^2 + 4) + 1)

        If p is a boundary point, the distance is infinity::

            sage: HyperbolicMethodsUHP.geod_dist_from_point(2, 2 + I, 5)
            +Infinity
        """
        # Here is the trick for computing distance to a geodesic:
        # find an isometry mapping the geodesic to the geodesic between
        # 0 and infinity (so corresponding to the line imag(z) = 0.
        # then any complex number is r exp(i*theta) in polar coordinates.
        # the mutual perpendicular between this point and imag(z) = 0
        # intersects imag(z) = 0 at ri.  So we calculate the distance
        # between r exp(i*theta) and ri after we transform the original
        # point.
        (bd_1, bd_2) = cls.boundary_points(start, end)
        if bd_1 + bd_2 != infinity:
            # Not a straight line
            # Map the endpoints to 0 and infinity and the midpoint to 1.
            T = cls._crossratio_matrix(bd_1, (bd_1 + bd_2)/2, bd_2)
        else:
            # Is a straight line
            # Map the endpoints to 0 and infinity and another endpoint
            # to 1
            T = cls._crossratio_matrix(bd_1, bd_1 + 1, bd_2)
        x = mobius_transform(T, p)
        return cls.point_dist(x, abs(x)*I)

    @classmethod
    def uncomplete(cls, start, end):
        r"""
        Return starting and ending points of a geodesic whose completion
        is the the geodesic starting at ``start`` and ending at ``end``.

        INPUT:

        - ``start`` -- a real number or infinity.
        - ``end`` -- a real number or infinity

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.uncomplete(0,1)
            [3/10*I + 1/10, 3/10*I + 9/10]
       """
        if start + end == infinity:
            p = min(real(start), real(end)) + I
        else:
            center = (real(start) + real(end))/Integer(2)
            radius = abs(real(end) - center)
            p = center + radius*I
        A = cls._to_std_geod(start, p, end).inverse()
        p1 = mobius_transform(A, I/Integer(3))
        p2 = mobius_transform(A, 3*I)
        return [p1, p2]

    @classmethod
    def angle(cls, start_1, end_1, start_2, end_2):
        r"""
        Return the angle  between any two given completed geodesics if
        they intersect.

        INPUT:

        -``other`` -- a hyperbolic geodesic in the UHP model

        OUTPUT:

        - The angle in radians between the two given geodesics

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: numerical_approx(HyperbolicMethodsUHP.angle(2, 4, 3, 3 + I))
            1.57079632679490

        If the geodesics are identical, return angle 0::

            sage: HyperbolicMethodsUHP.angle(2, 4, 2, 4)
            0
        """
        (p_1,p_2) = sorted(cls.boundary_points(start_1, end_1))
        (q_1,q_2) = sorted(cls.boundary_points(start_2, end_2))
        # if the geodesics are equal, the angle between them is 0
        if (abs(p_1 - q_1) < EPSILON \
                and abs(p_2 - q_2) < EPSILON):
            return 0
        elif p_2 != infinity: # geodesic not a straight line
            # So we send it to the geodesic with endpoints [0, oo]
            T = cls._crossratio_matrix(p_1, (p_1+p_2)/2, p_2)
        else:
            # geodesic is a straight line, so we send it to the geodesic
            # with endpoints [0,oo]
            T = cls._crossratio_matrix(p_1, p_1 +1, p_2)
        # b_1 and b_2 are the endpoints of the image of other
        b_1, b_2 = [mobius_transform(T, k) for k in [q_1, q_2]]
        # If other is now a straight line...
        if (b_1 == infinity or b_2 == infinity):
            # then since they intersect, they are equal
            return 0
        else:
            return real(arccos((b_1+b_2)/abs(b_2-b_1)))


####################
# Isometry Methods #
####################
    @classmethod
    def orientation_preserving(cls, M):
        r"""
        Return ``True`` if ``self`` is orientation preserving and ``False``
        otherwise.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: A = identity_matrix(2)
            sage: HyperbolicMethodsUHP.orientation_preserving(A)
            True
            sage: B = matrix(2,[0,1,1,0])
            sage: HyperbolicMethodsUHP.orientation_preserving(B)
            False
        """
        return bool(M.det() > 0)

    @classmethod
    def classification(cls, M):
        r"""
        Classify the hyperbolic isometry as elliptic, parabolic, or
        hyperbolic.

        A hyperbolic isometry fixes two points on the boundary of
        hyperbolic space, a parabolic isometry fixes one point on the
        boundary of hyperbolic space, and an elliptic isometry fixes
        no points.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.classification(identity_matrix(2))
            'identity'

            sage: HyperbolicMethodsUHP.classification(4*identity_matrix(2))
            'identity'

            sage: HyperbolicMethodsUHP.classification(matrix(2,[2,0,0,1/2]))
            'hyperbolic'

            sage: HyperbolicMethodsUHP.classification(matrix(2, [0, 3, -1/3, 6]))
            'hyperbolic'

            sage: HyperbolicMethodsUHP.classification(matrix(2,[1,1,0,1]))
            'parabolic'

            sage: HyperbolicMethodsUHP.classification(matrix(2,[-1,0,0,1]))
            'reflection'
        """
        A = M.n()
        A = A/(abs(A.det()).sqrt())
        tau = abs(A.trace())
        a = A.list()
        if A.det() > 0:
            tf = bool((a[0] - 1)**2 + a[1]**2 + a[2]**2 + (a[3] - 1)**2  <
                      EPSILON)
            tf = bool((a[0] - 1)**2 + a[1]**2 + a[2]**2 + (a[3] - 1)**2  <
                      EPSILON)
            if tf:
                return 'identity'
            if tau -  2 < -EPSILON:
                return 'elliptic'
            elif tau -2  > -EPSILON and tau -  2 < EPSILON:
                return 'parabolic'
            elif tau - 2  > EPSILON:
                return 'hyperbolic'
            else:
                raise ValueError("Something went wrong with classification!" +
                                 "Trace is " + str(A.trace()))
        else:  #The isometry reverses orientation.
            if tau < EPSILON:
                return 'reflection'
            else:
                return 'orientation-reversing hyperbolic'

    @classmethod
    def translation_length(cls, M):
        r"""
        For hyperbolic elements, return the translation length.
        Otherwise, raise a ValueError.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.translation_length(matrix(2,[2,0,0,1/2]))
            2*arccosh(5/4)

        ::

            sage: H = HyperbolicMethodsUHP.isometry_from_fixed_points(-1,1)
            sage: p = exp(i*7*pi/8)
            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import mobius_transform
            sage: Hp = mobius_transform(H, p)
            sage: bool((HyperbolicMethodsUHP.point_dist(p, Hp) - HyperbolicMethodsUHP.translation_length(H)) < 10**-9)
            True
        """
        d = sqrt(M.det()**2)
        tau = sqrt((M/sqrt(d)).trace()**2)
        if cls.classification(M) in ['hyperbolic',
                                     'oriention-reversing hyperbolic']:
            return 2*arccosh(tau/2)
            raise TypeError("Translation length is only defined for hyperbolic"
                            " transformations.")

    @classmethod
    def isometry_from_fixed_points(cls, repel, attract):
        r"""
        Given two fixed points ``repel`` and ``attract`` as complex
        numbers return a hyperbolic isometry with ``repel`` as repelling
        fixed point and ``attract`` as attracting fixed point.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.isometry_from_fixed_points(2 + I, 3 + I)
            Traceback (most recent call last):
            ...
            ValueError: fixed points of hyperbolic elements must be ideal

            sage: HyperbolicMethodsUHP.isometry_from_fixed_points(2, 0)
            [  -1    0]
            [-1/3 -1/3]
        """
        if imag(repel) + imag(attract) > EPSILON:
            raise ValueError("fixed points of hyperbolic elements must be ideal")
        repel = real(repel)
        attract = real(attract)
        if repel == infinity:
            A = cls._mobius_sending([infinity, attract, attract + 1],
                                    [infinity, attract, attract + 2])
        elif attract == infinity:
            A = cls._mobius_sending([repel, infinity, repel + 1],
                                    [repel, infinity, repel + 2])
        else:
            A = cls._mobius_sending([repel, attract, infinity],
                                    [repel, attract, max(repel, attract) + 1])
        return A


    @classmethod
    def fixed_point_set(cls, M):
        r"""
        Return the a list containing the fixed point set of
        orientation-preserving isometries.

        OUTPUT:

        - a list of hyperbolic points.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: H = matrix(2, [-2/3,-1/3,-1/3,-2/3])
            sage: (p1,p2) = HyperbolicMethodsUHP.fixed_point_set(H)
            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import mobius_transform
            sage: bool(mobius_transform(H, p1) == p1)
            True
            sage: bool(mobius_transform(H, p2) == p2)
            True

            sage: HyperbolicMethodsUHP.fixed_point_set(identity_matrix(2))
            Traceback (most recent call last):
            ...
            ValueError: the identity transformation fixes the entire hyperbolic plane
        """
        d = sqrt(M.det()**2)
        M = M/sqrt(d)
        tau = M.trace()**2
        M_cls = cls.classification(M)
        if M_cls == 'identity':
            raise ValueError("the identity transformation fixes the entire "
                             "hyperbolic plane")
        if M_cls == 'parabolic':
            if abs(M[1,0]) < EPSILON:
                return [infinity]
            else:
                # boundary point
                return [(M[0,0] - M[1,1])/(2*M[1,0])]
        elif M_cls=='elliptic':
            d = sqrt(tau - 4)
            return [(M[0,0] - M[1,1] + sign(M[1,0])*d)/(2*M[1,0])]
        elif M_cls == 'hyperbolic':
            if M[1,0]!= 0: #if the isometry doesn't fix infinity
                d = sqrt(tau - 4)
                p_1 = (M[0,0] - M[1,1]+d) / (2*M[1,0])
                p_2 = (M[0,0] - M[1,1]-d) / (2*M[1,0])
                return [p_1, p_2]
            else: #else, it fixes infinity.
                p_1 = M[0,1]/(M[1,1]-M[0,0])
                p_2 = infinity
                return [p_1, p_2]
        else:
            # raise NotImplementedError("fixed point set not implemented for"
            #                           " isometries of type {0}".format(M_cls))
            try:
                p, q = [M.eigenvectors_right()[k][1][0] for k in range(2)]
            except (IndexError):
                M = M.change_ring(RDF)
                p, q = [M.eigenvectors_right()[k][1][0] for k in range(2)]
            if p[1] == 0:
                p = infinity
            else:
                p = p[0]/p[1]
            if q[1] == 0:
                q = infinity
            else:
                q = q[0]/q[1]
            pts = [p, q]
            return [k for k in pts if imag(k) >= 0]


    @classmethod
    def repelling_fixed_point(cls, M):
        r"""
        For a hyperbolic isometry, return the attracting fixed point.
        Otherwise raise a ValueError.

        OUTPUT:

        - a hyperbolic point.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: A = matrix(2,[4,0,0,1/4])
            sage: HyperbolicMethodsUHP.repelling_fixed_point(A)
            0
        """
        if cls.classification(M) not in ['hyperbolic',
                                         'orientation-reversing hyperbolic']:
            raise ValueError("Repelling fixed point is defined only" +
                             "for hyperbolic isometries.")
        v = M.eigenmatrix_right()[1].column(1)
        if v[1] == 0:
            return infinity
        return v[0]/v[1]

    @classmethod
    def attracting_fixed_point(cls, M):
        r"""
        For a hyperbolic isometry, return the attracting fixed point.
        Otherwise raise a ValueError.

        OUTPUT:

        - a hyperbolic point.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: A = matrix(2,[4,0,0,1/4])
            sage: HyperbolicMethodsUHP.attracting_fixed_point(A)
            +Infinity
        """
        if cls.classification(M) not in \
                ['hyperbolic', 'orientation-reversing hyperbolic']:
            raise ValueError("Attracting fixed point is defined only" +
                             "for hyperbolic isometries.")
        v = M.eigenmatrix_right()[1].column(0)
        if v[1] == 0:
            return infinity
        return v[0]/v[1]

    @classmethod
    def random_isometry(cls, preserve_orientation = True, **kwargs):
        r"""
        Return a random isometry in the Upper Half Plane model.

        INPUT:

        - ``preserve_orientation`` -- if ``True`` return an orientation-preserving isometry.

        OUTPUT:

        - a hyperbolic isometry.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: A = HyperbolicMethodsUHP.random_isometry()
            sage: B = HyperbolicMethodsUHP.random_isometry(preserve_orientation=False)
            sage: B.det() < 0
            True
        """
        [a,b,c,d] = [RR.random_element() for k in range(4)]
        while abs(a*d - b*c) < EPSILON:
            [a,b,c,d] = [RR.random_element() for k in range(4)]
        M = matrix(RDF, 2,[a,b,c,d])
        M = M/(M.det()).abs().sqrt()
        if M.det() > 0:
            if preserve_orientation:
                return M
            else:
                return M*matrix(2,[0,1,1,0])
        else:
            if preserve_orientation:
                return M*matrix(2,[0,1,1,0])
            else:
                return M

###################
# Helping Methods #
###################


    @classmethod
    def _to_std_geod(cls, start, p, end):
        r"""
        Given the points of a geodesic in hyperbolic space, return the
        hyperbolic isometry that sends that geodesic to the geodesic
        through 0 and infinity that also sends the point ``p`` to ``I``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import mobius_transform
            sage: (p_1, p_2, p_3) = [HyperbolicMethodsUHP.random_point() for k in range(3)]
            sage: A = HyperbolicMethodsUHP._to_std_geod(p_1, p_2, p_3)
            sage: bool(abs(mobius_transform(A, p_1)) < 10**-9)
            True
            sage: bool(abs(mobius_transform(A, p_2) - I) < 10**-9)
            True
            sage: bool(mobius_transform(A, p_3) == infinity)
            True
        """
        B = matrix(2, [1, 0, 0, -I])
        return B*cls._crossratio_matrix(start, p, end)

    @classmethod
    def _crossratio_matrix(cls, p_0, p_1, p_2):
        r"""
        Given three points (the list `p`) in `\mathbb{CP}^{1}' in affine
        coordinates, return the linear fractional transformation taking
        the elements of `p` to `0`,`1`, and `\infty'.

        INPUT:

        - A list of three distinct elements of three distinct elements of `\mathbb{CP}^1` in affine coordinates.  That is, each element must be a complex number, `\infty`, or symbolic.

        OUTPUT:

        - An element of '\GL(2,\CC)`.


        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import mobius_transform
            sage: (p_1, p_2, p_3) = [HyperbolicMethodsUHP.random_point() for k in range(3)]
            sage: A = HyperbolicMethodsUHP._crossratio_matrix(p_1, p_2, p_3)
            sage: bool(abs(mobius_transform(A, p_1) < 10**-9))
            True
            sage: bool(abs(mobius_transform(A, p_2) - 1) < 10**-9)
            True
            sage: bool(mobius_transform(A, p_3) == infinity)
            True
            sage: (x,y,z) = var('x,y,z');  HyperbolicMethodsUHP._crossratio_matrix(x,y,z)
            [     y - z -x*(y - z)]
            [    -x + y  (x - y)*z]
        """
        if p_0 == infinity:
            return matrix(2,[0,-(p_1 - p_2), -1, p_2])
        elif p_1 == infinity:
            return matrix(2,[1, -p_0,1,-p_2])
        elif p_2 == infinity:
            return matrix(2,[1,-p_0, 0, p_1 - p_0])
        else:
            return matrix(2,[p_1 - p_2, (p_1 - p_2)*(-p_0), p_1 - p_0,
                          ( p_1 - p_0)*(-p_2)])


    @classmethod
    def _mobius_sending(cls, list1, list2):
        r"""
        Given two lists ``list1``, ``list2``, of three points each in
        `\mathbb{CP}^1`, return the linear fractional transformation
        taking the points in ``list1`` to the points in ``list2``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import mobius_transform
            sage: bool(abs(mobius_transform(HyperbolicMethodsUHP._mobius_sending([1,2,infinity],[3 - I, 5*I,-12]),1) - 3 + I) < 10^-4)
            True
            sage: bool(abs(mobius_transform(HyperbolicMethodsUHP._mobius_sending([1,2,infinity],[3 - I, 5*I,-12]),2) - 5*I) < 10^-4)
            True
            sage: bool(abs(mobius_transform(HyperbolicMethodsUHP._mobius_sending([1,2,infinity],[3 - I, 5*I,-12]),infinity) + 12) < 10^-4)
            True
        """
        if len(list1) != 3  or len(list2) != 3:
            raise TypeError("mobius_sending requires each list to be three"
                            "points long.")
        pl = list1 + list2
        z = pl[0:3]
        w = pl[3:6]
        A = cls._crossratio_matrix(z[0],z[1],z[2])
        B = cls._crossratio_matrix(w[0],w[1],w[2])
        return B.inverse()*A
