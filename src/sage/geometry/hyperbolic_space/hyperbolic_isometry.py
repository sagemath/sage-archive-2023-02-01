r"""
Hyperbolic Isometries

This module implements the abstract base class for isometries of
hyperbolic space of arbitrary dimension.  It also contains the
implementations for specific models of hyperbolic geometry.

The isometry groups of all implemented models are either matrix Lie
groups or are doubly covered by matrix Lie groups.  As such, the
isometry constructor takes a matrix as input.  However, since the
isometries themselves may not be matrices, quantities like the trace
and determinant are not directly accessible from this class.

AUTHORS:

- Greg Laun (2013): initial version

EXAMPLES:

We can construct isometries in the upper half plane model, abbreviated
UHP for convenience::

    sage: UHP.isometry(matrix(2,[1,2,3,4]))
    Isometry in UHP
    [1 2]
    [3 4]
    sage: A = UHP.isometry(matrix(2,[0,1,1,0]))
    sage: A.inverse()
    Isometry in UHP
    [0 1]
    [1 0]
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
from sage.misc.lazy_import import lazy_import
from sage.misc.lazy_attribute import lazy_attribute
lazy_import('sage.modules.free_module_element', 'vector')
from sage.rings.infinity import infinity
from sage.geometry.hyperbolic_space.hyperbolic_constants import EPSILON
from sage.misc.latex import latex
from sage.rings.all import CC
from sage.functions.other import real, imag

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_factory', 'HyperbolicAbstractFactory')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_methods',
            ['HyperbolicAbstractMethods', 'HyperbolicMethodsUHP'])

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_factory',
            ['HyperbolicFactoryUHP', 'HyperbolicFactoryPD',
             'HyperbolicFactoryKM', 'HyperbolicFactoryHM'])

class HyperbolicIsometry(SageObject):
    r"""
    Abstract base class for hyperbolic isometries.  This class should
    never be instantiated.

    INPUT:

    - ``A`` -- a matrix representing a hyperbolic isometry in the
      appropriate model

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
        sage: A = HyperbolicIsometryUHP(identity_matrix(2))
        sage: B = HyperbolicIsometryHM(identity_matrix(3))
    """
    HFactory = HyperbolicAbstractFactory
    HMethods = HyperbolicAbstractMethods

    #####################
    # "Private" Methods #
    #####################

    def __init__(self, A):
        r"""
        See `HyperbolicIsometry` for full documentation.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: HyperbolicIsometryUHP(matrix(2, [0,1,-1,0]))
            Isometry in UHP
            [ 0  1]
            [-1  0]
        """
        self._model = self.HFactory.get_model()
        self._model.isometry_test(A)
        self._matrix = A

    @lazy_attribute
    def _cached_matrix(self):
        r"""
        The representation of the current isometry used for
        calculations.  For example, if the current model uses the
        HyperbolicMethodsUHP class, then _cached_matrix will hold the
        SL(2,`\Bold{R}`) representation of self.matrix().

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: A = HyperbolicIsometryHM(identity_matrix(3))
            sage: A._cached_matrix
            [1 0]
            [0 1]
        """
        return self.model().isometry_to_model(self.matrix(),
                                              self.HMethods.model().short_name)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - a string

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: HyperbolicIsometryUHP(identity_matrix(2))
            Isometry in UHP
            [1 0]
            [0 1]
        """
        return "Isometry in {0}\n{1}".format(self.model_name(), self.matrix())

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: A = UHP.isometry(identity_matrix(2))
            sage: B = HM.isometry(identity_matrix(3))
            sage: latex(A)
            \pm \left(\begin{array}{rr}
            1 & 0 \\
            0 & 1
            \end{array}\right)

            sage: latex(B)
            \left(\begin{array}{rrr}
            1 & 0 & 0 \\
            0 & 1 & 0 \\
            0 & 0 & 1
            \end{array}\right)
        """
        if self.model().isometry_group_is_projective:
            return "\pm " + latex(self.matrix())
        else:
            return latex(self.matrix())

    def __eq__(self, other):
        r"""
        Return ``True`` if the isometries are the same and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: A = HyperbolicIsometryUHP(identity_matrix(2))
            sage: B = HyperbolicIsometryUHP(-identity_matrix(2))
            sage: A == B
            True
        """
        pos_matrix = bool(abs(self.matrix() - other.matrix()) < EPSILON)
        neg_matrix = bool(abs(self.matrix() + other.matrix()) < EPSILON)
        if self._model.isometry_group_is_projective:
            return self.model() == other.model() and (pos_matrix or neg_matrix)
        else:
            return self.model_name() == other.model_name() and pos_matrix

    def __pow__(self, n):
        r"""
        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: A = HyperbolicIsometryUHP(matrix(2,[3,1,2,1]))
            sage: A**3
            Isometry in UHP
            [41 15]
            [30 11]
        """
        return self.__class__(self.matrix()**n)

    def __mul__(self, other):
        r"""
        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: A = HyperbolicIsometryUHP(Matrix(2,[5,2,1,2]))
            sage: B = HyperbolicIsometryUHP(Matrix(2,[3,1,1,2]))
            sage: B*A
            Isometry in UHP
            [16  8]
            [ 7  6]
            sage: A = HyperbolicIsometryUHP(Matrix(2,[5,2,1,2]))
            sage: p = UHP.point(2 + I)
            sage: A*p
            Point in UHP 8/17*I + 53/17

            sage: g = UHP.geodesic(2 + I,4 + I)
            sage: A*g
            Geodesic in UHP from 8/17*I + 53/17 to 8/37*I + 137/37

            sage: A = diagonal_matrix([1, -1, 1])
            sage: A = HyperbolicIsometryHM(A)
            sage: A.orientation_preserving()
            False
            sage: p = HM.point((0, 1, sqrt(2)))
            sage: A*p
            Point in HM (0, -1, sqrt(2))
        """
        from sage.geometry.hyperbolic_space.hyperbolic_geodesic import HyperbolicGeodesic
        from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPoint
        if self.model_name() != other.model_name():
            raise TypeError("{0} and {1} are not in the same"
                            "model".format(self, other))
        if isinstance(other, HyperbolicIsometry):
            return type (self) (self.matrix()*other.matrix())
        elif isinstance(other, HyperbolicPoint):
            return self.HFactory.get_point(self.model().isometry_act_on_point(
                self.matrix(), other.coordinates()))
        elif isinstance(other, HyperbolicGeodesic):
            return self.HFactory.get_geodesic(self*other.start(), self*other.end())
        else:
            NotImplementedError("Multiplication is not defined between a "
                                "hyperbolic isometry and {0}".format(other))

    # def __call__ (self, other):
    #     r"""
    #     EXAMPLES::
    #
    #         sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
    #         sage: A = HyperbolicIsometryUHP(Matrix(2,[5,2,1,2]))
    #         sage: p = UHP.point(2 + I)
    #         sage: A(p)
    #         Point in UHP 8/17*I + 53/17.
    #
    #         sage: g = UHP.geodesic(2 + I,4 + I)
    #         sage: A (g)
    #         Geodesic in UHP from 8/17*I + 53/17 to 8/37*I + 137/37.
    #
    #         sage: A = diagonal_matrix([1, -1, 1])
    #         sage: A = HyperbolicIsometryHM(A)
    #         sage: A.orientation_preserving()
    #         False
    #         sage: p = HM.point((0, 1, sqrt(2)))
    #         sage: A(p)
    #         Point in HM (0, -1, sqrt(2)).
    #     """
    #     from sage.geometry.hyperbolic_space.hyperbolic_geodesic import HyperbolicGeodesic
    #     if self.model() != other.model():
    #         raise TypeError("{0} is not in the {1} model.".format(other, self.model_name()))
    #     if isinstance(other, HyperbolicGeodesic):
    #         return self.HFactory.get_geodesic(self(other.start()), self(other.end()))
    #     return self.HFactory.get_point(self.model().isometry_act_on_point(
    #         self.matrix(), other.coordinates()))

    #######################
    # Setters and Getters #
    #######################

    def matrix(self):
        r"""
        Return the matrix of the isometry.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: HyperbolicIsometryUHP(-identity_matrix(2)).matrix()
            [-1  0]
            [ 0 -1]
        """
        return self._matrix

    def inverse(self):
        r"""
        Return the inverse of the isometry ``self``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: A = HyperbolicIsometryUHP(matrix(2,[4,1,3,2]))
            sage: B = A.inverse()
            sage: A*B == HyperbolicIsometryUHP(identity_matrix(2))
            True
        """
        return self.HFactory.get_isometry(self.matrix().inverse())

    @classmethod
    def model(this):
        r"""
        Return the model to which the HyperbolicIsometry belongs.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: HyperbolicIsometryUHP(identity_matrix(2)).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelUHP'>

            sage: HyperbolicIsometryPD(identity_matrix(2)).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelPD'>

            sage: HyperbolicIsometryKM(identity_matrix(3)).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelKM'>

            sage: HyperbolicIsometryHM(identity_matrix(3)).model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelHM'>
        """
        return this.HFactory.get_model()

    @classmethod
    def model_name(this):
        r"""
        Return the short name of the hyperbolic model.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: HyperbolicIsometryUHP(identity_matrix(2)).model_name()
            'UHP'

            sage: HyperbolicIsometryPD(identity_matrix(2)).model_name()
            'PD'

            sage: HyperbolicIsometryKM(identity_matrix(3)).model_name()
            'KM'

            sage: HyperbolicIsometryHM(identity_matrix(3)).model_name()
            'HM'
        """
        return this.model().short_name

    def to_model(self, model_name):
        r"""
        Convert the current object to image in another model.

        INPUT:

        - ``model_name`` -- a string representing the image model.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: HyperbolicIsometryUHP(identity_matrix(2)).to_model('HM')
            Isometry in HM
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        from sage.geometry.hyperbolic_space.model_factory import ModelFactory
        factory = ModelFactory.find_factory(model_name)
        matrix = self.model().isometry_to_model(self.matrix(),  model_name)
        return factory.get_isometry(matrix)

    ###################
    # Boolean Methods #
    ###################

    def orientation_preserving(self):
        r"""
        Return ``True`` if ``self`` is orientation preserving and ``False``
        otherwise.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: A = HyperbolicIsometryUHP(identity_matrix(2))
            sage: A.orientation_preserving()
            True
            sage: B = HyperbolicIsometryUHP(matrix(2,[0,1,1,0]))
            sage: B.orientation_preserving()
            False
        """
        return self.HMethods.orientation_preserving(self._cached_matrix)

    ###################################
    # Methods implemented in HMethods #
    ###################################

    def classification(self):
        r"""
        Classify the hyperbolic isometry as elliptic, parabolic,
        hyperbolic or a reflection.

        A hyperbolic isometry fixes two points on the boundary of
        hyperbolic space, a parabolic isometry fixes one point on the
        boundary of hyperbolic space, and an elliptic isometry fixes no
        points.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: H = HyperbolicIsometryUHP(matrix(2,[2,0,0,1/2]))
            sage: H.classification()
            'hyperbolic'

            sage: P = HyperbolicIsometryUHP(matrix(2,[1,1,0,1]))
            sage: P.classification()
            'parabolic'

            sage: E = HyperbolicIsometryUHP(matrix(2,[-1,0,0,1]))
            sage: E.classification()
            'reflection'
        """
        return self.HMethods.classification(self._cached_matrix)

    def translation_length(self):
        r"""
        For hyperbolic elements, return the translation length;
        otherwise, raise a ``ValueError``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: H = HyperbolicIsometryUHP(matrix(2,[2,0,0,1/2]))
            sage: H.translation_length()
            2*arccosh(5/4)

        ::

            sage: f_1 = UHP.point(-1)
            sage: f_2 = UHP.point(1)
            sage: H = HyperbolicIsometryUHP.isometry_from_fixed_points(f_1, f_2)
            sage: p = UHP.point(exp(i*7*pi/8))
            sage: bool((p.dist(H*p) - H.translation_length()) < 10**-9)
            True
        """
        return self.HMethods.translation_length(self._cached_matrix)

    def axis(self, **graphics_options):
        r"""
        For a hyperbolic isometry, return the axis of the
        transformation; otherwise raise a ``ValueError``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: H = HyperbolicIsometryUHP(matrix(2,[2,0,0,1/2]))
            sage: H.axis()
            Geodesic in UHP from 0 to +Infinity

            It is an error to call this function on an isometry that is
            not hyperbolic::

            sage: P = HyperbolicIsometryUHP(matrix(2,[1,4,0,1]))
            sage: P.axis()
            Traceback (most recent call last):
            ...
            ValueError: the isometry is not hyperbolic: axis is undefined.
        """
        if self.classification() not in (
                ['hyperbolic', 'orientation-reversing hyperbolic']):
            raise ValueError("the isometry is not hyperbolic: axis is"
                             " undefined.")
        return self.HFactory.get_geodesic(*self.fixed_point_set())

    def fixed_point_set(self, **graphics_options):
        r"""
        Return the a list containing the fixed point set of orientation-
        preserving isometries.

        OUTPUT:

        - a list of hyperbolic points

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: H = HyperbolicIsometryUHP(matrix(2, [-2/3,-1/3,-1/3,-2/3]))
            sage: (p1,p2) = H.fixed_point_set()
            sage: H*p1 == p1
            True
            sage: H*p2 == p2
            True
            sage: A = HyperbolicIsometryUHP(matrix(2,[0,1,1,0]))
            sage: A.orientation_preserving()
            False
            sage: A.fixed_point_set()
            [Boundary point in UHP 1, Boundary point in UHP -1]

       ::

            sage: B = HyperbolicIsometryUHP(identity_matrix(2))
            sage: B.fixed_point_set()
            Traceback (most recent call last):
            ...
            ValueError: the identity transformation fixes the entire hyperbolic plane
        """
        pts = self.HMethods.fixed_point_set(self._cached_matrix)
        pts =  [self.HMethods.model().point_to_model(k, self.model_name()) for\
                         k in pts]
        return [self.HFactory.get_point(k, **graphics_options) for k in pts]

    def fixed_geodesic(self, **graphics_options):
        r"""
        If ``self`` is a reflection in a geodesic, return that geodesic.

        EXAMPLES::

            sage: A = UHP.isometry(matrix(2, [0, 1, 1, 0]))
            sage: A.fixed_geodesic()
            Geodesic in UHP from 1 to -1
        """
        fps = self.HMethods.fixed_point_set(self._cached_matrix)
        if len(fps) < 2:
            raise ValueError("Isometries of type"
                             " {0}".format(self.classification())
                             + " don't fix geodesics.")
        from sage.geometry.hyperbolic_space.model_factory import ModelFactory
        fact = ModelFactory.find_factory(self.HMethods.model_name())
        geod = fact.get_geodesic(fps[0], fps[1])
        return geod.to_model(self.model_name())

    def repelling_fixed_point(self, **graphics_options):
        r"""
        For a hyperbolic isometry, return the attracting fixed point;
        otherwise raise a ``ValueError``.

        OUTPUT:

        - a hyperbolic point

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: A = HyperbolicIsometryUHP(Matrix(2,[4,0,0,1/4]))
            sage: A.repelling_fixed_point()
            Boundary point in UHP 0
        """
        fp = self.HMethods.repelling_fixed_point(self._cached_matrix)
        fp = self.HMethods.model().point_to_model(fp, self.model_name())
        return self.HFactory.get_point(fp)

    def attracting_fixed_point(self, **graphics_options):
        r"""
        For a hyperbolic isometry, return the attracting fixed point;
        otherwise raise a `ValueError``.

        OUTPUT:

        - a hyperbolic point

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: A = HyperbolicIsometryUHP(Matrix(2,[4,0,0,1/4]))
            sage: A.attracting_fixed_point()
            Boundary point in UHP +Infinity
        """
        fp = self.HMethods.attracting_fixed_point(self._cached_matrix)
        fp = self.HMethods.model().point_to_model(fp, self.model_name())
        return self.HFactory.get_point(fp)

    @classmethod
    def isometry_from_fixed_points(cls, repel, attract):
        r"""
        Given two fixed points ``repel`` and ``attract`` as hyperbolic
        points return a hyperbolic isometry with ``repel`` as repelling
        fixed point and ``attract`` as attracting fixed point.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: p, q = [UHP.point(k) for k in [2 + I, 3 + I]]
            sage: HyperbolicIsometryUHP.isometry_from_fixed_points(p, q)
            Traceback (most recent call last):
            ...
            ValueError: fixed points of hyperbolic elements must be ideal

            sage: p, q = [UHP.point(k) for k in [2, 0]]
            sage: HyperbolicIsometryUHP.isometry_from_fixed_points(p, q)
            Isometry in UHP
            [  -1    0]
            [-1/3 -1/3]
        """
        try:
            A = cls.HMethods.isometry_from_fixed_points(repel._cached_coordinates,
                                                   attract._cached_coordinates)
            A = cls.HMethods.model().isometry_to_model(A, cls.model_name())
            return cls.HFactory.get_isometry(A)
        except(AttributeError):
            repel = cls.HFactory.get_point(repel)
            attract = cls.HFactory.get_point(attract)
            return cls.isometry_from_fixed_points(repel, attract)

    @classmethod
    def random_element(cls, preserve_orientation = True, **kwargs):
        r"""
        Return a random isometry in the Upper Half Plane model.

        INPUT:

        - ``preserve_orientation`` -- if ``True`` return an
          orientation-preserving isometry

        OUTPUT:

        - a hyperbolic isometry

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import *
            sage: A = HyperbolicIsometryUHP.random_element()
            sage: B = HyperbolicIsometryUHP.random_element(preserve_orientation=False)
            sage: B.orientation_preserving()
            False
        """
        A = cls.HMethods.random_isometry(preserve_orientation, **kwargs)
        A = cls.HMethods.model().isometry_to_model(A, cls.model_name())
        return cls.HFactory.get_isometry(A)


class HyperbolicIsometryUHP(HyperbolicIsometry):
    r"""
    Create a hyperbolic isometry in the UHP model.

    INPUT:

    - a matrix in `GL(2, \RR)`

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import HyperbolicIsometryUHP
        sage: A = HyperbolicIsometryUHP(identity_matrix(2))
    """

    HFactory = HyperbolicFactoryUHP
    HMethods = HyperbolicMethodsUHP

class HyperbolicIsometryPD(HyperbolicIsometry):
    r"""
    Create a hyperbolic isometry in the PD model.

    INPUT:

    - a matrix in `PU(1,1)`

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import HyperbolicIsometryPD
        sage: A = HyperbolicIsometryPD(identity_matrix(2))
    """
    HFactory = HyperbolicFactoryPD
    HMethods = HyperbolicMethodsUHP

class HyperbolicIsometryKM(HyperbolicIsometry):
    r"""
    Create a hyperbolic isometry in the KM model.

    INPUT:

    - a matrix in `SO(2,1)`

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import HyperbolicIsometryKM
        sage: A = HyperbolicIsometryKM(identity_matrix(3))
    """
    HFactory = HyperbolicFactoryKM
    HMethods = HyperbolicMethodsUHP

class HyperbolicIsometryHM(HyperbolicIsometry):
    r"""
    Create a hyperbolic isometry in the HM model.

    INPUT:

    - a matrix in `SO(2,1)`

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry import HyperbolicIsometryHM
        sage: A = HyperbolicIsometryHM(identity_matrix(3))
    """
    HFactory = HyperbolicFactoryHM
    HMethods = HyperbolicMethodsUHP

