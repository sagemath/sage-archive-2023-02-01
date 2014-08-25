# -*- coding: utf-8 -*-
r"""
Hyperbolic Models

In this module, a hyperbolic model is a collection of data that allow
the user to implement new models of hyperbolic space with minimal effort.
The data include facts about the underlying set (such as whether the
model is bounded), facts about the metric (such as whether the model is
conformal), facts about the isometry group (such as whether it is a
linear or projective group), and more.  Generally speaking, any data
or method that pertains to the model itself -- rather than the points,
geodesics, or isometries of the model -- is implemented in this module.

Abstractly, a model of hyperbolic space is a connected, simply connected
manifold equipped with a complete Riemannian metric of constant curvature
`-1`.  This module records information sufficient to enable computations
in hyperbolic space without explicitly specifying the underlying set or
its Riemannian metric.  Although, see the
`SageManifolds <http://sagemanifolds.obspm.fr/>`_ project if 
you would like to take this approach.

This module implements the abstract base class for a model of hyperbolic
space of arbitrary dimension.  It also contains the implementations of
specific models of hyperbolic geometry.

AUTHORS:

- Greg Laun (2013): Initial version.

EXAMPLES:

We illustrate how the classes in this module encode data by comparing
the upper half plane (UHP), Poincaré disk (PD) and hyperboloid (HM)
models.  First we import::

    sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP as U
    sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelPD as P
    sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelHM as H

We note that the UHP and PD models are bounded while the HM model is
not::

   sage: U.is_bounded() and P.is_bounded()
   True
   sage: H.is_bounded()
   False

The isometry groups of UHP and PD are projective, while that of HM is
linear:
    
    sage: U.isometry_group_is_projective
    True
    sage: H.isometry_group_is_projective
    False

The models are responsible for determining if the coordinates of points
and the matrix of linear maps are appropriate for constructing points
and isometries in hyperbolic space:
  
    sage: U.point_in_model(2 + I)
    True
    sage: U.point_in_model(2 - I)
    False
    sage: U.point_in_model(2)
    False
    sage: U.boundary_point_in_model(2)
    True

.. TODO::

    Implement a category for metric spaces.
"""

#***********************************************************************
#
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.misc.bindable_class import BindableClass
from sage.misc.lazy_import import lazy_import
from sage.functions.other import imag, real, sqrt
from sage.symbolic.constants import pi
from sage.functions.log import exp
from sage.functions.all import cos, sin, arccosh, arccos, sign
from sage.rings.all import CC, RR, RDF
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.symbolic.pynac import I
from sage.matrix.all import matrix
from sage.categories.homset import Hom
from sage.geometry.hyperbolic_space.hyperbolic_constants import EPSILON
from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPoint, HyperbolicPointUHP
from sage.geometry.hyperbolic_space.hyperbolic_isometry import (
            HyperbolicIsometryUHP, HyperbolicIsometryPD,
            HyperbolicIsometryKM, HyperbolicIsometryHM)
from sage.geometry.hyperbolic_space.hyperbolic_geodesic import (
            HyperbolicGeodesicUHP, HyperbolicGeodesicPD,
            HyperbolicGeodesicKM, HyperbolicGeodesicHM)
from sage.geometry.hyperbolic_space.hyperbolic_coercion import (
            CoercionUHPtoPD, CoercionUHPtoKM, CoercionUHPtoHM,
            CoercionPDtoUHP, CoercionPDtoKM, CoercionPDtoHM,
            CoercionKMtoUHP, CoercionKMtoPD, CoercionKMtoHM,
            CoercionHMtoUHP, CoercionHMtoPD, CoercionHMtoKM)

lazy_import('sage.modules.free_module_element', 'vector')

class HyperbolicModel(Parent, UniqueRepresentation, BindableClass):
    r"""
    Abstract base class for hyperbolic models.
    """
    #name = "abstract hyperbolic space"
    #short_name = "Abstract"
    #bounded = False
    #conformal = False
    #dimension = 0
    #isometry_group = None
    #isometry_group_is_projective = False
    #pt_conversion_dict = {}
    #isom_conversion_dict = {}

    def __init__(self, space, name, short_name, bounded, conformal,
                 dimension, isometry_group, isometry_group_is_projective):
        """
        Initialize ``self``.
        """
        self._name = name
        self._short_name = short_name
        self._bounded = bounded
        self._conformal = conformal
        self._dimension = dimension
        self._isometry_group = isometry_group
        self._isometry_group_is_projective = isometry_group_is_projective
        from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicModels
        Parent.__init__(self, category=HyperbolicModels(space))

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Hyperbolic plane in the {} model".format(self._name)

    def _element_constructor_(self, x, is_boundary=None, **graphics_options):
        """
        Construct an element of ``self``.
        """
        return self.get_point(x, is_boundary, **graphics_options)

    Element = HyperbolicPoint

    def name(self):
        """
        Return the name of this model.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.name()
            'Upper Half Plane Model'
        """
        return self._name

    def short_name(self):
        """
        Return the short name of this model.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.short_name()
            'UHP'
        """
        return self._short_name

    def is_bounded(self):
        """
        Return ``True`` if ``self`` is a bounded model.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.is_bounded()
            True
        """
        return self._bounded

    def is_conformal(self):
        """
        Return ``True`` if ``self`` is a conformal model.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.is_conformal()
            True
        """
        return self._conformal

    def is_isometry_group_projective(self):
        """
        Return ``True`` if the isometry group of ``self`` is projective.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.is_isometry_group_projective()
            True
        """
        return self._isometry_group_is_projective

    def point_in_model(self, p): #Abstract
        r"""
        Return ``True`` if the point is in the given model and ``False``
        otherwise.

        INPUT:

        - any object that can converted into a complex number

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: HyperbolicPlane.UHP.point_in_model(I)
            True
            sage: HyperbolicPlane.UHP.point_in_model(-I)
            False
        """
        return True

    def point_test(self, p): #Abstract
        r"""
        Test whether a point is in the model.  If the point is in the
        model, do nothing.  Otherwise, raise a ``ValueError``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: HyperbolicModelUHP.point_test(2 + I)
            sage: HyperbolicModelUHP.point_test(2 - I)
            Traceback (most recent call last):
            ...
            ValueError: -I + 2 is not a valid point in the UHP model
        """
        if not (cls.point_in_model(p) or cls.boundary_point_in_model(p)):
            error_string = "{0} is not a valid point in the {1} model"
            raise ValueError(error_string.format(p, cls.short_name))

    def boundary_point_in_model(self, p): #Abstract
        r"""
        Return ``True`` if the point is on the ideal boundary of hyperbolic
        space and ``False`` otherwise.

        INPUT:

        - any object that can converted into a complex number

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: HyperbolicPlane.UHP.boundary_point_in_model(I)
            False
        """
        return True

    def bdry_point_test(self, p): #Abstract
        r"""
        Test whether a point is in the model.  If the point is in the
        model, do nothing; otherwise raise a ``ValueError``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: HyperbolicModelUHP.bdry_point_test(2)
            sage: HyperbolicModelUHP.bdry_point_test(1 + I)
            Traceback (most recent call last):
            ...
            ValueError: I + 1 is not a valid boundary point in the UHP model
        """
        if not self._bounded or not cls.boundary_point_in_model(p):
            error_string = "{0} is not a valid boundary point in the {1} model"
            raise ValueError(error_string.format(p, self._short_name))

    def isometry_in_model(self, A): #Abstract
        r"""
        Return ``True`` if the input matrix represents an isometry of the
        given model and ``False`` otherwise.

        INPUT:

        - a matrix that represents an isometry in the appropriate model

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: HyperbolicPlane.UHP.isometry_in_model(identity_matrix(2))
            True

            sage: HyperbolicPlane.UHP.isometry_in_model(identity_matrix(3))
            False
        """
        return True

    def isometry_act_on_point(self, A, p): #Abstract
        r"""
        Given an isometry ``A`` and a point ``p`` in the current model,
        return image of ``p`` unduer the action `A \cdot p`.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: I2 = identity_matrix(2)
            sage: p = HyperbolicPlane.UHP.random_point().coordinates()
            sage: bool(norm(HyperbolicModelUHP.isometry_act_on_point(I2, p) - p) < 10**-9)
            True
        """
        return A * vector(p)

    def isometry_test(self, A): #Abstract
        r"""
        Test whether an isometry is in the model.

        If the isometry is in the model, do nothing. Otherwise, raise
        a ``ValueError``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: HyperbolicModelUHP.isometry_test(identity_matrix(2))
            sage: HyperbolicModelUHP.isometry_test(matrix(2, [I,1,2,1]))
            Traceback (most recent call last):
            ...
            ValueError:
            [I 1]
            [2 1] is not a valid isometry in the UHP model.
        """
        if not cls.isometry_in_model(A):
            error_string = "\n{0} is not a valid isometry in the {1} model."
            raise ValueError(error_string.format(A, cls.short_name))

    def get_point(self, coordinates, is_boundary=None, **graphics_options): #Abstract
        r"""
        Return a point in ``self``.

        Automatically determine the type of point to return given either
        (1) the coordinates of a point in the interior or ideal boundary
        of hyperbolic space or (2) a :class:`HyperbolicPoint` or
        :class:`HyperbolicBdryPoint` object.

        INPUT:

        - a point in hyperbolic space or on the ideal boundary

        OUTPUT:

        - a :class:`HyperbolicPoint`

        EXAMPLES:

        We can create an interior point via the coordinates::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: HyperbolicPlane().UHP().get_point(2*I)
            Point in UHP 2*I

        Or we can create a boundary point via the coordinates::

            sage: HyperbolicPlane().UHP().get_point(23)
            Boundary point in UHP 23

        Or we can create both types of points::

            sage: HyperbolicPlane().UHP().get_point(p)
            Point in UHP 2*I

            sage: HyperbolicPlane().UHP().get_point(q)
            Boundary point in UHP 23

            sage: HyperbolicPlane().UHP().get_point(12 - I)
            Traceback (most recent call last):
            ...
            ValueError: -I + 12 is neither an interior nor boundary point in the UHP model

        ::

            sage: HyperbolicPlane().UHP().get_point(2 + 3*I)
            Point in UHP 3*I + 2

            sage: HyperbolicPlane().PD().get_point(0)
            Point in PD 0

            sage: HyperbolicPlane().KM().get_point((0,0))
            Point in KM (0, 0)

            sage: HyperbolicPlane().HM().get_point((0,0,1))
            Point in HM (0, 0, 1)

            sage: p = HyperbolicPlane().UHP().get_point(I, color="red")
            sage: p.graphics_options()
            {'color': 'red'}

        ::

            sage: HyperbolicPlane().UHP().get_point(12)
            Boundary point in UHP 12

            sage: HyperbolicPlane().UHP().get_point(infinity)
            Boundary point in UHP +Infinity

            sage: HyperbolicPlane().PD().get_point(I)
            Boundary point in PD I

            sage: HyperbolicPlane().KM().get_point((0,-1))
            Boundary point in KM (0, -1)
        """

        if isinstance(coordinates, HyperbolicPoint):
            if coordinates.parent() is not self:
                coordinates = self(coordinates)
            coordinates.update_graphics(True, **graphics_options)
            return coordinates #both Point and BdryPoint

        if is_boundary is None:
            is_boundary = self.boundary_point_in_model(coordinates)
        return self.element_class(self, coordinates, is_boundary, **graphics_options)

    def get_geodesic(self, start, end=None, **graphics_options): #Abstract
        r"""
        Return a geodesic in the appropriate model.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_geodesic(I, 2*I)
            Geodesic in UHP from I to 2*I

            sage: HyperbolicPlane().PD().get_geodesic(0, I/2)
            Geodesic in PD from 0 to 1/2*I

            sage: HyperbolicPlane().KM().get_geodesic((1/2, 1/2), (0,0))
            Geodesic in KM from (1/2, 1/2) to (0, 0)

            sage: HyperbolicPlane().HM().get_geodesic((0,0,1), (1,0, sqrt(2)))
            Geodesic in HM from (0, 0, 1) to (1, 0, sqrt(2))
        """
        if end is None:
            if isinstance(start, HyperbolicGeodesic):
                G = start
                if G.model() is not self:
                    G = G.to_model(self)
                G.update_graphics(True, **graphics_options)
                return G
            raise ValueError("the start and end points must be specified")
        return self._Geodesic(self, self(start), self(end), **graphics_options)

    def get_isometry(self, A):
        r"""
        Return an isometry in the appropriate model given the matrix.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: HyperbolicFactoryUHP.get_isometry(identity_matrix(2))
            Isometry in UHP
            [1 0]
            [0 1]

            sage: HyperbolicFactoryPD.get_isometry(identity_matrix(2))
            Isometry in PD
            [1 0]
            [0 1]

            sage: HyperbolicFactoryKM.get_isometry(identity_matrix(3))
            Isometry in KM
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: HyperbolicFactoryHM.get_isometry(identity_matrix(3))
            Isometry in HM
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        if isinstance(A, HyperbolicIsometry):
            if A.model() is not self:
                return A.to_model(self)
            return A
        return self._Isometry(self, A)

    def random_element(self, **kwargs):
        r"""
        Return a random point in ``self``. 

        The points are uniformly distributed over the rectangle
        `[-10, 10] \times [0, 10 i]` in the upper half plane model.

        EXAMPLES::

            sage: p = HyperbolicPlane().UHP().random_element()
            sage: bool((p.coordinates().imag()) > 0)
            True

            sage: p = HyperbolicPointPD.random_element()
            sage: HyperbolicPlane.PD.point_in_model(p.coordinates())
            True

            sage: p = HyperbolicPointKM.random_element()
            sage: HyperbolicPlane.KM.point_in_model(p.coordinates())
            True

            sage: p = HyperbolicPointHM.random_element().coordinates()
            sage: bool((p[0]**2 + p[1]**2 - p[2]**2 - 1) < 10**-8)
            True
        """
        return self.random_point(**kwds)

    def random_point(self, **kwds):
        """
        Return a random point of ``self``.

        EXAMPLES::

            sage: p = HyperbolicPlane().UHP().random_point()
            sage: bool((p.coordinates().imag()) > 0)
            True

            sage: PD = HyperbolicPlane().PD()
            sage: p = PD.random_point()
            sage: PD.point_in_model(p.coordinates())
            True
        """
        UHP = self.realization_of().UHP()
        return self(UHP.random_point(**kwargs))

    def random_geodesic(self, **kwargs):
        r"""
        Return a random hyperbolic geodesic.

        EXAMPLES::

            sage: h = HyperbolicPlane().UHP().random_geodesic()
            sage: bool((h.endpoints()[0].coordinates()).imag() >= 0)
            True
        """
        UHP = self.realization_of().UHP()
        g_ends = [UHP.random_point(**kwargs) for k in range(2)]
        return self.get_geodesic(self(g_ends[0]), self(g_ends[1]))

    def random_isometry(self, preserve_orientation=True, **kwargs):
        r"""
        Return a random isometry in the model of ``self``.

        INPUT:

        - ``preserve_orientation`` -- if ``True`` return an
          orientation-preserving isometry

        OUTPUT:

        - a hyperbolic isometry

        EXAMPLES::

            sage: A = HyperbolicPlane().UHP().random_isometry()
            sage: A.orientation_preserving()
            True
            sage: B = HyperbolicPlane().UHP().random_isometry(preserve_orientation=False)
            sage: B.orientation_preserving()
            False
        """
        UHP = self.realization_of().UHP()
        A = UHP.random_isometry(preserve_orientation, **kwargs)
        return A.to_model(self)

class HyperbolicModelUHP(HyperbolicModel):
    r"""
    Upper Half Plane model.
    """
    Element = HyperbolicPointUHP
    _Geodesic = HyperbolicGeodesicUHP
    _Isometry = HyperbolicIsometryUHP

    def __init__(self, space):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: TestSuite(UHP).run()
        """
        HyperbolicModel.__init__(self, space,
            name="Upper Half Plane Model", short_name="UHP",
            bounded=True, conformal=True, dimension=2,
            isometry_group="PSL(2, \\RR)", isometry_group_is_projective=True)

    def _coerce_map_from_(self, X):
        """
        Return if the there is a coercion map from ``X`` to ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.has_coerce_map_from(HyperbolicPlane().PD())
            True
            sage: UHP.has_coerce_map_from(HyperbolicPlane().KM())
            True
            sage: UHP.has_coerce_map_from(HyperbolicPlane().HM())
            True
            sage: UHP.has_coerce_map_from(QQ)
            False
        """
        if isinstance(X, HyperbolicModelPD):
            return CoercionPDtoUHP(Hom(X, self))
        if isinstance(X, HyperbolicModelKM):
            return CoercionKMtoUHP(Hom(X, self))
        if isinstance(X, HyperbolicModelHM):
            return CoercionHMtoUHP(Hom(X, self))
        return super(HyperbolicModelUHP, self)._coerce_map_from_(X)

    def point_in_model(self, p): #UHP
        r"""
        Check whether a complex number lies in the open upper half plane.

        EXAMPLES::

            sage: HyperbolicPlane.UHP.point_in_model(1 + I)
            True
            sage: HyperbolicPlane.UHP.point_in_model(infinity)
            False
            sage: HyperbolicPlane.UHP.point_in_model(CC(infinity))
            False
            sage: HyperbolicPlane.UHP.point_in_model(RR(infinity))
            False
            sage: HyperbolicPlane.UHP.point_in_model(1)
            False
            sage: HyperbolicPlane.UHP.point_in_model(12)
            False
            sage: HyperbolicPlane.UHP.point_in_model(1 - I)
            False
            sage: HyperbolicPlane.UHP.point_in_model(-2*I)
            False
        """
        return bool(imag(CC(p)) > 0)

    def boundary_point_in_model(self, p): #UHP
        r"""
        Check whether a complex number is a real number or ``\infty``.
        In the ``UHP.model_name_name``, this is the ideal boundary of
        hyperbolic space.

        EXAMPLES::

            sage: HyperbolicPlane.UHP.boundary_point_in_model(1 + I)
            False
            sage: HyperbolicPlane.UHP.boundary_point_in_model(infinity)
            True
            sage: HyperbolicPlane.UHP.boundary_point_in_model(CC(infinity))
            True
            sage: HyperbolicPlane.UHP.boundary_point_in_model(RR(infinity))
            True
            sage: HyperbolicPlane.UHP.boundary_point_in_model(1)
            True
            sage: HyperbolicPlane.UHP.boundary_point_in_model(12)
            True
            sage: HyperbolicPlane.UHP.boundary_point_in_model(1 - I)
            False
            sage: HyperbolicPlane.UHP.boundary_point_in_model(-2*I)
            False
        """
        im = abs(imag(CC(p)).n())
        return bool( (im < EPSILON) or (p == infinity) )

    def isometry_act_on_point(self, A, p): #UHP
        r"""
        Given an isometry ``A`` and a point ``p`` in the current model,
        return image of ``p`` unduer the action `A \cdot p`.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: I2 = identity_matrix(2)
            sage: p = HyperbolicPlane.UHP.random_point().coordinates()
            sage: bool(norm(HyperbolicModelUHP.isometry_act_on_point(I2, p) - p) < 10**-9)
            True
        """
        return mobius_transform(A, p)

    def isometry_in_model(self, A): #UHP
        r"""
        Check that ``A`` acts as an isometry on the upper half plane.
        That is, ``A`` must be an invertible `2 \times 2` matrix with real
        entries.

        EXAMPLES::

            sage: A = matrix(2,[1,2,3,4])
            sage: HyperbolicPlane.UHP.isometry_in_model(A)
            True
            sage: B = matrix(2,[I,2,4,1])
            sage: HyperbolicPlane.UHP.isometry_in_model(B)
            False
        """
        return bool(A.ncols() == 2 and A.nrows() == 2 and
            sum([k in RR for k in A.list()]) == 4 and
            abs(A.det()) > -EPSILON)

    def point_to_model(self, coordinates, model_name): #UHP
        r"""
        Convert ``coordinates`` from the current model to the model
        specified in ``model_name``.

        INPUT:

        - ``coordinates`` -- the coordinates of a valid point in the
          current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES::

            sage: HyperbolicPlane.UHP.point_to_model(I, 'UHP')
            I
            sage: HyperbolicPlane.UHP.point_to_model(I, 'PD')
            0
            sage: HyperbolicPlane.UHP.point_to_model(3 + I, 'KM')
            (6/11, 9/11)
            sage: HyperbolicPlane.UHP.point_to_model(3 + I, 'HM')
            (3, 9/2, 11/2)

        It is an error to try to convert a boundary point to a model
        that doesn't support boundary points::

            sage: HyperbolicPlane.UHP.point_to_model(infinity, 'HM')
            Traceback (most recent call last):
            ...
            NotImplementedError: boundary points are not implemented for the HM model
        """
        p = coordinates
        if (cls.boundary_point_in_model(p) and not
                ModelFactory.find_model(model_name).is_bounded()):
            raise NotImplementedError("boundary points are not implemented for"
                                      " the {0} model".format(model_name))
        return cls.pt_conversion_dict[model_name](coordinates)

    def isometry_to_model(self, A, model_name): # UHP
        r"""
        Convert ``A`` from the current model to the model specified in
        ``model_name``.

        INPUT:

        - ``A`` -- a matrix in the current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelUHP
            sage: HyperbolicModelUHP.isometry_to_model(matrix(2,[0, 1, 1, 0]),'PD')
            [0 I]
            [I 0]
        """
        cls.isometry_test(A)
        if A.det() < 0 and model_name == 'PD':
            return cls.isom_conversion_dict[model_name](I * A)
        return cls.isom_conversion_dict[model_name](A)

    def get_background_graphic(self, **bdry_options): #UHP
        r"""
        Return a graphic object that makes the model easier to visualize.
        For the upper half space, the background object is the ideal boundary.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: circ = HyperbolicFactoryUHP.get_background_graphic()
        """
        from sage.plot.line import line
        bd_min = bdry_options.get('bd_min', -5)
        bd_max = bdry_options.get('bd_max', 5)
        return line(((bd_min, 0), (bd_max, 0)), color='black')

    #################
    # Point Methods #
    #################

    def point_dist(self, p1, p2):
        r"""
        Compute the distance between two points in the Upper Half Plane
        using the hyperbolic metric.

        EXAMPLES::

           sage: HyperbolicMethodsUHP.point_dist(4.0*I,I)
           1.38629436111989
        """
        cls.model().point_test(p1)
        cls.model().point_test(p2)
        num = (real(p2) - real(p1))**2 + (imag(p2) - imag(p1))**2
        denom = 2*imag(p1)*imag(p2)
        if denom == 0:
            return infinity
        return arccosh(1 + num/denom)

    def symmetry_in(self, p):
        r"""
        Return the involutary isometry fixing the given point.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.symmetry_in(3 + 2*I)
            [  3/2 -13/2]
            [  1/2  -3/2]
        """
        self.point_test(p)
        x, y = real(p), imag(p)
        if y > 0:
            return matrix(2,[x/y,-(x**2/y) - y,1/y,-(x/y)])

    def random_point(self, **kwargs):
        r"""
        Return a random point in the upper half plane. The points are
        uniformly distributed over the rectangle `[-10, 10] \times [0, 10i]`.

        EXAMPLES::

            sage: p = HyperbolicPlane().UHP().random_point()
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

    def boundary_points(self, p1, p2):
        r"""
        Given two points ``p1`` and ``p2`` in the hyperbolic plane,
        determine the endpoints of the complete hyperbolic geodesic
        through ``p1`` and ``p2``.

        INPUT:

        - ``p1``, ``p2`` -- points in the hyperbolic plane

        OUTPUT:

        - a list of boundary points

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.boundary_points(I, 2*I)
            [0, +Infinity]
            sage: HyperbolicMethodsUHP.boundary_points(1 + I, 2 + 4*I)
            [-sqrt(65) + 9, sqrt(65) + 9]
        """
        if p1 == p2:
            raise ValueError("{} and {} are not distinct".format(p1, p2))
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
        x,y = [real(k) for k in self.boundary_points(start, end)]
        if x == infinity:
            M = matrix(2,[1,-2*y,0,-1])
        elif y == infinity:
            M = matrix(2,[1,-2*x,0,-1])
        else:
            M = matrix(2,[(x+y)/(y-x),-2*x*y/(y-x),2/(y-x),-(x+y)/(y-x)])
        return M

    def common_perpendicular(self, start_1, end_1, start_2, end_2):
        r"""
        Return the unique hyperbolic geodesic perpendicular to two given
        geodesics, if such a geodesic exists; otherwise raise a
        ``ValueError``.

        INPUT:

        - ``other`` -- a hyperbolic geodesic in current model

        OUTPUT:

        - a hyperbolic geodesic

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.common_perpendicular(2, 3, 4, 5)
            [1/2*sqrt(3) + 7/2, -1/2*sqrt(3) + 7/2]

        It is an error to ask for the common perpendicular of two
        intersecting geodesics::

            sage: HyperbolicMethodsUHP.common_perpendicular(2, 4, 3, infinity)
            Traceback (most recent call last):
            ...
            ValueError: geodesics intersect; no common perpendicular exists
        """
        A = self.reflection_in(start_1, end_1)
        B = self.reflection_in(start_2, end_2)
        C = A*B
        if self.classification(C) != 'hyperbolic':
            raise ValueError("geodesics intersect; " +
                             "no common perpendicular exists")
        return self.fixed_point_set(C)

    def intersection(self, start_1, end_1, start_2, end_2):
        r"""
        Return the point of intersection of two complete geodesics
        (if such a point exists).

        INPUT:

        - ``other`` -- a hyperbolic geodesic in the current model

        OUTPUT:

        - a hyperbolic point

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.intersection(3, 5, 4, 7)
            [2/3*sqrt(-2) + 13/3]

        If the given geodesics do not intersect, the function raises an
        error::

            sage: HyperbolicMethodsUHP.intersection(4, 5, 5, 7)
            Traceback (most recent call last):
            ...
            ValueError: geodesics don't intersect

        If the given geodesics are identical, return that
        geodesic::

            sage: HyperbolicMethodsUHP.intersection(4 + I, 18*I, 4 + I, 18*I)
            [-1/8*sqrt(114985) - 307/8, 1/8*sqrt(114985) - 307/8]
        """
        start_1, end_1 = sorted(self.boundary_points(start_1, end_1))
        start_2, end_2 = sorted(self.boundary_points(start_2, end_2))
        if start_1 == start_2 and end_1 == end_2: # Unoriented geods are same
            return [start_1, end_1]
        if start_1 == start_2:
            return start_1
        elif end_1 == end_2:
            return end_1
        A = self.reflection_in(start_1, end_1)
        B = self.reflection_in(start_2, end_2)
        C = A*B
        if self.classification(C) in ['hyperbolic', 'parabolic']:
            raise ValueError("geodesics don't intersect")
        return self.fixed_point_set(C)

    def perpendicular_bisector(self, start, end):
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
        start, end = sorted((start, end))
        d = self.point_dist(start, end)/2
        end_1, end_2 = self.boundary_points(start, end)
        T = (matrix(2,[exp(d/2),0,0,exp(-d/2)])*
             matrix(2,[cos(pi/4),-sin(pi/4),sin(pi/4),cos(pi/4)]))
        S= self._to_std_geod(end_1, start, end_2)
        H = S.inverse()*T*S
        return [mobius_transform(H ,k) for k in [end_1, end_2]]

    def midpoint(self, start, end):
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
        half_dist = self.point_dist(start, end)/2
        end_1,end_2 = self.boundary_points(start, end)
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

    def geod_dist_from_point(self, start, end, p):
        r"""
        Return the hyperbolic distance from a given hyperbolic geodesic
        and a hyperbolic point.

        INPUT:

        - ``other`` -- a hyperbolic point in the same model

        OUTPUT:

        - the hyperbolic distance

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: HyperbolicMethodsUHP.geod_dist_from_point(2, 2 + I, I)
            arccosh(1/10*sqrt(5)*((sqrt(5) - 1)^2 + 4) + 1)

        If `p` is a boundary point, the distance is infinity::

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
        (bd_1, bd_2) = self.boundary_points(start, end)
        if bd_1 + bd_2 != infinity:
            # Not a straight line
            # Map the endpoints to 0 and infinity and the midpoint to 1.
            T = self._crossratio_matrix(bd_1, (bd_1 + bd_2)/2, bd_2)
        else:
            # Is a straight line
            # Map the endpoints to 0 and infinity and another endpoint
            # to 1
            T = self._crossratio_matrix(bd_1, bd_1 + 1, bd_2)
        x = mobius_transform(T, p)
        return self.point_dist(x, abs(x)*I)

    def uncomplete(self, start, end):
        r"""
        Return starting and ending points of a geodesic whose completion
        is the geodesic starting at ``start`` and ending at ``end``.

        INPUT:

        - ``start`` -- a real number or infinity
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
        A = self._to_std_geod(start, p, end).inverse()
        p1 = mobius_transform(A, I/Integer(3))
        p2 = mobius_transform(A, 3*I)
        return [p1, p2]

    def angle(self, start_1, end_1, start_2, end_2):
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
        (p_1,p_2) = sorted(self.boundary_points(start_1, end_1))
        (q_1,q_2) = sorted(self.boundary_points(start_2, end_2))
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

    ####################
    # Isometry Methods #
    ####################

    # Move to HyperbolicIsometryUHP
    def orientation_preserving(self, M):
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

    # Move to HyperbolicIsometryUHP
    def classification(self, M):
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
        A = A / (abs(A.det()).sqrt())
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
                raise ValueError("something went wrong with classification:" +
                                 " trace is {}".format(A.trace()))
        else:  #The isometry reverses orientation.
            if tau < EPSILON:
                return 'reflection'
            else:
                return 'orientation-reversing hyperbolic'

    # Move to HyperbolicIsometryUHP
    def translation_length(self, M):
        r"""
        For hyperbolic elements, return the translation length;
        otherwise, raise a ``ValueError``.

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
        tau = sqrt((M / sqrt(d)).trace()**2)
        if self.classification(M) in ['hyperbolic', 'oriention-reversing hyperbolic']:
            return 2*arccosh(tau/2)
        raise TypeError("translation length is only defined for hyperbolic"
                        " transformations")

    # Move to HyperbolicIsometryUHP
    def fixed_point_set(self, M):
        r"""
        Return the a list containing the fixed point set of
        orientation-preserving isometries.

        OUTPUT:

        - a list of hyperbolic points

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
        M = M / sqrt(d)
        tau = M.trace()**2
        M_cls = self.classification(M)
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

    # Move to HyperbolicIsometryUHP
    def repelling_fixed_point(self, M):
        r"""
        For a hyperbolic isometry, return the repelling fixed point;
        otherwise raise a ``ValueError``.

        OUTPUT:

        - a hyperbolic point

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: A = matrix(2,[4,0,0,1/4])
            sage: HyperbolicMethodsUHP.repelling_fixed_point(A)
            0
        """
        if self.classification(M) not in ['hyperbolic',
                                         'orientation-reversing hyperbolic']:
            raise ValueError("repelling fixed point is defined only" +
                             "for hyperbolic isometries")
        v = M.eigenmatrix_right()[1].column(1)
        if v[1] == 0:
            return infinity
        return v[0]/v[1]

    # Move to HyperbolicIsometryUHP
    def attracting_fixed_point(self, M):
        r"""
        For a hyperbolic isometry, return the attracting fixed point;
        otherwise raise a ``ValueError``.

        OUTPUT:

        - a hyperbolic point

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import HyperbolicMethodsUHP
            sage: A = matrix(2,[4,0,0,1/4])
            sage: HyperbolicMethodsUHP.attracting_fixed_point(A)
            +Infinity
        """
        if self.classification(M) not in \
                ['hyperbolic', 'orientation-reversing hyperbolic']:
            raise ValueError("Attracting fixed point is defined only" +
                             "for hyperbolic isometries.")
        v = M.eigenmatrix_right()[1].column(0)
        if v[1] == 0:
            return infinity
        return v[0]/v[1]

    def isometry_from_fixed_points(self, repel, attract):
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
            A = self._mobius_sending([infinity, attract, attract + 1],
                                     [infinity, attract, attract + 2])
        elif attract == infinity:
            A = self._mobius_sending([repel, infinity, repel + 1],
                                     [repel, infinity, repel + 2])
        else:
            A = self._mobius_sending([repel, attract, infinity],
                                     [repel, attract, max(repel, attract) + 1])
        return A

    def random_isometry(self, preserve_orientation=True, **kwargs):
        r"""
        Return a random isometry in the Upper Half Plane model.

        INPUT:

        - ``preserve_orientation`` -- if ``True`` return an
          orientation-preserving isometry

        OUTPUT:

        - a hyperbolic isometry

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
        M = M / (M.det()).abs().sqrt()
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

    def _to_std_geod(self, start, p, end):
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
        return B * self._crossratio_matrix(start, p, end)

    def _crossratio_matrix(self, p_0, p_1, p_2):
        r"""
        Given three points (the list `p`) in `\mathbb{CP}^{1}` in affine
        coordinates, return the linear fractional transformation taking
        the elements of `p` to `0`,`1`, and `\infty'.

        INPUT:

        - a list of three distinct elements of three distinct elements
          of `\mathbb{CP}^1` in affine coordinates; that is, each element
          must be a complex number, `\infty`, or symbolic.

        OUTPUT:

        - an element of `\GL(2,\CC)`

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

    def _mobius_sending(self, list1, list2):
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
        if len(list1) != 3 or len(list2) != 3:
            raise TypeError("mobius_sending requires each list to be three points long")
        pl = list1 + list2
        z = pl[0:3]
        w = pl[3:6]
        A = self._crossratio_matrix(z[0],z[1],z[2])
        B = self._crossratio_matrix(w[0],w[1],w[2])
        return B.inverse() * A

class HyperbolicModelPD(HyperbolicModel):
    r"""
    Poincaré Disk Model.
    """
    _Geodesic = HyperbolicGeodesicPD
    _Isometry = HyperbolicIsometryPD

    def __init__(self, space):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: TestSuite(PD).run()
        """
        HyperbolicModel.__init__(self, space,
            name="Poincare Disk Model", # u"Poincaré Disk Model"
            short_name="PD",
            bounded=True, conformal=True, dimension=2,
            isometry_group="PU(1, 1)", isometry_group_is_projective=True)

    def _coerce_map_from_(self, X):
        """
        Return if the there is a coercion map from ``X`` to ``self``.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: PD.has_coerce_map_from(HyperbolicPlane().UHP())
            True
            sage: PD.has_coerce_map_from(HyperbolicPlane().KM())
            True
            sage: PD.has_coerce_map_from(HyperbolicPlane().HM())
            True
            sage: PD.has_coerce_map_from(QQ)
            False
        """
        if isinstance(X, HyperbolicModelUHP):
            return CoercionUHPtoPD(Hom(X, self))
        if isinstance(X, HyperbolicModelKM):
            return CoercionKMtoPD(Hom(X, self))
        if isinstance(X, HyperbolicModelHM):
            return CoercionHMtoPD(Hom(X, self))
        return super(HyperbolicModelPD, self)._coerce_map_from_(X)

    def point_in_model(self, p): #PD
        r"""
        Check whether a complex number lies in the open unit disk.

        EXAMPLES::

            sage: HyperbolicPlane.PD.point_in_model(1.00)
            False

            sage: HyperbolicPlane.PD.point_in_model(1/2 + I/2)
            True

            sage: HyperbolicPlane.PD.point_in_model(1 + .2*I)
            False
        """
        return bool(abs(CC(p)) < 1)

    def boundary_point_in_model(self, p): #PD
        r"""
        Check whether a complex number lies in the open unit disk.

        EXAMPLES::

            sage: HyperbolicPlane.PD.boundary_point_in_model(1.00)
            True

            sage: HyperbolicPlane.PD.boundary_point_in_model(1/2 + I/2)
            False

            sage: HyperbolicPlane.PD.boundary_point_in_model(1 + .2*I)
            False
        """
        return bool(abs(abs(CC(p))- 1) < EPSILON)

    def isometry_act_on_point(self, A, p): #PD
        r"""
        Given an isometry ``A`` and a point ``p`` in the current model,
        return image of ``p`` unduer the action `A \cdot p`.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelPD
            sage: I2 = identity_matrix(2)
            sage: q = HyperbolicPlane.PD.random_point().coordinates()
            sage: bool(norm(HyperbolicModelPD.isometry_act_on_point(I2, q) - q) < 10**-9)
            True
        """
        _image = mobius_transform(A, p)
        if not PD_preserve_orientation(A):
            return mobius_transform(I*matrix(2,[0,1,1,0]), _image)
        return _image

    def isometry_in_model(self, A): #PD
        r"""
        Check if the given matrix ``A`` is in the group `U(1,1)`.

        EXAMPLES::

            sage: z = [CC.random_element() for k in range(2)]; z.sort(key=abs)
            sage: A = matrix(2,[z[1], z[0],z[0].conjugate(),z[1].conjugate()])
            sage: HyperbolicPlane.PD.isometry_in_model(A)
            True
        """
        # alpha = A[0][0]
        # beta = A[0][1]
        # Orientation preserving and reversing
        return PD_preserve_orientation(A) or PD_preserve_orientation(I*A)

    def point_to_model(self, coordinates, model_name): #PD
        r"""
        Convert ``coordinates`` from the current model to the model
        specified in ``model_name``.

        INPUT:

        - ``coordinates`` -- the coordinates of a valid point in the
          current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES::

            sage: HyperbolicPlane.PD.point_to_model(0, 'UHP')
            I

            sage: HyperbolicPlane.PD.point_to_model(I, 'UHP')
            +Infinity

            sage: HyperbolicPlane.PD.point_to_model(-I, 'UHP')
            0
        """
        if model_name == 'UHP' and coordinates == I:
            return infinity
        return super(HyperbolicModelPD, cls).point_to_model(coordinates, model_name)

    def isometry_to_model(self, A, model_name): #PD
        r"""
        Convert ``A`` from the current model to the model specified in
        ``model_name``.

        INPUT:

        - ``A`` -- a matrix in the current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES:

        We check that orientation-reversing isometries behave as they
        should::

            sage: HyperbolicPlane.PD.isometry_to_model(matrix(2,[0,I,I,0]),'UHP')
            [ 0 -1]
            [-1  0]
        """
        cls.isometry_test(A)
        # Check for orientation-reversing isometries.
        if (not PD_preserve_orientation(A) and model_name == 'UHP'):
            return cls.isom_conversion_dict[model_name](I*A)
        return cls.isom_conversion_dict[model_name](A)

    def get_background_graphic(self, **bdry_options): #PD
        r"""
        Return a graphic object that makes the model easier to visualize.
        For the Poincare disk, the background object is the ideal boundary.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: circ = HyperbolicFactoryPD.get_background_graphic()
        """
        from sage.plot.circle import circle
        return circle((0,0), 1, axes=False, color='black')

#####################################################################
## Klein disk model

class HyperbolicModelKM(HyperbolicModel):
    r"""
    Klein Model.
    """
    _Geodesic = HyperbolicGeodesicKM
    _Isometry = HyperbolicIsometryKM

    def __init__(self, space):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: KM = HyperbolicPlane().KM()
            sage: TestSuite(KM).run()
        """
        HyperbolicModel.__init__(self, space,
            name="Klein Disk Model", short_name="KM",
            bounded=True, conformal=False, dimension=2,
            isometry_group="PSO(2, 1)", isometry_group_is_projective=True)

    def _coerce_map_from_(self, X):
        """
        Return if the there is a coercion map from ``X`` to ``self``.

        EXAMPLES::

            sage: KM = HyperbolicPlane().UHP()
            sage: KM.has_coerce_map_from(HyperbolicPlane().UHP())
            True
            sage: KM.has_coerce_map_from(HyperbolicPlane().PD())
            True
            sage: KM.has_coerce_map_from(HyperbolicPlane().HM())
            True
            sage: KM.has_coerce_map_from(QQ)
            False
        """
        if isinstance(X, HyperbolicModelUHP):
            return CoercionUHPtoKM(Hom(X, self))
        if isinstance(X, HyperbolicModelPD):
            return CoercionPDtoKM(Hom(X, self))
        if isinstance(X, HyperbolicModelHM):
            return CoercionHMtoKM(Hom(X, self))
        return super(HyperbolicModelKM, self)._coerce_map_from_(X)

    def point_in_model(self, p): #KM
        r"""
        Check whether a point lies in the open unit disk.

        EXAMPLES::

            sage: HyperbolicPlane.KM.point_in_model((1,0))
            False

            sage: HyperbolicPlane.KM.point_in_model((1/2 , 1/2))
            True

            sage: HyperbolicPlane.KM.point_in_model((1 , .2))
            False
        """
        return len(p) == 2 and bool(p[0]**2 + p[1]**2 < 1)

    def boundary_point_in_model(self, p): #KM
        r"""
        Check whether a point lies in the unit circle, which corresponds
        to the ideal boundary of the hyperbolic plane in the Klein model.

        EXAMPLES::

            sage: HyperbolicPlane.KM.boundary_point_in_model((1,0))
            True

            sage: HyperbolicPlane.KM.boundary_point_in_model((1/2 , 1/2))
            False

            sage: HyperbolicPlane.KM.boundary_point_in_model((1 , .2))
            False
        """
        return len(p) == 2 and bool(abs(p[0]**2 + p[1]**2 - 1) < EPSILON)

    def isometry_act_on_point(self, A, p): #KM
        r"""
        Given an isometry ``A`` and a point ``p`` in the current model,
        return image of ``p`` unduer the action `A \cdot p`.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelKM
            sage: I3 = identity_matrix(3)
            sage: v = vector(HyperbolicPlane.KM.random_point().coordinates())
            sage: bool(norm(HyperbolicModelKM.isometry_act_on_point(I3, v) - v) < 10**-9)
            True
        """
        v = A*vector((list(p) + [1]))
        if v[2] == 0:
            return infinity
        return v[0:2]/v[2]

    def isometry_in_model(self, A): #KM
        r"""
        Check if the given matrix ``A`` is in the group `SO(2,1)`.

        EXAMPLES::

            sage: A = matrix(3,[1, 0, 0, 0, 17/8, 15/8, 0, 15/8, 17/8])
            sage: HyperbolicPlane.KM.isometry_in_model(A)
            True
        """
        from sage.geometry.hyperbolic_space.hyperbolic_constants import LORENTZ_GRAM
        return bool((A*LORENTZ_GRAM*A.transpose() - LORENTZ_GRAM).norm()**2 <
                EPSILON)

    def point_to_model(self, coordinates, model_name): #KM
        r"""
        Convert ``coordinates`` from the current model to the model
        specified in ``model_name``.

        INPUT:

        - ``coordinates`` -- the coordinates of a valid point in the
          current model
        - ``model_name`` -- a string denoting the model to be converted to

        OUTPUT:

        - the coordinates of a point in the ``short_name`` model

        EXAMPLES::

            sage: HyperbolicPlane.KM.point_to_model((0, 0), 'UHP')
            I

            sage: HyperbolicPlane.KM.point_to_model((0, 0), 'HM')
            (0, 0, 1)

            sage: HyperbolicPlane.KM.point_to_model((0,1), 'UHP')
            +Infinity
        """
        if model_name == 'UHP' and tuple(coordinates) == (0,1):
            return infinity
        return super(HyperbolicModelKM, cls).point_to_model(coordinates, model_name)

    def get_background_graphic(self, **bdry_options): #KM
        r"""
        Return a graphic object that makes the model easier to visualize.
        For the Klein model, the background object is the ideal boundary.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: circ = HyperbolicFactoryKM.get_background_graphic()
        """
        from sage.plot.circle import circle
        return circle((0,0), 1, axes=False, color='black')

#####################################################################
## Hyperboloid model

class HyperbolicModelHM(HyperbolicModel):
    r"""
    Hyperboloid Model.
    """
    _Geodesic = HyperbolicGeodesicHM
    _Isometry = HyperbolicIsometryHM

    def __init__(self, space):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: HM = HyperbolicPlane().HM()
            sage: TestSuite(HM).run()
        """
        HyperbolicModel.__init__(self, space,
            name="Hyperboloid Model", short_name="HM",
            bounded=False, conformal=True, dimension=2,
            isometry_group="SO(2, 1)", isometry_group_is_projective=True)

    def _coerce_map_from_(self, X):
        """
        Return if the there is a coercion map from ``X`` to ``self``.

        EXAMPLES::

            sage: HM = HyperbolicPlane().UHP()
            sage: HM.has_coerce_map_from(HyperbolicPlane().UHP())
            True
            sage: HM.has_coerce_map_from(HyperbolicPlane().PD())
            True
            sage: HM.has_coerce_map_from(HyperbolicPlane().KM())
            True
            sage: HM.has_coerce_map_from(QQ)
            False
        """
        if isinstance(X, HyperbolicModelUHP):
            return CoercionUHPtoHM(Hom(X, self))
        if isinstance(X, HyperbolicModelPD):
            return CoercionPDtoHM(Hom(X, self))
        if isinstance(X, HyperbolicModelKM):
            return CoercionKMtoHM(Hom(X, self))
        return super(HyperbolicModelHM, self)._coerce_map_from_(X)

    def point_in_model(self, p): #HM
        r"""
        Check whether a complex number lies in the hyperboloid.

        EXAMPLES::

            sage: HyperbolicPlane.HM.point_in_model((0,0,1))
            True

            sage: HyperbolicPlane.HM.point_in_model((1,0,sqrt(2)))
            True

            sage: HyperbolicPlane.HM.point_in_model((1,2,1))
            False
        """
        return len(p) == 3 and bool(p[0]**2 + p[1]**2 - p[2]**2 + 1 < EPSILON)

    def boundary_point_in_model(self, p):  #HM
        r"""
        Return ``False`` since the Hyperboloid model has no boundary points.

        EXAMPLES::

            sage: HyperbolicPlane.HM.boundary_point_in_model((0,0,1))
            False

            sage: HyperbolicPlane.HM.boundary_point_in_model((1,0,sqrt(2)))
            False

            sage: HyperbolicPlane.HM.boundary_point_in_model((1,2,1))
            False
        """
        return False

    def isometry_in_model(self, A):  #HM
        r"""
        Test that the matrix ``A`` is in the group `SO(2,1)^+`.

        EXAMPLES::

           sage: A = diagonal_matrix([1,1,-1])
           sage: HyperbolicPlane.HM.isometry_in_model(A)
           True
        """
        from sage.geometry.hyperbolic_space.hyperbolic_constants import LORENTZ_GRAM
        return bool((A*LORENTZ_GRAM*A.transpose() - LORENTZ_GRAM).norm()**2 < EPSILON)

    def get_background_graphic(self, **bdry_options): #HM
        r"""
        Return a graphic object that makes the model easier to visualize.
        For the hyperboloid model, the background object is the hyperboloid
        itself.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: circ = HyperbolicFactoryPD.get_background_graphic()
        """
        hyperboloid_opacity = bdry_options.get('hyperboloid_opacity', .1)
        z_height = bdry_options.get('z_height', 7.0)
        x_max = sqrt((z_height**2 - 1) / 2.0)
        from sage.plot.plot3d.all import plot3d
        from sage.all import var
        (x,y) = var('x,y')
        return plot3d((1 + x**2 + y**2).sqrt(), (x, -x_max, x_max),
                      (y,-x_max, x_max), opacity = hyperboloid_opacity, **bdry_options)

#####################################################################
## Helper functions

def PD_preserve_orientation(A):
    r"""
    For a PD isometry, determine if it preserves orientation.
    This test is more more involved than just checking the sign
    of the determinant, and it is used a few times in this file.

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_model import PD_preserve_orientation as orient
        sage: orient(matrix(2, [-I, 0, 0, I]))
        True
        sage: orient(matrix(2, [0, I, I, 0]))
        False
    """
    return bool(A[1][0] == A[0][1].conjugate() and A[1][1] == A[0][0].conjugate()
                and abs(A[0][0]) - abs(A[0][1]) != 0)

def mobius_transform(A, z):
    r"""
    Given a matrix ``A`` in `GL(2, \CC)` and a point ``z`` in the complex
    plane return the mobius transformation action of ``A`` on ``z``.

    INPUT:

    - ``A`` -- a `2 \times 2` invertible matrix over the complex numbers
    - ``z`` -- a complex number or infinity

    OUTPUT:

    - a complex number or infinity

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_methods import mobius_transform
        sage: mobius_transform(matrix(2,[1,2,3,4]),2 + I)
        2/109*I + 43/109
        sage: y = var('y')
        sage: mobius_transform(matrix(2,[1,0,0,1]),x + I*y)
        x + I*y

    The matrix must be square and `2 \times 2`::

        sage: mobius_transform(matrix([[3,1,2],[1,2,5]]),I)
        Traceback (most recent call last):
        ...
        TypeError: A must be an invertible 2x2 matrix over the complex numbers or a symbolic ring

        sage: mobius_transform(identity_matrix(3),I)
        Traceback (most recent call last):
        ...
        TypeError: A must be an invertible 2x2 matrix over the complex numbers or a symbolic ring

    The matrix can be symbolic or can be a matrix over the real
    or complex numbers, but must be invertible::

        sage: (a,b,c,d) = var('a,b,c,d');
        sage: mobius_transform(matrix(2,[a,b,c,d]),I)
        (I*a + b)/(I*c + d)

        sage: mobius_transform(matrix(2,[0,0,0,0]),I)
        Traceback (most recent call last):
        ...
        TypeError: A must be an invertible 2x2 matrix over the complex numbers or a symbolic ring
    """
    if A.ncols() == 2 and A.nrows() == 2 and A.det() != 0:
            (a,b,c,d) = A.list()
            if z == infinity:
                if c == 0:
                    return infinity
                return a/c
            if a*d - b*c < 0:
                w = z.conjugate() # Reverses orientation
            else:
                w = z
            if c*z + d == 0:
                return infinity
            else:
                return (a*w + b)/(c*w + d)
    else:
        raise TypeError("A must be an invertible 2x2 matrix over the"
                        " complex numbers or a symbolic ring")

