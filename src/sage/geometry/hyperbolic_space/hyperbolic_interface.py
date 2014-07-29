# -*- coding: utf-8 -*-
r"""
Interface to Hyperbolic Models

This module provides a convenient interface for interacting with models
of hyperbolic space as well as their points, geodesics, and isometries.

The primary point of this module is to allow the code that implements
hyperbolic space to be sufficiently decoupled while still providing a
convenient user experience.

The interfaces are by default given abbreviated names.  For example,
UHP (upper half plane model), PD (Poincaré disk model), KM (Klein disk
model), and HM (hyperboloid model).

AUTHORS:

- Greg Laun (2013): Initial version.

EXAMPLES::

    sage: UHP.point(2 + I)
    Point in UHP I + 2

    sage: PD.point(1/2 + I/2)
    Point in PD 1/2*I + 1/2
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
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_model',
            ['HyperbolicModel','HyperbolicModelUHP',
             'HyperbolicModelPD', 'HyperbolicModelHM', 'HyperbolicModelKM'])

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_factory',
            ['HyperbolicFactory', 'HyperbolicFactoryUHP',
             'HyperbolicFactoryPD', 'HyperbolicFactoryHM',
             'HyperbolicFactoryKM'])

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_point',
            ['HyperbolicPoint', 'HyperbolicPointUHP', 'HyperbolicPointPD',
             'HyperbolicPointHM', 'HyperbolicPointKM'])

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_bdry_point',
            ['HyperbolicBdryPointUHP', 'HyperbolicBdryPointPD',
             'HyperbolicBdryPointHM', 'HyperbolicBdryPointKM'])

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_geodesic',
            ['HyperbolicGeodesic', 'HyperbolicGeodesicUHP',
             'HyperbolicGeodesicPD', 'HyperbolicGeodesicHM',
             'HyperbolicGeodesicKM'])


lazy_import('sage.geometry.hyperbolic_space.hyperbolic_isometry',
            ['HyperbolicIsometry', 'HyperbolicIsometryUHP',
             'HyperbolicIsometryPD', 'HyperbolicIsometryHM',
             'HyperbolicIsometryKM'])


class HyperbolicUserInterface(UniqueRepresentation):
    r"""
    Abstract base class for hyperbolic interfaces.  These provide a user
    interface for interacting with models of hyperbolic geometry without
    having the interface dictate the class structure.
    """
    HModel = HyperbolicModel
    HFactory = HyperbolicFactory
    HPoint = HyperbolicPoint
    HIsometry = HyperbolicIsometry
    HGeodesic = HyperbolicGeodesic

    @classmethod
    def model_name(cls):
        r"""
        Return the full name of the hyperbolic model.

        EXAMPLES::

            sage: UHP.model_name()
            'Upper Half Plane Model'
            sage: PD.model_name()
            'Poincare Disk Model'
            sage: KM.model_name()
            'Klein Disk Model'
            sage: HM.model_name()
            'Hyperboloid Model'
        """
        return cls.HModel.name

    @classmethod
    def short_name(cls):
        r"""
        Return the short name of the hyperbolic model.

        EXAMPLES::

            sage: UHP.short_name()
            'UHP'
            sage: PD.short_name()
            'PD'
            sage: HM.short_name()
            'HM'
            sage: KM.short_name()
            'KM'
        """
        return cls.HModel.short_name

    @classmethod
    def is_bounded(cls):
        r"""
        Return ``True`` if the model is bounded and ``False`` otherwise.

        EXAMPLES::

            sage: UHP.is_bounded()
            True
            sage: PD.is_bounded()
            True
            sage: KM.is_bounded()
            True
            sage: HM.is_bounded()
            False
        """
        return cls.HModel.bounded

    @classmethod
    def point(cls, p, **kwargs):
        r"""
        Return a :class:`HyperbolicPoint` object in the current
        model with coordinates ``p``.

        EXAMPLES::

            sage: UHP.point(0)
            Boundary point in UHP 0

            sage: PD.point(I/2)
            Point in PD 1/2*I

            sage: KM.point((0,1))
            Boundary point in KM (0, 1)

            sage: HM.point((0,0,1))
            Point in HM (0, 0, 1)
        """
        return cls.HFactory.get_point(p, **kwargs)

    @classmethod
    def point_in_model(cls, p):
        r"""
        Return ``True`` if ``p`` gives the coordinates of a point in the
        interior of hyperbolic space in the model.

        EXAMPLES::

            sage: UHP.point_in_model(I)
            True
            sage: UHP.point_in_model(0) # Not interior point.
            False
            sage: HM.point_in_model((0,1, sqrt(2)))
            True
        """
        return cls.HModel.point_in_model(p)

    @classmethod
    def bdry_point_in_model(cls, p):
        r"""
        Return ``True`` if ``p`` gives the coordinates of a point on the
        ideal boundary of hyperbolic space in the current model.

        EXAMPLES::

            sage: UHP.bdry_point_in_model(0)
            True
            sage: UHP.bdry_point_in_model(I) # Not boundary point
            False
        """
        return cls.HModel.bdry_point_in_model(p)

    @classmethod
    def isometry_in_model(cls, A):
        r"""
        Return ``True`` if the matrix ``A`` acts isometrically on
        hyperbolic space in the current model.

        EXAMPLES::

            sage: A = matrix(2,[10,0,0,10]) # det(A) is not 1
            sage: UHP.isometry_in_model(A) # But A acts isometrically
            True
        """
        return cls.HModel.isometry_in_model(A)

    @classmethod
    def geodesic(cls, start, end, **kwargs):
        r"""
        Return an oriented :class:`HyperbolicGeodesic` object in the
        current model that starts at ``start`` and ends at ``end``.

        EXAMPLES::

            sage: UHP.geodesic(1, 0)
            Geodesic in UHP from 1 to 0

            sage: PD.geodesic(1, 0)
            Geodesic in PD from 1 to 0

            sage: KM.geodesic((0,1/2), (1/2, 0))
            Geodesic in KM from (0, 1/2) to (1/2, 0)

            sage: HM.geodesic((0,0,1), (0,1, sqrt(2)))
            Geodesic in HM from (0, 0, 1) to (0, 1, sqrt(2))
        """
        return cls.HFactory.get_geodesic(start, end, **kwargs)

    @classmethod
    def isometry(cls, A):
        r"""
        Return an :class:`HyperbolicIsometry` object in the current model
        that coincides with (in the case of linear isometry groups) or lifts
        to (in the case of projective isometry groups) the matrix ``A``.

        EXAMPLES::

            sage: UHP.isometry(identity_matrix(2))
            Isometry in UHP
            [1 0]
            [0 1]

            sage: PD.isometry(identity_matrix(2))
            Isometry in PD
            [1 0]
            [0 1]

            sage: KM.isometry(identity_matrix(3))
            Isometry in KM
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: HM.isometry(identity_matrix(3))
            Isometry in HM
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return cls.HFactory.get_isometry(A)

    @classmethod
    def random_point(cls, **kwargs):
        r"""
        Return a random point in the current model.

        EXAMPLES::

            sage: p = UHP.random_point()
        """
        return cls.HPoint.random_element(**kwargs)

    @classmethod
    def random_geodesic(cls, **kwargs):
        r"""
        Return a random geodesic in the current model.

        EXAMPLES::

            sage: p = UHP.random_geodesic()
        """
        return cls.HGeodesic.random_element(**kwargs)

    @classmethod
    def random_isometry(cls, **kwargs):
        r"""
        Return a random isometry in the current model.

        EXAMPLES::

            sage: p = UHP.random_isometry()
        """
        return cls.HIsometry.random_element(**kwargs)

    @classmethod
    def isometry_from_fixed_points(cls, p1, p2):
        r"""
        Given two ideal boundary points ``p1`` and ``p2``, return an
        isometry of hyperbolic type that fixes ``p1`` and ``p2``.

        INPUT:

        - ``p1``, ``p2`` -- points in the ideal boundary of hyperbolic
          space either as coordinates or as :class:`HyperbolicPoint`.

        OUTPUT:

        - a :class:`HyperbolicIsometry` in the current model whose
          classification is hyperbolic that fixes ``p1`` and ``p2``

        EXAMPLES::

            sage: UHP.isometry_from_fixed_points(0, 4)
            Isometry in UHP
            [  -1    0]
            [-1/5 -1/5]
            sage: UHP.isometry_from_fixed_points(UHP.point(0), UHP.point(4))
            Isometry in UHP
            [  -1    0]
            [-1/5 -1/5]
        """
        return cls.HIsometry.isometry_from_fixed_points(cls.point(p1),
                                                        cls.point(p2))
    @classmethod
    def dist(cls, a, b):
        r"""
        Return the hyperbolic distance between points ``a`` and ``b``.

        INPUT:

        - ``a`` -- a hyperbolic point
        - ``b`` -- a hyperbolic point

        EXAMPLES::

            sage: UHP.dist(UHP.point(I), UHP.point(2*I))
            arccosh(5/4)
            sage: UHP.dist(I, 2*I)
            arccosh(5/4)
        """
        try:
            return a.dist(b)
        except(AttributeError):
            return cls.point(a).dist(cls.point(b))

    @classmethod
    def point_to_model(cls, p, model):
        r"""
        Return the image of ``p`` in the model ``model``.

        INPUT:

        - ``p`` -- a point in the current model of hyperbolic space
          either as coordinates or as a :class:`HyperbolicPoint`

        - ``model`` -- the name of an implemented model of hyperbolic
          space of the same dimension

        EXAMPLES::

            sage: UHP.point_to_model(I, 'PD')
            0
            sage: PD.point_to_model(I, 'UHP')
            +Infinity
            sage: UHP.point_to_model(UHP.point(I), 'HM')
            (0, 0, 1)
        """
        if isinstance(p, HyperbolicPoint):
            p = p.coordinates()
        return cls.HModel.point_to_model(p, model)

    @classmethod
    def isometry_to_model(cls, A, model):
        r"""
        Return the image of ``A`` in the isometry group of the model
        given by ``model``.

        INPUT:

        - ``A`` -- a matrix acting isometrically on the current model
        - ``model`` -- the name of an implemented model of hyperbolic
          space of the same dimension

        EXAMPLES::

            sage: I2 = identity_matrix(2)
            sage: UHP.isometry_to_model(I2, 'HM')
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: UHP.isometry_to_model(UHP.isometry(I2), 'PD')
            [1 0]
            [0 1]
        """
        if isinstance(A, HyperbolicIsometry):
            A = A.matrix()
        return cls.HModel.isometry_to_model(A, model)


class UHP(HyperbolicUserInterface):
    r"""
    Hyperbolic interface for the UHP model.

    EXAMPLES::

        sage: UHP.point(I)
        Point in UHP I
    """
    HModel = HyperbolicModelUHP
    HFactory = HyperbolicFactoryUHP
    HPoint = HyperbolicPointUHP
    HIsometry = HyperbolicIsometryUHP
    HGeodesic = HyperbolicGeodesicUHP


class PD(HyperbolicUserInterface):
    r"""
    Hyperbolic interface for the PD model.

    EXAMPLES::

        sage: PD.point(I)
        Boundary point in PD I
    """
    HModel = HyperbolicModelPD
    HFactory = HyperbolicFactoryPD
    HPoint = HyperbolicPointPD
    HIsometry = HyperbolicIsometryPD
    HGeodesic = HyperbolicGeodesicPD


class KM(HyperbolicUserInterface):
    r"""
    Hyperbolic interface for the KM model.

    EXAMPLES::

        sage: KM.point((0,0))
        Point in KM (0, 0)
    """
    HModel = HyperbolicModelKM
    HFactory = HyperbolicFactoryKM
    HPoint = HyperbolicPointKM
    HIsometry = HyperbolicIsometryKM
    HGeodesic = HyperbolicGeodesicKM


class HM(HyperbolicUserInterface):
    r"""
    Hyperbolic interface for the HM model.

    EXAMPLES::

        sage: HM.point((0,0,1))
        Point in HM (0, 0, 1)
    """
    HModel = HyperbolicModelHM
    HFactory = HyperbolicFactoryHM
    HPoint = HyperbolicPointHM
    HIsometry = HyperbolicIsometryHM
    HGeodesic = HyperbolicGeodesicHM

