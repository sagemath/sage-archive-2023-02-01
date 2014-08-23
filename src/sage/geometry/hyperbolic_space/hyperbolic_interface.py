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

    sage: HyperbolicPlane.UHP.point(2 + I)
    Point in UHP I + 2

    sage: HyperbolicPlane.PD.point(1/2 + I/2)
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
from sage.structure.parent import Parent
from sage.categories.sets_cat import Sets
from sage.categories.realizations import Realizations, Category_realization_of_parent
from sage.misc.lazy_attribute import lazy_attribute
from sage.geometry.hyperbolic_space.hyperbolic_model import (
        HyperbolicModelUHP, HyperbolicModelPD,
        HyperbolicModelHM, HyperbolicModelKM)

def HyperbolicSpace(n):
    """
    Return ``n`` dimensional hyperbolic space.
    """
    if n == 2:
        return HyperbolicPlane()
    raise NotImplementedError("currently only implemented in dimension 2")

class HyperbolicPlane(Parent, UniqueRepresentation):
    """
    The hyperbolic plane `\mathbb{H}^2`.

    Here are the models currently implemented:

    - ``UHP`` -- upper half plane
    - ``PD`` -- Poincare disk
    - ``KM`` -- Klein disk
    - ``HM`` -- hyperboloid model
    """
    def __init__(self):
        """
        Initialize ``self``.
        """
        Parent.__init__(self, category=Sets().WithRealizations())

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Hyperbolic plane"

    def a_realization(self):
        """
        Return a realization of ``self``.
        """
        return self.UHP()

    UHP = HyperbolicModelUHP
    UpperHalfPlane = UHP

    PD = HyperbolicModelPD
    PoincareDisk = PD

    KM = HyperbolicModelKM
    KleinDisk = KM

    HM = HyperbolicModelHM
    Hyperboloid = HM

class HyperbolicModels(Category_realization_of_parent):
    r"""
    The category of hyperbolic models of hyperbolic space.
    """
    def __init__(self, base):
        r"""
        Initialize the hyperbolic models of hyperbolic space.

        INPUT:

        - ``base`` -- a hyperbolic space

        TESTS::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicModels
            sage: H = HyperbolicPlane()
            sage: models = HyperbolicModels(H)
            sage: H.UHP() in models
            True
        """
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicModels
            sage: H = HyperbolicPlane()
            sage: HyperbolicModels(H)
            Category of hyperbolic models of Hyperbolic plane
        """
        return "Category of hyperbolic models of {}".format(self.base())

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicModels
            sage: H = HyperbolicPlane()
            sage: models = HyperbolicModels(H)
            sage: models.super_categories()
            [Category of finite dimensional algebras with basis over Rational Field,
             Category of realizations of Descent algebra of 4 over Rational Field]
        """
        return [Sets(), Realizations(self.base())]

    class ParentMethods:
        @lazy_attribute
        def _computation_model(self):
            """
            Return the model in which to do computations by default.
            """
            return self.realization_of().UHP()

    class ElementMethods:
        pass

# TODO: Remove this class and move its doctests
class HyperbolicUserInterface(UniqueRepresentation):
    r"""
    Abstract base class for hyperbolic interfaces.  These provide a user
    interface for interacting with models of hyperbolic geometry without
    having the interface dictate the class structure.
    """
    @classmethod
    def model_name(cls):
        r"""
        Return the full name of the hyperbolic model.

        EXAMPLES::

            sage: HyperbolicPlane.UHP.model_name()
            'Upper Half Plane Model'
            sage: HyperbolicPlane.PD.model_name()
            'Poincare Disk Model'
            sage: HyperbolicPlane.KM.model_name()
            'Klein Disk Model'
            sage: HyperbolicPlane.HM.model_name()
            'Hyperboloid Model'
        """
        return cls.HModel.name

    @classmethod
    def short_name(cls):
        r"""
        Return the short name of the hyperbolic model.

        EXAMPLES::

            sage: HyperbolicPlane.UHP.short_name()
            'UHP'
            sage: HyperbolicPlane.PD.short_name()
            'PD'
            sage: HyperbolicPlane.HM.short_name()
            'HM'
            sage: HyperbolicPlane.KM.short_name()
            'KM'
        """
        return cls.HModel.short_name

    @classmethod
    def is_bounded(cls):
        r"""
        Return ``True`` if the model is bounded and ``False`` otherwise.

        EXAMPLES::

            sage: HyperbolicPlane.UHP.is_bounded()
            True
            sage: HyperbolicPlane.PD.is_bounded()
            True
            sage: HyperbolicPlane.KM.is_bounded()
            True
            sage: HyperbolicPlane.HM.is_bounded()
            False
        """
        return cls.HModel.bounded

    @classmethod
    def point(cls, p, **kwargs):
        r"""
        Return a :class:`HyperbolicPoint` object in the current
        model with coordinates ``p``.

        EXAMPLES::

            sage: HyperbolicPlane.UHP.point(0)
            Boundary point in UHP 0

            sage: HyperbolicPlane.PD.point(I/2)
            Point in PD 1/2*I

            sage: HyperbolicPlane.KM.point((0,1))
            Boundary point in KM (0, 1)

            sage: HyperbolicPlane.HM.point((0,0,1))
            Point in HM (0, 0, 1)
        """
        return cls.HFactory.get_point(p, **kwargs)

    @classmethod
    def point_in_model(cls, p):
        r"""
        Return ``True`` if ``p`` gives the coordinates of a point in the
        interior of hyperbolic space in the model.

        EXAMPLES::

            sage: HyperbolicPlane.UHP.point_in_model(I)
            True
            sage: HyperbolicPlane.UHP.point_in_model(0) # Not interior point.
            False
            sage: HyperbolicPlane.HM.point_in_model((0,1, sqrt(2)))
            True
        """
        return cls.HModel.point_in_model(p)

    @classmethod
    def bdry_point_in_model(cls, p):
        r"""
        Return ``True`` if ``p`` gives the coordinates of a point on the
        ideal boundary of hyperbolic space in the current model.

        EXAMPLES::

            sage: HyperbolicPlane.UHP.bdry_point_in_model(0)
            True
            sage: HyperbolicPlane.UHP.bdry_point_in_model(I) # Not boundary point
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
            sage: HyperbolicPlane.UHP.isometry_in_model(A) # But A acts isometrically
            True
        """
        return cls.HModel.isometry_in_model(A)

    @classmethod
    def geodesic(cls, start, end, **kwargs):
        r"""
        Return an oriented :class:`HyperbolicGeodesic` object in the
        current model that starts at ``start`` and ends at ``end``.

        EXAMPLES::

            sage: HyperbolicPlane.UHP.geodesic(1, 0)
            Geodesic in UHP from 1 to 0

            sage: HyperbolicPlane.PD.geodesic(1, 0)
            Geodesic in PD from 1 to 0

            sage: HyperbolicPlane.KM.geodesic((0,1/2), (1/2, 0))
            Geodesic in KM from (0, 1/2) to (1/2, 0)

            sage: HyperbolicPlane.HM.geodesic((0,0,1), (0,1, sqrt(2)))
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

            sage: HyperbolicPlane.UHP.isometry(identity_matrix(2))
            Isometry in UHP
            [1 0]
            [0 1]

            sage: HyperbolicPlane.PD.isometry(identity_matrix(2))
            Isometry in PD
            [1 0]
            [0 1]

            sage: HyperbolicPlane.KM.isometry(identity_matrix(3))
            Isometry in KM
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: HyperbolicPlane.HM.isometry(identity_matrix(3))
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

            sage: p = HyperbolicPlane.UHP.random_point()
        """
        return cls.HPoint.random_element(**kwargs)

    @classmethod
    def random_geodesic(cls, **kwargs):
        r"""
        Return a random geodesic in the current model.

        EXAMPLES::

            sage: p = HyperbolicPlane.UHP.random_geodesic()
        """
        return cls.HGeodesic.random_element(**kwargs)

    @classmethod
    def random_isometry(cls, **kwargs):
        r"""
        Return a random isometry in the current model.

        EXAMPLES::

            sage: p = HyperbolicPlane.UHP.random_isometry()
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

            sage: UHP = HyperbolicPlane.UHP
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

            sage: UHP = HyperbolicPlane.UHP
            sage: UHP.dist(UHP.point(I), UHP.point(2*I))
            arccosh(5/4)
            sage: HyperbolicPlane.UHP.dist(I, 2*I)
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

            sage: HyperbolicPlane.UHP.point_to_model(I, 'PD')
            0
            sage: HyperbolicPlane.PD.point_to_model(I, 'UHP')
            +Infinity
            sage: HyperbolicPlane.UHP.point_to_model(HyperbolicPlane.UHP.point(I), 'HM')
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

            sage: UHP = HyperbolicPlane.UHP
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
