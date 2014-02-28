r"""
AUTHORS:
- Greg Laun (2013): initial version


#***********************************************************************
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************
"""
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.lazy_import import lazy_import
lazy_import('sage.functions.other','sqrt')

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_model', 'HyperbolicModel')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_point', 'HyperbolicPoint')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_bdry_point', 'HyperbolicBdryPoint')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_isometry', 'HyperbolicIsometry')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_geodesic', 'HyperbolicGeodesic')

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_model', 'HyperbolicModelUHP')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_point', 'HyperbolicPointUHP')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_bdry_point', 'HyperbolicBdryPointUHP')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_isometry', 'HyperbolicIsometryUHP')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_geodesic', 'HyperbolicGeodesicUHP')

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_model', 'HyperbolicModelPD')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_point', 'HyperbolicPointPD')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_bdry_point', 'HyperbolicBdryPointPD')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_isometry', 'HyperbolicIsometryPD')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_geodesic', 'HyperbolicGeodesicPD')

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_model', 'HyperbolicModelKM')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_point', 'HyperbolicPointKM')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_bdry_point', 'HyperbolicBdryPointKM')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_isometry', 'HyperbolicIsometryKM')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_geodesic', 'HyperbolicGeodesicKM')

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_model', 'HyperbolicModelHM')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_point', 'HyperbolicPointHM')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_bdry_point', 'HyperbolicBdryPointHM')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_isometry', 'HyperbolicIsometryHM')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_geodesic', 'HyperbolicGeodesicHM')

class HyperbolicAbstractFactory (UniqueRepresentation):
    HModel = HyperbolicModel
    HPoint = HyperbolicPoint
    HBdryPoint = HyperbolicBdryPoint
    HIsometry = HyperbolicIsometry
    HGeodesic = HyperbolicGeodesic

    @classmethod
    def get_model(cls):
        r"""
        Return the current model as a class rather than a string.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: HyperbolicFactoryUHP.get_model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelUHP'>

            sage: HyperbolicFactoryPD.get_model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelPD'>

            sage: HyperbolicFactoryKM.get_model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelKM'>

            sage: HyperbolicFactoryHM.get_model()
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelHM'>

        """
        return cls.HModel

    @classmethod
    def get_interior_point(cls, coordinates, **graphics_options):
        r"""
        Return a point in the appropriate model given the coordinates.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: HyperbolicFactoryUHP.get_interior_point(2 + 3*I)
            Point in UHP 3*I + 2.

            sage: HyperbolicFactoryPD.get_interior_point(0)
            Point in PD 0.

            sage: HyperbolicFactoryKM.get_interior_point((0,0))
            Point in KM (0, 0).

            sage: HyperbolicFactoryHM.get_interior_point((0,0,1))
            Point in HM (0, 0, 1).

            sage: p = HyperbolicFactoryUHP.get_interior_point(I, color="red"); p.graphics_options()
            {'color': 'red'}
        """
        return cls.HPoint(coordinates, **graphics_options)

    @classmethod
    def get_bdry_point(cls, coordinates, **graphics_options):
        r"""
        Return a boundary point in the appropriate model given the
        coordinates.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: HyperbolicFactoryUHP.get_bdry_point(12)
            Boundary point in UHP 12.

            sage: HyperbolicFactoryUHP.get_bdry_point(infinity)
            Boundary point in UHP +Infinity.

            sage: HyperbolicFactoryPD.get_bdry_point(I)
            Boundary point in PD I.

            sage: HyperbolicFactoryKM.get_bdry_point((0,-1))
            Boundary point in KM (0, -1).
        """
        return cls.HBdryPoint(coordinates, **graphics_options)

    @classmethod
    def get_point(cls, coordinates, **graphics_options):
        r"""
        Automatically determine the type of point to return given either
        (1) the coordinates of a point in the interior or ideal boundary
        of hyperbolic space or (2) a HyperbolicPoint or
        HyperbolicBdryPoint object.

        INPUT:

        - a point in hyperbolic space or on the ideal boundary.

        OUTPUT:

        - A HyperbolicPoint or HyperbolicBdryPoint.

        EXAMPLES::

        We can create an interior point via the coordinates::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: p = HyperbolicFactoryUHP.get_point(2*I); p
            Point in UHP 2*I.

        Or we can create a boundary point via the coordinates::

            sage: q = HyperbolicFactoryUHP.get_point(23); q
            Boundary point in UHP 23.

        Or we can create both types of points by passing the
        HyperbolicPoint or HyperbolicBdryPoint object::

            sage: HyperbolicFactoryUHP.get_point(p)
            Point in UHP 2*I.

            sage: HyperbolicFactoryUHP.get_point(q)
            Boundary point in UHP 23.

            sage: HyperbolicFactoryUHP.get_point(12 - I)
            Traceback (most recent call last):
            ...
            ValueError: -I + 12 is neither an interior nor boundary point in the UHP model.
        """
        from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPoint
        if isinstance(coordinates, HyperbolicPoint):
            coordinates.update_graphics(True, **graphics_options)
            return coordinates #both Point and BdryPoint
        elif cls.HModel.point_in_model(coordinates):
            return cls.HPoint(coordinates, **graphics_options)
        elif cls.HModel.bdry_point_in_model(coordinates):
            return cls.HBdryPoint(coordinates, **graphics_options)
        else:
            e_1 = "{0} is neither an interior nor boundary".format(coordinates)
            e_2 = " point in the {0} model.".format(cls.get_model().short_name)
            raise ValueError(e_1 + e_2)

    @classmethod
    def get_geodesic(cls, start, end, **graphics_options):
        r"""
        Return a geodesic in the appropriate model.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: HyperbolicFactoryUHP.get_geodesic(I, 2*I)
            Geodesic in UHP from I to 2*I.

            sage: HyperbolicFactoryPD.get_geodesic(0, I/2)
            Geodesic in PD from 0 to 1/2*I.

            sage: HyperbolicFactoryKM.get_geodesic((1/2, 1/2), (0,0))
            Geodesic in KM from (1/2, 1/2) to (0, 0).

            sage: HyperbolicFactoryHM.get_geodesic((0,0,1), (1,0, sqrt(2)))
            Geodesic in HM from (0, 0, 1) to (1, 0, sqrt(2)).
        """
        return cls.HGeodesic(start, end, **graphics_options)

    @classmethod
    def get_isometry(cls, A):
        r"""
        Return an isometry in the appropriate model given the matrix.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: HyperbolicFactoryUHP.get_isometry(identity_matrix(2))
            Isometry in UHP
            [1 0]
            [0 1].

            sage: HyperbolicFactoryPD.get_isometry(identity_matrix(2))
            Isometry in PD
            [1 0]
            [0 1].

            sage: HyperbolicFactoryKM.get_isometry(identity_matrix(3))
            Isometry in KM
            [1 0 0]
            [0 1 0]
            [0 0 1].

            sage: HyperbolicFactoryHM.get_isometry(identity_matrix(3))
            Isometry in HM
            [1 0 0]
            [0 1 0]
            [0 0 1].
        """
        return cls.HIsometry(A)

    @classmethod
    def get_background_graphic(cls, **bdry_options):
        r"""
        Return a graphic object that makes the model easier to
        visualize.  For bounded models, such as the upper half space and
        the unit ball models, the background object is the ideal
        boundary.  For the hyperboloid model, the background object is
        the hyperboloid itself.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: circ = HyperbolicFactoryPD.get_background_graphic()
        """
        return None

class HyperbolicFactoryUHP (HyperbolicAbstractFactory, UniqueRepresentation):
    HModel = HyperbolicModelUHP
    HPoint = HyperbolicPointUHP
    HBdryPoint = HyperbolicBdryPointUHP
    HIsometry = HyperbolicIsometryUHP
    HGeodesic = HyperbolicGeodesicUHP

    @classmethod
    def get_background_graphic(cls, **bdry_options):
        r"""
        Return a graphic object that makes the model easier to visualize.
        For bounded models, such as the upper half space and the unit
        ball models, the background object is the ideal boundary.  For
        the hyperboloid model, the background object is the hyperboloid
        itself.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: circ = HyperbolicFactoryUHP.get_background_graphic()
        """
        from sage.plot.line import line
        bd_min = bdry_options.get('bd_min', -5)
        bd_max = bdry_options.get('bd_max', 5)
        return line(((bd_min, 0), (bd_max, 0)), color='black')

class HyperbolicFactoryPD (HyperbolicAbstractFactory, UniqueRepresentation):
    from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelPD
    HModel = HyperbolicModelPD
    HPoint = HyperbolicPointPD
    HBdryPoint = HyperbolicBdryPointPD
    HIsometry = HyperbolicIsometryPD
    HGeodesic = HyperbolicGeodesicPD

    @classmethod
    def get_background_graphic(cls, **bdry_options):
        r"""
        Return a graphic object that makes the model easier to visualize.
        For bounded models, such as the upper half space and the unit
        ball models, the background object is the ideal boundary.  For
        the hyperboloid model, the background object is the hyperboloid
        itself.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: circ = HyperbolicFactoryPD.get_background_graphic()
        """
        from sage.plot.circle import circle
        return circle((0,0),1, axes=False, color = 'black')


class HyperbolicFactoryKM (HyperbolicAbstractFactory, UniqueRepresentation):
    from sage.geometry.hyperbolic_space.hyperbolic_model import HyperbolicModelKM
    HModel = HyperbolicModelKM
    HPoint = HyperbolicPointKM
    HBdryPoint = HyperbolicBdryPointKM
    HIsometry = HyperbolicIsometryKM
    HGeodesic = HyperbolicGeodesicKM

    @classmethod
    def get_background_graphic(cls, **bdry_options):
        r"""
        Return a graphic object that makes the model easier to visualize.
        For bounded models, such as the upper half space and the unit
        ball models, the background object is the ideal boundary.  For
        the hyperboloid model, the background object is the hyperboloid
        itself.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: circ = HyperbolicFactoryKM.get_background_graphic()
        """
        from sage.plot.circle import circle
        return circle((0,0),1, axes=False, color = 'black')

class HyperbolicFactoryHM (HyperbolicAbstractFactory, UniqueRepresentation):
    HModel = HyperbolicModelHM
    HPoint = HyperbolicPointHM
    HIsometry = HyperbolicIsometryHM
    HGeodesic = HyperbolicGeodesicHM
    
    @classmethod
    def get_background_graphic(cls, **bdry_options):
        r"""
        Return a graphic object that makes the model easier to visualize.
        For bounded models, such as the upper half space and the unit
        ball models, the background object is the ideal boundary.  For
        the hyperboloid model, the background object is the hyperboloid
        itself.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_factory import *
            sage: circ = HyperbolicFactoryPD.get_background_graphic()
        """
        hyperboloid_opacity = bdry_options.get('hyperboloid_opacity', .1)
        z_height = bdry_options.get('z_height', 7.0)
        x_max = sqrt((z_height**2 - 1)/2.0)
        from sage.plot.plot3d.all import plot3d
        from sage.all import var
        (x,y) = var('x,y')
        return plot3d((1 + x**2 + y**2).sqrt(), (x, -x_max, x_max),\
                          (y,-x_max, x_max), opacity = hyperboloid_opacity, **bdry_options)
