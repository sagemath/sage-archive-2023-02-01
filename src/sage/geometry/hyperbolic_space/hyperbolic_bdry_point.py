"""
Ideal Boundary Points in Hyperbolic Space

This module implements the abstract base class for ideal points in 
hyperbolic space of arbitrary dimension.  It also contains the
implementations for specific models of hyperbolic geometry.

Note that not all models of hyperbolic space are bounded, meaning that
the ideal boundary is not the topological boundary of the set underlying
tho model.  For example, the unit disk model is bounded with boundary
given by the unit sphere.  The hyperboloid model is not bounded.

AUTHORS:

- Greg Laun (2013): initial version

EXAMPLES:

We can construct boundary points in the upper half plane model, 
abbreviated UHP for convenience::

    sage: UHP.point(3)
    Boundary point in UHP 3

Points on the boundary are infinitely far from interior points::

    sage: UHP.point(3).dist(UHP.point(I))
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

from sage.geometry.hyperbolic_space.hyperbolic_point import HyperbolicPoint
from sage.symbolic.pynac import I
from sage.rings.infinity import infinity
from sage.rings.all import RR
from sage.misc.lazy_import import lazy_import

lazy_import('sage.plot.point', 'point')
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_factory',
            ['HyperbolicAbstractFactory', 'HyperbolicFactoryHM',
             'HyperbolicFactoryUHP', 'HyperbolicFactoryPD',
             'HyperbolicFactoryKM'])
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_methods',
            ['HyperbolicAbstractMethods', 'HyperbolicMethodsUHP'])
lazy_import('sage.modules.free_module_element', 'vector')


class HyperbolicBdryPoint(HyperbolicPoint):
    r"""
    Abstract base class for points on the ideal boundary of hyperbolic
    space.  This class should never be instantiated.

    INPUT:

    - the coordinates of a hyperbolic boundary point in the appropriate model

    EXAMPLES:

    Note that the coordinate does not differentiate the different models::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import *
        sage: p = HyperbolicBdryPointUHP(1); p
        Boundary point in UHP 1

        sage: q = HyperbolicBdryPointPD(1); q
        Boundary point in PD 1

        sage: p == q
        False

        sage: bool(p.coordinates() == q.coordinates())
        True

    It is an error to specify a an interior point of hyperbolic space::

        sage: HyperbolicBdryPointUHP(0.2 + 0.3*I)
        Traceback (most recent call last):
        ...
        ValueError: 0.200000000000000 + 0.300000000000000*I is not a valid boundary point in the UHP model
    """
    def __init__(self, coordinates, **graphics_options):
        r"""
        See ``HyperbolicPoint`` for full documentation.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import *
            sage: HyperbolicBdryPointUHP(1)
            Boundary point in UHP 1
        """
        if not self.model().bounded:
            raise NotImplementedError(
                "{0} is not a bounded model; boundary"
                " points not implemented".format(self.model_name()))
        elif self.model().bdry_point_in_model(coordinates):
            if type(coordinates) == tuple:
                coordinates = vector(coordinates)
            self._coordinates = coordinates
        else:
            raise ValueError(
                "{0} is not a valid".format(coordinates) +
                " boundary point in the {0} model".format(self.model_name()))
        self._graphics_options = graphics_options

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import *
            sage: HyperbolicBdryPointUHP(infinity)
            Boundary point in UHP +Infinity

            sage: HyperbolicBdryPointPD(-1)
            Boundary point in PD -1

            sage: HyperbolicBdryPointKM((0, -1))
            Boundary point in KM (0, -1)
        """
        return "Boundary point in {0} {1}".format(self.model_name(),
                                                   self.coordinates())


class HyperbolicBdryPointUHP(HyperbolicBdryPoint):
    r"""
    Create a boundary point for the UHP model.

    INPUT:

    - the coordinates of a hyperbolic boundary point in the upper half plane

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import HyperbolicBdryPointUHP
        sage: q = HyperbolicBdryPointUHP(1); q
        Boundary point in UHP 1
    """
    HFactory = HyperbolicFactoryUHP
    HMethods = HyperbolicMethodsUHP

    def show(self, boundary=True, **options):
        r"""
        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import *
            sage: HyperbolicBdryPointUHP(0).show()
            sage: HyperbolicBdryPointUHP(infinity).show()
            Traceback (most recent call last):
            ...
            NotImplementedError: can't draw the point infinity
        """
        opts = dict([('axes', False), ('aspect_ratio',1)])
        opts.update(self.graphics_options())
        opts.update(options)
        from sage.misc.functional import numerical_approx
        p = self.coordinates()
        if p == infinity:
            raise NotImplementedError("can't draw the point infinity")
        p = numerical_approx(p)
        pic = point((p, 0), **opts)
        if boundary:
            bd_pic = self.HFactory.get_background_graphic(bd_min = p - 1,
                                                          bd_max = p + 1)
            pic = bd_pic + pic
        return pic


class HyperbolicBdryPointPD(HyperbolicBdryPoint):
    r"""
    Create a boundary point for the PD model.

    INPUT:

    - the coordinates of a hyperbolic boundary point in the Poincare disk

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import HyperbolicBdryPointPD
        sage: q = HyperbolicBdryPointPD(1); q
        Boundary point in PD 1
    """
    HFactory = HyperbolicFactoryPD
    HMethods = HyperbolicMethodsUHP


class HyperbolicBdryPointKM(HyperbolicBdryPoint):
    r"""
    Create a boundary point for the KM model.

    INPUT:

    - the coordinates of a hyperbolic boundary point in the Klein disk

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import HyperbolicBdryPointKM
        sage: q = HyperbolicBdryPointKM((1,0)); q
        Boundary point in KM (1, 0)
    """
    HFactory = HyperbolicFactoryKM
    HMethods = HyperbolicMethodsUHP

class HyperbolicBdryPointHM(HyperbolicBdryPoint):
    r"""
    A dummy class for the boundary points of the hyperboloid model.  The model
    is not bounded, so there are no boundary points.  The class is needed for
    compatibility reasons.

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_bdry_point import HyperbolicBdryPointHM
        sage: q = HyperbolicBdryPointHM((1,0,0)); q
        Traceback (most recent call last):
        ...
        NotImplementedError: HM is not a bounded model; boundary points not implemented
    """
    HFactory = HyperbolicFactoryHM
    HMethods = HyperbolicMethodsUHP

