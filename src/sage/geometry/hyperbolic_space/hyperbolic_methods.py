r"""
Hyperbolic Methods

This module should not be used directly by users.  It is provided for
developers of Sage.

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


class HyperbolicMethodsUHP(HyperbolicAbstractMethods):
    r"""
    Hyperbolic methods for the UHP model of hyperbolic space.
    """
    HModel = HyperbolicModelUHP

