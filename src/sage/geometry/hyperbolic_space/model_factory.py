r"""
Factory for Hyperbolic Models

AUTHORS:

- Greg Laun (2013): initial version
"""

#*****************************************************************************
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.lazy_import import lazy_import

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_model',
            ['HyperbolicModelUHP', 'HyperbolicModelPD',
             'HyperbolicModelHM', 'HyperbolicModelKM'])

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_factory',
             ['HyperbolicFactoryUHP', 'HyperbolicFactoryPD',
              'HyperbolicFactoryHM', 'HyperbolicFactoryKM'])

class ModelFactory(UniqueRepresentation):
    """
    Factory for creating the hyperbolic models.
    """
    @classmethod
    def find_model(cls, model_name):
        r"""
        Given the short name of a hyperbolic model, return that model.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.model_factory import ModelFactory
            sage: ModelFactory.find_model('UHP')
            <class 'sage.geometry.hyperbolic_space.hyperbolic_model.HyperbolicModelUHP'>
        """
        return {
            'UHP': HyperbolicModelUHP,
            'PD' : HyperbolicModelPD,
            'HM' : HyperbolicModelHM,
            'KM' : HyperbolicModelKM
            }[model_name]

    @classmethod
    def find_factory(cls, model_name):
        r"""
        Given the short name of a hyperbolic model, return the associated factory.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.model_factory import ModelFactory
            sage: ModelFactory.find_factory('UHP')
            <class 'sage.geometry.hyperbolic_space.hyperbolic_factory.HyperbolicFactoryUHP'>
        """
        return {
            'UHP' : HyperbolicFactoryUHP,
            'PD' : HyperbolicFactoryPD,
            'HM' : HyperbolicFactoryHM,
            'KM' : HyperbolicFactoryKM
            }[model_name]

