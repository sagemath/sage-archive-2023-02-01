r"""
Calculus method control

This module defines the class :class:`CalculusMethod` to govern the
calculus (symbolic or numerical) method used in manifolds.


AUTHORS:

- Marco Mancini, Eric Gourgoulhon, Michal Bejger (2017): initial version

"""

#******************************************************************************
#       Copyright (C) 2015 Marco Mancini <marco.mancini@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.structure.sage_object import SageObject
from sage.symbolic.ring import SR
from sage.manifolds.utilities import (simplify_chain_real,
                                      simplify_chain_generic,
                                      simplify_chain_generic_sympy,)

import sympy


# conversion functions
def _SR_to_Sympy(expression):
    return SR(expression)._sympy_()

def _Sympy_to_SR(expression):
    return SR(expression)


class CalculusMethod(SageObject):
    r"""
    Class for managing the calculus method used in manifolds.

    EXAMPLES:


    """
    _default = 'SR'  # default symbolic method
    _methods = ('SR','sympy') # implemented methods
    _tranf = {'SR': SR,'sympy' : _SR_to_Sympy} # translators

    def __init__(self,current=None,chart=None):
        r"""
        Construct the single instance of class .

        TEST::
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod()
            sage: calc_mth._repr_()
            'Possible symbolic methods:\n - SR (*) (default)\n - sympy'

            sage: calc_mth.set('sympy')
            sage: calc_mth._repr_()
            'Possible symbolic methods:\n - SR (default)\n - sympy (*)'
            sage: calc_mth.reset()

        Test of the singleton character::


        The test suite is passed::


        """

        self._current = self._default if current is None else current


        self._simplify = {'sympy': simplify_chain_generic_sympy }

        if chart.manifold().base_field_type() == 'real':
            self._simplify['SR'] = simplify_chain_real
        else:
            self._simplify['SR'] = simplify_chain_generic


    def set(self,method):
        r"""
        Set the current calculus method.

        TEST::
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod()
            sage: calc_mth.set('sympy')
            sage: calc_mth._repr_()
            'Possible symbolic methods:\n - SR (default)\n - sympy (*)'
            sage: calc_mth.reset()
        """

        self._test(method)
        self._current = method

    def _test(self,method):
        r"""
        Test if a calculus method is defined.

        TEST::
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod()
            sage: calc_mth._test('lala')
            Traceback (most recent call last):
            ...
            NotImplementedError: method lala not implemented yet.

        """

        if method not in self._methods:
            raise NotImplementedError("method "+method+" not " +
                                  "implemented yet.")

    def reset(self):
        r"""
        Set the current symbolic method to default one.


        TEST::
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod()
            sage: calc_mth.set('sympy')
            sage: calc_mth._repr_()
            'Possible symbolic methods:\n - SR (default)\n - sympy (*)'
            sage: calc_mth.reset()
            sage: calc_mth._repr_()
            'Possible symbolic methods:\n - SR (*) (default)\n - sympy'


        """
        self._current = self._default

    def _repr_(self):
        r"""
        String representation of the object.

        TEST::
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod()
            sage: calc_mth._repr_()
            'Possible symbolic methods:\n - SR (*) (default)\n - sympy'

        """
        resu = 'Possible symbolic methods:\n'
        for method in self._methods:
            resu += ' - {}'.format(method)
            if method == self._current: resu += ' (*)'
            if method == self._default: resu += ' (default)'
            resu += '\n'
        return resu[:-1]
