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
                                      simplify_chain_real_sympy,
                                      simplify_chain_generic_sympy,)

import sympy
from sympy.core.function import FunctionClass


# conversion functions
def _SR_to_Sympy(expression):
    #print('_SR_to_Sympy', expression)
    # first test if expression is yet Sympy
    # I have not found an elegant way to do that (MMancini)
    if isinstance(type(expression),FunctionClass):
            return expression


    return SR(expression)._sympy_()
    # except:

    #     if isinstance(type(expression),FunctionClass):
    #         return expression
    #     else:
    #         raise
    #         print (type(type(expression)),expression)
    #     return (_Sympy_to_SR(expression))._sympy_()

def _Sympy_to_SR(expression):
    try:
        return SR(expression)
    except:
        # If SR cannot transform a sympy expression is because it is a
        # sympy abstract function.
        a = expression._sage_()
        # As all sage objects have a ._sage_ operator, they have to be
        # catched
        if type(a) == type(expression):
            raise TypeError
        return a

class CalculusMethod(SageObject):
    r"""
    Class for managing the calculus method used in manifolds.

    EXAMPLES:


    """
    _default = 'SR'  # default symbolic method
    _methods = ('SR','sympy') # implemented methods
    _tranf = {'SR':  _Sympy_to_SR,'sympy' : _SR_to_Sympy} # translators

    def __init__(self,current=None,chart=None):
        r"""
        Construct the single instance of class .

        TEST::
            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X = M.chart('x y'); X
            Chart (M, (x, y))
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod(chart=X)
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


        self._simplify_dict = {}

        if chart is not None :
            self._bf_type = chart.manifold().base_field_type()
        else :
            self._bf_type = 'bo'

        if self._bf_type == 'real':
            self._simplify_dict['sympy'] = simplify_chain_real_sympy
            self._simplify_dict['SR'] = simplify_chain_real
        else:
            self._simplify_dict['sympy'] = simplify_chain_generic_sympy
            self._simplify_dict['SR'] = simplify_chain_generic



    def simplify(self,expression,method=None):
        if method is None : method = self._current
        return self._simplify_dict[method](expression)

    # def is_trivial_zero(self,expression,method=None):
    #     if method is None : method = self._current
    #     if method == 'SR':
    #         return expression.is_trivial_zero()
    #     elif method == 'sympy':
    #         return expression.is_zero


    def set(self,method):
        r"""
        Set the current calculus method.

        TEST::
            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X = M.chart('x y'); X
            Chart (M, (x, y))
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod(chart=X)
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
            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X = M.chart('x y'); X
            Chart (M, (x, y))
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod(chart=X)
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
            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X = M.chart('x y'); X
            Chart (M, (x, y))
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod(chart=X)
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
            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X = M.chart('x y'); X
            Chart (M, (x, y))
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod(chart=X)
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
