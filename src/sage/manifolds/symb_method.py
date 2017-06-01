r"""
Symbolic method control

This module defines the singleton class :class:`symb_method` to govern the
symbolic method used in manifolds.


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
from __future__ import absolute_import

from sage.structure.sage_object import SageObject
from sage.misc.fast_methods import Singleton
from sage.symbolic.ring import SR
import sympy


# conversion functions
def _SR_to_Sympy(expression):
    return SR(expression)._sympy_()

def _Sympy_to_SR(expression):
    return SR(expression)


class SymbMethod(Singleton, SageObject):
    r"""
    Singleton class for managing the symbolic method used in manifolds.

    EXAMPLES:


    """

    def __init__(self,current=None):
        r"""
        Construct the single instance of class Parallelism (singleton model).

        TEST::

            sage: symb = SymbMethod()
            sage: symb._repr_()
            'Possible symbolic methods:\n - SR (*) (default)\n - sympy'

            sage: symb.set('sympy')
            sage: symb._repr_()
            'Possible symbolic methods:\n - SR (default)\n - sympy (*)'
            sage: symb.reset()

        Test of the singleton character::


        The test suite is passed::


        """
        self._default = 'SR'  # default symbolic method
        self._methods = ('SR','sympy') # implemented methods
        self.tranf = {'SR': SR,'sympy' : _SR_to_Sympy} # translators
        self.current = self._default if current is None else current


    def set(self,method):
        self.test(method)
        self.current = method

    def test(self,method):
        if method not in self._methods:
            raise NotImplementedError("method "+method+" not " +
                                  "implemented yet.")

    def reset(self):
        r"""
        Set the current symbolic method to default one.


        TEST::
            sage: symb = SymbMethod()
            sage: symb.set('sympy')
            sage: symb._repr_()
            'Possible symbolic methods:\n - SR (default)\n - sympy (*)'
            sage: symb.reset()
            sage: symb._repr_()
            'Possible symbolic methods:\n - SR (*) (default)\n - sympy'


        """
        self.current = self._default

    def _repr_(self):
        r"""
        String representation of the object.

        TEST::
            sage: SymbMethod()._repr_()
            'Possible symbolic methods:\n - SR (*) (default)\n - sympy'

        """
        resu = 'Possible symbolic methods:\n'
        for method in self._methods:
            resu += ' - {}'.format(method)
            if method == self.current: resu += ' (*)'
            if method == self._default: resu += ' (default)'
            resu += '\n'
        return resu[:-1]
