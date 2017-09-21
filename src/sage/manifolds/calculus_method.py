r"""
Calculus method control

This module defines the class :class:`CalculusMethod` to govern the
calculus (symbolic or numerical) method used in manifolds.


AUTHORS:

- Marco Mancini, (2017): initial version

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


# conversion functions
def _SR_to_Sympy(expression):
    r"""
    Convert an expression from SR to ``sympy``. In the case that the
    expression is yet a ``sympy`` expression then it is returned

    INPUT:

    - ``expression`` -- SR or ``sympy`` symbolic expression

    OUTPUT:

    - ``expression`` -- ``sympy`` symbolic expression

    EXAMPLES::
        sage: from sage.manifolds.calculus_method import _SR_to_Sympy
        sage: a = x^2+sin(x)^2; a
        x^2 + sin(x)^2
        sage: type(a)
        <type 'sage.symbolic.expression.Expression'>
        sage: b = _SR_to_Sympy(a); b
        x**2 + sin(x)**2
        sage: type(b)
        <class 'sympy.core.add.Add'>

        The conversion of a ``sympy`` expression is not converted::

        sage: _SR_to_Sympy(b) == b
        True


    """
    from sympy.core.function import FunctionClass
    # first test if expression is yet Sympy
    # I have not found an elegant way to do that (MMancini)
    if isinstance(type(expression),FunctionClass):
            return expression

    return SR(expression)._sympy_()


def _Sympy_to_SR(expression):
    r"""
    Convert an expression from ``sympy`` to SR. In the case that the
    expression is yet a SR expression then it is returned

    INPUT:

    - ``expression`` -- ``sympy`` symbolic expression


    OUTPUT:

    - ``expression`` -- SR or ``sympy`` symbolic expression


    EXAMPLES::

        sage: from sage.manifolds.calculus_method import _Sympy_to_SR, _SR_to_Sympy
        sage: a = x^2+sin(x)^2;
        sage: type(a)
        <type 'sage.symbolic.expression.Expression'>
        sage: b = _SR_to_Sympy(a); b
        x**2 + sin(x)**2
        sage: type(b)
        <class 'sympy.core.add.Add'>
        sage: _Sympy_to_SR(b)==a
        x^2 + sin(x)^2 == x^2 + sin(x)^2
        sage: bool(_Sympy_to_SR(b)==a)
        True

    """
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
    Calculus method used in manifolds.
    The class store the possible calculus methods and permits to
    select some basic operations working on them.
    For the moment, only two calculus methods are implemented :
     - `Maxima`, indicated with `SR` (Symbolic Ring)
     - `Sympy`, indicated with 'sympy'.
    The current method is indicated with a `*`.


    EXAMPLES::

        sage: from sage.manifolds.calculus_method import CalculusMethod
        sage: calc_mth = CalculusMethod(bf_type='complex')
        sage: calc_mth._repr_()
        'Possible symbolic methods:\n - SR (*) (default)\n - sympy'

        sage: calc_mth.set('sympy')
        sage: print (calc_mth)
        Possible symbolic methods:
         - SR (default)
         - sympy (*)
        sage: calc_mth.reset()



    """
    _default = 'SR'  # default symbolic method
    _methods = ('SR','sympy') # implemented methods
    _tranf = {'SR':  _Sympy_to_SR,'sympy' : _SR_to_Sympy} # translators

    def __init__(self, current=None, bf_type='complex'):
        r"""
        Construct the single instance of the class.

        INPUT:

        - ``current`` -- Current symbolic method (default: `None`)

        - ``bf_type`` -- base field (default: `complex`)


        """

        self._current = self._default if current is None else current

        self._simplify_dict = {}

        # initiliazation of the dict for simplifications
        if bf_type == 'real':
            self._simplify_dict['sympy'] = simplify_chain_real_sympy
            self._simplify_dict['SR'] = simplify_chain_real
        else:
            self._simplify_dict['sympy'] = simplify_chain_generic_sympy
            self._simplify_dict['SR'] = simplify_chain_generic

        self._bf_type = bf_type


    def simplify(self, expression, method=None):
        r"""
        Simplification chain for calculus method.

        INPUT:

        - ``expression`` -- expression to simplify

        - ``method`` -- calculus method to use (default: `None`)


        EXAMPLES::

            sage: M = Manifold(2, 'M', field='complex', structure='topological')
            sage: X = M.chart('x y');
            sage: f = x^2+sin(x)^2+cos(x)^2
            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod(bf_type='real')
            sage: calc_mth.simplify(f)
            x^2 + 1

        Method cannot be mixed::

        sage: calc_mth.set('sympy')
            sage: calc_mth.simplify(f)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.symbolic.expression.Expression' object has no attribute 'combsimp'

        sage: calc_mth.simplify(f._sympy_())
        x**2 + 1

        or::

        sage: calc_mth.simplify(f,method="SR")
        x^2 + 1

        """

        if method is None : method = self._current
        return self._simplify_dict[method](expression)

    def is_trivial_zero(self, expression, method=None):
        r"""
        Chek if expression is trivially equal to zero without any
        simplification.

        INPUT:

        - ``expression`` -- expression

        - ``method`` -- calculus method to use (default: `None`)

        OUTPUT:

        - `True` is expression is trivially zero, `False` elsewhere.


        EXAMPLES::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod(bf_type='real')
            sage: f = sin(x)-sin(x)
            sage: calc_mth.is_trivial_zero(f)
            True
            sage: f = sin(x)^2 + cos(x)^2 - 1
            sage: calc_mth.is_trivial_zero(f)
            False


        """

        if method is None : method = self._current
        if method == 'SR':
            return expression.is_trivial_zero()
        elif method == 'sympy':
            # we have to test sympy is zero because it could be 'NoneType'
            if expression.is_zero:
                return True
            else:
                False


    def set(self,method):
        r"""
        Set the current calculus method.

        - ``method`` -- calculus method

        EXAMPLES::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod(bf_type='complex')
            sage: calc_mth.set('sympy')
            sage: calc_mth
            Possible symbolic methods:
             - SR (default)
             - sympy (*)
            sage: calc_mth.reset()
        """

        self._test(method)
        self._current = method

    def _test(self,method):
        r"""
        Test if a calculus method is defined.

        INPUT:

        - ``method`` -- calculus method

        TESTS::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod(bf_type='complex')
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


        EXAMPLES::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod(bf_type='complex')
            sage: calc_mth.set('sympy')
            sage: calc_mth
            Possible symbolic methods:
             - SR (default)
             - sympy (*)
            sage: calc_mth.reset()
            sage: calc_mth
             Possible symbolic methods:
             - SR (*) (default)
             - sympy

        """
        self._current = self._default

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: from sage.manifolds.calculus_method import CalculusMethod
            sage: calc_mth = CalculusMethod(bf_type='complex')
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
