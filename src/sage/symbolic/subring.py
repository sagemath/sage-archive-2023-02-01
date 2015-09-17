r"""
Subrings of the Symbolic Ring

AUTHORS:

- Daniel Krenn (2015)


Classes and Methods
===================
"""

#*****************************************************************************
# Copyright (C) 2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage
from ring import SymbolicRing, SR
from expression import Expression


class SymbolicSubringFactory(sage.structure.factory.UniqueFactory):
    r"""
    A factory creating a symbolic subring.

    EXAMPLES::

        sage: from sage.symbolic.subring import SymbolicSubring
        sage: V = var('a, b, c, r, s, t, x, y, z')
        sage: def var_in_subring(s, S):
        ....:     try:
        ....:         S(s)
        ....:         return True
        ....:     except ValueError:
        ....:         return False

    ::

        sage: A = SymbolicSubring(accepting_variables=(a, b, c)); A
        Symbolic Subring accepting variables a, b, c
        sage: tuple((v, var_in_subring(v, A)) for v in V)
        ((a, True), (b, True), (c, True),
         (r, False), (s, False), (t, False),
         (x, False), (y, False), (z, False))

    ::

        sage: R = SymbolicSubring(rejecting_variables=(r, s, t)); R
        Symbolic Subring rejecting variables r, s, t
        sage: tuple((v, var_in_subring(v, R)) for v in V)
        ((a, True), (b, True), (c, True),
         (r, False), (s, False), (t, False),
         (x, True), (y, True), (z, True))

    ::

        sage: C = SymbolicSubring(only_constants=True); C
        Symbolic Constants Subring
        sage: tuple((v, var_in_subring(v, C)) for v in V)
        ((a, False), (b, False), (c, False),
         (r, False), (s, False), (t, False),
         (x, False), (y, False), (z, False))

    TESTS::

        sage: SymbolicSubring(accepting_variables=tuple()) is C
        True

    ::

        sage: SymbolicSubring(rejecting_variables=tuple()) is SR
        True
    """
    def create_key_and_extra_args(
            self, accepting_variables=None, rejecting_variables=None,
            only_constants=False, **kwds):
        r"""
        Given the arguments and keyword, create a key that uniquely
        determines this object.
        """
        if accepting_variables is None and \
           rejecting_variables is None and \
           only_constants == False:
            raise ValueError('Cannot create a symbolic subring '
                             'since nothing specified.')
        if accepting_variables is not None and \
           rejecting_variables is not None or \
           rejecting_variables is not None and \
           only_constants == True or \
           only_constants == True and \
           accepting_variables is not None:
            raise ValueError('Cannot create a symbolic subring '
                             'since input is ambiguous.')

        if accepting_variables is not None:
            vars = tuple(accepting_variables)
            if vars:
                cls = SymbolicSubringAcceptingVars
            else:
                cls = SymbolicConstantsSubring
        elif rejecting_variables is not None:
            vars = tuple(rejecting_variables)
            cls = SymbolicSubringRejectingVars
        elif only_constants:
            vars = tuple()
            cls = SymbolicConstantsSubring

        vars = tuple(sorted(iter(SR(v) for v in vars), key=str))
        return (cls, vars), kwds


    def create_object(self, version, key, **kwds):
        r"""
        Create an object from the given arguments.
        """
        cls, vars = key
        if cls is SymbolicSubringRejectingVars and not vars:
            return SR
        return cls(vars, **kwds)


SymbolicSubring = SymbolicSubringFactory("SymbolicSubring")


class GenericSymbolicSubring(SymbolicRing):

    def __init__(self, vars):
        super(GenericSymbolicSubring, self).__init__()
        self._vars_ = set(vars)


    def _repr_variables_(self):
        return ', '.join(str(v) for v in sorted(self._vars_, key=str))


    def is_variable_valid(self, var):
        raise NotImplementedError('Not implemented in this abstract base class')


    def _element_constructor_(self, x):
        expression = super(GenericSymbolicSubring, self)._element_constructor_(x)
        if not all(self.is_variable_valid(var)
                   for var in expression.variables()):
            raise ValueError('%s is not contained in %s' % (x, self))


class SymbolicSubringAcceptingVars(GenericSymbolicSubring):

    def _repr_(self):
        return 'Symbolic Subring accepting variables %s' % \
            (self._repr_variables_())


    def is_variable_valid(self, var):
        return var in self._vars_


class SymbolicSubringRejectingVars(GenericSymbolicSubring):

    def _repr_(self):
        return 'Symbolic Subring rejecting variables %s' % \
            (self._repr_variables_())


    def is_variable_valid(self, var):
        return var not in self._vars_


class SymbolicConstantsSubring(GenericSymbolicSubring):

    def _repr_(self):
        return 'Symbolic Constants Subring'


    def is_variable_valid(self, var):
        return False
