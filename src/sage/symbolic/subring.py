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

from ring import SymbolicRing, SR
from expression import Expression


from sage.structure.factory import UniqueFactory
class SymbolicSubringFactory(UniqueFactory):
    r"""
    A factory creating a symbolic subring.

    INPUT:

    - ``accepting_variables`` (default: ``None``) -- a tuple or other
      iterable of variables. If specified, then a symbolic subring of expressions
      in only these variables is created.

    - ``rejecting_variables`` (default: ``None``) -- a tuple or other
      iterable of variables. If specified, then a symbolic subring of expressions
      in variables not distinct of these variables is created.

    - ``only_constants`` (default: ``False``) -- a boolean. If set,
      then a symbolic subring of constant expressions (i.e.,
      expressions without a variable) is created.

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
        Symbolic Subring accepting the variables a, b, c
        sage: tuple((v, var_in_subring(v, A)) for v in V)
        ((a, True), (b, True), (c, True),
         (r, False), (s, False), (t, False),
         (x, False), (y, False), (z, False))

    ::

        sage: R = SymbolicSubring(rejecting_variables=(r, s, t)); R
        Symbolic Subring rejecting the variables r, s, t
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

        See :class:`SymbolicSubringFactory` for details.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring.create_key_and_extra_args()
            Traceback (most recent call last):
            ...
            ValueError: Cannot create a symbolic subring since nothing specified.
            sage: SymbolicSubring.create_key_and_extra_args(
            ....:     accepting_variables=('a',), rejecting_variables=('r',))
            Traceback (most recent call last):
            ...
            ValueError: Cannot create a symbolic subring since input is ambiguous.
            sage: SymbolicSubring.create_key_and_extra_args(
            ....:     accepting_variables=('a',), only_constants=True)
            Traceback (most recent call last):
            ...
            ValueError: Cannot create a symbolic subring since input is ambiguous.
            sage: SymbolicSubring.create_key_and_extra_args(
            ....:     rejecting_variables=('r',), only_constants=True)
            Traceback (most recent call last):
            ...
            ValueError: Cannot create a symbolic subring since input is ambiguous.
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

        See :class:`SymbolicSubringFactory` for details.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(rejecting_variables=tuple()) is SR  # indirect doctest
            True
        """
        cls, vars = key
        if cls is SymbolicSubringRejectingVars and not vars:
            return SR
        return cls(vars, **kwds)


SymbolicSubring = SymbolicSubringFactory("SymbolicSubring")


class GenericSymbolicSubring(SymbolicRing):

    def __init__(self, vars):
        r"""
        An abstract base class for a symbolic subring.

        INPUT:

        - ``vars`` -- a tuple of symbolic variables.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',))  # indirect doctest
            Symbolic Subring accepting the variable a
            sage: SymbolicSubring(rejecting_variables=('r',))  # indirect doctest
            Symbolic Subring rejecting the variable r
            sage: SymbolicSubring(only_constants=True)  # indirect doctest
            Symbolic Constants Subring
            sage: SymbolicSubring(rejecting_variables=tuple())  # indirect doctest
            Symbolic Ring
        """
        super(GenericSymbolicSubring, self).__init__()
        self._vars_ = set(vars)


    def _repr_variables_(self):
        r"""
        Return a representation string of the variables.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=tuple())._repr_variables_()
            'no variable'
            sage: SymbolicSubring(accepting_variables=('a',))._repr_variables_()
            'the variable a'
            sage: SymbolicSubring(accepting_variables=('a', 'b'))._repr_variables_()
            'the variables a, b'
        """
        if not self._vars_:
            s = 'no variable'
        elif len(self._vars_) == 1:
            s = 'the variable '
        else:
            s = 'the variables '
        return s + ', '.join(str(v) for v in sorted(self._vars_, key=str))


    def is_variable_valid(self, variable):
        r"""
        Return whether the given ``variable`` is valid in this subring.

        INPUT:

        - ``variable`` -- a symbolic variable.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.symbolic.subring import GenericSymbolicSubring
            sage: GenericSymbolicSubring(vars=tuple()).is_variable_valid(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: Not implemented in this abstract base class
        """
        raise NotImplementedError('Not implemented in this abstract base class')


    def _element_constructor_(self, x):
        r"""
        Creates the element of this subring specified by the input ``x``.

        INPUT:

        - ``x`` -- an object.

        OUTPUT:

        An element of this symbolic subring.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: S = SymbolicSubring(accepting_variables=('a',))
            sage: S('a')  # indirect doctest
            sage: S('x')  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: x is not contained in Symbolic Subring accepting the variable a
        """
        expression = super(GenericSymbolicSubring, self)._element_constructor_(x)
        if not all(self.is_variable_valid(var)
                   for var in expression.variables()):
            raise ValueError('%s is not contained in %s' % (x, self))

    def _coerce_map_from_(self, P):
        r"""
        Return whether ``P`` coerces into this symbolic subring.

        INPUT:

        - ``P`` -- a parent.

        OUTPUT:

        A boolean or ``None``.

        TESTS::

            sage: from sage.symbolic.subring import GenericSymbolicSubring
            sage: GenericSymbolicSubring(vars=tuple()).has_coerce_map_from(SR)  # indirect doctest
            False
        """
        pass


from sage.categories.pushout import ConstructionFunctor
class GenericSymbolicSubringFunctor(ConstructionFunctor):
    r"""
    A base class for the functors constructing symbolic subrings.

    INPUT:

    - ``vars`` -- a tuple, set, or other iterable of symbolic variables.

    EXAMPLES::

        sage: from sage.symbolic.subring import SymbolicSubring
        sage: SymbolicSubring(only_constants=True).construction()[0]  # indirect doctest
        Subring<accepting no variable>

    .. SEEALSO::

        :class:`sage.categories.pushout.ConstructionFunctor`.
    """

    _functor_name = 'GenericSymbolicSubringFunctor'

    rank = 11

    # The symbolic subring construction returns an object admitting a
    # coercion map into the original, not vice versa.
    coercion_reversed = True

    _repr_type_ = 'generic'


    def __init__(self, vars):
        r"""
        See :class:`GenericSymbolicSubringFunctor` for details.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',)).construction()[0]  # indirect doctest
            Subring<accepting a>
        """
        self.vars = set(vars)
        from sage.categories.rings import Rings
        super(ConstructionFunctor, self).__init__(Rings(), Rings())


    def _repr_variables_(self):
        r"""
        Return a representation string of the variables.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F = SymbolicSubring(accepting_variables=('a',)).construction()[0]
            sage: F._repr_variables_()
            'a'
        """
        return ', '.join(str(v) for v in sorted(self.vars, key=str))


    def _repr_(self):
        r"""
        Return a representation string of this functor.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',))  # indirect doctest
            Symbolic Subring accepting the variable a
            sage: SymbolicSubring(rejecting_variables=('r',))  # indirect doctest
            Symbolic Subring rejecting the variable r
            sage: SymbolicSubring(only_constants=True)  # indirect doctest
            Symbolic Constants Subring
        """
        return 'Subring<%s%s%s>' % (
            self._repr_type_, ' ' if self._repr_type_ else '',
            self._repr_variables_() if self.vars else 'no variable')


    def merge(self, other):
        r"""
        Merge this functor with ``other`` if possible.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A functor or ``None``.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F = SymbolicSubring(accepting_variables=('a',)).construction()[0]
            sage: F.merge(F) is F
            True
        """
        if self == other:
            return self


    def __eq__(self, other):
        r"""
        Return if this functor is equal to ``other``.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F = SymbolicSubring(accepting_variables=('a',)).construction()[0]
            sage: F == F
            True
        """
        return type(self) == type(other) and self.vars == other.vars


    def __ne__(self, other):
        r"""
        Return if this functor is not equal to ``other``.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F = SymbolicSubring(accepting_variables=('a',)).construction()[0]
            sage: F != F
            False
        """
        return not self == other


class SymbolicSubringAcceptingVars(GenericSymbolicSubring):
    r"""
    The symbolic subring consisting of symbolic expressions in the given variables.
    """

    def _repr_(self):
        r"""
        Return a representation string of this symbolic subring.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',))  # indirect doctest
            Symbolic Subring accepting the variable a
        """
        return 'Symbolic Subring accepting %s' % \
            (self._repr_variables_())


    def is_variable_valid(self, variable):
        r"""
        Return whether the given ``variable`` is valid in this subring.

        INPUT:

        - ``variable`` -- a symbolic variable.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: S = SymbolicSubring(accepting_variables=('a',))
            sage: S.is_variable_valid('a')
            True
            sage: S.is_variable_valid('r')
            False
            sage: S.is_variable_valid('x')
            False
        """
        return SR(variable) in self._vars_


    def construction(self):
        r"""
        Return the functorial construction of this symbolic subring.

        OUTPUT:

        A tuple whose first entry is a construction functor and its second
        is the symbolic ring.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',)).construction()
            (Subring<accepting a>, Symbolic Ring)
        """
        return (SymbolicSubringAcceptingVarsFunctor(self._vars_), SR)
class SymbolicSubringRejectingVars(GenericSymbolicSubring):

    def _repr_(self):
        return 'Symbolic Subring rejecting variables %s' % \
            (self._repr_variables_())




class SymbolicConstantsSubring(GenericSymbolicSubring):
        return SR(variable) not in self._vars_


    def construction(self):
        r"""
        Return the functorial construction of this symbolic subring.

        OUTPUT:

        A tuple whose first entry is a construction functor and its second
        is the symbolic ring.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(rejecting_variables=('r',)).construction()
            (Subring<rejecting r>, Symbolic Ring)
        """
        return (SymbolicSubringRejectingVarsFunctor(self._vars_), SR)

    def _repr_(self):
        return 'Symbolic Constants Subring'


    def is_variable_valid(self, var):
        return False
