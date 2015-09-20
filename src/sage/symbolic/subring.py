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
            a
            sage: S('x')  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: x is not contained in Symbolic Subring accepting the variable a
        """
        expression = SR(x)
        if not all(self.is_variable_valid(var)
                   for var in expression.variables()):
            raise ValueError('%s is not contained in %s' % (x, self))
        return expression


    def _coerce_map_from_(self, P):
        r"""
        Return whether ``P`` coerces into this symbolic subring.

        INPUT:

        - ``P`` -- a parent.

        OUTPUT:

        A boolean or ``None``.

        TESTS::

            sage: from sage.symbolic.subring import GenericSymbolicSubring
            sage: GenericSymbolicSubring(vars=tuple()).has_coerce_map_from(SR)  # indirect doctest  # not tested see #19231
            False

        ::
            sage: from sage.symbolic.subring import SymbolicSubring
            sage: C = SymbolicSubring(only_constants=True)
            sage: C.has_coerce_map_from(ZZ)  # indirect doctest
            True
            sage: C.has_coerce_map_from(QQ)  # indirect doctest
            True
            sage: C.has_coerce_map_from(RR)  # indirect doctest
            True
            sage: C.has_coerce_map_from(RIF)  # indirect doctest
            True
            sage: C.has_coerce_map_from(CC)  # indirect doctest
            True
            sage: C.has_coerce_map_from(CIF)  # indirect doctest
            True
            sage: C.has_coerce_map_from(AA)  # indirect doctest
            True
            sage: C.has_coerce_map_from(QQbar)  # indirect doctest
            True
            sage: C.has_coerce_map_from(SR)  # indirect doctest  # not tested see #19231
            False
        """
        from sage.rings.real_mpfr import mpfr_prec_min
        from sage.rings.all import (ComplexField,
                                    RLF, CLF, AA, QQbar, InfinityRing)
        from sage.rings.real_mpfi import is_RealIntervalField
        from sage.rings.complex_interval_field import is_ComplexIntervalField

        if isinstance(P, type):
            return SR._coerce_map_from_(P)

        elif RLF.has_coerce_map_from(P) or \
             CLF.has_coerce_map_from(P) or \
             AA.has_coerce_map_from(P) or \
             QQbar.has_coerce_map_from(P):
            return True

        elif (P is InfinityRing
              or is_RealIntervalField(P) or is_ComplexIntervalField(P)):
            return True

        elif ComplexField(mpfr_prec_min()).has_coerce_map_from(P):
            return P not in (RLF, CLF, AA, QQbar)


    def __cmp__(self, other):
        """
        Compare two symbolic subrings.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: A = SymbolicSubring(accepting_variables=('a',))
            sage: B = SymbolicSubring(accepting_variables=('b',))
            sage: AB = SymbolicSubring(accepting_variables=('a', 'b'))
            sage: A == A
            True
            sage: A == B
            False
            sage: A == AB
            False
        """
        c = cmp(type(self), type(other))
        if c != 0:
            return c
        if self._vars_ == other._vars_:
            return 0
        return 1


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


    def _coerce_map_from_(self, P):
        r"""
        Return whether ``P`` coerces into this symbolic subring.

        INPUT:

        - ``P`` -- a parent.

        OUTPUT:

        A boolean or ``None``.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: A = SymbolicSubring(accepting_variables=('a',))
            sage: AB = SymbolicSubring(accepting_variables=('a', 'b'))
            sage: A.has_coerce_map_from(AB)  # indirect doctest
            False
            sage: AB.has_coerce_map_from(A)  # indirect doctest
            True
        """
        if isinstance(P, SymbolicSubringAcceptingVars):
            return self._vars_ >= P._vars_
        return super(SymbolicSubringAcceptingVars, self)._coerce_map_from_(P)


    def _an_element_(self):
        r"""
        Return an element of this symbolic subring.

        OUTPUT:

        A symbolic expression.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(accepting_variables=('a',)).an_element()
            a
        """
        return sorted(self._vars_, key=str)[0]


class SymbolicSubringAcceptingVarsFunctor(GenericSymbolicSubringFunctor):

    _functor_name = 'SymbolicSubringAcceptingVarsFunctor'

    _repr_type_ = 'accepting'


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
            sage: G = SymbolicSubring(rejecting_variables=('r',)).construction()[0]
            sage: F.merge(F) is F
            True
            sage: F.merge(G) is G
            True
        """
        if self == other:
            return self
        elif type(self) == type(other):
            return type(self)(self.vars | other.vars)
        elif isinstance(other, SymbolicSubringRejectingVarsFunctor):
            if not (self.vars & other.vars):
                return other


    def _apply_functor(self, R):
        """
        Apply this functor to the given symbolic ring `R`.

        INPUT:

        - ``R`` -- a symbolic ring.

        OUTPUT:

        A subring of ``R``.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F, R = SymbolicSubring(accepting_variables=('a',)).construction()
            sage: F(R)  # indirect doctest
            Symbolic Subring accepting the variable a

        TESTS::

            sage: F(F(R))
            Traceback (most recent call last):
            ...
            NotImplementedError: This functor can only be applied on the
            symbolic ring but Symbolic Subring accepting the variable a given.
        """
        if R is not SR:
            raise NotImplementedError('This functor can only be applied on '
                                      'the symbolic ring but %s given.' % (R,))
        return SymbolicSubring(accepting_variables=self.vars)


class SymbolicSubringRejectingVars(GenericSymbolicSubring):
    r"""
    The symbolic subring consisting of symbolic expressions whose variables
    are none of the given variables.
    """

    def _repr_(self):
        r"""
        Return a representation string of this symbolic subring.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(rejecting_variables=('r',))  # indirect doctest
            Symbolic Subring rejecting the variable r
        """
        return 'Symbolic Subring rejecting %s' % \
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
            sage: S = SymbolicSubring(rejecting_variables=('r',))
            sage: S.is_variable_valid('a')
            True
            sage: S.is_variable_valid('r')
            False
            sage: S.is_variable_valid('x')
            True
        """
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


    def _coerce_map_from_(self, P):
        r"""
        Return whether ``P`` coerces into this symbolic subring.

        INPUT:

        - ``P`` -- a parent.

        OUTPUT:

        A boolean or ``None``.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: R = SymbolicSubring(rejecting_variables=('r',))
            sage: RS = SymbolicSubring(rejecting_variables=('r', 's'))
            sage: RS.has_coerce_map_from(R)  # indirect doctest
            False
            sage: R.has_coerce_map_from(RS)  # indirect doctest
            True
        """
        if isinstance(P, SymbolicSubringRejectingVars):
            return self._vars_ <= P._vars_
        return super(SymbolicSubringRejectingVars, self)._coerce_map_from_(P)


    def _an_element_(self):
        r"""
        Return an element of this symbolic subring.

        OUTPUT:

        A symbolic expression.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(rejecting_variables=('r',)).an_element()
            some_variable
            sage: SymbolicSubring(rejecting_variables=('some_variable',)).an_element()
            some_some_variable
            sage: SymbolicSubring(rejecting_variables=('some_some_variable',)).an_element()
            some_variable
            sage: SymbolicSubring(rejecting_variables=('some_variable','some_some_variable')).an_element()
            some_some_some_variable
        """
        v = SR.an_element()
        while not self.is_variable_valid(v):
            v = SR('some_' + str(v))
        return v


class SymbolicSubringRejectingVarsFunctor(GenericSymbolicSubringFunctor):

    _functor_name = 'SymbolicSubringRejectingVarsFunctor'

    _repr_type_ = 'rejecting'


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
            sage: G = SymbolicSubring(rejecting_variables=('r',)).construction()[0]
            sage: G.merge(G) is G
            True
            sage: G.merge(F) is G
            True
        """
        if self == other:
            return self
        elif type(self) == type(other):
            return type(self)(self.vars & other.vars)
        elif isinstance(other, SymbolicSubringAcceptingVarsFunctor):
            if not (self.vars & other.vars):
                return self


    def _apply_functor(self, R):
        """
        Apply this functor to the given symbolic ring `R`.

        INPUT:

        - ``R`` -- a symbolic ring.

        OUTPUT:

        A subring of ``R``.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: F, R = SymbolicSubring(rejecting_variables=('r',)).construction()
            sage: F(R)  # indirect doctest
            Symbolic Subring rejecting the variable r

        TESTS::

            sage: F(F(R))
            Traceback (most recent call last):
            ...
            NotImplementedError: This functor can only be applied on the
            symbolic ring but Symbolic Subring rejecting the variable r given.
        """
        if R is not SR:
            raise NotImplementedError('This functor can only be applied on '
                                      'the symbolic ring but %s given.' % (R,))
        return SymbolicSubring(rejecting_variables=self.vars)


class SymbolicConstantsSubring(SymbolicSubringAcceptingVars):
    r"""
    The symbolic subring consisting of symbolic constants.
    """

    def _repr_(self):
        r"""
        Return a representation string of this symbolic subring.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(only_constants=True)  # indirect doctest
            Symbolic Constants Subring
        """
        return 'Symbolic Constants Subring'


    def is_variable_valid(self, variable):
        r"""
        Return whether the given ``variable`` is valid in this subring.

        INPUT:

        - ``variable`` -- a symbolic variable.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: S = SymbolicSubring(only_constants=True)
            sage: S.is_variable_valid('a')
            False
            sage: S.is_variable_valid('r')
            False
            sage: S.is_variable_valid('x')
            False
        """
        return False


    def _an_element_(self):
        r"""
        Return an element of this symbolic subring.

        OUTPUT:

        A symbolic expression.

        TESTS::

            sage: from sage.symbolic.subring import SymbolicSubring
            sage: SymbolicSubring(only_constants=True).an_element()
            I*pi*e
        """
        return SR('I') * SR('pi') * SR('e')
