r"""
(Asymptotic) Growth Groups

This module adds support for (asymptotic) growth groups. Such groups
are equipped with a partial order: the elements can be seen as
functions, and their behavior as the argument(s) get large (tend to
`\infty`) is compared.

Besides an abstract base class :class:`GenericGrowthGroup`, this module
contains concrete realizations of growth groups. At the moment there
is

- :class:`MonomialGrowthGroup` (whose elements are powers of a fixed symbol).

More complex growth groups can be constructed via cartesian products.

These growth groups are used behind the scenes when performing
calculations in an asymptotic ring (to be implemented).

AUTHORS:

- Benjamin Hackl (2015-01): initial version
- Daniel Krenn (2015-05-29): initial version and review
- Daniel Krenn (2015-06-02): cartesian products
- Benjamin Hackl (2015-07): growth group factory
- Benjamin Hackl (2015-08): exponential growth group, initial version

.. WARNING::

    As this code is experimental, warnings are thrown when a growth
    group is created for the first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup, GrowthGroup
        sage: GenericGrowthGroup(ZZ)
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        Growth Group Generic(ZZ)
        sage: GrowthGroup('x^ZZ * log(x)^ZZ')
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        Growth Group x^ZZ * log(x)^ZZ

.. NOTE::

    By using the following short notation for growth groups, their
    creation is very simple: *Monomial growth groups* (i.e. the
    group for powers of a fixed symbol;
    :class:`~sage.rings.asymptotic.growth_group.MonomialGrowthGroup`)
    are denoted as ``variable^base``, e.g. ``x^ZZ`` and ``y^QQ`` for
    the group of integer powers of `x`, and the group of rational
    powers of `y`, respectively.

    This also enables us to construct *logarithmic growth groups*,
    e.g. ``log(x)^ZZ``.

    Exponential growth groups, i.e. growth groups representing
    elements of the form `\operatorname{base}^\operatorname{variable}`
    are denoted as ``base^variable``. For example, ``QQ^x`` denotes
    the multiplicative group of exponential expressions `q^x`, where
    `q \in \mathbb{Q}^{\times}`.

EXAMPLES::

    sage: import sage.rings.asymptotic.growth_group as agg
    sage: G_x = agg.GrowthGroup('x^ZZ'); repr(G_x)
    'Growth Group x^ZZ'
    sage: G_xy = agg.GrowthGroup('x^ZZ * y^ZZ'); G_xy
    Growth Group x^ZZ * y^ZZ
    sage: G_xy.an_element()
    x * y
    sage: x = G_xy('x'); y = G_xy('y')
    sage: elem = x^21 * y^21; elem^2
    x^42 * y^42

A monomial growth group itself is totally ordered, all elements
are comparable. However, this does **not** hold for cartesian
products::

    sage: e1 = x^2 * y; e2 = x * y^2
    sage: e1 <= e2 or e2 <= e1
    False

In terms of uniqueness, we have the following behaviour::

    sage: agg.GrowthGroup('x^ZZ * y^ZZ') is agg.GrowthGroup('y^ZZ * x^ZZ')
    True

The above is ``True`` since the order of the factors does not play a role here; they use different variables. But when using the same variable, it plays a role::

    sage: agg.GrowthGroup('x^ZZ * log(x)^ZZ') is agg.GrowthGroup('log(x)^ZZ * x^ZZ')
    False

(Note that it is mathematically nonsense to make ``log(x)`` larger than ``x``.)

With the help of the short notation, even complicated growth groups
can be constructed easily::

    sage: G = agg.GrowthGroup('QQ^x * x^ZZ * log(x)^QQ * y^QQ')
    sage: G.an_element()
    (1/2)^x * x * log(x)^(1/2) * y^(1/2)
    sage: (x, y) = var('x y')
    sage: G(2^x * log(x) * y^(1/2)) * G(x^(-5) * 5^x * y^(1/3))
    10^x * x^(-5) * log(x) * y^(5/6)
"""

#*****************************************************************************
# Copyright (C) 2014--2015 Benjamin Hackl <benjamin.hackl@aau.at>
#               2014--2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage
from sage.misc.lazy_import import lazy_import
lazy_import('sage.rings.asymptotic.growth_group_cartesian', 'CartesianProductGrowthGroups')


def repr_short_to_parent(s):
    r"""
    Helper method for the growth group factory, which converts a short
    representation string to a parent.

    INPUT:

    A string.

    OUTPUT:

    A parent.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: agg.repr_short_to_parent('ZZ')
        Integer Ring
        sage: agg.repr_short_to_parent('QQ')
        Rational Field
        sage: agg.repr_short_to_parent('SR')
        Symbolic Ring

    TESTS::

        sage: agg.repr_short_to_parent('abcdef')
        Traceback (most recent call last):
        ...
        ValueError: Cannot create a parent out of 'abcdef'.
    """
    if s == 'ZZ':
        return sage.rings.integer_ring.ZZ
    elif s == 'QQ':
        return sage.rings.rational_field.QQ
    elif s == 'SR':
        return sage.symbolic.ring.SR
    else:
        raise ValueError("Cannot create a parent out of '%s'." % s)


def parent_to_repr_short(P):
    r"""
    Helper method which generates a short(er) representation string
    out of a parent.

    INPUT:

    A parent.

    OUTPUT:

    A string.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: agg.parent_to_repr_short(ZZ)
        'ZZ'
        sage: agg.parent_to_repr_short(QQ)
        'QQ'
        sage: agg.parent_to_repr_short(SR)
        'SR'
        sage: agg.parent_to_repr_short(ZZ[x])
        '(Univariate Polynomial Ring in x over Integer Ring)'
    """
    if P is sage.rings.integer_ring.ZZ:
        return 'ZZ'
    elif P is sage.rings.rational_field.QQ:
        return 'QQ'
    elif P is sage.symbolic.ring.SR:
        return 'SR'
    else:
        rep = repr(P)
        if ' ' in rep:
            rep = '(' + rep + ')'
        return rep


class Variable(sage.structure.unique_representation.CachedRepresentation,
               sage.structure.sage_object.SageObject):
    r"""
    A class managing the variable of a growth group.

    INPUT:

    - ``var`` -- an object whose representation string is used as the
      variable. It has to be a valid Python identifier. ``var`` can
      also be a tuple (or other iterable of such objects).

    - ``repr`` -- (default: ``None``) if specified, then this string
      will be displayed instead of ``var``. Use this to get
      e.g. ``log(x)^ZZ``: ``var`` is then used to specify the variable `x`.


    TESTS::

        sage: from sage.rings.asymptotic.growth_group import Variable
        sage: v = Variable('x'); repr(v), v.variable_names()
        ('x', ('x',))
        sage: v = Variable('x1'); repr(v), v.variable_names()
        ('x1', ('x1',))
        sage: v = Variable('x_42'); repr(v), v.variable_names()
        ('x_42', ('x_42',))
        sage: v = Variable(' x'); repr(v), v.variable_names()
        ('x', ('x',))
        sage: v = Variable('x '); repr(v), v.variable_names()
        ('x', ('x',))
        sage: v = Variable('x '); repr(v), v.variable_names()
        ('x', ('x',))
        sage: v = Variable(''); repr(v), v.variable_names()
        ('', ())

    ::

        sage: v = Variable(('x', 'y')); repr(v), v.variable_names()
        ('x, y', ('x', 'y'))
        sage: v = Variable(('x', 'log(y)')); repr(v), v.variable_names()
        ('x, log(y)', ('x', 'y'))
        sage: v = Variable(('x', 'log(x)')); repr(v), v.variable_names()
        Traceback (most recent call last):
        ...
        ValueError: Variable names ('x', 'x') are not pairwise distinct.

    ::

        sage: v = Variable('log(x)'); repr(v), v.variable_names()
        ('log(x)', ('x',))
        sage: v = Variable('log(log(x))'); repr(v), v.variable_names()
        ('log(log(x))', ('x',))

    ::

        sage: v = Variable('x', repr='log(x)'); repr(v), v.variable_names()
        ('log(x)', ('x',))
    """
    def __init__(self, var, repr=None):
        r"""
        See :class:`Variable` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('blub')
            blub
            sage: Variable('blub') is Variable('blub')
            True
        """
        from sage.symbolic.ring import isidentifier

        if not isinstance(var, (list, tuple)):
            var = (var,)
        var = tuple(str(v).strip() for v in var)

        if repr is None:
            var_bases = sum(iter(
                self.extract_variable_names(v)
                if not isidentifier(v) else (v,)
                for v in var), tuple())
            var_repr = ', '.join(var)
        else:
            for v in var:
                if not isidentifier(v):
                    raise ValueError("'%s' is not a valid name for a variable." % (v,))
            var_bases = var
            var_repr = str(repr).strip()

        if len(var_bases) != len(set(var_bases)):
            raise ValueError('Variable names %s are not pairwise distinct.' %
                             (var_bases,))
        self.var_bases = var_bases
        self.var_repr = var_repr


    def __hash__(self):
        r"""
        Return the hash if this variable.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: hash(Variable('blub'))  # random
            -123456789
        """
        return hash((self.var_repr,) + self.var_bases)


    def __eq__(self, other):
        r"""
        Compares if this variable equals ``other``.

        INPUT:

        - ``other`` -- another variable.

        OUTPUT:

        A boolean.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x') == Variable('x')
            True
            sage: Variable('x') == Variable('y')
            False
        """
        return self.var_repr == other.var_repr and self.var_bases == other.var_bases


    def __ne__(self, other):
        r"""
        Compares if this variable does not equal ``other``.

        INPUT:

        - ``other`` -- another variable.

        OUTPUT:

        A boolean.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x') != Variable('x')
            False
            sage: Variable('x') != Variable('y')
            True
        """
        return not self.__eq__(other)


    def _repr_(self):
        r"""
        Return a representation string of this variable.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('blub')  # indirect doctest
            blub
        """
        return self.var_repr


    def variable_names(self):
        r"""
        Return the names of the variables.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x').variable_names()
            ('x',)
            sage: Variable('log(x)').variable_names()
            ('x',)
        """
        return self.var_bases


    def is_monomial(self):
        r"""
        Returns if this is a monomial variable.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable('x').is_monomial()
            True
            sage: Variable('log(x)').is_monomial()
            False
        """
        return len(self.var_bases) == 1 and self.var_bases[0] == self.var_repr


    @staticmethod
    def extract_variable_names(s):
        r"""
        Finds the name of the variable for the given string.

        INPUT:

        - ``s`` -- a string.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import Variable
            sage: Variable.extract_variable_names('')
            ()
            sage: Variable.extract_variable_names('x')
            ('x',)
            sage: Variable.extract_variable_names('exp(x)')
            ('x',)
            sage: Variable.extract_variable_names('sin(cos(ln(x)))')
            ('x',)

        ::

            sage: Variable.extract_variable_names('log(77w)')
            Traceback (most recent call last):
            ....
            ValueError: '77w' is not a valid name for a variable.
            sage: Variable.extract_variable_names('log(x')
            Traceback (most recent call last):
            ....
            ValueError: Unbalanced parentheses in 'log(x'.
            sage: Variable.extract_variable_names('x)')
            Traceback (most recent call last):
            ....
            ValueError: Unbalanced parentheses in 'x)'.
            sage: Variable.extract_variable_names('log)x(')
            Traceback (most recent call last):
            ....
            ValueError: Unbalanced parentheses in 'log)x('.
            sage: Variable.extract_variable_names('log(x)+y')
            ('x', 'y')

        ::

            sage: Variable.extract_variable_names('a + b')
            ('a', 'b')
            sage: Variable.extract_variable_names('a+b')
            ('a', 'b')
            sage: Variable.extract_variable_names('a +b')
            ('a', 'b')
            sage: Variable.extract_variable_names('+a')
            ('a',)
            sage: Variable.extract_variable_names('a+')
            ('a',)
            sage: Variable.extract_variable_names('b!')
            ('b',)
            sage: Variable.extract_variable_names('-a')
            ('a',)
            sage: Variable.extract_variable_names('a*b')
            ('a', 'b')
            sage: Variable.extract_variable_names('2^q')
            ('q',)
            sage: Variable.extract_variable_names('77')
            ()
        """
        from sage.symbolic.ring import isidentifier
        import re
        numbers = re.compile(r"\d+$")
        vars = []

        def find_next_outer_parentheses(s):
            op = s.find('(')
            level = 1
            for i, c in enumerate(s[op+1:]):
                if c == ')':
                    level -= 1
                if c == '(':
                    level += 1
                if level == 0:
                    return op, op+i+1
            return op, -1

        def strip(s):
            s = s.strip()
            if not s:
                return

            # parentheses (...)
            # functions f(...)
            op, cl = find_next_outer_parentheses(s)
            if (op == -1) != (cl == -1) or op > cl:
                raise ValueError("Unbalanced parentheses in '%s'." % (s,))
            if cl != -1:
                strip(s[op+1:cl])
                strip(s[cl+1:])
                return

            # unary +a, a+, ...
            # binary a+b, a*b, ...
            for operator in ('**', '+', '-', '*', '/', '^', '!'):
                a, o, b = s.partition(operator)
                if o:
                    strip(a)
                    strip(b)
                    return

            # a number
            if numbers.match(s) is not None:
                return

            # else: a variable
            if not isidentifier(s):
                raise ValueError("'%s' is not a valid name for a variable." % (s,))
            vars.append(s)

        strip(s)
        return tuple(vars)


class GenericGrowthElement(sage.structure.element.MultiplicativeGroupElement):
    r"""
    An abstract implementation of a generic growth element.

    Growth elements form a group by multiplication, and (some of) the
    elements can be compared to each other, i.e., all elements form a
    poset.

    INPUT:

    - ``parent`` -- a :class:`GenericGrowthGroup`.

    - ``raw_element`` -- an element from the base of the parent.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G = agg.GenericGrowthGroup(ZZ)
        sage: g = agg.GenericGrowthElement(G, 42); g
        GenericGrowthElement(42)
        sage: g.parent()
        Growth Group Generic(ZZ)
        sage: G(raw_element=42) == g
        True
    """

    def __init__(self, parent, raw_element):
        r"""
        See :class:`GenericGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G(raw_element=42)
            GenericGrowthElement(42)

        TESTS::

            sage: G(raw_element=42).category()
            Category of elements of Growth Group Generic(ZZ)

        ::

            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G(raw_element=42).category()
            Category of elements of Growth Group Generic(ZZ)

        ::

            sage: agg.GenericGrowthElement(None, 0)
            Traceback (most recent call last):
            ...
            ValueError: The parent must be provided
        """
        if parent is None:
            raise ValueError('The parent must be provided')
        super(GenericGrowthElement, self).__init__(parent=parent)

        self._raw_element_ = parent.base()(raw_element)


    def _repr_(self):
        r"""
        A representation string for this abstract generic element.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: G = GenericGrowthGroup(ZZ)
            sage: G(raw_element=42)  # indirect doctest
            GenericGrowthElement(42)
            sage: H = GenericGrowthGroup(ZZ, 'h')
            sage: H(raw_element=42)  # indirect doctest
            GenericGrowthElement(42, h)
        """
        vars = ', '.join(self.parent()._var_.variable_names())
        if vars:
            vars = ', ' + vars
        return 'GenericGrowthElement(%s%s)' % (self._raw_element_, vars)


    def __hash__(self):
        r"""
        Return the hash of this element.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ);
            sage: hash(G(raw_element=42))  # random
            5656565656565656
        """
        return hash((self.parent(), self._raw_element_))


    def _mul_(self, other):
        r"""
        Abstract multiplication method for generic elements.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`.

        OUTPUT:

        A class:`GenericGrowthElement` representing the product with
        ``other``.

        .. NOTE::

            Inherited classes must override this.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: g = G.an_element()
            sage: g * g
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
        """
        raise NotImplementedError('Only implemented in concrete realizations.')


    def _div_(self, other):
        r"""
        Divide this growth element by another one.

        INPUT:

        - ``other`` -- an instance of :class:`GenericGrowthElement`.

        OUTPUT:

        An instance of :class:`GenericGrowthElement`.

        .. NOTE::

            This method is called by the coercion framework, thus, it can be
            assumed that this element is of the same type as ``other``.
            The output will be of the same type as well.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: e1 = P(raw_element=2)
            sage: e2 = e1._div_(P.gen()); e2
            x
            sage: e2 == e1 / P.gen()
            True
        """
        return self._mul_(~other)


    def __pow__(self, power):
        r"""
        Raises this growth element to the given ``power``.

        INPUT:

        - ``power`` -- a number. This can be anything that is a
          valid right hand side of ``*`` with elements of the
          parent's base.

        OUTPUT:

        The result of this exponentiation.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G.an_element()^7
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
        """
        raise NotImplementedError('Only implemented in concrete realizations.')


    def __eq__(self, other):
        r"""
        Return if this growth element is equal to ``other``.

        INPUT:

        - ``other`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

            The comparison of two elements with the same parent is done in
            :meth:`_eq_`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G.an_element() == G.an_element()
            True
            sage: G(raw_element=42) == G(raw_element=7)
            False

        ::

            sage: G_ZZ = agg.GenericGrowthGroup(ZZ)
            sage: G_QQ = agg.GenericGrowthGroup(QQ)
            sage: G_ZZ(raw_element=1) == G_QQ(raw_element=1)
            True

        ::

            sage: P_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P_ZZ.gen() == P_QQ.gen()
            True
            sage: ~P_ZZ.gen() == P_ZZ.gen()
            False
            sage: ~P_ZZ(1) == P_ZZ(1)
            True
        """
        from sage.structure.element import have_same_parent
        if have_same_parent(self, other):
            return self._eq_(other)

        from sage.structure.element import get_coercion_model
        import operator
        try:
            return get_coercion_model().bin_op(self, other, operator.eq)
        except TypeError:
            return False


    def _eq_(self, other):
        r"""
        Return if this :class:`GenericGrowthElement` is equal to ``other``.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function compares two instances of
            :class:`GenericGrowthElement`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: e1 = P(raw_element=1)
            sage: e1._eq_(P.gen())
            True
            sage: e2 = e1^4
            sage: e2 == e1^2 * e1 * e1
            True
            sage: e2 == e1
            False
        """
        return self._raw_element_ == other._raw_element_


    def __le__(self, other):
        r"""
        Return if this growth element is at most (less than or equal
        to) ``other``.

        INPUT:

        - ``other`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

            The comparison of two elements with the same parent is done in
            :meth:`_le_`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P_ZZ.gen() <= P_QQ.gen()^2
            True
            sage: ~P_ZZ.gen() <= P_ZZ.gen()
            True
        """
        from sage.structure.element import have_same_parent
        if have_same_parent(self, other):
            return self._le_(other)

        from sage.structure.element import get_coercion_model
        import operator
        try:
            return get_coercion_model().bin_op(self, other, operator.le)
        except TypeError:
            return False


    def _le_(self, other):
        r"""
        Return if this :class:`GenericGrowthElement` is at most (less
        than or equal to) ``other``.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function compares two instances of
            :class:`GenericGrowthElement`.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: e1 = G(raw_element=1); e2 = G(raw_element=2)
            sage: e1 <= e2  # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
        """
        raise NotImplementedError('Only implemented in concrete realizations.')


class GenericGrowthGroup(
        sage.structure.unique_representation.UniqueRepresentation,
        sage.structure.parent.Parent):
    r"""
    An abstract implementation for growth groups.

    INPUT:

    - ``base`` -- one of SageMath's parents, out of which the elements
      get their data (``raw_element``).

    - ``category`` -- (default: ``None``) the category of the newly
      created growth group. It has to be a subcategory of ``Join of
      Category of groups and Category of posets``. This is also the
      default category if ``None`` is specified.

    .. NOTE::

        This class should be derived to get concrete implementations.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G = agg.GenericGrowthGroup(ZZ); G
        Growth Group Generic(ZZ)

    .. SEEALSO::

        :class:`MonomialGrowthGroup`
    """
    # TODO: implement some sort of 'assume', where basic assumptions
    # for the variables can be stored. --> within the cartesian product

    # enable the category framework for elements
    Element = GenericGrowthElement


    @staticmethod
    def __classcall__(cls, base, var=None, category=None):
        r"""
        Normalizes the input in order to ensure a unique
        representation.

        For more information see :class:`MonomialGrowthGroup`.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P1 = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P2 = agg.MonomialGrowthGroup(ZZ, ZZ['x'].gen())
            sage: P3 = agg.MonomialGrowthGroup(ZZ, SR.var('x'))
            sage: P1 is P2 and P2 is P3
            True
            sage: P4 = agg.MonomialGrowthGroup(ZZ, buffer('xylophone', 0, 1))
            sage: P1 is P4
            True
            sage: P5 = agg.MonomialGrowthGroup(ZZ, 'x ')
            sage: P1 is P5
            True

        ::

            sage: L1 = agg.MonomialGrowthGroup(QQ, log(x))
            sage: L2 = agg.MonomialGrowthGroup(QQ, 'log(x)')
            sage: L1 is L2
            True
        """
        if var is None:
            var = Variable('')
        elif not isinstance(var, Variable):
            var = Variable(var)
        return super(GenericGrowthGroup, cls).__classcall__(
            cls, base, var, category)


    @sage.misc.superseded.experimental(trac_number=17601)
    def __init__(self, base, var, category=None):
        r"""
        See :class:`GenericGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GenericGrowthGroup(ZZ).category()
            Join of Category of groups and Category of posets

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'x')
            Growth Group x^ZZ
            sage: agg.MonomialGrowthGroup(QQ, SR.var('n'))
            Growth Group n^QQ
            sage: agg.MonomialGrowthGroup(ZZ, ZZ['y'].gen())
            Growth Group y^ZZ
            sage: agg.MonomialGrowthGroup(QQ, 'log(x)')
            Growth Group log(x)^QQ

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.ExponentialGrowthGroup(QQ, 'x')
            Growth Group QQ^x
            sage: agg.ExponentialGrowthGroup(SR, ZZ['y'].gen())
            Growth Group SR^y

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G.is_parent_of(G(raw_element=42))
            True
            sage: G2 = agg.GenericGrowthGroup(ZZ, category=FiniteGroups() & Posets())
            sage: G2.category()
            Join of Category of finite groups and Category of finite posets
            sage: G3 = agg.GenericGrowthGroup(ZZ, category=Rings())
            Traceback (most recent call last):
            ...
            ValueError: (Category of rings,) is not a subcategory of Join of Category of groups and Category of posets

        ::

            sage: G = agg.GenericGrowthGroup('42')
            Traceback (most recent call last):
            ...
            TypeError: 42 is not a valid base

        ::

            sage: agg.MonomialGrowthGroup('x', ZZ)
            Traceback (most recent call last):
            ...
            ValueError: 'Integer Ring' is not a valid name for a variable.
            sage: agg.MonomialGrowthGroup('x', 'y')
            Traceback (most recent call last):
            ...
            TypeError: x is not a valid base

        ::

            sage: agg.ExponentialGrowthGroup('x', ZZ)
            Traceback (most recent call last):
            ...
            ValueError: 'Integer Ring' is not a valid name for a variable.
            sage: agg.ExponentialGrowthGroup('x', 'y')
            Traceback (most recent call last):
            ...
            TypeError: x is not a valid base

        """
        if not isinstance(base, sage.structure.parent.Parent):
            raise TypeError('%s is not a valid base' % (base,))
        from sage.categories.groups import Groups
        from sage.categories.posets import Posets

        if category is None:
            category = Groups() & Posets()
        else:
            if not isinstance(category, tuple):
                category = (category,)
            if not any(cat.is_subcategory(Groups() & Posets()) for cat in
                       category):
                raise ValueError('%s is not a subcategory of %s'
                                 % (category, Groups() & Posets()))

        self._var_ = var
        super(GenericGrowthGroup, self).__init__(category=category,
                                                 base=base)


    def _repr_short_(self):
        r"""
        A short representation string of this abstract growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: GenericGrowthGroup(QQ)._repr_short_()
            'Generic(QQ)'
            sage: GenericGrowthGroup(QQ)
            Growth Group Generic(QQ)
            sage: GenericGrowthGroup(QQ, ('a', 'b'))
            Growth Group Generic(QQ, a, b)
        """
        vars = ', '.join(self._var_.variable_names())
        if vars:
            vars = ', ' + vars
        return 'Generic(%s%s)' % (parent_to_repr_short(self.base()), vars)


    def _repr_(self, condense=False):
        r"""
        A representations string of this growth group.

        INPUT:

        - ``condense`` -- (default: ``False``) if set, then a shorter
          output is returned, e.g. the prefix-string ``Growth Group``
          is not show in this case.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'x')  # indirect doctest
            Growth Group x^ZZ
            sage: agg.MonomialGrowthGroup(QQ, 'log(x)')  # indirect doctest
            Growth Group log(x)^QQ

        TESTS::

            sage: agg.MonomialGrowthGroup(QQ, 'log(x)')._repr_(condense=True)
            'log(x)^QQ'
        """
        pre = 'Growth Group ' if not condense else ''
        return '%s%s' % (pre, self._repr_short_())


    def __hash__(self):
        r"""
        Return the hash of this group.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: hash(agg.GenericGrowthGroup(ZZ))  # random
            4242424242424242

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: hash(P)  # random
            -1234567890123456789

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.ExponentialGrowthGroup(ZZ, 'x')
            sage: hash(P)  # random
            -1234567890123456789
        """
        return hash((self.__class__, self.base(), self._var_))


    def _an_element_(self):
        r"""
        Return an element of ``self``.

        INPUT:

        Nothing.

        OUTPUT:

        An element of ``self``.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ);
            sage: G.an_element()  # indirect doctest
            GenericGrowthElement(1)
        """
        return self.element_class(self, self.base().an_element())


    def le(self, left, right):
        r"""
        Return if the growth of ``left`` is at most (less than or
        equal to) the growth of ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: x = G.gen()
            sage: G.le(x, x^2)
            True
            sage: G.le(x^2, x)
            False
            sage: G.le(x^0, 1)
            True
        """
        return self(left) <= self(right)


    def _element_constructor_(self, data, raw_element=None):
        r"""
        Converts a given object to this growth group.

        INPUT:

        - ``data`` -- an object representing the element to be
          initialized.

        - ``raw_element`` -- (default: ``None``) if given, then this is
          directly passed to the element constructor (i.e., no conversion
          is performed).

        OUTPUT:

        An element of this growth group.

        .. NOTE::

            This method calls :meth:`_convert_`, which does the actual
            conversion from ``data``.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_ZZ = agg.GenericGrowthGroup(ZZ)
            sage: z = G_ZZ(raw_element=42); z
            GenericGrowthElement(42)
            sage: z is G_ZZ(z)
            True

        ::

            sage: G_QQ = agg.GenericGrowthGroup(QQ)
            sage: q = G_QQ(raw_element=42)
            sage: q is z
            False
            sage: G_ZZ(q)
            GenericGrowthElement(42)
            sage: G_QQ(z)
            GenericGrowthElement(42)
            sage: q is G_ZZ(q)
            False

        ::

            sage: G_ZZ()
            Traceback (most recent call last):
            ...
            ValueError: No input specified. Cannot continue.
            sage: G_ZZ('blub')
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert blub.
            sage: G_ZZ('x', raw_element=42)
            Traceback (most recent call last):
            ...
            ValueError: Input is ambigous: x as well as raw_element=42 are specified

        ::

            sage: x = agg.MonomialGrowthGroup(ZZ, 'x')(raw_element=1)
            sage: G_y = agg.MonomialGrowthGroup(ZZ, 'y')
            sage: G_y(x)
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert x.
        """
        if raw_element is None:
            if type(data) == self.element_class and data.parent() == self:
                return data
            elif isinstance(data, self.element_class):
                try:
                    if self._var_ != data.parent()._var_:
                        raise ValueError('Cannot convert %s.' % (data,))
                except AttributeError:
                    pass
                raw_element = data._raw_element_
            elif type(data) == int and data == 0:
                raise ValueError('No input specified. Cannot continue.')
            else:
                raw_element = self._convert_(data)
            if raw_element is None:
                raise ValueError('Cannot convert %s.' % (data,))
        elif type(data) != int or data != 0:
            raise ValueError('Input is ambigous: '
                             '%s as well as raw_element=%s '
                             'are specified' % (data, raw_element))

        return self.element_class(self, raw_element)


    def _convert_(self, data):
        r"""
        Converts ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        .. NOTE::

            This method always returns ``None`` in this abstract base
            class, and should be overridden in inherited class.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G._convert_('icecream') is None
            True
        """
        pass


    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this growth group.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: G_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: bool(G_ZZ.has_coerce_map_from(G_QQ))  # indirect doctest
            False
            sage: bool(G_QQ.has_coerce_map_from(G_ZZ))  # indirect doctest
            True

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P_x_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_x_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: bool(P_x_ZZ.has_coerce_map_from(P_x_QQ))  # indirect doctest
            False
            sage: bool(P_x_QQ.has_coerce_map_from(P_x_ZZ))  # indirect doctest
            True
            sage: P_y_ZZ = agg.MonomialGrowthGroup(ZZ, 'y')
            sage: bool(P_y_ZZ.has_coerce_map_from(P_x_ZZ))  # indirect doctest
            False
            sage: bool(P_x_ZZ.has_coerce_map_from(P_y_ZZ))  # indirect doctest
            False
            sage: bool(P_y_ZZ.has_coerce_map_from(P_x_QQ))  # indirect doctest
            False
            sage: bool(P_x_QQ.has_coerce_map_from(P_y_ZZ))  # indirect doctest
            False

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P_x_ZZ = agg.GrowthGroup('ZZ^x')
            sage: P_x_QQ = agg.GrowthGroup('QQ^x')
            sage: bool(P_x_ZZ.has_coerce_map_from(P_x_QQ))  # indirect doctest
            False
            sage: bool(P_x_QQ.has_coerce_map_from(P_x_ZZ))  # indirect doctest
            True
            sage: P_y_ZZ = agg.GrowthGroup('ZZ^y')
            sage: bool(P_y_ZZ.has_coerce_map_from(P_x_ZZ))  # indirect doctest
            False
            sage: bool(P_x_ZZ.has_coerce_map_from(P_y_ZZ))  # indirect doctest
            False
            sage: bool(P_y_ZZ.has_coerce_map_from(P_x_QQ))  # indirect doctest
            False
            sage: bool(P_x_QQ.has_coerce_map_from(P_y_ZZ))  # indirect doctest
            False
        """
        if isinstance(S, GenericGrowthGroup) and self._var_ == S._var_:
            if self.base().has_coerce_map_from(S.base()):
                return True


    def gens_monomial(self):
        r"""
        Return a generator of this growth group, in case one exists.

        INPUT:

        Nothing.

        OUTPUT:

        An empty tuple.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GenericGrowthGroup(ZZ).gens_monomial()
            ()

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GrowthGroup('ZZ^x').gens_monomial()
            ()
        """
        return tuple()


    def gens(self):
        r"""
        Return a tuple of all generators (as a group) of this growth
        group.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple whose entries are growth elements.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.gens()
            (x,)
            sage: agg.MonomialGrowthGroup(ZZ, 'log(x)').gens()
            (log(x),)
        """
        return (self(raw_element=self.base().one()),)


    def gen(self, n=0):
        r"""
        Return the `n`-th generator (as a group) of this growth group.

        INPUT:

        - ``n`` -- default: `0`.

        OUTPUT:

        A :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.gen()
            x

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('QQ^x')
            sage: P.gen()
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        """
        return self.gens()[n]


    def ngens(self):
        r"""
        Return the number of generators (as a group) of this growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A Python integer.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.ngens()
            1
            sage: agg.MonomialGrowthGroup(ZZ, 'log(x)').ngens()
            1

        ::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('QQ^x')
            sage: P.ngens()
            0
        """
        return len(self.gens())


    def variable_names(self):
        r"""
        Return the names of the variables.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GenericGrowthGroup(ZZ).variable_names()
            ()

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ').variable_names()
            ('x',)
            sage: GrowthGroup('log(x)^ZZ').variable_names()
            ('x',)

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('QQ^x').variable_names()
            ('x',)
            sage: GrowthGroup('QQ^(x*log(x))').variable_names()
            ('x',)
        """
        return self._var_.variable_names()


    CartesianProduct = CartesianProductGrowthGroups


from sage.categories.pushout import ConstructionFunctor
class AbstractGrowthGroupFunctor(ConstructionFunctor):
    r"""
    """
    _functor_name = 'AbstractGrowthGroup'
    rank = 13

    def __init__(self, var, domain):
        if var is None:
            var = Variable('')
        elif not isinstance(var, Variable):
            var = Variable(var)
        self.var = var
        super(ConstructionFunctor, self).__init__(
            domain, sage.categories.groups.Groups())


    def _repr_(self):
        return '%s[%s]' % (self._functor_name, self.var)


    def merge(self, other):
        if self == other:
            return self


    def __eq__(self, other):
        return type(self) == type(other) and self.var == other.var


    def __ne__(self, other):
        return not self.__eq__(other)


class MonomialGrowthElement(GenericGrowthElement):
    r"""
    An implementation of monomial growth elements.

    INPUT:

    - ``parent`` -- a :class:`GenericGrowthGroup`.

    - ``raw_element`` -- an element from the base ring of the parent.

      This ``raw_element`` is the exponent of the created monomial
      growth element.

    A monomial growth element represents a term of the type
    `\operatorname{variable}^{\operatorname{exponent}}`. The multiplication
    corresponds to the addition of the exponents.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
        sage: e1 = P(1); e1
        1
        sage: e2 = P(raw_element=2); e2
        x^2
        sage: e1 == e2
        False
        sage: P.le(e1, e2)
        True
        sage: P.le(e1, P.gen()) and P.le(P.gen(), e2)
        True
    """

    @property
    def exponent(self):
        r"""
        The exponent of this growth element.

        EXAMPLES:

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P(x^42).exponent
            42
        """
        return self._raw_element_


    def _repr_(self):
        r"""
        A representation string for this monomial growth element.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P(1)._repr_()
            '1'
            sage: P(x^5)  # indirect doctest
            x^5
            sage: P(x^(1/2))  # indirect doctest
            x^(1/2)

        TESTS::

            sage: P(x^-1)  # indirect doctest
            1/x
            sage: P(x^-42)  # indirect doctest
            x^(-42)
        """
        from sage.rings.integer_ring import ZZ

        var = repr(self.parent()._var_)
        if self.exponent == 0:
            return '1'
        elif self.exponent == 1:
            return var
        elif self.exponent == -1:
            return '1/' + var
        elif self.exponent in ZZ and self.exponent > 0:
            return var + '^' + str(self.exponent)
        else:
            return var + '^(' + str(self.exponent) + ')'


    def _mul_(self, other):
        r"""
        Multiply this monomial growth element with another.

        INPUT:

        - ``other`` -- a :class:`MonomialGrowthElement`

        OUTPUT:

        The product as a :class:`MonomialGrowthElement`.

        .. NOTE::

            Two monomial growth elements are multiplied by adding
            their exponents.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: a = P(x^2)
            sage: b = P(x^3)
            sage: c = a._mul_(b); c
            x^5
            sage: c == a * b
            True
            sage: a * b * a  # indirect doctest
            x^7
        """
        return self.parent()(raw_element=self.exponent + other.exponent)


    def __invert__(self):
        r"""
        Return the multiplicative inverse of this monomial growth element.

        INPUT:

        Nothing.

        OUTPUT:

        The multiplicative inverse as a :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: e1 = P(raw_element=2)
            sage: e2 = e1.__invert__(); e2
            x^(-2)
            sage: e2 == ~e1
            True
        """
        return self.parent()(raw_element=-self.exponent)


    def __pow__(self, power):
        r"""
        Raises this growth element to the given ``power``.

        INPUT:

        - ``power`` -- a number. This can be anything that is a
          valid right hand side of ``*`` with elements of the
          parent's base.

        OUTPUT:

        The result of this exponentiation, a :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: x = P.gen()
            sage: a = x^7; a
            x^7
            sage: a^(1/2)
            Traceback (most recent call last):
            ...
            ValueError: Growth Group x^ZZ disallows taking x^7 to the power of 1/2.
            sage: P = agg.MonomialGrowthGroup(QQ, 'x')
            sage: b = P.gen()^(7/2); b
            x^(7/2)
            sage: b^12
            x^42
        """
        new_exponent = self.exponent * power
        P = self.parent()
        if new_exponent in P.base():
            return P(raw_element=new_exponent)
        else:
            raise ValueError('%s disallows taking %s to the power '
                             'of %s.' % (P, self, power))


    def _le_(self, other):
        r"""
        Return if this :class:`MonomialGrowthElement` is at most
        (less than or equal to) ``other``.

        INPUT:

        - ``other`` -- a :class:`MonomialGrowthElement`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function compares two instances of
            :class:`MonomialGrowthElement`.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P_ZZ.gen() <= P_QQ.gen()^2  # indirect doctest
            True
        """
        return self.exponent <= other.exponent


class MonomialGrowthGroup(GenericGrowthGroup):
    r"""
    A growth group dealing with powers of a fixed object/symbol.

    The elements :class:`MonomialGrowthElement` of this group represent powers
    of a fixed base; the group law is the multiplication, which corresponds
    to the addition of the exponents of the monomials.

    INPUT:

    - ``base`` -- one of SageMath's parents, out of which the elements
      get their data (``raw_element``).

      As monomials are represented by this group, the elements in
      ``base`` are the exponents of these monomials.

    - ``var`` -- an object.

      The string representation of ``var`` acts as a base of the
      monomials represented by this group.

    - ``category`` -- (default: ``None``) the category of the newly
      created growth group. It has to be a subcategory of ``Join of
      Category of groups and Category of posets``. This is also the
      default category if ``None`` is specified.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: P = agg.MonomialGrowthGroup(ZZ, 'x'); P
        Growth Group x^ZZ
        sage: agg.MonomialGrowthGroup(ZZ, log(SR.var('y')))
        Growth Group log(y)^ZZ

    .. SEEALSO::

        :class:`GenericGrowthGroup`
    """

    # enable the category framework for elements
    Element = MonomialGrowthElement


    def _repr_short_(self):
        r"""
        A short representation string of this monomial growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'a')  # indirect doctest
            Growth Group a^ZZ


        TESTS::

            sage: agg.MonomialGrowthGroup(ZZ, 'a')._repr_short_()
            'a^ZZ'
            sage: agg.MonomialGrowthGroup(QQ, 'a')._repr_short_()
            'a^QQ'
            sage: agg.MonomialGrowthGroup(PolynomialRing(QQ, 'x'), 'a')._repr_short_()
            'a^(Univariate Polynomial Ring in x over Rational Field)'
        """
        return '%s^%s' % (self._var_, parent_to_repr_short(self.base()))


    def _convert_(self, data):
        r"""
        Converts ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P._convert_('icecream') is None
            True
            sage: P(1)  # indirect doctest
            1
            sage: P('x')  # indirect doctest
            x

        ::

            sage: P(x)  # indirect doctest
            x
            sage: P(x^-333)  # indirect doctest
            x^(-333)
            sage: P(log(x)^2)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert log(x)^2.

        ::

            sage: PR.<x> = ZZ[]; x.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: P(x^2)  # indirect doctest
            x^2

        ::

            sage: PSR.<x> = ZZ[[]]
            sage: P(x^42)  # indirect doctest
            x^42
            sage: P(x^12 + O(x^17))
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert x^12 + O(x^17).

        ::

            sage: R.<w,x> = ZZ[]
            sage: P(x^4242)  # indirect doctest
            x^4242
            sage: P(w^4242)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert w^4242.

        ::

            sage: PSR.<w,x> = ZZ[[]]
            sage: P(x^7)  # indirect doctest
            x^7
            sage: P(w^7)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert w^7.

        ::

            sage: P('x^7')
            x^7
            sage: P('1/x')
            1/x
            sage: P('x^(-2)')
            x^(-2)
            sage: P('x^-2')
            x^(-2)
        """
        if data == 1:
            return self.base().zero()
        var = repr(self._var_)
        if str(data) == var:
            return self.base().one()

        try:
            P = data.parent()
        except AttributeError:
            if var not in str(data):
                return  # this has to end here

            elif str(data) == '1/' + var:
                return self.base()(-1)
            elif str(data).startswith(var + '^'):
                return self.base()(str(data).replace(var + '^', '')
                                   .replace('(', '').replace(')', ''))
            else:
                return  # end of parsing


        from sage.symbolic.ring import SR
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        from sage.rings.polynomial.multi_polynomial_ring_generic import \
            MPolynomialRing_generic
        from sage.rings.power_series_ring import PowerSeriesRing_generic
        import operator
        if P is SR:
            if data.operator() == operator.pow:
                base, exponent = data.operands()
                if str(base) == var:
                    return exponent
        elif isinstance(P, (PolynomialRing_general, MPolynomialRing_generic)):
            if data.is_monomial() and len(data.variables()) == 1:
                if var == str(data.variables()[0]):
                    return data.degree()
        elif isinstance(P, PowerSeriesRing_generic):
            if hasattr(data, 'variables') and len(data.variables()) == 1:
                from sage.rings.integer_ring import ZZ
                if data.is_monomial() and data.precision_absolute() not in ZZ:
                    if var == str(data.variables()[0]):
                        return data.degree()
            elif var == str(data.variable()[0]):
                from sage.rings.integer_ring import ZZ
                if data.is_monomial() and data.precision_absolute() not in ZZ:
                    return data.degree()


    def gens_monomial(self):
        r"""
        Return a tuple containing generators of this growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple containing elements of this growth group.

        .. NOTE::

            If a :class:`MonomialGrowthGroup` models a logarithmic
            growth group (by having a variable name of the form
            ``log(...)``), an empty tuple is returned.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'x').gens_monomial()
            (x,)
            sage: agg.MonomialGrowthGroup(QQ, 'log(x)').gens_monomial()
            ()
        """
        if not self._var_.is_monomial():
            return tuple()
        return (self(raw_element=self.base().one()),)


    def construction(self):
        r"""
        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ').construction()
            (MonomialGrowthGroup[x], Integer Ring)
       """
        return MonomialGrowthGroupFunctor(self._var_), self.base()


class MonomialGrowthGroupFunctor(AbstractGrowthGroupFunctor):
    r"""
    """

    _functor_name = 'MonomialGrowthGroup'


    def __init__(self, var):
        super(MonomialGrowthGroupFunctor, self).__init__(var,
            sage.categories.commutative_additive_groups.CommutativeAdditiveGroups())


    def _apply_functor(self, base):
        return MonomialGrowthGroup(base, self.var)


class ExponentialGrowthElement(GenericGrowthElement):
    r"""
    An implementation of exponential growth elements.

    INPUT:

    - ``parent`` -- an :class:`ExponentialGrowthGroup`.

    - ``raw_element`` -- an element from the base ring of the parent.

      This ``raw_element`` is the base of the created exponential
      growth element.

    An exponential growth element represents a term of the type
    `\operatorname{base}^{\operatorname{variable}}`. The multiplication
    corresponds to the multiplication of the bases.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: P = agg.GrowthGroup('ZZ^x')
        sage: e1 = P(1); e1
        1
        sage: e2 = P(raw_element=2); e2
        2^x
        sage: e1 == e2
        False
        sage: P.le(e1, e2)
        True
        sage: P.le(e1, P(1)) and P.le(P(1), e2)
        True
    """

    @property
    def base(self):
        r"""
        The base of this exponential growth element.

        EXAMPLES:

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('ZZ^x')
            sage: P(42^x).base
            42
        """
        return self._raw_element_


    def _repr_(self):
        r"""
        A representation string for this exponential growth element.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('QQ^x')
            sage: P(1)._repr_()
            '1'
            sage: P(5^x)  # indirect doctest
            5^x
            sage: P((1/2)^x)  # indirect doctest
            (1/2)^x

        TESTS::

            sage: P((-1)^x)  # indirect doctest
            (-1)^x
        """
        from sage.rings.integer_ring import ZZ

        var = repr(self.parent()._var_)
        if self.base == 1:
            return '1'
        elif self.base in ZZ and self.base > 0 or str(self.base).startswith('sqrt'):
            return str(self.base) + '^' + var
        else:
            return '(' + str(self.base) + ')^' + var


    def _mul_(self, other):
        r"""
        Multiply this exponential growth element with another.

        INPUT:

        - ``other`` -- a :class:`ExponentialGrowthElement`

        OUTPUT:

        The product as a :class:`ExponentialGrowthElement`.

        .. NOTE::

            Two exponential growth elements are multiplied by
            multiplying their bases.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('ZZ^x')
            sage: a = P(2^x)
            sage: b = P(3^x)
            sage: c = a._mul_(b); c
            6^x
            sage: c == a * b
            True
            sage: a * b * a  # indirect doctest
            12^x
        """
        return self.parent()(raw_element=self.base * other.base)


    def __invert__(self):
        r"""
        Return the multiplicative inverse of this exponential growth element.

        INPUT:

        Nothing.

        OUTPUT:

        The multiplicative inverse as a :class:`ExponentialGrowthElement`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('ZZ^x')
            sage: e1 = P(raw_element=2)
            sage: e2 = e1.__invert__(); e2
            (1/2)^x
            sage: e2 == ~e1
            True
        """
        new_base = 1 / self.base
        try:
            return self.parent()(raw_element=new_base)
        except (ValueError, TypeError):
            new_parent = ExponentialGrowthGroup(new_base.parent(),
                                                self.parent()._var_)
            return new_parent(raw_element=new_base)


    def __pow__(self, power):
        r"""
        Takes this growth element to the given ``power``.

        INPUT:

        - ``power`` -- a number. This can anything that is valid to be
          on the right hand side of ``*`` with an elements of the
          parent's base.

        OUTPUT:

        The result of this exponentiation a :class:`ExponentialGrowthElement`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.GrowthGroup('ZZ^x')
            sage: a = P(7^x); a
            7^x
            sage: b = a^(1/2); b
            sqrt(7)^x
            sage: b.parent()
            Growth Group SR^x
            sage: b^12
            117649^x
        """
        new_base = self.base ** power
        try:
            return self.parent()(raw_element=new_base)
        except (ValueError, TypeError):
            pass

        new_parent = ExponentialGrowthGroup(new_base.parent(),
                                            self.parent()._var_)
        return new_parent(raw_element=new_base)


    def _le_(self, other):
        r"""
        Return if this :class:`ExponentialGrowthElement` is at most
        (less than or equal to) ``other``.

        INPUT:

        - ``other`` -- a :class:`ExponentialGrowthElement`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function compares two instances of
            :class:`ExponentialGrowthElement`.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P_ZZ = agg.GrowthGroup('ZZ^x')
            sage: P_SR = agg.GrowthGroup('SR^x')
            sage: P_ZZ(2^x) <= P_SR(sqrt(3)^x)^2  # indirect doctest
            True
        """
        return bool(abs(self.base) <= abs(other.base))


class ExponentialGrowthGroup(GenericGrowthGroup):
    r"""
    A growth group dealing with expressions involving a fixed
    variable/symbol as the exponent.

    The elements :class:`ExponentialGrowthElement` of this group
    represent exponential functions with bases from a fixed base
    ring; the group law is the multiplication.

    INPUT:

    - ``base`` -- one of SageMath's parents, out of which the elements
      get their data (``raw_element``).

      As exponential expressions are represented by this group,
      the elements in ``base`` are the bases of these exponentials.

    - ``var`` -- an object.

      The string representation of ``var`` acts as an exponent of the
      elements represented by this group.

    - ``category`` -- (default: ``None``) the category of the newly
      created growth group. It has to be a subcategory of ``Join of
      Category of groups and Category of posets``. This is also the
      default category if ``None`` is specified.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: P = agg.ExponentialGrowthGroup(QQ, 'x'); P
        Growth Group QQ^x

    .. SEEALSO::

        :class:`GenericGrowthGroup`
    """

    # enable the category framework for elements
    Element = ExponentialGrowthElement


    def _repr_short_(self):
        r"""
        A short representation string of this exponential growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.ExponentialGrowthGroup(QQ, 'a')  # indirect doctest
            Growth Group QQ^a


        TESTS::

            sage: agg.ExponentialGrowthGroup(QQ, 'a')._repr_short_()
            'QQ^a'
            sage: agg.ExponentialGrowthGroup(PolynomialRing(QQ, 'x'), 'a')._repr_short_()
            '(Univariate Polynomial Ring in x over Rational Field)^a'
        """
        return '%s^%s' % (parent_to_repr_short(self.base()), self._var_)


    def _convert_(self, data):
        r"""
        Converts given ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.ExponentialGrowthGroup(ZZ, 'x')
            sage: P._convert_('icecream') is None
            True
            sage: P(1)  # indirect doctest
            1
            sage: P('2^x')  # indirect doctest
            2^x

        ::

            sage: P(2^x)  # indirect doctest
            2^x
            sage: P((-333)^x)  # indirect doctest
            (-333)^x
            sage: P(0)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert 0.

        ::

            sage: P('7^x')
            7^x
            sage: P('(-2)^x')
            (-2)^x

        ::

            sage: P = agg.GrowthGroup('SR^x')
            sage: P(sqrt(3)^x)
            sqrt(3)^x
            sage: P((3^(1/3))^x)
            (3^(1/3))^x
        """
        if data == 1 or data == '1':
            return self.base().one()
        var = repr(self._var_)
        try:
            P = data.parent()
        except AttributeError:
            import re
            if var not in str(data):
                return  # this has to end here

            elif str(data).endswith('^' + var):
                return self.base()(str(data).replace('^' + var, '')
                                   .replace('(', '').replace(')', ''))
            else:
                return  # end of parsing


        from sage.symbolic.ring import SR
        import operator
        from sage.symbolic.operators import mul_vararg
        if P is SR:
            if data.operator() == operator.pow:
                base, exponent = data.operands()
                if str(exponent) == var:
                    return base
                elif exponent.operator() == mul_vararg:
                    return base ** (exponent / SR(var))


    def gens(self):
        r"""
        Return a tuple of all generators (as a group) of this growth
        group.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple whose entries are growth elements.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: E = GrowthGroup('ZZ^x')
            sage: E.gens()
            ()
        """
        return tuple()


    def construction(self):
        r"""
        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('QQ^x').construction()
            (ExponentialGrowthGroup[x], Rational Field)
        """
        return ExponentialGrowthGroupFunctor(self._var_), self.base()


class ExponentialGrowthGroupFunctor(AbstractGrowthGroupFunctor):
    r"""
    """

    _functor_name = 'ExponentialGrowthGroup'


    def __init__(self, var):
        super(ExponentialGrowthGroupFunctor, self).__init__(var,
            sage.categories.groups.Groups())


    def _apply_functor(self, base):
        return ExponentialGrowthGroup(base, self.var)


class GrowthGroupFactory(sage.structure.factory.UniqueFactory):
    r"""
    A factory creating asymptotic growth groups.

    INPUT:

    - ``specification`` -- a string.

    OUTPUT:

    An asymptotic growth group.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: agg.GrowthGroup('x^ZZ')
        Growth Group x^ZZ
        sage: agg.GrowthGroup('log(x)^QQ')
        Growth Group log(x)^QQ

    This factory can also be used to construct Cartesian products
    of growth groups::

        sage: agg.GrowthGroup('x^ZZ * y^ZZ')
        Growth Group x^ZZ * y^ZZ
        sage: agg.GrowthGroup('x^ZZ * log(x)^ZZ')
        Growth Group x^ZZ * log(x)^ZZ
        sage: agg.GrowthGroup('x^ZZ * log(x)^ZZ * y^QQ')
        Growth Group x^ZZ * log(x)^ZZ * y^QQ
        sage: agg.GrowthGroup('QQ^x * x^ZZ * y^QQ * QQ^z')
        Growth Group QQ^x * x^ZZ * y^QQ * QQ^z
    """
    def create_key_and_extra_args(self, specification, **kwds):
        r"""
        Given the arguments and keyword, create a key that uniquely
        determines this object.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GrowthGroup.create_key_and_extra_args('x^ZZ')
            (('x^ZZ',), {})
            sage: agg.GrowthGroup.create_key_and_extra_args('asdf')
            Traceback (most recent call last):
            ...
            ValueError: 'asdf' is not a valid string describing a growth group.
            sage: agg.GrowthGroup.create_key_and_extra_args('log(x)^ZZ * y^QQ')
            (('log(x)^ZZ', 'y^QQ'), {})
            sage: agg.GrowthGroup.create_key_and_extra_args('log(x)**ZZ * y**QQ')
            (('log(x)**ZZ', 'y**QQ'), {})
            sage: agg.GrowthGroup.create_key_and_extra_args('a^b * * c^d')
            Traceback (most recent call last):
            ...
            ValueError: 'a^b * * c^d' is invalid since a '*' follows a '*'
            sage: agg.GrowthGroup.create_key_and_extra_args('a^b * (c*d^e)')
            (('a^b', 'c*d^e'), {})
        """
        factors = list()
        balanced = True
        if specification and specification[0] == '*':
            raise ValueError("'%s' is invalid since it starts with a '*'." %
                             (specification,))
        for s in specification.split('*'):
            if not s:
                factors[-1] += '*'
                balanced = False
                continue
            if not s.strip():
                raise ValueError("'%s' is invalid since a '*' follows a '*'" %
                                 (specification,))
            if not balanced:
                s = factors.pop() + '*' + s
            balanced = s.count('(') == s.count(')')
            factors.append(s)

        def strip(s):
            s = s.strip()
            if not s:
                return s
            if s[0] == '(' and s[-1] == ')':
                s = s[1:-1]
            return s.strip()

        factors = tuple(strip(f) for f in factors)

        for f in factors:
            if '^' not in f and '**' not in f:
                raise ValueError("'%s' is not a valid string describing "
                                 "a growth group." % (f,))

        return factors, kwds


    def create_object(self, version, factors, **kwds):
        r"""
        Create an object from the given arguments.

        TESTS::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: agg.GrowthGroup('as^df')
            Traceback (most recent call last):
            ...
            ValueError: 'as^df' is not a valid string describing a growth group.
            sage: agg.GrowthGroup('x^y^z')
            Traceback (most recent call last):
            ...
            ValueError: Cannot decode x^y^z.
        """

        groups = []
        for factor in factors:
            b_and_e = factor.split('^')
            if len(b_and_e) != 2:
                raise ValueError('Cannot decode %s.' % (factor,))
            (b, e) = b_and_e

            try:
                # monomial growth group: 'var^base'
                groups.append(
                    MonomialGrowthGroup(repr_short_to_parent(e), b, **kwds))
                continue
            except (TypeError, ValueError):
                pass

            try:
                # exponential growth group: 'base^var'
                groups.append(
                    ExponentialGrowthGroup(repr_short_to_parent(b), e, **kwds))
                continue
            except (TypeError, ValueError):
                pass

            raise ValueError("'%s' is not a valid string describing "
                             "a growth group." % (factor,))

        if len(groups) == 1:
            return groups[0]

        from sage.categories.cartesian_product import cartesian_product
        return cartesian_product(groups)


GrowthGroup = GrowthGroupFactory("GrowthGroup")
