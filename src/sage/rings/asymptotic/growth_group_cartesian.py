r"""
Growth Groups as Cartesian Products

See :mod:`sage.rings.asymptotic.growth_group` for a description.

AUTHORS:

- Daniel Krenn (2015-06-02): cartesian products
- Benjamin Hackl (2015-07)

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

TESTS::

    sage: from sage.rings.asymptotic.growth_group import GrowthGroup
    sage: A = GrowthGroup('QQ^x * x^ZZ'); A
    Growth Group QQ^x * x^ZZ
    sage: A.construction()
    (The cartesian_product functorial construction,
     (Growth Group QQ^x, Growth Group x^ZZ))
    sage: A.construction()[1][0].construction()
    (ExponentialGrowthGroup[x], Rational Field)
    sage: A.construction()[1][1].construction()
    (MonomialGrowthGroup[x], Integer Ring)
    sage: B = GrowthGroup('x^ZZ * y^ZZ'); B
    Growth Group x^ZZ * y^ZZ
    sage: B.construction()
    (The cartesian_product functorial construction,
     (Growth Group x^ZZ, Growth Group y^ZZ))
    sage: C = GrowthGroup('x^ZZ * log(x)^ZZ * y^ZZ'); C
    Growth Group x^ZZ * log(x)^ZZ * y^ZZ
    sage: C.construction()
    (The cartesian_product functorial construction,
     (Growth Group x^ZZ * log(x)^ZZ, Growth Group y^ZZ))
    sage: C.construction()[1][0].construction()
    (The cartesian_product functorial construction,
     (Growth Group x^ZZ, Growth Group log(x)^ZZ))
    sage: C.construction()[1][1].construction()
    (MonomialGrowthGroup[y], Integer Ring)

::

    sage: cm = sage.structure.element.get_coercion_model()
    sage: D = GrowthGroup('QQ^x * x^QQ')
    sage: cm.common_parent(A, D)
    Growth Group QQ^x * x^QQ
    sage: E = GrowthGroup('ZZ^x * x^QQ')
    sage: cm.common_parent(A, E)
    Growth Group QQ^x * x^QQ

::

    sage: A.an_element()
    (1/2)^x*x
    sage: tuple(E.an_element())
    (1, x^(1/2))
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


class CartesianProductFactory(sage.structure.factory.UniqueFactory):
    r"""
    Creates various types of cartesian products depending on its input.

    INPUT:

    - ``growth_groups`` -- a tuple (or other iterable) of growth groups.

    - ``order`` -- (default: ``None``) if specified, then this order
      is taken for comparing two cartesian product elements. If ``order`` is
      ``None`` this is determined automatically.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: A = GrowthGroup('x^ZZ'); A
        Growth Group x^ZZ
        sage: B = GrowthGroup('log(x)^ZZ'); B
        Growth Group log(x)^ZZ
        sage: C = cartesian_product([A, B]); C
        Growth Group x^ZZ * log(x)^ZZ
        sage: C._le_ == C.le_lex
        True
        sage: D = GrowthGroup('y^ZZ'); D
        Growth Group y^ZZ
        sage: E = cartesian_product([A, D]); E
        Growth Group x^ZZ * y^ZZ
        sage: E._le_ == E.le_components
        True
        sage: F = cartesian_product([C, D]); F
        Growth Group x^ZZ * log(x)^ZZ * y^ZZ
        sage: F._le_ == F.le_components
        True
        sage: cartesian_product([A, E]); G
        Traceback (most recent call last):
        ...
        ValueError: Growth groups (Growth Group x^ZZ, Growth Group x^ZZ * y^ZZ)
        do not have pairwise disjoint variables.
        sage: cartesian_product([A, B, D])
        Growth Group x^ZZ * log(x)^ZZ * y^ZZ

    TESTS::

        sage: from sage.rings.asymptotic.growth_group_cartesian import CartesianProductFactory
        sage: CartesianProductFactory('factory')([A, B], category=Groups() & Posets())
        Growth Group x^ZZ * log(x)^ZZ
        sage: CartesianProductFactory('factory')([], category=Sets())
        Traceback (most recent call last):
        ...
        TypeError: Cannot create cartesian product without factors.
    """
    def create_key_and_extra_args(self, growth_groups, category, **kwds):
        r"""
        Given the arguments and keywords, create a key that uniquely
        determines this object.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group_cartesian import CartesianProductFactory
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: A = GrowthGroup('x^ZZ')
            sage: CartesianProductFactory('factory').create_key_and_extra_args(
            ....:     [A], category=Sets(), order='blub')
            (((Growth Group x^ZZ,), Category of sets), {'order': 'blub'})
        """
        return (tuple(growth_groups), category), kwds


    def create_object(self, version, args, **kwds):
        r"""
        Create an object from the given arguments.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: cartesian_product([GrowthGroup('x^ZZ')])  # indirect doctest
            Growth Group x^ZZ
        """
        growth_groups, category = args
        if not growth_groups:
            raise TypeError('Cannot create cartesian product without factors.')
        order = kwds.pop('order', None)
        if order is not None:
            return GenericProduct(growth_groups, category, order=order, **kwds)

        vg = tuple((g.variable_names(), g) for g in growth_groups)

        # check if all groups have a variable
        if not all(v for v, _ in vg):
            raise NotImplementedError('Growth groups %s have no variable.' %
                                      tuple(g for g in growth_groups
                                            if not g.variable_names()))

        # sort by variables
        from itertools import groupby, product
        vgs = tuple((v, tuple(gs)) for v, gs in
                    groupby(sorted(vg, key=lambda k: k[0]), key=lambda k: k[0]))

        # check if variables are pairwise disjoint
        for u, w in product(iter(v for v, _ in vgs), repeat=2):
            if u != w and not set(u).isdisjoint(set(w)):
                raise ValueError('Growth groups %s do not have pairwise disjoint '
                                 'variables.' % (growth_groups,))

        # build cartesian products
        u_groups = list()
        for _, gs in vgs:
            gs = tuple(g for _, g in gs)
            if len(gs) > 1:
                u_groups.append(UnivariateProduct(gs, category, **kwds))
            else:
                u_groups.append(gs[0])

        if len(u_groups) > 1:
            m_group = MultivariateProduct(tuple(u_groups), category, **kwds)
        else:
            m_group = u_groups[0]
        return m_group


CartesianProductGrowthGroups = CartesianProductFactory('CartesianProductGrowthGroups')


from sage.sets.cartesian_product import CartesianProductPosets
from growth_group import GenericGrowthGroup
class GenericProduct(CartesianProductPosets, GenericGrowthGroup):
    r"""
    A cartesian product of growth groups.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: P = agg.MonomialGrowthGroup(QQ, 'x')
        sage: L = agg.MonomialGrowthGroup(ZZ, 'log(x)')
        sage: C = cartesian_product([P, L], order='lex'); C
        Growth Group x^QQ * log(x)^ZZ
        sage: C.an_element()
        x^(1/2)*log(x)

    ::

        sage: Px = agg.MonomialGrowthGroup(QQ, 'x')
        sage: Lx = agg.MonomialGrowthGroup(ZZ, 'log(x)')
        sage: Cx = cartesian_product([Px, Lx], order='lex')
        sage: Py = agg.MonomialGrowthGroup(QQ, 'y')
        sage: C = cartesian_product([Cx, Py], order='components'); C
        Growth Group x^QQ * log(x)^ZZ * y^QQ
        sage: C.an_element()
        x^(1/2)*log(x)*y^(1/2)

    .. SEEALSO:

        :class:`~sage.sets.cartesian_product.CartesianProduct`,
        :class:`~sage.sets.cartesian_product.CartesianProductPosets`.
    """

    __classcall__ = CartesianProductPosets.__classcall__


    def __init__(self, sets, category, **kwds):
        r"""
        See :class:`GenericProduct` for details.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ * y^ZZ')  # indirect doctest
            Growth Group x^ZZ * y^ZZ
        """
        order = kwds.pop('order')
        CartesianProductPosets.__init__(self, sets, category, order, **kwds)

        vars = sum(iter(factor.variable_names()
                        for factor in self.cartesian_factors()),
                   tuple())
        from itertools import groupby
        from growth_group import Variable
        Vars = Variable(tuple(v for v, _ in groupby(vars)), repr=self._repr_short_())

        GenericGrowthGroup.__init__(self, sets[0], Vars, self.category(), **kwds)


    __hash__ = CartesianProductPosets.__hash__


    def some_elements(self):
        r"""
        Return some elements of this cartesian product of growth groups.

        See :class:`TestSuite` for a typical use case.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from itertools import islice
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('QQ^y * x^QQ * log(x)^ZZ')
            sage: tuple(islice(G.some_elements(), 10))
            (x^(1/2)*(1/2)^y,
             x^(-1/2)*log(x)*(-1/2)^y,
             x^2*log(x)^(-1)*2^y,
             x^(-2)*log(x)^2*(-2)^y,
             log(x)^(-2)*0^y,
             x*log(x)^3,
             x^(-1)*log(x)^(-3)*(-1)^y,
             x^42*log(x)^4*42^y,
             x^(2/3)*log(x)^(-4)*(2/3)^y,
             x^(-2/3)*log(x)^5*(-2/3)^y)
        """
        from itertools import izip
        return iter(
            self(c) for c in
            izip(*tuple(F.some_elements() for F in self.cartesian_factors())))


    def _element_constructor_(self, data):
        r"""
        Converts the given object to an element of this cartesian
        product.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ * y^ZZ')
            sage: G_log = GrowthGroup('x^ZZ * log(x)^ZZ * y^ZZ')

        Conversion from the symbolic ring works::

            sage: x,y = var('x y')
            sage: G(x^-3*y^2)
            x^(-3)*y^2
            sage: G(x^4), G(y^2)
            (x^4, y^2)
            sage: G(1)
            1

        Even more complex expressions can be parsed::

            sage: G_log(x^42*log(x)^-42*y^42)
            x^42*log(x)^(-42)*y^42

        TESTS::

            sage: G = GrowthGroup('x^ZZ * y^ZZ')
            sage: G('x'), G('y')
            (x, y)

        ::

            sage: G_log(log(x))
            log(x)

        ::

            sage: G(G.cartesian_factors()[0].gen())
            x

        ::

            sage: GrowthGroup('QQ^x * x^QQ')(['x^(1/2)'])
            x^(1/2)
            sage: l = GrowthGroup('x^ZZ * log(x)^ZZ')(['x', 'log(x)']); l
            x*log(x)
            sage: type(l)
            <class 'sage.rings.asymptotic.growth_group_cartesian.UnivariateProduct_with_category.element_class'>
            sage: GrowthGroup('QQ^x * x^QQ')(['2^log(x)'])
            Traceback (most recent call last):
            ...
            ValueError: ['2^log(x)'] is not in Growth Group QQ^x * x^QQ.
            > *previous* ValueError: 2^log(x) is not in any of the factors of
            Growth Group QQ^x * x^QQ
            sage: GrowthGroup('QQ^x * x^QQ')(['2^log(x)', 'x^55'])
            Traceback (most recent call last):
            ...
            ValueError: ['2^log(x)', 'x^55'] is not in Growth Group QQ^x * x^QQ.
            > *previous* ValueError: 2^log(x) is not in any of the factors of
            Growth Group QQ^x * x^QQ
        """
        def convert_factors(data, raw_data):
            try:
                return self._convert_factors_(data)
            except ValueError as e:
                from misc import combine_exceptions
                raise combine_exceptions(
                    ValueError('%s is not in %s.' % (raw_data, self)), e)

        if data == 1:
            return self.one()

        elif data is None:
            raise ValueError('%s cannot be converted.' % (data,))

        elif type(data) == self.element_class and data.parent() == self:
            return data

        elif isinstance(data, self.element_class):
            return self.element_class(self, data)

        elif isinstance(data, str):
            from misc import split_str_by_mul
            return convert_factors(split_str_by_mul(data), data)

        elif hasattr(data, 'parent'):
            P = data.parent()

            if P is self:
                return data

            elif P is sage.symbolic.ring.SR:
                from sage.symbolic.operators import mul_vararg
                if data.operator() == mul_vararg:
                    return convert_factors(data.operands(), data)

            # room for other parents (e.g. polynomial ring et al.)

        else:
            try:
                return super(GenericProduct, self)._element_constructor_(data)
            except ValueError:
                pass
            if isinstance(data, (tuple, list,
                                 sage.sets.cartesian_product.CartesianProduct.Element)):
                return convert_factors(tuple(data), data)

        return convert_factors((data,), data)


    _repr_ = GenericGrowthGroup._repr_


    def _repr_short_(self):
        r"""
        A short (shorter than :meth:`._repr_`) representation string
        for this cartesian product of growth groups.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: P = agg.MonomialGrowthGroup(QQ, 'x')
            sage: L = agg.MonomialGrowthGroup(ZZ, 'log(x)')
            sage: cartesian_product([P, L], order='lex')._repr_short_()
            'x^QQ * log(x)^ZZ'
        """
        return ' * '.join(S._repr_short_() for S in self.cartesian_factors())


    def _convert_factors_(self, factors):
        r"""
        Helper method. Try to convert some ``factors`` to an
        element of one of the cartesian factors and returns the product of
        all these factors.

        INPUT:

        - ``factors`` -- a tuple or other iterable.

        OUTPUT:

        An element of this cartesian product.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ * log(x)^QQ * y^QQ')
            sage: e1 = G._convert_factors_([x^2])
            sage: (e1, e1.parent())
            (x^2, Growth Group x^ZZ * log(x)^QQ * y^QQ)
        """
        from sage.misc.misc_c import prod

        def get_factor(data):
            for factor in self.cartesian_factors():
                try:
                    return factor, factor(data)
                except (ValueError, TypeError):
                    pass
            raise ValueError('%s is not in any of the factors of %s' % (data, self))

        return prod(list(self.cartesian_injection(*get_factor(f))
                         for f in factors))


    def cartesian_injection(self, factor, element):
        r"""
        Injects the given element into this cartesian product at the given factor.

        INPUT:

        - ``factor`` -- a growth group (a factor of this cartesian product).

        - ``element`` -- an element of ``factor``.

        OUTPUT:

        An element of this cartesian product.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ * y^QQ')
            sage: G.cartesian_injection(G.cartesian_factors()[1], 'y^7')
            y^7
        """
        return self(tuple((f.one() if f != factor else element)
                          for f in self.cartesian_factors()))


    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this growth group.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: A = GrowthGroup('QQ^x * x^QQ')
            sage: B = GrowthGroup('QQ^x * x^ZZ')
            sage: A.has_coerce_map_from(B)
            True
            sage: B.has_coerce_map_from(A)
            False
        """
        if CartesianProductPosets.has_coerce_map_from(self, S):
            return True

        elif isinstance(S, GenericProduct):
            factors = S.cartesian_factors()
        else:
            factors = (S,)

        if all(any(g.has_coerce_map_from(f) for g in self.cartesian_factors())
               for f in factors):
            return True


    def _pushout_(self, other):
        r"""
        Construct the pushout of this and the other growth group. This is called by
        :func:`sage.categories.pushout.pushout`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.categories.pushout import pushout
            sage: cm = sage.structure.element.get_coercion_model()
            sage: A = GrowthGroup('QQ^x * x^ZZ')
            sage: B = GrowthGroup('x^ZZ * log(x)^ZZ')
            sage: A._pushout_(B)
            Growth Group QQ^x * x^ZZ * log(x)^ZZ
            sage: pushout(A, B)
            Growth Group QQ^x * x^ZZ * log(x)^ZZ
            sage: cm.discover_coercion(A, B)
            ((map internal to coercion system -- copy before use)
             Conversion map:
               From: Growth Group QQ^x * x^ZZ
               To:   Growth Group QQ^x * x^ZZ * log(x)^ZZ,
             (map internal to coercion system -- copy before use)
             Conversion map:
               From: Growth Group x^ZZ * log(x)^ZZ
               To:   Growth Group QQ^x * x^ZZ * log(x)^ZZ)
            sage: cm.common_parent(A, B)
            Growth Group QQ^x * x^ZZ * log(x)^ZZ

        ::

            sage: C = GrowthGroup('QQ^x * x^QQ * y^ZZ')
            sage: D = GrowthGroup('x^ZZ * log(x)^QQ * QQ^z')
            sage: C._pushout_(D)
            Growth Group QQ^x * x^QQ * log(x)^QQ * y^ZZ * QQ^z
            sage: cm.common_parent(C, D)
            Growth Group QQ^x * x^QQ * log(x)^QQ * y^ZZ * QQ^z
            sage: A._pushout_(D)
            Growth Group QQ^x * x^ZZ * log(x)^QQ * QQ^z
            sage: cm.common_parent(A, D)
            Growth Group QQ^x * x^ZZ * log(x)^QQ * QQ^z
            sage: cm.common_parent(B, D)
            Growth Group x^ZZ * log(x)^QQ * QQ^z
            sage: cm.common_parent(A, C)
            Growth Group QQ^x * x^QQ * y^ZZ
            sage: E = GrowthGroup('log(x)^ZZ * y^ZZ')
            sage: cm.common_parent(A, E)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Growth Group QQ^x * x^ZZ' and 'Growth Group log(x)^ZZ * y^ZZ'

        ::

            sage: F = GrowthGroup('z^QQ')
            sage: pushout(C, F)
            Growth Group QQ^x * x^QQ * y^ZZ * z^QQ
        """
        from growth_group import GenericGrowthGroup, AbstractGrowthGroupFunctor
        from misc import merge_overlapping

        if isinstance(other, GenericProduct):
            Ofactors = other.cartesian_factors()
        elif isinstance(other, GenericGrowthGroup):
            Ofactors = (other,)
        elif (other.construction() is not None and
              isinstance(other.construction()[0], AbstractGrowthGroupFunctor)):
            Ofactors = (other,)
        else:
            return


        def pushout_univariate_factors(self, other, var, Sfactors, Ofactors):
            try:
                return merge_overlapping(
                    Sfactors, Ofactors,
                    lambda f: (type(f), f._var_.var_repr))
            except ValueError:
                pass

            cm = sage.structure.element.get_coercion_model()
            try:
                Z = cm.common_parent(*Sfactors+Ofactors)
                return (Z,), (Z,)
            except TypeError:
                pass

            def subfactors(F):
                for f in F:
                    if isinstance(f, GenericProduct):
                        for g in subfactors(f.cartesian_factors()):
                            yield g
                    else:
                        yield f

            try:
                return merge_overlapping(
                    tuple(subfactors(Sfactors)), tuple(subfactors(Ofactors)),
                    lambda f: (type(f), f._var_.var_repr))
            except ValueError:
                pass

            from sage.structure.coerce_exceptions import CoercionException
            raise CoercionException(
                'Cannot construct the pushout of %s and %s: The factors '
                'with variables %s are not overlapping, '
                'no common parent was found, and '
                'splitting the factors was unsuccessful.' % (self, other, var))


        class it:
            def __init__(self, it):
                self.it = it
                self.var = None
                self.factors = None
            def next(self):
                try:
                    self.var, factors = next(self.it)
                    self.factors = tuple(factors)
                except StopIteration:
                    self.var = None
                    self.factors = tuple()

        from itertools import groupby
        S = it(groupby(self.cartesian_factors(), key=lambda k: k.variable_names()))
        O = it(groupby(Ofactors, key=lambda k: k.variable_names()))

        newS = []
        newO = []

        S.next()
        O.next()
        while S.var is not None or O.var is not None:
            if S.var is not None and S.var < O.var:
                newS.extend(S.factors)
                newO.extend(S.factors)
                S.next()
            elif O.var is not None and S.var > O.var:
                newS.extend(O.factors)
                newO.extend(O.factors)
                O.next()
            else:
                SL, OL = pushout_univariate_factors(self, other, S.var,
                                                    S.factors, O.factors)
                newS.extend(SL)
                newO.extend(OL)
                S.next()
                O.next()

        assert(len(newS) == len(newO))

        if (len(self.cartesian_factors()) == len(newS) and
            len(other.cartesian_factors()) == len(newO)):
            # We had already all factors in each of the self and
            # other, thus splitting it in subproblems (one for
            # each factor) is the strategy to use. If a pushout is
            # possible :func:`sage.categories.pushout.pushout`
            # will manage this by itself.
            return

        from sage.categories.pushout import pushout
        from sage.categories.cartesian_product import cartesian_product
        return pushout(cartesian_product(newS), cartesian_product(newO))


    def gens_monomial(self):
        r"""
        Return a tuple containing generators of this growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple containing elements of this growth group.

        .. NOTE:

            This method calls the ``gens_monomial()`` method on the
            individual factors of this cartesian product and
            concatenates the respective outputs.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ * log(x)^ZZ * y^QQ * log(z)^ZZ')
            sage: G.gens_monomial()
            (x, y)

        TESTS::

            sage: all(g.parent() == G for g in G.gens_monomial())
            True
        """
        return sum(iter(
            tuple(self.cartesian_injection(factor, g) for g in factor.gens_monomial())
            for factor in self.cartesian_factors()),
                   tuple())


    def variable_names(self):
        r"""
        Return the names of the variables.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: GrowthGroup('x^ZZ * log(x)^ZZ * y^QQ * log(z)^ZZ').variable_names()
            ('x', 'y', 'z')
        """
        vars = sum(iter(factor.variable_names()
                        for factor in self.cartesian_factors()),
                   tuple())
        from itertools import groupby
        return tuple(v for v, _ in groupby(vars))


    class Element(CartesianProductPosets.Element):

        def _repr_(self):
            r"""
            A representation string for this cartesian product element.

            INPUT:

            Nothing.

            OUTPUT:

            A string.

            EXAMPLES::

                sage: import sage.rings.asymptotic.growth_group as agg
                sage: P = agg.MonomialGrowthGroup(QQ, 'x')
                sage: L = agg.MonomialGrowthGroup(ZZ, 'log(x)')
                sage: cartesian_product([P, L], order='lex').an_element()._repr_()
                'x^(1/2)*log(x)'
            """
            s = '*'.join(repr(v) for v in self.value if not v.is_one())
            if s == '':
                return '1'
            return s

        def log(self, base=None):
            r"""
            The logarithm of this element.

            INPUT:

            - ``base`` -- the base of the logarithm. If ``None``
              (default value) is used, the logarithm is the natural
              logarithm.

            OUTPUT:

            A growth element.

            EXAMPLES::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('x^ZZ * log(x)^ZZ')
                sage: x, = G.gens_monomial()
                sage: log(x)  # indirect doctest
                log(x)
                sage: log(x^5)
                Traceback (most recent call last):
                ...
                ValueError: The logarithm of x^5 cannot be constructed in Growth Group x^ZZ * log(x)^ZZ.

            ::

                sage: G = GrowthGroup('QQ^x * x^ZZ')
                sage: x, = G.gens_monomial()
                sage: el = x.rpow(2); el
                2^x
                sage: log(el)
                Traceback (most recent call last):
                ...
                ValueError: The logarithm of 2^x cannot be constructed in Growth Group QQ^x * x^ZZ.
                sage: log(el, base=2)
                x

            TESTS::

                sage: G = GrowthGroup("QQ['e']^x * x^ZZ")
                sage: x, = G.gens_monomial()
                sage: log(exp(x))
                x
            """
            lf = self.log_factor(base=base)
            if len(lf) == 1 and lf[0][1] == 1:
                return lf[0][0]
            if base:
                raise ValueError('The logarithm of %s with base %s cannot be '
                                 'constructed in %s.' % (self, base, self.parent()))

            raise ValueError('The logarithm of %s cannot be constructed in '
                             '%s.' % (self, self.parent()))

        def factor(self):
            r"""
            Return a factorization of this element into growth
            elements from atomic growth groups.

            INPUT:

            Nothing.

            OUTPUT:

            A list of growth elements.

            EXAMPLES::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('x^ZZ * log(x)^ZZ * y^ZZ')
                sage: x, y = G.gens_monomial()
                sage: x.factor()
                [x]
                sage: f = (x * y).factor(); f
                [x, y]
                sage: [factor.parent() for factor in f]
                [Growth Group x^ZZ, Growth Group y^ZZ]
                sage: f = (x * log(x)).factor(); f
                [x, log(x)]
                sage: [factor.parent() for factor in f]
                [Growth Group x^ZZ, Growth Group log(x)^ZZ]

            ::

                sage: G = GrowthGroup('x^ZZ * log(x)^ZZ * log(log(x))^ZZ * y^QQ')
                sage: x, y = G.gens_monomial()
                sage: f = (x * log(x) * y).factor(); f
                [x, log(x), y]
                sage: [factor.parent() for factor in f]
                [Growth Group x^ZZ, Growth Group log(x)^ZZ, Growth Group y^QQ]

            ::

                sage: G.one().factor()
                Traceback (most recent call last):
                ...
                ValueError: 1 does not have a factorization.
            """
            components = []
            for (component, factor) in zip(self.value,
                                           self.parent().cartesian_factors()):
                if component != factor.one():
                    if hasattr(component, 'factor'):
                        components = components + component.factor()
                    else:
                        components.append(component)

            if components:
                return components
            raise ValueError('%s does not have a factorization.' % (self,))

        def log_factor(self, base=None):
            r"""
            Return the logarithm of the factorization of this
            element.

            In particular, this function yields a list of pairs
            consisting of the growth and the corresponding
            coefficient resulting from taking the logarithm of every
            atomic factor of this element.

            INPUT:

            - ``base`` -- the base of the logarithm. If ``None``
              (default value) is used, the logarithm is the natural
              logarithm.

            OUTPUT:

            A list.

            .. SEEALSO::

                :meth:`factor`,
                :meth:`log`.

            EXAMPLES::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('QQ^x * x^ZZ * log(x)^ZZ * y^ZZ * log(y)^ZZ')
                sage: x, y = G.gens_monomial()
                sage: (x * y).log_factor()
                [[log(x), 1], [log(y), 1]]
                sage: (x^123).log_factor()
                [[log(x), 123]]
                sage: (G('2^x') * x^2).log_factor(base=2)
                [[x, 1], [log(x), 2/log(2)]]

            ::

                sage: G(1).log_factor()
                Traceback (most recent call last):
                ...
                ValueError: 1 does not have a factorization.

            ::

                sage: log(x).log_factor()
                Traceback (most recent call last):
                ...
                ValueError: Logarithm of log(x) cannot be constructed in Growth
                Group QQ^x * x^ZZ * log(x)^ZZ * y^ZZ * log(y)^ZZ.

            TESTS::

                sage: G = GrowthGroup("QQ['e']^x * x^ZZ * log(x)^ZZ")
                sage: x, = G.gens_monomial()
                sage: (exp(x) * x).log_factor()
                [[x, 1], [log(x), 1]]
            """
            P = self.parent()
            factors = self.factor()
            log_factors = []
            from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup, \
                MonomialGrowthGroup
            from sage.functions.log import log

            for factor in factors:
                try:
                    FP = factor.parent()
                    if isinstance(FP, MonomialGrowthGroup):
                        if base is None:
                            coef = factor._raw_element_
                        else:
                            coef = factor._raw_element_ / log(base)
                        growth = P('log(' + repr(FP._var_) + ')')
                        log_factors.append([growth, coef])

                    elif isinstance(FP, ExponentialGrowthGroup):
                        b = factor._raw_element_
                        if hasattr(b, 'is_monomial') and b.is_monomial():
                            if b.variable_name() == 'e' and not base:
                                coef = b.valuation()
                        else:
                            coef = log(b, base=base)

                        growth = P(repr(FP._var_))
                        log_factors.append([growth, coef])

                    else:
                        raise NotImplementedError('Taking the logarithm of %s '
                                                  'is not implemented.' %
                                                  (factor,))

                except (ValueError, TypeError):
                    raise ValueError('Logarithm of %s cannot be '
                                     'constructed in %s.' % (factor, P))

            return log_factors

        def rpow(self, base):
            r"""
            Take ``base`` to the power of this element.

            In other words, this is the exponential function with
            ``base`` as its base.

            INPUT:

            - ``base`` -- the base of the exponential function.

            OUTPUT:

            A growth element.

            EXAMPLES::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('QQ^x * x^ZZ')
                sage: x, = G.gens_monomial()
                sage: x.rpow(2)
                2^x
                sage: x.rpow(1/2)
                (1/2)^x

            ::

                sage: x.rpow(0)
                Traceback (most recent call last):
                ...
                ValueError: 0 is not an allowed base.
                sage: (x^2).rpow(2)
                Traceback (most recent call last):
                ...
                ValueError: Cannot construct 2^(x^2) in Growth Group QQ^x * x^ZZ.

            ::

                sage: G = GrowthGroup('QQ^(x*log(x)) * x^ZZ * log(x)^ZZ')
                sage: x, = G.gens_monomial()
                sage: (x * log(x)).rpow(2)
                2^(x*log(x))
            """
            P = self.parent()
            factors = self.factor()
            if base == 0:
                raise ValueError('%s is not an allowed base.' % (base,))

            from sage.rings.asymptotic.growth_group import MonomialGrowthGroup
            if len(factors) == 1:
                fp = factors[0].parent()
                if isinstance(fp, MonomialGrowthGroup) and repr(fp._var_).startswith('log('):
                    if factors[0]._raw_element_ == 1:
                        from sage.functions.log import log
                        new_elem = P(repr(fp._var_)[4:-1])
                        if base == 'e':
                            return new_elem
                        base_ring = new_elem.factor()[0].base_ring()
                        return new_elem ** log(base_ring(base))

            from sage.rings.asymptotic.growth_group import ExponentialGrowthGroup
            new_var = repr(self)
            if '*' in new_var or '^' in new_var:
                new_var = '(' + new_var + ')'

            def check_factor(factor):
                return new_var == repr(factor._var_) and \
                       isinstance(factor, ExponentialGrowthGroup)

            for cf in P.cartesian_factors():
                if hasattr(cf, 'cartesian_factors'):
                    for ccf in cf.cartesian_factors():
                        if check_factor(ccf):
                            return P(ccf(raw_element=ccf.base_ring()(base)))
                else:
                    if check_factor(cf):
                        return P(cf(raw_element=cf.base_ring()(base)))

            raise ValueError('Cannot construct %s^%s in %s.' % (base, new_var, P))


        def exp(self):
            r"""
            The exponential function of this element.

            INPUT:

            Nothing.

            OUTPUT:

            A growth element.

            EXAMPLES::

                sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                sage: G = GrowthGroup('x^ZZ * log(x)^ZZ * log(log(x))^ZZ')
                sage: x, = G.gens_monomial()
                sage: exp(log(x))
                x
                sage: exp(log(log(x)))
                log(x)

            ::

                sage: exp(x)
                Traceback (most recent call last):
                ...
                ValueError: Cannot construct e^x in Growth Group x^ZZ * log(x)^ZZ * log(log(x))^ZZ.
            """
            return self.rpow('e')


        def __invert__(self):
            r"""
            Return the multiplicative inverse of this cartesian product.

            OUTPUT:

            An growth element.

            .. NOTE::

                The result may live in a larger parent than we started with.

            TESTS::

                 sage: from sage.rings.asymptotic.growth_group import GrowthGroup
                 sage: G = GrowthGroup('ZZ^x * x^ZZ')
                 sage: g = G('2^x * x^3')
                 sage: (~g).parent()
                 Growth Group QQ^x * x^ZZ
            """
            new_element = tuple(~x for x in self.cartesian_factors())
            if all(n.parent() is x.parent()
                   for n, x in zip(new_element, self.cartesian_factors())):
                return self.parent()(new_element)
            else:
                from sage.categories.cartesian_product import cartesian_product
                new_parent = cartesian_product(
                    tuple(x.parent() for x in new_element))
                return new_parent(new_element)


    CartesianProduct = CartesianProductGrowthGroups


class UnivariateProduct(GenericProduct):
    def __init__(self, sets, category, **kwargs):
        r"""
        A cartesian product of growth groups with the same variables.

        TEST::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: type(GrowthGroup('x^ZZ * log(x)^ZZ'))  # indirect doctest
            <class 'sage.rings.asymptotic.growth_group_cartesian.UnivariateProduct_with_category'>
        """
        super(UnivariateProduct, self).__init__(
            sets, category, order='lex', **kwargs)


    CartesianProduct = CartesianProductGrowthGroups


class MultivariateProduct(GenericProduct):
    def __init__(self, sets, category, **kwargs):
        r"""
        A cartesian product of growth groups with the pairwise different variables.

        TEST::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: type(GrowthGroup('x^ZZ * y^ZZ'))  # indirect doctest
            <class 'sage.rings.asymptotic.growth_group_cartesian.MultivariateProduct_with_category'>
        """
        super(MultivariateProduct, self).__init__(
            sets, category, order='components', **kwargs)


    CartesianProduct = CartesianProductGrowthGroups
