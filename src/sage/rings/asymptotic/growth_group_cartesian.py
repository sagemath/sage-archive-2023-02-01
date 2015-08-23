r"""
Growth Groups as Cartesian Products

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


def merge_overlapping(A, B, key=None):
    r"""
    Merge the two overlapping tuples/lists.

    TESTS::

        sage: from sage.rings.asymptotic.growth_group_cartesian import merge_overlapping
        sage: def f(L, s):
        ....:     return list((ell, s) for ell in L)
        sage: key = lambda k: k[0]
        sage: merge_overlapping(f([0..3], 'a'), f([5..7], 'b'), key)
        Traceback (most recent call last):
        ...
        ValueError: Input does not have an overlap.
        sage: merge_overlapping(f([0..2], 'a'), f([4..7], 'b'), key)
        Traceback (most recent call last):
        ...
        ValueError: Input does not have an overlap.
        sage: merge_overlapping(f([4..7], 'a'), f([0..2], 'b'), key)
        Traceback (most recent call last):
        ...
        ValueError: Input does not have an overlap.
        sage: merge_overlapping(f([0..3], 'a'), f([3..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a'), (4, 'b')],
         [(0, 'a'), (1, 'a'), (2, 'a'), (3, 'b'), (4, 'b')])
        sage: merge_overlapping(f([3..4], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'b'), (2, 'b'), (3, 'a'), (4, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b'), (4, 'a')])
        sage: merge_overlapping(f([0..1], 'a'), f([0..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'b'), (3, 'b'), (4, 'b')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b'), (4, 'b')])
        sage: merge_overlapping(f([0..3], 'a'), f([0..1], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'a'), (3, 'a')])
        sage: merge_overlapping(f([0..3], 'a'), f([1..3], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'a'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: merge_overlapping(f([1..3], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: merge_overlapping(f([0..6], 'a'), f([3..4], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a'), (4, 'a'), (5, 'a'), (6, 'a')],
         [(0, 'a'), (1, 'a'), (2, 'a'), (3, 'b'), (4, 'b'), (5, 'a'), (6, 'a')])
        sage: merge_overlapping(f([0..3], 'a'), f([1..2], 'b'), key)
        ([(0, 'a'), (1, 'a'), (2, 'a'), (3, 'a')],
         [(0, 'a'), (1, 'b'), (2, 'b'), (3, 'a')])
        sage: merge_overlapping(f([1..2], 'a'), f([0..3], 'b'), key)
        ([(0, 'b'), (1, 'a'), (2, 'a'), (3, 'b')],
         [(0, 'b'), (1, 'b'), (2, 'b'), (3, 'b')])
        sage: merge_overlapping(f([1..3], 'a'), f([1..3], 'b'), key)
        ([(1, 'a'), (2, 'a'), (3, 'a')],
         [(1, 'b'), (2, 'b'), (3, 'b')])
    """
    if key is None:
        key = lambda k: k

    def find_overlapping_index(A, B):
        if len(B) > len(A) - 2:
            raise StopIteration
        matches = iter(i for i in xrange(1, len(A) - len(B))
                       if all(key(a) == key(b) for a, b in zip(A[i:i+len(B)], B)))
        return next(matches)

    def find_mergedoverlapping_index(A, B):
        """
        Return in index i where to merge two overlapping tuples/lists ``A`` and ``B``.

        Then ``A + B[i:]`` or ``A[:-i] + B`` are the merged tuples/lists.

        Adapted from http://stackoverflow.com/a/30056066/1052778.
        """
        matches = iter(i for i in xrange(min(len(A), len(B)), 0, -1)
                       if all(key(a) == key(b) for a, b in zip(A[-i:], B[:i])))
        return next(matches, 0)

    i = find_mergedoverlapping_index(A, B)
    if i > 0:
        return A + B[i:], A[:-i] + B

    i = find_mergedoverlapping_index(B, A)
    if i > 0:
        return B[:-i] + A, B + A[i:]

    try:
        i = find_overlapping_index(A, B)
    except StopIteration:
        pass
    else:
        return A, A[:i] + B + A[i+len(B):]

    try:
        i = find_overlapping_index(B, A)
    except StopIteration:
        pass
    else:
        return B[:i] + A + B[i+len(A):], B

    raise ValueError('Input does not have an overlap.')


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
        x^(1/2) * log(x)

    ::

        sage: Px = agg.MonomialGrowthGroup(QQ, 'x')
        sage: Lx = agg.MonomialGrowthGroup(ZZ, 'log(x)')
        sage: Cx = cartesian_product([Px, Lx], order='lex')
        sage: Py = agg.MonomialGrowthGroup(QQ, 'y')
        sage: C = cartesian_product([Cx, Py], order='components'); C
        Growth Group x^QQ * log(x)^ZZ * y^QQ
        sage: C.an_element()
        x^(1/2) * log(x) * y^(1/2)

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
            sage: G(x^-3 * y^2)
            x^(-3) * y^2
            sage: G(x^4), G(y^2)
            (x^4, y^2)
            sage: G(1)
            1

        Even more complex expressions can be parsed::

            sage: G_log(x^42 * log(x)^-42 * y^42)
            x^42 * log(x)^(-42) * y^42

        TESTS::

            sage: G = GrowthGroup('x^ZZ * y^ZZ')
            sage: G('x'), G('y')
            (x, y)

        ::

            sage: G_log(log(x))
            log(x)
        """
        if data == 1:
            return self.one()

        if data is None:
            raise ValueError('%s cannot be converted.' % (data,))

        if isinstance(data, list):
            try:
                obj = super(GenericProduct,
                            self)._element_constructor_(data)
                return obj
            except ValueError:
                return self.prod(self(elem) for elem in data)

        if hasattr(data, 'parent'):
            P = data.parent()
            if P is self:
                return data

            elif P is sage.symbolic.ring.SR:
                import operator
                from sage.symbolic.operators import mul_vararg
                op = data.operator()
                if op == operator.pow or data.is_symbol() \
                        or isinstance(op, sage.functions.log.Function_log):
                    return self(self._convert_to_factor_(data))
                elif op == mul_vararg:
                    return self(data.operands())
            # room for other parents (e.g. polynomial ring et al.)

        # try to convert the input to one of the factors
        data_conv = self._convert_to_factor_(data)
        if data_conv is not None:
            factors = self.cartesian_factors()
            return self([data_conv if factor == data_conv.parent() else 1 for
                         factor in factors])

        # final attempt: try parsing the representation string
        str_lst = str(data).replace(' ', '').split('*')
        return self(str_lst)


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


    def _convert_to_factor_(self, data):
        r"""
        Helper method. Try to convert some input ``data`` to an
        element of one of the cartesian factors of this product.

        INPUT:

        - ``data`` -- some input to be converted.

        OUTPUT:

        An element of an cartesian factor of this product,
        or ``None``.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ * log(x)^QQ * y^QQ')
            sage: e1 = G._convert_to_factor_(x^2)
            sage: (e1, e1.parent())
            (x^2, Growth Group x^ZZ * log(x)^QQ)
            sage: G._convert_to_factor_('asdf') is None
            True
        """
        for factor in self.cartesian_factors():
            try:
                if hasattr(factor, '_convert_to_factor_'):
                    return factor(factor._convert_to_factor_(data))
                return factor(data)
            except (ValueError, TypeError):
                continue


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


        #C = other.construction()[0]
        #if isinstance(C, PolynomialFunctor):
        #    return pushout(self, CartesianProductPolys((other,)))



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
        """
        return sum(iter(factor.gens_monomial()
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
                'x^(1/2) * log(x)'
            """
            s = ' * '.join(repr(v) for v in self.value if not v.is_one())
            if s == '':
                return '1'
            return s


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
