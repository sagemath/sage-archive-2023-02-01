# -*- coding: utf-8 -*-
r"""
Spaces of valuations

This module provides spaces of exponential pseudo-valuations on integral
domains. It currently, only provides support for such valuations if they are
discrete, i.e., their image is a discrete additive subgroup of the rational
numbers extended by `\infty`.

EXAMPLES::

    sage: from mac_lane import * # optional: standalone
    sage: pAdicValuation(QQ, 2).parent()
    Discrete valuations on Rational Field

AUTHORS:

- Julian Rueth (2016-10-14): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Fix doctests so they work in standalone mode (when invoked with sage -t, they run within the mac_lane/ directory)
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.getcwd())
    sys.path.append(os.path.dirname(os.getcwd()))

from sage.categories.homset import Homset
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.abstract_method import abstract_method
from sage.structure.unique_representation import UniqueRepresentation

class DiscretePseudoValuationSpace(UniqueRepresentation, Homset):
    r"""
    The space of discrete pseudo-valuations on ``domain``.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: H = DiscretePseudoValuationSpace(QQ)
        sage: pAdicValuation(QQ, 2) in H
        True

    TESTS::

        sage: TestSuite(H).run()

    """
    def __init__(self, domain):
        r"""
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: isinstance(pAdicValuation(QQ, 2).parent(), DiscretePseudoValuationSpace)
            True

        """
        from sage.categories.domains import Domains
        if domain not in Domains():
            raise ValueError("domain must be an integral domain")

        from value_group import DiscreteValuationCodomain
        # A valuation is a map from an additive semigroup to an additive semigroup, however, it
        # does not preserve that structure. It is therefore only a morphism in the category of sets.
        from sage.categories.all import Sets
        Homset.__init__(self, domain, DiscreteValuationCodomain(), category = Sets())

    @lazy_attribute
    def _abstract_element_class(self):
        r"""
        Return an abstract base class for all valuations in this space.

        This is used to extend every valuation with a number of generic methods
        that are independent of implementation details.

        Usually, extensions of this kind would be done by implementing an
        appropriate class ``MorphismMethods`` in the category of this homset.
        However, there is no category whose arrows are the valuations, so we
        need to move this magic down to the level of the actual homset.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: isinstance(pAdicValuation(QQ, 2), DiscretePseudoValuationSpace.ElementMethods) # indirect doctest
            True

        """
        class_name = "%s._abstract_element_class"%self.__class__.__name__
        from sage.structure.dynamic_class import dynamic_class
        return dynamic_class(class_name, (super(DiscretePseudoValuationSpace,self)._abstract_element_class, self.__class__.ElementMethods))

    def _an_element_(self):
        r"""
        Return a trivial valuation in this space.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: DiscretePseudoValuationSpace(QQ).an_element() # indirect doctest
            Trivial pseudo-valuation on Rational Field

        """
        from trivial_valuation import TrivialPseudoValuation
        return TrivialPseudoValuation(self.domain())

    def _repr_(self):
        r"""
        Return a printable representation of this space.
        
        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: DiscretePseudoValuationSpace(QQ) # indirect doctest
            Discrete pseudo-valuations on Rational Field

        """
        return "Discrete pseudo-valuations on %r"%(self.domain(),)

    def __contains__(self, x):
        r"""
        Return whether ``x`` is a valuation in this space.

        EXAMPLES:

            sage: from mac_lane import * # optional: standalone
            sage: H = DiscretePseudoValuationSpace(QQ)
            sage: H.an_element() in H
            True

        Elements of spaces which embed into this spaces are correctly handled::

            sage: pAdicValuation(QQ, 2) in H
            True

        """
        # override the logic from Homset with the original implementation for Parent
        # which entirely relies on a proper implementation of
        # _element_constructor_ and coercion maps
        from sage.structure.parent import Parent
        return Parent.__contains__(self, x)

    def __call__(self, x):
        r"""
        Create an element in this space from ``x``.

        EXAMPLES:

            sage: from mac_lane import * # optional: standalone
            sage: H = DiscretePseudoValuationSpace(QQ)
            sage: H(pAdicValuation(QQ, 2))
            2-adic valuation

        """
        # override the logic from Homset with the original implementation for Parent
        # which entirely relies on a proper implementation of
        # _element_constructor_ and coercion maps
        from sage.structure.parent import Parent
        return Parent.__call__(self, x)

    def _element_constructor_(self, x):
        r"""
        Create an element in this space from ``x``,

        EXAMPLES:

        We try to convert valuations defined on different domains by changing
        their base ring::

            sage: from mac_lane import * # optional: standalone
            sage: Z = DiscretePseudoValuationSpace(ZZ)
            sage: Q = DiscretePseudoValuationSpace(QQ)
            sage: v = pAdicValuation(ZZ, 2)
            sage: v in Q
            False
            sage: Q(v) in Q
            True
            sage: Q(v) in Z
            False
            sage: Z(Q(v)) in Z
            True

        We support coercions and conversions, even though they are not
        implemented here::

            sage: Z(v)
            2-adic valuation

        """
        if isinstance(x.parent(), DiscretePseudoValuationSpace):
            if x.domain() is not self.domain():
                try:
                    return self(x.change_ring(self.domain()))
                except NotImplementedError:
                    pass
            else:
                # If this is an element of a discrete pseudo-valuation space over the same domain,
                # then we treat it as an element of this space (see __contains__), i.e., we do not
                # actually change x.parent() to self here if it is, e.g., a discrete valuation space.
                # This might be surprising but is how facades work for example.
                return x
        raise ValueError("element can not be converted into the space of %r"%(self,))

    class ElementMethods:
        r"""
        Provides methods for discrete pseudo-valuations that are added
        automatically to valuations in this space.

        EXAMPLES:

        Here is an example of a method that is automagically added to a
        discrete valuation::

            sage: from mac_lane import * # optional: standalone
            sage: H = DiscretePseudoValuationSpace(QQ)
            sage: pAdicValuation(QQ, 2).is_discrete_pseudo_valuation() # indirect doctest
            True

        The methods will be provided even if the concrete types is not created
        with :meth:`__make_element_class__`::

            sage: from valuation import DiscretePseudoValuation
            sage: m = DiscretePseudoValuation(H)
            sage: m.parent() is H
            True
            sage: m.is_discrete_pseudo_valuation()
            True

        However, the category framework advises you to use inheritance::

            sage: m._test_category()
            Traceback (most recent call last):
            ...
            AssertionError: False is not true

        Using :meth:`__make_element_class__`, makes your concrete valuation
        inherit from this class::

            sage: m = H.__make_element_class__(DiscretePseudoValuation)(H)
            sage: m._test_category()

        """
        def is_discrete_pseudo_valuation(self):
            r"""
            Return whether this valuation is a discrete pseudo-valuation.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: pAdicValuation(QQ, 2).is_discrete_pseudo_valuation()
                True

            """
            return True

        @abstract_method
        def is_discrete_valuation(self):
            r"""
            Return whether this valuation is a discrete valuation, i.e.,
            whether it is a :meth:`is_discrete_pseudo_valuation` that only
            sends zero to `\infty`.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: pAdicValuation(QQ, 2).is_discrete_valuation()
                True
            
            """

        def is_trivial(self):
            r"""
            Return whether this valuation is trivial, i.e., whether it is
            constant `\infty` or constant zero for everything but the zero
            element.

            Subclasses need to override this method if they do not implement
            :meth:`uniformizer`.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: pAdicValuation(QQ, 7).is_trivial()
                False

            """
            from sage.rings.all import infinity
            if self(self.domain().one()) == infinity:
                # the constant infinity
                return True
            if self(self.uniformizer()) != 0:
                # not constant on the non-zero elements
                return False
            return True

        @abstract_method
        def uniformizer(self):
            r"""
            Return an element in :meth:`domain` which has positive valuation
            and generates the value group of this valuation.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: pAdicValuation(QQ, 11).uniformizer()
                11

            Trivial valuations have no uniformizer::

                sage: v = DiscretePseudoValuationSpace(QQ).an_element()
                sage: v.is_trivial()
                True
                sage: v.uniformizer()
                Traceback (most recent call last):
                ...
                ValueError: Trivial valuations do not define a uniformizing element
                
            """

        def value_group(self, **options):
            r"""
            Return the value group of this discrete pseudo-valuation, a
            discrete additive subgroup of the rational numbers.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: pAdicValuation(QQ, 2).value_group()
                Additive Abelian Group generated by 1

            A pseudo-valuation that is `\infty` everywhere, does not have a
            value group::

                sage: v = DiscretePseudoValuationSpace(QQ).an_element()
                sage: v.value_group()
                Traceback (most recent call last):
                ...
                ValueError: The trivial pseudo-valuation that is infinity everywhere does not have a value group.

            """
            from value_group import DiscreteValueGroup
            return DiscreteValueGroup(self(self.uniformizer()))

        def shift(self, x, s):
            r"""
            Return a modified version of ``x`` whose valuation is increased by
            ``s``.

            The element returned has essentially the same reduction, i.e., if
            ``x`` has valuation `v`, then the reduction of ``x`` in the residue
            ring of elements of valuation `\ge v` module elements of valuation
            `> v` is naturally the same as the reduction of ``shift(x, s)`` in
            the correspoding residue ring of elements of valuation `\ge v + s`.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(ZZ, 2)
                sage: v.shift(1, 10)
                1024
                sage: v.shift(1, -10)
                Traceback (most recent call last):
                ...
                TypeError: no conversion of this rational to integer

            """
            x = self.domain().coerce(x)
            from sage.rings.all import QQ, ZZ
            s = QQ.coerce(s)
            if s == 0:
                return x
            if s not in self.value_group():
                raise ValueError("s must be in the value group of this valuation")
            return self.domain()(x * (self.uniformizer() ** ZZ(s/self.value_group().gen())))

        @abstract_method
        def residue_ring(self):
            r"""
            Return the residue ring of this valuation, i.e., the elements of
            non-negative valuation module the elements of positive valuation.

            This is identical to :meth:`residue_field` when a residue field
            exists.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: pAdicValuation(QQ, 2).residue_ring()
                Finite Field of size 2
                sage: TrivialValuation(QQ).residue_ring()
                Rational Field

            Note that a residue ring always exists, even when a residue field
            may not::

                sage: TrivialPseudoValuation(QQ).residue_ring()
                Quotient of Rational Field by the ideal (1)
                sage: TrivialValuation(ZZ).residue_ring()
                Integer Ring
                sage: GaussValuation(ZZ['x'], pAdicValuation(ZZ, 2)).residue_ring()
                Univariate Polynomial Ring in x over Finite Field of size 2 (using NTL)


            """

        @abstract_method
        def reduce(self, x):
            r"""
            Return the image of ``x`` in the :meth:`residue_ring` of this
            valuation.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 2)
                sage: v.reduce(2)
                0
                sage: v.reduce(1)
                1
                sage: v.reduce(1/3)
                1
                sage: v.reduce(1/2)
                Traceback (most recent call last):
                ...
                ValueError: reduction is only defined for elements of non-negative valuation

            """

        @abstract_method
        def lift(self, X):
            r"""
            Return a lift of ``X`` in :meth:`domain` which reduces down to
            ``X`` again via :meth:`reduce`.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 2)
                sage: v.lift(v.residue_ring().one())
                1

            """

        def extension(self, ring):
            r"""
            Return the unique extension of this valuation to ``ring``.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(ZZ, 2)
                sage: w = v.extension(QQ)
                sage: w.domain()
                Rational Field

            """
            extensions = self.extensions(ring)
            assert(len(extensions))
            if len(extensions) > 1:
                raise ValueError("there is no unique extension of %r from %r to %r"%(self, self.domain(), ring))
            return extensions[0]

        def extensions(self, ring):
            r"""
            Return the extensions of this valuation to ``ring``.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(ZZ, 2)
                sage: v.extensions(QQ)
                [2-adic valuation]

            """
            if ring is self.domain():
                return [self]
            raise NotImplementedError("extending %r from %r to %r not implemented"%(self, self.domain(), ring))

        def restriction(self, ring):
            r"""
            Return the restriction of this valuation to ``ring``.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 2)
                sage: w = v.restriction(ZZ)
                sage: w.domain()
                Integer Ring

            """
            if ring is self.domain():
                return self
            raise NotImplementedError("restricting %r from %r to %r not implemented"%(self, self.domain(), ring))

        def change_ring(self, ring):
            r"""
            Return this valuation over ``ring``.

            Unlike :meth:`extension` or meth:`reduction`, this might not be
            completely sane mathematically. It is essentially a conversion of
            this valuation into another space of valuations.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 3)
                sage: v.change_ring(ZZ)
                3-adic valuation

            """
            if ring is self.domain():
                return self
            raise NotImplementedError("changing %r from %r to %r not implemented"%(self, self.domain(), ring))

        def _test_add(self, **options):
            r"""
            Check that the (strict) triangle equality is satisfied for the
            valuation of this ring.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(ZZ, 3)
                sage: v._test_add()

            """
            tester = self._tester(**options)
            S = tester.some_elements(self.domain().some_elements())
            from sage.categories.cartesian_product import cartesian_product
            for x,y in tester.some_elements(cartesian_product([S,S])):
                tester.assertGreaterEqual(self(x+y),min(self(x),self(y)))
                if self(x) != self(y):
                    tester.assertEqual(self(x+y),min(self(x),self(y)))

        def _test_infinite_zero(self, **options):
            r"""
            Check that zero is sent to infinity.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_infinite_zero()

            """
            tester = self._tester(**options)
            from sage.rings.all import infinity
            tester.assertEqual(self(self.domain().zero()), infinity)

        def _test_mul(self, **options):
            r"""
            Check that multiplication translates to addition of valuations.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_mul()

            """
            tester = self._tester(**options)
            S = list(self.domain().some_elements())
            from sage.categories.cartesian_product import cartesian_product
            for x,y in tester.some_elements(cartesian_product([S,S])):
                tester.assertEqual(self(x*y),self(x)+self(y))

        def _test_no_infinite_units(self, **options):
            r"""
            Checks that no units are sent to infinity.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_no_infinite_units()

            As multiplication translates to addition, pseudo-valuations which
            send a unit to infinity are necessarily trivial::

                sage: v = DiscretePseudoValuationSpace(QQ).an_element()
                sage: v(1)
                +Infinity
                sage: v.is_trivial()
                True

            """
            if not self.is_discrete_valuation() and self.is_trivial():
                return

            from sage.rings.all import infinity
            tester = self._tester(**options)
            for x in tester.some_elements(self.domain().some_elements()):
                if self(x) == infinity:
                    tester.assertFalse(x.is_unit())

        def _test_value_group(self, **options):
            r"""
            Check correctness of the value group.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_value_group()

            """
            from sage.rings.all import infinity
            tester = self._tester(**options)
            # check consistency of trivial valuations first
            if self.is_trivial():
                if self(self.domain().one()) == infinity:
                    # a trivial pseudo-valuation that sends everything to infinity
                    with tester.assertRaises(ValueError):
                        self.value_group()
                    return

            # check that all valuations are in the value group
            for x in tester.some_elements(self.domain().some_elements()):
                if self(x) != infinity:
                    tester.assertIn(self(x), self.value_group())

            if not self.is_trivial():
                # check that the uniformizer generates the value group
                tester.assertEqual(self.value_group().gen(), self(self.uniformizer()))

        def _test_shift(self, **options):
            r"""
            Check correctness of shifts.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_shift()

            """
            tester = self._tester(**options)

            if self.is_trivial():
                # trivial valuations can not perform non-trivial shifts
                return

            S = tester.some_elements(self.domain().some_elements())
            V = tester.some_elements(self.value_group().some_elements())
            from sage.categories.cartesian_product import cartesian_product
            for x, n in tester.some_elements(cartesian_product([S,V])):
                v = self(x)
                from sage.categories.fields import Fields
                if n < 0 and self.domain() not in Fields():
                    # note that shifting might not be possible in this case even if -n > v
                    continue
                y = self.shift(x, n)
                tester.assertIs(y.parent(), self.domain())
                tester.assertEqual(self(y), v + n)
                # shifts preserve reductions
                z = self.shift(y, -n)
                tester.assertEqual(self(z), v)
                if v >= 0:
                    tester.assertEqual(self.reduce(z), self.reduce(x))

        def _test_residue_ring(self, **options):
            r"""
            Check the correctness of residue fields.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_residue_ring()

            """
            tester = self._tester(**options)

            r = self.residue_ring()
            if r.zero() == r.one():
                # residue ring is the zero rng
                tester.assertGreater(self(1), 0)
                return

            c = self.residue_ring().characteristic()
            if c != 0:
                tester.assertGreater(self(c), 0)

        def _test_reduce(self, **options):
            r"""
            Check the correctness of reductions.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_reduce()

            """
            tester = self._tester(**options)

            for x in tester.some_elements(self.domain().some_elements()):
                if self(x) < 0:
                    with tester.assertRaises(ValueError):
                        self.reduce(x)
                    continue
                if self(x) == 0:
                    y = self.reduce(x)
                    tester.assertIn(y, self.residue_ring())
                    tester.assertNotEqual(y, 0)
                    if x.is_unit() and ~x in self.domain():
                        tester.assertTrue(y.is_unit())
                        tester.assertIn(~y, self.residue_ring())
                        tester.assertEqual(~y, self.reduce(self.domain()(~x)))
                if self(x) > 0:
                    tester.assertEqual(self.reduce(x), 0)

        def _test_lift(self, **options):
            r"""
            Check the correctness of lifts.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_lift()

            """
            tester = self._tester(**options)

            for X in tester.some_elements(self.residue_ring().some_elements()):
                x = self.lift(X)
                y = self.reduce(x)
                tester.assertEqual(X, y)
                if X != 0:
                    tester.assertEqual(self(x), 0)

        def _test_restriction(self, **options):
            r"""
            Check the correctness of reductions.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_restriction()

            """
            tester = self._tester(**options)

            tester.assertEqual(self.restriction(self.domain()), self)

        def _test_extension(self, **options):
            r"""
            Check the correctness of extensions.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_extension()

            """
            tester = self._tester(**options)

            tester.assertEqual(self.extension(self.domain()), self)
            tester.assertEqual(self.extensions(self.domain()), [self])

        def _test_change_ring(self, **options):
            r"""
            Check the correctness of :meth:`change_ring`.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_change_ring()

            """
            tester = self._tester(**options)

            tester.assertEqual(self.change_ring(self.domain()), self)


class DiscreteValuationSpace(DiscretePseudoValuationSpace):
    r"""
    The space of discrete valuations on ``domain``.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: H = DiscreteValuationSpace(QQ)
        sage: pAdicValuation(QQ, 2) in H
        True

    TESTS::

        sage: TestSuite(H).run()

    """
    def __init__(self, domain):
        r"""
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: isinstance(pAdicValuation(QQ, 2).parent(), DiscreteValuationSpace)
            True

        """
        DiscretePseudoValuationSpace.__init__(self, domain)

    def _element_constructor_(self, x):
        r"""
        Create an element in this space from ``x``.

        TESTS:

        Discrete pseudo-valuations that are actually valuations are contained
        in this space::

            sage: from mac_lane import * # optional: standalone
            sage: H = DiscretePseudoValuationSpace(QQ)
            sage: G = DiscreteValuationSpace(QQ)
            sage: v = pAdicValuation(QQ, 2)
            sage: v._set_parent(H)
            sage: v in G
            True
        
        Discrete pseudo-valuations that are discrete valuations can be
        converted into this space::

            sage: H = DiscreteValuationSpace(ZZ)
            sage: v = pAdicValuation(ZZ, 2)
            sage: v._set_parent(H)
            sage: v in G
            False
            sage: G(v) in G
            True

        """
        # We accept any discrete pseudo-valuation that claims to be a discrete valuation
        if isinstance(x.parent(), DiscretePseudoValuationSpace) and x.is_discrete_valuation():
            if x.domain() is not self.domain():
                # after we base-changed them into our domain
                return DiscretePseudoValuationSpace(self.domain())(x)
            # If this is a valuation in a discrete pseudo-valuation space over the same domain,
            # then we treat it as an element of this space (see __contains__), i.e., we do not
            # actually change x.parent() to self here if it is, e.g., a
            # discrete pseudo-valuation space.  This might be surprising but is
            # also how facades work for example.
            return x
        raise ValueError("element does not convert to a discrete valuation in %r"%(self,))

    def _repr_(self):
        r"""
        Return a printable representation of this space.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: DiscreteValuationSpace(QQ) # indirect doctest
            Discrete valuations on Rational Field

        """
        return "Discrete valuations on %r"%(self.domain(),)

    def _an_element_(self):
        r"""
        Return the trivial valuation in this space.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: DiscreteValuationSpace(QQ).an_element() # indirect doctest
            Trivial valuation on Rational Field

        """
        from trivial_valuation import TrivialValuation
        return TrivialValuation(self.domain())

    class ElementMethods(DiscretePseudoValuationSpace.ElementMethods):
        r"""
        Provides methods for discrete valuations that are added
        automatically to valuations in this space.

        EXAMPLES:

        Here is an example of a method that is automagically added to a
        discrete valuation::

            sage: from mac_lane import * # optional: standalone
            sage: pAdicValuation(QQ, 2).is_discrete_valuation() # indirect doctest
            True

        """
        def is_discrete_valuation(self):
            r"""
            Return whether this valuation is a discrete valuation.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: pAdicValuation(QQ, 2).is_discrete_valuation()
                True

            """
            return True

        def residue_field(self):
            r"""
            Return the residue field of this valuation, i.e., the field of
            fractions of the :meth:`residue_ring`, the elements of non-negative
            valuation module the elements of positive valuation.

            EXAMPLES::

                sage: from mac_lane import * # optional: standalone
                sage: pAdicValuation(QQ, 2).residue_field()
                Finite Field of size 2
                sage: TrivialValuation(QQ).residue_field()
                Rational Field

                sage: TrivialValuation(ZZ).residue_field()
                Rational Field
                sage: GaussValuation(ZZ['x'], pAdicValuation(ZZ, 2)).residue_field()
                Rational function field in x over Finite Field of size 2

            """
            ret = self.residue_ring()
            from sage.categories.fields import Fields
            if ret in Fields():
                return ret
            from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
            if is_PolynomialRing(ret):
                from sage.rings.function_field.all import FunctionField
                return FunctionField(ret.base_ring().fraction_field(), names=(ret.variable_name(),))
            return ret.fraction_field()

        def _test_no_infinite_nonzero(self, **options):
            r"""
            Check that only zero is sent to infinity.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_no_infinite_nonzero()

            """
            from sage.rings.all import infinity
            tester = self._tester(**options)
            for x in tester.some_elements(self.domain().some_elements()):
                if self(x) is infinity:
                    tester.assertEqual(x, 0)

        def _test_residue_field(self, **options):
            r"""
            Check the correctness of residue fields.

            TESTS::

                sage: from mac_lane import * # optional: standalone
                sage: v = pAdicValuation(QQ, 5)
                sage: v._test_residue_field()

            """
            tester = self._tester(**options)
            try:
                k = self.residue_field()
            except ValueError:
                from sage.categories.fields import Fields
                # a discrete valuation on a field has a residue field
                tester.assertFalse(self.domain() in Fields())
                return

            c = self.residue_field().characteristic()
            if c != 0:
                tester.assertGreater(self(c), 0)
