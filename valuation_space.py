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
        from trivial_valuation import TrivialDiscretePseudoValuation
        return self.__make_element_class__(TrivialDiscretePseudoValuation)(self.domain())

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
            return self(x.change_ring(self.domain()))
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
            return DiscreteValueGroup(self(self.uniformizer))

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

        # discrete valuations are discrete pseudo-valuations
        self.register_embedding(self.hom(lambda x:x, DiscretePseudoValuationSpace(domain)))

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
        # We try to convert discrete pseudo-valuations that claims to be a
        # valuation into this space
        if isinstance(x.parent(), DiscretePseudoValuationSpace) and x.is_discrete_valuation():
            # after we base-changed them into our domain
            return DiscretePseudoValuationSpace(self.domain())(x)
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
        from trivial_valuation import TrivialDiscreteValuation
        return self.__make_element_class__(TrivialDiscreteValuation)(self.domain())

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
