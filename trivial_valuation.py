# -*- coding: utf-8 -*-
r"""
Trivial valuations

This module provides trivial valuations which mainly exist to be able to
provide an implementation for
:meth:`DiscretePseudoValuationSpace._an_element_`.

EXAMPLES::

    sage: from mac_lane import * # optional: standalone
    sage: v = DiscreteValuationSpace(QQ).an_element(); v
    Trivial valuation on Rational Field
    sage: v(1)
    0

.. NOTE:

Note that the tests in this module do not create instances of valuations
directly since this gives the wrong inheritance structure on the resulting
objects::

    sage: v = TrivialDiscretePseudoValuation(QQ)
    sage: v._test_category()
    Traceback (most recent call last):
    ...
    AssertionError: False is not true

Instead, the valuations need to be created through the
``__make_element_class__`` of the containing space::

    sage: H = DiscretePseudoValuationSpace(QQ)
    sage: v = H.__make_element_class__(TrivialDiscretePseudoValuation)(QQ)
    sage: v._test_category()

Using ``an_element`` provides a much more readable shortcut for this::

    sage: v = H.an_element()
    sage: v._test_category()

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

from valuation import DiscretePseudoValuation

class TrivialDiscretePseudoValuation_base(DiscretePseudoValuation):
    r"""
    Base class for code shared by trivial valuations.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: v = DiscretePseudoValuationSpace(ZZ).an_element(); v # indirect doctest
        Trivial pseudo-valuation on Integer Ring

    TESTS::

        sage: TestSuite(v).run()

    """
    def __reduce__(self):
        r"""
        Return pickling information for this valuation.

        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscretePseudoValuationSpace(ZZ).an_element()
            sage: loads(dumps(v)) == v # indirect doctest
            True

        """
        return self.__class__, (self.domain(), )

    def _eq_(self, other):
        r"""
        Return whether this valuation is indistinguishable to the
        pseudo-valuation ``other``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscretePseudoValuationSpace(ZZ).an_element()
            sage: w = DiscretePseudoValuationSpace(ZZ).an_element()
            sage: v == w
            True

        """
        # other lives in the same space of valuations, say the space of
        # discrete valuations on the integers; therefore, it suffices to
        # compare the types since there are no paremeters to trivial valuations
        return type(self) == type(other)

    def uniformizer(self):
        r"""
        Return a uniformizing element for this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscretePseudoValuationSpace(QQ).an_element()
            sage: v.uniformizer()
            Traceback (most recent call last):
            ...
            ValueError: Trivial valuations do not define a uniformizing element

        """
        raise ValueError("Trivial valuations do not define a uniformizing element")

    def is_trivial(self):
        r"""
        Return whether this valuation is trivial.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscretePseudoValuationSpace(QQ).an_element()
            sage: v.is_trivial()
            True

        """
        return True

class TrivialDiscretePseudoValuation(TrivialDiscretePseudoValuation_base):
    r"""
    The trivial pseudo-valuation that is `\infty` everywhere.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: v = DiscretePseudoValuationSpace(QQ).an_element(); v
        Trivial pseudo-valuation on Rational Field

    TESTS::

        sage: TestSuite(v).run()

    """
    def __init__(self, domain):
        r"""
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscretePseudoValuationSpace(QQ).an_element()
            sage: isinstance(v, TrivialDiscretePseudoValuation)
            True
    
        """
        from .valuation_space import DiscretePseudoValuationSpace
        TrivialDiscretePseudoValuation_base.__init__(self, DiscretePseudoValuationSpace(domain))

    def is_discrete_valuation(self):
        r"""
        Return whether this is a discrete valuation.

        EXAMPLES:

        Returns ``False`` since this is only a pseudo-valuation::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscretePseudoValuationSpace(QQ).an_element()
            sage: v.is_discrete_valuation()
            False

        """
        return False

    def _call_(self, x):
        r"""
        Evaluate this valuation at ``x``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscretePseudoValuationSpace(QQ).an_element()
            sage: v(0)
            +Infinity
            sage: v(1)
            +Infinity

        """
        from sage.rings.all import infinity
        return infinity

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscretePseudoValuationSpace(QQ).an_element(); v # indirect doctest
            Trivial pseudo-valuation on Rational Field

        """
        return "Trivial pseudo-valuation on %r"%(self.domain(),)

    def value_group(self):
        r"""
        Return the value group of this valuation.

        EXAMPLES:

        A trivial discrete pseudo-valuation has no value group::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscretePseudoValuationSpace(QQ).an_element()
            sage: v.value_group()
            Traceback (most recent call last):
            ...
            ValueError: The trivial pseudo-valuation that is infinity everywhere does not have a value group.

        """
        raise ValueError("The trivial pseudo-valuation that is infinity everywhere does not have a value group.")

class TrivialDiscreteValuation(TrivialDiscretePseudoValuation_base):
    r"""
    The trivial valuation that is zero on non-zero elements.

    EXAMPLES::

        sage: from mac_lane import * # optional: standalone
        sage: v = DiscreteValuationSpace(QQ).an_element(); v # indirect doctest
        Trivial valuation on Rational Field

    TESTS::

        sage: TestSuite(v).run()

    """
    def __init__(self, domain):
        r"""
        TESTS::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscreteValuationSpace(QQ).an_element()
            sage: isinstance(v, TrivialDiscreteValuation)
            True
    
        """
        from .valuation_space import DiscreteValuationSpace
        TrivialDiscretePseudoValuation_base.__init__(self, DiscreteValuationSpace(domain))

    def _call_(self, x):
        r"""
        Evaluate this valuation at ``x``.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscreteValuationSpace(QQ).an_element()
            sage: v(0)
            +Infinity
            sage: v(1)
            0

        """
        from sage.rings.all import infinity
        return infinity if x == 0 else self.codomain().zero()

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscreteValuationSpace(QQ).an_element(); v # indirect doctest
            Trivial valuation on Rational Field

        """
        return "Trivial valuation on %r"%(self.domain(),)

    def value_group(self):
        r"""
        Return the value group of this valuation.

        EXAMPLES:

        A trivial discrete valuation has a trivial value group::

            sage: from mac_lane import * # optional: standalone
            sage: v = DiscreteValuationSpace(QQ).an_element()
            sage: v.value_group()
            Trivial Additive Abelian Group

        """
        from .value_group import DiscreteValueGroup
        return DiscreteValueGroup(0)
