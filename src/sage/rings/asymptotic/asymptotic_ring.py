r"""
Asymptotic Ring

This module provides a ring (called :class:`AsymptoticRing`) for
computations with :wikipedia:`asymptotic expansions <Asymptotic_expansion>`.


.. _asymptotic_ring_definition:

(Informal) Definition
=====================

An asymptotic expansion is a sum such as

.. MATH::

    5z^3 + 4z^2 + O(z)

as `z \to \infty` or

.. MATH::

    3x^{42}y^2 + 7x^3y^3 + O(x^2) + O(y)

as `x` and `y` tend to `\infty`. It is a truncated series (after a
finite number of terms), which approximates a function.

The summands of the asymptotic expansions are partially ordered. In
this module these summands are the following:

- Exact terms `c\cdot g` with a coefficient `c` and an element `g` of
  a growth group (:ref:`see below <asymptotic_ring_growth>`).

- `O`-terms `O(g)` (see :wikipedia:`Big O notation <Big_O_notation>`;
  also called *Bachmann--Landau notation*) for a growth group
  element `g` (:ref:`again see below <asymptotic_ring_growth>`).

See
:wikipedia:`the Wikipedia article on asymptotic expansions <Asymptotic_expansion>`
for more details.
Further examples of such elements can be found :ref:`here <asymptotic_ring_intro>`.


.. _asymptotic_ring_growth:

Growth Groups and Elements
--------------------------

The elements of a :doc:`growth group <growth_group>` are equipped with
a partial order and usually contain a variable. Examples---the order
is described below these examples---are

- elements of the form `z^q` for some integer or rational `q`
  (growth groups with :ref:`description strings <growth_group_description>`
  ``z^ZZ`` or ``z^QQ``),

- elements of the form `\log(z)^q` for some integer or rational `q`
  (growth groups ``log(z)^ZZ`` or ``log(z)^QQ``),

- elements of the form `a^z` for some
  rational `a` (growth group ``QQ^z``), or

- more sophisticated constructions like products
  `x^r \cdot \log(x)^s \cdot a^y \cdot y^q`
  (this corresponds to an element of the growth group
  ``x^QQ * log(x)^ZZ * QQ^y * y^QQ``).

The order in all these examples is induced by the magnitude of the
elements as `x`, `y`, or `z` (independently) tend to `\infty`. For
elements only using the variable `z` this means that `g_1 \leq g_2` if

.. MATH::

    \lim_{z\to\infty} \frac{g_1}{g_2} \leq 1.

.. NOTE::

    Asymptotic rings where the variable tend to some value distinct from
    `\infty` are not yet implemented.

To find out more about

- growth groups,

- on how they are created and

- about the above used *descriptions strings*

see the top of the module :doc:`growth group <growth_group>`.


.. WARNING::

    As this code is experimental, a warning is thrown when an
    asymptotic ring (or an associated structure) is created for the
    first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G = GrowthGroup('x^ZZ')
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
        sage: T = GenericTermMonoid(G, ZZ)
        sage: R.<x, y> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ)
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.


.. _asymptotic_ring_intro:

Introductory Examples
=====================

We start this series of examples by defining two asymptotic rings.


Two Rings
---------

A Univariate Asymptotic Ring
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, we construct the following (very simple) asymptotic ring in the variable `z`::

    sage: A.<z> = AsymptoticRing(growth_group='z^QQ', coefficient_ring=ZZ); A
    Asymptotic Ring <z^QQ> over Integer Ring

A typical element of this ring is
::

    sage: A.an_element()
    z^(3/2) + O(z^(1/2))

This element consists of two summands: the exact term with coefficient
`1` and growth `z^{3/2}` and the `O`-term `O(z^{1/2})`. Note that the
growth of `z^{3/2}` is larger than the growth of `z^{1/2}` as
`z\to\infty`, thus this expansion cannot be simplified (which would
be done automatically, see below).

Elements can be constructed via the generator `z` and the function
:func:`~sage.rings.big_oh.O`, for example

::

    sage: 4*z^2 + O(z)
    4*z^2 + O(z)

A Multivariate Asymptotic Ring
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, we construct a more sophisticated asymptotic ring in the
variables `x` and `y` by
::

    sage: B.<x, y> = AsymptoticRing(growth_group='x^QQ * log(x)^ZZ * QQ^y * y^QQ', coefficient_ring=QQ); B
    Asymptotic Ring <x^QQ * log(x)^ZZ * QQ^y * y^QQ> over Rational Field

Again, we can look at a typical (nontrivial) element::

    sage: B.an_element()
    1/8*x^(3/2)*log(x)^3*(1/8)^y*y^(3/2) + O(x^(1/2)*log(x)*(1/2)^y*y^(1/2))

Again, elements can be created using the generators `x` and `y`, as well as
the function :func:`~sage.rings.big_oh.O`::

    sage: log(x)*y/42 + O(1/2^y)
    1/42*log(x)*y + O((1/2)^y)

Arithmetical Operations
-----------------------

In this section we explain how to perform various arithmetical
operations with the elements of the asymptotic rings constructed
above.


The Ring Operations Plus and Times
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We start our calculations in the ring
::

    sage: A
    Asymptotic Ring <z^QQ> over Integer Ring

Of course, we can perform the usual ring operations `+` and `*`::

    sage: z^2 + 3*z*(1-z)
    -2*z^2 + 3*z
    sage: (3*z + 2)^3
    27*z^3 + 54*z^2 + 36*z + 8

In addition to that, special powers---our growth group ``z^QQ`` allows
the exponents to be out of `\QQ`---can also be computed::

    sage: (z^(5/2)+z^(1/7)) * z^(-1/5)
    z^(23/10) + z^(-2/35)

The central concepts of computations with asymptotic expansions is
that the `O`-notation can be used. For example, we have
::

    sage: z^3 + z^2 + z + O(z^2)
    z^3 + O(z^2)

where the result is simplified automatically. A more sophisticated example is
::

    sage: (z+2*z^2+3*z^3+4*z^4) * (O(z)+z^2)
    4*z^6 + O(z^5)


Division
^^^^^^^^

The asymptotic expansions support division. For example, we can
expand `1/(z-1)` to a geometric series::

    sage: 1 / (z-1)
    z^(-1) + z^(-2) + z^(-3) + z^(-4) + ... + z^(-20) + O(z^(-21))

A default precision (parameter ``default_prec`` of
:class:`AsymptoticRing`) is predefined. Thus, only the first `20`
summands are calculated. However, if we only want the first `5` exact
terms, we cut of the rest by using
::

    sage: (1 / (z-1)).truncate(5)
    z^(-1) + z^(-2) + z^(-3) + z^(-4) + z^(-5) + O(z^(-6))

or

::

    sage: 1 / (z-1) + O(z^(-6))
    z^(-1) + z^(-2) + z^(-3) + z^(-4) + z^(-5) + O(z^(-6))

Of course, we can work with more complicated expansions as well::

    sage: (4*z+1) / (z^3+z^2+z+O(z^0))
    4*z^(-2) - 3*z^(-3) - z^(-4) + O(z^(-5))

Not all elements are invertible, for instance,

::

    sage: 1 / O(z)
    Traceback (most recent call last):
    ...
    ZeroDivisionError: Cannot invert O(z).

is not invertible, since it includes `0`.


Powers, Expontials and Logarithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It works as simple as it can be; just use the usual operators ``^``,
``exp`` and ``log``. For example, we obtain the usual series expansion
of the logarithm
::

    sage: -log(1-1/z)
    z^(-1) + 1/2*z^(-2) + 1/3*z^(-3) + ... + O(z^(-21))

as `z \to \infty`.

Similarly, we can apply the exponential function of an asymptotic expansion::

    sage: exp(1/z)
    1 + z^(-1) + 1/2*z^(-2) + 1/6*z^(-3) + 1/24*z^(-4) + ... + O(z^(-20))

Arbitrary powers work as well; for example, we have
::

    sage: (1 + 1/z + O(1/z^5))^(1 + 1/z)
    1 + z^(-1) + z^(-2) + 1/2*z^(-3) + 1/3*z^(-4) + O(z^(-5))


Multivariate Arithmetic
^^^^^^^^^^^^^^^^^^^^^^^

Now let us move on to arithmetic in the multivariate ring

::

    sage: B
    Asymptotic Ring <x^QQ * log(x)^ZZ * QQ^y * y^QQ> over Rational Field

.. TODO::

    write this part


More Examples
=============


The mathematical constant e as a limit
--------------------------------------

The base of the natural logarithm `e` satisfies the equation

.. MATH::

    e = \lim_{n\to\infty} \left(1+\frac{1}{n}\right)^n

By using asymptotic expansions, we obtain the more precise result
::

    sage: E.<n> = AsymptoticRing(growth_group='n^ZZ', coefficient_ring=SR, default_prec=5); E
    Asymptotic Ring <n^ZZ> over Symbolic Ring
    sage: (1 + 1/n)^n
    e - 1/2*e*n^(-1) + 11/24*e*n^(-2) - 7/16*e*n^(-3) + 2447/5760*e*n^(-4) + O(n^(-5))


Selected Technical Details
==========================


Coercions and Functorial Constructions
--------------------------------------

The :class:`AsymptoticRing` fully supports
`coercion <../../../../coercion/index.html>`_. For example, the coefficient ring is automatically extended when needed::

    sage: A
    Asymptotic Ring <z^QQ> over Integer Ring
    sage: (z + 1/2).parent()
    Asymptotic Ring <z^QQ> over Rational Field

Here, the coefficient ring was extended to allow `1/2` as a
coefficent. Another example is
::

    sage: C.<c> = AsymptoticRing(growth_group='c^ZZ', coefficient_ring=ZZ['e'])
    sage: C.an_element()
    e^3*c^3 + O(c)
    sage: C.an_element() / 7
    1/7*e^3*c^3 + O(c)

Here the result's coefficient ring is the newly found
::

    sage: (C.an_element() / 7).parent()
    Asymptotic Ring <c^ZZ> over
    Univariate Polynomial Ring in e over Rational Field

Not only the coefficient ring can be extended, but the growth group as
well. For example, we can add/multiply elements of the asymptotic
rings ``A`` and ``C`` to get an expansion of new asymptotic ring::

    sage: r = c*z + c/2 + O(z); r
    c*z + 1/2*c + O(z)
    sage: r.parent()
    Asymptotic Ring <c^ZZ * z^QQ> over
    Univariate Polynomial Ring in e over Rational Field


Data Structures
---------------

The summands of an
:class:`asymptotic expansion <AsymptoticExpansion>` are wrapped
:doc:`growth group elements <growth_group>`.
This wrapping is done by the
:doc:`term monoid module <term_monoid>`.
However, inside an
:class:`asymptotic expansion <AsymptoticExpansion>` these summands
(terms) are stored together with their growth-relationship, i.e., each
summand knows its direct predecessors and successors. As a data
structure a special poset (namely a
:mod:`mutable poset <sage.data_structures.mutable_poset>`)
is used. We can have a look at this::

    sage: b = x^3*y + x^2*y + x*y^2 + O(x) + O(y)
    sage: print b.summands.repr_full(reverse=True)
    poset(x*y^2, x^3*y, x^2*y, O(x), O(y))
    +-- oo
    |   +-- no successors
    |   +-- predecessors:   x*y^2, x^3*y
    +-- x*y^2
    |   +-- successors:   oo
    |   +-- predecessors:   O(x), O(y)
    +-- x^3*y
    |   +-- successors:   oo
    |   +-- predecessors:   x^2*y
    +-- x^2*y
    |   +-- successors:   x^3*y
    |   +-- predecessors:   O(x), O(y)
    +-- O(x)
    |   +-- successors:   x*y^2, x^2*y
    |   +-- predecessors:   null
    +-- O(y)
    |   +-- successors:   x*y^2, x^2*y
    |   +-- predecessors:   null
    +-- null
    |   +-- successors:   O(x), O(y)
    |   +-- no predecessors


Various
=======

AUTHORS:

- Benjamin Hackl (2015)
- Daniel Krenn (2015)

ACKNOWLEDGEMENT:

- Benjamin Hackl, Clemens Heuberger and Daniel Krenn are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

- Benjamin Hackl is supported by the Google Summer of Code 2015.


Classes and Methods
===================
"""

# *****************************************************************************
# Copyright (C) 2015 Benjamin Hackl <benjamin.hackl@aau.at>
#               2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************

from sage.rings.ring import Algebra
from sage.structure.element import CommutativeAlgebraElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.superseded import experimental

class AsymptoticExpansion(CommutativeAlgebraElement):
    r"""
    Class for asymptotic expansions, i.e., the elements of an
    :class:`AsymptoticRing`.

    INPUT:

    - ``parent`` -- the parent of the asymptotic expansion.

    - ``summands`` -- the summands as a
      :class:`~sage.data_structures.mutable_poset.MutablePoset`, which
      represents the underlying structure.

    - ``simplify`` -- a boolean (default: ``True``). It controls
      automatic simplification (absorption) of the asymptotic expansion.

    - ``convert`` -- a boolean (default: ``True``). If set, then the
      ``summands`` are converted to the asymptotic ring (the parent of this
      expansion). If not, then the summands are taken as they are. In
      that case, the caller must ensure that the parent of the terms is
      set correctly.

    EXAMPLES:

    There are several ways to create asymptotic expansions; usually
    this is done by using the corresponding :class:`asymptotic rings <AsymptoticRing>`::

        sage: R_x.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); R_x
        Asymptotic Ring <x^QQ> over Rational Field
        sage: R_y.<y> = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=ZZ); R_y
        Asymptotic Ring <y^ZZ> over Integer Ring

    At this point, `x` and `y` are already asymptotic expansions::

        sage: type(x)
        <class 'sage.rings.asymptotic.asymptotic_ring.AsymptoticRing_with_category.element_class'>

    The usual ring operations, but allowing rational exponents (growth
    group ``x^QQ``) can be performed::

        sage: x^2 + 3*(x - x^(2/5))
        x^2 + 3*x - 3*x^(2/5)
        sage: (3*x^(1/3) + 2)^3
        27*x + 54*x^(2/3) + 36*x^(1/3) + 8

    One of the central ideas behind computing with asymptotic
    expansions is that the `O`-notation (see
    :wikipedia:`Big_O_notation`) can be used. For example, we have::

        sage: (x+2*x^2+3*x^3+4*x^4) * (O(x)+x^2)
        4*x^6 + O(x^5)

    In particular, :meth:`~sage.rings.big_oh.O` can be used to
    construct the asymptotic expansions. With the help of the
    :meth:`summands`, we can also have a look at the inner structure
    of an asymptotic expansion::

        sage: expr1 = x + 2*x^2 + 3*x^3 + 4*x^4; expr2 = O(x) + x^2
        sage: print(expr1.summands.repr_full())
        poset(x, 2*x^2, 3*x^3, 4*x^4)
        +-- null
        |   +-- no predecessors
        |   +-- successors:   x
        +-- x
        |   +-- predecessors:   null
        |   +-- successors:   2*x^2
        +-- 2*x^2
        |   +-- predecessors:   x
        |   +-- successors:   3*x^3
        +-- 3*x^3
        |   +-- predecessors:   2*x^2
        |   +-- successors:   4*x^4
        +-- 4*x^4
        |   +-- predecessors:   3*x^3
        |   +-- successors:   oo
        +-- oo
        |   +-- predecessors:   4*x^4
        |   +-- no successors
        sage: print(expr2.summands.repr_full())
        poset(O(x), x^2)
        +-- null
        |   +-- no predecessors
        |   +-- successors:   O(x)
        +-- O(x)
        |   +-- predecessors:   null
        |   +-- successors:   x^2
        +-- x^2
        |   +-- predecessors:   O(x)
        |   +-- successors:   oo
        +-- oo
        |   +-- predecessors:   x^2
        |   +-- no successors
        sage: print((expr1 * expr2).summands.repr_full())
        poset(O(x^5), 4*x^6)
        +-- null
        |   +-- no predecessors
        |   +-- successors:   O(x^5)
        +-- O(x^5)
        |   +-- predecessors:   null
        |   +-- successors:   4*x^6
        +-- 4*x^6
        |   +-- predecessors:   O(x^5)
        |   +-- successors:   oo
        +-- oo
        |   +-- predecessors:   4*x^6
        |   +-- no successors

    In addition to the monomial growth elements from above, we can
    also compute with logarithmic terms (simply by constructing the
    appropriate growth group)::

        sage: R_log = AsymptoticRing(growth_group='log(x)^QQ', coefficient_ring=QQ)
        sage: lx = R_log(log(SR.var('x')))
        sage: (O(lx) + lx^3)^4
        log(x)^12 + O(log(x)^10)

    .. SEEALSO::

        :doc:`growth_group`,
        :doc:`term_monoid`,
        :mod:`~sage.data_structures.mutable_poset`.
    """
    def __init__(self, parent, summands, simplify=True, convert=True):
        r"""
        See :class:`AsymptoticExpansion` for more information.

        TESTS::

            sage: R_x.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: R_y.<y> = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=ZZ)
            sage: R_x is R_y
            False
            sage: ex1 = x + 2*x^2 + 3*x^3 + 4*x^4 + 5*x^5
            sage: ex2 = x + O(R_x(1))
            sage: ex1 * ex2
            5*x^6 + O(x^5)

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: OT = TermMonoid('O', G, ZZ); ET = TermMonoid('exact', G, ZZ)
            sage: R = AsymptoticRing(G, ZZ)
            sage: lst = [ET(x, 1), ET(x^2, 2), OT(x^3), ET(x^4, 4)]
            sage: expr = R(lst, simplify=False); expr  # indirect doctest
            4*x^4 + O(x^3) + 2*x^2 + x
            sage: print expr.summands.repr_full()
            poset(x, 2*x^2, O(x^3), 4*x^4)
            +-- null
            |   +-- no predecessors
            |   +-- successors:   x
            +-- x
            |   +-- predecessors:   null
            |   +-- successors:   2*x^2
            +-- 2*x^2
            |   +-- predecessors:   x
            |   +-- successors:   O(x^3)
            +-- O(x^3)
            |   +-- predecessors:   2*x^2
            |   +-- successors:   4*x^4
            +-- 4*x^4
            |   +-- predecessors:   O(x^3)
            |   +-- successors:   oo
            +-- oo
            |   +-- predecessors:   4*x^4
            |   +-- no successors
            sage: expr._simplify_(); expr
            4*x^4 + O(x^3)
            sage: print expr.summands.repr_full()
            poset(O(x^3), 4*x^4)
            +-- null
            |   +-- no predecessors
            |   +-- successors:   O(x^3)
            +-- O(x^3)
            |   +-- predecessors:   null
            |   +-- successors:   4*x^4
            +-- 4*x^4
            |   +-- predecessors:   O(x^3)
            |   +-- successors:   oo
            +-- oo
            |   +-- predecessors:   4*x^4
            |   +-- no successors
            sage: R(lst, simplify=True) # indirect doctest
            4*x^4 + O(x^3)

        ::

            sage: R.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ)
            sage: e = R(x^2 + O(x))
            sage: from sage.rings.asymptotic.asymptotic_ring import AsymptoticExpansion
            sage: S = AsymptoticRing(growth_group='x^QQ', coefficient_ring=ZZ)
            sage: for s in AsymptoticExpansion(S, e.summands).summands.elements_topological():
            ....:     print s.parent()
            O-Term Monoid x^QQ with implicit coefficients in Integer Ring
            Exact Term Monoid x^QQ with coefficients in Integer Ring
            sage: for s in AsymptoticExpansion(S, e.summands,
            ....:         convert=False).summands.elements_topological():
            ....:     print s.parent()
            O-Term Monoid x^QQ with implicit coefficients in Rational Field
            Exact Term Monoid x^QQ with coefficients in Rational Field

        ::

            sage: AsymptoticExpansion(S, R(1/2).summands)
            Traceback (most recent call last):
            ...
            ValueError: Cannot include 1/2 with parent
            Exact Term Monoid x^QQ with coefficients in Rational Field in
            Asymptotic Ring <x^QQ> over Integer Ring
            > *previous* ValueError: 1/2 is not a coefficient in
            Exact Term Monoid x^QQ with coefficients in Integer Ring.
        """
        super(AsymptoticExpansion, self).__init__(parent=parent)

        from sage.data_structures.mutable_poset import MutablePoset
        if not isinstance(summands, MutablePoset):
            raise TypeError('Summands %s are not in a mutable poset as expected '
                            'when creating an element of %s.' % (summands, parent))

        if convert:
            from misc import combine_exceptions
            from term_monoid import TermMonoid
            def convert_terms(element):
                T = TermMonoid(term=element.parent(), asymptotic_ring=parent)
                try:
                    return T(element)
                except (ValueError, TypeError) as e:
                    raise combine_exceptions(
                        ValueError('Cannot include %s with parent %s in %s' %
                                   (element, element.parent(), parent)), e)
            new_summands = summands.copy()
            new_summands.map(convert_terms, topological=True, reverse=True)
            self._summands_ = new_summands
        else:
            self._summands_ = summands

        if simplify:
            self._simplify_()


    @property
    def summands(self):
        r"""
        The summands of this asymptotic expansion stored in the
        underlying data structure (a
        :class:`~sage.data_structures.mutable_poset.MutablePoset`).

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: expr = 7*x^12 + x^5 + O(x^3)
            sage: expr.summands
            poset(O(x^3), x^5, 7*x^12)

        .. SEEALSO::

            :class:`sage.data_structures.mutable_poset.MutablePoset`
        """
        return self._summands_


    def __hash__(self):
        r"""
        A hash value for this element.

        .. WARNING::

            This hash value uses the string representation and might not be
            always right.

        TESTS::

            sage: R_log = AsymptoticRing(growth_group='log(x)^QQ', coefficient_ring=QQ)
            sage: lx = R_log(log(SR.var('x')))
            sage: elt = (O(lx) + lx^3)^4
            sage: hash(elt) # random
            -4395085054568712393
        """
        return hash(str(self))


    def __nonzero__(self):
        r"""
        Return whether this asymptotic expansion is not identically zero.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        TESTS::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: bool(R(0))  # indirect doctest
            False
            sage: bool(x)  # indirect doctest
            True
            sage: bool(7*x^12 + x^5 + O(x^3))  # indirect doctest
            True
        """
        return bool(self._summands_)


    def __eq__(self, other):
        r"""
        Return whether this asymptotic expansion is equal to ``other``.

        INPUT:

        - ``other`` -- an object.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', QQ)
            sage: (1 + 2*x + 3*x^2) == (3*x^2 + 2*x + 1)  # indirect doctest
            True
            sage: O(x) == O(x)
            False

        TESTS::

            sage: x == None
            False

        ::

            sage: x == 'x'
            False
        """
        if other is None:
            return False
        try:
            return not bool(self - other)
        except (TypeError, ValueError):
            return False


    def __ne__(self, other):
        r"""
        Return whether this asymptotic expansion is not equal to ``other``.

        INPUT:

        - ``other`` -- an object.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', QQ)
            sage: (1 + 2*x + 3*x^2) != (3*x^2 + 2*x + 1)  # indirect doctest
            False
            sage: O(x) != O(x)
            True

        TESTS::

            sage: x != None
            True
        """
        return not self == other


    def has_same_summands(self, other):
        r"""
        Return whether this asymptotic expansion and ``other`` have the
        same summands.

        INPUT:

        - ``other`` -- an asymptotic expansion.

        OUTPUT:

        A boolean.

        .. NOTE::

            While for example ``O(x) == O(x)`` yields ``False``,
            these expansions *do* have the same summands and this method
            returns ``True``.

            Moreover, this method uses the coercion model in order to
            find a common parent for this asymptotic expansion and
            ``other``.

        EXAMPLES::

            sage: R_ZZ.<x_ZZ> = AsymptoticRing('x^ZZ', ZZ)
            sage: R_QQ.<x_QQ> = AsymptoticRing('x^ZZ', QQ)
            sage: sum(x_ZZ^k for k in range(5)) == sum(x_QQ^k for k in range(5))  # indirect doctest
            True
            sage: O(x_ZZ) == O(x_QQ)
            False

        TESTS::

            sage: x_ZZ.has_same_summands(None)
            False
        """
        if other is None:
            return False
        from sage.structure.element import have_same_parent
        if have_same_parent(self, other):
            return self._has_same_summands_(other)

        from sage.structure.element import get_coercion_model
        return get_coercion_model().bin_op(self, other,
                                           lambda self, other:
                                           self._has_same_summands_(other))


    def _has_same_summands_(self, other):
        r"""
        Return whether this :class:`AsymptoticExpansion` has the same
        summands as ``other``.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpansion`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method compares two :class:`AsymptoticExpansion`
            with the same parent.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', QQ)
            sage: O(x).has_same_summands(O(x))
            True
            sage: (1 + x + 2*x^2).has_same_summands(2*x^2 + O(x))  # indirect doctest
            False
        """
        if len(self.summands) != len(other.summands):
            return False
        from itertools import izip
        return all(s == o for s, o in
                   izip(self.summands.elements_topological(),
                        other.summands.elements_topological()))


    def _simplify_(self):
        r"""
        Simplify this asymptotic expansion.

        INPUT:

        Nothing.

        OUTPUT:

        Nothing, but modifies this asymptotic expansion.

        .. NOTE::

            This method is usually called during initialization of
            this asymptotic expansion.

        .. NOTE::

            This asymptotic expansion is simplified by letting
            `O`-terms that are included in this expansion absorb all
            terms with smaller growth.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import TermMonoid
            sage: G = GrowthGroup('x^ZZ')
            sage: OT = TermMonoid('O', G, ZZ); ET = TermMonoid('exact', G, ZZ)
            sage: R = AsymptoticRing(G, ZZ)
            sage: lst = [ET(x, 1), ET(x^2, 2), OT(x^3), ET(x^4, 4)]
            sage: expr = R(lst, simplify=False); expr  # indirect doctest
            4*x^4 + O(x^3) + 2*x^2 + x
            sage: expr._simplify_(); expr
            4*x^4 + O(x^3)
            sage: R(lst)  # indirect doctest
            4*x^4 + O(x^3)
        """
        self._summands_.merge(reverse=True)


    def _repr_(self):
        r"""
        A representation string for this asymptotic expansion.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: (5*x^2+12*x) * (x^3+O(x))  # indirect doctest
            5*x^5 + 12*x^4 + O(x^3)
            sage: (5*x^2-12*x) * (x^3+O(x))  # indirect doctest
            5*x^5 - 12*x^4 + O(x^3)
        """
        s = ' + '.join(repr(elem) for elem in
                       self.summands.elements_topological(reverse=True))
        s = s.replace('+ -', '- ')
        if not s:
            return '0'
        return s


    def _add_(self, other):
        r"""
        Add ``other`` to this asymptotic expansion.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpansion`.

        OUTPUT:

        The sum as an :class:`AsymptoticExpansion`.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: expr1 = x^123; expr2 = x^321
            sage: expr1._add_(expr2)
            x^321 + x^123
            sage: expr1 + expr2  # indirect doctest
            x^321 + x^123

        If an `O`-term is added to an asymptotic expansion, then
        the `O`-term absorbs everything it can::

            sage: x^123 + x^321 + O(x^555)  # indirect doctest
            O(x^555)

        TESTS::

            sage: x + O(x)
            O(x)
            sage: O(x) + x
            O(x)
        """
        return self.parent()(self.summands.union(other.summands),
                             simplify=True, convert=False)


    def _sub_(self, other):
        r"""
        Subtract ``other`` from this asymptotic expansion.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpansion`.

        OUTPUT:

        The difference as an :class:`AsymptoticExpansion`.

        .. NOTE::

            Subtraction of two asymptotic expansions is implemented
            by means of addition: `e_1 - e_2 = e_1 + (-1)\cdot e_2`.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: expr1 = x^123; expr2 = x^321
            sage: expr1 - expr2  # indirect doctest
            -x^321 + x^123
            sage: O(x) - O(x)
            O(x)
        """
        return self + self.parent().coefficient_ring(-1)*other


    def _mul_term_(self, term):
        r"""
        Helper method: multiply this asymptotic expansion by the
        asymptotic term ``term``.

        INPUT:

        - ``term`` -- an asymptotic term (see
          :doc:`term_monoid`).

        OUTPUT:

        The product as an :class:`AsymptoticExpansion`.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import OTermMonoid
            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: T = OTermMonoid(R.growth_group, ZZ)
            sage: expr = 10*x^2 + O(x)
            sage: t = T(R.growth_group.gen())
            sage: expr._mul_term_(t)
            O(x^3)
        """
        from term_monoid import ExactTerm
        simplify = not isinstance(term, ExactTerm)
        return self.parent()(self.summands.mapped(lambda element: term * element),
                             simplify=simplify, convert=False)


    def _mul_(self, other):
        r"""
        Multiply this asymptotic expansion by another asymptotic expansion ``other``.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpansion`.

        OUTPUT:

        The product as an :class:`AsymptoticExpansion`.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: ex1 = 5*x^12
            sage: ex2 = x^3 + O(x)
            sage: ex1 * ex2  # indirect doctest
            5*x^15 + O(x^13)

        .. TODO::

            The current implementation is the standard long
            multiplication. More efficient variants like Karatsuba
            multiplication, or methods that exploit the structure
            of the underlying poset shall be implemented at a later
            point.
        """
        return sum(self._mul_term_(term_other) for
                   term_other in other.summands.elements())


    def _rmul_(self, other):
        r"""
        Multiply this asymptotic expansion by an element ``other`` of its
        coefficient ring.

        INPUT:

        - ``other`` -- an element of the coefficient ring.

        OUTPUT:

        An :class:`AsymptoticExpansion`.

        TESTS::

            sage: A.<a> = AsymptoticRing(growth_group='QQ^a * a^QQ * log(a)^QQ', coefficient_ring=ZZ)
            sage: 2*a # indirect doctest
            2*a
        """
        if other.is_zero():
            return self.parent().zero()

        from term_monoid import TermMonoid
        E = TermMonoid('exact', asymptotic_ring=self.parent())
        e = E(self.parent().growth_group.one(), coefficient=other)
        return self._mul_term_(e)


    _lmul_ = _rmul_


    def _div_(self, other):
        r"""
        Divide this element through ``other``.

        INPUT:

        - ``other`` -- an asymptotic expansion.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', QQ, default_prec=5)
            sage: 1/x^42
            x^(-42)
            sage: (1 + 4*x) / (x + 2*x^2)
            2*x^(-1) - 1/2*x^(-2) + 1/4*x^(-3) - 1/8*x^(-4) + 1/16*x^(-5) + O(x^(-6))
            sage: x / O(x)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot invert O(x).
        """
        return self * ~other


    def __invert__(self, precision=None):
        r"""
        Return the multiplicative inverse of this element.

        INPUT:

        - ``precision`` -- the precision used for truncating the
          expansion. If ``None`` (default value) is used, the
          default precision of the parent is used.

        OUTPUT:

        An asymptotic expansion.

        .. WARNING::

            Due to truncation of infinite expansions, the element
            returned by this method might not fulfill
            ``el * ~el == 1``.

        .. TODO::

            As soon as `L`-terms are implemented, this
            implementation has to be adapted as well in order to
            yield correct results.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ, default_prec=4)
            sage: ~x
            x^(-1)
            sage: ~(x^42)
            x^(-42)
            sage: ex = ~(1 + x); ex
            x^(-1) - x^(-2) + x^(-3) - x^(-4) + O(x^(-5))
            sage: ex * (1+x)
            1 + O(x^(-4))
            sage: ~(1 + O(1/x))
            1 + O(x^(-1))

        TESTS::

            sage: A.<a> = AsymptoticRing(growth_group='a^ZZ', coefficient_ring=ZZ)
            sage: (1 / a).parent()
            Asymptotic Ring <a^ZZ> over Rational Field
            sage: (a / 2).parent()
            Asymptotic Ring <a^ZZ> over Rational Field

        ::

            sage: ~A(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Division by zero in 0.

        ::

            sage: B.<s, t> = AsymptoticRing(growth_group='s^ZZ * t^ZZ', coefficient_ring=QQ)
            sage: ~(s + t)
            Traceback (most recent call last):
            ...
            ValueError: Expansion s + t cannot be inverted since there are
            several maximal elements s, t.
        """
        if not self.summands:
            raise ZeroDivisionError('Division by zero in %s.' % (self,))

        elif len(self.summands) == 1:
            element = next(self.summands.elements())
            return self.parent()._create_element_in_extension_(
                ~element, element.parent())

        max_elem = tuple(self.summands.maximal_elements())
        if len(max_elem) != 1:
            raise ValueError('Expansion %s cannot be inverted since there '
                             'are several maximal elements %s.' %
                             (self, ', '.join(str(e) for e in
                                              sorted(max_elem, key=str))))
        max_elem = max_elem[0]

        imax_elem = ~max_elem
        if imax_elem.parent() is max_elem.parent():
            new_self = self
        else:
            new_self = self.parent()._create_element_in_extension_(
                imax_elem, max_elem.parent()).parent()(self)

        one = new_self.parent().one()
        geom = one - new_self._mul_term_(imax_elem)

        expanding = True
        result = one
        while expanding:
            new_result = (geom*result + one).truncate(precision=precision)
            if new_result.has_same_summands(result):
                expanding = False
            result = new_result
        return result._mul_term_(imax_elem)


    invert = __invert__


    def truncate(self, precision=None):
        r"""
        Truncate this asymptotic expansion.

        INPUT:

        - ``precision`` -- a positive integer or ``None``. Number of
          summands that are kept. If ``None`` (default value) is
          given, then ``default_prec`` from the parent is used.

        OUTPUT:

        An asymptotic expansion.

        .. NOTE::

            For example, truncating an asymptotic expansion with
            ``precision=20`` does not yield an expansion with exactly 20
            summands! Rather than that, it keeps the 20 summands
            with the largest growth, and adds appropriate
            `O`-Terms.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', QQ)
            sage: ex = sum(x^k for k in range(5)); ex
            x^4 + x^3 + x^2 + x + 1
            sage: ex.truncate(precision=2)
            x^4 + x^3 + O(x^2)
            sage: ex.truncate(precision=0)
            O(x^4)
            sage: ex.truncate()
            x^4 + x^3 + x^2 + x + 1
        """
        if precision is None:
            precision = self.parent().default_prec

        if len(self.summands) <= precision:
            return self

        summands = self.summands.copy()
        from term_monoid import TermMonoid
        def convert_terms(element):
            if convert_terms.count < precision:
                convert_terms.count += 1
                return element
            T = TermMonoid(term='O', asymptotic_ring=self.parent())
            return T(element)
        convert_terms.count = 0
        summands.map(convert_terms, topological=True, reverse=True)
        return self.parent()(summands, simplify=True, convert=False)


    def __pow__(self, exponent, precision=None):
        r"""
        Calculate the power of this asymptotic expansion to the given ``exponent``.

        INPUT:

        - ``exponent`` -- an element.

        - ``precision`` -- the precision used for truncating the
          expansion. If ``None`` (default value) is used, the
          default precision of the parent is used.

        OUTPUT:

        An asymptotic expansion.

        TESTS::

            sage: R_QQ.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ)
            sage: x^(1/7)
            x^(1/7)
            sage: R_ZZ.<y> = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=ZZ)
            sage: y^(1/7)
            y^(1/7)
            sage: (y^(1/7)).parent()
            Asymptotic Ring <y^QQ> over Rational Field
            sage: (x^(1/2) + O(x^0))^15
            x^(15/2) + O(x^7)
            sage: (y^2 + O(y))^(1/2)  # not tested, see #19316
            y + O(1)
            sage: (y^2 + O(y))^(-2)
            y^(-4) + O(y^(-5))

        ::

            sage: B.<z> = AsymptoticRing(growth_group='z^QQ * log(z)^QQ', coefficient_ring=QQ)
            sage: (z^2 + O(z))^(1/2)
            z + O(1)

        ::

            sage: A.<x> = AsymptoticRing('QQ^x * x^SR * log(x)^ZZ', QQ)
            sage: x * 2^x
            2^x*x
            sage: 5^x * 2^x
            10^x
            sage: 2^log(x)
            x^(log(2))
            sage: 2^(x + 1/x)
            2^x + log(2)*2^x*x^(-1) + 1/2*log(2)^2*2^x*x^(-2) + ... + O(2^x*x^(-20))
            sage: _.parent()
            Asymptotic Ring <QQ^x * x^SR * log(x)^QQ> over Symbolic Ring

        See :trac:`19110`::

            sage: O(x)^(-1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot take O(x) to exponent -1.
            > *previous* ZeroDivisionError: rational division by zero

        ::

            sage: B.<z> = AsymptoticRing(growth_group='z^QQ * log(z)^QQ', coefficient_ring=QQ, default_prec=5)
            sage: z^(1+1/z)
            z + log(z) + 1/2*z^(-1)*log(z)^2 + 1/6*z^(-2)*log(z)^3 +
            1/24*z^(-3)*log(z)^4 + O(z^(-4)*log(z)^5)

        ::

            sage: B(0)^(-7)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot take 0 to the negative exponent -7.
            sage: B(0)^SR.var('a')
            Traceback (most recent call last):
            ...
            NotImplementedError: Taking 0 to the exponent a not implemented.

        ::

            sage: C.<s, t> = AsymptoticRing(growth_group='s^QQ * t^QQ', coefficient_ring=QQ)
            sage: (s + t)^s
            Traceback (most recent call last):
            ...
            ValueError: Cannot take s + t to the exponent s.
            > *previous* ValueError: log(s + t) cannot be constructed since
            there are several maximal elements s, t.
        """
        if not self.summands:
            if exponent == 0:
                return self.parent().one()
            elif exponent > 0:
                return self.parent().zero()
            elif exponent < 0:
                raise ZeroDivisionError('Cannot take %s to the negative exponent %s.' %
                                        (self, exponent))
            else:
                raise NotImplementedError('Taking %s to the exponent %s not implemented.' %
                                          (self, exponent))

        elif len(self.summands) == 1:
            element = next(self.summands.elements())
            if isinstance(exponent, AsymptoticExpansion) and element.is_constant():
                return exponent.rpow(base=element.coefficient, precision=precision)
            try:
                return self.parent()._create_element_in_extension_(
                    element ** exponent, element.parent())
            except (ArithmeticError, TypeError, ValueError):
                if not isinstance(exponent, AsymptoticExpansion):
                    raise

        from sage.rings.integer_ring import ZZ
        try:
            exponent = ZZ(exponent)
        except (TypeError, ValueError):
            pass
        else:
            return super(AsymptoticExpansion, self).__pow__(exponent)

        try:
            return (exponent * self.log(precision=precision)).exp(precision=precision)
        except (TypeError, ValueError, ZeroDivisionError) as e:
            from misc import combine_exceptions
            raise combine_exceptions(
                ValueError('Cannot take %s to the exponent %s.' % (self, exponent)), e)


    pow = __pow__


    def O(self):
        r"""
        Convert all terms in this asymptotic expansion to `O`-terms.

        INPUT:

        Nothing.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: AR.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: O(x)
            O(x)
            sage: type(O(x))
            <class 'sage.rings.asymptotic.asymptotic_ring.AsymptoticRing_with_category.element_class'>
            sage: expr = 42*x^42 + x^10 + O(x^2); expr
            42*x^42 + x^10 + O(x^2)
            sage: expr.O()
            O(x^42)
            sage: O(AR(0))
            0
            sage: (2*x).O()
            O(x)

        .. SEEALSO::

            :func:`sage.rings.power_series_ring.PowerSeriesRing`,
            :func:`sage.rings.laurent_series_ring.LaurentSeriesRing`.
        """
        return sum(self.parent().create_summand('O', growth=element)
                   for element in self.summands.maximal_elements())


    def log(self, base=None, precision=None):
        r"""
        The logarithm of this asymptotic expansion.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.

        - ``precision`` -- the precision used for truncating the
          expansion. If ``None`` (default value) is used, the
          default precision of the parent is used.

        OUTPUT:

        An asymptotic expansion.

        .. NOTE::

            Computing the logarithm of an asymptotic expansion
            is possible if and only if there is exactly one maximal
            summand in the expansion.

        ALGORITHM:

        If the expansion has more than one summand,
        the asymptotic expansion for `\log(1+t)` as `t` tends to `0`
        is used.

        .. TODO::

            As soon as `L`-terms are implemented, this
            implementation has to be adapted as well in order to
            yield correct results.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ * log(x)^ZZ', coefficient_ring=QQ)
            sage: log(x)
            log(x)
            sage: log(x^2)
            2*log(x)
            sage: log(x-1)
            log(x) - x^(-1) - 1/2*x^(-2) - 1/3*x^(-3) - ... + O(x^(-21))

        TESTS::

            sage: log(R(1))
            0
            sage: log(R(0))
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot compute log(0) in
            Asymptotic Ring <x^ZZ * log(x)^ZZ> over Rational Field.
            sage: C.<s, t> = AsymptoticRing(growth_group='s^ZZ * t^ZZ', coefficient_ring=QQ)
            sage: log(s + t)
            Traceback (most recent call last):
            ...
            ValueError: log(s + t) cannot be constructed since there are
            several maximal elements s, t.
        """
        P = self.parent()

        if not self.summands:
            raise ArithmeticError('Cannot compute log(0) in %s.' % (self.parent(),))

        elif len(self.summands) == 1:
            if self.is_one():
                return P.zero()
            element = next(self.summands.elements())
            return sum(P._create_element_in_extension_(l, element.parent())
                       for l in element.log_term(base=base))

        max_elem = tuple(self.summands.maximal_elements())
        if len(max_elem) != 1:
            raise ValueError('log(%s) cannot be constructed since there '
                             'are several maximal elements %s.' %
                             (self, ', '.join(str(e) for e in
                                              sorted(max_elem, key=str))))
        max_elem = max_elem[0]

        imax_elem = ~max_elem
        if imax_elem.parent() is max_elem.parent():
            new_self = self
        else:
            new_self = P._create_element_in_extension_(
                imax_elem, max_elem.parent()).parent()(self)

        one = new_self.parent().one()
        geom = one - new_self._mul_term_(imax_elem)

        from sage.rings.integer_ring import ZZ
        expanding = True
        result = -geom
        geom_k = geom
        k = ZZ(1)
        while expanding:
            k += ZZ(1)
            geom_k *= geom
            new_result = (result - geom_k * ~k).truncate(precision=precision)
            if new_result.has_same_summands(result):
                expanding = False
            result = new_result

        result += new_self.parent()(max_elem).log()
        if base:
            from sage.functions.log import log
            result = result / log(base)
        return result


    def is_little_o_of_one(self):
        r"""
        Return whether this expansion is of order `o(1)`.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: A.<x> = AsymptoticRing('x^ZZ * log(x)^ZZ', QQ)
            sage: (x^4 * log(x)^(-2) + x^(-4) * log(x)^2).is_little_o_of_one()
            False
            sage: (x^(-1) * log(x)^1234 + x^(-2) + O(x^(-3))).is_little_o_of_one()
            True
            sage: (log(x) - log(x-1)).is_little_o_of_one()
            True

        ::

            sage: A.<x, y> = AsymptoticRing('x^QQ * y^QQ * log(y)^ZZ', QQ)
            sage: (x^(-1/16) * y^32 + x^32 * y^(-1/16)).is_little_o_of_one()
            False
            sage: (x^(-1) * y^(-3) + x^(-3) * y^(-1)).is_little_o_of_one()
            True
            sage: (x^(-1) * y / log(y)).is_little_o_of_one()
            False
            sage: (log(y-1)/log(y) - 1).is_little_o_of_one()
            True
        """
        return all(term.is_little_o_of_one() for term in self.summands.maximal_elements())


    def rpow(self, base, precision=None):
        r"""
        Return the power of ``base`` to this asymptotic expansion.

        INPUT:

        - ``base`` -- an element or ``'e'``.

        - ``precision`` -- the precision used for truncating the
          expansion. If ``None`` (default value) is used, the
          default precision of the parent is used.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: A.<x> = AsymptoticRing('x^ZZ', QQ)
            sage: (1/x).rpow('e', precision=5)
            1 + x^(-1) + 1/2*x^(-2) + 1/6*x^(-3) + 1/24*x^(-4) + O(x^(-5))

        TESTS::

            sage: x.rpow(SR.var('y'))
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot construct y^x in Growth Group x^ZZ
            > *previous* TypeError: unsupported operand parent(s) for '*':
            'Growth Group x^ZZ' and 'Growth Group SR^x'
        """
        if isinstance(base, AsymptoticExpansion):
            return base.__pow__(self, precision=precision)

        P = self.parent()

        # first: remove terms from a copy of this term such that a
        # term in o(1) remains

        expr_o = self.summands.copy()
        large_terms = []
        for term in self.summands.elements_topological():
            if not term.is_little_o_of_one():
                large_terms.append(term)
                expr_o.remove(term.growth)

        expr_o = P(expr_o)

        # next: try to take the exponential function of the large elements

        try:
            large_result = P.prod(
                P._create_element_in_extension_(term.rpow(base),
                                              term.parent())
                for term in large_terms)
        except (TypeError, ValueError) as e:
            from misc import combine_exceptions
            raise combine_exceptions(
                ValueError('Cannot construct the power of %s to the '
                           'exponent %s in %s.' %
                           (base, self, self.parent())), e)

        # then: expand expr_o

        if not expr_o:
            return large_result


        if base == 'e':
            geom = expr_o
        else:
            from sage.functions.log import log
            geom = expr_o * log(base)
        P = geom.parent()

        expanding = True
        result = P.one()
        geom_k = P.one()
        from sage.rings.integer_ring import ZZ
        k = ZZ(0)
        while expanding:
            k += ZZ(1)
            geom_k *= geom
            new_result = (result + geom_k / k.factorial()).truncate(precision=precision)
            if new_result.has_same_summands(result):
                expanding = False
            result = new_result

        return result * large_result


    def exp(self, precision=None):
        r"""
        Return the exponential of (i.e., the power of `e` to) this asymptotic expansion.

        INPUT:

        - ``precision`` -- the precision used for truncating the
          expansion. If ``None`` (default value) is used, the
          default precision of the parent is used.

        OUTPUT:

        An asymptotic expansion.

        .. NOTE::

            The exponential function of this expansion can only be
            computed exactly if the respective growth element can be
            constructed in the underlying growth group.

        ALGORITHM:

        If the corresponding growth can be constructed, return
        the exact exponential function. Otherwise, if this term
        is `o(1)`, try to expand the series and truncate
        according to the given precision.

        .. TODO::

            As soon as `L`-terms are implemented, this
            implementation has to be adapted as well in order to
            yield correct results.

        EXAMPLES::

            sage: A.<x> = AsymptoticRing('(e^x)^ZZ * x^ZZ * log(x)^ZZ', SR)
            sage: exp(x)
            e^x
            sage: exp(2*x)
            (e^x)^2
            sage: exp(x + log(x))
            e^x*x

        ::

            sage: (x^(-1)).exp(precision=7)
            1 + x^(-1) + 1/2*x^(-2) + 1/6*x^(-3) + ... + O(x^(-7))

        TESTS::

            sage: A.<x> = AsymptoticRing('(e^x)^ZZ * x^QQ * log(x)^QQ', SR)
            sage: exp(log(x))
            x
            sage: log(exp(x))
            x

        ::

            sage: exp(x+1)
            e*e^x
        """
        return self.rpow('e', precision=precision)


    def substitute(self, rules=None, domain=None, **kwds):
        r"""
        Substitute the given ``rules`` in this asymptotic expansion.

        INPUT:

        - ``rules`` -- a dictionary.

        - ``kwds`` -- keyword arguments will be added to the
          substitution ``rules``.

        - ``domain`` -- (default: ``None``) a parent. The neutral
          elements `0` and `1` (rules for the keys ``'_zero_'`` and
          ``'_one_'``, see note box below) are taken out of this
          domain. If ``None``, then this is determined automatically.

        OUTPUT:

        An object.

        .. NOTE::

          The neutral element of the asymptotic ring is replaced by
          the value to the key ``'_zero_'``; the neutral element of
          the growth group is replaced by the value to the key
          ``'_one_'``.

        EXAMPLES::

            sage: A.<x> = AsymptoticRing(growth_group='(e^x)^QQ * x^ZZ * log(x)^ZZ', coefficient_ring=QQ, default_prec=5)

        ::

            sage: (e^x * x^2 + log(x)).subs(x=SR('s'))
            s^2*e^s + log(s)
            sage: _.parent()
            Symbolic Ring

        ::

            sage: (x^3 + x + log(x)).subs(x=x+5).truncate(5)
            x^3 + 15*x^2 + 76*x + log(x) + 130 + O(x^(-1))
            sage: _.parent()
            Asymptotic Ring <(e^x)^QQ * x^ZZ * log(x)^ZZ> over Rational Field

        ::

            sage: (e^x * x^2 + log(x)).subs(x=2*x)
            4*(e^x)^2*x^2 + log(x) + log(2)
            sage: _.parent()
            Asymptotic Ring <(e^x)^QQ * x^QQ * log(x)^QQ> over Symbolic Ring

        ::

            sage: (x^2 + log(x)).subs(x=4*x+2).truncate(5)
            16*x^2 + 16*x + log(x) + log(4) + 4 + 1/2*x^(-1) + O(x^(-2))
            sage: _.parent()
            Asymptotic Ring <(e^x)^QQ * x^ZZ * log(x)^ZZ> over Symbolic Ring

        ::

            sage: (e^x * x^2 + log(x)).subs(x=RIF(pi))
            229.534211738584?
            sage: _.parent()
            Real Interval Field with 53 bits of precision

        .. SEEALSO::

            :meth:`sage.symbolic.expression.Expression.subs`

        TESTS::

            sage: x.subs({'y': -1})
            Traceback (most recent call last):
            ...
            ValueError: Cannot substitute y in x since it is not a generator of
            Asymptotic Ring <(e^x)^QQ * x^ZZ * log(x)^ZZ> over Rational Field.
            sage: B.<u, v, w> = AsymptoticRing(growth_group='u^QQ * v^QQ * w^QQ', coefficient_ring=QQ)
            sage: (1/u).subs({'u': 0})
            Traceback (most recent call last):
            ...
            TypeError: Cannot apply the substitution rules {u: 0} on u^(-1) in
            Asymptotic Ring <u^QQ * v^QQ * w^QQ> over Rational Field.
            > *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Asymptotic Ring <u^QQ * v^QQ * w^QQ> over Rational Field.
            >> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Exact Term Monoid u^QQ * v^QQ * w^QQ with coefficients in Rational Field.
            >...> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Growth Group u^QQ * v^QQ * w^QQ.
            >...> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Growth Group u^QQ.
            >...> *previous* ZeroDivisionError: rational division by zero
            sage: (1/u).subs({'u': 0, 'v': SR.var('v')})
            Traceback (most recent call last):
            ...
            TypeError: Cannot apply the substitution rules {u: 0, v: v} on u^(-1) in
            Asymptotic Ring <u^QQ * v^QQ * w^QQ> over Rational Field.
            > *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Asymptotic Ring <u^QQ * v^QQ * w^QQ> over Rational Field.
            >> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Exact Term Monoid u^QQ * v^QQ * w^QQ with coefficients in Rational Field.
            >...> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Growth Group u^QQ * v^QQ * w^QQ.
            >...> *previous* ZeroDivisionError: Cannot substitute in u^(-1) in
            Growth Group u^QQ.
            >...> *previous* ZeroDivisionError: rational division by zero

        ::

            sage: u.subs({u: 0, 'v': SR.var('v')})
            0
            sage: v.subs({u: 0, 'v': SR.var('v')})
            v
            sage: _.parent()
            Symbolic Ring

        ::

            sage: u.subs({SR.var('u'): -1})
            Traceback (most recent call last):
            ...
            TypeError: Cannot substitute u in u since it is neither an
            asymptotic expansion nor a string
            (but a <type 'sage.symbolic.expression.Expression'>).

        ::

            sage: u.subs({u: 1, 'u': 1})
            1
            sage: u.subs({u: 1}, u=1)
            1
            sage: u.subs({u: 1, 'u': 2})
            Traceback (most recent call last):
            ...
            ValueError: Cannot substitute in u: duplicate key u.
            sage: u.subs({u: 1}, u=3)
            Traceback (most recent call last):
            ...
            ValueError: Cannot substitute in u: duplicate key u.
        """
        # check if nothing to do
        if not rules and not kwds:
            return self

        # init and process keyword arguments
        gens = self.parent().gens()
        locals = kwds or dict()

        # update with rules
        if isinstance(rules, dict):
            for k, v in rules.iteritems():
                if not isinstance(k, str) and k not in gens:
                    raise TypeError('Cannot substitute %s in %s '
                                    'since it is neither an '
                                    'asymptotic expansion '
                                    'nor a string (but a %s).' %
                                    (k, self, type(k)))
                k = str(k)
                if k in locals and locals[k] != v:
                    raise ValueError('Cannot substitute in %s: '
                                     'duplicate key %s.' % (self, k))
                locals[k] = v
        elif rules is not None:
            raise TypeError('Substitution rules %s have to be a dictionary.' %
                            (rules,))

        # fill up missing rules
        for g in gens:
            locals.setdefault(str(g), g)

        # check if all keys are generators
        gens_str = tuple(str(g) for g in gens)
        for k in locals:
            if str(k) not in gens_str:
                raise ValueError('Cannot substitute %s in %s '
                                 'since it is not a generator of %s.' %
                                 (k, self, self.parent()))

        # determine 0 and 1
        if domain is None and \
               ('_zero_' not in locals or '_one_' not in locals):
            P = self.parent()
            for g in gens:
                G = locals[str(g)].parent()
                if G is not P:
                    domain = G
                    break
            else:
                domain = P
        locals.setdefault('_zero_', domain.zero())
        locals.setdefault('_one_', domain.one())

        # do the actual substitution
        try:
            return self._substitute_(locals)
        except (ArithmeticError, TypeError, ValueError) as e:
            from misc import combine_exceptions
            rules = '{' + ', '.join(
                '%s: %s' % (k, v)
                for k, v in sorted(locals.iteritems(),
                                   key=lambda k: str(k[0]))
                if not k.startswith('_') and
                not any(k == str(g) and v is g for g in gens)) + '}'
            raise combine_exceptions(
                TypeError('Cannot apply the substitution rules %s on %s '
                          'in %s.' % (rules, self, self.parent())), e)


    subs = substitute


    def _substitute_(self, rules):
        r"""
        Substitute the given ``rules`` in this asymptotic expansion.

        INPUT:

        - ``rules`` -- a dictionary.
          The neutral element of the asymptotic ring is replaced by the value
          to key ``'_zero_'``.

        OUTPUT:

        An object.

        TESTS::

            sage: A.<z> = AsymptoticRing(growth_group='z^QQ', coefficient_ring=QQ)
            sage: z._substitute_({'z': SR.var('a')})
            a
            sage: _.parent()
            Symbolic Ring
            sage: A(0)._substitute_({'_zero_': 'zero'})
            'zero'
            sage: (1/z)._substitute_({'z': 4})
            1/4
            sage: _.parent()
            Rational Field
            sage: (1/z)._substitute_({'z': 0})
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot substitute in z^(-1) in
            Asymptotic Ring <z^QQ> over Rational Field.
            > *previous* ZeroDivisionError: Cannot substitute in z^(-1) in
            Exact Term Monoid z^QQ with coefficients in Rational Field.
            >> *previous* ZeroDivisionError: Cannot substitute in z^(-1) in
            Growth Group z^QQ.
            >...> *previous* ZeroDivisionError: rational division by zero
        """
        if not self.summands:
            return rules['_zero_']
        from sage.symbolic.operators import add_vararg
        try:
            return add_vararg(
                *tuple(s._substitute_(rules)
                       for s in self.summands.elements_topological()))
        except (ArithmeticError, TypeError, ValueError) as e:
            from misc import substitute_raise_exception
            substitute_raise_exception(self, e)


    def symbolic_expression(self, R=None):
        r"""
        Return this asymptotic expansion as a symbolic expression.

        INPUT:

        - ``R`` -- (a subring of) the symbolic ring or ``None``.
          The output is will be an element of ``R``. If ``None``,
          then the symbolic ring is used.

        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: A.<x, y, z> = AsymptoticRing(growth_group='x^ZZ * y^QQ * log(y)^QQ * QQ^z * z^QQ', coefficient_ring=QQ)
            sage: SR(A.an_element())  # indirect doctest
            1/8*(1/8)^z*x^3*y^(3/2)*z^(3/2)*log(y)^(3/2) +
            Order((1/2)^z*x*sqrt(y)*sqrt(z)*sqrt(log(y)))

        TESTS::

            sage: a = A.an_element(); a
            1/8*x^3*y^(3/2)*log(y)^(3/2)*(1/8)^z*z^(3/2) +
            O(x*y^(1/2)*log(y)^(1/2)*(1/2)^z*z^(1/2))
            sage: a.symbolic_expression()
            1/8*(1/8)^z*x^3*y^(3/2)*z^(3/2)*log(y)^(3/2) +
            Order((1/2)^z*x*sqrt(y)*sqrt(z)*sqrt(log(y)))
            sage: _.parent()
            Symbolic Ring

        ::

            sage: from sage.symbolic.ring import SymbolicRing
            sage: class MySymbolicRing(SymbolicRing):
            ....:     pass
            sage: mySR = MySymbolicRing()
            sage: a.symbolic_expression(mySR).parent() is mySR
            True
        """
        if R is None:
            from sage.symbolic.ring import SR
            R = SR

        return self.substitute(dict((g, R(R.var(str(g))))
                                    for g in self.parent().gens()),
                               domain=R)


    _symbolic_ = symbolic_expression  # will be used by SR._element_constructor_


class AsymptoticRing(Algebra, UniqueRepresentation):
    r"""
    A ring consisting of :class:`asymptotic expansions <AsymptoticExpansion>`.

    INPUT:

    - ``growth_group`` -- either a partially ordered group (see
      :doc:`growth_group`) or a string
      describing such a growth group (see
      :class:`~sage.rings.asymptotic.growth_group.GrowthGroupFactory`).

    - ``coefficient_ring`` -- the ring which contains the
      coefficients of the expansions.

    - ``default_prec`` -- a positive integer. This is the number of
      summands that are kept before truncating an infinite series.

    - ``category`` -- the category of the parent can be specified
      in order to broaden the base structure. It has to be a
      subcategory of ``Category of rings``. This is also the default
      category if ``None`` is specified.

    EXAMPLES:

    We begin with the construction of an asymptotic ring in various
    ways. First, we simply pass a string specifying the underlying
    growth group::

        sage: R1_x.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); R1_x
        Asymptotic Ring <x^QQ> over Rational Field
        sage: x
        x

    This is equivalent to the following code, which explicitly
    specifies the underlying growth group::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: G_QQ = GrowthGroup('x^QQ')
        sage: R2_x.<x> = AsymptoticRing(growth_group=G_QQ, coefficient_ring=QQ); R2_x
        Asymptotic Ring <x^QQ> over Rational Field

    Of course, the coefficient ring of the asymptotic ring and the
    base ring of the underlying growth group do not need to
    coincide::

        sage: R_ZZ_x.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=ZZ); R_ZZ_x
        Asymptotic Ring <x^QQ> over Integer Ring

    Note, we can also create and use logarithmic growth groups::

        sage: R_log = AsymptoticRing(growth_group='log(x)^ZZ', coefficient_ring=QQ); R_log
        Asymptotic Ring <log(x)^ZZ> over Rational Field

    Other growth groups are available. See :doc:`asymptotic_ring` for
    more examples.

    Below there are some technical details.

    According to the conventions for parents, uniqueness is ensured::

        sage: R1_x is R2_x
        True

    Furthermore, the coercion framework is also involved. Coercion
    between two asymptotic rings is possible (given that the
    underlying growth groups and coefficient rings are chosen
    appropriately)::

        sage: R1_x.has_coerce_map_from(R_ZZ_x)
        True

    Additionally, for the sake of convenience, the coefficient ring
    also coerces into the asymptotic ring (representing constant
    quantities)::

        sage: R1_x.has_coerce_map_from(QQ)
        True

    TESTS::

        sage: from sage.rings.asymptotic.asymptotic_ring import AsymptoticRing as AR_class
        sage: class AR(AR_class):
        ....:     class Element(AR_class.Element):
        ....:         __eq__ = AR_class.Element.has_same_summands
        sage: A = AR(growth_group='z^QQ', coefficient_ring=QQ)
        sage: from itertools import islice
        sage: TestSuite(A).run(  # not tested  # long time  # see #19424
        ....:     verbose=True,
        ....:     elements=tuple(islice(A.some_elements(), 10)),
        ....:     skip=('_test_some_elements',  # to many elements
        ....:           '_test_distributivity'))  # due to cancellations: O(z) != O(z^2)
    """

    # enable the category framework for elements
    Element = AsymptoticExpansion


    @staticmethod
    def __classcall__(cls, growth_group=None, coefficient_ring=None,
                      names=None, category=None, default_prec=None):
        r"""
        Normalizes the input in order to ensure a unique
        representation of the parent.

        For more information see :class:`AsymptoticRing`.

        EXAMPLES:

        ``__classcall__`` unifies the input to the constructor of
        :class:`AsymptoticRing` such that the instances generated
        are unique. Also, this enables the use of the generation
        framework::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: MG = GrowthGroup('x^ZZ')
            sage: AR1 = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: AR2.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR1 is AR2
            True

        The bracket notation can only be used if the growth group
        has a generator::

            sage: AR.<lx> = AsymptoticRing(growth_group='log(x)^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError:  Growth Group log(x)^ZZ does not provide any
            generators but name 'lx' given.

        The names of the generators have to agree with the names used in
        the growth group except for univariate rings::

            sage: A.<icecream> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ); A
            Asymptotic Ring <x^ZZ> over Integer Ring
            sage: icecream
            x
            sage: A.<x, y> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ); A
            Asymptotic Ring <x^ZZ * y^ZZ> over Integer Ring
            sage: A.<y, x> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Names 'y', 'x' do not coincide with
            generators 'x', 'y' of Growth Group x^ZZ * y^ZZ.
            sage: A.<a, b> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Names 'a', 'b' do not coincide with
            generators 'x', 'y' of Growth Group x^ZZ * y^ZZ.
            sage: A.<x, b> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Names 'x', 'b' do not coincide with
            generators 'x', 'y' of Growth Group x^ZZ * y^ZZ.
            sage: A.<x> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Name 'x' do not coincide with
            generators 'x', 'y' of Growth Group x^ZZ * y^ZZ.
            sage: A.<x, y, z> = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Names 'x', 'y', 'z' do not coincide with
            generators 'x', 'y' of Growth Group x^ZZ * y^ZZ.

        TESTS::

            sage: AsymptoticRing(growth_group=None, coefficient_ring=ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Growth group not specified. Cannot continue.
            sage: AsymptoticRing(growth_group='x^ZZ', coefficient_ring=None)
            Traceback (most recent call last):
            ...
            ValueError: Coefficient ring not specified. Cannot continue.
            sage: AsymptoticRing(growth_group='x^ZZ', coefficient_ring='icecream')
            Traceback (most recent call last):
            ...
            ValueError: icecream is not a ring. Cannot continue.
        """
        from sage.categories.sets_cat import Sets
        from sage.categories.rings import Rings

        Sets_parent_class = Sets().parent_class
        while issubclass(cls, Sets_parent_class):
            cls = cls.__base__

        if isinstance(growth_group, str):
            from growth_group import GrowthGroup
            growth_group = GrowthGroup(growth_group)

        if growth_group is None:
            raise ValueError('Growth group not specified. Cannot continue.')

        if coefficient_ring is None:
            raise ValueError('Coefficient ring not specified. Cannot continue.')
        if coefficient_ring not in Rings():
            raise ValueError('%s is not a ring. Cannot continue.' % (coefficient_ring,))

        strgens = tuple(str(g) for g in growth_group.gens_monomial())
        def format_names(N):
            return ('s ' if len(N) != 1 else ' ') + ', '.join("'%s'" % n for n in N)
        if names and not strgens:
            raise ValueError('%s does not provide any generators but name%s given.' %
                             (growth_group, format_names(names)))
        elif names is not None and len(names) == 1 and len(strgens) == 1:
            pass
        elif names is not None and names != strgens:
            raise ValueError('Name%s do not coincide with generator%s of %s.' %
                             (format_names(names), format_names(strgens), growth_group))

        if category is None:
            from sage.categories.commutative_algebras import CommutativeAlgebras
            from sage.categories.rings import Rings
            category = CommutativeAlgebras(Rings())

        if default_prec is None:
            from sage.misc.defaults import series_precision
            default_prec = series_precision()

        return super(AsymptoticRing,
                     cls).__classcall__(cls, growth_group, coefficient_ring,
                                        category=category,
                                        default_prec=default_prec)


    @experimental(trac_number=17601)
    def __init__(self, growth_group, coefficient_ring, category, default_prec):
        r"""
        See :class:`AsymptoticRing` for more information.

        TESTS::

            sage: R1 = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ); R1
            Asymptotic Ring <x^ZZ> over Integer Ring
            sage: R2.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); R2
            Asymptotic Ring <x^QQ> over Rational Field
            sage: R1 is R2
            False

        ::

            sage: R3 = AsymptoticRing('x^ZZ')
            Traceback (most recent call last):
            ...
            ValueError: Coefficient ring not specified. Cannot continue.
        """
        self._coefficient_ring_ = coefficient_ring
        self._growth_group_ = growth_group
        self._default_prec_ = default_prec
        super(AsymptoticRing, self).__init__(base_ring=coefficient_ring,
                                             category=category)


    @property
    def growth_group(self):
        r"""
        The growth group of this asymptotic ring.

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.growth_group
            Growth Group x^ZZ

        .. SEEALSO::

            :doc:`growth_group`
        """
        return self._growth_group_


    @property
    def coefficient_ring(self):
        r"""
        The coefficient ring of this asymptotic ring.

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.coefficient_ring
            Integer Ring
        """
        return self._coefficient_ring_


    @property
    def default_prec(self):
        r"""
        The default precision of this asymptotic ring.

        This is the parameter used to determine how many summands
        are kept before truncating an infinite series (which occur
        when inverting asymptotic expansions).

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.default_prec
            20
            sage: AR = AsymptoticRing('x^ZZ', ZZ, default_prec=123)
            sage: AR.default_prec
            123
        """
        return self._default_prec_


    def change_parameter(self, **kwds):
        r"""
        Return an asymptotic ring with a change in one or more of the given parameters.

        INPUT:

        - ``growth_group`` -- (default: ``None``) the new growth group.

        - ``coefficient_ring`` -- (default: ``None``) the new coefficient ring.

        - ``category`` -- (default: ``None``) the new category.

        - ``default_prec`` -- (default: ``None``) the new default precision.

        OUTPUT:

        An asymptotic ring.

        EXAMPLES::

            sage: A = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: A.change_parameter(coefficient_ring=QQ)
            Asymptotic Ring <x^ZZ> over Rational Field

        TESTS::

            sage: A.change_parameter(coefficient_ring=ZZ) is A
            True
        """
        parameters = ('growth_group', 'coefficient_ring', 'default_prec')
        values = dict()
        for parameter in parameters:
            values[parameter] = kwds.get(parameter, getattr(self, parameter))
        values['category'] = self.category()
        if isinstance(values['growth_group'], str):
            from growth_group import GrowthGroup
            values['growth_group'] = GrowthGroup(values['growth_group'])
        if all(values[parameter] is getattr(self, parameter)
               for parameter in parameters) and values['category'] is self.category():
            return self
        from misc import underlying_class
        return underlying_class(self)(**values)


    @staticmethod
    def _create_empty_summands_():
        r"""
        Create an empty data structure suitable for storing and working
        with summands.

        INPUT:

        Nothing.

        OUTPUT:

        A :class:`~sage.data_structures.mutable_poset.MutablePoset`.

        TESTS::

            sage: AsymptoticRing._create_empty_summands_()
            poset()
        """
        from sage.data_structures.mutable_poset import MutablePoset
        from term_monoid import can_absorb, absorption
        return MutablePoset(key=lambda element: element.growth,
                            can_merge=can_absorb,
                            merge=absorption)


    def _create_element_in_extension_(self, term, old_term_parent=None):
        r"""
        Create an element in an extension of this asymptotic ring which
        is chosen according to the input.

        INPUT:

        - ``term`` -- the element data.

        - ``old_term_parent`` -- the parent of ``term`` is compared to this
          parent. If both are the same or ``old_parent`` is ``None``,
          then the result is an expansion in this (``self``) asymptotic ring.

        OUTPUT:

        An element.

        EXAMPLES::

            sage: A = AsymptoticRing('z^ZZ', ZZ)
            sage: a = next(A.an_element().summands.elements_topological())
            sage: B = AsymptoticRing('z^QQ', QQ)
            sage: b = next(B.an_element().summands.elements_topological())
            sage: c = A._create_element_in_extension_(a, a.parent())
            sage: next(c.summands.elements_topological()).parent()
            O-Term Monoid z^ZZ with implicit coefficients in Integer Ring
            sage: c = A._create_element_in_extension_(b, a.parent())
            sage: next(c.summands.elements_topological()).parent()
            O-Term Monoid z^QQ with implicit coefficients in Rational Field

        TESTS::

            sage: c = A._create_element_in_extension_(b, None)
            sage: next(c.summands.elements_topological()).parent()
            O-Term Monoid z^QQ with implicit coefficients in Rational Field
        """
        if old_term_parent is None or term.parent() is old_term_parent:
            parent = self
        else:
            # Insert an 'if' here once terms can have different
            # coefficient rings, as this will be for L-terms.
            parent = self.change_parameter(
                growth_group=term.parent().growth_group,
                coefficient_ring=term.parent().coefficient_ring)
        return parent(term, simplify=False, convert=False)


    def _element_constructor_(self, data, simplify=True, convert=True):
        r"""
        Convert a given object to this asymptotic ring.

        INPUT:

        - ``data`` -- an object representing the element to be
          initialized.

        - ``simplify`` -- (default: ``True``) if set, then the constructed
          element is simplified (terms are absorbed) automatically.

        - ``convert`` -- (default: ``True``) passed on to the element
          constructor.  If set, then the ``summands`` are converted to
          the asymptotic ring (the parent of this expansion). If not,
          then the summands are taken as they are. In that case, the
          caller must ensure that the parent of the terms is set
          correctly.

        OUTPUT:

        An element of this asymptotic ring.

        TESTS::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR(5)  # indirect doctest
            5
            sage: AR(3*x^2)  # indirect doctest
            3*x^2
            sage: x = ZZ['x'].gen(); x.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: AR(x)
            x
            sage: y = ZZ['y'].gen(); AR(y)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Polynomial y is not in
            Asymptotic Ring <x^ZZ> over Integer Ring
            > *previous* ValueError: Growth y is not in
            Exact Term Monoid x^ZZ with coefficients in Integer Ring.
            >> *previous* ValueError: y is not in Growth Group x^ZZ.

        ::

            sage: A = AsymptoticRing(growth_group='p^ZZ', coefficient_ring=QQ)
            sage: P.<p> = QQ[]
            sage: A(p)  # indirect doctest
            p
            sage: A(p^11)  # indirect doctest
            p^11
            sage: A(2*p^11)  # indirect doctest
            2*p^11
            sage: A(3*p^4 + 7/3*p - 8)  # indirect doctest
            3*p^4 + 7/3*p - 8

        ::

            sage: S = AsymptoticRing(growth_group='x^ZZ * y^ZZ', coefficient_ring=QQ)
            sage: var('x, y')
            (x, y)
            sage: S(x + y)  # indirect doctest
            x + y
            sage: S(2*x - 4*x*y^6)  # indirect doctest
            -4*x*y^6 + 2*x

        ::

            sage: A.<a,b> = AsymptoticRing('a^ZZ * b^ZZ', QQ)
            sage: 1/a
            a^(-1)

        ::

            sage: P.<a, b, c> = ZZ[]
            sage: A(a + b)
            a + b
            sage: A(a + c)
            Traceback (most recent call last):
            ...
            ValueError: Polynomial a + c is not in
            Asymptotic Ring <a^ZZ * b^ZZ> over Rational Field
            > *previous* ValueError: Growth c is not in
            Exact Term Monoid a^ZZ * b^ZZ with coefficients in Rational Field.
            >> *previous* ValueError: c is not in Growth Group a^ZZ * b^ZZ.
            >...> *previous* ValueError: c is not in any of the factors of
            Growth Group a^ZZ * b^ZZ

        ::

            sage: M = AsymptoticRing('m^ZZ', ZZ)
            sage: N = AsymptoticRing('n^ZZ', QQ)
            sage: N(M.an_element())  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot include m^3 with parent
            Exact Term Monoid m^ZZ with coefficients in Integer Ring
            in Asymptotic Ring <n^ZZ> over Rational Field
            > *previous* ValueError: m^3 is not in Growth Group n^ZZ

        ::

            sage: M([1])  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: Not all list entries of [1] are asymptotic terms,
            so cannot create an asymptotic expansion in
            Asymptotic Ring <m^ZZ> over Integer Ring.
            sage: M(SR.var('a') + 1)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Symbolic expression a + 1 is not in
            Asymptotic Ring <m^ZZ> over Integer Ring.
            > *previous* ValueError: a is not in
            Exact Term Monoid m^ZZ with coefficients in Integer Ring.
            >> *previous* ValueError: Factor a of a is neither a coefficient
            (in Integer Ring) nor growth (in Growth Group m^ZZ).
        """
        from sage.data_structures.mutable_poset import MutablePoset
        if isinstance(data, MutablePoset):
            return self.element_class(self, data, simplify=simplify, convert=convert)

        if type(data) == self.element_class and data.parent() == self:
            return data

        if isinstance(data, AsymptoticExpansion):
            return self.element_class(self, data.summands,
                                      simplify=simplify, convert=convert)

        from term_monoid import GenericTerm
        if isinstance(data, GenericTerm):
            data = (data,)

        if isinstance(data, (list, tuple)):
            if not all(isinstance(elem, GenericTerm) for elem in data):
                raise TypeError('Not all list entries of %s '
                                'are asymptotic terms, so cannot create an '
                                'asymptotic expansion in %s.' % (data, self))
            summands = AsymptoticRing._create_empty_summands_()
            summands.union_update(data)
            return self.element_class(self, summands,
                                      simplify=simplify, convert=convert)

        if not data:
            summands = AsymptoticRing._create_empty_summands_()
            return self.element_class(self, summands,
                                      simplify=simplify, convert=False)

        try:
            P = data.parent()
        except AttributeError:
            return self.create_summand('exact', data)

        from misc import combine_exceptions
        from sage.symbolic.ring import SymbolicRing
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        from sage.rings.polynomial.multi_polynomial_ring_generic import is_MPolynomialRing
        from sage.rings.power_series_ring import is_PowerSeriesRing

        if isinstance(P, SymbolicRing):
            from sage.symbolic.operators import add_vararg
            if data.operator() == add_vararg:
                summands = []
                for summand in data.operands():
                    # TODO: check if summand is an O-Term here
                    # (see #19425, #19426)
                    try:
                        summands.append(self.create_summand('exact', summand))
                    except ValueError as e:
                        raise combine_exceptions(
                            ValueError('Symbolic expression %s is not in %s.' %
                                       (data, self)), e)
                return sum(summands)

        elif is_PolynomialRing(P):
            p = P.gen()
            try:
                return sum(self.create_summand('exact', growth=p**i,
                                               coefficient=c)
                           for i, c in enumerate(data))
            except ValueError as e:
                raise combine_exceptions(
                    ValueError('Polynomial %s is not in %s' % (data, self)), e)

        elif is_MPolynomialRing(P):
            try:
                return sum(self.create_summand('exact', growth=g, coefficient=c)
                           for c, g in iter(data))
            except ValueError as e:
                raise combine_exceptions(
                    ValueError('Polynomial %s is not in %s' % (data, self)), e)

        elif is_PowerSeriesRing(P):
            raise NotImplementedError(
                'Cannot convert %s from the %s to an asymptotic expansion '
                'in %s, since growths at other points than +oo are not yet '
                'supported.' % (data, P, self))
            # Delete lines above as soon as we can deal with growths
            # other than the that at going to +oo.
            p = P.gen()
            try:
                result = self(data.polynomial())
            except ValueError as e:
                raise combine_exceptions(
                    ValueError('Powerseries %s is not in %s' % (data, self)), e)
            prec = data.precision_absolute()
            if prec < sage.rings.infinity.PlusInfinity():
                try:
                    result += self.create_summand('O', growth=p**prec)
                except ValueError as e:
                    raise combine_exceptions(
                        ValueError('Powerseries %s is not in %s' %
                                   (data, self)), e)
            return result

        return self.create_summand('exact', data)


    def _coerce_map_from_(self, R):
        r"""
        Return whether ``R`` coerces into this asymptotic ring.

        INPUT:

        - ``R`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            There are two possible cases: either ``R`` coerces in the
            :meth:`coefficient_ring` of this asymptotic ring, or ``R``
            itself is an asymptotic ring, where both the
            :meth:`growth_group` and the :meth:`coefficient_ring` coerce into
            the :meth:`growth_group` and the :meth:`coefficient_ring` of this
            asymptotic ring, respectively.

        TESTS::

            sage: AR_ZZ = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ); AR_ZZ
            Asymptotic Ring <x^ZZ> over Integer Ring
            sage: x_ZZ = AR_ZZ.gen()
            sage: AR_QQ = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); AR_QQ
            Asymptotic Ring <x^QQ> over Rational Field
            sage: x_QQ = AR_QQ.gen()
            sage: AR_QQ.has_coerce_map_from(AR_ZZ)  # indirect doctest
            True
            sage: x_ZZ * x_QQ
            x^2

        ::

            sage: AR_QQ.has_coerce_map_from(QQ)
            True
            sage: AR_QQ.has_coerce_map_from(ZZ)
            True
        """
        from sage.data_structures.mutable_poset import MutablePoset
        if R == MutablePoset:
            return
        if self.coefficient_ring.has_coerce_map_from(R):
            return True
        if self.growth_group.has_coerce_map_from(R):
            return True
        elif isinstance(R, AsymptoticRing):
            if self.growth_group.has_coerce_map_from(R.growth_group) and \
                    self.coefficient_ring.has_coerce_map_from(R.coefficient_ring):
                return True


    def _repr_(self):
        r"""
        A representation string of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ',
            ....:                     coefficient_ring=ZZ)
            sage: repr(AR)  # indirect doctest
            'Asymptotic Ring <x^ZZ> over Integer Ring'
        """
        try:
            G = '<' + self.growth_group._repr_(condense=True) + '>'
        except TypeError:
            G = repr(self.growth_group)
        return 'Asymptotic Ring %s over %s' % (G, self.coefficient_ring)


    def _an_element_(self):
        r"""
        Return an element of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        An :class:`AsymptoticExpansion`.

        EXAMPLES::

            sage: AsymptoticRing(growth_group='z^QQ', coefficient_ring=ZZ).an_element()
            z^(3/2) + O(z^(1/2))
            sage: AsymptoticRing(growth_group='z^ZZ', coefficient_ring=QQ).an_element()
            1/8*z^3 + O(z)
            sage: AsymptoticRing(growth_group='z^QQ', coefficient_ring=QQ).an_element()
            1/8*z^(3/2) + O(z^(1/2))
        """
        from term_monoid import TermMonoid
        E = TermMonoid('exact', asymptotic_ring=self)
        O = TermMonoid('O', asymptotic_ring=self)
        return self(E.an_element(), simplify=False, convert=False)**3 + \
            self(O.an_element(), simplify=False, convert=False)


    def some_elements(self):
        r"""
        Return some elements of this term monoid.

        See :class:`TestSuite` for a typical use case.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from itertools import islice
            sage: A = AsymptoticRing(growth_group='z^QQ', coefficient_ring=ZZ)
            sage: tuple(islice(A.some_elements(), 10))
            (z^(3/2) + O(z^(1/2)),
             O(z^(1/2)),
             z^(3/2) + O(z^(-1/2)),
             -z^(3/2) + O(z^(1/2)),
             O(z^(-1/2)),
             O(z^2),
             z^6 + O(z^(1/2)),
             -z^(3/2) + O(z^(-1/2)),
             O(z^2),
             z^(3/2) + O(z^(-2)))
        """
        from sage.misc.mrange import cantor_product
        from term_monoid import TermMonoid
        E = TermMonoid('exact', asymptotic_ring=self)
        O = TermMonoid('O', asymptotic_ring=self)
        return iter(self(e, simplify=False, convert=False)**3 +
                    self(o, simplify=False, convert=False)
                    for e, o in cantor_product(
                            E.some_elements(), O.some_elements()))


    def gens(self):
        r"""
        Return a tuple with generators of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple of asymptotic expansions.

        .. NOTE::

            Generators do not necessarily exist. This depends on the
            underlying growth group. For example,
            :class:`monomial growth groups <sage.rings.asymptotic.growth_group.MonomialGrowthGroup>`
            have a generator, and exponential growth groups
            do not.

        EXAMPLES::

            sage: AR.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.gens()
            (x,)
            sage: B.<y,z> = AsymptoticRing(growth_group='y^ZZ * z^ZZ', coefficient_ring=QQ)
            sage: B.gens()
            (y, z)
        """
        return tuple(self.create_summand('exact',
                                         growth=g,
                                         coefficient=self.coefficient_ring(1))
                     for g in self.growth_group.gens_monomial())


    def gen(self, n=0):
        r"""
        Return the ``n``-th generator of this asymptotic ring.

        INPUT:

        - ``n`` -- (default: `0`) a non-negative integer.

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: R.gen()
            x
        """
        return self.gens()[n]


    def ngens(self):
        r"""
        Return the number of generators of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: AR.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.ngens()
            1
        """
        return len(self.growth_group.gens_monomial())


    def create_summand(self, type, data=None, **kwds):
        r"""
        Create a simple asymptotic expansion consisting of a single
        summand.

        INPUT:

        - ``type`` -- 'O' or 'exact'.

        - ``data`` -- the element out of which a summand has to be created.

        - ``growth`` -- an element of the :meth:`growth_group`.

        - ``coefficient`` -- an element of the :meth:`coefficient_ring`.

        .. NOTE::

            Either ``growth`` and ``coefficient`` or ``data`` have to
            be specified.

        OUTPUT:

        An asymptotic expansion.

        .. NOTE::

            This method calls the factory :class:`TermMonoid
            <sage.rings.asymptotic.term_monoid.TermMonoidFactory>`
            with the appropriate arguments.

        EXAMPLES::

            sage: R = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: R.create_summand('O', x^2)
            O(x^2)
            sage: R.create_summand('exact', growth=x^456, coefficient=123)
            123*x^456
            sage: R.create_summand('exact', data=12*x^13)
            12*x^13

        TESTS::

            sage: R.create_summand('exact', data='12*x^13')
            12*x^13
            sage: R.create_summand('exact', data='x^13 * 12')
            12*x^13
            sage: R.create_summand('exact', data='x^13')
            x^13
            sage: R.create_summand('exact', data='12')
            12
            sage: R.create_summand('exact', data=12)
            12

        ::

            sage: R.create_summand('O', growth=42*x^2, coefficient=1)
            Traceback (most recent call last):
            ...
            ValueError: Growth 42*x^2 is not in O-Term Monoid x^ZZ with implicit coefficients in Integer Ring.
            > *previous* ValueError: 42*x^2 is not in Growth Group x^ZZ.

        ::

            sage: AR.<z> = AsymptoticRing('z^QQ', QQ)
            sage: AR.create_summand('exact', growth='z^2')
            Traceback (most recent call last):
            ...
            TypeError: Cannot create exact term: only 'growth' but
            no 'coefficient' specified.
        """
        from term_monoid import TermMonoid
        TM = TermMonoid(type, asymptotic_ring=self)

        if data is None:
            try:
                data = kwds.pop('growth')
            except KeyError:
                raise TypeError("Neither 'data' nor 'growth' are specified.")
            if type == 'exact' and kwds.get('coefficient') is None:
                raise TypeError("Cannot create exact term: only 'growth' "
                                "but no 'coefficient' specified.")


        if type == 'exact' and kwds.get('coefficient') == 0:
            return self.zero()

        return self(TM(data, **kwds), simplify=False, convert=False)


    def variable_names(self):
        r"""
        Return the names of the variables.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: A = AsymptoticRing(growth_group='x^ZZ * QQ^y', coefficient_ring=QQ)
            sage: A.variable_names()
            ('x', 'y')
        """
        return self.growth_group.variable_names()


    def construction(self):
        r"""
        Return the construction of this asymptotic ring.

        OUTPUT:

        A pair whose first entry is an
        :class:`asymptotic ring construction functor <AsymptoticRingFunctor>`
        and its second entry the coefficient ring.

        EXAMPLES::

            sage: A = AsymptoticRing(growth_group='x^ZZ * QQ^y', coefficient_ring=QQ)
            sage: A.construction()
            (AsymptoticRing<x^ZZ * QQ^y>, Rational Field)

        .. SEEALSO::

            :doc:`asymptotic_ring`,
            :class:`AsymptoticRing`,
            :class:`AsymptoticRingFunctor`.
        """
        return AsymptoticRingFunctor(self.growth_group), self.coefficient_ring


from sage.categories.pushout import ConstructionFunctor
class AsymptoticRingFunctor(ConstructionFunctor):
    r"""
    A :class:`construction functor <sage.categories.pushout.ConstructionFunctor>`
    for :class:`asymptotic rings <AsymptoticRing>`.

    INPUT:

    - ``growth_group`` -- a partially ordered group (see
      :class:`AsymptoticRing` or
      :doc:`growth_group` for details).

    EXAMPLES::

        sage: AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ).construction()  # indirect doctest
        (AsymptoticRing<x^ZZ>, Rational Field)

    .. SEEALSO::

        :doc:`asymptotic_ring`,
        :class:`AsymptoticRing`,
        :class:`sage.rings.asymptotic.growth_group.AbstractGrowthGroupFunctor`,
        :class:`sage.rings.asymptotic.growth_group.ExponentialGrowthGroupFunctor`,
        :class:`sage.rings.asymptotic.growth_group.MonomialGrowthGroupFunctor`,
        :class:`sage.categories.pushout.ConstructionFunctor`.

    TESTS::

        sage: X = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
        sage: Y = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=QQ)
        sage: cm = sage.structure.element.get_coercion_model()
        sage: cm.record_exceptions()
        sage: cm.common_parent(X, Y)
        Asymptotic Ring <x^ZZ * y^ZZ> over Rational Field
        sage: sage.structure.element.coercion_traceback()  # not tested

    ::

        sage: from sage.categories.pushout import pushout
        sage: pushout(AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ), QQ)
        Asymptotic Ring <x^ZZ> over Rational Field
    """

    rank = 13


    def __init__(self, growth_group):
        r"""
        See :class:`AsymptoticRingFunctor` for details.

        TESTS::

            sage: from sage.rings.asymptotic.asymptotic_ring import AsymptoticRingFunctor
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: AsymptoticRingFunctor(GrowthGroup('x^ZZ'))
            AsymptoticRing<x^ZZ>
        """
        self.growth_group = growth_group

        from sage.categories.rings import Rings
        super(ConstructionFunctor, self).__init__(
            Rings(), Rings())


    def _repr_(self):
        r"""
        Return a representation string of this functor.

        OUTPUT:

        A string.

        TESTS::

            sage: AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ).construction()[0]  # indirect doctest
            AsymptoticRing<x^ZZ>
        """
        return 'AsymptoticRing<%s>' % (self.growth_group._repr_(condense=True),)


    def _apply_functor(self, coefficient_ring):
        r"""
        Apply this functor to the given ``coefficient_ring``.

        INPUT:

        - ``base`` - anything :class:`~sage.rings.asymptotic.growth_group.MonomialGrowthGroup` accepts.

        OUTPUT:

        An :class:`AsymptoticRing`.

        EXAMPLES::

            sage: A = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
            sage: F, C = A.construction()
            sage: F(C)  # indirect doctest
            Asymptotic Ring <x^ZZ> over Rational Field
        """
        return AsymptoticRing(growth_group=self.growth_group,
                              coefficient_ring=coefficient_ring)


    def merge(self, other):
        r"""
        Merge this functor with ``other`` if possible.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A functor or ``None``.

        EXAMPLES::

            sage: X = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
            sage: Y = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=QQ)
            sage: F_X = X.construction()[0]
            sage: F_Y = Y.construction()[0]
            sage: F_X.merge(F_X)
            AsymptoticRing<x^ZZ>
            sage: F_X.merge(F_Y)
            AsymptoticRing<x^ZZ * y^ZZ>
        """
        if self == other:
            return self

        if isinstance(other, AsymptoticRingFunctor):
            from sage.structure.element import get_coercion_model
            cm = get_coercion_model()
            try:
                G = cm.common_parent(self.growth_group, other.growth_group)
            except TypeError:
                pass
            else:
                return AsymptoticRingFunctor(G)


    def __eq__(self, other):
        r"""
        Return whether this functor is equal to ``other``.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: X = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
            sage: Y = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=QQ)
            sage: F_X = X.construction()[0]
            sage: F_Y = Y.construction()[0]
            sage: F_X == F_X
            True
            sage: F_X == F_Y
            False
        """
        return type(self) == type(other) and \
            self.growth_group == other.growth_group


    def __ne__(self, other):
        r"""
        Return whether this functor is not equal to ``other``.

        INPUT:

        - ``other`` -- a functor.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: X = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=QQ)
            sage: Y = AsymptoticRing(growth_group='y^ZZ', coefficient_ring=QQ)
            sage: F_X = X.construction()[0]
            sage: F_Y = Y.construction()[0]
            sage: F_X != F_X
            False
            sage: F_X != F_Y
            True
        """
        return not self.__eq__(other)
