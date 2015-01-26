r"""
Mutable Poset

This module provides a class representing a finite partially ordered
set (poset) for the purpose of being used as a data structure. Thus
the here introduced posets are mutable, i.e., elements can be added and
removed from a poset at any time.

To get in touch with Sage's "usual" posets, start with the page
:mod:`Posets <sage.combinat.posets.__init__>` in the reference manual.


.. _mutable_poset_intro:

Introduction
============


.. _mutable_poset_examples:

Examples
========

First Steps
-----------

We start by creating an empty poset. This is simply done by

::

    sage: from sage.data_structures.mutable_poset import MutablePoset as MP
    sage: P = MP()
    sage: P
    poset()

A poset should contain elements, thus let us add them with

::

    sage: P.add(42)
    sage: P.add(7)
    sage: P.add(13)
    sage: P.add(3)

Let us look at the poset again::

    sage: P
    poset(3, 7, 13, 42)

We see that they elements are sorted using `\leq` which exists on the
integers `\ZZ`. Since this is even a total order, we could have used
the more efficient :class:`MutableToset`. (Note that also other data
structures are suitable for totally ordered sets.)


A less boring Example
---------------------

Let us continue with a less boring example. We define the class

::

    sage: class T(tuple):
    ....:     def __le__(left, right):
    ....:         return all(l <= r for l, r in zip(left, right))

It is equipped with a `\leq`-operation which makes `a \leq b` if all
entries of `a` are at most `b`. For example, we have

::

    sage: a = T((1,1))
    sage: b = T((2,1))
    sage: c = T((1,2))
    sage: a <= b, a <= c, b <= c
    (True, True, False)

The last comparison gives ``False``, since the first entries give `2 \leq 1`.

Now, let us add such elements to a poset::

    sage: Q = MP()
    sage: Q.add(T((1, 1)))
    sage: Q.add(T((3, 3)))
    sage: Q.add(T((4, 1)))
    sage: Q.add(T((3, 2)))
    sage: Q.add(T((2, 3)))
    sage: Q.add(T((2, 2)))
    sage: Q
    poset((1, 1), (2, 2), (2, 3), (3, 2), (3, 3), (4, 1))

In the representation above, the elements are sorted topologically,
smallest first. This does not show (directly) more structural
information. We can overcome this and display a "wiring layout" by
typing::

    sage: print Q.repr_full(reverse=True)
    poset((3, 3), (2, 3), (3, 2), (2, 2), (4, 1), (1, 1))
    +-- oo
    |   +-- no successors
    |   +-- predecessors:   (3, 3), (4, 1)
    +-- (3, 3)
    |   +-- successors:   oo
    |   +-- predecessors:   (2, 3), (3, 2)
    +-- (2, 3)
    |   +-- successors:   (3, 3)
    |   +-- predecessors:   (2, 2)
    +-- (3, 2)
    |   +-- successors:   (3, 3)
    |   +-- predecessors:   (2, 2)
    +-- (2, 2)
    |   +-- successors:   (2, 3), (3, 2)
    |   +-- predecessors:   (1, 1)
    +-- (4, 1)
    |   +-- successors:   oo
    |   +-- predecessors:   (1, 1)
    +-- (1, 1)
    |   +-- successors:   (2, 2), (4, 1)
    |   +-- predecessors:   null
    +-- null
    |   +-- successors:   (1, 1)
    |   +-- no predecessors

Note that we use ``reverse=True`` to let the elements appear from
largest (on the top) to smallest (on the bottom).

If you look at the output above, you'll see two additional elements,
namely ``oo`` (`\infty`) and ``null`` (`\emptyset`). So what are these
strange animals? The answer is simple and maybe you can guess it
already. The `\infty`-element is larger than every other element,
therefore a successor of the maximal elements in the poset. Similarly,
the `\emptyset`-element is smaller than any other element, therefore a
predecessor of the poset's minimal elements. Both do not have to scare
us; they are just there and sometimes useful.


AUTHORS:

- Daniel Krenn (2015-01-21): initial version

ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the Austrian Science Fund (FWF): P 24644-N26.

Classes and their Methods
=========================
"""
#*****************************************************************************
# Copyright (C) 2015 Daniel Krenn <dev@danielkrenn.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************

import sage


# *****************************************************************************


class MutablePosetElement(sage.structure.sage_object.SageObject):
    r"""
    An element of a mutable poset.
    """
    def __init__(self, poset, value):
        r"""
        See :class:`MutablePosetElement` for details.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: MutablePosetElement(P, (1, 2))
            (1, 2)
        """
        self._poset_ = poset
        self._value_ = value
        self._predecessors_ = set()
        self._successors_ = set()


    @property
    def poset(self):
        r"""
        The poset to which the element belongs.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: e.poset is P
            True
        """
        return self._poset_


    @property
    def value(self):
        r"""
        The value of the element.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: e.value
            (1, 2)
        """
        return self._value_


    @property
    def key(self):
        r"""
        The key of the element.

        The value of the element is converted by the poset to the key.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: P = MP()
            sage: e = MutablePosetElement(P, (1, 2))
            sage: e.key
            (1, 2)
            sage: Q = MP(key=lambda k: k[0])
            sage: f = MutablePosetElement(Q, (1, 2))
            sage: f.key
            1
        """
        return self.poset.get_key(self._value_)


    def predecessors(self, reverse=False):
        r"""
        Return the predecessors of the element.

        INPUT:

        - ``reverse`` -- (default: ``False``) if set, then returns
          successors instead.

        OUTPUT:

        A set.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: e.predecessors()
            set()
        """
        if reverse:
            return self._successors_
        return self._predecessors_


    def successors(self, reverse=False):
        r"""
        Return the successors of the element.

        INPUT:

        - ``reverse`` -- (default: ``False``) if set, then returns
          predecessors instead.

        OUTPUT:

        A set.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: e.successors()
            set()
        """
        if reverse:
            return self._predecessors_
        return self._successors_


    def is_special(self):
        r"""

        Return if the element is either the null-element, i.e., the
        element smaller than any possible other element or the
        infinity-element, i.e., the element larger than any possible
        other element.

        INPUT:

        Nothing.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.null.is_special()
            True
            sage: P.oo.is_special()
            True
        """
        return self.value is None


    def is_null(self, reverse=False):
        r"""
        Return if the element is the null-element, i.e., the element
        smaller than any possible other element.

        INPUT:

        - ``reverse`` -- (default: ``False``) if set, then returns
          if the element is the largest possible.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.null.is_null()
            True
            sage: P.oo.is_null()
            False
        """
        return self.value is None and not self.predecessors(reverse)


    def is_oo(self, reverse=False):
        r"""
        Return if the element is the infinity-element, i.e., the element
        larger than any possible other element.

        INPUT:

        - ``reverse`` -- (default: ``False``) if set, then returns
          if the element is the smallest possible.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.null.is_oo()
            False
            sage: P.oo.is_oo()
            True
        """
        return self.value is None and not self.successors(reverse)


    def __repr__(self):
        r"""
        Return the representation of the element.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        This methods usually returns the representation string of its
        :meth:`value`. The only exception is if this value is
        ``None``. In this case either ``'null'`` or ``'oo'`` is
        returned depending in the nonexistence of predecessors and
        sucessors respectively.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: repr(MutablePosetElement(P, (1, 2)))  # indirect doctest
            '(1, 2)'
            sage: repr(P.null)  # indirect doctest
            'null'
            sage: repr(P.oo)  # indirect doctest
            'oo'
        """
        if self.value is None:
            if not self.predecessors():
                return 'null'
            if not self.successors():
                return 'oo'
        return repr(self.value)


    def __hash__(self):
        r"""
        Return the hash of the element.

        INPUT:

        Nothing.

        OUTPUT:

        A hash value.

        This returns the hash value of the key of this element.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: hash(MutablePosetElement(P, (1, 2))) == hash((1, 2))
            True
        """
        return hash(self.key)


    def le(left, right, reverse=False):
        r"""
        Return if ``left`` is less or equal to ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        - ``reverse`` -- (default: ``False``) if set, then return
          ``right <= left`` instead.

        OUTPUT:

        ``True`` or ``False``.

        This methods usually returns if the keys of the given
        elements are less or equal. The only exception is if the
        value is ``None``. In this case the elements are considered as
        special elements: If it has no predecessors, then it is
        interpreted as an element smaller than any other, if it has no
        successors, then as larger than any other.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: z = P.null
            sage: oo = P.oo
            sage: z <= e
            True
            sage: e <= oo
            True
            sage: z <= oo
            True
            sage: oo <= z
            False
            sage: oo <= e
            False
            sage: e <= z
            False
            sage: z <= z
            True
            sage: oo <= oo
            True
            sage: e <= e
            True

        ::

            sage: z.le(e, reverse=True)
            False
            sage: e.le(oo, reverse=True)
            False
            sage: z.le(oo, reverse=True)
            False
            sage: oo.le(z, reverse=True)
            True
            sage: oo.le(e, reverse=True)
            True
            sage: e.le(z, reverse=True)
            True
            sage: z.le(z, reverse=True)
            True
            sage: oo.le(oo, reverse=True)
            True
            sage: e.le(e, reverse=True)
            True
        """
        if reverse:
            left, right = (right, left)

        if left.value is None:
            if not left.predecessors():
                # null on the left
                return True
            else:
                # oo on the left
                if right.value is None:
                    # null or oo on the right
                    return not right.successors()
                else:
                    # not null, not oo on the right
                    return False
        if right.value is None:
            if not right.successors():
                # oo on the right
                return True
            else:
                # null on the right
                if left.value is None:
                    # null or oo on the left
                    return not left.predecessors()
                else:
                    # not null, not oo on the right
                    return False
        return left.key <= right.key


    __le__ = le


    def eq(left, right):
        r"""
        Return if ``left`` is equal to ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        ``True`` or ``False``.

        This method compares the keys of the elements.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: z = P.null
            sage: oo = P.oo
            sage: z == z
            True
            sage: oo == oo
            True
            sage: e == e
            True
            sage: z == e
            False
            sage: e == oo
            False
            sage: oo == z
            False
        """
        if left.value is None and right.value is None:
            return left.is_null() == right.is_null()
        return left.key == right.key


    __eq__ = eq


    def _copy_all_linked_(self, memo, poset):
        r"""
        Helper function for :meth:`MutablePoset.copy`.

        INPUT:

        - ``memo`` -- a dictionary which assigns to the id of ``self``
          a copy of ``self``.

        - ``poset`` -- the poset to which the newly created element belongs.

        OUTPUT:

        A new element.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: Q = MP()
            sage: memo = {}
            sage: z = P.null._copy_all_linked_(memo, Q)
            sage: z.poset is Q
            True
            sage: oo = z.successors().pop()
            sage: oo == P.oo
            True
        """
        try:
            return memo[id(self)]
        except KeyError:
            pass

        new = self.__class__(poset, self.value)
        memo[id(self)] = new

        for reverse in (False, True):
            for e in self.successors(reverse):
                new.successors(reverse).add(e._copy_all_linked_(memo, poset))

        return new


    def _search_covers_(self, covers, element, reverse=False):
        r"""
        Helper function for :meth:`covers`.

        INPUT:

        - ``covers`` -- a set which finally contains all covers.

        - ``element`` -- the element for which to find the covering elements.

        - ``reverse`` -- (default: ``False``) if not set, then find
          the lower covers, otherwise find the upper covers.

        OUTPUT:

        ``True`` or ``False``.

        Note that ``False`` is returned if we do not have
        ``self <= element``.

       TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1, 1)))
            (1, 1, 1)
            sage: P.add(T((1, 3, 1)))
            (1, 3, 1)
            sage: P.add(T((2, 1, 2)))
            (2, 1, 2)
            sage: P.add(T((4, 4, 2)))
            (4, 4, 2)
            sage: P.add(T((1, 2, 2)))
            (1, 2, 2)
            sage: e = P.add(T((2, 2, 2))); e
            (2, 2, 2)
            sage: covers = set()
            sage: P.null._search_covers_(covers, e)
            True
            sage: sorted(covers, key=lambda c: repr(c.value))
            [(1, 2, 2), (2, 1, 2)]
        """
        if not self.le(element, reverse) or self == element:
            return False
        if not any([e._search_covers_(covers, element, reverse)
                    for e in self.successors(reverse)]):
            covers.add(self)
        return True


    def covers(self, element, reverse=False):
        r"""
        Return the covers of the given element (considering only
        elements which originate from ``self``).

        INPUT:

        - ``element`` -- the element for which to find the covering elements.

        - ``reverse`` -- (default: ``False``) if not set, then find
          the lower covers, otherwise find the upper covers.

        OUTPUT:

        A set of the covers.

        Suppose ``reverse`` is ``False``. This method returns all the
        lower covers of the given ``element``, i.e., elements in the
        poset, which are at most the given element and maximal with
        this property. Only elements which are (not necessarily
        direct) successors of ``self`` are considered.

        If ``reverse`` is ``True``, then the reverse direction is
        taken, i.e., in the text above replace lower covers by upper
        covers, maximal by minimal, and successors by predecessors.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            (1, 1)
            sage: P.add(T((1, 3)))
            (1, 3)
            sage: P.add(T((2, 1)))
            (2, 1)
            sage: P.add(T((4, 4)))
            (4, 4)
            sage: P.add(T((1, 2)))
            (1, 2)
            sage: e = P.add(T((2, 2))); e
            (2, 2)
            sage: sorted(P.null.covers(e),
            ....:        key=lambda c: repr(c.value))
            [(1, 2), (2, 1)]
            sage: sorted(P.oo.covers(e, reverse=True),
            ....:        key=lambda c: repr(c.value))
            [(4, 4)]
        """
        covers = set()
        self._search_covers_(covers, element, reverse)
        return covers


    def _iter_depth_first_visit_(self, marked, reverse=False, key=None):
        r"""
        Helper function for :meth:`iter_depth_first`.

        INPUT:

        - ``marked`` -- a set in which marked elements are stored.

        - ``reverse`` -- (default: ``False``) -- if set, reverses the order.

        - ``key`` -- (default: ``None``) a function used for sorting
          the successors. If this is ``None``, no sorting occurrs.

        OUTPUT:

        An iterator.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(42)
            42
            sage: P.add(5)
            5
            sage: marked = set()
            sage: list(P.oo._iter_depth_first_visit_(marked, True))
            [oo, 42, 5, null]
        """
        if self in marked:
            return
        marked.add(self)
        yield self
        S = self.successors(reverse)
        if key is not None:
            S = sorted(S, key=key)
        for element in S:
            for e in element._iter_depth_first_visit_(marked, reverse, key):
                yield e


    def iter_depth_first(self, reverse=False, key=None):
        r"""
        Iterates over all elements in depth first order.

        INPUT:

        - ``reverse`` -- (default: ``False``) -- if set, reverses the
          order, i.e., ``False`` starts at bottom (`\emptyset`),
          ``True`` starts at top (`\infty`).

        - ``key`` -- (default: ``None``) a function used for sorting
          the direct successors of an element (used in case of a
          tie). If this is ``None``, no sorting occurrs.

        OUTPUT:

        An iterator.

        ALGORITHM:

        See :wikipedia:`Depth-first_search`.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            (1, 1)
            sage: P.add(T((1, 3)))
            (1, 3)
            sage: P.add(T((2, 1)))
            (2, 1)
            sage: P.add(T((4, 4)))
            (4, 4)
            sage: P.add(T((1, 2)))
            (1, 2)
            sage: P.add(T((2, 2)))
            (2, 2)
            sage: list(P.null.iter_depth_first(reverse=False,
            ....:                                key=lambda c: repr(c)))
            [null, (1, 1), (1, 2), (1, 3), (4, 4), oo, (2, 2), (2, 1)]
            sage: list(P.oo.iter_depth_first(reverse=True,
            ....:                                key=lambda c: repr(c)))
            [oo, (4, 4), (1, 3), (1, 2), (1, 1), null, (2, 2), (2, 1)]
        """
        marked = set()
        return self._iter_depth_first_visit_(marked, reverse, key)


    def _iter_topological_visit_(self, marked, reverse=False, key=None):
        r"""
        Helper function for :meth:`iter_topological`.

        INPUT:

        - ``marked`` -- a set in which marked elements are stored.

        - ``reverse`` -- (default: ``False``) -- if set, reverses the order.

        - ``key`` -- (default: ``None``) a function used for sorting
          the successors. If this is ``None``, no sorting occurrs.

        OUTPUT:

        An iterator.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(42)
            42
            sage: P.add(5)
            5
            sage: marked = set()
            sage: list(P.null._iter_topological_visit_(marked, True))
            [oo, 42, 5, null]
        """
        if self in marked:
            return
        marked.add(self)
        S = self.predecessors(reverse)
        if key is not None:
            S = sorted(S, key=key)
        for element in S:
            for e in element._iter_topological_visit_(marked, reverse, key):
                yield e
        yield self


    def iter_topological(self, reverse=False, key=None):
        r"""
        Iterates over all elements in topological order.

        INPUT:

        - ``reverse`` -- (default: ``False``) -- if set, reverses the
          order, i.e., ``False`` gives smallest elements first,
          ``True`` gives largest first.

        - ``key`` -- (default: ``None``) a function used for sorting
          the direct successors of an element (used in case of a
          tie). If this is ``None``, no sorting occurrs.

        OUTPUT:

        An iterator.

        ALGORITHM:

        Here a simplified version of the algorithm found in [T1976]_
        and [CLRS2001]_ is used. See also
        :wikipedia:`Topological_sorting`.

        .. [T1976] Robert E. Tarjan, *Edge-disjoint spanning trees and
           depth-first search*, Acta Informatica 6 (2), 1976, 171-185,
           :doi:`10.1007/BF00268499`.

        .. [CLRS2001] Thomas H. Cormen, Charles E. Leiserson, Ronald
           L. Rivest and Clifford Stein, *Section 22.4: Topological
           sort*, Introduction to Algorithms (2nd ed.), MIT Press and
           McGraw-Hill, 2001, 549-552, ISBN 0-262-03293-7.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            (1, 1)
            sage: P.add(T((1, 3)))
            (1, 3)
            sage: P.add(T((2, 1)))
            (2, 1)
            sage: P.add(T((4, 4)))
            (4, 4)
            sage: P.add(T((1, 2)))
            (1, 2)
            sage: P.add(T((2, 2)))
            (2, 2)

        ::

            sage: for e in P.elements_topological(include_special=True,
            ....:                                 reverse=True):
            ....:     print e
            ....:     print list(e.iter_topological(reverse=True,
            ....:                                   key=lambda c: repr(c)))
            oo
            [oo]
            (4, 4)
            [oo, (4, 4)]
            (1, 3)
            [oo, (4, 4), (1, 3)]
            (2, 2)
            [oo, (4, 4), (2, 2)]
            (1, 2)
            [oo, (4, 4), (1, 3), (2, 2), (1, 2)]
            (2, 1)
            [oo, (4, 4), (2, 2), (2, 1)]
            (1, 1)
            [oo, (4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1)]
            null
            [oo, (4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1), null]

        ::

            sage: for e in P.elements_topological(include_special=True,
            ....:                                 reverse=True):
            ....:     print e
            ....:     print list(e.iter_topological(reverse=False,
            ....:                                   key=lambda c: repr(c)))
            oo
            [null, (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (4, 4), oo]
            (4, 4)
            [null, (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (4, 4)]
            (1, 3)
            [null, (1, 1), (1, 2), (1, 3)]
            (2, 2)
            [null, (1, 1), (1, 2), (2, 1), (2, 2)]
            (1, 2)
            [null, (1, 1), (1, 2)]
            (2, 1)
            [null, (1, 1), (2, 1)]
            (1, 1)
            [null, (1, 1)]
            null
            [null]
        """
        marked = set()
        return self._iter_topological_visit_(marked, reverse, key)


# *****************************************************************************


def sorted_set_by_tuple(S, T):
    r"""
    Return an iterator over ``S`` respecting the order given by ``T``.

    INPUT:

    - ``S`` -- a set (or something which supports containment test).

    - ``T`` -- a tuple (or other iterable).

    OUTPUT:

    An iterator.

    In the iterator all elements of ``T``, which are also in ``S``
    appear. The order given by ``T`` is kept.

    EXAMPLES::

        sage: from sage.data_structures.mutable_poset import sorted_set_by_tuple
        sage: tuple(sorted_set_by_tuple({3, 4, 6}, (5, 4, 1, 2, 3, 6)))
        (4, 3, 6)
    """
    return iter(ell for ell in T if ell in S)


# *****************************************************************************


class MutablePoset(sage.structure.sage_object.SageObject):
    r"""
    A mutable poset (partially ordered set) as data structure.
    """
    def __init__(self, data=None, key=None):
        r"""
        See :class:`MutablePoset` for details.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: MP()
            poset()
        """

        if data is not None:
            raise NotImplementedError

        self._null_ = MutablePosetElement(self, None)
        self._oo_ = MutablePosetElement(self, None)
        self._null_.successors().add(self._oo_)
        self._oo_.predecessors().add(self._null_)
        self._elements_ = {}

        if key is None:
            self._key_ = lambda k: k
        else:
            self._key_ = key


    @property
    def null(self):
        r"""
        The element `\emptyset` which is smaller than any other element.

        EXAMPLES:

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: z = P.null; z
            null
            sage: z.is_null()
            True
        """
        return self._null_


    @property
    def oo(self):
        r"""
        The element `\infty` which is larger than any other element.

        EXAMPLES:

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: oo = P.oo; oo
            oo
            sage: oo.is_oo()
            True
        """
        return self._oo_


    def element(self, key):
        r"""
        Return the element corresponding to ``key``.

        INPUT:

        ``key`` -- the key of an object.

        OUTPUT:

        An instance of :class:`MutablePosetElement`.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P.add(42)
            sage: e = P.element(42); e
            42
            sage: type(e)
            <class 'sage.data_structures.mutable_poset.MutablePosetElement'>
        """
        return self._elements_[key]


    def get_key(self, value):
        r"""
        Return the key corresponding to ``value``.

        INPUT:

        - ``value`` -- an object.

        OUTPUT:

        An object (the key of ``value``).

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: P = MP()
            sage: P.get_key(None) is None
            True
            sage: P.get_key((1, 2))
            (1, 2)
            sage: Q = MP(key=lambda k: k[0])
            sage: Q.get_key((1, 2))
            1
        """
        if value is None:
            return None
        return self._key_(value)


    def copy(self):
        r"""
        Creates a shallow copy.

        INPUT:

        Nothing.

        OUTPUT:

        A poset with the same content as ``self``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            (1, 1)
            sage: P.add(T((1, 3)))
            (1, 3)
            sage: P.add(T((2, 1)))
            (2, 1)
            sage: P.add(T((4, 4)))
            (4, 4)
            sage: P.add(T((1, 2)))
            (1, 2)
            sage: Q = P.copy()
            sage: P.repr_full() == Q.repr_full()
            True
        """
        new = self.__class__()
        memo = {}
        new._null_ = self._null_._copy_all_linked_(memo, new)
        new._oo_ = memo[id(self._oo_)]
        new._elements_ = dict((f.key, f) for f in
                              iter(memo[id(e)]
                                   for e in self._elements_.itervalues()))
        return new


    __copy__ = copy


    def elements(self, include_special=False, reverse=False):
        r"""
        Return an iterator over all elements.

        INPUT:

        - ``include_special`` -- (default: ``False``) if set, then
          including a smallest element (`\emptyset`) and a largest element
          (`\infty`).

        - ``reverse`` -- (default: ``False``) if set, the order is
          reversed. This only affects the elements `\emptyset` and `\infty`.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: tuple(P.elements())
            ()
            sage: tuple(P.elements(include_special=True))
            (null, oo)
            sage: tuple(P.elements(include_special=True, reverse=True))
            (oo, null)
        """
        if include_special:
            yield self.null if not reverse else self.oo
        for e in self._elements_.itervalues():
            yield e
        if include_special:
            yield self.oo if not reverse else self.null


    def elements_topological(self, include_special=False,
                             reverse=False, key=None):
        r"""
        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            (1, 1)
            sage: P.add(T((1, 3)))
            (1, 3)
            sage: P.add(T((2, 1)))
            (2, 1)
            sage: P.add(T((4, 4)))
            (4, 4)
            sage: P.add(T((1, 2)))
            (1, 2)
            sage: P.add(T((2, 2)))
            (2, 2)
            sage: list(P.elements_topological())
            [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (4, 4)]
            sage: list(P.elements_topological(reverse=True))
            [(4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1)]
            sage: list(P.elements_topological(include_special=True))
            [null, (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (4, 4), oo]
            sage: list(P.elements_topological(
            ....:     include_special=True, reverse=True))
            [oo, (4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1), null]
        """
        if key is None:
            key = lambda c: repr(c)
        element = self.oo if not reverse else self.null
        return iter(e for e in element.iter_topological(reverse, key)
                    if include_special or not e.is_special())


    def repr(self, include_special=False, reverse=False):
        r"""
        Return a representation of the poset.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: print MP().repr()
            poset()
        """
        s = 'poset('
        s += ', '.join(repr(element) for element in
                       self.elements_topological(include_special, reverse))
        s += ')'
        return s


    def repr_full(self, reverse=False):
        r"""
        Return a representation with ordering details of the poset.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: print MP().repr_full(reverse=True)
            poset()
            +-- oo
            |   +-- no successors
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   oo
            |   +-- no predecessors
        """
        sortedelements = tuple(
            self.elements_topological(include_special=True, reverse=reverse))
        strings = [self.repr(include_special=False, reverse=reverse)]
        for element in sortedelements:
            strings.append('+-- ' + repr(element))
            for rev in (not reverse, reverse):
                what = 'successors' if not rev else 'predecessors'
                if element.successors(rev):
                    s = '|   +-- ' + what + ':   '
                    s += ', '.join(repr(e) for e in
                                   sorted_set_by_tuple(element.successors(rev),
                                                       sortedelements))
                else:
                    s = '|   +-- no ' + what
                strings.append(s)
        return '\n'.join(strings)


    __repr__ = repr


    def __contains__(self, key):
        r"""
        Tests if ``key`` is encapsulated by one of the poset's elements.

        INPUT:

        - ``key`` -- an object.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            (1, 1)
            sage: T((1, 1)) in P
            True
            sage: T((1, 2)) in P
            False
        """
        return key in self._elements_


    def add(self, value):
        r"""
        Add the given object as element to the poset.

        INPUT:

        - ``value`` -- an object (hashable and supporting comparison
          with the operator ``<=``.

        OUTPUT:

        The created (or already existing) element.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            (1, 1)
            sage: P.add(T((1, 3)))
            (1, 3)
            sage: P.add(T((2, 1)))
            (2, 1)
            sage: P.add(T((4, 4)))
            (4, 4)
            sage: P.add(T((1, 2)))
            (1, 2)
            sage: print P.repr_full(reverse=True)
            poset((4, 4), (1, 3), (1, 2), (2, 1), (1, 1))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4)
            +-- (4, 4)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3), (2, 1)
            +-- (1, 3)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2)
            +-- (1, 2)
            |   +-- successors:   (1, 3)
            |   +-- predecessors: (1, 1)
            +-- (2, 1)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 1)
            +-- (1, 1)
            |   +-- successors:   (1, 2), (2, 1)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 1)
            |   +-- no predecessors
            sage: P.add(T((2, 2)))
            (2, 2)
            sage: print P.repr_full(reverse=True)
            poset((4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4)
            +-- (4, 4)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3), (2, 2)
            +-- (1, 3)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2)
            +-- (2, 2)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2), (2, 1)
            +-- (1, 2)
            |   +-- successors:   (1, 3), (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (2, 1)
            |   +-- successors:   (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (1, 1)
            |   +-- successors:   (1, 2), (2, 1)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 1)
            |   +-- no predecessors

        When adding an element which is already in the poset, the
        existing one is returned::

            sage: e = T((2, 2))
            sage: f = P.add(e).value; f
            (2, 2)
            sage: e == f, e is f
            (True, False)

        TESTS::

            sage: R = MP(key=lambda k: T(k[2:3]))
            sage: R.add(T((1, 1, 42)))
            (1, 1, 42)
            sage: R.add(T((1, 3, 42)))
            (1, 1, 42)
            sage: R.add(T((2, 1, 7)))
            (2, 1, 7)
            sage: R.add(T((4, 4, 42)))
            (1, 1, 42)
            sage: R.add(T((1, 2, 7)))
            (2, 1, 7)
            sage: R.add(T((2, 2, 7)))
            (2, 1, 7)
            sage: print R.repr_full(reverse=True)
            poset((1, 1, 42), (2, 1, 7))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (1, 1, 42)
            +-- (1, 1, 42)
            |   +-- successors:   oo
            |   +-- predecessors: (2, 1, 7)
            +-- (2, 1, 7)
            |   +-- successors:   (1, 1, 42)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (2, 1, 7)
            |   +-- no predecessors
        """
        if value is None:
            raise ValueError('None is not allowed as value.')
        key = self.get_key(value)

        try:
            return self._elements_[key]
        except KeyError:
            pass

        new = MutablePosetElement(self, value)
        smaller = self.null.covers(new, reverse=False)
        larger = self.oo.covers(new, reverse=True)

        for reverse in (False, True):
            sm = smaller if not reverse else larger
            la = larger if not reverse else smaller
            for element in sm:
                for e in element.successors(reverse).intersection(la):
                    e.predecessors(reverse).remove(element)
                    element.successors(reverse).remove(e)
                new.predecessors(reverse).add(element)
                element.successors(reverse).add(new)

        self._elements_[key] = new
        return new


    def remove(self, key):
        r"""
        Remove the given object from the poset.

        INPUT:

        - ``value`` -- an object.

        OUTPUT:

        Nothing.

        If the element is not a member, raise a ``KeyError``.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add(T((1, 1)))
            (1, 1)
            sage: P.add(T((1, 3)))
            (1, 3)
            sage: P.add(T((2, 1)))
            (2, 1)
            sage: P.add(T((4, 4)))
            (4, 4)
            sage: P.add(T((1, 2)))
            (1, 2)
            sage: P.add(T((2, 2)))
            (2, 2)
            sage: print P.repr_full(reverse=True)
            poset((4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4)
            +-- (4, 4)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3), (2, 2)
            +-- (1, 3)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2)
            +-- (2, 2)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2), (2, 1)
            +-- (1, 2)
            |   +-- successors:   (1, 3), (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (2, 1)
            |   +-- successors:   (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (1, 1)
            |   +-- successors:   (1, 2), (2, 1)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 1)
            |   +-- no predecessors
            sage: P.remove(T((1, 2)))
            sage: print P.repr_full(reverse=True)
            poset((4, 4), (1, 3), (2, 2), (2, 1), (1, 1))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4)
            +-- (4, 4)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3), (2, 2)
            +-- (1, 3)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 1)
            +-- (2, 2)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (2, 1)
            +-- (2, 1)
            |   +-- successors:   (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (1, 1)
            |   +-- successors:   (1, 3), (2, 1)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 1)
            |   +-- no predecessors

        TESTS::

            sage: Q = MP(key=lambda k: T(k[0:2]))
            sage: Q.add(T((1, 1, 42)))
            (1, 1, 42)
            sage: Q.add(T((1, 3, 42)))
            (1, 3, 42)
            sage: Q.add(T((2, 1, 7)))
            (2, 1, 7)
            sage: Q.add(T((4, 4, 42)))
            (4, 4, 42)
            sage: Q.add(T((1, 2, 7)))
            (1, 2, 7)
            sage: Q.add(T((2, 2, 7)))
            (2, 2, 7)
            sage: print Q.repr_full(reverse=True)
            poset((4, 4, 42), (1, 3, 42), (2, 2, 7),
                  (1, 2, 7), (2, 1, 7), (1, 1, 42))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4, 42)
            +-- (4, 4, 42)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3, 42), (2, 2, 7)
            +-- (1, 3, 42)
            |   +-- successors:   (4, 4, 42)
            |   +-- predecessors: (1, 2, 7)
            +-- (2, 2, 7)
            |   +-- successors:   (4, 4, 42)
            |   +-- predecessors: (1, 2, 7), (2, 1, 7)
            +-- (1, 2, 7)
            |   +-- successors:   (1, 3, 42), (2, 2, 7)
            |   +-- predecessors: (1, 1, 42)
            +-- (2, 1, 7)
            |   +-- successors:   (2, 2, 7)
            |   +-- predecessors: (1, 1, 42)
            +-- (1, 1, 42)
            |   +-- successors:   (1, 2, 7), (2, 1, 7)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 1, 42)
            |   +-- no predecessors
            sage: Q.remove((1,1))
            sage: print Q.repr_full(reverse=True)
            poset((4, 4, 42), (1, 3, 42), (2, 2, 7), (1, 2, 7), (2, 1, 7))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4, 42)
            +-- (4, 4, 42)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3, 42), (2, 2, 7)
            +-- (1, 3, 42)
            |   +-- successors:   (4, 4, 42)
            |   +-- predecessors: (1, 2, 7)
            +-- (2, 2, 7)
            |   +-- successors:   (4, 4, 42)
            |   +-- predecessors: (1, 2, 7), (2, 1, 7)
            +-- (1, 2, 7)
            |   +-- successors:   (1, 3, 42), (2, 2, 7)
            |   +-- predecessors: null
            +-- (2, 1, 7)
            |   +-- successors:   (2, 2, 7)
            |   +-- predecessors: null
            +-- null
            |   +-- successors:   (1, 2, 7), (2, 1, 7)
            |   +-- no predecessors
        """
        if key is None:
            raise ValueError('None is not allowed as key.')

        try:
            element = self._elements_[key]
        except KeyError:
            raise KeyError('Key %s is not contained in this poset.' % (key,))

        for reverse in (False, True):
            for p in element.predecessors(reverse):
                S = p.successors(reverse)
                S.remove(element)
                D = set(s for s in p.iter_depth_first(reverse)
                        if s in element.successors(reverse))
                S.update(element.successors(reverse))
                S.difference_update(D)
        del self._elements_[key]


# *****************************************************************************


class MutableTosetElement(MutablePosetElement):
    pass


# *****************************************************************************


class MutableToset(MutablePoset):
    r"""
    A mutable toset (totally ordered set) as data structure.
    """
    pass


# *****************************************************************************
