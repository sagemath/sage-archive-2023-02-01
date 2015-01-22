r"""
Asymptotic Ring
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

        INPUT:

        Nothing.

        OUTPUT:

        A poset.

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

        INPUT:

        Nothing.

        OUTPUT:

        An object.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: e.value
            (1, 2)
        """
        return self._value_


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

        Return if the element is either the zero-element, i.e., the
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
            sage: P._zero_.is_special()
            True
            sage: P._oo_.is_special()
            True
        """
        return self.value is None


    def is_zero(self, reverse=False):
        r"""
        Return if the element is the zero-element, i.e., the element
        smaller than any possible other element.

        INPUT:

        - ``reverse`` -- (default: ``False``) if set, then returns
          if the element is the largest possible.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: P._zero_.is_zero()
            True
            sage: P._oo_.is_zero()
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
            sage: P._zero_.is_oo()
            False
            sage: P._oo_.is_oo()
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
        ``None``. In this case either ``'zero'`` or ``'oo'`` is
        returned depending in the nonexistence of predecessors and
        sucessors respectively.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: repr(MutablePosetElement(P, (1, 2)))  # indirect doctest
            '(1, 2)'
            sage: repr(P._zero_)  # indirect doctest
            'zero'
            sage: repr(P._oo_)  # indirect doctest
            'oo'
        """
        if self.value is None:
            if not self.predecessors():
                return 'zero'
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

        This returns the hash value of the element.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: hash(MutablePosetElement(P, (1, 2))) == hash((1, 2))
            True
        """
        return hash(self.value)


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

        This methods usually returns if the values of the given
        elements are less or equal. The only exception is if this
        value is ``None``. In this case the elements are considered as
        special elements: If it has no predecessors, then it is
        interpreted as an element smaller than any other, if it has no
        successors, then as larger than any other.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: z = P._zero_
            sage: oo = P._oo_
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
                # zero on the left
                return True
            else:
                # oo on the left
                if right.value is None:
                    # zero or oo on the right
                    return not right.successors()
                else:
                    # not zero, not oo on the right
                    return False
        if right.value is None:
            if not right.successors():
                # oo on the right
                return True
            else:
                # zero on the right
                if left.value is None:
                    # zero or oo on the left
                    return not left.predecessors()
                else:
                    # not zero, not oo on the right
                    return False
        return left.value <= right.value


    __le__ = le


    def eq(left, right):
        r"""
        Return if ``left`` is equal to ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        ``True`` or ``False``.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: from sage.data_structures.mutable_poset import MutablePosetElement
            sage: e = MutablePosetElement(P, (1, 2))
            sage: z = P._zero_
            sage: oo = P._oo_
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
            return left.is_zero() == right.is_zero()
        return left.value == right.value


    __eq__ = eq


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
            sage: P.add_element(T((1, 1, 1)))
            (1, 1, 1)
            sage: P.add_element(T((1, 3, 1)))
            (1, 3, 1)
            sage: P.add_element(T((2, 1, 2)))
            (2, 1, 2)
            sage: P.add_element(T((4, 4, 2)))
            (4, 4, 2)
            sage: P.add_element(T((1, 2, 2)))
            (1, 2, 2)
            sage: e = P.add_element(T((2, 2, 2))); e
            (2, 2, 2)
            sage: covers = set()
            sage: P._zero_._search_covers_(covers, e)
            True
            sage: covers
            {(2, 1, 2), (1, 2, 2)}
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
            sage: P.add_element(T((1, 1)))
            (1, 1)
            sage: P.add_element(T((1, 3)))
            (1, 3)
            sage: P.add_element(T((2, 1)))
            (2, 1)
            sage: P.add_element(T((4, 4)))
            (4, 4)
            sage: P.add_element(T((1, 2)))
            (1, 2)
            sage: e = P.add_element(T((2, 2))); e
            (2, 2)
            sage: P._zero_.covers(e)
            {(2, 1), (1, 2)}
            sage: P._oo_.covers(e, reverse=True)
            {(4, 4)}
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
            sage: P.add_element(42)
            42
            sage: P.add_element(5)
            5
            sage: marked = set()
            sage: list(P._oo_._iter_depth_first_visit_(marked, True))
            [oo, 42, 5, zero]
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
          order, i.e., ``False`` starts at bottom (`0`),
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
            sage: P.add_element(T((1, 1)))
            (1, 1)
            sage: P.add_element(T((1, 3)))
            (1, 3)
            sage: P.add_element(T((2, 1)))
            (2, 1)
            sage: P.add_element(T((4, 4)))
            (4, 4)
            sage: P.add_element(T((1, 2)))
            (1, 2)
            sage: P.add_element(T((2, 2)))
            (2, 2)
            sage: list(P._zero_.iter_depth_first(reverse=False,
            ....:                                key=lambda c: repr(c)))
            [zero, (1, 1), (1, 2), (1, 3), (4, 4), oo, (2, 2), (2, 1)]
            sage: list(P._oo_.iter_depth_first(reverse=True,
            ....:                                key=lambda c: repr(c)))
            [oo, (4, 4), (1, 3), (1, 2), (1, 1), zero, (2, 2), (2, 1)]
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
            sage: P.add_element(42)
            42
            sage: P.add_element(5)
            5
            sage: marked = set()
            sage: list(P._zero_._iter_topological_visit_(marked, True))
            [oo, 42, 5, zero]
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
        and [CLRS]_ is used. See also
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
            sage: P.add_element(T((1, 1)))
            (1, 1)
            sage: P.add_element(T((1, 3)))
            (1, 3)
            sage: P.add_element(T((2, 1)))
            (2, 1)
            sage: P.add_element(T((4, 4)))
            (4, 4)
            sage: P.add_element(T((1, 2)))
            (1, 2)
            sage: P.add_element(T((2, 2)))
            (2, 2)

        ::

            sage: for e in P.elements_topological(reverse=True):
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
            zero
            [oo, (4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1), zero]

        ::

            sage: for e in P.elements_topological(reverse=True):
            ....:     print e
            ....:     print list(e.iter_topological(reverse=False,
            ....:                                   key=lambda c: repr(c)))
            oo
            [zero, (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (4, 4), oo]
            (4, 4)
            [zero, (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (4, 4)]
            (1, 3)
            [zero, (1, 1), (1, 2), (1, 3)]
            (2, 2)
            [zero, (1, 1), (1, 2), (2, 1), (2, 2)]
            (1, 2)
            [zero, (1, 1), (1, 2)]
            (2, 1)
            [zero, (1, 1), (2, 1)]
            (1, 1)
            [zero, (1, 1)]
            zero
            [zero]
        """
        marked = set()
        return self._iter_topological_visit_(marked, reverse, key)


# *****************************************************************************


def _sort_set_by_tuple_iter_(S, T):
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

        sage: from sage.data_structures.mutable_poset import _sort_set_by_tuple_iter_
        sage: tuple(_sort_set_by_tuple_iter_({3, 4, 6}, (5, 4, 1, 2, 3, 6)))
        (4, 3, 6)
    """
    return iter(ell for ell in T if ell in S)


# *****************************************************************************


class MutablePoset(sage.structure.sage_object.SageObject):
    r"""
    A mutable poset.
    """
    def __init__(self, data=None):
        r"""
        See :class:`MutablePoset` for details.

        TESTS::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: MP()
            poset()
        """

        if data is not None:
            raise NotImplementedError

        self._zero_ = MutablePosetElement(self, None)
        self._oo_ = MutablePosetElement(self, None)
        self._zero_.successors().add(self._oo_)
        self._oo_.predecessors().add(self._zero_)
        self._elements_ = {}



    def elements(self, include_special=False, reverse=False):
        r"""
        Return an iterator over all elements.

        INPUT:

        - ``include_special`` -- (default: ``False``) if set, then
          including a smallest element (`0`) and a largest element
          (`\infty`).

        - ``reverse`` -- (default: ``False``) if set, the order is
          reversed. This only affects the elements `0` and `\infty`.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: P = MP()
            sage: tuple(P.elements())
            ()
            sage: tuple(P.elements(include_special=True))
            (zero, oo)
            sage: tuple(P.elements(include_special=True, reverse=True))
            (oo, zero)
        """
        if include_special:
            yield self._zero_ if not reverse else self._oo_
        for e in self._elements_.itervalues():
            yield e
        if include_special:
            yield self._oo_ if not reverse else self._zero_


    def elements_topological(self, reverse=False, key=None):
        r"""
        EXAMPLES::

            sage: from sage.data_structures.mutable_poset import MutablePoset as MP
            sage: class T(tuple):
            ....:     def __le__(left, right):
            ....:         return all(l <= r for l, r in zip(left, right))
            sage: P = MP()
            sage: P.add_element(T((1, 1)))
            (1, 1)
            sage: P.add_element(T((1, 3)))
            (1, 3)
            sage: P.add_element(T((2, 1)))
            (2, 1)
            sage: P.add_element(T((4, 4)))
            (4, 4)
            sage: P.add_element(T((1, 2)))
            (1, 2)
            sage: P.add_element(T((2, 2)))
            (2, 2)
            sage: list(P.elements_topological())
            [zero, (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (4, 4), oo]
            sage: list(P.elements_topological(reverse=True))
            [oo, (4, 4), (1, 3), (2, 2), (1, 2), (2, 1), (1, 1), zero]
        """
        if key is None:
            key = lambda c: repr(c)
        element = self._oo_ if not reverse else self._zero_
        return element.iter_topological(reverse, key)


    def repr(self):
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
        s += ', '.join(repr(element) for element in self.elements())
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
            |   +-- predecessors: zero
            +-- zero
            |   +-- successors:   oo
            |   +-- no predecessors
        """
        sortedelements = tuple(
            self.elements(include_special=True, reverse=reverse))
        strings = [self.repr()]
        for element in sortedelements:
            s = '+-- ' + repr(element) + '\n'
            if element.successors():
                s += '|   +-- successors:   '
                s += ', '.join(repr(e) for e in
                               _sort_set_by_tuple_iter_(element.successors(),
                                                        sortedelements))
            else:
                s += '|   +-- no successors'
            s += '\n'
            if element.predecessors():
                s += '|   +-- predecessors: '
                s += ', '.join(repr(e) for e in
                               _sort_set_by_tuple_iter_(element.predecessors(),
                                                        sortedelements))
            else:
                s += '|   +-- no predecessors'
            strings.append(s)
        return '\n'.join(strings)


    __repr__ = repr



    def add_element(self, value):
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
            sage: P.add_element(T((1, 1)))
            (1, 1)
            sage: P.add_element(T((1, 3)))
            (1, 3)
            sage: P.add_element(T((2, 1)))
            (2, 1)
            sage: P.add_element(T((4, 4)))
            (4, 4)
            sage: P.add_element(T((1, 2)))
            (1, 2)
            sage: print P.repr_full(reverse=True)
            poset((1, 2), (1, 3), (1, 1), (2, 1), (4, 4))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4)
            +-- (1, 2)
            |   +-- successors:   (1, 3)
            |   +-- predecessors: (1, 1)
            +-- (1, 3)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2)
            +-- (1, 1)
            |   +-- successors:   (1, 2), (2, 1)
            |   +-- predecessors: zero
            +-- (2, 1)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 1)
            +-- (4, 4)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3), (2, 1)
            +-- zero
            |   +-- successors:   (1, 1)
            |   +-- no predecessors
            sage: P.add_element(T((2, 2)))
            (2, 2)
            sage: print P.repr_full(reverse=True)
            poset((1, 2), (1, 3), (4, 4), (2, 1), (2, 2), (1, 1))
            +-- oo
            |   +-- no successors
            |   +-- predecessors: (4, 4)
            +-- (1, 2)
            |   +-- successors:   (1, 3), (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (1, 3)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2)
            +-- (4, 4)
            |   +-- successors:   oo
            |   +-- predecessors: (1, 3), (2, 2)
            +-- (2, 1)
            |   +-- successors:   (2, 2)
            |   +-- predecessors: (1, 1)
            +-- (2, 2)
            |   +-- successors:   (4, 4)
            |   +-- predecessors: (1, 2), (2, 1)
            +-- (1, 1)
            |   +-- successors:   (1, 2), (2, 1)
            |   +-- predecessors: zero
            +-- zero
            |   +-- successors:   (1, 1)
            |   +-- no predecessors

        When adding an element which is already in the poset, the
        existing one is returned::

            sage: e = T((2, 2))
            sage: f = P.add_element(e).value; f
            (2, 2)
            sage: e == f, e is f
            (True, False)
        """
        try:
            return self._elements_[value]
        except KeyError:
            pass
        if value is None:
            raise ValueError('None is not allowed as value.')

        new = MutablePosetElement(self, value)
        smaller = self._zero_.covers(new, reverse=False)
        larger = self._oo_.covers(new, reverse=True)

        for reverse in (False, True):
            sm = smaller if not reverse else larger
            la = larger if not reverse else smaller
            for element in sm:
                for e in element.successors(reverse).intersection(la):
                    e.predecessors(reverse).remove(element)
                    element.successors(reverse).remove(e)
                new.predecessors(reverse).add(element)
                element.successors(reverse).add(new)

        self._elements_[value] = new
        return new


# *****************************************************************************
