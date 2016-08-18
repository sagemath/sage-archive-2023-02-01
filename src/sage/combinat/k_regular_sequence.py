r"""
`k`-regular Sequences

An introduction and formal definition of `k`-regular sequences can be
found, for example, on the :wikipedia:`k-regular_sequence` or in
[AS2003]_.


Examples
========

Binary sum of digits
--------------------

::

    sage: Seq2 = kRegularSequenceSpace(2, ZZ)
    sage: S = Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
    ....:          left=vector([0, 1]), right=vector([1, 0]))
    sage: S
    2-regular sequence 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, ...
    sage: all(S[n] == sum(n.digits(2)) for n in srange(10))
    True

Number of odd entries in Pascal's triangle
------------------------------------------

::

    sage: @cached_function
    ....: def u(n):
    ....:     if n <= 1:
    ....:         return n
    ....:     return 2*u(floor(n/2)) + u(ceil(n/2))
    sage: tuple(u(n) for n in srange(10))
    (0, 1, 3, 5, 9, 11, 15, 19, 27, 29)

    sage: U = Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
    ....:          left=vector([0, 1]), right=vector([1, 0]), transpose=True)
    sage: all(U[n] == u(n) for n in srange(30))
    True


Various
=======

.. SEEALSO::

    :doc:`sage/rings/cfinite_sequence`,
    :doc:`sage/combinat/binary_recurrence_sequences`.

REFERENCES:

.. [AS2003] Jean-Paul Allouche, Jeffrey Shallit,
   *Automatic Sequences: Theory, Applications, Generalizations*,
   Cambridge University Press, 2003.

AUTHORS:

- Daniel Krenn (2016)

ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the
  Austrian Science Fund (FWF): P 24644-N26.


Classes and Methods
===================
"""
#*****************************************************************************
#       Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.recognizable_series import RecognizableSeries
from sage.combinat.recognizable_series import RecognizableSeriesSpace
from sage.misc.cachefunc import cached_method
from sage.structure.element import Element

class kRegularSequence(RecognizableSeries):

    def __init__(self, parent, mu, left=None, right=None):
        r"""
        A `k`-regular sequence.

        INPUT:

        - ``parent`` -- an instance of :class:`kRegularSequenceSpace`.

        - ``mu`` -- a tuple or other iterable of square matrices,
          all of which have the same dimension.

        - ``left`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the left to the matrix product. If ``None``, then this
          multiplication is skipped.

        - ``right`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the left to the matrix product. If ``None``, then this
          multiplication is skipped.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:      vector([0, 1]), vector([1, 0]),
            ....:      transpose=True)
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...

        .. SEEALSO::

            :doc:`k-regular sequence <k_regular_sequence>`,
            :class:`kRegularSequenceSpace`.
        """
        super(kRegularSequence, self).__init__(
            parent=parent, mu=mu, left=left, right=right)


    def _repr_(self):
        r"""
        Return a representation string of this `k`-regular sequence.

        OUTPUT:

        A string

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: s = Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:           vector([0, 1]), vector([1, 0]), transpose=True)
            sage: repr(s)  # indirect doctest
            '2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...'
        """
        from sage.misc.lazy_list import lazy_list_formatter
        return lazy_list_formatter(
            self,
            name='{}-regular sequence'.format(self.parent().k),
            opening_delimiter='', closing_delimiter='',
            preview=10)



    def __getitem__(self, n):
        r"""
        Return the `n`th entry of this sequence.

        INPUT:

        - ``n`` -- a nonnegative integer.

        OUTPUT:

        An element of the universe of the sequence.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: S = Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:          left=vector([0, 1]), right=vector([1, 0]))
            sage: S[7]
            3

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: W = Seq2.indices()
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Seq2((M0, M1))
            sage: S._mu_of_word_(W(0.digits(2))) == M0
            True
            sage: S._mu_of_word_(W(1.digits(2))) == M1
            True
            sage: S._mu_of_word_(W(3.digits(2))) == M1^2
            True

        ::

            sage: S._product_of_mu_(-1)
            Traceback (most recent call last):
            ...
            ValueError: m=-1 is not a nonnegative integer.
        """
        from sage.rings.integer_ring import ZZ
        n = ZZ(n)
        W = self.parent().indices()
        w = W(n.digits(self.parent().k))
        return super(kRegularSequence, self).__getitem__(w)


    def __iter__(self):
        r"""
        Return an iterator.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: S = Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:          left=vector([0, 1]), right=vector([1, 0]))
            sage: from itertools import islice
            sage: tuple(islice(S, 10))
             (0, 1, 1, 2, 1, 2, 2, 3, 1, 2)

        TESTS::

            sage: it = iter(S)
            sage: iter(it) is it
            True
            sage: iter(S) is not it
            True
        """
        from itertools import count
        return iter(self[n] for n in count())


class kRegularSequenceSpace(RecognizableSeriesSpace):

    Element = kRegularSequence

    @classmethod
    def __normalize__(cls, k, coefficients=None, **kwds):
        r"""
        Normalizes the input in order to ensure a unique
        representation.

        For more information see :class:`kRegularSequenceSpace`.

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2.category()
            Category of sets

        ::

            sage: Seq2 is kRegularSequenceSpace(2)
            True
        """
        from sage.arith.srange import srange
        from sage.rings.integer_ring import ZZ
        if coefficients is None:
            coefficients = ZZ
        nargs = super(kRegularSequenceSpace, cls).__normalize__(
            coefficients, alphabet=srange(k), **kwds)
        return (k,) + nargs


    def __init__(self, k, coefficients, indices, category):
        r"""
        The space of `k`-regular Sequences over the given ``universe``.

        INPUT:

        - ``k`` -- an integer at least `2` specifying the base.

        - ``coefficients`` -- a (semi-)ring.

        - ``category`` -- (default: ``None``) the category of the
          sequence space. If ``None``, then the category of
          :class:`~sage.categories.sets_cat.Sets` is used.

        EXAMPLES::

            sage: kRegularSequenceSpace(2, ZZ)
            Space of 2-regular sequences over Integer Ring
            sage: kRegularSequenceSpace(3, ZZ)
            Space of 3-regular sequences over Integer Ring

        .. SEEALSO::

            :doc:`k-regular sequence <k_regular_sequence>`,
            :class:`kRegularSequence`.
        """
        self.k = k
        super(kRegularSequenceSpace, self).__init__(coefficients, indices, category)


    def _repr_(self):
        r"""
        Return a representation string of this `k`-regular sequence space.

        OUTPUT:

        A string

        TESTS::

            sage: repr(kRegularSequenceSpace(2, ZZ))  # indirect doctest
            'Space of 2-regular sequences over Integer Ring'
        """
        return 'Space of {}-regular sequences over {}'.format(self.k, self.base())

