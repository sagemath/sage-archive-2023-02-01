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
    ....:          initial=vector([0, 1]), selection=vector([1, 0]))
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
    ....:          initial=vector([0, 1]), selection=vector([1, 0]), transpose=True)
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

from sage.misc.cachefunc import cached_method
from sage.structure.element import Element


class kRegularSequence(Element):

    def __init__(self, parent, matrices, initial=None, selection=None,
                 output_function=None, transpose=False):
        r"""
        A `k`-regular sequence.

        INPUT:

        - ``parent`` -- an instance of :class:`kRegularSequenceSpace`.

        - ``matrices`` -- a tuple or other iterable of square matrices,
          all of which have the same dimension.

        - ``initial`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the left to the matrix product. If ``None``, then this
          multiplication is skipped.

        - ``selection`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the left to the matrix product. If ``None``, then this
          multiplication is skipped.

        - ``output_function`` -- (default: ``None``) a function, which is
          applied after evaluating the sequence. This may be used to
          extract the value of a 1x1 matrix.

        - ``transpose`` -- (default: ``False``) a boolean. If set, then
          each of the ``matrices``. Additionally the vectors ``initial``
          and ``selection`` are switched and (if possible)
          transposed as well.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:      vector([0, 1]), vector([1, 0]),
            ....:      transpose=True)
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...

        Using an output function::

            sage: Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:      Matrix([[0, 1]]), Matrix([[1], [0]]),
            ....:      lambda o: o[0, 0], transpose=True)
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...
        """
        super(kRegularSequence, self).__init__(parent=parent)

        def tr(M):
            try:
                return M.transpose() if transpose else M
            except AttributeError:
                return M

        self.matrices = tuple(tr(M) for M in matrices)
        self.k = len(self.matrices)
        self.d = self.matrices[0].nrows()
        if not all(M.dimensions() == (self.d, self.d) for M in self.matrices):
            raise ValueError  # TODO

        if not transpose:
            self.initial = initial
            self.selection = selection
        else:
            self.initial = tr(selection)
            self.selection = tr(initial)

        if output_function is None:
            self.output_function = lambda o: o
        else:
            self.output_function = output_function


    def _repr_(self):
        r"""
        Return a representation string of this `k`-regular sequence

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
            [self[n] for n in xrange(11)],  # once slicing works, use self
            name='{}-regular sequence'.format(self.parent().k),
            opening_delimiter='', closing_delimiter='',
            preview=10)


    def info(self):
        r"""
        Displays the matrices of the `k`-linear representation, the initial
        vector and the selection vector.

        OUTPUT:

        Nothing; printing to standard output.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:      initial=vector([0, 1]), selection=vector([1, 0])).info()
            matrices:
            (
            [1 0]  [ 0 -1]
            [0 1], [ 1  2]
            )
            initial:
            (0, 1)
            selection:
            (1, 0)
        """
        from sys import displayhook
        print('matrices:')
        displayhook(self.matrices)
        print('initial:')
        displayhook(self.initial)
        print('selection:')
        displayhook(self.selection)


    @cached_method
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
            ....:          initial=vector([0, 1]), selection=vector([1, 0]))
            sage: S[7]
            3
        """
        result = self._product_of_matrices_(n)
        if self.initial is not None:
            result = self.initial * result
        if self.selection is not None:
            result = result * self.selection
        return self.output_function(result)


    @cached_method
    def _product_of_matrices_(self, m):
        r"""
        Return the product of matrices according to the `k`-ary
        digit expansion of `m`.

        INPUT:

        - ``m`` -- a nonnegative integer.

        OUTPUT:

        A matrix.

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Seq2((M0, M1))
            sage: S._product_of_matrices_(0) == M0
            True
            sage: S._product_of_matrices_(1) == M1
            True
            sage: S._product_of_matrices_(3) == M1^2
            True

        ::

            sage: S._product_of_matrices_(-1)
            Traceback (most recent call last):
            ...
            ValueError: m=-1 is not a nonnegative integer.
         """
        k = self.parent().k
        if m < 0:
            raise ValueError('m={} is not a nonnegative integer.'.format(m))
        if 0 <= m < k:
            return self.matrices[m]
        n = m // k
        r = m - n*k
        return self.matrices[r] * self._product_of_matrices_(n)


from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

class kRegularSequenceSpace(UniqueRepresentation, Parent):

    Element = kRegularSequence

    def __init__(self, k, base, category=None):
        r"""
        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2
            Space of 2-regular sequences over Integer Ring
            sage: Seq2.category()
            Category of sets
        """
        from sage.categories.sets_cat import Sets
        self.k = k
        super(kRegularSequenceSpace, self).__init__(
            category=category or Sets(), base=base)


    def _repr_(self):
        return 'Space of {}-regular sequences over {}'.format(self.k, self.base())

