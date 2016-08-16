r"""
Recognizable Series


Various
=======

REFERENCES:

.. [BR2010] Jean Berstel, Christophe Reutenauer,
   *Noncommutative Rational Series With Applications*,
   Cambridge, 2010.

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


class PrefixClosedSet(object):

    def __init__(self, alphabet=None, words=None):
        r"""
        A prefix-closed set.

        Creation of this prefix-closed sets is interactive
        iteratively.

        INPUT:

        - ``alphabet`` -- finite words over this ``alphabet``
          will be created.

        - ``words`` -- specify the finite words directly
          (instead of via ``alphabet``).

        EXAMPLES::

            sage: from sage.combinat.recognizable_series import PrefixClosedSet
            sage: P = PrefixClosedSet(alphabet=[0, 1]); P
            [word: ]

        See :meth:`populate_interactive` for further examples.

        TESTS::

            sage: P = PrefixClosedSet(
            ....:         words=Words([0, 1], infinite=False)); P
            [word: ]
        """
        if alphabet is not None:
            from sage.combinat.words.words import Words
            self.words = Words(alphabet, infinite=False)
        else:
            self.words = words
        self.elements = [self.words([])]


    def __repr__(self):
        r"""
        A representation string of this prefix-closed set

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.combinat.recognizable_series import PrefixClosedSet
            sage: P = PrefixClosedSet(alphabet=[0, 1])
            sage: repr(P)  # indirect doctest
            '[word: ]'
        """
        return repr(self.elements)


    def add(self, w, check=True):
        r"""
        Add a word to this prefix-closed set.

        INPUT:

        - ``w`` -- a word.

        - ``check`` -- (default: ``True``) if set, then it is verified
          whether all proper prefixes of ``w`` are already in this
          prefix-closed set.

        OUTPUT:

        Nothing, but a
        :python:`RuntimeError<library/exceptions.html#exceptions.ValueError>`
        is raised if the check fails.

        EXAMPLES::

            sage: from sage.combinat.recognizable_series import PrefixClosedSet
            sage: P = PrefixClosedSet(alphabet=[0, 1])
            sage: W = P.words
            sage: P.add(W([0])); P
            [word: , word: 0]
            sage: P.add(W([0, 1])); P
            [word: , word: 0, word: 01]
            sage: P.add(W([1, 1]))
            Traceback (most recent call last):
            ...
            ValueError: Cannot add as not all prefixes of 11 are included yet.
        """
        if check and any(p not in self.elements
                         for p in w.prefixes_iterator()
                         if p != w):
            raise ValueError('Cannot add as not all prefixes of '
                             '{} are included yet.'.format(w))
        self.elements.append(w)


    def populate_interactive(self):
        r"""
        Return an iterator over possible new elements.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from sage.combinat.recognizable_series import PrefixClosedSet
            sage: P = PrefixClosedSet(alphabet=[0, 1]); P
            [word: ]
            sage: for n, p in enumerate(P.populate_interactive()):
            ....:     print('{}?'.format(p))
            ....:     if n in (0, 2, 3, 5):
            ....:         P.add(p)
            ....:         print('...added')
            0?
            ...added
            1?
            00?
            ...added
            01?
            ...added
            000?
            001?
            ...added
            010?
            011?
            0010?
            0011?
            sage: P.elements
            [word: , word: 0, word: 00, word: 01, word: 001]
        """
        n = 0
        it = self.words.iterate_by_length(1)
        while n < len(self.elements):
            try:
                nn = next(it)
                yield self.elements[n] + nn #next(it)
            except StopIteration:
                n += 1
                it = self.words.iterate_by_length(1)


    def prefix_set(self):
        r"""
        Return the prefix set corresponding to this prefix
        closed set.

        OUTPUT:

        A list.

        EXAMPLES::

            sage: from sage.combinat.recognizable_series import PrefixClosedSet
            sage: P = PrefixClosedSet(alphabet=[0, 1]); P
            [word: ]
            sage: for n, p in enumerate(P.populate_interactive()):
            ....:     if n in (0, 1, 2, 4, 6):
            ....:         P.add(p)
            sage: P
            [word: , word: 0, word: 1, word: 00, word: 10, word: 000]
            sage: P.prefix_set()
            [word: 01, word: 11, word: 001, word: 100,
             word: 101, word: 0000, word: 0001]
        """
        return [p + a
                for p in self.elements
                for a in self.words.iterate_by_length(1)
                if p + a not in self.elements]


class RecognizableSeries(Element):

    def __init__(self, parent, mu, left=None, right=None):
        r"""
        A recognizable series.

        - ``parent`` -- an instance of :class:`RecognizableSeriesSpace`.

        - ``mu`` -- a family of square matrices, all of which have the
          same dimension. The indices of this family are the alphabet.
          ``mu`` may be a list or tuple of the same cardinality as the
          alphabet as well. See :meth:`mu` for more details.

        - ``left`` -- (default: ``None``) a vector.  When evaluating a
          coefficient, this vector is multiplied from the left to the
          matrix obtained from :meth:`mu` applying on a word. If
          ``None``, then this multiplication is skipped.
          See :meth:`left` for more details.

        - ``right`` -- (default: ``None``) a vector.  When evaluating a
          coefficient, this vector is multiplied from the right to the
          matrix obtained from :meth:`mu` applying on a word. If
          ``None``, then this multiplication is skipped.
          See :meth:`right` for more details.

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:     vector([0, 1]), vector([1, 0]),
            ....:     transpose=True)
            [1] + 3*[01] + [10] + 5*[11] + 9*[001] + 3*[010] + ...

        .. SEEALSO::

            :doc:`recognizable series <recognizable_series>`,
            :class:`RecognizableSeriesSpace`.

        TESTS::

            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Rec(mu=(M0, M1))  # not tested TODO
        """
        super(RecognizableSeries, self).__init__(parent=parent)

        from sage.sets.family import Family

        A = self.parent().alphabet()
        if isinstance(mu, (list, tuple)):
            if len(A) != len(mu):
                raise ValueError('TODO')
            mu = dict(zip(A, mu))
        mu = Family(mu)
        # TODO some check that alphabet is correct
        #if A != self.mu:
        #    raise ValueError('TODO')

        #self.k = len(self.matrices)
        #self.d = self.matrices[0].nrows()
        #if not all(M.dimensions() == (self.d, self.d) for M in self.matrices):
        #    raise ValueError  # TODO

        self._left_ = left
        self._mu_ = mu
        self._right_ = right


    @property
    def mu(self):
        r"""
        When evaluating a coefficient, this is applied on each letter
        of a word; the result is a matrix.
        This extends :meth:`mu` to words over the parent's
        :meth:`~RecognizableSeriesSpace.alphabet`.
        """
        return self._mu_


    @property
    def left(self):
        r"""
        When evaluating a coefficient, this vector is multiplied from
        the left to the matrix obtained from :meth:`mu` applied on a
        word.
        """
        return self._left_


    @property
    def right(self):
        r"""
        When evaluating a coefficient, this vector is multiplied from
        the right to the matrix obtained from :meth:`mu` applied on a
        word.
        """
        return self._right_


    def _repr_(self, latex=False):
        r"""
        A representation string for this recognizable series.

        INPUT:

        - ``latex`` -- (default: ``False``) a boolean. If set, then
          LaTeX-output is returned.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:         vector([0, 1]), vector([1, 0]), transpose=True)
            sage: repr(S)  # indirect doctest
            '[1] + 3*[01] + [10] + 5*[11] + 9*[001] + 3*[010] + ...'
        """
        if not self:
            return '0'

        from itertools import islice

        if latex:
            from sage.misc.latex import latex as latex_repr
            f = latex_repr
            times = ' '
        else:
            f = repr
            times = '*'

        def summand(c, w):
            if c == 1:
                return f(w)
            return f(c) + times + f(w)

        s = ' + '.join(summand(c, w)
                       for c, w in islice(self, 10))
        s = s.replace('+ -', '- ')
        return s + ' + ...'


    def _latex_(self):
        r"""
        A LaTeX-representation string for this recognizable series.

        OUTPUT:

        A string.

        TESTS::

        """
        return self._repr_(latex=True)


    @cached_method
    def __getitem__(self, w):
        r"""
        Return the coefficient to word `w` of this series.

        INPUT:

        - ``w`` -- a word over the parent's
          :meth:`~RecognizableSeriesSpace.alphabet`.

        OUTPUT:

        An element in the parent's
        :meth:`~RecognizableSeriesSpace.coefficients`.

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: W = Rec.indices()
            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:          left=vector([0, 1]), right=vector([1, 0]))
            sage: S[W(7.digits(2))]
            3
        """
        result = self._mu_of_word_(w)
        if self.left is not None:
            result = self.left * result
        if self.right is not None:
            result = result * self.right
        return result


    @cached_method
    def _mu_of_empty_word_(self):
        r"""
        Return :meth:`mu` applied on the empty word.

        OUTPUT:

        A matrix.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: W = Rec.indices()
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Rec({W([0]): M0, W([1]): M1})
            sage: S._mu_of_empty_word_()
            [1 0]
            [0 1]
            sage: I = Matrix([[1, 0], [0, 1]])
            sage: T = Rec({W([]): I, W([0]): M0, W([1]): M1})
            sage: T._mu_of_empty_word_()
            [1 0]
            [0 1]
            sage: _ is I
            True
        """
        eps = self.parent().indices()()
        try:
            return self.mu[eps]
        except KeyError:
            return next(iter(self.mu)).parent().one()


    @cached_method
    def _mu_of_word_(self, w):
        r"""
        Return :meth:`mu` applied on the word `w`.

        INPUT:

        - ``w`` -- a word over the parent's
          :meth:`~RecognizableSeriesSpace.alphabet`.

        OUTPUT:

        A matrix.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: W = Rec.indices()
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Rec((M0, M1))
            sage: S._mu_of_word_(W([0])) == M0
            True
            sage: S._mu_of_word_(W([1])) == M1
            True
            sage: S._mu_of_word_(W(3.digits(2))) == M1^2
            True

        ::

            sage: S._mu_of_word_(-1)
            Traceback (most recent call last):
            ...
            ValueError: m=-1 is not a nonnegative integer.
        """
        if w not in self.parent().indices():
            raise ValueError('TODO')
        from sage.misc.misc_c import prod
        return prod(tuple(self.mu[a] for a in w), z=self._mu_of_empty_word_())


    def __iter__(self):
        r"""
        Return an iterator over pairs ``(coefficient, basis(index))``.

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:         left=vector([0, 1]), right=vector([1, 0]))
            sage: from itertools import islice
            sage: tuple(islice(S, 10))
            ((1, [1]),
             (1, [01]),
             (1, [10]),
             (2, [11]),
             (1, [001]),
             (1, [010]),
             (2, [011]),
             (1, [100]),
             (2, [101]),
             (2, [110]))

        TESTS::

            sage: it = iter(S)
            sage: iter(it) is it
            True
            sage: iter(S) is not it
            True
        """
        if not self:
            return iter([])
        A = self.parent()._algebra_
        B = A.basis()
        return iter((self[w], B[w])
                    for w in self.parent().indices() if self[w] != 0)


    def is_trivial_zero(self):
        r"""
        Return whether this recognizable series is trivially equal to
        zero (without any :meth:`minimization <minmized>`).

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:     left=vector([0, 1]), right=vector([1, 0])).is_trivial_zero()
            False
            sage: Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:     left=vector([0, 0]), right=vector([1, 0])).is_trivial_zero()
            True
            sage: Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:     left=vector([0, 1]), right=vector([0, 0])).is_trivial_zero()
            True

        The following two differ in the coefficient of the empty word::

            sage: Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:     left=vector([0, 1]), right=vector([1, 0])).is_trivial_zero()
            True
            sage: Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:     left=vector([1, 1]), right=vector([1, 1])).is_trivial_zero()
            False

        TESTS::

            sage: Rec.zero().is_trivial_zero()
            True

        The following are nonzero as the coefficient of the empty word
        is nonzero::

            sage: Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:     left=None, right=vector([1, 0])).is_trivial_zero()
            False
            sage: Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:     left=None, right=None).is_trivial_zero()
            False
            sage: Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:     left=vector([0, 1]), right=None).is_trivial_zero()
            False

        The following is zero, but not trivially zero::

            sage: S = Rec((Matrix([[1, 0], [0, 0]]), Matrix([[1, 0], [0, 0]])),
            ....:         left=vector([0, 1]), right=vector([1, 0]))
            sage: S.is_trivial_zero()
            False
            sage: S.is_zero()
            True
        """
        return (self.left is not None and not self.left) or \
            (self.right is not None and not self.right) or \
            (all(not self.mu[a] for a in self.parent().alphabet()) and
             not self[self.parent().indices()()])



    def __nonzero__(self):
        r"""
        Return whether this recognizable series is nonzero.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: bool(Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:          left=vector([0, 1]), right=vector([1, 0])))
            False
            sage: bool(Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:          left=vector([0, 1]), right=vector([1, 0])))
            False
            sage: bool(Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:          left=vector([0, 0]), right=vector([1, 0])))
            False
            sage: bool(Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:          left=vector([0, 1]), right=vector([0, 0])))
            False

        ::

            sage: S = Rec((Matrix([[1, 0], [0, 0]]), Matrix([[1, 0], [0, 0]])),
            ....:         left=vector([0, 1]), right=vector([1, 0]))
            sage: bool(S)
            False
            sage: S  # not tested
        """
        if self.is_trivial_zero():
            return False
        try:
            M = self.minimized()
        except ValueError:
            pass
        else:
            if M.is_trivial_zero():
                return False
        return True


    def transposed(self):
        r"""
        Return the transposed series.

        OUTPUT:

        A :class:`RecognizableSeries`.

        Each of the ``matrices`` is transposed. Additionally
        the vectors ``left`` and ``right`` are switched.

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:         vector([0, 1]), vector([1, 0]), transpose=True)
            sage: S
            [1] + 3*[01] + [10] + 5*[11] + 9*[001] + 3*[010]
                + 15*[011] + [100] + 11*[101] + 5*[110] + ...
            sage: S.mu[0], S.mu[1], S.left, S.right
            (
            [3 0]  [ 0  1]
            [6 1], [-6  5], (1, 0), (0, 1)
            )
            sage: T = S.transposed()
            sage: T
            [1] + [01] + 3*[10] + 5*[11] + [001] + 3*[010]
                + 5*[011] + 9*[100] + 11*[101] + 15*[110] + ...
            sage: T.mu[0], T.mu[1], T.left, T.right
            (
            [3 6]  [ 0 -6]
            [0 1], [ 1  5], (0, 1), (1, 0)
            )
        """
        def tr(M):
            try:
                return M.transpose()
            except AttributeError:
                return M
        return self.parent()(self.mu.map(tr),
                             left=tr(self.right),
                             right=tr(self.left))


    def minimized(self):
        r"""

        EXAMPLES::

            sage: from itertools import islice, izip
            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])

            sage: S = Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:         vector([0, 1]), vector([1, 0]),
            ....:         transpose=True)
            sage: S
            [1] + 3*[01] + [10] + 5*[11] + 9*[001] + 3*[010]
                + 15*[011] + [100] + 11*[101] + 5*[110] + ...
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            (
            [3 0]  [ 0  1]
            [6 1], [-6  5], (1, 0), (0, 1)
            )
            sage: all(c == d and v == w
            ....:     for (c, v), (d, w) in islice(izip(iter(S), iter(M)), 20))
            True

            sage: S = Rec((Matrix([[2, 0], [1, 1]]), Matrix([[2, 0], [2, 1]])),
            ....:         vector([1, 0]), vector([1, 1]))
            sage: S
            [] + 2*[0] + 2*[1] + 4*[00] + 4*[01] + 4*[10] + 4*[11]
               + 8*[000] + 8*[001] + 8*[010] + ...
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([2], [2], (1), (1))
            sage: all(c == d and v == w
            ....:     for (c, v), (d, w) in islice(izip(iter(S), iter(M)), 20))
            True
        """
        if self.left is None or self.right is None:
            raise ValueError("Cannot minmize as 'left' or 'right' is None.")
        return self._minimized_right_()._minimized_left_()


    def _minimized_right_(self):
        return self.transposed()._minimized_left_().transposed()


    def _minimized_left_(self):
        r"""

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[0, 0], [0, 0]]), Matrix([[0, 0], [0, 0]])),
            ....:         vector([1, 1]), vector([1, 1]))
            sage: M = S._minimized_left_()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([0], [0], (1), (2))
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([0], [0], (1), (2))

        ::

            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:         vector([1, -1]), vector([1, 1]))._minimized_left_()
            sage: S.mu[0], S.mu[1], S.left, S.right
            ([1], [1], (1), (0))
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([], [], (), ())

            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:         vector([1, 1]), vector([1, -1]))
            sage: M = S._minimized_left_()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([1], [1], (1), (0))
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([], [], (), ())

            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[1, 0], [0, 1]])),
            ....:         left=vector([0, 1]), right=vector([1, 0]))
            sage: M = S._minimized_left_()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([1], [1], (1), (0))
            sage: M = S.minimized()
            sage: M.mu[0], M.mu[1], M.left, M.right
            ([], [], (), ())
        """
        from sage.matrix.constructor import Matrix
        from sage.modules.free_module_element import vector
        from sage.rings.integer_ring import ZZ

        pcs = PrefixClosedSet(self.parent().indices())
        left = self.left * self._mu_of_word_(pcs.elements[0])
        if left.is_zero():
            return self.parent().zero()
        Left = [left]
        for p in pcs.populate_interactive():
            left = self.left * self._mu_of_word_(p)
            try:
                Matrix(Left).solve_left(left)
            except ValueError:
                # no solution found
                pcs.add(p)
                Left.append(left)
        P = pcs.elements
        C = pcs.prefix_set()

        ML = Matrix(Left)

        def alpha(c):
            return ML.solve_left(self.left * self._mu_of_word_(c))

        def mu_prime_entry(a, p, q, iq):
            c = p + a
            if c == q:
                return ZZ(1)
            elif c in C:
                return alpha(c)[iq]
            else:
                return ZZ(0)

        mu_prime = []
        for a in self.parent().alphabet():
            a = self.parent().indices()([a])
            M = [[mu_prime_entry(a, p, q, iq) for iq, q in enumerate(P)]
                 for p in P]
            mu_prime.append(Matrix(M))

        left_prime = vector([ZZ(1)] + (len(P)-1)*[ZZ(0)])
        if self.right is None:
            right_prime = None
        else:
            right_prime = vector(self[p] for p in P)

        return self.parent().element_class(
            self.parent(), mu_prime, left_prime, right_prime)


from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

class RecognizableSeriesSpace(UniqueRepresentation, Parent):

    Element = RecognizableSeries


    @staticmethod
    def __classcall__(cls,
                      coefficients=None, alphabet=None, 
                      indices=None, algebra=None,
                      category=None):
        r"""
        Normalizes the input in order to ensure a unique
        representation.

        For more information see :class:`ReconizableSeriesSpace`.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec.category()
            Category of set algebras over Integer Ring

        """
        if algebra is None:
            if indices is None:
                if alphabet is None:
                    raise ValueError('TODO')
                from sage.combinat.words.words import Words
                indices = Words(alphabet, infinite=False)
                from sage.combinat.words.word_options import WordOptions
                WordOptions(identifier='')
            if coefficients is None:
                raise ValueError('TODO')
            algebra = indices.algebra(coefficients)
            def key_to_cmp(a, b, key):
                return (key(a) > key(b)) - (key(a) < key(b))
            algebra.print_options(
                prefix='',
                generator_cmp=lambda a, b: key_to_cmp(a, b, key=lambda k: (len(k), k)))

        from sage.categories.sets_cat import Sets
        category = category or algebra.category()

        return super(RecognizableSeriesSpace, cls).__classcall__(
            cls, algebra, category)


    def __init__(self, algebra, category):
        r"""
        The space of recognizable series on the given alphabet and
        with the given coefficients.

        INPUT:

        - ``coefficients`` -- a (semi-)ring.

        - ``alphabet`` -- a tuple, list or totally ordered set.
          If specified, then the ``indices`` are the
          finite words over this ``alphabet``.

        - ``indices`` -- a SageMath parent.

        - ``algebra`` -- a SageMath parent.
          If specified, then ``coefficients``
          and ``indices`` are determined by this ``algebra``.

        - ``category`` -- (default: ``None``) the category of this
          space. If ``None``, then the category of the ``algebra``
          is used.

        EXAMPLES::

            sage: S1 = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S1
            Space of recognizable series on {0, 1} with coefficients in Integer Ring
            sage: S2 = RecognizableSeriesSpace(coefficients=ZZ, alphabet=[0, 1])
            sage: S2
            Space of recognizable series on {0, 1} with coefficients in Integer Ring
            sage: S3 = RecognizableSeriesSpace(ZZ, indices=Words([0, 1], infinite=False))
            sage: S3
            Space of recognizable series on {0, 1} with coefficients in Integer Ring
            sage: S4 = RecognizableSeriesSpace(algebra=Words([0, 1], infinite=False).algebra(ZZ))
            sage: S4
            Space of recognizable series on {0, 1} with coefficients in Integer Ring
            sage: S1 is S2 is S3 is S4
            True

        .. SEEALSO::

            :doc:`recognizable series <recognizable_series>`,
            :class:`RecognizableSeries`.
        """
        self._algebra_ = algebra
        super(RecognizableSeriesSpace, self).__init__(
            category=category, base=algebra.base())


    def alphabet(self):
        r"""
        Return the alphabet of this recognizable series space.

        OUTPUT:

        A totally ordered set.

        EXAMPLES::

            sage: RecognizableSeriesSpace(ZZ, [0, 1]).alphabet()
            {0, 1}

        TESTS::

            sage: type(RecognizableSeriesSpace(ZZ, [0, 1]).alphabet())
            <class 'sage.sets.totally_ordered_finite_set.TotallyOrderedFiniteSet_with_category'>
        """
        return self._algebra_.indices().alphabet()


    def indices(self):
        r"""
        Return the indices of the recognizable series.

        OUTPUT:

        The set of finite words over the alphabet.

        EXAMPLES::

            sage: RecognizableSeriesSpace(ZZ, [0, 1]).indices()
            Finite words over {0, 1}
        """
        return self._algebra_.indices()


    def coefficients(self):
        r"""
        Return the coefficients of this recognizable series space.

        OUTPUT:

        A (semi-)ring.

        EXAMPLES::

            sage: RecognizableSeriesSpace(ZZ, [0, 1]).coefficients()
            Integer Ring
        """
        return self.base()


    def _repr_(self):
        r"""
        Return a representation string of this recognizable sequence
        space.

        OUTPUT:

        A string.

        TESTS::

            sage: repr(RecognizableSeriesSpace(ZZ, [0, 1]))  # indirect doctest
            'Space of recognizable series on {0, 1} with coefficients in Integer Ring'
        """
        return 'Space of recognizable series on {} ' \
               'with coefficients in {}'.format(self.alphabet(),
                                                self.coefficients())


    def _element_constructor_(self, mu,
                              left=None, right=None,
                              transpose=False):
        r"""

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec.zero()
            0
            sage: type(_)
            <class 'sage.combinat.recognizable_series.RecognizableSeriesSpace_with_category.element_class'>
        """
        if isinstance(mu, int) and mu == 0:
            from sage.matrix.constructor import Matrix
            from sage.modules.free_module_element import vector
            from sage.sets.family import Family

            return self.element_class(
                self, Family(self.alphabet(), lambda a: Matrix()),
                vector([]), vector([]))

        element = self.element_class(self, mu, left, right)
        if transpose:
            return element.transposed()
        else:
            return element

