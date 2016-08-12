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


class RecognizableSeries(Element):

    def __init__(self, parent, mu, left=None, right=None, transpose=False):
        r"""
        A recognizable series.

        TODO-INPUT:

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

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:     vector([0, 1]), vector([1, 0]),
            ....:     transpose=True)
            [1] + 3*[01] + [10] + 5*[11] + 9*[001] + 3*[010] + ...

        Using an output function::

            sage: Rec = RecognizableSeriesSpace(
            ....:           ZZ, [0, 1], output_function=lambda o: o[0, 0])
            sage: Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:     Matrix([[0, 1]]), Matrix([[1], [0]]),
            ....:     transpose=True)
            [1] + 3*[01] + [10] + 5*[11] + 9*[001] + 3*[010] + ...

        .. SEEALSO::

            :doc:`k-regular sequence <k_regular_sequence>`,
            :class:`kRegularSequenceSpace`.

        TESTS::

            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Rec(mu=(M0, M1))  # not tested
        """
        super(RecognizableSeries, self).__init__(parent=parent)

        from sage.sets.family import Family

        def tr(M):
            try:
                return M.transpose()
            except AttributeError:
                return M

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

        if not transpose:
            self.left = left
            self.mu = mu
            self.right = right
        else:
            self.left = tr(right)
            self.mu = mu.map(tr)
            self.right = tr(left)


    def _repr_(self):
        r"""
        Return a representation string of this recognizable series.

        OUTPUT:

        A string

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: S = Rec((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:         vector([0, 1]), vector([1, 0]), transpose=True)
            sage: repr(S)  # indirect doctest
            '[1] + 3*[01] + [10] + 5*[11] + 9*[001] + 3*[010] + ...'
        """
        from itertools import islice
        A = self.parent()._algebra_
        B = A.basis()
        return repr(sum((c * w for c, w in islice(self, 10)),
                        A.zero())) + ' + ...'


    @cached_method
    def __getitem__(self, w):
        r"""
        Return the coefficient to word `w` of this sequence.

        INPUT:

        - ``n`` -- a nonnegative integer.

        OUTPUT:

        An element of the universe of the sequence.

        EXAMPLES::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: W = Rec.indices()
            sage: S = Rec((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:          left=vector([0, 1]), right=vector([1, 0]))
            sage: S[W(7.digits(2))]
            3
        """
        result = self._product_of_matrices_(w)
        if self.left is not None:
            result = self.left * result
        if self.right is not None:
            result = result * self.right
        return self.parent()._output_function_(result)


    @cached_method
    def _mu_empty_(self):
        eps = self.parent().indices()()
        try:
            return self.mu[eps]
        except KeyError:
            return next(iter(self.mu.values())).parent().one()


    @cached_method
    def _product_of_matrices_(self, w):
        r"""
        Return the product of matrices according to the `k`-ary
        digit expansion of `m`.

        INPUT:

        - ``m`` -- a nonnegative integer.

        OUTPUT:

        A matrix.

        TESTS::

            sage: Rec = RecognizableSeriesSpace(ZZ, [0, 1])
            sage: W = Rec.indices()
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Rec((M0, M1))
            sage: S._product_of_matrices_(W([0])) == M0
            True
            sage: S._product_of_matrices_(W([1])) == M1
            True
            sage: S._product_of_matrices_(W(3.digits(2))) == M1^2
            True

        ::

            sage: S._product_of_matrices_(-1)
            Traceback (most recent call last):
            ...
            ValueError: m=-1 is not a nonnegative integer.
        """
        if w not in self.parent().indices():
            raise ValueError('TODO')
        from sage.misc.misc_c import prod
        return prod(tuple(self.mu[a] for a in w), z=self._mu_empty_())


    def __iter__(self):
        r"""
        Return an iterator.

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
        # TODO self == 0: return iter()
        A = self.parent()._algebra_
        B = A.basis()
        return iter((self[w], B[w])
                    for w in self.parent().indices() if self[w] != 0)


from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

class RecognizableSeriesSpace(UniqueRepresentation, Parent):

    Element = RecognizableSeries

    @staticmethod
    def __classcall__(cls,
                      coefficients=None, alphabet=None, 
                      indices=None, algebra=None,
                      output_function=None, category=None):
        r"""
        Normalizes the input in order to ensure a unique
        representation.

        For more information see :class:`kRegularSequenceSpace`.

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
            cls, algebra, output_function, category)


    def __init__(self, algebra, output_function, category):
        r"""
        The space of `k`-regular Sequences over the given ``universe``.

        INPUT:

        - ``k`` -- an integer at least `2` specifying the base.

        - ``universe`` -- (default: ``None``) a object (e.g. a SageMath parent)
          in which the entries of a sequence live.
          If ``None``, then the integer ring `\ZZ` is used.

        - ``category`` -- (default: ``None``) the category of the
          sequence space. If ``None``, then the category of
          :class:`~sage.categories.sets_cat.Sets` is used.

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

            :doc:`k-regular sequence <k_regular_sequence>`,
            :class:`kRegularSequence`.
        """
        self._algebra_ = algebra
        if output_function is None:
            self._output_function_ = lambda o: o
        else:
            self._output_function_ = output_function
        super(RecognizableSeriesSpace, self).__init__(
            category=category, base=algebra.base())


    def alphabet(self):
        return self._algebra_.indices().alphabet()

    def indices(self):
        return self._algebra_.indices()

    def coefficients(self):
        return self.base()

    def _repr_(self):
        r"""
        Return a representation string of this `k`-regular sequence space.

        OUTPUT:

        A string

        TESTS::

            sage: repr(RecognizableSeriesSpace(ZZ, [0, 1]))  # indirect doctest
            'Space of recognizable series on {0, 1} with coefficients in Integer Ring'
        """
        return 'Space of recognizable series on {} ' \
               'with coefficients in {}'.format(self.alphabet(),
                                                 self.coefficients())

