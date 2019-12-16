r"""
`k`-regular Sequences

An introduction and formal definition of `k`-regular sequences can be
found, for example, on the :wikipedia:`k-regular_sequence` or in
[AS2003]_.


.. WARNING::

    As this code is experimental, warnings are thrown when a
    `k`-regular sequence space is created for the first time in a
    session (see :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: Seq2 = kRegularSequenceSpace(2, ZZ)
        doctest:...: FutureWarning: This class/method/function is
        marked as experimental. It, its functionality or its interface
        might change without a formal deprecation.
        See http://trac.sagemath.org/21202 for details.


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

    :mod:`recognizable series <sage.combinat.recognizable_series>`,
    :mod:`sage.rings.cfinite_sequence`,
    :mod:`sage.combinat.binary_recurrence_sequences`.

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
from __future__ import absolute_import

from .recognizable_series import RecognizableSeries
from .recognizable_series import RecognizableSeriesSpace
from sage.misc.cachefunc import cached_method


class kRegularSequence(RecognizableSeries):

    def __init__(self, parent, mu, left=None, right=None):
        r"""
        A `k`-regular sequence.

        INPUT:

        - ``parent`` -- an instance of :class:`kRegularSequenceSpace`.

        - ``mu`` -- a family of square matrices, all of which have the
          same dimension. The indices of this family are `0,...,k-1`.
          ``mu`` may be a list or tuple of cardinality `k`
          as well. See
          :meth:`~sage.combinat.recognizable_series.RecognizableSeries.mu`
          for more details.

        - ``left`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the left to the matrix product. If ``None``, then this
          multiplication is skipped.

        - ``right`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the right to the matrix product. If ``None``, then this
          multiplication is skipped.

        When created via the parent :class:`kRegularSequenceSpace`, then
        the following option is available.

        - ``transpose`` -- (default: ``False``) a boolean. If set, then
          each of the matrices in
          :meth:`mu <sage.combinat.recognizable_series.RecognizableSeries.mu>`
          is transposed. Additionally the vectors
          :meth:`left <sage.combinat.recognizable_series.RecognizableSeries.left>`
          and
          :meth:`right <sage.combinat.recognizable_series.RecognizableSeries.right>`
          are switched.
          (This is done by calling :meth:`~sage.combinat.recognizable_series.RecognizableSeries.transposed`.)

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

        A string.

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


    @cached_method
    def __getitem__(self, n, **kwds):
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

            sage: S[-1]
            Traceback (most recent call last):
            ...
            OverflowError: can't convert negative value to unsigned char

        ::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: W = Seq2.indices()
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Seq2((M0, M1), [0, 1], [1, 1])
            sage: S._mu_of_word_(W(0.digits(2))) == M0
            True
            sage: S._mu_of_word_(W(1.digits(2))) == M1
            True
            sage: S._mu_of_word_(W(3.digits(2))) == M1^2
            True
        """
        return self.coefficient_of_word(self.parent()._n_to_index_(n), **kwds)


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
    r"""
    The space of `k`-regular Sequences over the given ``coefficients``.

    INPUT:

    - ``k`` -- an integer at least `2` specifying the base.

    - ``coefficients`` -- a (semi-)ring. If not specified (``None``),
      then the integer ring is used.

    - ``category`` -- (default: ``None``) the category of this
      space.

    EXAMPLES::

        sage: kRegularSequenceSpace(2, ZZ)
        Space of 2-regular sequences over Integer Ring
        sage: kRegularSequenceSpace(3, ZZ)
        Space of 3-regular sequences over Integer Ring

    .. SEEALSO::

        :doc:`k-regular sequence <k_regular_sequence>`,
        :class:`kRegularSequence`.
    """

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


    def __init__(self, k, *args, **kwds):
        r"""
        See :class:`kRegularSequenceSpace` for details.

        INPUT:

        - ``k`` -- an integer at least `2` specifying the base.

        Other input arguments are passed on to
        :meth:`~sage.combinat.recognizable_series.RecognizableSeriesSpace.__init__`.

        TESTS::

            sage: kRegularSequenceSpace(2, ZZ)
            Space of 2-regular sequences over Integer Ring
            sage: kRegularSequenceSpace(3, ZZ)
            Space of 3-regular sequences over Integer Ring

        .. SEEALSO::

            :doc:`k-regular sequence <k_regular_sequence>`,
            :class:`kRegularSequence`.
        """
        self.k = k
        super(kRegularSequenceSpace, self).__init__(*args, **kwds)


    def _repr_(self):
        r"""
        Return a representation string of this `k`-regular sequence space.

        OUTPUT:

        A string.

        TESTS::

            sage: repr(kRegularSequenceSpace(2, ZZ))  # indirect doctest
            'Space of 2-regular sequences over Integer Ring'
        """
        return 'Space of {}-regular sequences over {}'.format(self.k, self.base())


    def _n_to_index_(self, n):
        r"""
        Convert `n` to an index usable by the underlying
        recognizable series.

        INPUT:

        - ``n`` -- a nonnegative integer.

        OUTPUT:

        A word.

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2._n_to_index_(6)
            word: 011
            sage: Seq2._n_to_index_(-1)
            Traceback (most recent call last):
            ...
            OverflowError: can't convert negative value to unsigned char
        """
        from sage.rings.integer_ring import ZZ
        n = ZZ(n)
        W = self.indices()
        return W(n.digits(self.k))


    def _parse_recursions_(self, equations, function, var, n0=0):
        r"""Parse recursion equations as admissible in :meth:`~.recursions`.

        INPUT:

        - ``equations`` -- see :meth:`~recursions`.

        - ``function`` -- see :meth:`~recursions`.

        - ``var`` -- see :meth:`~recursions`.

        OUTPUT:

        A namedtuple ``recursion_rules`` containing all the
        significant information obtained from ``equations``.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: Seq2._parse_recursions_([
            ....:     f(4*n) == f(2*n) + 2*f(2*n + 1) + 3*f(2*n - 2),
            ....:     f(4*n + 1) == 4*f(2*n) + 5*f(2*n + 1) + 6*f(2*n - 2),
            ....:     f(4*n + 2) == 7*f(2*n) + 8*f(2*n + 1) + 9*f(2*n - 2),
            ....:     f(4*n + 3) == 10*f(2*n) + 11*f(2*n + 1) + 12*f(2*n - 2),
            ....:     f(0) == 1, f(1) == 2, f(2) == 1], f, n, 42)
            recursion_rules(M=2, m=1, l=-2, u=1, ll=-6, uu=3, dim=11,
            coeffs={(0, 1): 2, (0, 0): 1, (3, 1): 11, (3, 0): 10, (2, -2): 9,
            (2, 1): 8, (2, 0): 7, (3, -2): 12, (0, -2): 3, (1, 0): 4, (1, -2): 6,
            (1, 1): 5}, start_values={0: 1, 1: 2, 2: 1}, n0=42)

        Stern--Brocot Sequence::

            sage: Seq2._parse_recursions_([f(2*n) == f(n),
            ....:    f(2*n + 1) == f(n) + f(n + 1), f(0) == 0,
            ....:    f(1) == 1, f(2) == 1], f, n)
            recursion_rules(M=1, m=0, l=0, u=1, ll=0, uu=2, dim=3,
            coeffs={(1, 0): 1, (0, 0): 1, (1, 1): 1},
            start_values={0: 0, 1: 1, 2: 1}, n0=0)

        TESTS:

            The following tests check that the equations are well-formed::

                sage: Seq2._parse_recursions_([], f, n)
                Traceback (most recent call last):
                ...
                ValueError: List of recursion equations is empty.

            ::

                sage: Seq2._parse_recursions_([f(4*n + 1)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: f(4*n + 1) is not an equation with ==.

            ::

                sage: Seq2._parse_recursions_([f(2*n) + 1 == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: f(2*n) + 1 is not an evaluation of f.

            ::

                sage: Seq2._parse_recursions_([f(2*n, 5) == 3], f, n)
                Traceback (most recent call last):
                ...
                ValueError: f(2*n, 5) does not have one argument.

            ::

                sage: Seq2._parse_recursions_([f(1/n + 1) == f(n)], f, n)
                Traceback (most recent call last):
                ....:
                ValueError: 1/n + 1 is not a polynomial in n.

            ::

                sage: Seq2._parse_recursions_([f(4*n^2) == f(2*n^2)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 4*n^2 is not a polynomial of degree smaller 2.

            ::

                sage: Seq2._parse_recursions_([f(3*n + 1) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 3 is not a power of 2.

            ::

                sage: Seq2._parse_recursions_([f(n + 1) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1 is less than 2.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(n), f(2*n) == 0], f, n)
                Traceback (most recent call last):
                ...
                ValueError: There are more than one recursions for f(2*n).

            ::

                sage: Seq2._parse_recursions_([f(2*n + 2) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 2 is not smaller than 2.

            ::

                sage: Seq2._parse_recursions_([f(2*n - 1) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: -1 is smaller than 0.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == 2*n], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 2*n does not contain f.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == 1/2*f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1/2 is not a valid coefficient.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == 1/f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1/f(n) is not a valid right hand side.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == 1/f(n) + 2*f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1/f(n) is not a valid summand.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(n + 1/2)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: n + 1/2 does not have integer coefficients.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(1/2*n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1/2*n does not have integer coefficients.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(n^2 + 1)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: n^2 + 1 does not have degree 1.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(1)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1 does not have degree 1.

            ::

                sage: Seq2._parse_recursions_([f(4*n) == f(2*n) + f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1 does not equal 2.

            ::

                sage: Seq2._parse_recursions_([f(4*n) == f(2*n), f(4*n + 1) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1 does not equal 2.

            ::

                sage: Seq2._parse_recursions_([f(4*n) == f(3*n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 3 is not a power of 2.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(4*n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 4 is not smaller than 2.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(2*n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 2 is not smaller than 2.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: Recursions for [f(2*n + 1)] are missing.

            ::

                sage: Seq2._parse_recursions_([f(4*n) == f(n), f(4*n + 3) == 0], f, n)
                Traceback (most recent call last):
                ...
                ValueError: Recursions for [f(4*n + 1), f(4*n + 2)] are missing.

            Finally, also for the zero-sequence the output is as expected::

                sage: Seq2._parse_recursions_([f(2*n) == 0, f(2*n + 1) == 0], f, n)
                recursion_rules(M=1, m=0, l=0, u=0, ll=0, uu=0, dim=1,
                coeffs={}, start_values={}, n0=0)
        """
        from collections import namedtuple

        from sage.arith.srange import srange
        from sage.functions.log import log
        from sage.functions.other import ceil, floor
        from sage.rings.integer_ring import ZZ
        from sage.symbolic.operators import add_vararg, mul_vararg, operator

        k = self.k
        base_ring = self.base()
        indices_right = [0]
        coeffs = {}
        start_values = {}
        remainders = []

        def _parse_multiplication_(op):
            if op.operator() != mul_vararg or len(op.operands()) != 2:
                raise ValueError("")
            operands = op.operands()
            if operands[1].operator() == function:
                return [operands[0], operands[1]]
            elif operands[0].operator() == function:
                return [operands[1], operands[0]]
            else:
                raise ValueError('%s does not contain %s.' % (op, function))

        def _parse_one_summand_(summand):
            if summand.operator() == mul_vararg:
                coeff, op = _parse_multiplication_(summand)
            elif summand.operator() == function:
                coeff, op = 1, summand
            else:
                raise ValueError('%s is not a valid summand.' % (summand,))
            try:
                poly = ZZ[var](op.operands()[0])
            except TypeError:
                raise ValueError('%s does not have integer coefficients.'
                                 % (op.operands()[0],))
            if poly.degree() != 1:
                raise ValueError("%s does not have degree 1."
                                 % (poly,))
            d, base_power_m = list(poly)
            m = log(base_power_m, base=k)
            return [coeff, m, d]

        if not equations:
            raise ValueError("List of recursion equations is empty.")

        for eq in equations:
            if eq.operator() != operator.eq:
                raise ValueError("%s is not an equation with ==."  % eq)
            left_side, right_side = eq.operands()
            if left_side.operator() != function:
                raise ValueError("%s is not an evaluation of %s."
                                 % (left_side, function))
            if  len(left_side.operands()) != 1:
                raise ValueError("%s does not have one argument." %
                                 (left_side,))
            try:
                polynomial_left = base_ring[var](left_side.operands()[0])
            except Exception:
                raise ValueError("%s is not a polynomial "
                                 "in %s." % (left_side.operands()[0], var))
            if polynomial_left.degree()  > 1:
                raise ValueError("%s is not a polynomial of degree smaller 2."
                                 % (polynomial_left,))
            if polynomial_left in base_ring and right_side in base_ring:
                start_values.update({polynomial_left: right_side})
            else:
                poly_left = ZZ[var](left_side.operands()[0])
                [r, base_power_M] = list(poly_left)
                M_new = log(base_power_M, base=k)
                try:
                    if M != log(base_power_M, base=k):
                        raise ValueError("%s does not equal %s."
                                         % (base_power_M, k^M))
                except NameError:
                    M = M_new
                    if M not in ZZ:
                        raise ValueError("%s is not a power of %s."
                                         % (base_power_M, k))
                    if M < 1:
                        raise ValueError("%s is less than %s."
                                         % (base_power_M, k))
                if r in remainders:
                    raise ValueError("There are more than one recursions for %s."
                                     % (left_side,))
                if r not in ZZ:
                    raise ValueError("%s is not an integer." % (r,))
                if r >= k**M:
                    raise ValueError("%s is not smaller than %s." % (r, k**M))
                if r < 0:
                    raise ValueError("%s is smaller than 0." % (r,))
                remainders.append(r)

                if right_side != 0:
                    if (len(right_side.operands()) == 1 and right_side.operator() == function
                        or right_side.operator() == mul_vararg and len(right_side.operands()) == 2):
                        summands = [right_side]
                    elif right_side.operator() == add_vararg:
                        summands = right_side.operands()
                    else: # check this again
                        raise ValueError("%s is not a valid right hand side."
                                         % (right_side,))
                    for summand in summands:
                        coeff, new_m, d = _parse_one_summand_(summand)
                        if coeff not in base_ring:
                            raise ValueError("%s is not a valid coefficient."
                                             % (coeff,))
                        try:
                            if m != new_m:
                                raise ValueError("%s does not equal %s."
                                                 % (k**new_m, k**m))
                        except NameError:
                            m = new_m
                            if m not in ZZ:
                                raise ValueError("%s is not a power of %s."
                                                 % (k**m, k))
                            if M <= m:
                                raise ValueError("%s is not smaller than %s."
                                                 % (k**m, k**M))

                        indices_right.append(d)
                        coeffs.update({(r, d): coeff})

        remainders.sort()
        if remainders != srange(k**M):
            missing_equations = [function(k**M*var + r)
                                 for r in srange(k**M)
                                 if r not in remainders]
            raise ValueError("Recursions for %s are missing."
                             % missing_equations)

        if not coeffs:
            m = M - 1

        l = min(indices_right)
        u = max(indices_right)
        ll = (floor((l*k**(M-m) - k**M + 1)/(k**(M-m) - 1)) + 1)*(l < 0)
        uu = max([ceil((u*k**(M-m) + k**M - k**m)/(k**(M-m) - 1)) - 1, k**m - 1])
        dim = (k**M - 1)/(k - 1) + (M - m)*(uu - ll - k**m + 1)

        recursion_rules = namedtuple('recursion_rules',
                                     ['M', 'm', 'l', 'u',
                                      'll', 'uu', 'dim',
                                      'coeffs', 'start_values', 'n0'])

        return recursion_rules(M=M, m=m, l=l, u=u, ll=ll, uu=uu, dim=dim,
                               coeffs=coeffs, start_values=start_values, n0=n0)


    def _get_matrix_from_recursions_(self, recursion_rules, rem, function, var):
        r"""
        Construct the matrix for remainder ``rem`` of the linear
        representation of the sequence induced by ``recursion_rules``.

        INPUT:

        - ``recursion_rules`` -- A namedtuple generated by
          :meth:`~_parse_recursions_`.

        - ``rem`` -- An integer between ``0`` and ``k - 1``.

        OUTPUT:

        A matrix.

        EXAMPLES:

        The following example illustrates how the coefficients in the
        right hand sides of the recursions correspond to the entries of
        the matrices. ::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: rules = Seq2._parse_recursions_([
            ....:     f(8*n) == 1*f(2*n) + 2*f(2*n + 1) + 3*f(2*n + 2),
            ....:     f(8*n + 1) == 4*f(2*n) + 5*f(2*n + 1) + 6*f(2*n + 2),
            ....:     f(8*n + 2) == 7*f(2*n) + 8*f(2*n + 1) + 9*f(2*n + 2),
            ....:     f(8*n + 3) == 10*f(2*n) + 11*f(2*n + 1) + 12*f(2*n + 2),
            ....:     f(8*n + 4) == 13*f(2*n) + 14*f(2*n + 1) + 15*f(2*n + 2),
            ....:     f(8*n + 5) == 16*f(2*n) + 17*f(2*n + 1) + 18*f(2*n + 2),
            ....:     f(8*n + 6) == 19*f(2*n) + 20*f(2*n + 1) + 21*f(2*n + 2),
            ....:     f(8*n + 7) == 22*f(2*n) + 23*f(2*n + 1) + 24*f(2*n + 2),],
            ....:     f, n)
            sage: Seq2._get_matrix_from_recursions_(rules, 0, f, n)
            [ 0  1  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  1  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  1  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  1  0  0]
            [ 0  1  2  3  0  0  0  0  0  0  0  0  0]
            [ 0  4  5  6  0  0  0  0  0  0  0  0  0]
            [ 0  7  8  9  0  0  0  0  0  0  0  0  0]
            [ 0 10 11 12  0  0  0  0  0  0  0  0  0]
            [ 0 13 14 15  0  0  0  0  0  0  0  0  0]
            [ 0 16 17 18  0  0  0  0  0  0  0  0  0]
            [ 0 19 20 21  0  0  0  0  0  0  0  0  0]
            sage: Seq2._get_matrix_from_recursions_(rules, 1, f, n)
            [ 0  0  1  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  1  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  1  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  1  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  1]
            [ 0 13 14 15  0  0  0  0  0  0  0  0  0]
            [ 0 16 17 18  0  0  0  0  0  0  0  0  0]
            [ 0 19 20 21  0  0  0  0  0  0  0  0  0]
            [ 0 22 23 24  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  1  2  3  0  0  0  0  0  0  0]
            [ 0  0  0  4  5  6  0  0  0  0  0  0  0]
            [ 0  0  0  7  8  9  0  0  0  0  0  0  0]

        Stern--Brocot Sequence::

            sage: SB_rules = Seq2._parse_recursions_([
            ....:     f(2*n) == f(n), f(2*n + 1) == f(n) + f(n + 1),
            ....:     f(0) == 0, f(1) == 1], f, n)
            sage: Seq2._get_matrix_from_recursions_(SB_rules, 0, f, n)
            [1 0 0]
            [1 1 0]
            [0 1 0]
            sage: Seq2._get_matrix_from_recursions_(SB_rules, 1, f, n)
            [1 1 0]
            [0 1 0]
            [0 1 1]

        """
        from sage.arith.srange import srange
        from sage.matrix.constructor import Matrix
        from sage.matrix.matrix_space import MatrixSpace
        from sage.matrix.special import block_matrix, zero_matrix
        from sage.modules.free_module_element import vector

        k = self.k
        M = recursion_rules.M
        m = recursion_rules.m
        l = recursion_rules.l
        u = recursion_rules.u
        ll = recursion_rules.ll
        uu = recursion_rules.uu
        n_start = recursion_rules.n_start
        dim = recursion_rules.dim
        coeffs = recursion_rules.coeffs
        start_values = recursion_rules.start_values
        n0 = recursion_rules.n0

        mat = []
        current_shift = 0

        for base_power in srange(m - 1):
            current_shift += k**base_power
            for d in srange(k**base_power):
                row = dim*[0]
                dd = k**base_power*rem + d
                index = current_shift + dd
                row[index] = 1
                mat.append(row)

        if m > 0:
            current_shift += k**(m-1)
            final_shift = current_shift - ll
            for d in srange(k**(m-1)):
                row = dim*[0]
                dd = k**(m-1)*rem + d
                index = current_shift + dd - ll
                row[index] = 1
                mat.append(row)
        else:
            final_shift = -ll

        for base_power in srange(m, M - 1):
            current_shift += k**base_power - k**m + uu - 2*ll + 1
            for d in srange(ll, k**base_power - k**m + uu + 1):
                row = dim*[0]
                dd = k**base_power*rem + d
                index = current_shift + dd
                row[index] = 1
                mat.append(row)

        for d in srange(ll, k**(M-1) - k**m + uu + 1):
            dd, r = d.quo_rem(k**M)
            rr = k**(M-1)*rem + r
            row = dim*[0]
            if rr < k**M:
                for i in range(l, u + 1):
                    try:
                        row[k**m*dd + i + final_shift] = coeffs[(rr, i)]
                    except KeyError:
                        pass
            else:
                for i in srange(l, u + 1):
                    try:
                        row[(dd + 1)*k**m + i + final_shift] = coeffs[(rr - k**M, i)]
                    except KeyError:
                        pass
            mat.append(row)

        mat = Matrix(mat)

        if n0 == 0:
            return mat
        else:
            arguments = [k**j*var + d for j in srange(m) for d in srange(k**j)] + \
                        [k**j*var + d for j in srange(m, M) for d in srange(ll, k**j - k**m + uu + 1)]
            W = []
            for i in srange(n0):
                v_eval_i = []
                v_eval_ki_plus_r = []
                for a in arguments:
                    try:
                        temp = a.substitute(var==i)
                        v_eval_i.append(start_values[temp])
                        temp = a.substitute(var==k*i+rem)
                        v_eval_ki_plus_r.append(start_values[temp])
                    except KeyError:
                        raise ValueError('Start value %s is missing.'
                                         % (function(temp),))
                W.append(list(vector(v_eval_ki_plus_r) - mat*vector(v_eval_i)))

            J = []
            for i in srange(n0):
                J.append([int(i >= rem and i % k == rem and j*k == i - rem) for j in srange(n0)])

            Mat = MatrixSpace(self.base_ring(), dim + n0, dim + n0)
            return Mat(block_matrix([[mat, Matrix(W).transpose()],
                                     [zero_matrix(n0, dim), Matrix(J)]]))


    def _get_left_from_recursions_(self, dim):
        r"""
        """
        from sage.modules.free_module_element import vector

        return vector([1] + (dim - 1)*[0])


    def _get_right_from_recursions_(self, recursion_rules, function):
        r"""
        """
        from sage.arith.srange import srange
        from sage.modules.free_module_element import vector

        base = self.k
        M = recursion_rules.M
        m = recursion_rules.m
        ll = recursion_rules.ll
        uu = recursion_rules.uu
        start_values = recursion_rules.start_values
        n0 = recursion_rules.n0
        right = []

        for j in srange(m):
            for d in srange(base**j):
                try:
                    right.append(start_values[d])
                except KeyError:
                    raise ValueError('Start value %s is missing.'
                                     % (function(d),))

        for j in srange(m, M):
            for d in srange(ll, base**j - base**m + uu + 1):
                try:
                    right.append(start_values[d])
                except KeyError:
                    raise ValueError('Start value %s is missing.'
                                     % (function(d),))

        if n0 >= 1:
            right = right + [1] + (n0 - 1)*[0]

        return vector(right)


    def recursions(self, equations, function, var, n0=0, minimize=False):
        r"""
        Construct a `k`-regular sequence that fulfills the recursions
        given in ``equations``.

        INPUT:

        - ``equations`` -- A list of equations where the elements have
          either the form

          - ``f(k^M * n + r) == sum(f(k^m * n + k) for k in srange(l, u + 1)``
            for some integers ``0 <= r < k^M`` and ``M > m >= 0``
            and some ``l <= u`` -- valid for all non-negative
            integers ``n`` -- and there is an equation of this form
            for all ``r``

          or the form

          - ``f(k) == t`` for some integer ``k`` and some ``t``.

        - ``function`` -- symbolic function ``f`` occuring in the equations.

        - ``var`` -- symbolic variable ``n`` occuring in the equations.

        - ``minimize`` -- (default: ``False``) a boolean. If ``True``,
          the linear representation of the resulting sequence is
          minimized (see :meth:`~minimized`).

        OUTPUT:

        A :class:`kRegularSequence`.

        EXAMPLES:

        Stern--Brocot Sequence::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: Seq2.recursions([
            ....:     f(2*n) == f(n), f(2*n + 1) == f(n) + f(n + 1),
            ....:     f(0) == 0, f(1) == 1, f(2) == 1], f, n)
            2-regular sequence 0, 1, 1, 2, 1, 3, 2, 3, 1, 4, ...

        Number of Odd Entries in Pascal's Triangle::

            sage: Seq2.recursions([
            ....:     f(2*n) == 3*f(n), f(2*n + 1) == 2*f(n) + f(n + 1),
            ....:     f(0) == 0, f(1) == 1, f(2) == 3, f(3) == 5], f, n)
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...

        Number of Unbordered Factors in the Thue--Morse Sequence::

            sage: Seq2.recursions([
            ....:     f(8*n) == 2*f(4*n), f(8*n + 1) == f(4*n + 1),
            ....:     f(8*n + 2) == f(4*n + 1) + f(4*n + 3),
            ....:     f(8*n + 3) == -f(4*n + 1) + f(4*n + 2),
            ....:     f(8*n + 4) == 2*f(4*n + 2), f(8*n + 5) == f(4*n + 3),
            ....:     f(8*n + 6) == -f(4*n + 1) + f(4*n + 2) + f(4*n + 3),
            ....:     f(8*n + 7) == 2*f(4*n + 1) + f(4*n + 3),
            ....:     f(0) == 1, f(1) == 2, f(2) == 2, f(3) == 4, f(4) == 2,
            ....:     f(5) == 4, f(6) == 6, f(7) == 0, f(8) == 4, f(9) == 4,
            ....:     f(10) == 4, f(11) == 4, f(12) == 12, f(13) == 0, f(14)==4,
            ....:     f(15) == 4, f(16) == 8, f(17) == 4, f(18) == 8, f(19) == 0,
            ....:     f(20) == 8, f(21) == 4, f(22) == 4, f(23) == 8, f(24) == 24,
            ....:     f(25) == 0, f(26) == 4, f(27) == 4, f(28) == 8, f(29) == 4,
            ....:     f(30) == 8, f(31) == 4, f(32) == 16, f(33) == 4], f, n, 3)
            2-regular sequence 1, 2, 2, 4, 2, 4, 6, 0, 4, 4, ...

        Number of Non-Zero Elements in the a Generalized Pascal's Triangle [TODO: reference]::

            sage: Seq2 = kRegularSequenceSpace(2, QQ)
            sage: Seq2.recursions([
            ....:     f(4*n) == 5/3*f(2*n) - 1/3*f(2*n + 1),
            ....:     f(4*n + 1) == 4/3*f(2*n) + 1/3*f(2*n + 1),
            ....:     f(4*n + 2) == 1/3*f(2*n) + 4/3*f(2*n + 1),
            ....:     f(4*n + 3) == -1/3*f(2*n) + 5/3*f(2*n + 1),
            ....:     f(0) == 1, f(1) == 2, f(2) == 3, f(3) == 3], f, n)
            2-regular sequence 1, 2, 3, 3, 4, 5, 5, 4, 5, 7, ...
        """
        from sage.arith.srange import srange

        k = self.k
        mu = []

        recursion_rules = self._parse_recursions_(equations, function, var, n0)

        for rem in srange(k):
            mu.append(self._get_matrix_from_recursions_(recursion_rules, rem, function, var))

        seq = self(mu, self._get_left_from_recursions_(recursion_rules.dim + n0),
                   self._get_right_from_recursions_(recursion_rules, function))

        if minimize:
            return seq.minimized()
        else:
            return seq
