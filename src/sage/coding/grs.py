r"""
Generalized Reed-Solomon Code

Given `n` different evaluation points `\alpha_1, \dots, \alpha_n` from some
finite field `F`, and `n` column multipliers `\beta_1, \dots, \beta_n`, the
corresponding GRS code of dimension `k` is the set:

.. math::

    \{ (\beta_1 f(\alpha_1), \ldots, \beta_n f(\alpha_n)  \mid  f \in F[x], \deg f < k \}

This file contains the following elements:

    - :class:`GeneralizedReedSolomonCode`, the class for GRS codes
    - :class:`GRSEvaluationVectorEncoder`, an encoder with a vectorial message space
    - :class:`GRSEvaluationPolynomialEncoder`, an encoder with a polynomial message space
"""

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import matrix, diagonal_matrix
from sage.modules.free_module_element import vector
from sage.misc.cachefunc import cached_method
from copy import copy
from linear_code import AbstractLinearCode
from encoder import Encoder
from sage.misc.misc_c import prod

class GeneralizedReedSolomonCode(AbstractLinearCode):
    r"""
    Representation of a Generalized Reed-Solomon code.

    INPUT:

    - ``evaluation_points`` -- A list of evaluation points in a finite field F

    - ``dimension`` -- The dimension of the code

    - ``column_multipliers`` -- (default: ``None``) List of column multipliers in F for this code.
      All column multipliers are set to 1 if default value is kept.

    EXAMPLES:

    We construct a GRS code with a manually built support, without specifying column multipliers::

        sage: F = GF(7)
        sage: support = [F(i) for i in range(1,7)]
        sage: C = codes.GeneralizedReedSolomonCode(support,3)
        sage: C
        [6, 3, 4] Generalized Reed-Solomon Code over Finite Field of size 7


    We construct a GRS code without specifying column multipliers::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: C
        [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59

    It is also possible to specify the column multipliers::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, F.list()[1:n+1])
        sage: C
        [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, evaluation_points, dimension, column_multipliers=None):
        r"""
        TESTS:

        If the evaluation points are not from a finite field, it raises an error::

            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(list()[:n], k)
            Traceback (most recent call last):
            ...
            ValueError: Evaluation points must be in a finite field

        If the column multipliers are not from a finite field, or not in the same
        finite field as the evaluation points, it raises an error::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, list()[1:n+1])
            Traceback (most recent call last):
            ...
            ValueError: Column multipliers must be in a finite field

            sage: F2 = GF(61)
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, F2.list()[1:n+1])
            Traceback (most recent call last):
            ...
            ValueError: Column multipliers and evaluation points must be in the same field

        The number of column multipliers is checked as well::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, F.list()[1:n])
            Traceback (most recent call last):
            ...
            ValueError: There must be exactly 40 column multipliers

        It is not allowed to have 0 as a column multiplier::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, F.list()[:n])
            Traceback (most recent call last):
            ...
            ValueError: All column multipliers must be non-zero

        And all the evaluation points must be different::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode([F.one()]*n, k)
            Traceback (most recent call last):
            ...
            ValueError: All evaluation points must be different
        """
        F = vector(evaluation_points).base_ring()
        if F.is_finite() == False:
            raise ValueError("Evaluation points must be in a finite field")
        super(GeneralizedReedSolomonCode, self).__init__(F, len(evaluation_points),\
                "EvaluationVector", "Syndrome")
        self._dimension = dimension
        self._evaluation_points = copy(evaluation_points)

        if column_multipliers is None:
            self._column_multipliers = [self.base_field().one()] * self._length
        else:
            Fc = vector(column_multipliers).base_ring()
            if Fc.is_finite() == False:
                raise ValueError("Column multipliers must be in a finite field")
            elif Fc != self.base_field():
                raise ValueError("Column multipliers and evaluation points\
                        must be in the same field")
            self._column_multipliers = copy(column_multipliers)
        if len(self._column_multipliers) != self._length:
            raise ValueError("There must be exactly %s column multipliers"\
                    % self._length)
        if 0 in self._column_multipliers:
            raise ValueError("All column multipliers must be non-zero")
        if len(self._evaluation_points) != len(set(self._evaluation_points)):
            raise ValueError("All evaluation points must be different")

    def __eq__(self, other):
        r"""
        Tests equality between Generalized Reed-Solomon Code objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C1 = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C2 = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C1.__eq__(C2)
            True
        """
        return isinstance(other, GeneralizedReedSolomonCode) \
                and self.base_field() == other.base_field() \
                and self.length() == other.length() \
                and self.dimension() == other.dimension() \
                and self.evaluation_points() == other.evaluation_points() \
                and self.column_multipliers() == other.column_multipliers()

    def __ne__(self, other):
        r"""
        Tests inequality of Generalized Reed-Solomon Code objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C1 = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C2 = codes.GeneralizedReedSolomonCode(F.list()[:n], k + 1)
            sage: C1.__ne__(C2)
            True
        """
        return not self.__eq__(other)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C
            [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "[%s, %s, %s] Generalized Reed-Solomon Code over %s"\
                % (self.length(), self.dimension(),\
                self.minimum_distance(), self.base_field())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: latex(C)
            [40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "[%s, %s, %s] \\textnormal{ Generalized Reed-Solomon Code over } %s"\
                % (self.length(), self.dimension() ,self.minimum_distance(),\
                self.base_field()._latex_())

    def minimum_distance(self):
        r"""
        Returns the minimum distance of ``self``. Since a GRS code is always MDS,
        this always returns ``C.length() - C.dimension() + 1``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.minimum_distance()
            29
        """
        return self.length() - self.dimension() + 1

    def evaluation_points(self):
        r"""
        Returns the list of evaluation points of ``self``.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.evaluation_points()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        return self._evaluation_points

    def column_multipliers(self):
        r"""
        Returns the list of column multipliers of ``self``.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.column_multipliers()
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        return self._column_multipliers

    @cached_method
    def multipliers_product(self):
        r"""
        Returns the list of products of the j-th column multiplier of ``self`` with
        the j-th column multiplier of ``self``'s parity check matrix, for j between 0 and
        the length of ``self``.

        AUTHORS:

        This function is taken from codinglib (https://bitbucket.org/jsrn/codinglib/)
        and was written by Johan Nielsen.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.multipliers_product()
            [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        """
        a = self.evaluation_points()
        return [ 1/prod([ a[i] - a[h] for h in range(0, len(a)) if h != i ])
                    for i in range(0,len(a)) ]

    @cached_method
    def parity_column_multipliers(self):
        r"""
        Returns the list of column multipliers of ``self``'s parity check matrix.

        AUTHORS:

        This function is taken from codinglib (https://bitbucket.org/jsrn/codinglib/)
        and was written by Johan Nielsen.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.parity_column_multipliers()
            [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        """
        n = self.length()
        col_mults = self.column_multipliers()
        etas = self.multipliers_product()
        return [ etas[i]/col_mults[i] for i in range(0,n) ]

    @cached_method
    def parity_check_matrix(self):
        r"""
        Returns the parity check matrix of ``self``.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.parity_check_matrix()
            [10  9  8  7  6  5  4  3  2  1]
            [ 0  9  5 10  2  3  2 10  5  9]
            [ 0  9 10  8  8  4  1  4  7  4]
            [ 0  9  9  2 10  9  6  6  1  3]
            [ 0  9  7  6  7  1  3  9  8  5]
        """
        F = self.base_ring()
        n = self.length()
        d = self.minimum_distance()
        alphas = self.evaluation_points()
        return matrix(F, d-1, n, lambda i,j : alphas[j] ** i) *\
                diagonal_matrix(F, self.parity_column_multipliers())









class GRSEvaluationVectorEncoder(Encoder):
    r"""
    An encoder which can encode vectors into codewords.

    INPUT:

    - ``code`` -- The associated code of this encoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
        sage: E
        Evaluation vector-style encoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: E
            Evaluation vector-style encoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        super(GRSEvaluationVectorEncoder, self).__init__(code)
        self._R = code.base_field()['x']

    def __eq__(self, other):
        r"""
        Tests equality between GRSEvaluationVectorEncoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D1 = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: D2 = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: D1.__eq__(D2)
            True
        """
        return isinstance(other, GRSEvaluationVectorEncoder) \
                and self.code() == other.code()

    def __ne__(self, other):
        r"""
        Tests difference between GRSEvaluationVectorEncoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C1 = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C2 = codes.GeneralizedReedSolomonCode(F.list()[:n], k + 1)
            sage: D1 = codes.encoders.GRSEvaluationVectorEncoder(C1)
            sage: D2 = codes.encoders.GRSEvaluationVectorEncoder(C2)
            sage: D1.__ne__(D2)
            True
        """
        return not self.__eq__(other)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: E
            Evaluation vector-style encoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Evaluation vector-style encoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: latex(E)
            \textnormal{Evaluation vector-style encoder for the }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Evaluation vector-style encoder for the }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of ``self``

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: E.generator_matrix()
            [1 1 1 1 1 1 1 1 1 1]
            [0 1 2 3 4 5 6 7 8 9]
            [0 1 4 9 5 3 3 5 9 4]
            [0 1 8 5 9 4 7 2 6 3]
            [0 1 5 4 3 9 9 3 4 5]
        """
        base_field = self.code().base_field()
        dimension = self.code().dimension()
        length = self.code().length()
        alphas = self.code().evaluation_points()
        col_mults = self.code().column_multipliers()
        return matrix(base_field, dimension, length, lambda i,j : col_mults[j]*alphas[j]**i)

    def unencode_nocheck(self, c):
        r"""
        Returns the message corresponding to ``c``.
        Does not check if ``c`` belongs to the code.

        INPUT:

        - ``c`` -- A vector with the same length as the code

        OUTPUT:

        - An element of the message space

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10 , 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: c = vector(F, (10, 3, 9, 6, 5, 6, 9, 3, 10, 8))
            sage: E.unencode_nocheck(c)
            (10, 3, 1, 0, 0)
        """
        C = self.code()
        alphas = self.code().evaluation_points()
        col_mults = self.code().column_multipliers()
        length = self.code().length()
        dimension = self.code().dimension()

        c = [c[i]/col_mults[i] for i in range(length)]
        points = [(alphas[i], c[i]) for i in range(dimension)]

        Pc = self._R.lagrange_polynomial(points).list()
        Pc = Pc + [self.code().base_field().zero()]*(dimension - len(Pc))

        m = vector(self.code().base_field(), Pc)
        return m










class GRSEvaluationPolynomialEncoder(Encoder):
    r"""
    An encoder which can encode polynomials into codewords.

    INPUT:

    - ``code`` -- The associated code of this encoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: E = codes.encoders.GRSEvaluationPolynomialEncoder(C)
        sage: E
        Evaluation polynomial-style encoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: E
            Evaluation polynomial-style encoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        super(GRSEvaluationPolynomialEncoder, self).__init__(code)
        self._R = code.base_field()['x']

    def __eq__(self, other):
        r"""
        Tests equality between GRSEvaluationPolynomialEncoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D1 = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: D2 = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: D1.__eq__(D2)
            True
        """
        return isinstance(other, GRSEvaluationPolynomialEncoder) \
                and self.code() == other.code()

    def __ne__(self, other):
        r"""
        Tests difference between GRSEvaluationPolynomialEncoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C1 = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C2 = codes.GeneralizedReedSolomonCode(F.list()[:n], k + 1)
            sage: D1 = codes.encoders.GRSEvaluationPolynomialEncoder(C1)
            sage: D2 = codes.encoders.GRSEvaluationPolynomialEncoder(C2)
            sage: D1.__ne__(D2)
            True
        """
        return not self.__eq__(other)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: E
            Evaluation polynomial-style encoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Evaluation polynomial-style encoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: latex(E)
            \textnormal{Evaluation polynomial-style encoder for the }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Evaluation polynomial-style encoder for the }%s" % self.code()._latex_()

    def encode(self, p):
        r"""
        Transforms ``p`` into an element of the associated code of ``self``.

        INPUT:

        - ``p`` -- A polynomial from ``self`` message space

        OUTPUT:

        - A codeword in associated code of ``self``

        EXAMPLES::

            sage: F = GF(11)
            sage: K.<x>=F[]
            sage: n, k = 10 , 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: p = x^2 + 3*x + 10
            sage: E.encode(p)
            (10, 3, 9, 6, 5, 6, 9, 3, 10, 8)
        """
        alphas    = self.code().evaluation_points()
        col_mults = self.code().column_multipliers()
        field     = self.code().base_ring()
        length    = self.code().length()
        u = vector(field, [col_mults[i]*p(alphas[i]) for i in range(length)])
        return u

    def unencode_nocheck(self, c):
        r"""
        Returns the message corresponding to ``c``.
        Does not check if ``c`` belongs to the code.

        INPUT:

        - ``c`` -- A vector with the same length as the code

        OUTPUT:

        - An element of the message space

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10 , 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: c = vector(F, (10, 3, 9, 6, 5, 6, 9, 3, 10, 8))
            sage: E.unencode_nocheck(c)
            x^2 + 3*x + 10
        """

        alphas = self.code().evaluation_points()
        col_mults = self.code().column_multipliers()
        length = self.code().length()
        dimension = self.code().dimension()

        c = [c[i]/col_mults[i] for i in range(length)]
        points = [(alphas[i], c[i]) for i in range(dimension)]

        Pc = self._R.lagrange_polynomial(points)
        return Pc

    def message_space(self):
        r"""
        Returns the message space of ``self``

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10 , 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: E.message_space()
            Univariate Polynomial Ring in x over Finite Field of size 11
        """
        return self._R
