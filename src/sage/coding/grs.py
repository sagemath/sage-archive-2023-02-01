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
    - :class:`GRSBerlekampWelchDecoder`, a decoder based on the Berlekamp-Welch decoding algorithm
    - :class:`GRSGaoDecoder`, a decoder based on the Gao decoding algorithm
    - :class:`GRSErrorErasureDecoder`, a decoder able to correct both errors and erasures
    - :class:`GRSKeyEquationSyndromeDecoder`, a decoder using key equation decoding
      based on syndrome polynomial

REFERENCES:

    .. [N13] Johan S. R. Nielsen, List Decoding of Algebraic Codes, 2013
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
from sage.rings.finite_rings.constructor import GF
from sage.categories.cartesian_product import cartesian_product
from sage.modules.free_module_element import vector
from sage.modules.free_module import VectorSpace
from sage.rings.integer import Integer
from sage.misc.cachefunc import cached_method
from copy import copy
from linear_code import (AbstractLinearCode,
                         LinearCodeSyndromeDecoder,
                         LinearCodeNearestNeighborDecoder)
from encoder import Encoder
from decoder import Decoder, DecodingError
from sage.rings.arith import xgcd
from sage.misc.misc_c import prod
from sage.functions.other import binomial, floor, sqrt
from sage.calculus.var import var
from sage.misc.functional import symbolic_sum
from sage.rings.integer_ring import ZZ
from gs_interpolation import *#gs_construct_Q_linalg as gs_construct_Q_linalg
from gs_rootfinding import *

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
        super(GeneralizedReedSolomonCode, self).__init__(F, len(evaluation_points), "EvaluationVector", "Syndrome")
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

            This function is taken from codinglib [Nielsen]_

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

            This function is taken from codinglib [Nielsen]_

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
        return [ etas[i]/col_mults[i] for i in range(n) ]

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

    def dual_code(self):
        r"""
        Returns the dual code of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: Cd = C.dual_code()

        Dual code of the dual code is the original code::

            sage: C == Cd.dual_code()
            True
        """
        col_mults = self.parity_column_multipliers()
        return GeneralizedReedSolomonCode(self.evaluation_points(), self.length() - self.dimension(), col_mults)

    def covering_radius(self):
        r"""
        Returns the covering radius of ``self``.

        The covering radius of a linear code `C` is the smallest
        number `r` s.t. any element of
        the ambient space of `C` is at most at distance `r` to `C`.

        As Reed-Solomon codes are MDS codes, their covering radius
        is always `d-1`, where `d` is the minimum distance.

        This is a custom method for GRS codes which is faster than
        generic implementation.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.covering_radius()
            5
        """
        return self.length() - self.dimension()

    @cached_method
    def weight_distribution(self):
        r"""
        Returns the weight distribution of ``self``.

        The weight distribution is returned as a list, where the
        `i-th` entry corresponds to the number of words of weight `i` in the code.

        This is a custom method for GRS codes which is faster than
        generic implementation.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.weight_distribution()
            [1, 0, 0, 0, 0, 0, 2100, 6000, 29250, 61500, 62200]
        """
        d = self.minimum_distance()
        n = self.length()
        q = self.base_ring().order()
        s = var('s')
        wd = [1] + [0] * (d - 1)
        for i in range(d, n+1):
            tmp = binomial(n, i) * (q - 1)
            wd.append(tmp * symbolic_sum(binomial(i-1, s) * (-1) ** s * q ** (i - d - s), s, 0, i-d))
        return wd

    def weight_enumerator(self):
        r"""
        Returns the weight enumerator of ``self``.

        This is a custom method for GRS codes which is faster than
        generic implementation.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.weight_enumerator()
            62200*x^10 + 61500*x^9 + 29250*x^8 + 6000*x^7 + 2100*x^6 + 1
        """
        PolRing = ZZ['x']
        x = PolRing.gen()
        s = var('s')
        wd = self.weight_distribution()
        d = self.minimum_distance()
        n = self.length()
        w_en = PolRing(1)
        for i in range(n + 1 - d):
            w_en += wd[i + d] * x ** (i + d)
        return w_en








####################### encoders ###############################

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









####################### decoders ###############################

class GRSBerlekampWelchDecoder(Decoder):
    r"""
    A decoder based on the Berlekamp-Welch decoding algorithm.

    This algorithm recovers the error locator polynomial  by solving a linear system.
    See [HJ04-BW] for details.

    REFERENCES:

    .. [HJ04-BW] Tom Hoiholdt and Joirn Justesen, A Course In Error-Correcting Codes,
       EMS, 2004, pp 51-52

    INPUT:

    - ``code`` -- A code associated to this decoder

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
        sage: D
        Berlekamp-Welch decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: D
            Berlekamp-Welch decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        super(GRSBerlekampWelchDecoder, self).__init__(code, code.ambient_space(),\
            "EvaluationPolynomial")

    def __eq__(self, other):
        r"""
        Tests equality between GRSBerlekampWelchDecoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D1 = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: D2 = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: D1.__eq__(D2)
            True
        """
        return isinstance(other, GRSBerlekampWelchDecoder) \
                and self.code() == other.code()\
                and self.input_space() == other.input_space()

    def __ne__(self, other):
        r"""
        Tests difference between GRSBerlekampWelchDecoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C1 = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C2 = codes.GeneralizedReedSolomonCode(F.list()[:n], k + 1)
            sage: D1 = codes.decoders.GRSBerlekampWelchDecoder(C1)
            sage: D2 = codes.decoders.GRSBerlekampWelchDecoder(C2)
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
            sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: D
            Berlekamp-Welch decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Berlekamp-Welch decoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: latex(D)
            \textnormal{Berlekamp Welch decoder for the }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Berlekamp Welch decoder for the }%s"\
                % self.code()._latex_()

    def decode_to_message(self, r):
        r"""
        Decodes ``r`` to an element in message space of ``self``.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self`` message space

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: r = vector(F, (8, 2, 6, 10, 6, 10, 7, 6, 7, 1))
            sage: D.decode_to_message(r)
            x^4 + 7*x^3 + 10*x^2 + 9*x + 8

        If we try to decode a word with too many errors, it returns
        an exception::

            sage: r[0] = r[1] = r[2] = 3
            sage: D.decode_to_message(r)
            Traceback (most recent call last):
            ...
            DecodingError
        """
        C = self.code()
        col_mults = C.column_multipliers()
        if r in C:
            return self.connected_encoder().unencode_nocheck(r)

        r = [r[i]/col_mults[i] for i in range(0, C.length())]

        t  = (C.minimum_distance()-1) // 2
        l0 = C.length()-1-t
        l1 = C.length()-1-t-(C.dimension()-1)
        S  = matrix(C.base_field(), C.length(), l0+l1+2, lambda i,j :\
                (C.evaluation_points()[i])**j if j<(l0+1)\
                else r[i]*(C.evaluation_points()[i])**(j-(l0+1)))
        S  = S.right_kernel()
        S  = S.basis_matrix().row(0)
        R = C.base_field()['x']

        Q0 = R(S.list_from_positions(xrange(0, l0+1)))
        Q1 = R(S.list_from_positions(xrange(l0+1 , l0+l1+2)))

        if not Q1.divides(Q0):
            raise DecodingError()
        f = (-Q0)//Q1

        return f

    def decoding_radius(self):
        r"""
        Returns maximal number of errors that ``self`` can decode.

        OUTPUT:

        - the number of errors as an integer

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSBerlekampWelchDecoder(C)
            sage: D.decoding_radius()
            14
        """
        return (self.code().minimum_distance()-1)//2











class GRSGaoDecoder(Decoder):
    r"""
    A decoder based on the Gao decoding algorithm.

    Gao decoding algorithm uses early terminated extended Euclidean algorithm
    to find the error locator polynomial. See [G02] for details.

    REFERENCES:

    .. [G02] Shuhong Gao, A new algorithm for decoding Reed-Solomon Codes, January 31, 2002

    INPUT:

    - ``code`` -- The associated code of this decoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: D = codes.decoders.GRSGaoDecoder(C)
        sage: D
        Gao decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: D
            Gao decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        super(GRSGaoDecoder, self).__init__(code, code.ambient_space(),\
                "EvaluationPolynomial")

    def __eq__(self, other):
        r"""
        Tests equality of GRSGaoDecoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D1 = codes.decoders.GRSGaoDecoder(C)
            sage: D2 = codes.decoders.GRSGaoDecoder(C)
            sage: D1.__eq__(D2)
            True
        """
        return isinstance(other, GRSGaoDecoder) \
                and self.code() == other.code()\
                and self.input_space() == other.input_space()

    def __ne__(self, other):
       r"""
       Tests difference of GRSGaoDecoder objects.

       EXAMPLES::

           sage: F = GF(59)
           sage: n, k = 40, 12
           sage: C1 = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
           sage: C2 = codes.GeneralizedReedSolomonCode(F.list()[:n], k + 1)
           sage: D1 = codes.decoders.GRSGaoDecoder(C1)
           sage: D2 = codes.decoders.GRSGaoDecoder(C2)
           sage: D1.__ne__(D2)
           True
       """
       return not self.__eq__(other)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: D
            Gao decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Gao decoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: latex(D)
            \textnormal{Gao decoder for the }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Gao decoder for the }%s" % self.code()._latex_()

    @cached_method
    def precompute(self, PolRing):
        r"""
        Return the unique monic polynomial vanishing on the evaluation points.
        Helper function for internal purposes.

        INPUT:

        - ``PolRing`` -- polynomial ring of the output

        OUTPUT:

        - a polynomial over ``PolRing``

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: P = PolynomialRing(F,'x')
            sage: D.precompute(P)
            x^10 + 10*x^9 + x^8 + 10*x^7 + x^6 + 10*x^5 + x^4 + 10*x^3 + x^2 + 10*x
        """
        alphas = self.code().evaluation_points()
        G = 1
        x = PolRing.gens()[0]
        for i in range(0, self.code().length()):
            G = G*(x-self.code().evaluation_points()[i])
        return G

    def partial_xgcd(self, a, b, PolRing):
        r"""
        Returns the greatest common divisor of ``a`` and ``b``.

        The computation stops whenever the degree of the xgcd falls below
        `\frac{d+k}{2}`, where `d` is the dimension of ``self.code()`` and
        `k` its dimension.

        This is a helper function, used in :meth:`decode_to_message`.

        INPUT:

        - ``a, b`` -- polynomials over ``PolRing``

        - ``PolRing`` -- polynomial ring of the output

        OUTPUT:

        - a tuple of polynomials ``(r, s)`` where ``r`` is to the xgcd of
          ``a`` and ``b`` and ``s`` is the Bezout coefficient of ``a``.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: P = PolynomialRing(F,'x')
            sage: x = P.parameter()
            sage: a = 5*x^2 + 9*x + 8
            sage: b = 10*x^2 + 3*x + 5
            sage: D.partial_xgcd(a, b, P)
            (10*x^2 + 3*x + 5, 1)
        """
        stop = floor(self.code().dimension() + self.code().length()) / 2
        s = PolRing.one()
        prev_s = PolRing.zero()

        r = b
        prev_r = a
        while(r.degree() >= stop):
            q = prev_r.quo_rem(r)[0]
            (prev_r, r) = (r, prev_r - q * r)
            (prev_s, s) = (s, prev_s - q * s)

        return (r, s)

    def decode_to_message(self, r):
        r"""
        Decodes ``r`` to an element in message space of ``self``

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self`` message space

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: r = vector(F, (8, 2, 6, 10, 6, 10, 7, 6, 7, 1))
            sage: D.decode_to_message(r)
            x^4 + 7*x^3 + 10*x^2 + 9*x + 8

        If we try to decode a word with too many errors, it returns
        an exception::

            sage: r[0] = r[1] = r[2] = 3
            sage: D.decode_to_message(r)
            Traceback (most recent call last):
            ...
            DecodingError
        """
        C = self.code()
        alphas = C.evaluation_points()
        col_mults = C.column_multipliers()
        PolRing = C.base_field()['x']
        G = self.precompute(PolRing)

        if r in C:
            return self.connected_encoder().unencode_nocheck(r)

        points = [(alphas[i], r[i]/col_mults[i]) for i in \
                range(0, C.length())]
        R = PolRing.lagrange_polynomial(points)

        (Q1, Q0) = self.partial_xgcd(G, R, PolRing)

        if not Q0.divides(Q1):
            raise DecodingError()
        h = Q1//Q0

        return h

    def decoding_radius(self):
        r"""
        Return maximal number of errors that ``self`` can decode

        OUTPUT:

        - the number of errors as an integer

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSGaoDecoder(C)
            sage: D.decoding_radius()
            14
        """
        return (self.code().minimum_distance()-1)//2










class GRSErrorErasureDecoder(Decoder):
    r"""
    Construct a decoder for GRS codes. This decoder is able to correct both errors
    and erasures within its decoding radius.

    INPUT:

    - ``code`` -- The associated code of this decoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: D = codes.decoders.GRSErrorErasureDecoder(C)
        sage: D
        Error-Erasure decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSErrorErasureDecoder(C)
            sage: D
            Error-Erasure decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        input_space = cartesian_product([code.ambient_space(),\
                VectorSpace(GF(2), code.ambient_space().dimension())])
        super(GRSErrorErasureDecoder, self).__init__(code, input_space, "EvaluationVector")

    def __eq__(self, other):
        r"""
        Test equality of GRSErrorErasureDecoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D1 = codes.decoders.GRSErrorErasureDecoder(C)
            sage: D2 = codes.decoders.GRSErrorErasureDecoder(C)
            sage: D1.__eq__(D2)
            True
        """
        return isinstance(other, GRSErrorErasureDecoder) \
                and self.code() == other.code()

    def __ne__(self, other):
        r"""
        Tests difference of GRSErrorErasureDecoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C1 = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C2 = codes.GeneralizedReedSolomonCode(F.list()[:n], k + 1)
            sage: D1 = codes.decoders.GRSErrorErasureDecoder(C1)
            sage: D2 = codes.decoders.GRSErrorErasureDecoder(C2)
            sage: D1.__ne__(D2)
            True
        """
        return not self.__eq__(other)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSErrorErasureDecoder(C)
            sage: D
            Error-Erasure decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Error-Erasure decoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSErrorErasureDecoder(C)
            sage: latex(D)
            \textnormal{Error-Erasure decoder for the }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Error-Erasure decoder for the }%s"\
                % self.code()._latex_()

    def decode_to_message(self, word_and_erasure_vector):
        r"""
        Decode ``word_and_erasure_vector`` to an element in message space
        of ``self``

        INPUT:

        - ``word_and_erasure_vector`` -- a pair of vectors, where
          first element is a codeword of ``self`` and second element
          is a vector of GF(2) containing erasure positions

        OUTPUT:

        - a vector of ``self`` message space

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSErrorErasureDecoder(C)
            sage: r = vector(F, (8, 2, 6, 10, 6, 10, 7, 6, 7, 1))
            sage: e = vector(GF(2), (0, 0, 0, 1, 0, 0, 0, 1, 0, 0))
            sage: w_e = (r, e)
            sage: D.decode_to_message(w_e)
            (8, 9, 10, 7, 1)

        If we try to decode a word with too many erasures, it returns
        an exception::

            sage: e = vector(GF(2), (1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
            sage: w_e = (r, e)
            sage: D.decode_to_message(w_e)
            Traceback (most recent call last):
            ...
            DecodingError: Too many erasures in the received word
        """

        word, erasure_vector = word_and_erasure_vector
        if erasure_vector.hamming_weight() >= self.code().minimum_distance():
            raise DecodingError("Too many erasures in the received word")

        shorten_word = vector(self.code().base_ring(), [word[i] for i in range(len(word))\
                if erasure_vector[i]!=1])
        C1_length = len(shorten_word)
        C1_evaluation_points = [self.code().evaluation_points()[i] for i in\
                range(self.code().length()) if erasure_vector[i]!=1]
        C1_column_multipliers = [self.code().column_multipliers()[i] for i in\
                range(self.code().length()) if erasure_vector[i]!=1]
        C1 = GeneralizedReedSolomonCode(C1_evaluation_points,\
                self.code().dimension(), C1_column_multipliers)
        return C1.decode_to_message(shorten_word)

    def decoding_radius(self, number_erasures):
        r"""
        Return maximal number of errors that ``self`` can decode according
        to how many erasures it receives

        INPUT:

        - ``number_erasures`` -- the number of erasures when we try to decode

        OUTPUT:

        - the number of errors as an integer

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: D = codes.decoders.GRSErrorErasureDecoder(C)
            sage: D.decoding_radius(5)
            12

        If we receive too many erasures, it returns an exception as codeword will
        be impossible to decode::

            sage: D.decoding_radius(30)
            Traceback (most recent call last):
            ...
            ValueError: The number of erasures exceed decoding capability
        """
        diff = self.code().minimum_distance() - number_erasures
        if diff <= 0:
            raise ValueError("The number of erasures exceed decoding capability")
        else :
            return diff // 2








class GRSKeyEquationSyndromeDecoder(Decoder):
    r"""
    Key equation decoding based on the syndrome polynomial.

    This decoder uses early terminated extended euclidean algorithm
    to solve the key equations, as described in [R].

    REFERENCES:

        .. [R06] Introduction to Coding Theory, Ron Roth, Cambridge University Press, 2006

    INPUT:

    - ``code`` -- The associated code of this decoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
        sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
        sage: D
        Key equation decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: D
            Key equation decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        if (code.base_field())(0) in code.evaluation_points():
            raise ValueError("Impossible to decode a GRS code which contains 0 amongst its evaluation points")
        super(GRSKeyEquationSyndromeDecoder, self).__init__(code, code.ambient_space(),\
                "EvaluationVector")

    def __eq__(self, other):
        r"""
        Tests equality of GRSKeyEquationSyndromeDecoder objects.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D1 = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: D2 = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: D1.__eq__(D2)
            True
        """
        return isinstance(other, GRSKeyEquationSyndromeDecoder) \
                and self.code() == other.code()\
                and self.input_space() == other.input_space()

    def __ne__(self, other):
       r"""
       Tests difference of GRSKeyEquationSyndromeDecoder objects.

       EXAMPLES::

           sage: F = GF(59)
           sage: n, k = 40, 12
           sage: C1 = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
           sage: C2 = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k + 1)
           sage: D1 = codes.decoders.GRSKeyEquationSyndromeDecoder(C1)
           sage: D2 = codes.decoders.GRSKeyEquationSyndromeDecoder(C2)
           sage: D1.__ne__(D2)
           True
       """
       return not self.__eq__(other)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: D
            Key equation decoder for the [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Key equation decoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: latex(D)
            \textnormal{Key equation decoder for the }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Key equation decoder for the }%s" % self.code()._latex_()

    def partial_xgcd(self, a, b, PolRing):
        r"""
        Returns the greatest common divisor of ``a`` and ``b``.

        The computation stops whenever the degree of the xgcd falls below
        the degree of ``t``, with ``t`` the Euler coefficient of ``b``.

        This is a helper function, used in :meth:`decode_to_message`.

        INPUT:

        - ``a, b`` -- polynomials over ``PolRing``

        - ``PolRing`` -- polynomial ring of the output

        OUTPUT:

        - a tuple of polynomials ``(r, s)`` where ``r`` is to the xgcd of
          ``a`` and ``b`` and ``s`` is the Bezout coefficient of ``a``.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: P = PolynomialRing(F,'x')
            sage: x = P.parameter()
            sage: a = 5*x^2 + 9*x + 8
            sage: b = 10*x^2 + 3*x + 5
            sage: D.partial_xgcd(a, b, P)
            (5, 8*x + 10)
        """
        stop = (self.code().length() - self.code().dimension()) / 2

        prev_t = PolRing.zero()
        t = PolRing.one()

        prev_r = a
        r = b

        while(r.degree() >= t.degree()):
            q = prev_r.quo_rem(r)[0]
            prev_r, r = r, prev_r - q * r
            prev_t, t = t, prev_t - q * t

        return (r, t)

    def syndrome(self, r):
        r"""
        Returns the coefficients of the syndrome polynomial of ``r``.

        This is a helper function, used in :meth:`decode_to_message`.

        INPUT:

        - ``r`` -- a vector of the ambient space of ``self.code()``

        OUTPUT:

        - a list

        EXAMPLES::


            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: r = vector(F, (8, 2, 6, 10, 6, 10, 7, 6, 7, 2))
            sage: D.syndrome(r)
            [1, 10, 1, 10, 1]
        """
        C = self.code()
        F = C.base_ring()
        S = []
        col_mults = C.parity_column_multipliers()
        alphas = C.evaluation_points()

        for l in range(C.minimum_distance() - 1):
            Sl = F.zero()
            for j in range(C.length()):
                Sl += r[j] * col_mults[j] * (alphas[j] ** l)
            S.append(Sl)

        return S

    def forney_formula(self, error_evaluator, error_locator):
        r"""
        Returns the error vector computed through Forney's formula.

        This is a helper function, used in :meth:`decode_to_message`.

        INPUT:

        - ``error_evaluator``, ``error_locator`` -- two polynomials

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: R.<x> = F[]
            sage: evaluator, locator = R(10), R([10, 10])
            sage: D.forney_formula(evaluator, locator)
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
        """
        C = self.code()
        alphas = C.evaluation_points()
        col_mults = C.parity_column_multipliers()
        ELPp = error_locator.derivative()
        e = []

        for i in range(C.length()):
            alpha_inv = 1/alphas[i]
            if error_locator(alpha_inv) == 0:
                e.append(-alphas[i]/col_mults[i] * error_evaluator(alpha_inv)/ELPp(alpha_inv))
            else:
                e.append(0)

        return vector(C.base_ring(), e)

    def decode_to_code(self, r):
        r"""
        Decodes ``r`` to an element in ``self.code()``.

        INPUT:

        - ``r`` -- a vector of the ambient space of ``self.code()``

        OUTPUT:

        - a vector of ``self.code()``

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: r = vector(F, (8, 2, 6, 10, 6, 10, 7, 6, 7, 2))
            sage: D.decode_to_code(r)
            (8, 2, 6, 10, 6, 10, 7, 6, 7, 1)
        """
        C = self.code()
        F = C.base_field()
        PolRing = C.base_field()['x']
        x = PolRing.gen()

        if r in C:
            return r

        S = PolRing(self.syndrome(r))
        a = x ** (C.minimum_distance() - 1)

        (EEP, ELP) = self.partial_xgcd(a, S, PolRing)

        e = self.forney_formula(EEP, ELP)
        return r - e

    def decoding_radius(self):
        r"""
        Return maximal number of errors that ``self`` can decode

        OUTPUT:

        - the number of errors as an integer

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
            sage: D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
            sage: D.decoding_radius()
            14
        """
        return (self.code().minimum_distance()-1)//2

IMPOSSIBLE_PARAMS = "Impossible parameters for the Guruswami-Sudan algorithm"

def _gs_satisfactory(tau, s, l, C = None, n_k = None):
    r"""
    Returns whether input parameters satisfy the governing equation of
    Guruswami-Sudan.

    See [N13]_ page 49, definition 3.3 and proposition 3.4 for details.

    INPUT:

    - ``tau`` -- an integer, number of errrors one expects Guruswami-Sudan algorithm
      to correct
    - ``s`` -- an integer, multiplicity parameter of Guruswami-Sudan algorithm
    - ``l`` -- an integer, list size parameter
    - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
    - ``n_k`` -- (default: ``None``) a tuple of integers, respectively the
      length and the dimension of the :class:`GeneralizedReedSolomonCode`

    ..NOTE::

        One has to provide either ``C`` or ``(n, k)``. If none or both are
        given, an exception will be raised.

    EXAMPLES::

        sage: from sage.coding.grs import _gs_satisfactory as gs_sat
        sage: tau, s, l = 97, 1, 2
        sage: n, k = 250, 70
        sage: gs_sat(tau, s, l, n_k = (n, k))
        True

    One can also pass a GRS code::

        sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
        sage: gs_sat(tau, s, l, C = C)
        True

    Another example where ``s`` and ``l`` does not satisfy the equation::

        sage: tau, s, l = 118, 47, 80
        sage: gs_sat(tau, s, l, n_k = (n, k))
        False

    If one provides both ``C`` and ``n_k`` an exception is returned::

        sage: tau, s, l = 97, 1, 2
        sage: n, k = 250, 70
        sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
        sage: gs_sat(tau, s, l, C = C, n_k = (n, k))
        Traceback (most recent call last):
        ...
        ValueError: Please provide only the code or its length and dimension

    Same if one provides none of these::

        sage: gs_sat(tau, s, l)
        Traceback (most recent call last):
        ...
        ValueError: Please provide either the code or its length and dimension
    """
    if C is not None and n_k is not None:
        raise ValueError("Please provide only the code or its length and dimension")
    elif C is None and n_k is None:
        raise ValueError("Please provide either the code or its length and dimension")
    elif C is not None:
        n, k = C.length(), C.dimension()
    elif n_k is not None and not isinstance(n_k, tuple):
        raise ValueError("n_k has to be a tuple")
    elif n_k is not None:
        n, k = n_k[0], n_k[1]
    return l > 0 and s > 0 and n * s * (s+1) < (l+1) * (2*s*(n-tau) - (k-1) * l)

def _s_l_from_tau(tau, C = None, n_k = None):
    r"""
    Returns suitable ``s`` and ``l`` according to input parameters.

    See [N13]_ pages 53-54, proposition 3.11 for details.

    INPUT:

    - ``tau`` -- an integer, number of errrors one expects Guruswami-Sudan algorithm
      to correct
    - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
    - ``n_k`` -- (default: ``None``) a tuple of integers, respectively the
      length and the dimension of the :class:`GeneralizedReedSolomonCode`

    OUTPUT:

    - ``(s, l)`` -- a couple of integers, where:
        - ``s`` is the multiplicity parameter of Guruswami-Sudan algorithm and
        - ``l`` is the list size parameter

    ..NOTE::

        One has to provide either ``C`` or ``(n, k)``. If none or both are
        given, an exception will be raised.

    EXAMPLES::

        sage: from sage.coding.grs import _s_l_from_tau as s_l
        sage: tau = 97
        sage: n, k = 250, 70
        sage: s_l(tau, n_k = (n, k))
        (2, 3)

    Same one with a GRS code::

        sage: C = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 70)
        sage: s_l(tau, C = C)
        (2, 3)

    Another one with a bigger ``tau``::

        sage: s_l(118, C = C)
        (47, 89)
    """
    if C is not None and n_k is not None:
        raise ValueError("Please provide only the code or its length and dimension")
    elif C is None and n_k is None:
        raise ValueError("Please provide either the code or its length and dimension")
    elif C is not None:
        n, k = C.length(), C.dimension()
    elif n_k is not None and not isinstance(n_k, tuple):
        raise ValueError("n_k has to be a tuple")
    elif n_k is not None:
        n, k = n_k[0], n_k[1]

    w = k - 1
    atau = n - tau
    smin = tau * w / (atau ** 2 - n * w)
    s = floor(1 + smin)
    D = (s - smin) * (atau ** 2 - n * w) * s + (w**2) /4
    l = floor(atau / w * s + 0.5 - sqrt(D)/w)
    assert _gs_satisfactory(tau,s,l, n_k = (n, k)) , IMPOSSIBLE_PARAMS
    return (s, l)

def _ligt(x):
    r"""
    Returns the least integer greater than ``x``.

    EXAMPLES::

        sage: from sage.coding.grs import _ligt
        sage: _ligt(41)
        42

    It works with any type of numbers (not only integers)::

        sage: _ligt(41.041)
        42
    """
    return floor(x+1)

def _gilt(x):
    r"""
    Returns the greatest integer smaller than ``x``.

    EXAMPLES::

        sage: from sage.coding.grs import _gilt
        sage: _gilt(43)
        42

    It works with any type of numbers (not only integers)::

        sage: _gilt(43.041)
        43
    """
    if x in ZZ:
        return Integer(x-1)
    else:
        return floor(x)

def _solve2deg_int(a,b,c):
    r"""
    Returns the greatest integer range `[i1, i2]` such that
    `i1 > x1` and `i2 < x2` where `x1`,`x2` are the two zeroes of the equation in `x`:
    `ax^2+bx+c=0`.

    If there is no real solution to the equation, it returns an empty, range with negative coefficients.

    INPUT:

    - ``a``, ``b`` and ``c`` -- coefficients of a second degree equation, ``a`` being the coefficient of
      the higher degree term.

    EXAMPLES::

        sage: from sage.coding.grs import _solve2deg_int
        sage: _solve2deg_int(1, -5, 1)
        (1, 4)

    If there is no real solution::

        sage: _solve2deg_int(50, 5, 42)
        (-2, -1)
    """
    D = b**2 - 4*a*c
    if D < 0:
        return (-2,-1)
    sD = float(sqrt(D))
    minx, maxx = (-b-sD)/2.0/a , (-b+sD)/2.0/a
    mini, maxi = (_ligt(minx), _gilt(maxx))
    if mini > maxi:
        return (-2,-1)
    else:
        return (mini,maxi)

def _find_minimal_satisfiable(f, startn=1, contiguous=True):
    r"""
    Returns the minimal integral ``n``, `n > 0` such that ``f(n) == True``.

    If the interval for which `f` is true is contiguous and open
    towards infinity, a logarithmic algorithm is used, otherwise linear.
    `startn` can be given as a hint to a value that might be true.

    INPUT:

    - ``f`` -- a function

    - ``startn`` -- (default: ``1``) the starting point of the algorithm.

    - ``contiguous`` -- (default: ``True``) boolean describing the contiguousity of ``f``'s
      interval

    EXAMPLES::

        sage: from sage.coding.grs import _find_minimal_satisfiable
        sage: def f(x):
        ....:    return None if x > 10 or x == 1 else x + 1

        sage: _find_minimal_satisfiable(f)
        2
    """
    if not contiguous:
        n = startn
        if f(n):
            while f(n) and n > 0:
                n = n - 1
            return n + 1
        else:
            while not f(n):
                n = n + 1
            return n
    else:
        maxn = startn
        minn = 1
        # Keep doubling n to find one that works and then binary
        while not f(maxn):
            minn = maxn + 1
            maxn *= 2
        while minn < maxn:
            tryn = minn + floor((maxn - minn) * 0.5)
            if f(tryn):
                maxn = tryn
            else:
                minn = tryn + 1
        return maxn

def best_s_l_from_tau(tau, C = None, n_k = None):
    r"""
    Returns the best ``s`` and ``l`` possible according to input parameters.

    INPUT:

    - ``tau`` -- an integer, number of errrors one expects Guruswami-Sudan algorithm
      to correct
    - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
    - ``n_k`` -- (default: ``None``) a tuple of integers, respectively the
      length and the dimension of the :class:`GeneralizedReedSolomonCode`

    OUTPUT:

    - ``(s, l)`` -- a couple of integers, where:
        - ``s`` is the multiplicity parameter of Guruswami-Sudan algorithm and
        - ``l`` is the list size parameter

    ..NOTE::

        One has to provide either ``C`` or ``(n, k)``. If none or both are
        given, an exception will be raised.

    EXAMPLES::

        sage: from sage.coding.grs import best_s_l_from_tau
        sage: tau, n, k = 97, 250, 70
        sage: best_s_l_from_tau(tau, n_k = (n, k))
        (1, 2)

    Another one with a bigger tau::

        sage: tau, n, k = 118, 250, 70
        sage: best_s_l_from_tau(tau, n_k = (n, k))
        (47, 89)
    """
    if C is not None and n_k is not None:
        raise ValueError("Please provide only the code or its length and dimension")
    elif C is None and n_k is None:
        raise ValueError("Please provide either the code or its length and dimension")
    elif C is not None:
        n, k = C.length(), C.dimension()
    elif n_k is not None and not isinstance(n_k, tuple):
        raise ValueError("n_k has to be a tuple")
    elif n_k is not None:
        n, k = n_k[0], n_k[1]
    (firsts, firstl) = _s_l_from_tau(tau, n_k = (n, k))
    def try_l(l):
        (mins,maxs) = _solve2deg_int(n, n-2*(l+1)*(n-tau), (k-1)*l*(l+1))
        if maxs > 0 and maxs >= mins:
            return max(1, mins)
        else:
            return None
    l = _find_minimal_satisfiable(try_l, firstl)
    s = try_l(l)
    assert _gs_satisfactory(tau, s, l, n_k = (n, k)) , IMPOSSIBLE_PARAMS
    return (s, l)

def list_decoding_range(n, d, q=None):
    r"""
    Returns the minimal and maximal number of errors correctable by a
    Johnson-distance list decoder beyond half the minimal distance.

    INPUT:

    - ``n`` -- an integer, the length of the code
    - ``d`` -- an integer, the minimum distance of the code
    - ``q`` -- (default: ``None``) ????????

    EXAMPLES::

        sage: from sage.coding.grs import list_decoding_range
        sage: list_decoding_range(250, 181)
        (91, 118)
    """
    if q is None:
        return (_ligt((d-1)/2), _gilt(n - sqrt(n*(n-d))))
    else:
        f = (q-1.)/q
        return (_ligt((d-1)/2), _gilt(f*(n-sqrt(n*(n-d/f)))))

def guruswami_sudan_decoding_radius(C = None, n_k = None, l = None, s = None):
    r"""
    Returns the maximal decoding radius of the Guruswami-Sudan decoder.

    If ``s`` is set but ``l`` is not it will return the best decoding radius using this ``s``
    alongside with the required ``l``, same for ``l``.

    INPUT:

    - ``C`` -- (default: ``None``) a :class:`GeneralizedReedSolomonCode`
    - ``n_k`` -- (default: ``None``) a tuple of integers, respectively the
      length and the dimension of the :class:`GeneralizedReedSolomonCode`
    - ``s`` -- (default: ``None``) an integer, the multiplicity parameter of Guruswami-Sudan algorithm
    - ``l`` -- (default: ``None``) an integer, the list size parameter

    ..NOTE::

        One has to provide either ``C`` or ``(n, k)``. If none or both are
        given, an exception will be raised.

    EXAMPLES::

        sage: from sage.coding.grs import guruswami_sudan_decoding_radius
        sage: n, k = 250, 70
        sage: guruswami_sudan_decoding_radius(n_k = (n, k))
        (118, (47, 89))
    """
    if C is not None and n_k is not None:
        raise ValueError("Please provide only the code or its length and dimension")
    elif C is None and n_k is None:
        raise ValueError("Please provide either the code or its length and dimension")
    elif C is not None:
        n, k = C.length(), C.dimension()
    elif n_k is not None and not isinstance(n_k, tuple):
        raise ValueError("n_k has to be a tuple")
    elif n_k is not None:
        n, k = n_k[0], n_k[1]

    def get_tau(s,l):
        if s<=0 or l<=0:
            return -1
        return _gilt(n - n/2*(s+1)/(l+1) - (k-1)/2*l/s)
    if l==None and s==None:
        tau = list_decoding_range(n,n-k+1)[1]
        return (tau, best_s_l_from_tau(tau, n_k = (n, k)))
    if l!=None and s!=None:
        return (get_tau(s,l), (s,l))
    if s!= None:
        # maximising tau under condition
        # n*(s+1 choose 2) < (ell+1)*s*(n-tau) - (ell+1 choose 2)*(k-1)
        # knowing n and s, we can just minimise
        # ( n*(s+1 choose 2) + (ell+1 choose 2)*(k-1) )/(ell+1)
        # Differentiating and setting to zero yields ell best choice:
        lmax = sqrt(n*s*(s+1)/(k-1)) - 1
        #the best integral value will be
        (l,tau) = find_integral_max(lmax, lambda l: get_tau(s,l))
        assert _gs_satisfactory(tau,s,l, n_k = (n, k)), IMPOSSIBLE_PARAMS
        #Note that we have not proven that this ell is minimial in integral
        #sense! It just seems that this most often happens
        return (tau,(s,l))
    if l!= None:
        # Acquired similarly to when restricting s
        smax = sqrt((k-1)/n*l*(l+1))
        (s,tau) = find_integral_max(smax, lambda s: get_tau(s,l))
        assert _gs_satisfactory(tau,s,l, n_k = (n, k)), IMPOSSIBLE_PARAMS
        return (get_tau(s,l), (s,l))

def poly2list(p, len):
    """Convert the polynomial p into a list of coefficients of length len"""
    return list(p) + [0]*max(0, len-p.degree()-1)

class GRSGuruswamiSudanDecoder(Decoder):

    from gs_interpolation import gs_construct_Q_linalg as gs_construct_Q_linalg
    from gs_rootfinding import *

    def __init__(self, C, tau = None, params = None, Qfinder = None, root_finder = None):
        n, k = C.length(), C.dimension()
        if tau:
            self.tau = tau
            self.s, self.ell = best_s_l_from_tau(tau, n_k = (n, k))
        elif params:
            self.s = params[0]
            self.ell = params[1]
            (self.tau,_) = guruswami_sudan_decoding_radius(n, k, s=self.s, l=self.ell)
        else:
            raise Exception("Specify either tau or params")
        #TODO: Precompute stuff for the various Qfinders
        self.Qfinder = Qfinder if Qfinder else gs_construct_Q_linalg
        self.root_finder = root_finder if root_finder else rootfind_bivariate
        super(GRSGuruswamiSudanDecoder, self).__init__(C, C.ambient_space(), C._default_encoder_name)

    def decode_to_code(self, r):
        C = self.code()
        n,k,d,alphas,colmults = C.length(), C.dimension(), C.minimum_distance(),\
                C.evaluation_points(), C.column_multipliers()
        ## SETUP INTERPOLATION PROBLEM
        wy = k-1
        points = [ (alphas[i], r[i]/colmults[i]) for i in range(0,len(alphas)) ]
        ## SOLVE INTERPOLATION
        Q = self.Qfinder(points, self.tau, (self.s,self.ell), wy)
        ## EXAMINE THE FACTORS AND CONVERT TO CODEWORDS
        #factors = self.root_finder(Q, maxd=None, algorithm = "roth_ruckenstein")
        factors = rootfind_bivariate(Q, algorithm = "roth_ruckenstein")
        if not factors:
            return None
        return [ poly2list(f, k) for f in factors ]
####################### registration ###############################

GeneralizedReedSolomonCode._registered_encoders["EvaluationVector"] = GRSEvaluationVectorEncoder
GeneralizedReedSolomonCode._registered_encoders["EvaluationPolynomial"] = GRSEvaluationPolynomialEncoder

GeneralizedReedSolomonCode._registered_decoders["BerlekampWelch"] = GRSBerlekampWelchDecoder
GRSBerlekampWelchDecoder._decoder_type = {"hard-decision", "unique", "always-succeed"}
GeneralizedReedSolomonCode._registered_decoders["Gao"] = GRSGaoDecoder
GRSGaoDecoder._decoder_type = {"hard-decision", "unique", "always-succeed"}
GeneralizedReedSolomonCode._registered_decoders["ErrorErasure"] = GRSErrorErasureDecoder
GRSErrorErasureDecoder._decoder_type = {"error-erasure", "unique", "always-succeed"}
GeneralizedReedSolomonCode._registered_decoders["KeyEquationSyndrome"] = GRSKeyEquationSyndromeDecoder
GRSKeyEquationSyndromeDecoder._decoder_type = {"hard-decision", "unique", "always-succeed"}

GeneralizedReedSolomonCode._registered_decoders["Syndrome"] = LinearCodeSyndromeDecoder
GeneralizedReedSolomonCode._registered_decoders["NearestNeighbor"] = LinearCodeNearestNeighborDecoder
