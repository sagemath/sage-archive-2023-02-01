r"""
Generalized Reed-Solomon code

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
from linear_code import (AbstractLinearCode,
                         LinearCodeSyndromeDecoder,
                         LinearCodeNearestNeighborDecoder)
from encoder import Encoder
from sage.misc.misc_c import prod
from sage.functions.other import binomial
from sage.calculus.var import var
from sage.misc.functional import symbolic_sum
from sage.rings.integer_ring import ZZ

class GeneralizedReedSolomonCode(AbstractLinearCode):
    r"""
    Representation of a Generalized Reed-Solomon code.

    INPUT:

    - ``evaluation_points`` -- A list of distinct elements of some finite field `F`.

    - ``dimension`` -- The dimension of the resulting code.

    - ``column_multipliers`` -- (default: ``None``) List of non-zero elements of `F`.
      All column multipliers are set to 1 if default value is kept.

    EXAMPLES:

    A Reed-Solomon code can be constructed by taking all non-zero elements of
    the field as evaluation points, and specifying no column multipliers::

        sage: F = GF(7)
        sage: evalpts = [F(i) for i in range(1,7)]
        sage: C = codes.GeneralizedReedSolomonCode(evalpts,3)
        sage: C
        [6, 3, 4] Generalized Reed-Solomon Code over Finite Field of size 7

    More generally, the following is a GRS code where the evaluation points are
    a subset of the field and includes zero::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: C
        [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59

    It is also possible to specify the column multipliers::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: colmults = F.list()[1:n+1]
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, colmults)
        sage: C
        [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, evaluation_points, dimension, column_multipliers=None):
        r"""
        TESTS:

        If the evaluation points are not from a finite field, it raises an error::

            sage: C = codes.GeneralizedReedSolomonCode([1,2,3], 1)
            Traceback (most recent call last):
            ...
            ValueError: Evaluation points must be in a finite field (and Integer Ring is not one)

        If the evaluation points are not from the same finite field, it raises an error::

            sage: F2, F3 = GF(2) , GF(3)
            sage: C = codes.GeneralizedReedSolomonCode([F2.zero(),F2.one(),F3(2)], 1)
            Traceback (most recent call last):
            ...
            ValueError: Failed converting all evaluation points to the same field (unable to find a common ring for all elements)

        If the column multipliers cannot be converted into the finite are not from a finite field, or cannot be not in the same
        finite field as the evaluation points, it raises an error::

            sage: F = GF(59)
            sage: F2 = GF(61)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, [.3]*n )
            Traceback (most recent call last):
            ...
            ValueError: Failed converting all evaluation points and column multipliers to the same field (unable to find a common ring for all elements)

            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, F2.list()[1:n+1])
            Traceback (most recent call last):
            ...
            ValueError: Failed converting all evaluation points and column multipliers to the same field (unable to find a common ring for all elements)

        The number of column multipliers is checked as well::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, F.list()[1:n])
            Traceback (most recent call last):
            ...
            ValueError: There must be the same number of evaluation points as column multipliers

        It is not allowed to have 0 as a column multiplier::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k, F.list()[:n])
            Traceback (most recent call last):
            ...
            ValueError: All column multipliers must be non-zero

        And all the evaluation points must be different. Note that they should
        be different after converting into the same field::

            sage: F = GF(5)
            sage: C = codes.GeneralizedReedSolomonCode([ F(0), 1, 2, 3, 5 ], 3)
            Traceback (most recent call last):
            ...
            ValueError: All evaluation points must be different

        The dimension is not allowed to exceed the length::

            sage: F = GF(59)
            sage: n, k = 40, 100
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            Traceback (most recent call last):
            ...
            ValueError: The dimension must be a positive integer at most the length of the code.
        """
        if column_multipliers:
            if len(evaluation_points) != len(column_multipliers):
                raise ValueError("There must be the same number of evaluation points as column multipliers");
            try:
                common_points = vector(list(evaluation_points) + list(column_multipliers))
                F = common_points.base_ring()
                self._evaluation_points = common_points[:len(evaluation_points)]
                self._column_multipliers = common_points[len(evaluation_points):]
            except (TypeError, ValueError) as e:
                raise ValueError("Failed converting all evaluation points and column multipliers to the same field (%s)" % e.message)
        else:
            try:
                self._evaluation_points = vector(evaluation_points)
                F = self._evaluation_points.base_ring()
                self._column_multipliers = vector(F, [F.one()] * len(self._evaluation_points))
            except (TypeError, ValueError) as e:
                raise ValueError("Failed converting all evaluation points to the same field (%s)" % e.message)

        if F.is_finite() == False or F.is_field() == False:
            raise ValueError("Evaluation points must be in a finite field (and %s is not one)" % F)
        super(GeneralizedReedSolomonCode, self).__init__(F, \
                len(self._evaluation_points), "EvaluationVector", "Syndrome")

        if dimension not in ZZ or dimension > self._length or dimension < 1:
            raise ValueError("The dimension must be a positive integer at most the length of the code.")
        self._dimension = dimension

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
        Returns the minimum distance of ``self``. Since a GRS code is always
        Maximum-Distance-Separable (MDS), this returns ``C.length() -
        C.dimension() + 1``.

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
        Returns the evaluation points of ``self`` as a vector.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.evaluation_points()
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
        """
        return self._evaluation_points

    def column_multipliers(self):
        r"""
        Returns the column multipliers of ``self`` as a vector.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10, 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.column_multipliers()
            (1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
        """
        return self._column_multipliers

    @cached_method
    def multipliers_product(self):
        r"""
        Returns the component-wise product of the column multipliers of ``self``
        with the column multipliers of the dual GRS code.

        This is a simple Cramer's rule-like expression on the evaluation points
        of ``self``. Recall that the column multipliers of the dual GRS code is
        also the column multipliers of the parity check matrix of ``self``.

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
        return self.dual_code().generator_matrix()

    @cached_method
    def dual_code(self):
        r"""
        Returns the dual code of ``self``, which is also a GRS code.

        EXAMPLES::

            sage: F =  GF(59)
            sage: colmults = [ F.random_element() for i in range(40) ]
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:40], 12, colmults)
            sage: Cd = C.dual_code(); Cd
            [40, 28, 13] Generalized Reed-Solomon Code over Finite Field of size 59

        The dual code of the dual code is the original code::

            sage: C == Cd.dual_code()
            True
        """
        col_mults = self.parity_column_multipliers()
        return GeneralizedReedSolomonCode(self.evaluation_points(), self.length() - self.dimension(), col_mults)

    def covering_radius(self):
        r"""
        Returns the covering radius of ``self``.

        The covering radius of a linear code `C` is the smallest
        number `r` s.t. any element of the ambient space of `C` is at most at
        distance `r` to `C`.

        As GRS codes are Maximum Distance Separable codes (MDS), their covering
        radius is always `d-1`, where `d` is the minimum distance. This is
        opposed to random linear codes where the covering radius is
        computationally hard to determine.

        EXAMPLES::

            sage: F = GF(2^8, 'a')
            sage: n, k = 256, 100
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: C.covering_radius()
            156
        """
        return self.length() - self.dimension()

    @cached_method
    def weight_distribution(self):
        r"""
        Returns the list whose `i`'th entry is the number of words of weight `i`
        in ``self``.

        Computing the weight distribution for a GRS code is very fast. Note that
        for random linear codes, it is computationally hard.

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
        Returns the polynomial whose coefficient to `x^i` is the number of codewords of weight `i` in ``self``.

        Computing the weight enumerator for a GRS code is very fast. Note that
        for random linear codes, it is computationally hard.

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








class GRSEvaluationVectorEncoder(Encoder):
    r"""
    Encoder for Generalized Reed-Solomon codes which encodes vectors into codewords.

    INPUT:

    - ``code`` -- The associated code of this encoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
        sage: E
        Evaluation vector-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59

    Actually, we can construct the encoder from ``C`` directly::

        sage: E = C.encoder("EvaluationVector")
        sage: E
        Evaluation vector-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: E
            Evaluation vector-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
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
            sage: D1 is D2
            False
        """
        return isinstance(other, GRSEvaluationVectorEncoder) \
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: E
            Evaluation vector-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Evaluation vector-style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationVectorEncoder(C)
            sage: latex(E)
            \textnormal{Evaluation vector-style encoder for }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Evaluation vector-style encoder for }%s" % self.code()._latex_()

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
        C = self.code()
        alphas = C.evaluation_points()
        col_mults = C.column_multipliers()
        return matrix(C.base_field(), C.dimension(), C.length(), lambda i,j: col_mults[j] * alphas[j]**i)








class GRSEvaluationPolynomialEncoder(Encoder):
    r"""
    Encoder for Generalized Reed-Solomon codes which uses evaluation of
    polynomials to obtain codewords.

    INPUT:

    - ``code`` -- The associated code of this encoder.

    EXAMPLES::

        sage: F = GF(59)
        sage: n, k = 40, 12
        sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
        sage: E = codes.encoders.GRSEvaluationPolynomialEncoder(C)
        sage: E
        Evaluation polynomial-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        sage: E.message_space()
        Univariate Polynomial Ring in x over Finite Field of size 59

    Actually, we can construct the encoder from ``C`` directly::

        sage: E = C.encoder("EvaluationPolynomial")
        sage: E
        Evaluation polynomial-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = codes.encoders.GRSEvaluationPolynomialEncoder(C)
            sage: E
            Evaluation polynomial-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
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
            sage: D1 is D2
            False
            sage: D1.__eq__(D2)
            True
        """
        return isinstance(other, GRSEvaluationPolynomialEncoder) \
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = C.encoder("EvaluationPolynomial")
            sage: E
            Evaluation polynomial-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        return "Evaluation polynomial-style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(59)
            sage: n, k = 40, 12
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = C.encoder("EvaluationPolynomial")
            sage: latex(E)
            \textnormal{Evaluation polynomial-style encoder for }[40, 12, 29] \textnormal{ Generalized Reed-Solomon Code over } \Bold{F}_{59}
        """
        return "\\textnormal{Evaluation polynomial-style encoder for }%s" % self.code()._latex_()

    def encode(self, p):
        r"""
        Transforms the polynomial ``p`` into a codeword of :meth:`code`.

        INPUT:

        - ``p`` -- A polynomial from the message space of ``self`` of degree
          less than ``self.code().dimension()``.

        OUTPUT:

        - A codeword in associated code of ``self``

        EXAMPLES::

            sage: F = GF(11)
            sage: Fx.<x> = F[]
            sage: n, k = 10 , 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = C.encoder("EvaluationPolynomial")
            sage: p = x^2 + 3*x + 10
            sage: c = E.encode(p); c
            (10, 3, 9, 6, 5, 6, 9, 3, 10, 8)
            sage: c in C
            True

        If a polynomial of too high degree is given, an error is raised::

            sage: p = x^10
            sage: E.encode(p)
            Traceback (most recent call last):
            ...
            ValueError: The polynomial to encode must have degree at most 4

        If ``p`` is not an element of the proper polynomial ring, an error is raised::

            sage: Qy.<y> = QQ[]
            sage: p = y^2 + 1
            sage: E.encode(p)
            Traceback (most recent call last):
            ...
            ValueError: The value to encode must be in Univariate Polynomial Ring in x over Finite Field of size 11
        """
        M = self.message_space()
        if p not in M:
            raise ValueError("The value to encode must be in %s" % M)
        C = self.code()
        if p.degree() >= C.dimension():
            raise ValueError("The polynomial to encode must have degree at most %s" % (C.dimension() - 1))
        alphas    = C.evaluation_points()
        col_mults = C.column_multipliers()
        c = vector(C.base_ring(), [col_mults[i]*p(alphas[i]) for i in range(C.length())])
        return c

    def unencode_nocheck(self, c):
        r"""
        Returns the message corresponding to the codeword ``c``.

        Use this method with caution: it does not check if ``c``
        belongs to the code, and if this is not the case, the output is
        unspecified. Instead, use :meth:`unencode`.

        INPUT:

        - ``c`` -- A codeword of :meth:`code`.

        OUTPUT:

        - An polynomial of degree less than ``self.code().dimension()``.

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10 , 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = C.encoder("EvaluationPolynomial")
            sage: c = vector(F, (10, 3, 9, 6, 5, 6, 9, 3, 10, 8))
            sage: c in C
            True
            sage: p = E.unencode_nocheck(c); p
            x^2 + 3*x + 10
            sage: E.encode(p) == c
            True

        Note that no error is thrown if ``c`` is not a codeword, and that the
        result is undefined::

            sage: c = vector(F, (11, 3, 9, 6, 5, 6, 9, 3, 10, 8))
            sage: c in C
            False
            sage: p = E.unencode_nocheck(c); p
            6*x^4 + 6*x^3 + 2*x^2
            sage: E.encode(p) == c
            False

        """
        C = self.code()
        alphas    = C.evaluation_points()
        col_mults = C.column_multipliers()

        c = [c[i]/col_mults[i] for i in range(C.length())]
        points = [(alphas[i], c[i]) for i in range(C.dimension())]

        Pc = self._R.lagrange_polynomial(points)
        return Pc

    def message_space(self):
        r"""
        Returns the message space of ``self``

        EXAMPLES::

            sage: F = GF(11)
            sage: n, k = 10 , 5
            sage: C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
            sage: E = C.encoder("EvaluationPolynomial")
            sage: E.message_space()
            Univariate Polynomial Ring in x over Finite Field of size 11
        """
        return self._R

####################### registration ###############################

GeneralizedReedSolomonCode._registered_encoders["EvaluationVector"] = GRSEvaluationVectorEncoder
GeneralizedReedSolomonCode._registered_encoders["EvaluationPolynomial"] = GRSEvaluationPolynomialEncoder
GeneralizedReedSolomonCode._registered_decoders["Syndrome"] = LinearCodeSyndromeDecoder
GeneralizedReedSolomonCode._registered_decoders["NearestNeighbor"] = LinearCodeNearestNeighborDecoder
