r"""
Golay code

Given the parmeters, alphabet ["binary"/GF(2), "ternary"/GF(3)], and extended ["true"/"false"],
returns the corresponding Golay code.

This file contains the following elements:

    - :class:`GolayCode`, the class for Golay codes
    - :class:`GolayGeneratorMatrixEncoder`, an encoder that uses the generator matrix
    - :class:`GolaySyndromeDecoder`, a decoder which corrects errors using the syndrome
"""

#*****************************************************************************
#       Copyright (C) 2016 Arpit Merchant <arpitdm@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.categories.cartesian_product import cartesian_product
from sage.modules.free_module_element import vector
from sage.modules.free_module import VectorSpace
from sage.rings.integer import Integer
from .linear_code import (AbstractLinearCode,
                         LinearCodeSyndromeDecoder,
                         LinearCodeNearestNeighborDecoder)
from .encoder import Encoder
from .decoder import Decoder, DecodingError
from sage.rings.arith import xgcd
from sage.misc.misc_c import prod
from sage.functions.other import binomial, floor, sqrt
from sage.calculus.var import var
from sage.misc.functional import symbolic_sum
from sage.rings.integer_ring import ZZ

class GolayCode(AbstractLinearCode):
    r"""
    Representation of a Golay Code.

    INPUT:

    - ``alphabet`` -- Either ``"binary"``, ``"ternary"``, ``GF(2)`` or ``GF(3)``.
      The created Golay code will be defined over the specified alphabet.

    - ``extended`` -- (default: ``True``) if set to ``True``, creates an extended Golay
      code.

    EXAMPLES::

    A Golay Code can be constructed by specifying the alphabet and type parameters.

        sage: C = GolayCode(binary, true)
        sage: C
        [24, 12, 8] Extended Binary Golay Code over Finite Field of size 2.

        sage: C = GolayCode(binary, false)
        sage: C
        [23, 12, 7] Binary Golay Code over Finite Field of size 2.

        sage: C = GolayCode(ternary, true)
        sage: C
        [12, 6, 6] Extended Ternary Golay Code over Finite Field of size 3.

        sage: C = GolayCode(ternary, false)
        sage: C
        [11, 6, 5] Ternary Golay Code over Finite Field of size 3.
    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, alphabet, extended=True):
        r"""
        TESTS:

        If the alphabet is not binary or ternary, it raises an error:
            sage: C = codes.GolayCode(quaternary, true)
            Traceback (most recent call last):
            ...
            ValueError: Finite field must be set to 'binary' or 'ternary'

        If the parameter ``extended`` is not "true" or "false", it raises an error:
            sage: C = codes.GolayCode(binary, yes)
            Traceback (most recent call last):
            ...
            ValueError: Extension parameter of the code must be set to 'true' or 'false'
        """
        if alphabet != "binary" or alphabet != "ternary":
            raise ValueError("Finite field must be set to 'binary' or 'ternary'")
        if extended != "true" or extended != "false":
            raise ValueError("Extension parameter of the code must be set to 'true' or 'false'")

        if alphabet == "binary":
            F = GF(2)
            self._dimension = 12
            if extended == "true":
                length = 24
            elif extended == "false":
                length = 23
        elif alphabet == "ternary":
            F = GF(3)
            self._dimension = 6
            if extended == "true":
                length = 12
            elif extended == "false":
                length = 11
        super(GolayCode, self).__init__(F, length, "GolayGeneratorMatrixEncoder", "GolaySyndromeDecoder")

    def __eq__(self, other):
        r"""
        Test equality between Golay Code objects.

        EXAMPLES::

            sage: C1 = codes.GolayCode(binary, true)
            sage: C2 = codes.GolayCode(binary, true)
            sage: C1.__eq__(C2)
            True

            sage: C1 = codes.GolayCode(binary, true)
            sage: C2 = codes.GolayCode(ternary, true)
            sage: C1.__eq__(C2)
            False
        """
        return isinstance(other, GolayCode) \
                and self.base_field() == other.base_field() \
                and self.length() == other.length() \

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.GolayCode(binary, true)
            sage: C
            [24, 12, 8] Golay Code over Finite Field of size 2
        """
        return "[%s, %s, %s] Golay Code over %s"\
                % (self.length(), self.dimension(),\
                self.minimum_distance(), self.base_field())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.GolayCode(binary, true)
            sage: latex(C)
            [24, 12, 8] \textnormal{ Golay Code over } \Bold{F}_{2}
        """
        return "[%s, %s, %s] \\textnormal{ Golay Code over } %s"\
                % (self.length(), self.dimension() ,self.minimum_distance(),
                self.base_field()._latex_())

    def dual_code(self):
        r"""
        Return the dual code of ``self``.

        As Golay Codes are self-dual codes, it just returns ``self``.

        EXAMPLES::

            sage: C = codes.GolayCode(binary, true)
            sage: Cd = C.dual_code(); Cd
            [24, 12, 8] Golay Code over Finite Field of size 2
        """
        return self

    def minimum_distance(self):
        r"""
        Return the minimum distance of ``self``.

        The minimum distance of Golay codes is already known,
        and is thus returned immediately without computing anything.

        EXAMPLES::

            sage: C = codes.GolayCode(binary, true)
            sage: C.minimum_distance()
            8
        """
        n = self.length()
        if n == 24:
            return 8
        elif n == 23:
            return 7
        elif n == 12:
            return 6
        elif n == 11:
            return 5

    def covering_radius(self):
        r"""
        Return the covering radius of ``self``.

        The covering radius of a linear code `C` is the smallest
        integer `r` s.t. any element of the ambient space of `C` is at most at
        distance `r` to `C`.

        The covering radii of all Golay codes are known, and are thus returned
        by this method without performing any computation

        EXAMPLES::

            sage: C = codes.GolayCode("ternary", False)
            sage: C.covering_radius()
            2
        """
        n = self.length()
        if n == 23:
            return 3
        elif n == 24:
            return 4
        elif n = 11:
            return 2
        elif n == 12:
            return 3

    def weight_distribution(self):
        r"""
        Return the list whose `i`'th entry is the number of words of weight `i`
        in ``self``.

        The weight distribution of all Golay codes are known, and are thus returned
        by this method without performing any computation
        MWS (67, 69)

        EXAMPLES::

            sage: C = codes.GolayCode(GF(3))
            sage: C.weight_distribution()
            [1, 0, 0, 0, 0, 0, 264, 0, 0, 440, 0, 0, 24]
        """
        n = self.length()
        if n == 23:
            return ([1]+[0]*6+[253]+[506]+[0]*2+[1288]*[1288]+[0]*2+[506]
                    +[253]+[0]*6+[1])
        if n == 24:
            return ([1]+[0]*7+[759]+[0]*3+[2576]+[0]*3+[759]+[0]*7+[1])
        if n == 11:
            return [1]*[0]*4+[132]*2+[0]+[330]+[110]+[0]+[24]
        if n == 12:
            return [1]+[0]*5+[264]+[0]*2+[440]+[0]*2+[24]

    def weight_enumerator(self):
        r"""
        Return the polynomial whose coefficient to `x^i` is the number of
        codewords of weight `i` in ``self``.

        The weight distribution of all Golay codes are known, and are thus returned
        by this method without performing any computation

        EXAMPLES::

            sage: C = codes.GolayCode(GF(3))
            sage: C.weight_enumerator()
            1 + 264*x**6 + 440*x**9 + 24*x**12
        """
        R = PolynomialRing(ZZ, "x")
        x = R.gen()
        n = self.length()
        if n == 23:
            return (1 + 253*x**7 + 506*x**8 + 1288*x**11 +
                1288*x**12 + 506*x**15 + 253*x**16 + x**23)
        if n == 24:
            return 1 + 759*x**8 + 2576*x**12 + 759*x**16 + x**24
        if n == 11:
            return 1 + 132*x**5 + 132*x**6 + 330*x**8 + 110*x**9 + 24*x**11
        if n == 12:
            return 1 + 264*x**6 + 440*x**9 + 24*x**12








####################### encoders ###############################


class GolayGeneratorMatrixEncoder(Encoder):
    r"""
    Encoder for Golay codes which encodes vectors into codewords.

    INPUT:

    - ``code`` -- The associated code of this encoder.

    EXAMPLES::

        sage: C = codes.GolayCode(binary, true)
        sage: E = codes.encoders.GolayGeneratorMatrixEncoder(C)
        sage: E
        Generator matrix-style encoder for [24, 12, 8] Golay Code over Finite Field of size 2

    Alternatively, we can construct the encoder from ``C`` directly::

        sage: E = C.encoder("GeneratorMatrix")
        sage: E
        Generator matrix-style encoder for [24, 12, 8] Golay Code over Finite Field of size 2
    """
    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: C = codes.GolayCode(binary, true)
            sage: E = codes.encoders.GolayGeneratorMatrixEncoder(C)
            sage: E
            Generator matrix-style encoder for [24, 12, 8] Golay Code over Finite Field of size 2
        """
        super(GolayGeneratorMatrixEncoder, self).__init__(code)
        self._R = code.base_field()['x']

    def __eq__(self, other):
        r"""
        Tests equality between GolayGeneratorMatrixEncoder objects.

        EXAMPLES::

            sage: C = codes.GolayCode(binary, true)
            sage: E1 = codes.encoders.GolayGeneratorMatrixEncoder(C)
            sage: E2 = codes.encoders.GolayGeneratorMatrixEncoder(C)
            sage: E1.__eq__(E2)
            True
        """
        return isinstance(other, GolayGeneratorMatrixEncoder) \
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.GolayCode(binary, true)
            sage: E = codes.encoders.GolayGeneratorMatrixEncoder(C)
            sage: E
            Generator matrix-style encoder for [24, 12, 8] Golay Code over Finite Field of size 2
        """
        return "Generator matrix-style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.GolayCode(binary, true)
            sage: E = codes.encoders.GolayGeneratorMatrixEncoder(C)
            sage: latex(E)
            \textnormal{Generator matrix-style encoder for }[24, 12, 8] \textnormal{ Golay Code over } \Bold{F}_{2}
        """
        return "\\textnormal{Generator matrix-style encoder for }%s" % self.code()._latex_()

    def generator_matrix(self):
        r"""
        Returns a generator matrix of ``self``

        EXAMPLES::

            sage: C = codes.GolayCode(binary, true)
            sage: E = codes.encoders.GolayGeneratorMatrixEncoder(C)
            sage: E.generator_matrix()
            [1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1]
            [0 1 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 1 0 0 0 1 0]
            [0 0 1 0 0 0 0 0 0 0 0 0 1 1 0 1 1 1 0 0 0 1 0 1]
            [0 0 0 1 0 0 0 0 0 0 0 0 1 0 1 1 1 0 0 0 1 0 1 1]
            [0 0 0 0 1 0 0 0 0 0 0 0 1 1 1 1 0 0 0 1 0 1 1 0]
            [0 0 0 0 0 1 0 0 0 0 0 0 1 1 1 0 0 0 1 0 1 1 0 1]
            [0 0 0 0 0 0 1 0 0 0 0 0 1 1 0 0 0 1 0 1 1 0 1 1]
            [0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 1 1 0 1 1 1]
            [0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 1 0 1 1 0 1 1 1 0]
            [0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 1 1 0 1 1 1 0 0]
            [0 0 0 0 0 0 0 0 0 0 1 0 1 1 0 1 1 0 1 1 1 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 1 0 1 1 1 0 0 0 1]
        """
        C = self.code()
        if C.base_field() == GF(2):
        G = matrix(GF(2),
            [[1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
             [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1]])



