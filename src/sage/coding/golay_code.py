r"""
Golay code

Given the parmeters, alphabet ["binary"/GF(2), "ternary"/GF(3)], and extended ["true"/"false"],
returns the corresponding Golay code.

This file contains the following elements:

    - :class:`GolayCode`, the class for Golay codes
    - :class:`GolayGeneratorMatrixEncoder`, an encoder that uses the generator matrix
"""

#*****************************************************************************
#       Copyright (C) 2016 Arpit Merchant <arpitdm@gmail.com>
#                     2016 David Lucas    <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from .linear_code import (AbstractLinearCode,
                          LinearCodeGeneratorMatrixEncoder)

class GolayCode(AbstractLinearCode):
    r"""
    Representation of a Golay Code.

    INPUT:

    - ``base_field`` -- The base field over which the code is defined.
      Can only be ``GF(2)`` or ``GF(3)``.

    - ``extended`` -- (default: ``True``) if set to ``True``, creates an extended Golay
      code.

    EXAMPLES::

        sage: codes.GolayCode(GF(2))
        [24, 6] extended Golay code over GF(2)

    Another example with the perfect binary Golay code::

        sage: codes.GolayCode(GF(2), False)
        [23, 6] Golay code over GF(2)
    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, extended=True):
        r"""
        TESTS:

        If ``base_field`` is not ``GF(2)`` or ``GF(3)``, it raises an error:
            sage: C = codes.GolayCode(ZZ, true)
            Traceback (most recent call last):
            ...
            ValueError: finite_field must be either GF(2) or GF(3)
        """
        if base_field not in [GF(2), GF(3)]:
            raise ValueError("finite_field must be either GF(2) or GF(3)")
        if extended not in [True, False]:
            raise ValueError("extension must be either True or False")

        if base_field is GF(2):
            length = 23
            self._dimension = 12
        else:
            length = 11
            self._dimension = 6
        if extended:
            length += 1
        super(GolayCode, self).__init__(base_field, length, "GeneratorMatrix", "Syndrome")

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

        If ``self`` is an extended Golay code, ``self`` is returned.
        Otherwise, it returns the output of
        :meth:`sage.coding.linear_code.AbstractLinearCode.dual_code`

        EXAMPLES::

            sage: C = codes.GolayCode(binary, true)
            sage: Cd = C.dual_code(); Cd
            [24, 12, 8] Golay Code over Finite Field of size 2
        """
        n = self.length()
        if n % 2 == 0:
            return self
        return super(GolayCode, self).dual_code()

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
        elif n == 11:
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

    def generator_matrix(self):
        r"""
        Return a generator matrix of ``self``

        Generator matrices of all Golay codes are known, and are thus returned
        by this method without performing any computation

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
        n = self.length()
        if n == 23:
            G = matrix(GF(2),
            [[1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
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
        elif n == 24:
            G = matrix(GF(2),
            [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1],
             [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0],
             [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1],
             [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0],
             [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1],
             [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1],
             [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1],
             [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]])
        elif n == 11:
            G = matrix(GF(3),
            [[2, 0, 1, 2, 1, 1, 0, 0, 0, 0, 0],
             [0, 2, 0, 1, 2, 1, 1, 0, 0, 0, 0],
             [0, 0, 2, 0, 1, 2, 1, 1, 0, 0, 0],
             [0, 0, 0, 2, 0, 1, 2, 1, 1, 0, 0],
             [0, 0, 0, 0, 2, 0, 1, 2, 1, 1, 0],
             [0, 0, 0, 0, 0, 2, 0, 1, 2, 1, 1]])
        else:
            G = matrix(GF(3),
            [[1, 0, 0, 0, 0, 0, 2, 0, 1, 2, 1, 2],
             [0, 1, 0, 0, 0, 0, 1, 2, 2, 2, 1, 0],
             [0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1],
             [0, 0, 0, 1, 0, 0, 1, 1, 0, 2, 2, 2],
             [0, 0, 0, 0, 1, 0, 2, 1, 2, 2, 0, 1],
             [0, 0, 0, 0, 0, 1, 0, 2, 1, 2, 2, 1]])
        return G

    def parity_check_matrix(self):
        r"""
        Return the parity check matrix of ``self``.

        The parity check matrix of a linear code `C` corresponds to the
        generator matrix of the dual code of `C`.

        Parity check matrices of all Golay codes are known, and are thus returned
        by this method without performing any computation.

        EXAMPLES::

            sage: C = codes.GolayCode(GF(3), False)
            sage: C.parity_check_matrix()
            [1 0 0 0 0 1 2 2 2 1 0]
            [0 1 0 0 0 0 1 2 2 2 1]
            [0 0 1 0 0 2 1 2 0 1 2]
            [0 0 0 1 0 1 1 0 1 1 1]
            [0 0 0 0 1 2 2 2 1 0 1]
        """
        n = self.length()
        if n == 23:
            H = matrix(GF(2),
                [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0],
                 [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1],
                 [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0],
                 [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1],
                 [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1],
                 [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1],
                 [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1]])
        elif n == 11:
            H = matrix(GF(3),
                [[1, 0, 0, 0, 0, 1, 2, 2, 2, 1, 0],
                 [0, 1, 0, 0, 0, 0, 1, 2, 2, 2, 1],
                 [0, 0, 1, 0, 0, 2, 1, 2, 0, 1, 2],
                 [0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1],
                 [0, 0, 0, 0, 1, 2, 2, 2, 1, 0, 1]])
        else:
            H = self.generator_matrix()
        return H




####################### registration ###############################

GolayCode._registered_encoders["GeneratorMatrix"] = LinearCodeGeneratorMatrixEncoder
