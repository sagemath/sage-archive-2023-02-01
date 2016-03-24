r"""
BCH Code

Let `F = GF(q)` and `\Phi` be the splitting field of
`x^{n} - 1` over`F`, with `n` a positive integer.
Let `\alpha` be an element of multiplicative order `n` in `\Phi`.

A BCH code consists of all codewords `c(x) \in F_{n}[x]` such that
`c(\alpha^{a}) = 0`, for `a = b, b + l, b + 2\times l, \dots, b + (\delta - 2) \times l`,
with `b`, `\delta` integers such that `b < \delta` and `0 < \delta \leq n`.
"""

#*****************************************************************************
#       Copyright (C) 2016 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from linear_code import AbstractLinearCode
from cyclic_code import CyclicCode
from grs import GeneralizedReedSolomonCode
from encoder import Encoder
from decoder import Decoder, DecodingError
from sage.modules.free_module_element import vector
from sage.misc.misc_c import prod
from sage.rings.integer import Integer
from sage.rings.ring import Field
from copy import copy

class BCHCode(CyclicCode):
    r"""
    Representation of a BCH code seen as a cyclic code.

    INPUT:

    - ``base_field`` -- the base field for this code

    - ``length`` -- the length of the code

    - ``starting_point`` -- the first element to add in the defining set

    - ``delta`` -- the ending point for the elements in the defining set

    - ``jump_size`` -- (default: ``1``) the jump size between two elements of the defining set

    - ``b`` -- (default: ``0``) is exactly the same as ``starting_point``. It is only here
      for retro-compatibility purposes with the old signature of `BCHCode` and will be removed soon.

    EXAMPLES::

        sage: C = codes.BCHCode(GF(2), 15, 1, 7)
        sage: C
        [15, 5] BCH Code over Finite Field of size 2 with x^10 + x^8 + x^5 + x^4 + x^2 + x + 1 as generator polynomial

        sage: C = codes.BCHCode(GF(2), 15, 1, 4, 3)
        sage: C
        [15, 7] BCH Code over Finite Field of size 2 with x^8 + x^7 + x^5 + x^4 + x^3 + x + 1
        as generator polynomial
    """

    def __init__(self, base_field, length, starting_point, delta, jump_size = 1, b = 0):
        """
        TESTS:

        ``delta`` must be between 2 and ``length`` (inclusive), otherwise an exception
        will be raised::

            sage: C = codes.BCHCode(GF(2), 15, 1, 16)
            Traceback (most recent call last):
            ...
            ValueError: delta must belong to [2, n]
        """
        if not (delta <= length and delta > 1):
            raise ValueError("delta must belong to [2, n]")
        if isinstance(base_field, (Integer, int)) and isinstance(starting_point, Field):
            from sage.misc.superseded import deprecation
            deprecation(42042, "codes.BCHCode(n, delta, F, b=0) is now deprecated. Please use the new signature instead.")
            delta = copy(length)
            length = copy(base_field)
            F = copy(base_field)
        if not isinstance(base_field, Field):
            raise ValueError("base_field has to be a finite field")
        elif not base_field.is_finite():
            raise ValueError("base_field has to be a finite field")

        D = []
        point = copy(starting_point)
        for i in range(0, delta - 1):
            D.append(point)
            point = (point + jump_size) % length

        try:
            super(BCHCode, self).__init__(field = base_field, length = length, D = D)
        except ValueError, e:
            raise e
        self._default_decoder_name = "UnderlyingGRS"
        self._jump_size = jump_size
        self._starting_point = b
        self._delta = delta

    def __eq__(self, other):
        r"""
        Tests equality between BCH Code objects.

        EXAMPLES::

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: C1 = codes.BCHCode(F, n, 1, 2)
            sage: C2 = codes.BCHCode(F, n, 1, 2)
            sage: C1 == C2
            True
        """
        return isinstance(other, BCHCode) \
                and self.length() == other.length() \
                and self.jump_size() == other.jump_size() \
                and self.starting_point() == other.starting_point() \

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 1, 7)
            sage: C
            [15, 5] BCH Code over Finite Field of size 2 with x^10 + x^8 + x^5 + x^4 + x^2 + x + 1 as generator polynomial
        """
        return "[%s, %s] BCH Code over %s with %s as generator polynomial"\
                % (self.length(), self.dimension(),\
                self.base_field(), self.generator_polynomial())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 1, 7)
            sage: latex(C)
            [15, 5] \textnormal{ BCH Code over } \Bold{F}_{2} \textnormal{ with } x^{10} + x^{8} + x^{5} + x^{4} + x^{2} + x + 1 \textnormal{ as generator polynomial}
        """
        return "[%s, %s] \\textnormal{ BCH Code over } %s \\textnormal{ with } %s \\textnormal{ as generator polynomial}"\
                % (self.length(), self.dimension(),\
                self.base_field()._latex_(), self.generator_polynomial()._latex_())

    def jump_size(self):
        r"""
        Returns the jump size between two consecutive elements of the defining set of ``self``.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 1, 4, 2)
            sage: C.jump_size()
            2
        """
        return self._jump_size

    def starting_point(self):
        r"""
        Returns the starting point which was used to compute the elements in
        the defining set of ``self``.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 1, 4, 2)
            sage: C.starting_point
            1
        """
        return self._starting_point

    def bch_to_grs(self):
        r"""
        Returns the underlying GRS code from which ``self`` was derived.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 1, 2, 2)
            sage: C.bch_to_grs()
            [15, 13, 3] Generalized Reed-Solomon Code over Finite Field in b of size 2^4
        """
        l = self.jump_size()
        alpha = self.root_of_unity()
        b = self.starting_point()
        n = self.length()

        grs_dim = n - self.bch_bound(arithmetic = True) + 1
        evals = []
        pcm = []
        for i in range(1, n + 1):
            evals.append(alpha ** (l * (i - 1)))
            pcm.append(alpha ** (b * (i - 1)))


        multipliers_product = [1/prod([evals[i] - evals[h] for h in range(len(evals)) if h != i])
                    for i in range(len(evals))]
        column_multipliers = [multipliers_product[i]/pcm[i] for i in range(n)]

        return GeneralizedReedSolomonCode(evals, grs_dim, column_multipliers)
