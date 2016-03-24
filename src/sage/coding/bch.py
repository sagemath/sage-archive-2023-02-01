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

class BCHCode(CyclicCode):
    r"""
    Representation of a BCH code seen as a cyclic code.

    INPUT:

    - ``F`` -- the base field for this code

    - ``n`` -- the length of the code

    - ``b`` -- the starting point for the elements in the defining set

    - ``delta`` -- the ending point for the elements in the defining set

    - ``l`` -- (default: ``1``) the jump size between two elements of the defining set

    EXAMPLES::

        sage: C = codes.BCHCode(GF(2), 15, 1, 7)
        sage: C
        [15, 5] BCH Code over Finite Field of size 2 with x^10 + x^8 + x^5 + x^4 + x^2 + x + 1 as generator polynomial

        sage: C = codes.BCHCode(GF(2), 15, 1, 4, 3)
        sage: C
        [15, 7] BCH Code over Finite Field of size 2 with x^8 + x^7 + x^5 + x^4 + x^3 + x + 1
        as generator polynomial
    """

    def __init__(self, F, n, b, delta, l = 1):

        if not (delta <= n and delta > 1):
            raise ValueError("delta must belong to [2, n]")

        D = []
        d = b
        for i in range(0, delta - 1):
            D.append(d)
            d = (d + l) % n

        try:
            super(BCHCode, self).__init__(field = F, length = n, D = D)
        except ValueError, e:
            raise e
        self._default_decoder_name = "UnderlyingGRS"
        self._jump_size = l
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

    def __ne__(self, other):
        r"""
        Tests inequality of BCH Code objects.

        EXAMPLES::

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: C1 = codes.BCHCode(F, 15, 1, 2)
            sage: C2 = codes.BCHCode(F, 13, 1, 2)
            sage: C1 != C2
            True
        """
        return not self.__eq__(other)

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
