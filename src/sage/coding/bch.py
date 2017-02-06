r"""
BCH Code

Let `F = GF(q)` and `\Phi` be the splitting field of `x^{n} - 1` over `F`,
with `n` a positive integer. Let also `\alpha` be an element of multiplicative
order `n` in `\Phi`. Finally, let `b, \delta, \ell` be integers such that
`0 \le b \leq n`, `2 \le \delta \leq n` and `\alpha^\ell` generates the
multiplicative group `\Phi^{\times}`.

A BCH code over `F` with designed distance `\delta` is a cyclic code whose
codewords `c(x) \in F[x]` satisfy `c(\alpha^{a}) = 0`, for all integers `a` in
the arithmetic sequence
`b, b + \ell, b + 2 \times \ell, \dots, b + (\delta - 2) \times \ell`.

TESTS:

This class uses the following experimental feature:
:class:`sage.coding.relative_finite_field_extension.RelativeFiniteFieldExtension`.
This test block is here only to trigger the experimental warning so it does not
interferes with doctests::

    sage: from sage.coding.relative_finite_field_extension import *
    sage: Fqm.<aa> = GF(16)
    sage: Fq.<a> = GF(4)
    sage: RelativeFiniteFieldExtension(Fqm, Fq)
    doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
    See http://trac.sagemath.org/20284 for details.
    Relative field extension between Finite Field in aa of size 2^4 and Finite Field in a of size 2^2
"""

# *****************************************************************************
#       Copyright (C) 2016 David Lucas <david.lucas@inria.fr>
#                     2017 Julien Lavauzelle <julien.lavauzelle@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from .cyclic_code import CyclicCode
from .grs import GeneralizedReedSolomonCode
from .decoder import Decoder, DecodingError
from sage.modules.free_module_element import vector
from sage.misc.misc_c import prod
from sage.categories.fields import Fields
from sage.rings.integer_ring import ZZ
from sage.arith.all import gcd
from sage.rings.all import Zmod
from copy import copy


class BCHCode(CyclicCode):
    r"""
    Representation of a BCH code seen as a cyclic code.

    INPUT:

    - ``base_field`` -- the base field for this code

    - ``length`` -- the length of the code

    - ``designed_distance`` -- the resulting minimum distance of the code

    - ``primitive_root`` -- (default: ``None``) the primitive root to use when
      creating the set of roots for the generating polynomial over the
      splitting field. It has to be of multiplicative order ``length`` over
      this field. If the splitting field is not ``field``, it also has to be a
      polynomial in ``zx``, where ``x`` is the degree of the extension field.
      For instance, over ``GF(16)``, it has to be a polynomial in ``z4``.

    - ``offset`` -- (default: ``0``) the first element in the defining set

    - ``jump_size`` -- (default: ``1``) the jump size between two elements of
      the defining set. It must be coprime with the multiplicative order of
      ``primitive_root``.

    - ``b`` -- (default: ``0``) is exactly the same as ``offset``. It is only
      here for retro-compatibility purposes with the old signature of `BCHCode`
      and will be removed soon.

    EXAMPLES::

        sage: C = codes.BCHCode(GF(2), 15, 7, offset = 1)
        sage: C
        [15, 5] BCH Code over GF(2) with designed distance 7
        sage: C.generator_polynomial()
        x^10 + x^8 + x^5 + x^4 + x^2 + x + 1

        sage: C = codes.BCHCode(GF(2), 15, 4, offset = 1, jump_size = 8)
        sage: C
        [15, 7] BCH Code over GF(2) with designed distance 4
        sage: C.generator_polynomial()
        x^8 + x^7 + x^6 + x^4 + 1
    """

    def __init__(self, base_field, length, designed_distance,
                 primitive_root=None, offset=0, jump_size=1, b=0):
        """
        TESTS:

        ``designed_distance`` must be between 2 and ``length`` (inclusive),
        otherwise an exception is raised::

            sage: C = codes.BCHCode(GF(2), 15, 16)
            Traceback (most recent call last):
            ...
            ValueError: designed_distance must belong to [2, n]
        """
        if not (designed_distance <= length and designed_distance > 1):
            raise ValueError("designed_distance must belong to [2, n]")
        if base_field in ZZ and designed_distance in Fields:
            from sage.misc.superseded import deprecation
            deprecation(20335, "codes.BCHCode(n, designed_distance, F, b=0) is now deprecated. Please use the new signature instead.")
            (length, designed_distance, base_field) = (base_field, length, designed_distance)
            offset = b
        if base_field not in Fields or not base_field.is_finite():
            raise ValueError("base_field has to be a finite field")

        q = base_field.cardinality()
        s = Zmod(length)(q).multiplicative_order()
        if gcd(jump_size, q ** s - 1) != 1:
            raise ValueError("jump_size must be coprime with the order of "
                             "the multicative group of the splitting field")

        D = [(offset + jump_size * i) % length
             for i in range(designed_distance - 1)]

        try:
            super(BCHCode, self).__init__(field=base_field, length=length,
                                          D=D, primitive_root=primitive_root)
        except ValueError, e:
            raise e
        self._default_decoder_name = "UnderlyingGRS"
        self._jump_size = jump_size
        self._offset = offset
        self._designed_distance = designed_distance

    def __eq__(self, other):
        r"""
        Tests equality between BCH Code objects.

        EXAMPLES::

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: C1 = codes.BCHCode(F, n, 2)
            sage: C2 = codes.BCHCode(F, n, 2)
            sage: C1 == C2
            True
        """
        return (isinstance(other, BCHCode) and
                self.length() == other.length() and
                self.jump_size() == other.jump_size() and
                self.offset() == self.offset() and
                self.primitive_root() == self.primitive_root())

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 7, offset = 1)
            sage: C
            [15, 5] BCH Code over GF(2) with designed distance 7
        """
        return "[%s, %s] BCH Code over GF(%s) with designed distance %d"\
                % (self.length(), self.dimension(),
                   self.base_field().cardinality(), self.designed_distance())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 7, offset = 1)
            sage: latex(C)
            [15, 5] \textnormal{ BCH Code over } \Bold{F}_{2} \textnormal{ with designed distance } 7
        """
        return ("[%s, %s] \\textnormal{ BCH Code over } %s \\textnormal{ with designed distance } %s"
                % (self.length(), self.dimension(),
                self.base_field()._latex_(), self.designed_distance()))

    def jump_size(self):
        r"""
        Returns the jump size between two consecutive elements of the defining
        set of ``self``.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 4, jump_size = 2)
            sage: C.jump_size()
            2
        """
        return self._jump_size

    def offset(self):
        r"""
        Returns the offset which was used to compute the elements in
        the defining set of ``self``.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 4, offset = 1)
            sage: C.offset()
            1
        """
        return self._offset

    def designed_distance(self):
        r"""
        Returns the designed distance of ``self``.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 4, offset = 1)
            sage: C.designed_distance()
            4
        """
        return self._designed_distance

    def bch_to_grs(self):
        r"""
        Returns the underlying GRS code from which ``self`` was derived.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 3)
            sage: C.bch_to_grs()
            [15, 13, 3] Generalized Reed-Solomon Code over GF(16)
        """
        l = self.jump_size()
        b = self.offset()
        n = self.length()
        designed_distance = self.designed_distance()
        grs_dim = n - designed_distance + 1

        alpha = self.primitive_root()
        alpha_l = alpha ** l
        alpha_b = alpha ** b
        evals = [alpha_l ** i for i in range(n)]
        pcm = [alpha_b ** i for i in range(n)]

        multipliers_product = [1/prod([evals[i] - evals[h] for h in range(n) if h != i]) for i in range(n)]
        column_multipliers = [multipliers_product[i]/pcm[i] for i in range(n)]

        return GeneralizedReedSolomonCode(evals, grs_dim, column_multipliers)


class BCHUnderlyingGRSDecoder(Decoder):
    r"""
    A decoder which decodes through the underlying
    :class:`sage.coding.grs.GeneralizedReedSolomonCode` code of the provided
    BCH code.

    INPUT:

    - ``code`` -- The associated code of this decoder.

    - ``grs_decoder`` -- The string name of the decoder to use over the
      underlying GRS code

    - ``**kwargs`` -- All extra arguments are forwarded to the GRS decoder
    """

    def __init__(self, code, grs_decoder="KeyEquationSyndrome", **kwargs):
        r"""

        EXAMPLES::

            sage: C = codes.BCHCode(GF(4, 'a'), 15, 3, jump_size=2, offset=1)
            sage: D = codes.decoders.BCHUnderlyingGRSDecoder(C)
            sage: D
            Decoder through the underlying GRS code of [15, 11] BCH Code over GF(4) with designed distance 3
        """

        self._grs_code = code.bch_to_grs()
        self._grs_decoder = self._grs_code.decoder(grs_decoder, **kwargs)
        self._decoder_type = copy(self._decoder_type)
        self._decoder_type.remove("dynamic")
        self._decoder_type = self._grs_decoder.decoder_type()
        super(BCHUnderlyingGRSDecoder, self).__init__(
            code, code.ambient_space(), "Vector")

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(4, 'a'), 15, 3, jump_size=2, offset=1)
            sage: D = codes.decoders.BCHUnderlyingGRSDecoder(C)
            sage: D
            Decoder through the underlying GRS code of [15, 11] BCH Code over GF(4) with designed distance 3
        """
        return "Decoder through the underlying GRS code of %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(4, 'a'), 15, 3, jump_size=2, offset=1)
            sage: D = codes.decoders.BCHUnderlyingGRSDecoder(C)
            sage: latex(D)
            \textnormal{Decoder through the underlying GRS code of } [15, 11] \textnormal{ BCH Code over } \Bold{F}_{2^{2}} \textnormal{ with designed distance } 3
        """
        return ("\\textnormal{Decoder through the underlying GRS code of } %s"
                % self.code()._latex_())

    def grs_code(self):
        r"""
        Returns the underlying GRS code of
        :meth:`sage.coding.decoder.Decoder.code`.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 3)
            sage: D = codes.decoders.BCHUnderlyingGRSDecoder(C)
            sage: D.grs_code()
            [15, 13, 3] Generalized Reed-Solomon Code over GF(16)
        """
        return self._grs_code

    def grs_decoder(self):
        r"""
        Returns the decoder used to decode words of :meth:`grs_code`.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(4, 'a'), 15, 3, jump_size=2, offset=1)
            sage: D = codes.decoders.BCHUnderlyingGRSDecoder(C)
            sage: D.grs_decoder()
            Key equation decoder for [15, 13, 3] Generalized Reed-Solomon Code over GF(16)
        """
        return self._grs_decoder

    def bch_word_to_grs(self, c):
        r"""
        Returns ``c`` converted as a codeword of :meth:`grs_code`.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(2), 15, 3)
            sage: D = codes.decoders.BCHUnderlyingGRSDecoder(C)
            sage: c = C.random_element()
            sage: y = D.bch_word_to_grs(c)
            sage: y.parent()
            Vector space of dimension 15 over Finite Field in z4 of size 2^4
            sage: y in D.grs_code()
            True
        """
        mapping = self.code().field_embedding().embedding()
        a = [mapping(i) for i in c]
        return vector(a)

    def grs_word_to_bch(self, c):
        r"""
        Returns ``c`` converted as a codeword of
        :meth:`sage.coding.decoder.Decoder.code`.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(4, 'a'), 15, 3, jump_size=2, offset=1)
            sage: D = codes.decoders.BCHUnderlyingGRSDecoder(C)
            sage: Cgrs = D.grs_code()
            sage: Fgrs = Cgrs.base_field()
            sage: b = Fgrs.gen()
            sage: c = vector(Fgrs, [0, b^2 + b, 1, b^2 + b, 0, 1, 1, 1, b^2 + b, 0, 0, b^2 + b + 1, b^2 + b, 0, 1])
            sage: D.grs_word_to_bch(c)
            (0, a, 1, a, 0, 1, 1, 1, a, 0, 0, a + 1, a, 0, 1)
        """
        C = self.code()
        FE = C.field_embedding()
        a = map(FE.cast_into_relative_field, c)
        return vector(a)

    def decode_to_code(self, y):
        r"""
        Decodes ``y`` to a codeword in :meth:`sage.coding.decoder.Decoder.code`.

        EXAMPLES::

            sage: F = GF(4, 'a')
            sage: a = F.gen()
            sage: C = codes.BCHCode(F, 15, 3, jump_size=2, offset=1)
            sage: D = codes.decoders.BCHUnderlyingGRSDecoder(C)
            sage: y = vector(F, [a, a + 1, 1, a + 1, 1, a, a + 1, a + 1, 0, 1, a + 1, 1, 1, 1, a])
            sage: D.decode_to_code(y)
            (a, a + 1, 1, a + 1, 1, a, a + 1, a + 1, 0, 1, a + 1, 1, 1, 1, a)
            sage: D.decode_to_code(y) in C
            True
        """
        D = self.grs_decoder()
        ygrs = self.bch_word_to_grs(y)
        try:
            cgrs = D.decode_to_code(ygrs)
        except DecodingError, e:
            raise e
        if "list-decoder" in D.decoder_type():
            l = []
            for c in cgrs:
                l.append(self.grs_word_to_bch(c))
            return l
        return self.grs_word_to_bch(cgrs)

    def decoding_radius(self):
        r"""
        Returns maximal number of errors that ``self`` can decode.

        EXAMPLES::

            sage: C = codes.BCHCode(GF(4, 'a'), 15, 3, jump_size=2, offset=1)
            sage: D = codes.decoders.BCHUnderlyingGRSDecoder(C)
            sage: D.decoding_radius()
            1
        """
        return (self.code().bch_bound(arithmetic=True)[0] - 1) // 2


####################### registration ###############################

BCHCode._registered_decoders["UnderlyingGRS"] = BCHUnderlyingGRSDecoder
BCHUnderlyingGRSDecoder._decoder_type = {"dynamic"}
