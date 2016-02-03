r"""
Subfield subcode

Let `C` be a `[n, k]` code over `\GF(q^t)`.
Let `Cs = \{c \in C | \forall i, c_i \in \GF(q)\}`, `c_i` being the `i`-th
coordinate of `c`.

`Cs` is called the subfield subcode of `C` over `\GF(q)`
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

from linear_code import (AbstractLinearCode,
                         LinearCodeSyndromeDecoder,
                         LinearCodeNearestNeighborDecoder)
from encoder import Encoder
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.functions.all import log

class SubfieldSubcode(AbstractLinearCode):
    r"""
    Representation of a subfield subcode

    EXAMPLES::

        sage: C = codes.RandomLinearCode(7, 3, GF(8, 'a'))
        sage: Cs = codes.SubfieldSubcode(C, 4)
        sage: Cs
        Subfield subcode of order 4 coming from Linear code of length 7, dimension 3 over Finite Field in a of size 2^3
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, original_code, subfield_order):
        r"""
        TESTS:

        ``subfield_order`` has to divide the order of ``original_code``'s base field,
        otherwise an error is raised::

            sage: C = codes.RandomLinearCode(7, 3, GF(8, 'a'))
            sage: Cs = codes.SubfieldSubcode(C, 3)
            Traceback (most recent call last):
            ...
            ValueError: subfield_order must divide the order of original_code's base field

        """
        if not isinstance(original_code, AbstractLinearCode):
            raise ValueError("original_code must be a linear code")
        if not isinstance(subfield_order, (int, Integer)):
            raise ValueError("subfield_order must be a Python int or a Sage integer")
        subfield_order = Integer(subfield_order)
        if not subfield_order.divides(original_code.base_field().order()):
            raise ValueError("subfield_order must divide the order of original_code's base field")
        self._original_code = original_code
        self._subfield_order = subfield_order
        if subfield_order.is_prime():
            F = GF(subfield_order)
        else:
            F = GF(subfield_order, 'x')
        super(SubfieldSubcode, self).__init__(F, original_code.length(), "ParityCheck", "Syndrome")

    def __eq__(self, other):
        r"""
        Tests equality between Subfield Subcode objects.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(7, 3, GF(8, 'a'))
            sage: Cs1 = codes.SubfieldSubcode(C, 4)
            sage: Cs2 = codes.SubfieldSubcode(C, 4)
            sage: Cs1 == Cs2
            True
        """
        return isinstance(other, SubfieldSubcode) \
                and self.original_code() == other.original_code()\
                and self.base_field().order() == other.base_field.order()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(7, 3, GF(8, 'a'))
            sage: Cs = codes.SubfieldSubcode(C, 4)
            sage: Cs
            Subfield subcode of order 4 coming from Linear code of length 7, dimension 3 over Finite Field in a of size 2^3
        """
        return "Subfield subcode of order %s coming from %s"\
                % (self.base_field().order(), self.original_code())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(7, 3, GF(8, 'a'))
            sage: Cs = codes.SubfieldSubcode(C, 4)
            sage: Cs

        """
        return "\\textnormal{Subfield subcode of order %s coming from }%s"\
                % (self.base_field().order(), self.original_code())

    def dimension(self):
        r"""
        Returns the dimension of ``self``.

        """
        return self.generator_matrix().nrows()

    def dimension_upper_bound(self):
        r"""
        Returns an upper bound for the dimension of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(7, 3, GF(8, 'a'))
            sage: Cs = codes.SubfieldSubcode(C, 4)
            sage: Cs.dimension_upper_bound()
            3
        """
        return self.original_code().dimension()

    def dimension_lower_bound(self):
        r"""
        Returns a lower bound for the dimension of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(7, 3, GF(8, 'a'))
            sage: Cs = codes.SubfieldSubcode(C, 4)
            sage: Cs.dimension_lower_bound()
            1
        """
        C = self._original_code()
        n = C.length()
        k = C.dimension()
        F = C.base_field()
        t = log(F.order() // self.base_field().order(), F.characteristic())
        return n - t*(n-k)

    def original_code(self):
        r"""
        Returns the original code of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(7, 3, GF(8, 'a'))
            sage: Cs = codes.SubfieldSubcode(C, 4)
            sage: Cs.original_code()
            Linear code of length 7, dimension 3 over Finite Field in a of size 2^3
        """
        return self._original_code

    @cached_method
    def parity_check_matrix(self):
        r"""
        Returns a parity check matrix of ``self``.

        """
        raise NotImplementedError








#Purely TEMPORARY, will be integrated with trac #19930

class SubfieldSubcodeParityCheckEncoder(Encoder):
    r"""
    Encoder based on :meth:`parity_check_matrix` for Linear codes.

    It constructs the generator matrix through the parity check matrix.

    INPUT:

    - ``code`` -- The associated code of this encoder.
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeParityCheckEncoder(C)
            sage: E
            Parity check matrix-based encoder for the Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        super(LinearCodeParityCheckEncoder, self).__init__(code)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeParityCheckEncoder(C)
            sage: E
            Parity check matrix-based encoder for the Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        return "Parity check matrix-based encoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeParityCheckEncoder(C)
            sage: latex(E)
            \textnormal{Parity check matrix-based encoder for the }[7, 4]\textnormal{ Linear code over }\Bold{F}_{2}
        """
        return "\\textnormal{Parity check matrix-based encoder for the }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeParityCheckEncoder(C)
            sage: E.generator_matrix()
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1]
        """
        return self.code().parity_check_matrix().right_kernel_matrix()


####################### registration ###############################

SubfieldSubcode._registered_encoders["ParityCheck"] = SubfieldSubcodeParityCheckEncoder
SubfieldSubcode._registered_decoders["Syndrome"] = LinearCodeSyndromeDecoder
SubfieldSubcode._registered_decoders["NearestNeighbor"] = LinearCodeNearestNeighborDecoder
