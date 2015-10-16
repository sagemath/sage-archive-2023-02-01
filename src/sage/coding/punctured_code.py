r"""
Punctured code

Let `C` be a linear code. Let `C_i` be the set of all words of `C` with the `i`-th coordinate being
removed. `C_i` is the punctured code of `C` on the `i`-th position.
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

from linear_code import AbstractLinearCode
from encoder import Encoder
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.modules.free_module import VectorSpace

def puncture(v, points, code):
    r"""
    Returns v punctured as the positions listed in ``points``.

    INPUT:

    - ``v`` -- a vector

    - ``points`` -- a list of integers

    - ``code`` -- the code in which ``v`` lives

    EXAMPLES::

        sage: C = codes.RandomLinearCode(11, 5, GF(7))
        sage: Cp = codes.PuncturedCode(C, 3)
        sage: v = vector(GF(7), (2,3,0,2,1,5,1,5,6,5,3))
        sage: sage.coding.punctured_code.puncture(v, Cp.punctured_positions(), Cp)
        (2, 3, 0, 1, 5, 1, 5, 6, 5, 3)
    """
    S = code.ambient_space()
    vl = v.list()
    v_final = []
    start = 0
    for i in points:
        v_final += vl[start:i]
        start = i + 1
    v_final += vl[start:len(vl)]
    return S(v_final)


class PuncturedCode(AbstractLinearCode):
    r"""
    Representation of a punctured code.

    - ``C`` -- A linear code

    - ``positions`` -- the positions where ``C`` will be punctured. It can be either an Sage integer
      or a Python int if you need to puncture only one position, or a list of positions to puncture.
      If the same position is passed several times, it will be considered only once.

    EXAMPLES::

        sage: C = codes.RandomLinearCode(11, 5, GF(7))
        sage: Cp = codes.PuncturedCode(C, 3)
        sage: Cp
        Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3]

        sage: Cp = codes.PuncturedCode(C, [3, 5])
        sage: Cp
        Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3, 5]
    """
    _registered_encoders = {}

    def __init__(self, C, positions):
        r"""
        TESTS::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 13)
            Traceback (most recent call last):
            ...
            ValueError: Positions to puncture must be positive integers smaller than the length of the provided code
        """
        if not isinstance(positions, (Integer, int, tuple, list)):
            raise TypeError("positions must be either a Sage Integer, a Python int, a tuple or a list")
        if isinstance(positions, (Integer, int)):
            positions = [positions]
        if not isinstance(C, AbstractLinearCode):
            raise ValueError("Provided code must be a linear code")
        if not all (i in range(0, C.length()) for i in positions):
            raise ValueError("Positions to puncture must be positive integers smaller than the length of the provided code")
        unique_positions = set()
        for i in positions:
            unique_positions.add(i)
        positions = []
        for i in unique_positions:
            positions.append(i)
        super(PuncturedCode, self).__init__(C.base_ring(), C.length() - len(positions), \
                "PuncturedMatrix")
        positions.sort()
        self._original_code = C
        self._positions = positions

    def __eq__(self, other):
        r"""
        Tests equality between two Punctured codes.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp1 = codes.PuncturedCode(C, 2)
            sage: Cp2 = codes.PuncturedCode(C, 2)
            sage: Cp1 == Cp2
            True
        """
        return isinstance(other, PuncturedCode) \
                and self.punctured_positions() == other.punctured_positions() \
                and self.original_code() == other.original_code()

    def __ne__(self, other):
        r"""
        Tests inequality between two Punctured codes.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp1 = codes.PuncturedCode(C, 2)
            sage: Cp2 = codes.PuncturedCode(C, 3)
            sage: Cp1 != Cp2
            True
        """
        return not self.__eq__(other)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp
            Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3]
        """
        return "Punctured code coming from %s punctured on position(s) %s"\
                % (self.original_code(), self.punctured_positions())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: latex(Cp)
            \textnormal{Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) } [3]
        """
        return "\\textnormal{Punctured code coming from %s punctured on position(s) } %s"\
                % (self.original_code(), self.punctured_positions())

    def punctured_positions(self):
        r"""
        Returns the list of positions which were punctured on the original code

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp.punctured_positions()
            [3]
        """
        return self._positions

    def original_code(self):
        r"""
        Returns the linear code which was punctured to get ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp.original_code()
            Linear code of length 11, dimension 5 over Finite Field of size 7
        """
        return self._original_code

    def dimension(self):
        r"""
        Returns the dimension of ``self``

        EXAMPLES::

            sage: set_random_seed(42)
            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp.dimension()
            5
        """
        if hasattr(self, '_dimension'):
            return self._dimension
        self._dimension = self.generator_matrix().rank()
        return self._dimension

    def random_element(self, *args, **kwds):
        r"""
        Returns a random codeword of ``self``.

        This methods avoids computation of ``self``'s :meth:`sage.coding.linear_code.generator_matrix`
        by puncturing the result of ``self``'s :meth:`original_code`'s
        :meth:`sage.coding.linear_code.random_element`.

        INPUT:

        - ``agrs``, ``kwds`` - extra positional arguments passed to
          :meth:`sage.modules.free_module.random_element`.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: set_random_seed(10)
            sage: Cp.random_element()
            (2, 0, 1, 3, 3, 3, 2, 6, 0, 5)
        """
        C_original = self.original_code()
        m = (C_original.base_ring() ** C_original.dimension()).random_element()
        c = C_original.encode(m)
        return puncture(c, self.punctured_positions(), self)

    def encode(self, m, original_encode=False, encoder_name=None, **kwargs):
        r"""
        Transforms an element of the message space into an element of the code.

        INPUT:

        - ``m`` -- a vector of the message space of the code.

        - ``original_encode`` -- (default: ``False``) if this is set to ``True``,
          ``m`` will be encoded using an Encoder of ``self``'s :meth:`original_code`.
          This allow to avoid the computation of a generator matrix for ``self``.

        - ``encoder_name`` -- (default: ``None``) Name of the encoder which will be used
          to encode ``word``. The default encoder of ``self`` will be used if
          default value is kept

        OUTPUT:

        - an element of ``self``

        EXAMPLES::

           sage: M = matrix(GF(7), [[1, 0, 0, 0, 3, 4, 6], [0, 1, 0, 6, 1, 6, 4], [0, 0, 1, 5, 2, 2, 4]])
           sage: C_original = LinearCode(M)
           sage: Cp = codes.PuncturedCode(C_original, 2)
           sage: m = vector(GF(7), [1, 3, 5])
           sage: Cp.encode(m)
            (1, 3, 5, 5, 0, 2)
        """
        if original_encode:
            c = self.original_code().encode(m, encoder_name, **kwargs)
            return puncture(c, self.punctured_positions, self)
        return self.encoder(encoder_name, **kwargs).encode(m)










class PuncturedCodePuncturedMatrixEncoder(Encoder):
    r"""
    Encoder using original code generator matrix to compute the punctured code's one.

    INPUT:

    - ``code`` -- The associated code of this encoder.
    """

    def __init__(self, code):
        r"""
        EXAMPLES::
            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: E
            Punctured matrix-based encoder for the Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3]
        """
        super(PuncturedCodePuncturedMatrixEncoder, self).__init__(code)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: E
            Punctured matrix-based encoder for the Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3]
        """
        return "Punctured matrix-based encoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: latex(E)
            \textnormal{Punctured matrix-based encoder for the }\textnormal{Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) } [3]
        """
        return "\\textnormal{Punctured matrix-based encoder for the }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.

        EXAMPLES::

            sage: set_random_seed(10)
            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: E.generator_matrix()
            [1 0 0 0 0 5 2 6 0 6]
            [0 1 0 0 0 5 2 2 1 1]
            [0 0 1 0 0 6 2 4 0 4]
            [0 0 0 1 0 0 6 3 3 3]
            [0 0 0 0 1 0 1 3 4 3]
        """
        C = self.code().original_code()
        pos = self.code().punctured_positions()
        M = C.generator_matrix()
        G = M.delete_columns(pos)
        G = G.echelon_form()
        delete = []
        cpt = 0
        for i in G.rows():
            if i.is_zero():
                delete.append(cpt)
            cpt += 1
        return G.delete_rows(delete)

####################### registration ###############################

PuncturedCode._registered_encoders["PuncturedMatrix"] = PuncturedCodePuncturedMatrixEncoder
