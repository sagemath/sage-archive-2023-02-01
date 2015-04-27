r"""
Linear Code Encoders
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

from encoder import Encoder
from sage.misc.cachefunc import cached_method

class EncoderLinearCodeGeneratorMatrix(Encoder):
    r"""
    Encoder based on generator_matrix for Linear codes.

    The only purpose of this encoder is to include the existing code
    into the new Encoder and Decoder structure.

    INPUT:

    - ``code`` -- The associated code of this encoder.
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = EncoderLinearCodeGeneratorMatrix(C)
            sage: E
            Generator matrix-based encoder for the Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        super(EncoderLinearCodeGeneratorMatrix, self).__init__(code)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = EncoderLinearCodeGeneratorMatrix(C)
            sage: E
            Generator matrix-based encoder for the Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        return "Generator matrix-based encoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = EncoderLinearCodeGeneratorMatrix(C)
            sage: latex(E)
            \textnormal{Generator matrix-based encoder for the }[7, 4]\textnormal{ Linear code over }\Bold{F}_{2}
        """
        return "\\textnormal{Generator matrix-based encoder for the }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = EncoderLinearCodeGeneratorMatrix(C)
            sage: E.generator_matrix()
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
        """
        if hasattr(self.code(), "_generator_matrix"):
            return self.code()._generator_matrix
        else:
            return self.code().generator_matrix()
