r"""
Encoder

Representation of a bijection between a message space and a code.
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

from sage.modules.free_module_element import vector
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject

class Encoder(SageObject):
    r"""
    Abstract top-class for :class:`Encoder` objects.

    Every encoder class should inherit from this abstract class.

    To implement an encoder, you need to:

    - inherit from :class:`Encoder`

    - call ``Encoder.__init__`` in the subclass constructor.
      Example: ``super(SubclassName, self).__init__(code)``.
      By doing that, your subclass will have its ``code`` parameter initialized.
      You need of course to complete the constructor by adding any additional parameter
      needed to describe properly the code defined in the subclass.

    Then, if the message space is a vector space, default implementations of :meth:`encode` and
    :meth:`unencode_nocheck` methods are provided. These implementations rely on :meth:`generator_matrix`
    which you need to override to use the default implementations.

    If the message space is not of the form `F^k`, where `F` is a finite field,
    you cannot have a generator matrix.
    In that case, you need to override :meth:`encode` and :meth:`unencode_nocheck`.

    Equality methods (``__eq__`` and ``__ne__``) might be useful for encoding in advanced
    codes constructions (like concatenated codes). If provided default implementation of
    these methods is not enough for your subclass, you are strongly encouraged to override
    them.

    As :class:`Encoder` is not designed to be instanciated, it does not have any representation
    methods. You should implement ``_repr_`` and ``_latex_`` methods in the sublclass.

    REFERENCES:

    .. [Nielsen] Johan S. R. Nielsen, (https://bitbucket.org/jsrn/codinglib/)
    """

    def __init__(self, code):
        r"""
        Initializes mandatory parameters for an :class:`Encoder` object.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every linear code. An abstract
        encoder object should never be created.

        INPUT:

        - ``code`` -- the associated code of ``self``

        EXAMPLES:

        We first create a new :class:`Encoder` subclass::

            sage: class EncoderExample(sage.coding.encoder.Encoder):
            ....:   def __init__(self, code):
            ....:       super(EncoderExample, self).__init__(code)

        We now create a member of our newly made class::

            sage: G = Matrix(GF(2), [[1, 0, 0, 1], [0, 1, 1, 1]])
            sage: C = LinearCode(G)
            sage: E = EncoderExample(C)

        We can check its parameters::

            sage: E.code()
            Linear code of length 4, dimension 2 over Finite Field of size 2
        """
        self._code = code

    def encode(self, word):
        r"""
        Transforms an element of the message space into a codeword.

        This is a default implementation which assumes that the message
        space of the encoder is `F^k`, where `F` is
        :meth:`sage.coding.linear_code.AbstractLinearCode.base_field`
        and ``k`` is :meth:`sage.coding.linear_code.AbstractLinearCode.dimension`.
        If this is not the case, this method should be overwritten by the subclass.

        INPUT:

        - ``word`` -- a vector of the message space of the code

        OUTPUT:

        - a vector of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (0, 1, 1, 0))
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E.encode(word)
            (1, 1, 0, 0, 1, 1, 0)

        If ``word`` is not in the message space of ``self``, it will return an exception::

            sage: word = random_vector(GF(7), 4)
            sage: E.encode(word)
            Traceback (most recent call last):
            ...
            ValueError: Vector to encode must be in a Vector space of dimension 4 over Finite Field of size 2
        """
        M = self.message_space()
        if word not in M:
            raise ValueError("Vector to encode must be in a %s" % M)
        return vector(word) * self.generator_matrix()

    def unencode(self, c, nocheck=False):
        r"""
        Returns the message corresponding to ``c``.

        INPUT:

        - ``c`` -- a vector of the same length as ``self`` over the
          base field of ``self``

        - ``nocheck`` -- (default: ``False``) checks if ``c`` is in ``self``. If this is set
          to ``True``, the return value of this method is not guaranteed to be correct.

        OUTPUT:

        - a vector of the message space of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: c = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E.unencode(c)
            (0, 1, 1, 0)
        """
        if nocheck == False:
            if c not in self.code():
                raise EncodingError("Given word is not in the code")
            else:
                return self.unencode_nocheck(c)
        else:
            return self.unencode_nocheck(c)

    @cached_method
    def _unencoder_matrix(self):
        r"""
        Finds an information set for the matrix ``G`` returned by :meth:`generator_matrix`,
        and returns the inverse of that submatrix of ``G``.

        AUTHORS:

            This function is taken from codinglib [Nielsen]_

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = C.encoder()
            sage: E._unencoder_matrix()
            [0 0 1 1]
            [0 1 0 1]
            [1 1 1 0]
            [0 1 1 1]
        """
        Gt = self.generator_matrix().matrix_from_columns(self.code().information_set())
        return Gt.inverse()

    def unencode_nocheck(self, c):
        r"""
        Returns the message corresponding to ``c``.

        When ``c`` is not a codeword, the output is unspecified.

        AUTHORS:

            This function is taken from codinglib [Nielsen]_

        INPUT:

        - ``c`` -- a vector of the same length as ``self`` over the
          base field of ``self``

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: c = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: c in C
            True
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E.unencode_nocheck(c)
            (0, 1, 1, 0)

        We take a vector that does not belong to C::

            sage: c = vector(GF(2), (1, 1, 0, 0, 1, 1, 1))
            sage: c in C
            False
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E.unencode_nocheck(c)
            (0, 1, 1, 0)
            sage: m = vector(GF(2), (0, 1, 1, 0))
            sage: c1 = E.encode(m)
            sage: c == c1
            False
        """
        U = self._unencoder_matrix()
        info_set = self.code().information_set()
        cc = vector( c[i] for i in info_set )
        return cc * U

    def code(self):
        r"""
        Returns the code in which :meth:`encode` has its output.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = C.encoder()
            sage: E.code()
            Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        return self._code

    def message_space(self):
        r"""
        Returns the ambient space of allowed input to :meth:`encode`.
        Note that :meth:`encode` is possibly a partial function over
        the ambient space.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = C.encoder()
            sage: E.message_space()
            Vector space of dimension 4 over Finite Field of size 2
        """
        return self.code().base_field()**(self.code().dimension())

    @abstract_method(optional = True)
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.

        This is an abstract method and it should be implemented separately.
        Reimplementing this for each subclass of :class:`Encoder` is not mandatory
        (as a generator matrix only makes sense when the message space is of the `F^k`,
        where `F` is the base field of :meth:`code`.)
        """

class EncodingError(Exception):
    r"""
    Special exception class to indicate an error during encoding or unencoding.
    """
    pass
