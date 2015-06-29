r"""
Decoder

Representation of an error-correction algorithm for a code.

AUTHORS:

- David Joyner (2009-02-01): initial version
- David Lucas (2015-06-29): abstract class version

"""
#*****************************************************************************
#       Copyright (C) 2009 David Joyner <wdjoyner@gmail.com>
#                     2015 David Lucas <david.lucas@inria.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.structure.sage_object import SageObject

class Decoder(SageObject):
    r"""
    Abstract top-class for :class:`Decoder` objects.

    Every decoder class should inherit from this abstract class.

    This class provides:

    - ``code``, the associated code of the decoder

    - ``input_space``, the space of every input for the decoder

    - ``connected_encoder_name``, the name of the encoder used to recover the original message
      from the word to decode

    - ``decoder_type``, a set of keywords which describe the behaviour of the decoder.

    - some methods for decoder objects

    To implement an decoder, you need to:

    - inherit from :class:`Decoder`

    - call :class:`Decoder`'s :meth:`__init__` in the subclass constructor.
      Example: ``super(SubclassName, self).__init__(code, input_space,
      connected_encoder_name, decoder_type)``.
      By doing that, your subclass will have all the parameters described above initialized.
      You need of course to complete the constructor by adding any additional parameter
      needed to describe properly the decoder defined in the subclass.

    Then, you need to override one of decoding methods, either :meth:`decode_to_code` or
    :meth:`decode_to_message`. You also need to override the method :meth:`decoding_radius`.

    Equality methods (``__eq__`` and ``__ne__``) might be useful for decoding in advanced
    codes constructions (like concatenated codes). It is thus strongly encouraged to override
    them.

    As :class:`Decoder` is not designed to be implemented, it does not have any representation
    methods. You should implement ``_repr_`` and ``_latex_`` methods in the sublclass.

    .. NOTE::

        For consistency on decoders, please follow this convention on names for subclasses:
        for a new decoder named ``DecName``, for code family ``CodeFam``, call it
        ``CodeFamDecNameDecoder``.

    """

    def __init__(self, code, input_space, connected_encoder_name,\
        decoder_type):
        r"""
        Initializes mandatory parameters for :class:`Decoder` objects.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every decoder. An abstract
        decoder object should never be created.

        INPUT:

        - ``code`` -- the associated code of ``self``

        - ``input_space`` -- the input space of ``self``

        - ``connected_encoder_name`` -- the associated encoder, which will be
          used by ``self`` to recover elements from the message space

        - ``decoder_type`` -- a set of types of ``self``. Describes the
          behaviour of ``self``.


        EXAMPLES:

        We first create a new :class:`Decoder` subclass::

            sage: class DecoderExample(sage.coding.decoder.Decoder):
            ....:   def __init__(self, code):
            ....:       in_space = code.base_field()
            ....:       connected_enc = "GeneratorMatrix"
            ....:       decoder_type = {"unique", "always-succeed", "is_example"}
            ....:       super(DecoderExample, self).__init__(code, in_space, connected_enc, decoder_type)

        We now create a member of our brand new class::

            sage: G = Matrix(GF(2), [[1, 0, 0, 1], [0, 1, 1, 1]])
            sage: C = LinearCode(G)
            sage: D = DecoderExample(C)

        We can check its parameters::

            sage: D.input_space()
            Finite Field of size 2
            sage: D.decoder_type()
            {'always-succeed', 'is_example', 'unique'}
            sage: D.connected_encoder()
            Generator matrix-based encoder for the Linear code of length 4, dimension 2 over Finite Field of size 2
            sage: D.code()
            Linear code of length 4, dimension 2 over Finite Field of size 2
        """
        self._code = code
        self._input_space = input_space
        self._connected_encoder_name = connected_encoder_name
        self._decoder_type = decoder_type

    def decoder_type(self):
        r"""
        Returns the set of types of ``self``. These types describe the nature of ``self``
        and its decoding algorithm.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1, 0, 0, 1], [0, 1, 1, 1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: D.decoder_type()
            {'always-succeed', 'complete', 'hard-decision', 'unique'}
        """
        return self._decoder_type

    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``r`` and returns a codeword.

        This is a default implementation which assumes that the method
        :meth:`decode_to_message()` has been implemented, else it returns an exception.

        INPUT:

        - ``r`` -- a element of the input space of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: word in C
            True
            sage: w_err = word + vector(GF(2), (1, 0, 0, 0, 0, 0, 0))
            sage: w_err in C
            False
            sage: D = C.decoder()
            sage: D.decode_to_code(w_err)
            (1, 1, 0, 0, 1, 1, 0)
        """
        if hasattr(self, "defaulting_decode_to_message"):
            raise NotImplementedError
        else:
            word = self.decode_to_message(r)
            return self.connected_encoder().encode(word)

    def connected_encoder(self):
        r"""
        Returns the connected encoder of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: D.connected_encoder()
            Generator matrix-based encoder for the Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        return self.code().encoder(name=self._connected_encoder_name)

    def decode_to_message(self, r):
        r"""
        Decodes ``r`` to the message space of ``self.code()``.

        This is a default implementation, which assumes that the method
        :meth:`decode_to_code()` has been implemented, else it returns an exception.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: w_err = word + vector(GF(2), (1, 0, 0, 0, 0, 0, 0))
            sage: D = C.decoder()
            sage: D.decode_to_message(w_err)
            (0, 1, 1, 0)
        """
        self.defaulting_decode_to_message = True
        return self.code().unencode(self.decode_to_code(r))

    def code(self):
        r"""
        Returns the code in which ``self.decode_to_code()`` has its output.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: D.code()
            Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        return self._code

    def message_space(self):
        r"""
        Returns the message space of connected encoder of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: D.message_space()
            Vector space of dimension 4 over Finite Field of size 2
        """
        return self.connected_encoder().message_space()

    def input_space(self):
        r"""
        Returns the input space of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: D.input_space()
            Vector space of dimension 7 over Finite Field of size 2
        """
        if hasattr(self, "_input_space"):
            return self._input_space
        else:
            raise NotImplementedError("Decoder does not have an _input_space parameter")

    @abstract_method(optional = True)
    def decoding_radius(self, **kwargs):
        r"""
        Returns the maximal number of errors that ``self`` is able to correct.

        This is an abstract method and should be implemented in subclasses.
        """
        raise NotImplementedError

class DecodingFailure(Exception):
    r"""
    Special exception class to indicate a failure during decoding.
    """
    pass
"""
Decoding methods for linear error-correcting codes.

Methods implemented:

 * nearest neighbor
 * syndrome

AUTHOR:

- David Joyner (2009-02-01): initial version

.. TODO::

    Add lots more methods!
"""
#*****************************************************************************
#       Copyright (C) 2009 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.decorators import rename_keyword

def syndrome(C, v):
    """
    The vector v represents a received word, so should
    be in the same ambient space V as C. Returns the
    elements in V (including v) which belong to the
    syndrome of v (ie, the coset v+C, sorted by weight).

    EXAMPLES::

        sage: C = codes.HammingCode(2,GF(3)); C
        Linear code of length 4, dimension 2 over Finite Field of size 3
        sage: V = VectorSpace(GF(3), 4)
        sage: v = V([0, 2, 0, 1])
        sage: from sage.coding.decoder import syndrome
        sage: syndrome(C, v)
         [(0, 0, 1, 0), (0, 2, 0, 1), (2, 0, 0, 2), (1, 1, 0, 0), (2, 2, 2, 0), (1, 0, 2, 1), (0, 1, 2, 2), (1, 2, 1, 2), (2, 1, 1, 1)]

    """
    V = C.ambient_space()
    if not isinstance(v, list):
        v = v.list()
    v = V(v)
    coset = [[c + v, (c + v).hamming_weight()] for c in C]
    return [x[0] for x in sorted(coset, key=lambda x: x[1])]

def coset_leader(C, v):
    """
    The vector v represents a received word, so should
    be in the same ambient space V as C. Returns an
    element of the syndrome of v of lowest weight.

    EXAMPLES::

        sage: C = codes.HammingCode(2,GF(3)); C
        Linear code of length 4, dimension 2 over Finite Field of size 3
        sage: V = VectorSpace(GF(3), 4)
        sage: v = V([0, 2, 0, 1])
        sage: from sage.coding.decoder import coset_leader
        sage: coset_leader(C, v)
        ((0, 0, 1, 0), 1)
        sage: coset_leader(C, v)[0]-v in C
        True

    """
    coset = [[c + v, (c + v).hamming_weight()] for c in C]
    wts = [x[1] for x in coset]
    min_wt = min(wts)
    s = C[0]  # initializing
    w = v.hamming_weight()  # initializing
    for x in coset:
        if x[1] == min_wt:
            w = x[1]
            s = x[0]
            break
    return s, w

@rename_keyword(deprecation=6094, method="algorithm")
def decode(C, v, algorithm="syndrome"):
    """
    The vector v represents a received word, so should
    be in the same ambient space V as C. Returns an
    element in C which is closest to v in the Hamming
    metric.

    Methods implemented include "nearest neighbor" (essentially
    a brute force search) and "syndrome".

    EXAMPLES::

        sage: C = codes.HammingCode(2,GF(3))
        sage: V = VectorSpace(GF(3), 4)
        sage: v = V([0, 2, 0, 1])
        sage: v in C
        False
        sage: from sage.coding.decoder import decode
        sage: c = decode(C, v);c
        (0, 2, 2, 1)
        sage: c in C
        True
        sage: c = decode(C, v, algorithm="nearest neighbor");c
        (0, 2, 2, 1)
        sage: C = codes.HammingCode(3,GF(3)); C
        Linear code of length 13, dimension 10 over Finite Field of size 3
        sage: V = VectorSpace(GF(3), 13)
        sage: v = V([2]+[0]*12)
        sage: decode(C, v)  # long time (9s on sage.math, 2011)
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    """
    V = C.ambient_space()
    if not isinstance(v, list):
        v = v.list()
    v = V(v)
    if algorithm == "nearest neighbor":
        diffs = [[c - v, (c - v).hamming_weight()] for c in C]
        diffs.sort(key=lambda x: x[1])
        return diffs[0][0] + v
    if algorithm == "syndrome":
        return -V(syndrome(C, v)[0]) + v
