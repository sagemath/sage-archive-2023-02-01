r"""
Decoders for AG codes

This module implements decoders for evaluation and differential AG codes.

The implemented algorithm for unique decoding of AG codes, named K, is from
[LBO2014]_ and [Lee2016]_.

EXAMPLES::

    sage: F.<a> = GF(9)
    sage: A2.<x,y> = AffineSpace(F, 2)
    sage: C = Curve(y^3 + y - x^4)
    sage: Q, = C.places_at_infinity()
    sage: O = C(0,0).place()
    sage: pls = C.places()
    sage: pls.remove(Q)
    sage: pls.remove(O)
    sage: G = -O + 18*Q
    sage: code = codes.EvaluationAGCode(pls, G)
    sage: code
    [26, 15] evaluation AG code over GF(9)
    sage: decoder = code.decoder('K')
    sage: tau = decoder.decoding_radius()
    sage: tau
    4

The ``decoder`` is now ready for correcting vectors received from a noisy
channel::

    sage: channel = channels.StaticErrorRateChannel(code.ambient_space(), tau)
    sage: message_space = decoder.message_space()
    sage: message = message_space.random_element()
    sage: encoder = decoder.connected_encoder()
    sage: sent_codeword = encoder.encode(message)
    sage: received_vector = channel(sent_codeword)
    sage: (received_vector - sent_codeword).hamming_weight()
    4
    sage: decoder.decode_to_code(received_vector) == sent_codeword
    True
    sage: decoder.decode_to_message(received_vector) == message
    True

AUTHORS:

- Kwankyu Lee (2019-03): initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Kwankyu Lee <kwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport cython

from sage.rings.all import PolynomialRing
from sage.rings.infinity import infinity
from sage.rings.function_field.all import FunctionField

from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix

from .encoder import Encoder
from .decoder import Decoder, DecodingError


class EvaluationAGCodeEncoder(Encoder):
    """
    Encoder of an evaluation AG code

    INPUT:

    - ``code`` -- an evaluation AG code

    - ``decoder`` -- a decoder of the code

    EXAMPLES::

        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: F = C.function_field()
        sage: pls = F.places()
        sage: p = C([0,0])
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: code = codes.EvaluationAGCode(D, G)
        sage: dec = code.decoder('K', Q)
        sage: enc = dec.connected_encoder()
        sage: enc
        Encoder for [8, 5] evaluation AG code over GF(4)

    Constructing a decoder for an AG code takes much time. So we save the
    decoder (and the connected encoder) by::

        sage: save(dec, SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')  # not tested

    and load it for later examples.
    """
    def __init__(self, code, decoder=None):
        """
        Initialize.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: enc = dec.connected_encoder()
            sage: TestSuite(enc).run(skip='_test_pickling')
        """
        super().__init__(code)

        if decoder is None:
            decoder = code.decoder('K')

        self._encode = decoder._encode
        self._unencode = decoder._decode

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: enc = dec.connected_encoder()
            sage: {enc: 1}
            {Encoder for [8, 5] evaluation AG code over GF(4): 1}
        """
        return hash((self.code(), self._encode))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: code = codes.EvaluationAGCode(D, G)
            sage: dec1 = code.decoder('K', Q)
            sage: enc1 = dec1.connected_encoder()
            sage: dec2 = code.decoder('K', Q)
            sage: enc2 = dec2.connected_encoder()
            sage: enc1 == enc2
            True
        """
        if self is other:
            return True

        if not isinstance(other, EvaluationAGCodeEncoder):
            return False

        return self.code() == other.code() and self._encode == other._encode

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: enc = dec.connected_encoder()
            sage: enc
            Encoder for [8, 5] evaluation AG code over GF(4)
        """
        return "Encoder for {}" .format(self.code())

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: enc = dec.connected_encoder()
            sage: latex(enc)
            \text{Encoder for }[8, 5]\text{ evaluation AG code over }\Bold{F}_{2^{2}}
        """
        return r"\text{{Encoder for }}{}".format(self.code()._latex_())

    def encode(self, message):
        """
        Return the codeword encoded from the message.

        INPUT:

        - ``message`` -- a vector in the message space

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: enc = dec.connected_encoder()
            sage: msg = enc.message_space().random_element()
            sage: codeword = enc.encode(msg)
            sage: enc.unencode(codeword) == msg
            True
        """
        return self._encode(message)

    def unencode_nocheck(self, codeword):
        """
        Return the message unencoded from ``codeword``.

        INPUT:

        - ``codeword`` -- a vector in the code

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: enc = dec.connected_encoder()
            sage: msg = enc.message_space().random_element()
            sage: codeword = enc.encode(msg)
            sage: enc.unencode(codeword) in enc.message_space()  # indirect doctest
            True
        """
        return self._unencode(codeword)


class DifferentialAGCodeEncoder(Encoder):
    """
    Encoder of a differential AG code

    INPUT:

    - ``code`` -- a differential AG code

    - ``decoder`` -- a decoder of the code

    EXAMPLES::

        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: F = C.function_field()
        sage: pls = F.places()
        sage: p = C([0,0])
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: code = codes.DifferentialAGCode(D, G)
        sage: dec = code.decoder('K', Q)
        sage: enc = dec.connected_encoder(); enc
        Encoder for [8, 3] differential AG code over GF(4)

    Constructing a decoder for an AG code takes much time. So we save the
    decoder (and the connected encoder) by::

        sage: save(dec, SAGE_SRC + '/sage/coding/tests/differential_ag_code')  # not tested

    and load it for later examples.
    """
    def __init__(self, code, decoder=None):
        """
        Initialize.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: enc = dec.connected_encoder()
            sage: TestSuite(enc).run(skip='_test_pickling')
        """
        super().__init__(code)

        if decoder is None:
            decoder = code.decoder('K')

        self._encode = decoder._encode
        self._unencode = decoder._decode

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: enc = dec.connected_encoder()
            sage: {enc: 1}
            {Encoder for [8, 3] differential AG code over GF(4): 1}
        """
        return hash((self.code(), self._encode))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: code = codes.DifferentialAGCode(D, G)
            sage: dec1 = code.decoder('K', Q)
            sage: enc1 = dec1.connected_encoder()
            sage: dec2 = code.decoder('K', Q)
            sage: enc2 = dec2.connected_encoder()
            sage: enc1 == enc2
            True
        """
        if self is other:
            return True

        if not isinstance(other, DifferentialAGCodeEncoder):
            return False

        return self.code() == other.code() and self._encode == other._encode

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: enc = dec.connected_encoder()
            sage: enc
            Encoder for [8, 3] differential AG code over GF(4)
        """
        return "Encoder for {}".format(self.code())

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: enc = dec.connected_encoder()
            sage: latex(enc)
            \text{Encoder for }[8, 3]\text{ differential AG code over }\Bold{F}_{2^{2}}
        """
        return r"\text{{Encoder for }}{}".format(self.code()._latex_())

    def encode(self, message):
        """
        Return the codeword encoded from the message.

        INPUT:

        - ``message`` -- a vector in the message space

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: enc = dec.connected_encoder()
            sage: msg = enc.message_space().random_element()
            sage: codeword = enc.encode(msg)
            sage: enc.unencode(codeword) == msg
            True
        """
        return self._encode(message)

    def unencode_nocheck(self, codeword):
        """
        Return the message unencoded from ``codeword``.

        INPUT:

        - ``codeword`` -- a vector in the code

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: enc = dec.connected_encoder()
            sage: msg = enc.message_space().random_element()
            sage: codeword = enc.encode(msg)
            sage: enc.unencode(codeword) in enc.message_space()  # indirect doctest
            True
        """
        return self._unencode(codeword)


class EvaluationAGCodeUniqueDecoder(Decoder):
    """
    Unique decoder for evaluation AG codes.

    INPUT:

    - ``code`` -- an evaluation AG code

    - ``Q`` -- (optional) a place, not one of the places supporting the code

    - ``basis`` -- (optional) a basis of the space of functions to evaluate

    - ``verbose`` -- if ``True``, verbose information is printed

    EXAMPLES::

        sage: k.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(k, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: F = C.function_field()
        sage: pls = F.places()
        sage: p = C(0,0)
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: code = codes.EvaluationAGCode(D, G)
        sage: dec = code.decoder('K', Q)
        sage: enc = dec.connected_encoder()
        sage: chan = channels.StaticErrorRateChannel(code.ambient_space(), 1)
        sage: rv = chan.transmit(code.random_element())
        sage: enc.encode(dec.decode_to_message(rv)) in code
        True

    If ``basis`` is given, that defines the associated evaluation encoding map::

        sage: basis = tuple(G.basis_function_space())
        sage: dec2 = code.decoder('K', Q, basis)
        sage: enc2 = dec2.connected_encoder()
        sage: f = basis[0]
        sage: cw = vector(f.evaluate(p) for p in D)
        sage: enc2.unencode(cw)
        (1, 0, 0, 0, 0)
        sage: enc2.encode(_) == cw
        True
        sage: f = basis[1]
        sage: cw = vector(f.evaluate(p) for p in D)
        sage: enc2.unencode(cw)
        (0, 1, 0, 0, 0)
        sage: enc2.encode(_) == cw
        True

    The default ``basis`` is given by ``code.basis_functions()``.
    """
    _decoder_type = {'always-succeed'}

    def __init__(self, code, Q=None, basis=None, verbose=False):
        """
        Initialize.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: TestSuite(dec).run()
        """
        if not code.dimension() > 0:
            raise ValueError("no decoder for degenerate codes")

        F = code.base_function_field()
        K = F.constant_base_field()

        if Q is None:
            # try to get a rational place not in the support of the AG code
            deg = 1
            for p in F.places(deg):
                if p not in code._pls:
                    Q = p
                    break
            if Q is None:  # if none, then take a nonrational place
                while Q is None:
                    deg += 1
                    Q = F.get_place(deg)
        elif Q in code._pls:
            raise ValueError("Q is one of the places defining the code")

        if verbose:
            print('auxiliary place: {} of degree {}'.format(Q, Q.degree()))

        super().__init__(code, code.ambient_space(), connected_encoder_name='evaluation')

        if Q.degree() > 1:
            circuit = _EvaluationAGCodeDecoder_K_extension(code._pls, code._G, Q,
                                                           verbose=verbose)
        else:
            circuit = _EvaluationAGCodeDecoder_K(code._pls, code._G, Q,
                                                 verbose=verbose)

        if basis is None:
            basis = code._basis_functions

        C = matrix([circuit.decode(vector(K, [f.evaluate(p) for p in code._pls]))
                    for f in basis])

        self._extension = Q.degree() > 1
        self._K = K
        self._basis = tuple(basis)
        self._C = C
        self._Cinv = C.inverse()

        self._circuit = circuit
        self._info = circuit.info
        self._Q = Q

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: {dec: 1}
            {Unique decoder for [8, 5] evaluation AG code over GF(4): 1}
        """
        return hash((self.code(), self._Q))

    def __eq__(self, other):
        """
        Check whether ``other`` equals ``self``.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: c1 = codes.EvaluationAGCode(D, 5*Q)
            sage: c2 = codes.EvaluationAGCode(D, 5*Q)
            sage: c1 == c2
            True
        """
        if self is other:
            return True

        if not isinstance(other, type(self)):
            return False

        return (self.code() == other.code() and self._Q == other._Q
                and self._basis == other._basis)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: dec
            Unique decoder for [8, 5] evaluation AG code over GF(4)
        """
        return "Unique decoder for {}".format(self.code())

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: latex(dec)
            \text{Unique decoder for }[8, 5]\text{ evaluation AG code over }\Bold{F}_{2^{2}}
        """
        return r"\text{{Unique decoder for }}{}".format(self.code()._latex_())

    def _encode(self, message):
        r"""
        Return the codeword encoded from ``message``.

        INPUT:

        - ``message`` -- a vector to be encoded to a codeword

        TESTS:

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: enc = dec.connected_encoder()
            sage: msg = enc.message_space().random_element()
            sage: dec._decode(dec._encode(msg)) == msg
            True
        """
        K = self._K
        C = self._C
        circuit = self._circuit

        if self._extension:
            internal_message = circuit.lift(vector(K, message)) * C
            return circuit.pull_back(circuit.encode(internal_message))
        else:
            return circuit.encode(vector(K, message) * C)

    def _decode(self, vector, **kwargs):
        r"""
        Return the message decoded from ``vector``.

        INPUT:

        - ``vector`` -- a vector to be decoded to a message

        TESTS:

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: code = dec.code()
            sage: cw = code.random_element()
            sage: dec._encode(dec._decode(cw)) == cw
            True
        """
        Cinv = self._Cinv
        circuit = self._circuit

        if self._extension:
            internal_message = circuit.decode(circuit.lift(vector), **kwargs) * Cinv
            return circuit.pull_back(internal_message)
        else:
            return circuit.decode(vector, **kwargs) * Cinv

    def connected_encoder(self, *args, **kwargs):
        r"""
        Returns the connected encoder for this decoder.

        INPUT:

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to the
          constructor of the connected encoder

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: dec.connected_encoder()
            Encoder for [8, 5] evaluation AG code over GF(4)
        """
        return self.code().encoder(self._connected_encoder_name, self, *args, **kwargs)

    def decoding_radius(self):
        r"""
        Return the decoding radius of the decoder.

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: dec.decoding_radius()
            1
        """
        return self._info['decoding_radius']

    def decode_to_message(self, received_vector, **kwargs):
        r"""
        Return the message decoded from ``received_vector``.

        INPUT:

        - ``received_vector`` -- a vector in the ambient space of the code

        - ``verbose`` -- boolean; if ``True``, verbose information on the decoding process
          is printed

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: enc = dec.connected_encoder()
            sage: code = dec.code()
            sage: chan = channels.StaticErrorRateChannel(code.ambient_space(), 1)
            sage: rv = chan.transmit(code.random_element())
            sage: msg = dec.decode_to_message(rv)
            sage: cw = enc.encode(msg)
            sage: (cw - rv).hamming_weight() == 1
            True
        """
        return self._decode(received_vector, **kwargs)

    def decode_to_code(self, received_vector, **kwargs):
        r"""
        Return the codeword decoded from ``received_vector``.

        INPUT:

        - ``received_vector`` -- a vector in the ambient space of the code

        - ``verbose`` -- boolean; if ``True``, verbose information on the decoding process
          is printed

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/evaluation_ag_code')
            sage: code = dec.code()
            sage: chan = channels.StaticErrorRateChannel(code.ambient_space(), 1)
            sage: rv = chan.transmit(code.random_element())
            sage: cw = dec.decode_to_code(rv)
            sage: (cw - rv).hamming_weight() == 1
            True
        """
        return self._encode(self._decode(received_vector, **kwargs))


class DifferentialAGCodeUniqueDecoder(Decoder):
    """
    Unique decoder for a differential AG codes.

    INPUT:

    - ``code`` -- an evaluation AG code

    - ``Q`` -- (optional) a place, not one of the places supporting the code

    - ``basis`` -- (optional) a basis of the space of differentials to take residues

    - ``verbose`` -- if ``True``, verbose information is printed

    EXAMPLES::

        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: F = C.function_field()
        sage: pls = F.places()
        sage: p = C(0,0)
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: code = codes.DifferentialAGCode(D, G)
        sage: dec = code.decoder('K', Q)
        sage: enc = dec.connected_encoder()
        sage: chan = channels.StaticErrorRateChannel(code.ambient_space(), 2)
        sage: rv = chan.transmit(code.random_element())
        sage: enc.encode(dec.decode_to_message(rv)) in code
        True

    If ``basis`` is given, that defines the associated residue encoding map::

        sage: basis = tuple((G - sum(D)).basis_differential_space())
        sage: dec2 = code.decoder('K', Q, basis)
        sage: enc2 = dec2.connected_encoder()
        sage: w = basis[0]
        sage: cw = vector(w.residue(p) for p in D)
        sage: enc2.unencode(cw)
        (1, 0, 0)
        sage: enc2.encode(_) == cw
        True
        sage: w = basis[1]
        sage: cw = vector(w.residue(p) for p in D)
        sage: enc2.unencode(cw)
        (0, 1, 0)
        sage: enc2.encode(_) == cw
        True

    The default ``basis`` is given by ``code.basis_differentials()``.

    Saving the decoder for later examples and tests::

        sage: save(dec, DOT_SAGE + 'temp/dec')

    The saved object will be removed after the final example.
    """
    _decoder_type = {'always-succeed'}

    def __init__(self, code, Q=None, basis=None, verbose=False):
        """
        Initialize.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: TestSuite(dec).run()
        """
        if not code.dimension() > 0:
            raise ValueError("no decoder for degenerate codes")

        F = code.base_function_field()
        K = F.constant_base_field()

        if Q is None:
            # try to get a rational place not in the support of the AG code
            deg = 1
            for p in F.places(deg):
                if p not in code._pls:
                    Q = p
                    break
            if Q is None:  # then take a nonrational place
                while Q is None:
                    deg += 1
                    Q = F.get_place(deg)
        elif Q in code._pls:
            raise ValueError("Q is one of the places defining the code")

        if verbose:
            print('auxiliary place: {} of degree {}'.format(Q, Q.degree()))

        super().__init__(code, code.ambient_space(), connected_encoder_name='residue')

        if Q.degree() > 1:
            circuit = _DifferentialAGCodeDecoder_K_extension(code._pls, code._G, Q,
                                                             verbose=verbose)
        else:
            circuit = _DifferentialAGCodeDecoder_K(code._pls, code._G, Q,
                                                   verbose=verbose)

        if basis is None:
            basis = code._basis_differentials

        C = matrix([circuit.decode(vector(K, [b.residue(p) for p in code._pls]))
                    for b in basis])

        self._extension = Q.degree() > 1
        self._K = K
        self._basis = tuple(basis)
        self._C = C
        self._Cinv = C.inverse()

        self._circuit = circuit
        self._info = circuit.info
        self._Q = Q

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: {dec: 1}
            {Unique decoder for [8, 3] differential AG code over GF(4): 1}
        """
        return hash((self.code(), self._Q))

    def __eq__(self, other):
        """
        Check whether ``other`` equals ``self``.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: c1 = codes.DifferentialAGCode(D, 5*Q)
            sage: c2 = codes.DifferentialAGCode(D, 5*Q)
            sage: c1 == c2
            True
        """
        if self is other:
            return True

        if not isinstance(other, type(self)):
            return False

        return (self.code() == other.code() and self._Q == other._Q
                and self._basis == other._basis)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: dec
            Unique decoder for [8, 3] differential AG code over GF(4)
        """
        return "Unique decoder for {}".format(self.code())

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: latex(dec)
            \text{Unique decoder for }[8, 3]\text{ differential AG code over }\Bold{F}_{2^{2}}
        """
        return r"\text{{Unique decoder for }}{}".format(self.code()._latex_())

    def _encode(self, message):
        r"""
        Return the codeword encoded from ``message``.

        INPUT:

        - ``message`` -- a vector to be encoded to a codeword

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: enc = dec.connected_encoder()
            sage: msg = enc.message_space().random_element()
            sage: dec._decode(dec._encode(msg)) == msg
            True
        """
        K = self._K
        C = self._C
        circuit = self._circuit

        if self._extension:
            internal_message = circuit.lift(vector(K, message)) * C
            return circuit.pull_back(circuit.encode(internal_message))
        else:
            return circuit.encode(vector(K, message) * C)

    def _decode(self, vector, **kwargs):
        r"""
        Return the message decoded from ``vector``.

        INPUT:

        - ``vector`` -- a vector to be decoded to a message

        TESTS::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: code = dec.code()
            sage: cw = code.random_element()
            sage: dec._encode(dec._decode(cw)) == cw
            True
        """
        Cinv = self._Cinv
        circuit = self._circuit

        if self._extension:
            internal_message = circuit.decode(circuit.lift(vector), **kwargs) * Cinv
            return circuit.pull_back(internal_message)
        else:
            return circuit.decode(vector, **kwargs) * Cinv

    def connected_encoder(self, *args, **kwargs):
        r"""
        Return the connected encoder for this decoder.

        INPUT:

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to the
          constructor of the connected encoder

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: dec.connected_encoder()
            Encoder for [8, 3] differential AG code over GF(4)
        """
        return self.code().encoder(self._connected_encoder_name, self, *args, **kwargs)

    def decoding_radius(self):
        r"""
        Return the decoding radius of the decoder.

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: dec.decoding_radius()
            2
        """
        return self._info['decoding_radius']

    def decode_to_message(self, received_vector, **kwargs):
        r"""
        Return the message decoded from ``received_vector``.

        INPUT:

        - ``received_vector`` -- a vector in the ambient space of the code

        - ``verbose`` -- boolean; if ``True``, verbose information on the decoding process
          is printed

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: enc = dec.connected_encoder()
            sage: code = dec.code()
            sage: chan = channels.StaticErrorRateChannel(code.ambient_space(), 2)
            sage: rv = chan.transmit(code.random_element())
            sage: msg = dec.decode_to_message(rv)
            sage: cw = enc.encode(msg)
            sage: (cw - rv).hamming_weight() == 2
            True
        """
        return self._decode(received_vector, **kwargs)

    def decode_to_code(self, received_vector, **kwargs):
        r"""
        Return the codeword decoded from ``received_vector``.

        INPUT:

        - ``received_vector`` -- a vector in the ambient space of the code

        - ``verbose`` -- boolean; if ``True``, verbose information on the decoding process
          is printed

        EXAMPLES::

            sage: dec = load(SAGE_SRC + '/sage/coding/tests/differential_ag_code')
            sage: enc = dec.connected_encoder()
            sage: code = dec.code()
            sage: chan = channels.StaticErrorRateChannel(code.ambient_space(), 2)
            sage: rv = chan.transmit(code.random_element())
            sage: cw = dec.decode_to_code(rv)
            sage: (cw - rv).hamming_weight() == 2
            True
        """
        return self._encode(self._decode(received_vector, **kwargs))


cdef class _Decoder_K(object):
    """
    Common base class for both differential and evaluation AG code decoder K.
    """
    cdef bint is_differential
    cdef int code_length, designed_distance, gamma, s0, tau
    cdef list code_basis, message_index, hvecs, eta_vecs
    cdef list dR, dRbar
    cdef object mul_mat, coeff_mat
    cdef object W, x

    cdef readonly dict info

    def encode(self, message):
        """
        Encode ``message`` to a codeword.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _EvaluationAGCodeDecoder_K
            sage: circuit = _EvaluationAGCodeDecoder_K(D, G, Q)
            sage: rv = vector([0, 0, 0, a, 0, a, a + 1, 0])
            sage: msg = circuit.decode(rv)
            sage: circuit.decode(circuit.encode(msg)) == msg
            True
        """
        code_basis = self.code_basis
        message_index = self.message_index
        return vector(sum([message[i]*code_basis[i] for i in range(len(message_index))]))

    def degree(self, f):
        """
        Return the degree of polynomial ``f``

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _EvaluationAGCodeDecoder_K
            sage: circuit = _EvaluationAGCodeDecoder_K(D, G, Q)
            sage: circuit.degree(0)
            -Infinity
        """
        if f.is_zero():
            return -infinity
        else:
            return f.degree()

    def exponents(self, int s):
        """
        Return the exponents of the monomial with weighted degree s.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _EvaluationAGCodeDecoder_K
            sage: circuit = _EvaluationAGCodeDecoder_K(D, G, Q)
            sage: circuit.exponents(11)
            (8, 1)
        """
        cdef int i, d

        gamma = self.gamma
        dRbar = self.dRbar  # dWbar for differential AG code

        i = s % gamma
        d = dRbar[i]

        return (s - d) // gamma, i

    def substitution(self, vec, w, k, i):
        """
        Substitute z with (z + w*phi_s).

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q = p.place()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _EvaluationAGCodeDecoder_K
            sage: circuit = _EvaluationAGCodeDecoder_K(D, G, Q)
            sage: W.<x> = F[]
            sage: circuit.substitution(vector([0, a*x^2 + a*x, x + 1, 0]), a, 1, 1)
            (0, 0, x + 1, 0)
        """
        gamma = self.gamma
        mul_mat = self.mul_mat
        W = self.W
        x = self.x

        a = vec[gamma:]
        b = vec[:gamma]
        c = (w * x**k) * sum([a[j] * mul_mat[j][i] for j in range(gamma)])
        return vector(W, (c + b).list() + a.list())

    def decode(self, received_vector, bint verbose=False,
               bint detect_decoding_failure=True,
               bint detect_Q_polynomial=True):
        """
        Return the coefficients that maps to the corrected codeword from the received vector.

        INPUT:

        - ``received_vector`` -- a received vector in the ambient space of the
          code

        - ``verbose`` -- boolean; if ``True``, verbose information is printed

        - ``detect_decoding_failure`` -- boolean; if ``True``, early failure
          detection is activated

        - ``detect_Q_polynomial`` -- boolean; if ``True``, a Q-polynomial is
          detected for fast decoding

        If decoding fails for some reason, ``DecodingError`` is raised. The
        message contained in the exception indicates the type of the decoding
        failure.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _EvaluationAGCodeDecoder_K
            sage: circuit = _EvaluationAGCodeDecoder_K(D, G, Q)
            sage: rv = vector(F, [1, a, 1, a + 1, a + 1, a + 1, 1, a + 1])
            sage: circuit.decode(rv)
            (1, 0, a + 1, a + 1, a)
        """
        cdef int s, sk, si, i
        cdef int k, ip, count, delta, wlt
        cdef bint found_Q

        code_length = self.code_length
        designed_distance = self.designed_distance

        gamma = self.gamma
        dR = self.dR
        dRbar = self.dRbar  # dWbar for differential AG code

        hvecs = self.hvecs
        eta_vecs = self.eta_vecs
        mul_mat = self.mul_mat
        coeff_mat = self.coeff_mat

        message_index = self.message_index
        code_basis = self.code_basis

        s0 = self.s0
        tau = self.tau

        W = self.W
        x = self.x

        # auxiliary functions
        degree = self.degree
        exponents = self.exponents
        substitution = self.substitution

        K = W.base_ring()

        if verbose:
            width = 7 * (K.degree() + 2)
            # auxiliary function for verbose printing
            def vprint_g(g, s):
                if verbose > 1:
                    print(g)
                else:
                    print('[', end='')
                    for i in reversed(range(gamma)):
                        t = g[gamma + i]
                        wd = gamma * degree(t) + dR[i] + s
                        s1 = '{} y{}z'.format(0 if t == 0 else t.lt(), i)
                        if t != 0:
                            s2 = '{:<4}'.format('({})'.format(wd))
                        else:
                            s2 = '(-) '
                        print(('{:>' + str(width) + '} ').format(s1 + s2), end='')
                    for i in reversed(range(gamma)):
                        t = g[i]
                        wd = gamma * degree(t) + dRbar[i]
                        s1 = '{} w{}' if self.is_differential else '{} Y{}'
                        s1 = s1.format(0 if t == 0 else t.lt(), i)
                        if t != 0:
                            s2 = '{:<4}'.format('({})'.format(wd))
                        else:
                            s2 = '(-) '
                        print(('{:>' + str(width) + '} ').format(s1 + s2), end='')
                    print(']')

        message = []

        # construct the initial generators of the interpolation module
        hvec = sum([received_vector[i] * hvecs[i] for i in range(code_length)])

        # weighted degree of hvec
        wd_hvec = max([gamma * degree(hvec[i]) + dRbar[i] for i in range(gamma)])

        if wd_hvec <= 0:
            if verbose:
                print("no error")

            for s in message_index:
                sk, si = exponents(s)
                message.append(hvec[si][sk])
        else:
            mat = []
            for i in range(gamma):
                row = vector(eta_vecs[i].list() + [W.zero() for j in range(gamma)])
                mat.append(row)
            for i in range(gamma):
                std = [W.zero() for j in range(gamma)]
                std[i] = W.one()
                row = vector(sum([-hvec[j] * mul_mat[i][j] for j in range(gamma)]).list() + std)
                mat.append(row)

            gbmat = matrix(mat)

            nu = []
            for i in range(gamma):
                nu.append(gbmat[i, i].lc())

            found_Q = False
            s = wd_hvec

            while s >= s0:
                if verbose:
                    print("# s = {}".format(s))
                    print("generators (leading terms):")
                    for i in reversed(range(gamma)):
                        g = gbmat[gamma + i]
                        print("F{} ".format(i), end='')
                        vprint_g(g, s)
                    for i in reversed(range(gamma)):
                        g = gbmat[i]
                        print("G{} ".format(i), end='')
                        vprint_g(g, s)

                sk, si = exponents(s)
                delta = 0
                mu = []
                i_k = []
                i_prime = []
                i_value = []
                i_count = []
                voting_value = []
                voting_count = []
                for i in range(gamma):
                    # detect decoding failure
                    dlt = degree(gbmat[gamma + i, gamma + i])
                    delta += dlt
                    if detect_decoding_failure and delta > tau:
                        # more errors than tau; declare failure
                        if verbose:
                            print("detected decoding failure")
                        raise DecodingError("more errors than decoding radius")

                    # detect Q-polynomial
                    wlt = gamma * dlt + dR[i]
                    if detect_Q_polynomial and wlt + s + tau < designed_distance:
                        found_Q = True
                        posQ = (s, i)
                        break

                    k, ip = exponents(wlt + s)
                    count = degree(gbmat[ip, ip]) - k
                    i_k.append(k)
                    i_prime.append(ip)
                    i_count.append(count)

                if found_Q:
                    break

                if s > 0 or sk < 0:  # not s in message_index
                    for i in range(gamma):
                        k = i_k[i]
                        ip = i_prime[i]

                        if k < 0:
                            value = K.zero()
                        else:
                            value = -gbmat[gamma + i, ip][k]

                        mu.append(1)
                        i_value.append(value)
                    winner = 0
                else:
                    for i in range(gamma):
                        k = i_k[i]
                        ip = i_prime[i]

                        mui = gbmat[gamma + i, gamma + i].lc() * coeff_mat[i, si]
                        value = -gbmat[gamma + i, ip][k] / mui

                        mu.append(mui)
                        i_value.append(value)

                        cbar = max(i_count[i], 0)
                        try:
                            pos = voting_value.index(value)
                            voting_count[pos] += cbar
                        except ValueError:
                            voting_value.append(value)
                            voting_count.append(cbar)

                    # voting
                    c = -1
                    for i in range(len(voting_value)):
                        if c < voting_count[i]:
                            c = voting_count[i]
                            winner = voting_value[i]

                if verbose:
                    print("i_prime:", i_prime)
                    print("i_count:", i_count)
                    print("i_value:", i_value)

                    if s <= 0 and sk >= 0:  # s in message_index
                        print("voting:", list(zip(voting_value, voting_count)))

                for i in range(gamma):
                    row_i = gbmat[gamma + i]
                    row_ip = gbmat[i_prime[i]]
                    if winner != 0:
                        row_i = substitution(row_i, winner, sk, si)
                        row_ip = substitution(row_ip, winner, sk, si)
                    if i_value[i] == winner:
                        nrow_ip = row_ip
                        nrow_i = row_i
                    else:
                        nnu = mu[i] * (winner - i_value[i])
                        if i_count[i] > 0:
                            nrow_ip = row_i
                            nrow_i = x**i_count[i] * row_i - nnu / nu[i_prime[i]] * row_ip
                            nu[i_prime[i]] = nnu
                        else:
                            nrow_ip = row_ip
                            nrow_i = row_i - nnu / nu[i_prime[i]] * x**(-i_count[i]) * row_ip
                    gbmat[i_prime[i]] = nrow_ip
                    gbmat[gamma + i] = nrow_i

                if s <= 0 and sk >= 0:  # s in message_index
                    if verbose:
                        print("message symbol:", winner)
                    message.insert(0, winner)

                s -= 1

            if found_Q:
                s, i = posQ
                if verbose:
                    print("found a Q-polynomial at s = {}, F{}".format(s, i))
                dlt = gamma * degree(gbmat[gamma + i,gamma + i]) + dR[i]

                while s >= s0:
                    if verbose:
                        print("# s = {}".format(s))
                        print("F{} ".format(i), end='')
                        vprint_g(gbmat[gamma + i], s)
                    sk, si = exponents(s)
                    if s <= 0 and sk >= 0:  # s in message_index
                        k, ip = exponents(dlt + s)
                        mui = gbmat[gamma + i,gamma + i].lc() * coeff_mat[i, si]
                        value = -gbmat[gamma + i, ip][k] / mui
                        if not value.is_zero():
                            gbmat[gamma + i] = substitution(gbmat[gamma + i], value, sk, si)
                        if verbose:
                            print("message symbol:", value)
                        message.insert(0, value)
                    s -= 1

                for j in range(gamma):
                    if not gbmat[gamma + i, j].is_zero():
                        if verbose:
                            print("detected decoding failure at division")
                        raise DecodingError("decoding failed")

        return vector(K, message)


@cython.auto_pickle(True)
cdef class _EvaluationAGCodeDecoder_K(_Decoder_K):
    """
    Unique decoding algorithm K for evaluation AG codes.

    INPUT:

    - ``pls`` -- a list of places of a function field

    - ``G`` -- a divisor of the function field

    - ``Q`` -- a rational place not in ``pls``

    - ``verbose`` -- if ``True``, verbose information is printed.

    TESTS::

        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: pls = C.places()
        sage: p = C([0,0])
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: from sage.coding.ag_code_decoders import _EvaluationAGCodeDecoder_K
        sage: circuit = _EvaluationAGCodeDecoder_K(D, G, Q)
        sage: rv = vector([a, 0, 0, a, 1, 1, a + 1, 0])
        sage: cw = circuit.encode(circuit.decode(rv))
        sage: rv - cw
        (a + 1, 0, 0, 0, 0, 0, 0, 0)
    """
    def __init__(self, pls, G, Q, verbose=False):
        """
        Initialize.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _EvaluationAGCodeDecoder_K
            sage: circuit = _EvaluationAGCodeDecoder_K(D, G, Q)
            sage: TestSuite(circuit).run(skip='_test_pickling')
        """
        cdef int i, s, n, r, d

        D = sum(pls)
        F = D.parent().function_field()
        K = F.constant_base_field()
        W = PolynomialRing(K, name='x')  # working polynomial ring
        x = W.gen()

        # length of the code
        code_length = len(pls)

        # compute gamma
        gamma = 1
        while True:
            if Q.divisor(gamma).dimension() > 1:
                break
            gamma += 1

        # compute xR
        for f in Q.divisor(gamma).basis_function_space():
            if f.valuation(Q) == -gamma:
                xR = f
                break

        # apery R
        dR = [0 for i in range(gamma)]
        yR = [None for i in range(gamma)]
        s = 0
        n = 0
        while n < gamma:
            g = 0
            for b in Q.divisor(s).basis_function_space():
                if b.valuation(Q) == -s:
                    g = b
                    break
            r = s % gamma
            if g != 0 and not yR[r]:
                dR[r] = s
                yR[r] = g
                n += 1
            s += 1

        # gaps of L
        gaps = set()
        for d in dR:
            gaps.update([d - gamma*(i+1) for i in range(d // gamma)])
        gaps = list(gaps)

        # genus of L
        genus = len(gaps)

        # apery Rbar
        dRbar = [0 for i in range(gamma)]
        yRbar = [None for i in range(gamma)]
        s = -G.degree()
        n = 0
        while n < gamma:
            B = (Q.divisor(s) + G).basis_function_space()
            g = 0
            for b in B:
                if b.valuation(Q) + G.multiplicity(Q) == -s:
                    g = b
                    break
            r = s % gamma
            if g != 0 and not yRbar[r]:
                dRbar[r] = s
                yRbar[r] = g
                n += 1
            s += 1

        if verbose:
            print("gamma:", gamma)
            print("x = {}".format(xR))
            print("Apery system of R")
            for i in range(gamma):
                print(" {}: {}, y{} = {}".format(i, dR[i], i, yR[i]))
            print("Apery system of Rbar")
            for i in range(gamma):
                print(" {}: {}, Y{} = {}".format(i, dRbar[i], i, yRbar[i]))

        # ev map for the monomial whose weighted degree is s
        evxR = vector(K, [xR.evaluate(p) for p in pls])
        evyRbar = [vector(K, [yRbar[i].evaluate(p) for p in pls]) for i in range(gamma)]

        self.is_differential = False
        self.code_length = code_length
        self.designed_distance = code_length - G.degree()
        self.gamma = gamma
        self.dR = dR
        self.dRbar = dRbar
        self.W = W
        self.x = x

        degree = self.degree
        exponents = self.exponents

        def monomial(s):
            sk, si = exponents(s)
            return xR**sk * yRbar[si]

        def ev_mon(s):
            sk, si = exponents(s)
            return vector([evxR[i]**sk * evyRbar[si][i] for i in range(code_length)])

        def next(s):  # integer next to s in Rbar
            while True:
                s += 1
                if exponents(s)[0] >= 0:
                    return s

        # minimum of nongaps of Rbar
        s0 = next(-G.degree() - 1)

        # basis of the code ev(L(G))
        message_index = []
        code_basis = []
        s = s0
        v = ev_mon(s)
        V = v.parent()
        while s <= 0:
            if not V.are_linearly_dependent(code_basis + [v]):
                message_index.append(s)
                code_basis.append(v)
            s = next(s)
            v = ev_mon(s)

        # compute a basis of J and h-functions via FGLM algorithm
        def get_eta_basis():
            cdef int num, s

            basis = [None for i in range(gamma)]
            s = s0
            mat = matrix(ev_mon(s))
            delta = [exponents(s)]
            num = 0
            while num < gamma:
                s = next(s)
                e = exponents(s)
                if basis[e[1]] is None:
                    v = ev_mon(s)
                    try:
                        sol = mat.solve_left(v)
                        gen = [W.zero() for i in range(gamma)]
                        for i in range(len(delta)):
                            gen[delta[i][1]] += -sol[i] * x**delta[i][0]
                        gen[e[1]] += x**e[0]
                        basis[e[1]] = vector(gen)
                        num += 1
                    except ValueError:
                        mat = matrix(list(mat) + [v])
                        delta.append(e)

            vecs = [None for i in range(code_length)]
            matinv = mat.inverse()
            for i in range(code_length):
                h = [W.zero() for k in range(gamma)]
                for j in range(code_length):
                    h[delta[j][1]] += matinv[i,j]*x**delta[j][0]
                vecs[i] = vector(h)

            return basis, vecs

        eta_vecs, hvecs = get_eta_basis()

        if verbose:
            print("message indices:", message_index)
            print("eta basis:", eta_vecs)
            print("Lagrange polynomials")
            for i in range(code_length):
                print("h{} = {}".format(i, hvecs[i]))

        # Lee-O'Sullivan bound
        def nu(s):
            cdef int m = 0
            for i in range(gamma):
                e = exponents(s + dR[i])
                m += max(0, degree(eta_vecs[e[1]][e[1]]) - e[0])
            return m

        nus = [nu(s) for s in message_index]
        dLO = min(nus)
        tau = (dLO - 1) // 2

        if verbose:
            print("gaps:", gaps)
            print("genus:", genus)
            print("nus:", nus)
            print("dLO:", dLO)
            print("tau:", tau)

        # the vector form corresponding to f in Rbar
        def vec_form(f):
            r = f
            l = [W.zero() for i in range(gamma)]
            while r != 0:
                s = -r.valuation(Q) - G.multiplicity(Q)
                e = exponents(s)
                mon = xR**e[0] * yRbar[e[1]]
                c = (r / mon).evaluate(Q)
                l[e[1]] += c*x**e[0]
                r -= c*mon
            return vector(l)

        # the matrix of the leading coefficient of y_i*ybar_j and the product
        def get_mul_mat():
            cm = matrix.zero(K, gamma, gamma)
            vm = [[None for j in range(gamma)] for i in range(gamma)]
            for i in range(gamma):
                for j in range(gamma):
                    f = yR[i] * yRbar[j]
                    v = vec_form(f)
                    e = exponents(dR[i] + dRbar[j])
                    cm[i,j] = v[e[1]][e[0]]
                    vm[i][j] = v
            return vm, cm

        mul_mat, coeff_mat = get_mul_mat()

        if verbose:
            print("multiplication table")
            for i in range(gamma):
                for j in range(gamma):
                    print("y{} * Y{}:".format(i, j), mul_mat[i][j])
            print("coefficient array")
            print(coeff_mat)

        self.code_basis = code_basis
        self.message_index = message_index
        self.hvecs = hvecs
        self.eta_vecs = eta_vecs
        self.mul_mat = mul_mat
        self.coeff_mat = coeff_mat
        self.s0 = s0
        self.tau = tau

        info = {}
        info['designed_distance'] = dLO
        info['decoding_radius'] = tau

        self.info = info


@cython.auto_pickle(True)
cdef class _DifferentialAGCodeDecoder_K(_Decoder_K):
    """
    Unique decoding algorithm K for differential AG codes.

    INPUT:

    - ``pls`` -- a list of places of a function field

    - ``G`` -- a divisor of the function field

    - ``Q`` -- a rational place not in ``pls``

    - ``verbose`` -- if ``True``, verbose information is printed

    TESTS::

        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: pls = C.places()
        sage: p = C([0,0])
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: from sage.coding.ag_code_decoders import _DifferentialAGCodeDecoder_K
        sage: circuit = _DifferentialAGCodeDecoder_K(D, G, Q)
        sage: rv = vector([1, a, 1, a, 1, a, a, a + 1])
        sage: cw = circuit.encode(circuit.decode(rv))
        sage: rv - cw
        (0, 0, 0, a + 1, 1, 0, 0, 0)
    """
    def __init__(self, pls, G, Q, verbose=False):
        """
        Initialize.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _DifferentialAGCodeDecoder_K
            sage: circuit = _DifferentialAGCodeDecoder_K(D, G, Q)
            sage: TestSuite(circuit).run(skip='_test_pickling')
        """
        D = sum(pls)
        F = D.parent().function_field()
        K = F.constant_base_field()
        W = PolynomialRing(K, name='x')  # working polynomial ring
        x = W.gen()

        # length of the code
        code_length = len(pls)

        # compute gamma
        gamma = 1
        while True:
            if Q.divisor(gamma).dimension() > 1:
                break
            gamma += 1

        # compute xR
        for xR in Q.divisor(gamma).basis_function_space():
            if xR.valuation(Q) == -gamma:
                break

        # apery R
        dR = [0 for i in range(gamma)]
        yR = [None for i in range(gamma)]
        s = 0
        n = 0
        while n < gamma:
            g = 0
            for b in Q.divisor(s).basis_function_space():
                if b.valuation(Q) == -s:
                    g = b
                    break
            r = s % gamma
            if g != 0 and not yR[r]:
                dR[r] = s
                yR[r] = g
                n += 1
            s += 1

        # gaps of L
        gaps = set()
        for d in dR:
            gaps.update([d - gamma*(i + 1) for i in range(d // gamma)])
        gaps = list(gaps)

        # genus of L
        genus = len(gaps)

        # apery Wbar
        dWbar = [0 for i in range(gamma)]
        wWbar = [None for i in range(gamma)]
        s = -code_length + G.degree() - 2 * genus + 2
        n = 0
        while n < gamma:
            B = (-D + G - Q.divisor(s)).basis_differential_space()
            g = 0
            for b in B:
                if b.valuation(Q) == G.multiplicity(Q) - s:
                    g = b
                    break
            r = s % gamma
            if g != 0 and not wWbar[r]:
                dWbar[r] = s
                wWbar[r] = g
                n += 1
            s += 1

        if verbose:
            print("gamma:", gamma)
            print("x = {}".format(xR))
            print("Apery system of R")
            for i in range(gamma):
                print(" {}: {}, y{} = {}".format(i, dR[i], i, yR[i]))
            print("Apery system of Wbar")
            for i in range(gamma):
                print(" {}: {}, w{} = {}".format(i, dWbar[i], i, wWbar[i]))

        # res map for the monomial whose weighted degree is s
        evxR = vector(K, [xR.evaluate(p) for p in pls])
        reswWbar = [vector(K, [wWbar[i].residue(p) for p in pls]) for i in range(gamma)]

        self.is_differential = True
        self.code_length = code_length
        self.designed_distance = G.degree() - 2 * genus + 2
        self.gamma = gamma
        self.dR = dR
        self.dRbar = dWbar
        self.W = W
        self.x = x

        degree = self.degree
        exponents = self.exponents

        def monomial(s):
            sk, si = exponents(s)
            return xR**sk * wWbar[si]

        def res_mon(s):
            sk, si = exponents(s)
            return vector([evxR[i]**sk * reswWbar[si][i] for i in range(code_length)])

        def next(s):  # integer next to s in Wbar
            while True:
                s += 1
                if exponents(s)[0] >= 0:
                    return s

        # minimum of nongaps of Wbar
        s0 = next(-code_length + G.degree() - 2*genus + 1)

        # basis of the code res(Omega(G))
        message_index = []
        code_basis = []
        s = s0
        v = res_mon(s)
        V = v.parent()
        while s <= 0:
            if not V.are_linearly_dependent(code_basis + [v]):
                message_index.append(s)
                code_basis.append(v)
            s = next(s)
            v = res_mon(s)

        # compute a basis of J and h-functions via FGLM algorithm
        def get_eta_basis():
            basis = [None for i in range(gamma)]
            s = s0
            mat = matrix(res_mon(s))
            delta = [exponents(s)]
            num = 0
            while num < gamma:
                s = next(s)
                e = exponents(s)
                if basis[e[1]] is None:
                    v = res_mon(s)
                    try:
                        sol = mat.solve_left(v)
                        gen = [W.zero() for i in range(gamma)]
                        for i in range(len(delta)):
                            gen[delta[i][1]] += -sol[i]*x**delta[i][0]
                        gen[e[1]] += x**e[0]
                        basis[e[1]] = vector(gen)
                        num += 1
                    except ValueError:
                        mat = matrix(list(mat) + [v])
                        delta.append(e)

            vecs = [None for i in range(code_length)]
            matinv = mat.inverse()
            for i in range(code_length):
                h = [W.zero() for k in range(gamma)]
                for j in range(code_length):
                    h[delta[j][1]] += matinv[i,j]*x**delta[j][0]
                vecs[i] = vector(h)

            return basis, vecs

        eta_vecs, hvecs = get_eta_basis()

        if verbose:
            print("message indices:", message_index)
            print("eta basis:", eta_vecs)
            print("Lagrange polynomials")
            for i in range(code_length):
                print("h{} = {}".format(i, hvecs[i]))

        # Lee-O'Sullivan bound
        def nu(s):
            m = 0
            for i in range(gamma):
                e = exponents(s + dR[i])
                m += max(0, degree(eta_vecs[e[1]][e[1]]) - e[0])
            return m

        nus = [nu(s) for s in message_index]
        dLO = min(nus)
        tau = (dLO - 1) // 2

        if verbose:
            print("gaps:", gaps)
            print("genus:", genus)
            print("nus:", nus)
            print("dLO:", dLO)
            print("tau:", tau)

        # the vector form corresponding to f in Wbar
        def vec_form(f):
            r = f
            l = [W.zero() for i in range(gamma)]
            while r != 0:
                s = -r.valuation(Q) + G.valuation(Q)
                e = exponents(s)
                mon = xR**e[0] * wWbar[e[1]]
                c = (r / mon).evaluate(Q)
                l[e[1]] += c * x**e[0]
                r -= c*mon
            return vector(l)

        # the matrix of the leading coefficient of y_i*w_j and the product
        def get_mul_mat():
            cm = matrix.zero(K, gamma, gamma)
            vm = [[None for j in range(gamma)] for i in range(gamma)]
            for i in range(gamma):
                for j in range(gamma):
                    f = yR[i] * wWbar[j]
                    v = vec_form(f)
                    e = exponents(dR[i] + dWbar[j])
                    cm[i,j] = v[e[1]][e[0]]
                    vm[i][j] = v
            return vm, cm

        mul_mat, coeff_mat = get_mul_mat()

        if verbose:
            print("multiplication table")
            for i in range(gamma):
                for j in range(gamma):
                    print("y{} * w{}:".format(i, j), mul_mat[i][j])
            print("coefficient array")
            print(coeff_mat)

        self.code_basis = code_basis
        self.message_index = message_index
        self.hvecs = hvecs
        self.eta_vecs = eta_vecs
        self.mul_mat = mul_mat
        self.coeff_mat = coeff_mat
        self.s0 = s0
        self.tau = tau

        info = {}
        info['designed_distance'] = dLO
        info['decoding_radius'] = tau

        self.info = info


cdef class _Decoder_K_extension(object):
    """
    Unique decoding algorithm K for AG codes via constant field extension.

    INPUT:

    - ``pls`` -- a list of places of a function field

    - ``G`` -- a divisor of the function field

    - ``Q`` -- a nonrational place

    - ``verbose`` -- if ``True``, verbose information is printed

    TESTS::

        sage: A.<x,y> = AffineSpace(GF(4), 2)
        sage: C = Curve(y^2 + y - x^3)
        sage: pls = C.places()
        sage: F = C.function_field()
        sage: G = 1*F.get_place(4)
        sage: code = codes.EvaluationAGCode(pls, G)
        sage: code.decoder('K')
        Unique decoder for [9, 4] evaluation AG code over GF(4)

    ::

        sage: P.<x,y> = ProjectiveSpace(GF(4), 1)
        sage: C = Curve(P)
        sage: pls = C.places()
        sage: len(pls)
        5
        sage: F = C.function_field()
        sage: G = F.get_place(2).divisor()
        sage: code = codes.EvaluationAGCode(pls, G)
        sage: code.decoder('K')
        Unique decoder for [5, 3] evaluation AG code over GF(4)

    """
    cdef object _embedK, _K
    cdef object decoder_ext
    cdef type decoder_cls
    cdef readonly dict info

    def __init__(self, pls, G, Q, verbose=False):
        """
        Initialize.

        TESTS::

            sage: A.<x,y> = AffineSpace(GF(4), 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: F = C.function_field()
            sage: G = 1*F.get_place(4)
            sage: code = codes.EvaluationAGCode(pls, G)
            sage: decoder = code.decoder('K')
            sage: TestSuite(decoder).run(skip='_test_pickling')
        """
        F = G.parent().function_field()
        K = F.constant_base_field()
        F_base = F.base_field()

        K_ext = K.extension(Q.degree())

        if verbose:
            print('extended constant field:', K_ext)

        F_ext_base = FunctionField(K_ext, F_base.variable_name())

        if F.degree() > 1:
            # construct constant field extension F_ext of F
            def_poly = F.polynomial().base_extend(F_ext_base)
            F_ext = F_ext_base.extension(def_poly, names=def_poly.variable_name())
        else: # rational function field
            F_ext = F_ext_base

        O_ext = F_ext.maximal_order()
        Oinf_ext = F_ext.maximal_order_infinite()

        # embedding of F into F_ext
        embedK = K_ext.coerce_map_from(K)
        embedF_base = F_base.hom(F_ext_base.gen(), embedK)

        if F.degree() > 1:
            embedF = F.hom(F_ext.gen(), embedF_base)
        else:
            embedF = embedF_base

        self._embedK = embedK
        self._K = K

        Div_ext = F_ext.divisor_group()

        def conorm_prime_divisor(pl):
            """
            Conorm map for prime divisors.
            """
            if pl.is_infinite_place():
                ideal = Oinf_ext.ideal([embedF(g) for g in pl.prime_ideal().gens()])
            else:
                ideal = O_ext.ideal([embedF(g) for g in pl.prime_ideal().gens()])
            return ideal.divisor()

        def conorm(d):
            """
            Conorm map from F to F_ext.
            """
            c = Div_ext.zero()
            for pl, mul in d.list():
                c += mul * conorm_prime_divisor(pl)
            return c

        def lift_place(pl):
            """
            Get a place of F_ext lying above pl.
            """
            return conorm_prime_divisor(pl).support()[0]

        pls_ext = [lift_place(pl) for pl in pls]
        G_ext = conorm(G)
        Q_ext = lift_place(Q)

        self.decoder_ext = self.decoder_cls(pls_ext, G_ext, Q_ext, verbose=verbose)
        self.info = self.decoder_ext.info

    def lift(self, v):
        """
        Lift a vector over the base field to a vector over the extension field.

        TESTS::

            sage: A.<x,y> = AffineSpace(GF(4), 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: F = C.function_field()
            sage: G = 1*F.get_place(4)
            sage: code = codes.EvaluationAGCode(pls, G)
            sage: decoder = code.decoder('K')
            sage: lift = decoder._circuit.lift
            sage: pull_back = decoder._circuit.pull_back
            sage: v = code.random_element()
            sage: pull_back(lift(v)) == v
            True
        """
        embedK = self._embedK
        return vector(embedK(e) for e in v)

    def pull_back(self, v):
        """
        Pull back a vector over the extension field to a vector over the base field.

        TESTS::

            sage: A.<x,y> = AffineSpace(GF(4), 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: F = C.function_field()
            sage: G = 1*F.get_place(4)
            sage: code = codes.DifferentialAGCode(pls, G)
            sage: decoder = code.decoder('K')
            sage: lift = decoder._circuit.lift
            sage: pull_back = decoder._circuit.pull_back
            sage: v = code.random_element()
            sage: pull_back(lift(v)) == v
            True
        """
        K = self._K
        return vector(K(e) for e in v)

    def encode(self, message, **kwargs):
        """
        Encode ``message`` to a codeword.

        TESTS::

            sage: A.<x,y> = AffineSpace(GF(4), 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: F = C.function_field()
            sage: G = 1*F.get_place(4)
            sage: code = codes.EvaluationAGCode(pls, G)
            sage: decoder = code.decoder('K')
            sage: cw = code.random_element()
            sage: circuit = decoder._circuit
            sage: circuit.encode(circuit.decode(cw)) == cw
            True
        """
        return self.decoder_ext.encode(message, **kwargs)

    def decode(self, received_vector, **kwargs):
        """
        Decode received vector to a message.

        INPUT:

        - ``received_vector`` -- a vector in the ambient space of the code

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(GF(4), 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: F = C.function_field()
            sage: G = 1*F.get_place(4)
            sage: code = codes.DifferentialAGCode(pls, G)
            sage: decoder = code.decoder('K')
            sage: cw = code.random_element()
            sage: circuit = decoder._circuit
            sage: circuit.encode(circuit.decode(cw)) == cw
            True
        """
        return self.decoder_ext.decode(received_vector, **kwargs)


@cython.auto_pickle(True)
cdef class _EvaluationAGCodeDecoder_K_extension(_Decoder_K_extension):
    """
    Unique decoding algorithm K for AG codes via constant field extension.

    INPUT:

    - ``pls`` -- a list of places of a function field

    - ``G`` -- a divisor of the function field

    - ``Q`` -- a nonrational place

    - ``verbose`` -- if ``True``, verbose information is printed

    EXAMPLES::

        sage: A.<x,y> = AffineSpace(GF(4), 2)
        sage: C = Curve(y^2 + y - x^3)
        sage: pls = C.places()
        sage: F = C.function_field()
        sage: G = 1*F.get_place(4)
        sage: code = codes.EvaluationAGCode(pls, G)
        sage: code.decoder('K')
        Unique decoder for [9, 4] evaluation AG code over GF(4)
    """
    def __init__(self, pls, G, Q, verbose=False):
        """
        Initialize.

        TESTS::

            sage: A.<x,y> = AffineSpace(GF(4), 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: F = C.function_field()
            sage: G = 1*F.get_place(3)
            sage: code = codes.EvaluationAGCode(pls, G)
            sage: decoder = code.decoder('K')
            sage: TestSuite(decoder).run(skip='_test_pickling')
        """
        self.decoder_cls = _EvaluationAGCodeDecoder_K

        super().__init__(pls, G, Q, verbose=verbose)


@cython.auto_pickle(True)
cdef class _DifferentialAGCodeDecoder_K_extension(_Decoder_K_extension):
    """
    Unique decoding algorithm K for AG codes via constant field extension.

    INPUT:

    - ``pls`` -- a list of places of a function field

    - ``G`` -- a divisor of the function field

    - ``Q`` -- a nonrational place

    - ``verbose`` -- if ``True``, verbose information is printed

    EXAMPLES::

        sage: A.<x,y> = AffineSpace(GF(4), 2)
        sage: C = Curve(y^2 + y - x^3)
        sage: pls = C.places()
        sage: F = C.function_field()
        sage: G = 1*F.get_place(4)
        sage: code = codes.DifferentialAGCode(pls, G)
        sage: code.decoder('K')
        Unique decoder for [9, 5] differential AG code over GF(4)
    """
    def __init__(self, pls, G, Q, verbose=False):
        """
        Initialize.

        TESTS::

            sage: A.<x,y> = AffineSpace(GF(4), 2)
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: F = C.function_field()
            sage: G = 1*F.get_place(3)
            sage: code = codes.DifferentialAGCode(pls, G)
            sage: decoder = code.decoder('K')
            sage: TestSuite(decoder).run(skip='_test_pickling')
        """
        self.decoder_cls = _DifferentialAGCodeDecoder_K

        super().__init__(pls, G, Q, verbose=verbose)
