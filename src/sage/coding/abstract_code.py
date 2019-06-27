r"""
Generic structures for codes

Class supporting methods available for any type of code (linear, non-linear) and
over any metric (Hamming, rank).

Any class inheriting from AbstractCode can use the encode/decode framework.
"""
#TODO: imports
from sage.modules.module import Module
<<<<<<< HEAD
=======
from sage.categories.modules import Modules
from sage.misc.cachefunc import cached_method
from copy import copy
from .encoder import Encoder
from .decoder import Decoder, DecodingError
from sage.rings.integer import Integer

import inspect
from sage.misc.sageinspect import sage_getargspec

from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF

from sage.structure.parent import Parent
>>>>>>> 4f44acd853... Fixed some dependencies. Category still set up wrong.

#TODO: credits?

def _explain_constructor(cl):
    r"""
    Internal function for use error messages when constructing encoders and decoders.

    EXAMPLES::

<<<<<<< HEAD
        sage: from sage.coding.linear_code import _explain_constructor, LinearCodeSyndromeDecoder
=======
        sage: from sage.coding.linear_code import LinearCodeSyndromeDecoder
        sage: from sage.coding.abstract_code import _explain_constructor
>>>>>>> 4f44acd853... Fixed some dependencies. Category still set up wrong.
        sage: cl = LinearCodeSyndromeDecoder
        sage: _explain_constructor(cl)
        "The constructor requires no arguments.\nIt takes the optional arguments ['maximum_error_weight'].\nSee the documentation of sage.coding.linear_code.LinearCodeSyndromeDecoder for more details."

        sage: from sage.coding.information_set_decoder import LinearCodeInformationSetDecoder
        sage: cl = LinearCodeInformationSetDecoder
        sage: _explain_constructor(cl)
        "The constructor requires the arguments ['number_errors'].\nIt takes the optional arguments ['algorithm'].\nIt accepts unspecified arguments as well.\nSee the documentation of sage.coding.information_set_decoder.LinearCodeInformationSetDecoder for more details."
    """
    if inspect.isclass(cl):
        argspec = sage_getargspec(cl.__init__)
        skip = 2 # skip the self and code arguments
    else:
        # Not a class, assume it's a factory function posing as a class
        argspec = sage_getargspec(cl)
        skip = 1 # skip code argument
    if argspec.defaults:
        args = argspec.args[skip:-len(argspec.defaults)]
        kwargs = argspec.args[-len(argspec.defaults):]
        opts = "It takes the optional arguments {}.".format(kwargs)
    else:
        args = argspec.args[skip:]
        opts = "It takes no optional arguments."
    if args:
        reqs = "The constructor requires the arguments {}.".format(args)
    else:
        reqs = "The constructor requires no arguments."
    if argspec.varargs or argspec.keywords:
        var = "It accepts unspecified arguments as well.\n"
    else:
        var = ""
    return("{}\n{}\n{}See the documentation of {}.{} for more details."\
            .format(reqs, opts, var, cl.__module__, cl.__name__))
<<<<<<< HEAD


class AbstractCode(Module):

    def __init__(self, base_ring, length, default_encoder_name,
                 default_decoder_name, metric='Hamming'):
=======

class AbstractCode(Module):

    def __init__(self, base_ring, length, metric='Hamming'):
>>>>>>> 4f44acd853... Fixed some dependencies. Category still set up wrong.
        """
        """
        _registered_encoders = {}
        _registered_decoders = {}

<<<<<<< HEAD
        self._base_ring = base_ring
        self._length = length
        self._metric = metric

=======
        if not isinstance(length, (int, Integer)):
            raise ValueError("length must be a Python int or a Sage Integer")
        if length <= 0:
            raise ValueError("length must be a non-zero positive integer")
        if not base_ring.is_ring():
            raise ValueError("'base_ring' must be a ring (and {} is not one)".format(base_ring))

        self._length = length
        self._metric = metric

        #self._default_decoder_name = default_decoder_name
        #self._default_encoder_name = default_encoder_name

        cat = Modules(base_ring).FiniteDimensional()
        Parent.__init__(self, base=base_ring, facade=True, category=cat)

>>>>>>> 4f44acd853... Fixed some dependencies. Category still set up wrong.
    def __getstate__(self):
        """
        Used for pickling codes.

        TESTS::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: '_registered_encoders' in C.__getstate__()
            True
        """
        d = super(AbstractCode, self).__getstate__()
        d['_registered_encoders'] = self._registered_encoders
        d['_registered_decoders'] = self._registered_decoders
        return d

    def _repr_(self):
        r"""
        Return an error message requiring to override ``_repr_`` in ``self``.

        As one has to implement specific representation methods (`_repr_` and `_latex_`)
        when writing a new code class which inherits from :class:`AbstractLinearCode`,
        the generic call to `_repr_` has to fail.

        EXAMPLES:

        This was taken from :trac:`20899` (and thus ensures this method fixes what was
        described in this ticket).

        We create a new code class, its dedicated encoder
        and set appropriate parameters::

            sage: from sage.coding.linear_code import AbstractLinearCode
            sage: from sage.coding.encoder import Encoder
            sage: class MyCode(AbstractLinearCode):
            ....:    _registered_encoders = {}
            ....:    _registered_decoders = {}
            ....:    def __init__(self):
            ....:        super(MyCode, self).__init__(GF(5), 10, "Monkey", "Syndrome")
            ....:        self._dimension = 2

            sage: class MonkeyEncoder(Encoder):
            ....:    def __init__(self, C):
            ....:        super(MonkeyEncoder, self).__init__(C)
            ....:    @cached_method
            ....:    def generator_matrix(self):
            ....:        return matrix(GF(5), 2, 10, [ [1]*5 + [0]*5, [0]*5 + [1]*5 ])
            sage: MyCode._registered_encoders["Monkey"] = MonkeyEncoder
            sage: MyCode._registered_decoders["Syndrome"] = codes.decoders.LinearCodeSyndromeDecoder

        We check we get a sensible error message while asking for a string
        representation of an instance of our new class:

            sage: C = MyCode()
            sage: C #random
            Traceback (most recent call last):
            ...
            RuntimeError: Please override _repr_ in the implementation of <class '__main__.MyCode_with_category'>
        """
        raise RuntimeError("Please override _repr_ in the implementation of {}".format(self.parent()))

    def _latex_(self):
        r"""
        Return an error message requiring to override ``_latex_`` in ``self``.

        As one has to implement specific representation methods (`_repr_` and `_latex_`)
        when writing a new code class which inherits from :class:`AbstractLinearCode`,
        the generic call to `_latex_` has to fail.

        EXAMPLES:

        This was taken from :trac:`20899` (and thus ensures this method fixes what was
        described in this ticket).

        We create a new code class, its dedicated encoder
        and set appropriate parameters::

            sage: from sage.coding.linear_code import AbstractLinearCode
            sage: from sage.coding.encoder import Encoder
            sage: class MyCode(AbstractLinearCode):
            ....:    _registered_encoders = {}
            ....:    _registered_decoders = {}
            ....:    def __init__(self):
            ....:        super(MyCode, self).__init__(GF(5), 10, "Monkey", "Syndrome")
            ....:        self._dimension = 2

            sage: class MonkeyEncoder(Encoder):
            ....:    def __init__(self, C):
            ....:        super(MonkeyEncoder, self).__init__(C)
            ....:    @cached_method
            ....:    def generator_matrix(self):
            ....:        return matrix(GF(5), 2, 10, [ [1]*5 + [0]*5, [0]*5 + [1]*5 ])
            sage: MyCode._registered_encoders["Monkey"] = MonkeyEncoder
            sage: MyCode._registered_decoders["Syndrome"] = codes.decoders.LinearCodeSyndromeDecoder

        We check we get a sensible error message while asking for a string
        representation of an instance of our new class:

            sage: C = MyCode()
            sage: latex(C)
            Traceback (most recent call last):
            ...
            RuntimeError: Please override _latex_ in the implementation of <class '__main__.MyCode_with_category'>
        """
        raise RuntimeError("Please override _latex_ in the implementation of {}".format(self.parent()))

    def list(self):
        r"""
        Return a list of all elements of this linear code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: Clist = C.list()
            sage: Clist[5]; Clist[5] in C
            (1, 0, 1, 0, 1, 0, 1)
            True
        """
        return [x for x in self]

<<<<<<< HEAD
    def base_ring(self):
        r"""
        Return the base ring of ``self``.

        EXAMPLES::

            sage: G  = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)
            sage: C.base_field()
            Finite Field of size 2
        """
        return self.base_ring()
=======
    #def base_ring(self):
    #    r"""
    #    Return the base ring of ``self``.
    #    """
    #    return self.base_ring()
>>>>>>> 4f44acd853... Fixed some dependencies. Category still set up wrong.

    def length(self):
        r"""
        Returns the length of this code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.length()
            7
        """
        return self._length

    def metric(self):
        """
        Return the metric of ``self``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
<<<<<<< HEAD
            sage: C.metric
            Hamming
=======
            sage: C.metric()
            'Hamming'
>>>>>>> 4f44acd853... Fixed some dependencies. Category still set up wrong.
        """
        return self._metric

###################### Encoding-Decoding #######################################

    def add_decoder(self, name, decoder):
        r"""
        Adds an decoder to the list of registered decoders of ``self``.

        .. NOTE::

            This method only adds ``decoder`` to ``self``, and not to any member of the class
            of ``self``. To know how to add an :class:`sage.coding.decoder.Decoder`, please refer
            to the documentation of :class:`AbstractLinearCode`.

        INPUT:

        - ``name`` -- the string name for the decoder

        - ``decoder`` -- the class name of the decoder

        EXAMPLES:

        First of all, we create a (very basic) new decoder::

            sage: class MyDecoder(sage.coding.decoder.Decoder):
            ....:   def __init__(self, code):
            ....:       super(MyDecoder, self).__init__(code)
            ....:   def _repr_(self):
            ....:       return "MyDecoder decoder with associated code %s" % self.code()

        We now create a new code::

            sage: C = codes.HammingCode(GF(2), 3)

        We can add our new decoder to the list of available decoders of C::

            sage: C.add_decoder("MyDecoder", MyDecoder)
            sage: sorted(C.decoders_available())
            ['InformationSet', 'MyDecoder', 'NearestNeighbor', 'Syndrome']

        We can verify that any new code will not know MyDecoder::

            sage: C2 = codes.HammingCode(GF(2), 3)
            sage: sorted(C2.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']

        TESTS:

        It is impossible to use a name which is in the dictionary of available decoders::

            sage: C.add_decoder("Syndrome", MyDecoder)
            Traceback (most recent call last):
            ...
            ValueError: There is already a registered decoder with this name
        """
        if self._registered_decoders == self.__class__._registered_decoders:
            self._registered_decoders = copy(self._registered_decoders)
            reg_dec = self._registered_decoders
            if name in reg_dec:
                raise ValueError("There is already a registered decoder with this name")
            reg_dec[name] = decoder
        else:
            if name in self._registered_decoders:
                raise ValueError("There is already a registered decoder with this name")
            reg_dec[name] = decoder

    def add_encoder(self, name, encoder):
        r"""
        Adds an encoder to the list of registered encoders of ``self``.

        .. NOTE::

            This method only adds ``encoder`` to ``self``, and not to any member of the class
            of ``self``. To know how to add an :class:`sage.coding.encoder.Encoder`, please refer
            to the documentation of :class:`AbstractLinearCode`.

        INPUT:

        - ``name`` -- the string name for the encoder

        - ``encoder`` -- the class name of the encoder

        EXAMPLES:

        First of all, we create a (very basic) new encoder::

            sage: class MyEncoder(sage.coding.encoder.Encoder):
            ....:   def __init__(self, code):
            ....:       super(MyEncoder, self).__init__(code)
            ....:   def _repr_(self):
            ....:       return "MyEncoder encoder with associated code %s" % self.code()

        We now create a new code::

            sage: C = codes.HammingCode(GF(2), 3)

        We can add our new encoder to the list of available encoders of C::

            sage: C.add_encoder("MyEncoder", MyEncoder)
            sage: sorted(C.encoders_available())
            ['MyEncoder', 'Systematic']

        We can verify that any new code will not know MyEncoder::

            sage: C2 = codes.HammingCode(GF(2), 3)
            sage: sorted(C2.encoders_available())
            ['Systematic']

        TESTS:

        It is impossible to use a name which is in the dictionary of available encoders::

            sage: C.add_encoder("Systematic", MyEncoder)
            Traceback (most recent call last):
            ...
            ValueError: There is already a registered encoder with this name
        """
        if self._registered_encoders == self.__class__._registered_encoders:
            self._registered_encoders = copy(self._registered_encoders)
            reg_enc = self._registered_encoders
            if name in reg_enc:
                raise ValueError("There is already a registered encoder with this name")
            reg_enc[name] = encoder
        else:
            if name in self._registered_encoders:
                raise ValueError("There is already a registered encoder with this name")
            reg_enc[name] = encoder

    def decode_to_code(self, word, decoder_name=None, *args, **kwargs):
        r"""
        Corrects the errors in ``word`` and returns a codeword.

        INPUT:

        - ``word`` -- a vector of the same length as ``self`` over
          the base field of ``self``

        - ``decoder_name`` -- (default: ``None``) Name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to :meth:`decoder`

        OUTPUT:

        - A vector of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: w_err = word + vector(GF(2), (1, 0, 0, 0, 0, 0, 0))
            sage: C.decode_to_code(w_err)
            (1, 1, 0, 0, 1, 1, 0)

        It is possible to manually choose the decoder amongst the list of the available ones::

            sage: sorted(C.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: C.decode_to_code(w_err, 'NearestNeighbor')
            (1, 1, 0, 0, 1, 1, 0)
        """
        D = self.decoder(decoder_name, *args, **kwargs)
        return D.decode_to_code(word)

    def decode_to_message(self, word, decoder_name=None, *args, **kwargs):
        r"""
        Correct the errors in word and decodes it to the message space.

        INPUT:

        - ``word`` -- a vector of the same length as ``self`` over the
          base field of ``self``

        - ``decoder_name`` -- (default: ``None``) Name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to :meth:`decoder`

        OUTPUT:

        - A vector of the message space of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: C.decode_to_message(word)
            (0, 1, 1, 0)

        It is possible to manually choose the decoder amongst the list of the available ones::

            sage: sorted(C.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: C.decode_to_message(word, 'NearestNeighbor')
            (0, 1, 1, 0)
        """
        return self.unencode(self.decode_to_code(word, decoder_name, *args, **kwargs), **kwargs)

    @cached_method
    def decoder(self, decoder_name=None, *args, **kwargs):
        r"""
        Return a decoder of ``self``.

        INPUT:

        - ``decoder_name`` -- (default: ``None``) name of the decoder which will be
          returned. The default decoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments will be forwarded to the constructor of the decoder
          that will be returned by this method

        OUTPUT:

        - a decoder object

        Besides creating the decoder and returning it, this method also stores
        the decoder in a cache. With this behaviour, each decoder will be created
        at most one time for ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.decoder()
            Syndrome decoder for [7, 4] linear code over GF(2) handling errors of weight up to 1


        If the name of a decoder which is not known by ``self`` is passed,
        an exception will be raised::

            sage: sorted(C.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: C.decoder('Try')
            Traceback (most recent call last):
            ...
            ValueError: There is no Decoder named 'Try'. The known Decoders are: ['InformationSet', 'NearestNeighbor', 'Syndrome']

        Some decoders take extra arguments. If the user forgets to supply these,
        the error message attempts to be helpful::

            sage: C.decoder('InformationSet')
            Traceback (most recent call last):
            ...
            ValueError: Constructing the InformationSet decoder failed, possibly due to missing or incorrect parameters.
            The constructor requires the arguments ['number_errors'].
            It takes the optional arguments ['algorithm'].
            It accepts unspecified arguments as well.
            See the documentation of sage.coding.information_set_decoder.LinearCodeInformationSetDecoder for more details.

        """
        if decoder_name is None:
            decoder_name = self._default_decoder_name
        if decoder_name in self._registered_decoders:
            decClass = self._registered_decoders[decoder_name]
            try:
                return decClass(self, *args, **kwargs)
            except TypeError:
                raise ValueError(
                        "Constructing the {0} decoder failed, possibly due "
                        "to missing or incorrect parameters.\n{1}".format(
                            decoder_name, _explain_constructor(decClass)))
        else:
            raise ValueError(
                    "There is no Decoder named '{0}'. The known Decoders are: "
                    "{1}".format(decoder_name, self.decoders_available()))

    def decoders_available(self, classes=False):
        r"""
        Returns a list of the available decoders' names for ``self``.

        INPUT:

        - ``classes`` -- (default: ``False``) if ``classes`` is set to ``True``,
          return instead a ``dict`` mapping available decoder name to the
          associated decoder class.

        OUTPUT: a list of strings, or a `dict` mapping strings to classes.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.decoders_available()
            ['InformationSet', 'NearestNeighbor', 'Syndrome']

            sage: dictionary = C.decoders_available(True)
            sage: sorted(dictionary.keys())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: dictionary['NearestNeighbor']
            <class 'sage.coding.linear_code.LinearCodeNearestNeighborDecoder'>
        """
        if classes:
            return copy(self._registered_decoders)

        return sorted(self._registered_decoders)

    def encode(self, word, encoder_name=None, *args, **kwargs):
        r"""
        Transforms an element of a message space into a codeword.

        INPUT:

        - ``word`` -- a vector of a message space of the code.

        - ``encoder_name`` -- (default: ``None``) Name of the encoder which will be used
          to encode ``word``. The default encoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        .. NOTE::

            The default encoder always has `F^{k}` as message space, with `k` the dimension
            of ``self`` and `F` the base ring of ``self``.

        One can use the following shortcut to encode a word ::

            C(word)

        OUTPUT:

        - a vector of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector((0, 1, 1, 0))
            sage: C.encode(word)
            (1, 1, 0, 0, 1, 1, 0)
            sage: C(word)
            (1, 1, 0, 0, 1, 1, 0)

        It is possible to manually choose the encoder amongst the list of the available ones::

            sage: sorted(C.encoders_available())
            ['GeneratorMatrix', 'Systematic']
            sage: word = vector((0, 1, 1, 0))
            sage: C.encode(word, 'GeneratorMatrix')
            (1, 1, 0, 0, 1, 1, 0)
        """
        E = self.encoder(encoder_name, *args, **kwargs)
        return E.encode(word)

    @cached_method
    def encoder(self, encoder_name=None, *args, **kwargs):
        r"""
        Returns an encoder of ``self``.

        The returned encoder provided by this method is cached.

        This methods creates a new instance of the encoder subclass designated by ``encoder_name``.
        While it is also possible to do the same by directly calling the subclass' constructor,
        it is strongly advised to use this method to take advantage of the caching mechanism.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          returned. The default encoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to the constructor of the encoder
          this method will return.

        OUTPUT:

        - an Encoder object.

        .. NOTE::

            The default encoder always has `F^{k}` as message space, with `k` the dimension
            of ``self`` and `F` the base ring of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.encoder()
            Generator matrix-based encoder for [7, 4] linear code over GF(2)

        We check that the returned encoder is cached::

            sage: C.encoder.is_in_cache()
            True

        If the name of an encoder which is not known by ``self`` is passed,
        an exception will be raised::

            sage: sorted(C.encoders_available())
            ['GeneratorMatrix', 'Systematic']
            sage: C.encoder('NonExistingEncoder')
            Traceback (most recent call last):
            ...
            ValueError: There is no Encoder named 'NonExistingEncoder'. The known Encoders are: ['GeneratorMatrix', 'Systematic']

        Some encoders take extra arguments. If the user incorrectly supplies
        these, the error message attempts to be helpful::

            sage: C.encoder('Systematic', strange_parameter=True)
            Traceback (most recent call last):
            ...
            ValueError: Constructing the Systematic encoder failed, possibly due to missing or incorrect parameters.
            The constructor requires no arguments.
            It takes the optional arguments ['systematic_positions'].
            See the documentation of sage.coding.linear_code.LinearCodeSystematicEncoder for more details.
        """
        if encoder_name is None:
            encoder_name = self._default_encoder_name
        if encoder_name in self._registered_encoders:
            encClass = self._registered_encoders[encoder_name]
            try:
                return encClass(self, *args, **kwargs)
            except TypeError:
                raise ValueError(
                        "Constructing the {0} encoder failed, possibly due "
                        "to missing or incorrect parameters.\n{1}".format(
                            encoder_name, _explain_constructor(encClass)))
        else:
            raise ValueError(
                    "There is no Encoder named '{0}'. The known Encoders are: "
                    "{1}".format(encoder_name, self.encoders_available()))

    def encoders_available(self, classes=False):
        r"""
        Returns a list of the available encoders' names for ``self``.

        INPUT:

        - ``classes`` -- (default: ``False``) if ``classes`` is set to ``True``,
          return instead a ``dict`` mapping available encoder name to the
          associated encoder class.

        OUTPUT: a list of strings, or a `dict` mapping strings to classes.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.encoders_available()
            ['GeneratorMatrix', 'Systematic']
            sage: dictionary = C.encoders_available(True)
            sage: sorted(dictionary.items())
            [('GeneratorMatrix', <class 'sage.coding.linear_code.LinearCodeGeneratorMatrixEncoder'>),
             ('Systematic', <class 'sage.coding.linear_code.LinearCodeSystematicEncoder'>)]
        """
        if classes:
            return copy(self._registered_encoders)

        return sorted(self._registered_encoders)

    def unencode(self, c, encoder_name=None, nocheck=False, **kwargs):
        r"""
        Returns the message corresponding to ``c``.

        This is the inverse of :meth:`encode`.

        INPUT:

        - ``c`` -- a codeword of ``self``.

        - ``encoder_name`` -- (default: ``None``) name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``nocheck`` -- (default: ``False``) checks if ``c`` is in ``self``. You might set
          this to ``True`` to disable the check for saving computation. Note that if ``c`` is
          not in ``self`` and ``nocheck = True``, then the output of :meth:`unencode` is
          not defined (except that it will be in the message space of ``self``).

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        OUTPUT:

        - an element of the message space of ``encoder_name`` of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: c = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: C.unencode(c)
            (0, 1, 1, 0)
        """
        E = self.encoder(encoder_name, **kwargs)
        return E.unencode(c, nocheck)

    def random_element(self, *args, **kwds):
        """
        Returns a random codeword; passes other positional and keyword
        arguments to ``random_element()`` method of vector space.

        OUTPUT:

        - Random element of the vector space of this code

        EXAMPLES::

            sage: C = codes.HammingCode(GF(4,'a'), 3)
            sage: C.random_element() # random test
            (1, 0, 0, a + 1, 1, a, a, a + 1, a + 1, 1, 1, 0, a + 1, a, 0, a, a, 0, a, a, 1)

        Passes extra positional or keyword arguments through::

            sage: C.random_element(prob=.5, distribution='1/n') # random test
            (1, 0, a, 0, 0, 0, 0, a + 1, 0, 0, 0, 0, 0, 0, 0, 0, a + 1, a + 1, 1, 0, 0)

        TESTS:

        Test that the codeword returned is immutable (see :trac:`16469`)::

            sage: c = C.random_element()
            sage: c.is_immutable()
            True

        Test that codeword returned has the same parent as any non-random codeword
        (see :trac:`19653`)::

            sage: C = codes.random_linear_code(GF(16, 'a'), 10, 4)
            sage: c1 = C.random_element()
            sage: c2 = C[1]
            sage: c1.parent() == c2.parent()
            True
        """
        E = self.encoder()
        M = E.message_space()
        m = M.random_element(*args, **kwds)
        c = E.encode(m)
        c.set_immutable()
        return c
