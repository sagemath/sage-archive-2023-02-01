# -*- coding: utf-8 -*-
r"""
Generic structures for linear codes over the rank metric

TESTS::

This module uses the following experimental feature:
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


Rank Metric
===========
In coding theory, the most common metric is the Hamming metric, where distance
between two codewords is  given by the number of positions in which they differ.
An alternative to this is the rank metric. Take two fields, `F_q` and `F_{q^m}`,
and define a code `C` to be a set of vectors of length `n` with entries from
`F_{q^m}`. Let `c` be a codeword. We can represent it as an `m \times n` matrix
`M` over `F_q`.

A detailed description on the relationship between the two representations can
be found in :meth:`sage.coding.linear_rank_metric.to_matrix_representation`
and :meth:`sage.coding.linear_rank_metric.from_matrix_representation`.

We can define a metric using the rank of the matrix representation of the
codewords. A distance between two codewords `a, b` is the rank of the matrix
representation of `a - b`. A weight of a codeword `c` is the rank of the matrix
representation of `c`.

One of the main uses of rank metric is in the McEliece cryptosystems, where it
significantly decreases the key sizes.


Linear Rank Metric Code and Gabidulin Codes
============================================
The class :class:`sage.coding.linear_rank_metric.LinearRankMetricCode` is a
representative of an unstructured generic example of linear rank metric codes.

Gabidulin codes are the main family of structured linear rank metric codes
studied at the moment. These codes are similar to Reed-Solomon codes.


``AbstractLinearCode``
----------------------

This is a base class designed to contain methods, features and parameters
shared by every linear rank metric code. For instance, generic algorithms for
computing the minimum distance, etc. Many of these algorithms are slow,
e.g. exponential in the code length. It also contains methods for swapping
between vector and matrix representation of elements. However, using the
encoder/decoder framework requires the vector representation of codewords.

``AbstractLinearCode`` is an abstract class for linear rank metric codes,
so any linear code rank metric class should inherit from this class.
Also ``AbstractLinearCode`` should never itself be instantiated.

See :class:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode`
for details and examples.


``LinearRankMetricCode``
--------------

This class is used to represent arbitrary and unstructured linear rank metric
codes. It mostly relies directly on generic methods provided by
``AbstractLinearRankMetricCode``, which means that basic operations on the code
(e.g. computation of the minimum distance) will use slow algorithms.

A ``LinearRankMetricCode`` is instantiated by providing a generator matrix::

    sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
    sage: C = codes.LinearRankMetricCode(GF(4), G)
    sage: C
    [3, 2] linear rank metric code over GF(64)
    sage: C.generator_matrix()
    [1 1 0]
    [0 0 1]
    sage: c = vector(GF(64), (1, 1, 1))
    sage: c in C
    True

Further references
------------------
Read more about rank metric and Gabidulin codes on https://en.wikipedia.org/wiki/Rank_error-correcting_code


AUTHORS:

- Marketa Slukova (2019-07): initial version

TESTS::

    sage: MS = MatrixSpace(GF(2),4,7)
    sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
    sage: C = LinearCode(G)
    sage: C == loads(dumps(C))
    True
"""
from sage.coding.abstract_code import AbstractCode
from sage.rings.integer import Integer
from sage.coding.relative_finite_field_extension import *
from sage.misc.cachefunc import cached_method
from .encoder import Encoder
from .decoder import Decoder
from sage.categories.fields import Fields
from sage.categories.modules import Modules
from sage.structure.parent import Parent
from sage.matrix.constructor import Matrix
from sage.modules.free_module import VectorSpace
from copy import copy
from sage.structure.element import is_Matrix
from sage.rings.integer_ring import ZZ

def to_matrix_representation(base_field, sub_field, v):
    r"""
    Returns a matrix representation of ``v`` over ``sub_field``.

    Let `(b_1, b_2, \ldots, b_m)`, `b_i \in GF(q^m)`, be a base of `GF(q^m)` as
    a vector space over `GF(q)`. Take an element `x \in GF(q^m)`. We can write
    `x` as `x = u_1 b_1 + u_2 b_2 + \ldots u_m b_m`, where `u_i \in GF(q)`. This
    way we can represent an element from `GF(q^m)` as a vector of length `m`
    over `GF(q)`.

    Given a vector ``v`` of length `n` with entries from ``base_field``, we can
    represent each entry as a vector of length `m`, yielding an `m \times n`
    matrix over ``sub_field``.

    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import to_matrix_representation
        sage: x = GF(64).gen()
        sage: a = vector(GF(64), (x + 1, x + 1, 1))
        sage: to_matrix_representation(GF(64), GF(4), a)
        [1 1 1]
        [1 1 0]
        [0 0 0]
    """
    if not v.is_vector():
        raise TypeError("Input must be a vector")
    n = v.length()
    FE = RelativeFiniteFieldExtension(base_field, sub_field)
    m = base_field.degree()//sub_field.degree()
    g = Matrix(sub_field, m, n, lambda i,j: FE.relative_field_representation(v[j])[i])
    return g

def from_matrix_representation(base_field, sub_field, word):
    r"""
    Returns a vector representation of ``word`` over ``base_field``.

    Given an `m \times n` matrix over ``sub_field``, we can represent each of
    its columns as an element of ``base_field``, yielding a vector of length `n`
    over ``base_field``.


    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import from_matrix_representation
        sage: m = Matrix(GF(4), [[1, 1, 1], [1, 1, 0], [0, 0, 0]])
        sage: from_matrix_representation(GF(64), GF(4), m)
        (z6 + 1, z6 + 1, 1)
    """
    if not is_Matrix(word):
        raise TypeError("Input must be a matrix")
    FE = RelativeFiniteFieldExtension(base_field, sub_field)
    v = []
    for i in range(word.ncols()):
        v.append(FE.absolute_field_representation(word.column(i)))
    return vector(v)

def rank_weight(base_field, sub_field, c):
    """
    Returns the rank of ``c`` as a matrix over ``sub_field``.

    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import rank_weight
        sage: x = GF(64).gen()
        sage: a = vector(GF(64), (x + 1, x + 1, 1))
        sage: rank_weight(GF(64), GF(4), a)
        2
    """
    if c.is_vector():
        c = to_matrix_representation(base_field, sub_field, c)
    return c.rank()

def rank_distance(base_field, sub_field, a, b):
    """
    Returns the rank of ``a`` - ``b`` as a matrix over ``sub_field``.

    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import rank_distance
        sage: x = GF(64).gen()
        sage: a = vector(GF(64), (x + 1, x + 1, 1))
        sage: b = vector(GF(64), (1, 0, 0))
        sage: rank_distance(GF(64), GF(4), a, b)
        2
    """
    if a.is_vector():
        a = to_matrix_representation(base_field, sub_field, a)
    if b.is_vector():
        b = to_matrix_representation(base_field, sub_field, b)
    if not (a.nrows() == b.nrows()) and (a.ncols() == b.ncols()):
        raise ValueError("The dimensions of {} and {} have to be identical".format(a, b))
    if not (a.base_ring() == b.base_ring()):
        raise ValueError("The base field of {} and {} has to be the same".format(a, b))
    return (a - b).rank()


class AbstractLinearRankMetricCode(AbstractCode):
    r"""
    Abstract class for linear rank metric codes.

    This class contains methods that can be used on families of linear rank
    metric codes. Every linear rank metric code class should inherit from this
    abstract class.

    This class is intended for codes which are linear over the ``base_field``.

    Codewords of rank metric codes have two representations. They can either be
    written as a vector of length `n` over `GF(q^m)`, or an `m \times n` matrix
    over `GF(q)`. The current implementation of linear rank metric codes
    supports only the vector representation. This means that to use the
    encoder/decoder framework, one has to work with vectors. However, one can
    always get the matrix representation using the
    :meth:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.to_matrix` method. To go
    back to a vector, use the
    :meth:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.from_matrix` method.

    Instructions on how to make a new family of rank metric codes is analogous
    to making a new family of linear codes over the Hamming metric, instructions
    for which are in :class:`sage.coding.linear_code.AbstractLinearCode`.
    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, sub_field, length, dimension, \
            default_encoder_name, default_decoder_name, field_extension=None):
        """
        Initializes mandatory parameters that every linear rank metric code has.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every linear rank metric code.
        The class :class:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode`
        should never be directly instantiated.

        INPUT:

        - ``base_field`` -- the base field of ``self``

        - ``sub_field`` -- the sub field of ``self``

        - ``length`` -- the length of ``self`` (a Python int or a Sage Integer),
          must be > 0 and at most the degree of the field extension

        - ``dimension`` -- the dimension of ``self``

        - ``default_encoder_name`` -- the name of the default encoder of ``self``

        - ``default_decoder_name`` -- the name of the default decoder of ``self``

        - ``field_extension`` -- representation of the elements of the relative
          extension of `base_field` over `sub_field` (default: ``None``)

        EXAMPLES:

        The following example demonstrates how to use subclass
        `AbstractLinearRankMetricCode` for representing a new family of rank
        metric codes. The example family is non-sensical::

            sage: from sage.coding.linear_rank_metric import AbstractLinearRankMetricCode
            sage: class MyRankMetricCode(AbstractLinearRankMetricCode):
            ....:   def __init__(self, base_field, sub_field, length, dimension, generator_matrix):
            ....:       sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.__init__(self, base_field, sub_field, length, dimension, "GeneratorMatrix", "NearestNeighbor")
            ....:       self._generator_matrix = generator_matrix
            ....:   def __iter__(self):
            ....:       from sage.modules.finite_submodule_iter import \
            ....:                                               FiniteFieldsubspace_iterator
            ....:       return FiniteFieldsubspace_iterator(self.generator_matrix(), immutable=True)
            ....:   def generator_matrix(self):
            ....:       return self._generator_matrix
            ....:   def _repr_(self):
            ....:       return "[%d, %d] dummy rank metric code over GF(%s)" % (self.length(), self.dimension(), self.base_field().cardinality())

        We now instantiate a member of our newly made code family::

            sage: generator_matrix = matrix(GF(8), 3, 3,
            ....:                           {(i,i):1 for i in range(3)})
            sage: C = MyRankMetricCode(GF(8), GF(2), 3, 3, generator_matrix)

        We can check its existence and parameters::

            sage: C
            [3, 3] dummy rank metric code over GF(8)

        We can check that it is truly a part of the framework category::

            sage: C.parent()
            <class '__main__.MyRankMetricCode_with_category'>
            sage: C.category()
            Category of finite dimensional vector spaces with basis over Finite Field in z3 of size 2^3

        And any method that works on rank metric linear codes works for our new dummy code::

            sage: C.minimum_distance()
            1
            sage: C.metric()
            'rank'

        TESTS:

        If the dimension is neither a Python int nor a Sage Integer, it will
        raise a exception::

            sage: C = MyRankMetricCode(GF(8), GF(2), 3, 3.0, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: dimension must be a Python int or a Sage Integer

        If ``base_field`` is not a field, an error is raised::

            sage: C = MyRankMetricCode(ZZ, GF(2), 3, 3, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: 'base_field' must be a field (and Integer Ring is not one)

        If ``sub_field`` is not a field, an error is raised::

            sage: C = MyRankMetricCode(GF(8), ZZ, 3, 3, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: 'sub_field' must be a field (and Integer Ring is not one)

        If the name of the default decoder is not known by the class, it will raise
        a exception::

            sage: from sage.coding.linear_rank_metric import AbstractLinearRankMetricCode
            sage: class MyRankMetricCode2(AbstractLinearRankMetricCode):
            ....:   def __init__(self, base_field, sub_field, length, dimension, generator_matrix):
            ....:       sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.__init__(self, base_field, sub_field, length, dimension, "GeneratorMatrix", "Fail")
            ....:       self._generator_matrix = generator_matrix
            ....:   def __iter__(self):
            ....:       from sage.modules.finite_submodule_iter import \
            ....:                                               FiniteFieldsubspace_iterator
            ....:       return FiniteFieldsubspace_iterator(self.generator_matrix(), immutable=True)
            ....:   def generator_matrix(self):
            ....:       return self._generator_matrix
            ....:   def _repr_(self):
            ....:       return "[%d, %d] dummy rank metric code over GF(%s)" % (self.length(), self.dimension(), self.base_field().cardinality())

            sage: C = MyRankMetricCode2(GF(8), GF(2), 3, 3, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: You must set a valid decoder as default decoder for this code, by filling in the dictionary of registered decoders

        If the name of the default encoder is not known by the class, it will raise
        an exception::

            sage: from sage.coding.linear_rank_metric import AbstractLinearRankMetricCode
            sage: class MyRankMetricCode2(AbstractLinearRankMetricCode):
            ....:   def __init__(self, base_field, sub_field, length, dimension, generator_matrix):
            ....:       sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.__init__(self, base_field, sub_field, length, dimension, "Fail", "NearestNeighbor")
            ....:       self._generator_matrix = generator_matrix
            ....:   def __iter__(self):
            ....:       from sage.modules.finite_submodule_iter import \
            ....:                                               FiniteFieldsubspace_iterator
            ....:       return FiniteFieldsubspace_iterator(self.generator_matrix(), immutable=True)
            ....:   def generator_matrix(self):
            ....:       return self._generator_matrix
            ....:   def _repr_(self):
            ....:       return "[%d, %d] dummy rank metric code over GF(%s)" % (self.length(), self.dimension(), self.base_field().cardinality())

            sage: C = MyRankMetricCode2(GF(8), GF(2), 3, 3, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: You must set a valid encoder as default encoder for this code, by filling in the dictionary of registered encoders

        If ``length`` is bigger than the degree of the extension, an error is
        raised::

            sage: C = MyRankMetricCode(GF(64), GF(4), 4, 3, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: 'length' can be at most the degree of the extension, 3

        """
        #TODO: if field_extension is provided, then what? how to check?
        #TODO: check linearity over the big field?

        self._registered_encoders["GeneratorMatrix"] = LinearRankMetricCodeGeneratorMatrixEncoder
        self._registered_decoders["NearestNeighbor"] = LinearRankMetricCodeNearestNeighborDecoder
        self._registered_encoders["Systematic"] = LinearRankMetricCodeSystematicEncoder

        if not isinstance(dimension, (int, Integer)):
            raise ValueError("dimension must be a Python int or a Sage Integer")
        if not base_field.is_field():
            raise ValueError("'base_field' must be a field (and {} is not one)".format(base_field))
        if not sub_field.is_field():
            raise ValueError("'sub_field' must be a field (and {} is not one)".format(sub_field))
        if not field_extension:
            field_extension = RelativeFiniteFieldExtension(base_field, sub_field)
        if not default_encoder_name in self._registered_encoders:
            raise ValueError("You must set a valid encoder as default encoder for this code, by filling in the dictionary of registered encoders")
        if not default_decoder_name in self._registered_decoders:
            raise ValueError("You must set a valid decoder as default decoder for this code, by filling in the dictionary of registered decoders")
        m = base_field.degree() // sub_field.degree()
        if length > m:
            raise ValueError("'length' can be at most the degree of the extension, {}".format(m))


        self._base_field = base_field
        self._sub_field = sub_field
        self._length = length
        self._dimension = dimension
        self._field_extension = field_extension
        self._extension_degree = m

        super(AbstractLinearRankMetricCode, self).__init__(length, "GeneratorMatrix",
        "NearestNeighbor", "rank")
        cat = Modules(base_field).FiniteDimensional().WithBasis().Finite()
        Parent.__init__(self, base=base_field, facade=False, category=cat)

    def base_field(self):
        """
        Returns the base field of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.base_field()
            Finite Field in z6 of size 2^6
        """
        return self._base_field

    def sub_field(self):
        """
        Returns the sub field of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.sub_field()
            Finite Field in z2 of size 2^2
        """
        return self._sub_field

    def dimension(self):
        """
        Returns the dimension of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.dimension()
            2
        """
        return self._dimension

    def cardinality(self):
        r"""
        Return the size of this code.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.cardinality()
            4096
        """
        return self.base_field().order()**self.dimension()

    __len__ = cardinality

    def field_extension(self):
        """
        Returns the field extension of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.field_extension()
            Relative field extension between Finite Field in z6 of size 2^6 and Finite Field in z2 of size 2^2
        """
        return self._field_extension

    def extension_degree(self):
        r"""
        Returns `m`, the degree of the field extension of ``self``.

        Let ``base_field`` be `GF(q^m)` and ``sub_field`` be `GF(q)`. Then this
        function returns `m`.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.extension_degree()
            3
        """

        return self._extension_degree

    def distance(self, left, right):
        """
        Returns the rank of the matrix of ``left`` - ``right``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: x = GF(64).gen()
            sage: a = vector(GF(64), (x + 1, x + 1, 1))
            sage: b = vector(GF(64), (1, 0, 0))
            sage: C.distance(a, b)
            2
        """
        return rank_distance(self._base_field, self._sub_field, left, right)

    def minimum_distance(self):
        r"""
        Returns the minimum distance of ``self``.

        This algorithm simply iterates over all the elements of the code and
        returns the minimum weight.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.minimum_distance()
            1
        """
        d = float('inf')
        for c in self:
            if c == self.zero():
                continue
            d = min(self.weight(c), d)
        return d

    def weight(self, word):
        """
        Returns the weight of the code word, i.e. its rank.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: x = GF(64).gen()
            sage: a = vector(GF(64), (x + 1, x + 1, 1))
            sage: C.weight(a)
            2
        """
        return rank_weight(self._base_field, self._sub_field, word)

    def to_matrix(self, word):
        """
        Returns the matrix representation of a word.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: x = GF(64).gen()
            sage: a = vector(GF(64), (x + 1, x + 1, 1))
            sage: C.to_matrix(a)
            [1 1 1]
            [1 1 0]
            [0 0 0]
        """
        return to_matrix_representation(self._base_field, self._sub_field, word)

    def from_matrix(self, word):
        """
        Returns the vector representation of a word.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: x = GF(64).gen()
            sage: m = Matrix(GF(4), [[1, 1, 1], [1, 1, 0], [0, 0, 0]])
            sage: C.from_matrix(m)
            (z6 + 1, z6 + 1, 1)
        """
        return from_matrix_representation(self._base_field, self._sub_field, word)

    def ambient_space(self):
        r"""
        Returns the ambient vector space of `self`.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.ambient_space()
            Vector space of dimension 3 over Finite Field in z6 of size 2^6
        """
        return VectorSpace(self.base_field(),self.length())

    @cached_method
    def zero(self):
        r"""
        Returns the zero vector of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.zero()
            (0, 0, 0)
        """
        return self.ambient_space().zero()

    def generator_matrix(self, encoder_name=None, **kwargs):
        r"""
        Returns a generator matrix of ``self``.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          used to compute the generator matrix. The default encoder of ``self``
          will be used if default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.generator_matrix()
            [1 1 0]
            [0 0 1]
        """
        E = self.encoder(encoder_name, **kwargs)
        return E.generator_matrix()

    @cached_method
    def parity_check_matrix(self):
        r"""
        Returns the parity check matrix of ``self``.

        The parity check matrix of a linear rank metric code `C` corresponds to
        the generator matrix of the dual code of `C`.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.parity_check_matrix()
            [1 1 0]
        """
        G = self.generator_matrix()
        H = G.right_kernel()
        M = H.basis_matrix()
        M.set_immutable()
        return M

    def syndrome(self, r):
        r"""
        Returns the syndrome of ``r``.

        The syndrome of ``r`` is the result of `H \times r` where `H` is
        the parity check matrix of ``self``. If ``r`` belongs to ``self``,
        its syndrome equals to the zero vector.

        INPUT:

        - ``r`` -- a vector of the same length as ``self``

        OUTPUT:

        - a column vector

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: x = GF(64).gen()
            sage: a = vector(GF(64), (x + 1, x + 1, 1))
            sage: C.syndrome(a)
            (0)
        """
        return self.parity_check_matrix()*r

    def __contains__(self, v):
        r"""
        Returns True if `v` can be coerced into `self`. Otherwise, returns False.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: x = GF(64).gen()
            sage: a = vector(GF(64), (x + 1, x + 1, 1))
            sage: a in C
            True
        """
        if not v in self.ambient_space() or len(v) != self.length():
            return False
        return self.syndrome(v) == 0

    @cached_method
    def information_set(self):
        """
        Return an information set of the code.

        Return value of this method is cached.

        A set of column positions of a generator matrix of a code
        is called an information set if the corresponding columns
        form a square matrix of full rank.

        OUTPUT:

        - Information set of a systematic generator matrix of the code.

        See :meth:`sage.coding.linear_code.AbstractLinearCode.information_set`
        for more information.
        """
        return self.encoder("Systematic").systematic_positions()

    def is_information_set(self, positions):
        """
        Return whether the given positions form an information set.

        INPUT:

        - A list of positions, i.e. integers in the range 0 to `n-1` where `n`
          is the length of `self`.

        OUTPUT:

        - A boolean indicating whether the positions form an information set.

        See :meth:`sage.coding.linear_code.AbstractLinearCode.is_information_set`
        for more information.
        """
        try:
            self.encoder("Systematic", systematic_positions=tuple(positions))
            return True
        except ValueError:
            return False


class LinearRankMetricCode(AbstractLinearRankMetricCode):
    r"""
    Linear rank metric codes over a finite field, represented using a
    generator matrix.

    This class should be used for arbitrary and unstructured linear rank metric
    codes. This means that basic operations on the code, such as the computation
    of the minimum distance, will use generic, slow algorithms.

    If you are looking for constructing a code from a more specific family, see
    if the family has been implemented by investigating ``codes.<tab>``. These
    more specific classes use properties particular to that family to allow
    faster algorithms, and could also have family-specific methods.

    EXAMPLES::

        sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
        sage: C = codes.LinearRankMetricCode(GF(4), G)
        sage: C
        [3, 2] linear rank metric code over GF(64)
        sage: C.base_field()
        Finite Field in z6 of size 2^6
        sage: C.sub_field()
        Finite Field in z2 of size 2^2
        sage: C.length()
        3
        sage: C.dimension()
        2
        sage: C[2]
        (z6, z6, 0)
    """

    def __init__(self, sub_field, generator, field_extension=None):
        r"""
        See the docstring for :meth:`LinearRankMetricCode`.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G) # indirect doctest
            sage: C
            [3, 2] linear rank metric code over GF(64)

        INPUT:

        - ``sub_field`` -- the sub field of ``self``

        - ``generator`` -- a generator matrix over the ``base_field`` with
          dimension `k \times n`, where `k` is the dimension of the code and
          `n` its length.

        - ``field_extension`` -- representation of the elements of the relative
          extension of `base_field` over `sub_field` (default: ``None``)
        """
        base_field = generator.base_ring()
        if not base_field.is_field():
            raise ValueError("'generator' must be defined on a field (not a ring)")

        try:
            basis = None
            if hasattr(generator,"nrows"): # generator matrix case
                if generator.rank() < generator.nrows():
                    basis = generator.row_space().basis()
            else:
                basis = generator.basis() # vector space etc. case
            if not basis is None:
                from sage.matrix.constructor import matrix
                generator = matrix(base_field, basis)
                if generator.nrows() == 0:
                    raise ValueError("this linear code contains no non-zero vector")
        except AttributeError:
            # Assume input is an AbstractLinearRankMetricCode, extract its generator matrix
            generator = generator.generator_matrix()

        self._generator_matrix = generator
        self._dimension = generator.rank()
        self._length = generator.ncols()
        super(LinearRankMetricCode, self).__init__(base_field, sub_field, self._length, self._dimension, "GeneratorMatrix", "NearestNeighbor", field_extension)

    def __iter__(self):
        """
        Return an iterator over the elements of this linear rank metric code.

        EXAMPLES::

            sage: G = Matrix(GF(8), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(2), G)
            sage: [list(c) for c in C if C.weight(c) < 2][:5]
            [[0, 0, 0], [1, 1, 0], [z3, z3, 0], [z3 + 1, z3 + 1, 0], [z3^2, z3^2, 0]]
        """
        from sage.modules.finite_submodule_iter import \
                                                FiniteFieldsubspace_iterator
        return FiniteFieldsubspace_iterator(self.generator_matrix(), immutable=True)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C
            [3, 2] linear rank metric code over GF(64)
        """
        R = self.base_field()
        if R in Fields():
            return "[%s, %s] linear rank metric code over GF(%s)"%(self.length(), self.dimension(), R.cardinality())
        else:
            return "[%s, %s] linear rank metric code over %s"%(self.length(), self.dimension(), R)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: latex(C)
            [3, 2]\textnormal{ Linear rank metric code over }\Bold{F}_{2^{6}}
        """
        return "[%s, %s]\\textnormal{ Linear rank metric code over }%s"\
                % (self.length(), self.dimension(), self.base_field()._latex_())

    def __hash__(self):
        r"""
        Returns the hash value of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: hash(C) #random
            -2463654037231890597
        """
        Str = str(self)
        G = self.generator_matrix()
        return hash((Str, G)) ^ hash(Str) ^ hash(G)

    def generator_matrix(self, encoder_name=None, **kwargs):
        r"""
        Returns a generator matrix of ``self``.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          used to compute the generator matrix. ``self._generator_matrix``
          will be returned if default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: C.generator_matrix()
            [1 1 0]
            [0 0 1]
        """
        if encoder_name is None or encoder_name is 'GeneratorMatrix':
            g = self._generator_matrix
        else:
            g = super(LinearRankMetricCode, self).generator_matrix(encoder_name, **kwargs)
        g.set_immutable()
        return g


####################### encoders ###############################
class LinearRankMetricCodeGeneratorMatrixEncoder(Encoder):
    r"""
    Encoder based on generator_matrix for linear rank metric codes.

    This is the default encoder of a generic linear rank metric code, and should
    never be used for other codes than :class:`LinearRankMetricCode`.
    """

    def __init__(self, code):
        r"""
        INPUT:

        - ``code`` -- The associated :class:`LinearRankMetricCode` of this encoder.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: E = codes.encoders.LinearRankMetricCodeGeneratorMatrixEncoder(C)
            sage: E
            Generator matrix-based encoder for [3, 2] linear rank metric code over GF(64)
        """
        super(LinearRankMetricCodeGeneratorMatrixEncoder, self).__init__(code)

    def __eq__(self, other):
        r"""
        Tests equality between LinearRankMetricCodeGeneratorMatrixEncoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: E1 = codes.encoders.LinearRankMetricCodeGeneratorMatrixEncoder(C)
            sage: E2 = codes.encoders.LinearRankMetricCodeGeneratorMatrixEncoder(C)
            sage: E1 == E2
            True
        """
        return isinstance(other, LinearRankMetricCodeGeneratorMatrixEncoder)\
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: E = codes.encoders.LinearRankMetricCodeGeneratorMatrixEncoder(C)
            sage: E
            Generator matrix-based encoder for [3, 2] linear rank metric code over GF(64)
        """
        return "Generator matrix-based encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: E = codes.encoders.LinearRankMetricCodeGeneratorMatrixEncoder(C)
            sage: latex(E)
            \textnormal{Generator matrix-based encoder for }[3, 2]\textnormal{ Linear rank metric code over }\Bold{F}_{2^{6}}
        """
        return "\\textnormal{Generator matrix-based encoder for }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: E = codes.encoders.LinearRankMetricCodeGeneratorMatrixEncoder(C)
            sage: E.generator_matrix()
            [1 1 0]
            [0 0 1]
        """
        g = self.code().generator_matrix()
        g.set_immutable()
        return g

class LinearRankMetricCodeSystematicEncoder(Encoder):
    r"""
    Encoder based on a generator matrix in systematic form for Linear codes.

    To encode an element of its message space, this encoder first builds a
    generator matrix in systematic form. What is called systematic form here
    is the reduced row echelon form of a matrix, which is not necessarily
    `[I \vert H]`, where `I` is the identity block and `H` the parity block.
    One can refer to :meth:`sage.coding.linear_code.LinearCodeSystematicEncoder.generator_matrix`
    for a concrete example.
    Once such a matrix has been computed, it is used to encode any message
    into a codeword.

    This encoder can also serve as the default encoder of a code defined by a
    parity check matrix: if the :class:`LinearRankMetricCodeSystematicEncoder`
    detects that it is the default encoder, it computes a generator matrix as the
    reduced row echelon form of the right kernel of the parity check matrix.

    For more information on how to use this encoder see
    :meth:`sage.coding.linear_code.LinearCodeSystematicEncoder`

    EXAMPLES:

    The following demonstrates the basic usage of :class:`LinearCodeSystematicEncoder`::

            sage: G = Matrix(GF(256), [[1,1,1,0,0,0,0,0],\
                                     [1,0,0,1,1,0,0,0],\
                                     [0,1,0,1,0,1,0,0],\
                                     [1,1,0,1,0,0,1,1]])
            sage: C = codes.LinearRankMetricCode(GF(2), G)
            sage: E = codes.encoders.LinearRankMetricCodeSystematicEncoder(C)
            sage: E.generator_matrix()
            [1 0 0 0 0 1 1 1]
            [0 1 0 0 1 0 1 1]
            [0 0 1 0 1 1 0 0]
            [0 0 0 1 1 1 1 1]
            sage: E2 = codes.encoders.LinearRankMetricCodeSystematicEncoder(C, systematic_positions=[5,4,3,2])
            sage: E2.generator_matrix()
            [1 0 0 0 0 1 1 1]
            [0 1 0 0 1 0 1 1]
            [1 1 0 1 0 0 1 1]
            [1 1 1 0 0 0 0 0]
    """

    def __init__(self, code, systematic_positions=None):
        r"""

        INPUT:

        - ``code`` -- The associated code of this encoder.

        - ``systematic_positions`` -- (default: ``None``) the positions in codewords that
          should correspond to the message symbols. A list of `k` distinct integers in
          the range 0 to `n-1` where `n` is the length of the code and `k` its
          dimension. The 0th symbol of a message will then be at position
          ``systematic_positions[0]``, the 1st index at position
          ``systematic_positions[1]``, etc. A ``ValueError`` is raised at
          construction time if the supplied indices do not form an information set.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: E = codes.encoders.LinearRankMetricCodeSystematicEncoder(C)
            sage: E
            Systematic encoder for [3, 2] linear rank metric code over GF(64)
        """
        super(LinearRankMetricCodeSystematicEncoder, self).__init__(code)
        self._systematic_positions = tuple(systematic_positions) if systematic_positions else None
        if systematic_positions:
            # Test that systematic_positions consists of integers in the right
            # range. We test that len(systematic_positions) = code.dimension()
            # in self.generator_matrix() to avoid possible infinite recursion.
            if (not all( e in ZZ and e >= 0 and e < code.length() for e in systematic_positions)) \
               or len(systematic_positions) != len(set(systematic_positions)):
                raise ValueError("systematic positions must be a tuple of distinct integers in the range 0 to n-1 where n is the length of the code")
            # Test that the systematic positions are an information set
            self.generator_matrix()


    def __eq__(self, other):
        r"""
        Tests equality between LinearCodeSystematicEncoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: E1 = codes.encoders.LinearRankMetricCodeSystematicEncoder(C)
            sage: E2 = codes.encoders.LinearRankMetricCodeSystematicEncoder(C)
            sage: E1 == E2
            True
        """
        return isinstance(other, LinearRankMetricCodeSystematicEncoder)\
                and self.code() == other.code()\
                and self.systematic_positions() == other.systematic_positions()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: E = codes.encoders.LinearRankMetricCodeSystematicEncoder(C)
            sage: E
            Systematic encoder for [3, 2] linear rank metric code over GF(64)
        """
        return "Systematic encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: E = codes.encoders.LinearRankMetricCodeSystematicEncoder(C)
            sage: latex(E)
            \textnormal{Systematic encoder for }[3, 2]\textnormal{ Linear rank metric code over }\Bold{F}_{2^{6}}
        """
        return "\\textnormal{Systematic encoder for }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix in systematic form of the associated code of ``self``.

        Systematic form here means that a subsets of the columns of the matrix
        forms the identity matrix.

        .. NOTE::

            The matrix returned by this method will not necessarily be `[I \vert H]`, where `I`
            is the identity block and `H` the parity block. If one wants to know which columns
            create the identity block, one can call :meth:`systematic_positions`

        EXAMPLES::

            sage: G = Matrix(GF(256), [[1,1,1,0,0,0,0],\
                                     [1,0,0,1,1,0,0],\
                                     [0,1,0,1,0,1,0],\
                                     [1,1,0,1,0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(2), G)
            sage: E = codes.encoders.LinearRankMetricCodeSystematicEncoder(C)
            sage: E.generator_matrix()
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1]
        """
        C = self.code()
        # This if statement detects if this encoder is itself the default encoder.
        # In this case, attempt building the generator matrix from the parity
        # check matrix
        if hasattr(self, "_use_pc_matrix"):
            if self._use_pc_matrix == 1:
                self._use_pc_matrix = 2
                return C.parity_check_matrix().right_kernel_matrix()
            else:
                raise ValueError("a parity check matrix must be specified if LinearRankMetricCodeSystematicEncoder is the default encoder")
        else:
            self._use_pc_matrix = 1
            M = copy(C.generator_matrix())
        if not self._systematic_positions:
            M.echelonize()
        else:
            k = M.nrows() # it is important that k is *not* computed as C.dimension() to avoid possible cyclic dependency
            if len(self._systematic_positions) != k:
                raise ValueError("systematic_positions must be a tuple of length equal to the dimension of the code")
            # Permute the columns of M and bring to reduced row echelon formb
            perm = self.systematic_permutation()
            M.permute_columns(perm)
            M.echelonize()
            if M[:,:k].is_singular():
                raise ValueError("systematic_positions are not an information set")
            M.permute_columns(perm.inverse())
        M.set_immutable()
        return M

    def systematic_permutation(self):
        r"""
        Returns a permutation which would take the systematic positions into [0,..,k-1]

        EXAMPLES::

            sage: G = Matrix(GF(256), [[1,0,0,0,1,1,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(2), G)
            sage: E = codes.encoders.LinearRankMetricCodeSystematicEncoder(C)
            sage: E.systematic_positions()
            (0, 1, 6)
            sage: E.systematic_permutation()
            [1, 2, 7, 3, 4, 5, 6]
        """
        n = self.code().length()
        systematic_positions = self.systematic_positions()
        k = len(systematic_positions)
        lp = [ None ]*n
        for (i,j) in zip(range(k), systematic_positions):
            lp[i] = j
        j = k
        set_sys_pos = set(systematic_positions)
        for i in range(n):
            if not i in set_sys_pos:
                lp[j] = i
                j += 1
        from sage.combinat.permutation import Permutation
        return Permutation([1 + e for e in lp])

    def systematic_positions(self):
        r"""
        Returns a tuple containing the indices of the columns which form an
        identity matrix when the generator matrix is in systematic form.

        EXAMPLES::

            sage: G = Matrix(GF(256), [[1,1,1,0,0,0,0],\
                                     [1,0,0,1,1,0,0],\
                                     [0,1,0,1,0,1,0],\
                                     [1,1,0,1,0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(2), G)
            sage: E = codes.encoders.LinearRankMetricCodeSystematicEncoder(C)
            sage: E.systematic_positions()
            (0, 1, 2, 3)
        """
        return self._systematic_positions if self._systematic_positions else self.generator_matrix().pivots()



####################### decoders ###############################
class LinearRankMetricCodeNearestNeighborDecoder(Decoder):
    r"""
    Construct a decoder for Linear Rank Metric Codes.

    This decoder will decode to the nearest codeword found.
    """

    def __init__(self, code):
        r"""

        INPUT:

        - ``code`` -- A code associated to this decoder

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D
            Nearest neighbor decoder for [3, 2] linear rank metric code over GF(64)
        """
        super(LinearRankMetricCodeNearestNeighborDecoder, self).__init__(code, code.ambient_space(), \
                code._default_encoder_name)

    def __eq__(self, other):
        r"""
        Tests equality between LinearRankMetricCodeNearestNeighborDecoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: D1 = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D2 = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D1 == D2
            True
        """
        return isinstance(other, LinearRankMetricCodeNearestNeighborDecoder)\
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D
            Nearest neighbor decoder for [3, 2] linear rank metric code over GF(64)
        """
        return "Nearest neighbor decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: latex(D)
            \textnormal{Nearest neighbor decoder for }[3, 2]\textnormal{ Linear rank metric code over }\Bold{F}_{2^{6}}
        """
        return "\\textnormal{Nearest neighbor decoder for }%s" % self.code()._latex_()

    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``word`` and returns a codeword.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self``'s message space
        """
        c_min = self.code().zero()
        h_min = self.code().weight(r)
        for c in self.code():
            if self.code().weight(c-r) < h_min:
                h_min = self.code().weight(c-r)
                c_min = c
        c_min.set_immutable()
        return c_min

    def decoding_radius(self):
        r"""
        Return maximal number of errors ``self`` can decode.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(GF(4), G)
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D.decoding_radius()
            0
        """
        return (self.code().minimum_distance()-1) // 2
