#TODO: explain what rank metric is, how it works and what it is used for
#TODO: add generator matrix encoder
#TODO: add nearest neighbout decoder
#TODO: add minimum distance method

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

def to_matrix_representation(base_field, sub_field, v):
    """
    """
    if not v.is_vector():
        raise TypeError("Input must be a vector")
    n = v.length()
    FE = RelativeFiniteFieldExtension(base_field, sub_field)
    m = base_field.degree()//sub_field.degree()
    g = Matrix(sub_field, m, n, lambda i,j: FE.relative_field_representation(v[j])[i])
    return g

def from_matrix_representation(base_field, sub_field, m):
    """
    """
    if not m.is_matrix():
        raise TypeError("Input must be a matrix")
    FE = RelativeFiniteFieldExtension(base_field, sub_field)
    v = []
    for i in range(m.ncols()):
        v.append(FE.absolute_field_representation(m.column(i)))
    return vector(v)

def rank_weight(base_field, sub_field, c):
    """
    """
    if c.is_vector():
        c = to_matrix_representation(base_field, sub_field, c)
    return c.rank()

def rank_distance(base_field, sub_field, a, b):
    """
    """
    if a.is_vector():
        a = to_matrix_representation(base_field, sub_field, a)
    if b.is_vector():
        b = to_matrix_representation(base_field, sub_field, b)
    return (a - b).rank()


class AbstractLinearRankMetricCode(AbstractCode):
    """
    Abstract class for linear rank metric codes.

    This class contains methods that can be used on families of linear rank
    metric codes. Every linear rank metric code class should inherit from this
    abstract class.

    This class is intended for codes which are linear over the ``base_field``.

    TODO: example of how to make a class of linear rank metric codes
    TODO: talk about representation
    TODO: add experimental trigger
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

        - ``length`` -- the length of ``self`` (a Python int or a Sage Integer,
          must be > 0)

        - ``dimension`` -- the dimension of ``self``

        - ``field_extension`` -- representation of the elements of the relative
          extension of `base_field` over `sub_field` (default: ``None``)

        - ``default_encoder_name`` -- the name of the default encoder of ``self``

        - ``default_decoder_name`` -- the name of the default decoder of ``self``
        """
        #TODO: check that sub_field is a sub_field of base_field
        #TODO: if field_extension is provided, then what? how to check?
        #TODO: check linearity over the big field
        #TODO: check category framework with Dima

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

        self._base_field = base_field
        self._sub_field = sub_field
        self._length = length
        self._dimension = dimension
        self._field_extension = field_extension

        super(AbstractLinearRankMetricCode, self).__init__(length, "GeneratorMatrix",
        "NearestNeighbor", "rank")
        cat = Modules(base_field).FiniteDimensional().WithBasis().Finite()
        Parent.__init__(self, base=base_field, facade=False, category=cat)

    def base_field(self):
        """
        Returns the base field of ``self``.
        """
        return self._base_field

    def sub_field(self):
        """
        Returns the sub field of ``self``.
        """
        return self._sub_field

    def dimension(self):
        """
        Returns the dimension of ``self``.
        """
        return self._dimension

    def cardinality(self):
        r"""
        Return the size of this code.
        """
        return self.base_field().order()**self.dimension()

    __len__ = cardinality

    def field_extension(self):
        """
        Returns the field extension of ``self``.
        """
        return self._field_extension

    def distance(self, left, right):
        """
        Returns the rank of the matrix of ``left`` - ``right``.
        """
        return rank_distance(self._base_field, self._sub_field, left, right)

    def minimum_distance(self):
        r"""
        Return an error requiring to override ``minimum_distance`` in ``self``.

        There is currently no general algorithm calculating the minimum distance
        of linear rank metric codes. One has to implement the specific method
        when writing a new code class which inherits from
        :class:`AbstractLinearRankMetricCode`.
        The generic call to ``minimum_distance`` has to fail.
        """

        #TODO: check that the formatting self.parent() works
        raise RuntimeError("Please override minimum_distance in the implementation of {}".format(self.parent()))

    def weight(self, word):
        """
        Returns the weight of the code word - its rank.
        """
        return rank_weight(self._base_field, self._sub_field, word)

    def to_matrix(self, word):
        """
        Returns the matrix representation of a word.
        """
        return to_matrix_representation(self._base_field, self._sub_field, word)

    def from_matrix(self, word):
        """
        Returns the vector representation of a word.
        """
        return from_matrix_representation(self._base_field, self._sub_field, word)

    def ambient_space(self):
        r"""
        Returns the ambient vector space of `self`.
        """
        return VectorSpace(self.base_field(),self.length())

    @cached_method
    def zero(self):
        r"""
        Returns the zero vector of ``self``.
        """
        return self.ambient_space().zero()


class LinearRankMetricCode(AbstractLinearRankMetricCode):
    r"""
    Linear rank metric codes over a finite field, represented using a
    generator matrix.

    This class should be used for arbitrary and unstructured linear rank metric
    codes. This means that basic operations on the code, such as the computation
    of the minimum distance, will use generic, slow algorithms.

    If you are looking for constructing a code from a more specific family, see
    if the family has been implemented by investigating `codes.<tab>`. These
    more specific classes use properties particular to that family to allow
    faster algorithms, and could also have family-specific methods.
    """

    def __init__(self, sub_field, generator, field_extension=None):
        r"""
        See the docstring for :meth:`LinearRankMetricCode`.

        INPUT:

        - ``generator`` -- a generator matrix over the ``base_field`` with
          dimension `k \times n`, where `k` is the dimension of the code and
          `n` its length.

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
            # Assume input is an AbstractLinearCode, extract its generator matrix
            generator = generator.generator_matrix()

        self._generator_matrix = generator
        self._dimension = generator.rank()
        self._length = generator.ncols()
        super(LinearRankMetricCode, self).__init__(base_field, sub_field, self._length, self._dimension, "GeneratorMatrix", "NearestNeighbor")

    def __iter__(self):
        """
        Return an iterator over the elements of this linear code.
        """
        from sage.modules.finite_submodule_iter import \
                                                FiniteFieldsubspace_iterator
        return FiniteFieldsubspace_iterator(self.generator_matrix(), immutable=True)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.
        """
        R = self.base_field()
        if R in Fields():
            return "[%s, %s] linear rank metric code over GF(%s)"%(self.length(), self.dimension(), R.cardinality())
        else:
            return "[%s, %s] linear rank metric code over %s"%(self.length(), self.dimension(), R)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.
        """
        return "[%s, %s]\\textnormal{ Linear rank metric code over }%s"\
                % (self.length(), self.dimension(), self.base_field()._latex_())

    def __hash__(self):
        r"""
        Returns the hash value of ``self``.
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
        """
        if encoder_name is None or encoder_name is 'GeneratorMatrix':
            g = self._generator_matrix
        else:
            g = super(LinearRankMetricCode, self).generator_matrix(encoder_name, **kwargs)
        g.set_immutable()
        return g

    @cached_method
    def parity_check_matrix(self):
        r"""
        Returns the parity check matrix of ``self``.

        The parity check matrix of a linear rank metric code `C` corresponds to
        the generator matrix of the dual code of `C`.
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
        """
        return self.parity_check_matrix()*r

    def __contains__(self, v):
        r"""
        Returns True if `v` can be coerced into `self`. Otherwise, returns False.
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
        """
        try:
            self.encoder("Systematic", systematic_positions=tuple(positions))
            return True
        except ValueError:
            return False

####################### encoders ###############################
class LinearRankMetricCodeGeneratorMatrixEncoder(Encoder):
    r"""
    Encoder based on generator_matrix for linear rank metric codes.

    This is the default encoder of a generic linear rank metriccode, and should
    never be used for other codes than :class:`LinearRankMetricCode`.

    INPUT:

    - ``code`` -- The associated :class:`LinearRankMetricCode` of this encoder.
    """

    def __init__(self, code):
        r"""
        """
        super(LinearRankMetricCodeGeneratorMatrixEncoder, self).__init__(code)

    def __eq__(self, other):
        r"""
        Tests equality between LinearRankMetricCodeGeneratorMatrixEncoder objects.
        """
        return isinstance(other, LinearRankMetricCodeGeneratorMatrixEncoder)\
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.
        """
        return "Generator matrix-based encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.
        """
        return "\\textnormal{Generator matrix-based encoder for }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.
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
    One can refer to :meth:`LinearCodeSystematicEncoder.generator_matrix`
    for a concrete example.
    Once such a matrix has been computed, it is used to encode any message
    into a codeword.

    This encoder can also serve as the default encoder of a code defined by a
    parity check matrix: if the :class:`LinearCodeSystematicEncoder` detects
    that it is the default encoder, it computes a generator matrix as the
    reduced row echelon form of the right kernel of the parity check matrix.

    INPUT:

    - ``code`` -- The associated code of this encoder.

    - ``systematic_positions`` -- (default: ``None``) the positions in codewords that
      should correspond to the message symbols. A list of `k` distinct integers in
      the range 0 to `n-1` where `n` is the length of the code and `k` its
      dimension. The 0th symbol of a message will then be at position
      ``systematic_positions[0]``, the 1st index at position
      ``systematic_positions[1]``, etc. A ``ValueError`` is raised at
      construction time if the supplied indices do not form an information set.
    """

    def __init__(self, code, systematic_positions=None):
        r"""
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
        """
        return isinstance(other, LinearRankMetricCodeSystematicEncoder)\
                and self.code() == other.code()\
                and self.systematic_positions() == other.systematic_positions()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        return "Systematic encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.
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
        """
        return self._systematic_positions if self._systematic_positions else self.generator_matrix().pivots()

####################### decoders ###############################
class LinearRankMetricCodeNearestNeighborDecoder(Decoder):
    r"""
    Construct a decoder for Linear Rank Metric Codes. This decoder will decode
    to the nearest codeword found.
    """

    def __init__(self, code):
        r"""
        INPUT:

        - ``code`` -- A code associated to this decoder
        """
        super(LinearRankMetricCodeNearestNeighborDecoder, self).__init__(code, code.ambient_space(), \
                code._default_encoder_name)

    def __eq__(self, other):
        r"""
        Tests equality between LinearRankMetricCodeNearestNeighborDecoder objects.
        """
        return isinstance(other, LinearRankMetricCodeNearestNeighborDecoder)\
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.
        """
        return "Nearest neighbor decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.
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
