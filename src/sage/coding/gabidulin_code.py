r"""
Gabidulin Code

This module provides the :class:`~sage.coding.gabidulin.GabidulinCode`, which constructs
Gabidulin Codes that are the rank metric equivalent of Reed Solomon codes and are
defined as the evaluation codes of degree-restricted skew polynomials.

This module also provides :class:`~sage.coding.gabidulin.GabidulinPolynomialEvaluationEncoder`,
an encoder with a skew polynomial message space and :class:`~sage.coding.gabidulin.GabidulinGeneratorMatrixEncoder`,
an encoder based on the generator matrix. It also provides a decoder
:class:`~sage.coding.gabidulin.GabidulinGaoDecoder` which corrects errors using
the Gao algorithm in the rank metric.

AUTHOR:

- Arpit Merchant (2016-08-16)

"""

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.modules.free_module import VectorSpace
from sage.rings.integer import Integer
from encoder import Encoder
from decoder import Decoder, DecodingError
from sage.rings.integer_ring import ZZ
from sage.functions.other import floor
from sage.coding.relative_finite_field_extension import *
from sage.coding.rank_metric import *

class GabidulinCode(AbstractRankMetricCode):
    """
    A Gabidulin Code.

    DEFINITION:

    A linear Gabidulin Code Gab[n, k] over `F_{q^m}` of length `n` (at most
    `m`) and dimension `k` (at most `n`) is the set of all codewords, that
    are the evaluation of a `q`-degree restricted skew polynomial `f(x)`
    belonging to the skew polynomial constructed over the base ring `F_{q^m}`
    and the twisting homomorphism `\sigma`.

    .. math::

        \{ \text{Gab[n, k]} = \big\{ (f(g_0) f(g_1) ... f(g_{n-1})) = f(\textbf{g}) : \text{deg}_{q}f(x) < k \big\} \}

    where the fixed evaluation points `g_0, g_1,..., g_{n-1}` are linearly
    independent over `F_{q^m}`.

    EXAMPLES:

        sage: from sage.coding.gabidulin import *
        sage: Fqm.<aa> = GF(2^9)
        sage: Fq.<a> = GF(8)
        sage: Frob = Fqm.frobenius_endomorphism()
        sage: S.<x> = Fqm['x', Frob]
        sage: C = GabidulinCode(Fqm, Fq, 2, 2, Frob); C
        [2, 2, 1] Linear Gabidulin Code over Finite Field in aa of size 2^9
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, sub_field, length, dimension, \
            twisting_homomorphism, field_extension=None, evaluation_points=None):
        """
        Representation of a Gabidulin Code.

        INPUT:

        - ``base_field`` -- finite field of order `q^m` where `q` is a prime power
          and `m` is an integer

        - ``sub_field`` -- finite field of order `q` which is a subfield of the
          ``base_field``

        - ``length`` -- length of the resulting code

        - ``dimension`` -- dimension of the resulting code

        - ``twisting_homomorphism`` -- homomorphism of the underlying skew polynomial
          ring, the message space of the resulting code

        - ``field_extension`` -- representation of the elements of the relative
          extension of `base_field` over `sub_field` (default: ``None``)

        - ``evaluation_points`` -- list of elements `g_0, g_1,...,g_{n-1}` of the
          `base_field` that are linearly independent over the `sub_field` (default:
          ``None``)

        TESTS:

        A Gabidulin Code can be constructed in the following way:

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: C
            [2, 2, 1] Linear Gabidulin Code over Finite Field in aa of size 2^4

        If the `base_field` is not a finite field, an error is raised:

            sage: Fqm.<aa> = RR
            sage: Fq.<a> = GF(4)
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob); C
            Traceback (most recent call last):
            ...
            ValueError: absolute_field has to be a finite field

        If the `sub_field` is not a finite field, an error is raised:

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = RR
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob); C
            Traceback (most recent call last):
            ...
            ValueError: relative_field has to be a finite field

        If the `sub_field` is not a subfield of `base_field`, an error is
        raised:

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(8)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob); C
            Traceback (most recent call last):
            ...
            ValueError: relative_field has to be a subfield of absolute_field

        If the `length` is not at most power of the `base_field`, an error
        is raised:

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 3, 2, Frob); C
            Traceback (most recent call last):
            ...
            ValueError: length of the code must be a positive integer less than or equal to the base_field_power which is 2

        If the `dimension` is not at most `length`, an error is raised:

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 3, Frob); C
            Traceback (most recent call last):
            ...
            ValueError: dimension of the code must be a positive integer less than or equal to its length which is 2

        If the number of evaluation points is not equal to the length
        of the code, an error is raised:

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5)
            sage: evals = [ aa^i for i in range(21) ]
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob, FE, evals); C
            Traceback (most recent call last):
            ...
            ValueError: the number of evaluation points should be equal to the length of the code

        If evaluation points are not linearly independent over the ``base_field``,
        an error is raised:

            sage: evals = [ aa*i for i in range(2) ]
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob, FE, evals); C
            Traceback (most recent call last):
            ...
            ValueError: the evaluation points provided are not linearly independent
        """
        super(GabidulinCode, self).__init__(base_field, sub_field, \
                length, dimension, "PolynomialEvaluation", "Gao", field_extension)
        if not field_extension:
            field_extension = self.field_extension()

        m = base_field.degree()/sub_field.degree()
        self._m = m
        if not length <= m or length not in ZZ or length < 1:
            raise ValueError("length of the code must be a positive integer less than or equal to the base_field_power which is %d" % m )
        if not dimension <= length or dimension not in ZZ or dimension < 1:
            raise ValueError("dimension of the code must be a positive integer less than or equal to its length which is %d" % length )
        V = VectorSpace(sub_field, m)
        self._vector_space = V
        S = base_field['x', twisting_homomorphism]
        self._message_space = S

        if evaluation_points is None:
            evaluation_points = field_extension.absolute_field_basis()[:length]
        else:
            if not len(evaluation_points) == length:
                raise ValueError("the number of evaluation points should be equal to the length of the code")
            for i in range(length):
                if not evaluation_points[i] in base_field:
                    raise ValueError("evaluation point does not belong to absolute field")
            basis = [field_extension.relative_field_representation(evaluation_points[i]) for i in range(length)]
            if V.linear_dependence(basis):
                raise ValueError("the evaluation points provided are not linearly independent")
        self._evaluation_points = evaluation_points

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob); C
            [2, 2, 1] Linear Gabidulin Code over Finite Field in aa of size 2^4
        """
        return "[%s, %s, %s] Linear Gabidulin Code over %s" \
                % (self.length(), self.dimension(),
                self.minimum_distance(), self.base_field())

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob);
            sage: latex(C)
            [2, 2, 1] \textnormal{ Linear Gabidulin Code over } \Bold{F}_{2^{4}}
        """
        return "[%s, %s, %s] \\textnormal{ Linear Gabidulin Code over } %s"\
                % (self.length(), self.dimension() ,self.minimum_distance(),
                self.base_field()._latex_())

    def __eq__(self, other):
        """
        Tests equality between Gabidulin Code objects.

        INPUT:

        - ``other`` -- another Gabidulin Code object

        OUTPUT:

        Return ``True`` or ``False``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C1 = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: C2 = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: C1.__eq__(C2)
            True

            sage: Fqmm.<aa> = GF(64)
            sage: Frob = Fqmm.frobenius_endomorphism()
            sage: C3 = codes.GabidulinCode(Fqmm, Fq, 2, 2, Frob)
            sage: C3.__eq__(C2)
            False
        """
        return isinstance(other, GabidulinCode) \
                and self.field_extension().absolute_field() == other.field_extension().absolute_field() \
                and self.field_extension().relative_field() == other.field_extension().relative_field() \
                and self.length() == other.length() \
                and self.dimension() == other.dimension() \
                and self.evaluation_points() == other.evaluation_points() \

    def minimum_distance(self):
        """
        Return the minimum distance of ``self``. Since Gabidulin Codes are
        Maximum-Distance-Separable (MDS), this returns ``self.length() -
        self.dimension() + 1``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 20, 15, Frob)
            sage: C.minimum_distance()
            6
        """
        return self.length() - self.dimension() + 1

    # parity_evaluation_points needs to be fixed. dual_code and parity_check_matrix depend on this.
    def parity_evaluation_points(self):
        eval_pts = self.evaluation_points()
        n = self.length()
        k = self.dimension()
        Fqm = self.base_field()
        q = self.field_extension().relative_field().order()
        coefficient_matrix = matrix(Fqm, n - 1, n, lambda i,j: pow(eval_pts[j], pow(q, -n + k + 1 + i))) #rewrite using sigma
        solution_space = coefficient_matrix.right_kernel() #these two lines need to be replaced
        parity_eval_pts = solution_space.random_element()
        return parity_eval_pts

    def dual_code(self):
        parity_eval_pts = self.parity_evaluation_points()
        return GabidulinCode(self.field_extension().absolute_field(),
                self.field_extension().relative_field(),
                self.length(), self.length() - self.dimension(), parity_eval_pts)

    def parity_check_matrix(self):
        E = GabidulinGeneratorMatrixEncoder(self.dual_code())
        return E.generator_matrix()

    def m(self): #rename
        return self._m

    def vector_space(self):
        """
        Return the vector space formed by the ``base_field`` over
        the ``sub_field`` of ``self``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 20, 15, Frob)
            sage: C.vector_space()
            Vector space of dimension 20 over Finite Field of size 5
        """
        return self._vector_space

    def generator_matrix(self):
        """
        Return the generator matrix of ``self``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: C.generator_matrix()
            [   1   aa]
            [   1 aa^2]
        """
        E = GabidulinGeneratorMatrixEncoder(self)
        return E.generator_matrix()

    def evaluation_points(self):
        """
        Return the evaluation points of ``self``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5^4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 4, 4, Frob)
            sage: C.evaluation_points()
            [1, aa, aa^2, aa^3]
        """
        return self._evaluation_points

    def message_space(self):
        """
        Return the message space of ``self``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5^4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 4, 4, Frob)
            sage: C.message_space()
            Skew Polynomial Ring in x over Finite Field in aa of size 5^20 twisted by aa |--> aa^5
        """
        return self._message_space

    def random_element(self):
        """
        Return a random element of ``self``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: C.random_element()
            (aa^3 + aa^2, aa^3)
        """
        E = self.encoder()
        M = self.message_space()
        p = M.random_element()
        while p.degree() >= self.dimension():
            p = M.random_element()
        c = E.encode(p)
        c.set_immutable()
        return c

####################### encoders ###############################


####################### encoders ###############################


class GabidulinPolynomialEvaluationEncoder(Encoder):

    def __init__(self, code):
        """
        This method constructs the encoder for Gabidulin Codes which
        evaluates the skew polynomial at the `n` evaluation points to
        form a codeword of length `n`.

        INPUT:

        - ``code`` -- The associated code of this encoder.

        EXAMPLES::

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: E
            Polynomial evaluation style encoder for [2, 2, 1] Linear Gabidulin Code over Finite Field in aa of size 2^4

        Alternatively, we can construct the encoder from ``C`` directly::

            sage: E = C.encoder("PolynomialEvaluation")
            sage: E
            Polynomial evaluation style encoder for [2, 2, 1] Linear Gabidulin Code over Finite Field in aa of size 2^4
        """
        if not isinstance(code, GabidulinCode):
            raise ValueError("code must be a Gabidulin code")
        super(GabidulinPolynomialEvaluationEncoder, self).__init__(code)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES:

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5^4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 4, 4, Frob)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C); E
            Polynomial evaluation style encoder for [4, 4, 1] Linear Gabidulin Code over Finite Field in aa of size 5^20
        """
        return "Polynomial evaluation style encoder for %s" % self.code()

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES:

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5^4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 4, 4, Frob)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: latex(E)
            \textnormal{Polynomial evaluation style encoder for } [4, 4, 1] \textnormal{ Linear Gabidulin Code over } \Bold{F}_{5^{20}}
        """
        return "\\textnormal{Polynomial evaluation style encoder for } %s" % self.code()._latex_()

    def __eq__(self, other):
        """
        Tests equality between Gabidulin Polynomial Evaluation
        Encoder objects.

        INPUT:

        - ``other`` -- another Gabidulin Polynomial Evaluation Encoder

        OUTPUT:

        Return ``True`` or ``False``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C1 = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: E1 = codes.encoders.GabidulinPolynomialEvaluationEncoder(C1)
            sage: C2 = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: E2 = codes.encoders.GabidulinPolynomialEvaluationEncoder(C2)
            sage: E1.__eq__(E2)
            True

            sage: Fqmm.<aa> = GF(64)
            sage: Frob = Fqmm.frobenius_endomorphism()
            sage: C3 = codes.GabidulinCode(Fqmm, Fq, 2, 2, Frob)
            sage: E3 = codes.encoders.GabidulinPolynomialEvaluationEncoder(C3)
            sage: E3.__eq__(E2)
            False
        """
        return isinstance(other, GabidulinPolynomialEvaluationEncoder) \
                and self.code() == other.code()

    def encode(self, p, form="vector"):
        """
        Transform the polynomial `p` into a codeword of the associated
        Gabidulin Code of ``self``.

        INPUT:

        - ``p`` -- skew polynomial belonging to the message space of
          the associated code of ``self``

        - ``form`` -- type parameter taking strings "vector" or "matrix"
          as values and converts codeword into the respective form
          (default: "vector")

        OUTPUT:

        The codeword corresponding to `p` in the vector or matrix form.

        EXAMPLES:

            sage: Fqm.<aa> = GF(2^9)
            sage: Fq.<a> = GF(2^3)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', Frob]
            sage: p = (aa^6 + aa^2 + aa + 1)*x + aa^7 + aa^5 + aa^4 + aa^2
            sage: codeword_vector = E.encode(p, "vector"); codeword_vector
            (aa^7 + aa^6 + aa^5 + aa^4 + aa + 1, aa^6 + aa^5 + aa^4 + aa^2)
            sage: codeword_matrix = E.encode(p, "matrix"); codeword_matrix
            [    a^2       0]
            [    a^2       0]
            [a^2 + a   a + 1]

        TESTS:

        If the skew polynomial to encode has degree greater than or equal
        to the dimension of the code, an error is raised::

            sage: t = aa^4*x^2 + aa
            sage: codeword_vector = E.encode(t, "vector"); codeword_vector
            Traceback (most recent call last):
            ...
            ValueError: the skew polynomial to encode must have degree at most 1

        The skew polynomial to encode, `p`, must belong to the message
        space of the code. Otherwise, an error is raised::

            sage: Fqmm = GF(2^12)
            sage: frob = Fqmm.frobenius_endomorphism()
            sage: S.<x> = Fqmm['x', frob]
            sage: q = S.random_element()
            sage: codeword_vector = E.encode(q, "vector"); codeword_vector
            Traceback (most recent call last):
            ...
            ValueError: the message to encode must be in Skew Polynomial Ring in x over Finite Field in aa of size 2^9 twisted by aa |--> aa^2
        """
        C = self.code()
        M = C.message_space()
        if p not in M:
            raise ValueError("the message to encode must be in %s" % M)
        if p.degree() >= C.dimension():
            raise ValueError("the skew polynomial to encode must have degree at most %s" % (C.dimension() - 1))
        eval_pts = C.evaluation_points()
        codeword = p.multi_point_evaluation(eval_pts)
        if form == "vector":
            return vector(codeword)
        elif form == "matrix":
            return to_matrix_representation(C, vector(codeword))
        else:
            return ValueError("the argument 'form' takes only either 'vector' or 'matrix' as valid input")

    def unencode_nocheck(self, c):
        """
        Returns the message corresponding to the codeword ``c``.

        Use this method with caution: it does not check if ``c``
        belongs to the code, and if this is not the case, the output is
        unspecified.

        INPUT:

        - ``c`` -- A codeword of :meth:`code`

        OUTPUT:

        A skew polynomial of degree less than ``self.code().dimension()``.

        EXAMPLES:

            sage: Fqm.<aa> = GF(2^9)
            sage: Fq.<a> = GF(2^3)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', Frob]
            sage: p = (aa^6 + aa^4)*x + aa^2 + aa
            sage: codeword_vector = E.encode(p, "vector")
            sage: E.unencode_nocheck(codeword_vector)
            (aa^6 + aa^4)*x + aa^2 + aa
        """
        C = self.code()
        eval_pts = C.evaluation_points()
        values = [c[i] for i in range(len(c))]
        p = C.message_space().interpolation_polynomial(eval_pts, values)
        return p

    def message_space(self):
        """
        Return the message space of the associated code of ``self``.

        EXAMPLES:

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5^4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 4, 4, Frob)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: E.message_space()
            Skew Polynomial Ring in x over Finite Field in aa of size 5^20 twisted by aa |--> aa^5
        """
        return self.code().message_space()

class GabidulinGeneratorMatrixEncoder(Encoder):

    def __init__(self, code):
        """
        This method constructs the generator matrix encoder for
        Gabidulin Codes.

        INPUT:

        - ``code`` -- The associated code of this encoder.

        EXAMPLES::

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: E = codes.encoders.GabidulinGeneratorMatrixEncoder(C)
            sage: E
            Generator matrix style encoder for [2, 2, 1] Linear Gabidulin Code over Finite Field in aa of size 2^4

        Alternatively, we can construct the encoder from ``C`` directly::

            sage: E = C.encoder("GeneratorMatrix")
            sage: E
            Generator matrix style encoder for [2, 2, 1] Linear Gabidulin Code over Finite Field in aa of size 2^4
        """
        if not isinstance(code, GabidulinCode):
            raise ValueError("code must be a Gabidulin code")
        super(GabidulinGeneratorMatrixEncoder, self).__init__(code)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES:

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5^4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 4, 4, Frob)
            sage: E = codes.encoders.GabidulinGeneratorMatrixEncoder(C); E
            Generator matrix style encoder for [4, 4, 1] Linear Gabidulin Code over Finite Field in aa of size 5^20
        """
        return "Generator matrix style encoder for %s" % self.code()

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES:

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5^4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 4, 4, Frob)
            sage: E = codes.encoders.GabidulinGeneratorMatrixEncoder(C)
            sage: latex(E)
            \textnormal{Generator matrix style encoder for } [4, 4, 1] \textnormal{ Linear Gabidulin Code over } \Bold{F}_{5^{20}}
        """
        return "\\textnormal{Generator matrix style encoder for } %s" % self.code()._latex_()

    def __eq__(self, other):
        """
        Tests equality between Gabidulin Generator Matrix
        Encoder objects.

        INPUT:

        - ``other`` -- another Gabidulin Generator Matrix Encoder

        OUTPUT:

        Return ``True`` or ``False``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C1 = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: E1 = codes.encoders.GabidulinGeneratorMatrixEncoder(C1)
            sage: C2 = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: E2 = codes.encoders.GabidulinGeneratorMatrixEncoder(C2)
            sage: E1.__eq__(E2)
            True

            sage: Fqmm.<aa> = GF(64)
            sage: Frob = Fqmm.frobenius_endomorphism()
            sage: C3 = codes.GabidulinCode(Fqmm, Fq, 2, 2, Frob)
            sage: E3 = codes.encoders.GabidulinGeneratorMatrixEncoder(C3)
            sage: E3.__eq__(E2)
            False
        """
        return isinstance(other, GabidulinGeneratorMatrixEncoder) \
                and self.code() == other.code()

    def generator_matrix(self):
        """
        Return the generator matrix of ``self``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(2^9)
            sage: Fq.<a> = GF(2^3)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 3, 3, Frob)
            sage: C.generator_matrix()
            [   1   aa aa^2]
            [   1 aa^2 aa^4]
            [   1 aa^4 aa^8]
        """
        C = self.code()
        eval_pts = C.evaluation_points()
        k = C.dimension()
        sigma = C.message_space().twist_map()
        create_matrix_elements = lambda A,k,f: reduce(lambda L,x: [x] + \
                map(lambda l: map(f,l), L), [A]*k, [])
        return matrix(C.base_field(), C.dimension(), C.length(), \
                create_matrix_elements(eval_pts, C.dimension(), sigma))


####################### decoders ###############################


####################### decoders ###############################


class GabidulinGaoDecoder(Decoder):

    def __init__(self, code):
        """
        This method constructs the Gao style decoder for
        Gabidulin Codes.

        INPUT:

        - ``code`` -- The associated code of this encoder.

        EXAMPLES::

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: D
            Gao decoder for [2, 2, 1] Linear Gabidulin Code over Finite Field in aa of size 2^4

        Alternatively, we can construct the encoder from ``C`` directly::

            sage: D = C.decoder("Gao")
            sage: D
            Gao decoder for [2, 2, 1] Linear Gabidulin Code over Finite Field in aa of size 2^4
        """
        if not isinstance(code, GabidulinCode):
            raise ValueError("code has to be a Gabidulin Code")
        super(GabidulinGaoDecoder, self).__init__(code, code.message_space(), "PolynomialEvaluation")

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES:

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5^4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 4, 4, Frob)
            sage: D = codes.decoders.GabidulinGaoDecoder(C); D
            Gao decoder for [4, 4, 1] Linear Gabidulin Code over Finite Field in aa of size 5^20
        """
        return "Gao decoder for %s" % self.code()

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES:

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5^4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 4, 4, Frob)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: latex(D)
            \textnormal{Gao decoder for } [4, 4, 1] \textnormal{ Linear Gabidulin Code over } \Bold{F}_{5^{20}}
        """
        return "\\textnormal{Gao decoder for } %s" % self.code()._latex_()

    def __eq__(self, other):
        """
        Tests equality between Gabidulin Gao Decoder objects.

        INPUT:

        - ``other`` -- another Gabidulin Gao Decoder

        OUTPUT:

        Return ``True`` or ``False``.

        EXAMPLES::

            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C1 = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: D1 = codes.decoders.GabidulinGaoDecoder(C1)
            sage: C2 = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: D2 = codes.decoders.GabidulinGaoDecoder(C2)
            sage: D1.__eq__(D2)
            True

            sage: Fqmm.<aa> = GF(64)
            sage: Frob = Fqmm.frobenius_endomorphism()
            sage: C3 = codes.GabidulinCode(Fqmm, Fq, 2, 2, Frob)
            sage: D3 = codes.decoders.GabidulinGaoDecoder(C3)
            sage: D3.__eq__(D2)
            False
        """
        return isinstance(other, GabidulinGaoDecoder) \
            and self.code() == other.code()

    def _partial_xgcd(self, a, b, d_stop):
        """
        Compute the partial gcd of `a` and `b` using the right linearized
        extended Euclidean algorithm up to the `d_stop` iterations. This
        is a private method for internal use only.

        INPUT:

        - ``a`` -- a skew polynomial

        - ``b`` -- another skew polynomial

        - ``d_stop`` -- the number of iterations for which the algorithm
          is to be run

        OUTPUT:

        - ``r_c`` -- right linearized remainder of `a` and `b`

        - ``u_c`` -- right linearized quotient of `a` and `b`

        EXAMPLES:

            sage: Fqm.<aa> = GF(2^9)
            sage: Fq.<a> = GF(2^3)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', Frob]
            sage: p = (aa^6 + aa^4)*x + aa^2 + aa
            sage: codeword_vector = E.encode(p, "vector")
            sage: r = D.decode_to_message(codeword_vector) #indirect_doctest
            sage: r
            (aa^6 + aa^4)*x + aa^2 + aa
        """
        C = self.code()
        S = C.message_space()
        if (a not in S) or (b not in S):
            raise ValueError("both the input polynomials must belong to %s" % S)
        if a.degree() < b.degree():
            raise ValueError("degree of first polynomial must be greater than or equal to degree of second polynomial")
        r_p = a
        r_c = b
        u_p = S.zero()
        u_c = S.one()
        v_p = u_c
        v_c = u_p

        while r_c.degree() >= d_stop:
            (q, r_c), r_p = r_p.right_quo_rem(r_c), r_c
            u_c, u_p = u_p - q*u_c, u_c
            v_c, v_p = v_p - q*v_c, v_c
        return r_c, u_c

    def _decode_to_code_and_message(self, r):
        """
        Return the decoded codeword and message (skew polynomial)
        corresponding to the received codeword `r`. This is a
        private method for internal use only.

        INPUT:

        - ``r`` -- received codeword

        OUTPUT:

        The decoded codeword and decoded message corresponding to
        the received codeword `r`.

        EXAMPLES:

            sage: Fqm.<aa> = GF(2^9)
            sage: Fq.<a> = GF(2^3)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', Frob]
            sage: p = (aa^6 + aa^4)*x + aa^2 + aa
            sage: codeword_vector = E.encode(p, "vector")
            sage: r = D.decode_to_message(codeword_vector) #indirect doctest
            sage: r
            (aa^6 + aa^4)*x + aa^2 + aa
        """
        C = self.code()
        length = len(r)
        if not length <= C.m() or length < 1:
            raise ValueError("length of the received code must be a positive integer \
                    less than or equal to the base_field_power which is %d" % m )
        eval_pts = C.evaluation_points()
        S = C.message_space()

        if length == C.dimension() or r in C:
            return r, self.connected_encoder().unencode_nocheck(r)

        R = S.interpolation_polynomial(eval_pts, list(r))
        r_out, u_out = self._partial_xgcd(S.minimal_vanishing_polynomial(eval_pts), \
                R, floor((C.length() + C.dimension())//2))
        quo, rem = r_out.left_quo_rem(u_out)
        if not rem.is_zero():
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        if quo not in S:
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        c = self.connected_encoder().encode(quo)
        if rank_weight(C, (c-r)) > self.decoding_radius():
            raise DecodingError("Decoding failed because the number of errors exceeded the decoding radius")
        return c, quo

    def decode_to_code(self, r):
        """
        Return the decoded codeword corresponding to the
        received word `r`.

        INPUT:

        - ``r`` -- received codeword

        OUTPUT:

        The decoded codeword corresponding to the received codeword.

        EXAMPLES:

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 6, 4, Frob)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', Frob]
            sage: p = 4*aa^11*x^2 + aa^4*x + 2*aa
            sage: codeword_vector = E.encode(p, "vector")
            sage: r = D.decode_to_code(codeword_vector)
            sage: r
            (aa^6 + aa^4)*x + aa^2 + aa
        """
        return self._decode_to_code_and_message(r)[0]

    def decode_to_message(self, r):
        """
        Return the skew polynomial (message) corresponding to the
        received word `r`.

        INPUT:

        - ``r`` -- received codeword

        OUTPUT:

        The message corresponding to the received codeword.

        EXAMPLES:

            sage: Fqm.<aa> = GF(2^9)
            sage: Fq.<a> = GF(2^3)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 2, 2, Frob)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: E = codes.encoders.GabidulinPolynomialEvaluationEncoder(C)
            sage: S.<x> = Fqm['x', Frob]
            sage: p = (aa^6 + aa^4)*x + aa^2 + aa
            sage: codeword_vector = E.encode(p, "vector")
            sage: r = D.decode_to_message(codeword_vector)
            sage: r
            (aa^6 + aa^4)*x + aa^2 + aa
        """
        return self._decode_to_code_and_message(r)[1]

    def decoding_radius(self):
        """
        Return the decoding radius of the Gabidulin Gao Decoder.

        EXAMPLES:

            sage: Fqm.<aa> = GF(5^20)
            sage: Fq.<a> = GF(5)
            sage: Frob = Fqm.frobenius_endomorphism()
            sage: C = codes.GabidulinCode(Fqm, Fq, 20, 4, Frob)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: D.decoding_radius()
            8
        """
        return (self.code().minimum_distance()-1)//2

############################## registration ####################################

GabidulinCode._registered_encoders["PolynomialEvaluation"] = GabidulinPolynomialEvaluationEncoder
GabidulinCode._registered_encoders["GeneratorMatrix"] = GabidulinGeneratorMatrixEncoder

GabidulinCode._registered_decoders["Gao"] = GabidulinGaoDecoder
