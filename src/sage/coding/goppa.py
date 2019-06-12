from sage.coding.linear_code import AbstractLinearCode
from sage.coding.encoder import Encoder
from sage.coding.decoder import Decoder
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.modules.free_module_element import vector
from sage.coding.all import codes

def _columnize(element):
    v = vector(element)
    return v.column()

class GoppaCode(AbstractLinearCode):
    """Implementation of Goppa codes, a generalization of narrow-sense BCH codes.

    INPUTS:

    -generating_pol -- a monic polynomial with coefficients in a finite field F_p^m,
    the code is defined over F_p, p must be a prime number

    -defining_set -- a set of elements of F_p^m that are not roots of g,
    its cardinality is the length of the code

    EXAMPLES::

        sage: from sage.coding.goppa import GoppaCode
        sage: F8=GF(2**6)
        sage: RR.<xx> = F8[]
        sage: g = xx^9+1
        sage: L = [a for a in F8.list() if g(a)!=0]
        sage: C = GoppaCode(g, L)
        sage: C
        [55, 16] Goppa code
    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, generating_pol, defining_set):

        self._field = generating_pol.base_ring()
        self._length = len(defining_set)
        self._generating_pol = generating_pol
        self._defining_set = defining_set

        super(GoppaCode, self).__init__(self._field, self._length, "GoppaEncoder", "Syndrome")

        if not generating_pol.is_monic():
            raise ValueError("generating_pol must be monic")
        F = self._field
        if (not F.is_field() or not F.is_finite()):
            raise ValueError("generating_pol must be defined over a finite field")
        for a in defining_set:
            if generating_pol(a) == 0:
                raise ValueError("Defining elements cannot be roots of generating polynomial")



    def _repr_(self):
        """Representation of a Goppa code.

        EXAMPLES::

            sage: from sage.coding.goppa import GoppaCode
            sage: F8=GF(2**3)
            sage: RR.<xx> = F8[]
            sage: g = xx^2+xx+1
            sage: L = [a for a in F8.list() if g(a)!=0]
            sage: C = GoppaCode(g, L)
            sage: C._repr_()
            '[8, 2] Goppa code'

        """
        return "[%s, %s] Goppa code" %(self.length(), self.dimension())

    def _latex_(self):
        """Returns latex representation of self.

        EXAMPLES::

            sage: from sage.coding.goppa import GoppaCode
            sage: F8=GF(2**3)
            sage: RR.<xx> = F8[]
            sage: g = xx^2+xx+1
            sage: L = [a for a in F8.list() if g(a)!=0]
            sage: C = GoppaCode(g, L)
            sage: C._latex_()
            '\\textnormal{Goppa code of length } 8'

        """

        return ("\\textnormal{Goppa code of length } %s" % self.length())

    def __eq__(self, other):
        """Equality check between GoppaCode objects.

        EXAMPLES::

            sage: from sage.coding.goppa import GoppaCode
            sage: F8=GF(2**3)
            sage: RR.<xx> = F8[]
            sage: g = xx^2+xx+1
            sage: L = [a for a in F8.list() if g(a)!=0]
            sage: C = GoppaCode(g, L)
            sage: D = GoppaCode(g, L)
            sage: D.__eq__(C)
            True

        """
        return (isinstance(other, GoppaCode)
           and self.length() == other.length()
           and self._generating_pol == other._generating_pol
           and self._defining_set == other._defining_set)

    def parity_check_matrix(self):
        """Returns a parity check matrix for 'self'.

        The element in row t, column i is h[i]*(D[i]**t), where:

        -h[i] is the inverse of g(D[i])
        -D[i] is the ith element of the defining set

        In the resulting d * n matrix we interpret each entry as an m-column vector
        and return a dm * n matrix.


        EXAMPLES::

            sage: from sage.coding.goppa import GoppaCode
            sage: F8=GF(2**3)
            sage: RR.<xx> = F8[]
            sage: g = xx^2+xx+1
            sage: L = [a for a in F8.list() if g(a)!=0]
            sage: C = GoppaCode(g, L)
            sage: C
            [8, 2] Goppa code
            sage: C.parity_check_matrix()
            [1 0 0 0 0 0 0 1]
            [0 0 1 0 1 1 1 0]
            [0 1 1 1 0 0 1 0]
            [0 1 1 1 1 1 1 1]
            [0 1 0 1 1 0 1 0]
            [0 0 1 1 1 1 0 0]

        """

        g = self._generating_pol
        F = g.base_ring()
        m = self.base_field().degree()
        n = self._length
        d = g.degree()
        alpha = F.primitive_element()

        D = self._defining_set
        h = [(g(D[i]).inverse_of_unit()) for i in range(n)]

        #assemble top row
        M = _columnize(alpha)
        for i in range(n):
            v = _columnize(h[i])
            M = M.augment(v)
        M = M.delete_columns([0])
        old = M

        for t in range(1,d):
            #assemble row
            M = _columnize(alpha)
            for i in range(n):
                v = _columnize(h[i]*(D[i]**t))
                M = M.augment(v)
            M = M.delete_columns([0])
            new = M
            old = old.stack(new)

        return old

    def _parity_check_matrix_Vandermonde(self):
        """Alternative way to generate parity check matrix for Goppa codes using
        Vandermonde matrix.

        EXAMPLES::

            sage: from sage.coding.goppa import GoppaCode
            sage: F8=GF(2**3)
            sage: RR.<xx> = F8[]
            sage: g = xx^2+xx+1
            sage: L = [a for a in F8.list() if g(a)!=0]
            sage: C = GoppaCode(g, L)
            sage: C
            [8, 2] Goppa code
            sage: C._parity_check_matrix_Vandermonde() == C.parity_check_matrix()
            True

        """
        L = self._defining_set
        g = self._generating_pol
        t = g.degree()

        from sage.matrix.constructor import vandermonde, matrix, diagonal_matrix, block_matrix

        V = matrix.vandermonde(L)
        V = V.transpose()
        GL = [g(i) for i in L]
        GLI = [j.inverse_of_unit() for j in GL]
        D = diagonal_matrix(GLI)
        VF = matrix([V.row(i) for i in range(t)])
        H = VF*D

        matrices = [matrix([vector(i) for i in H.row(j)]) for j in range(t)]
        matrices = [m.transpose() for m in matrices]

        m = block_matrix(t, 1, matrices)
        return m

    def distance_bound(self):
        """Estimates a lower bound for the minimun distance of the code using the degree of g.
        Min distance is guarenteed to be larger or equal to this bound.

        EXAMPLES::

            sage: from sage.coding.goppa import GoppaCode
            sage: F8=GF(2**3)
            sage: RR.<xx> = F8[]
            sage: g = xx^2+xx+1
            sage: L = [a for a in F8.list() if g(a)!=0]
            sage: C = GoppaCode(g, L)
            sage: C
            [8, 2] Goppa code
            sage: C.distance_bound()
            3
            sage: C.minimum_distance()
            5

        """
        return 1 + (self._generating_pol).degree()


class GoppaCodeEncoder(Encoder):

    def __init__(self, code):
        super(GoppaCodeEncoder, self).__init__(code)

    def _repr_(self):
        return "Goppa encoder for the %s" % self.code()

    def _latex_(self):
        return "\textnormal{Goppa encoder for the } %s" % self.code()

    def __eq__(self, other):
        return (isinstance(other, GoppaCodeEncoder)
           and self.code() == other.code())

    def generator_matrix(self):
        """Returns a generator matrix for 'self'.

        Dimension of resulting matrix is d * n?


        EXAMPLES::

            sage: from sage.coding.goppa import GoppaCode
            sage: F8=GF(2**3)
            sage: RR.<xx> = F8[]
            sage: g = xx^2+xx+1
            sage: L = [a for a in F8.list() if g(a)!=0]
            sage: C = GoppaCode(g, L)
            sage: C
            [8, 2] Goppa code
            sage: C.generator_matrix()
            [1 0 0 1 0 1 1 1]
            [0 1 1 1 1 1 1 0]

        """
        c = self.code()
        pmat = c.parity_check_matrix()
        aux = codes.from_parity_check_matrix(pmat)
        return aux.generator_matrix()

GoppaCode._registered_encoders["GoppaEncoder"] = GoppaCodeEncoder
