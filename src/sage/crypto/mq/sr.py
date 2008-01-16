r"""
Small Scale Variants of the AES (SR) Polynomial System Generator.

We support polynomial system generation over $GF(2)$ and
$GF(2^e)$. Also, we support both the specification of SR as given in
the paper and we support a variant of SR* which is AES.

AUTHORS:
    -- Martin Albrecht <malb@informatik.uni-bremen.de> (2007-09) initial version

EXAMPLES:
    sage: sr = mq.SR(1,1,1,4)

    $n$ is the number of rounds, $r$ the number of rows in the state
    array, $c$ the number of columns in the state array. $e$ the
    degree of the underlying field.

    sage: sr.n, sr.r, sr.c, sr.e
    (1, 1, 1, 4)

    By default variables are ordered reverse to as they appear., e.g.:

    sage: sr.R
    Multivariate Polynomial Ring in k100, k101, k102, k103, x100, x101, x102, x103, w100, w101, w102, w103, s000, s001, s002, s003, k000, k001, k002, k003 over Finite Field in a of size 2^4

    For SR(1,1,1,4) the ShiftRows matrix isn't that interresting:

    sage: sr.ShiftRows
    [1 0 0 0]
    [0 1 0 0]
    [0 0 1 0]
    [0 0 0 1]

    Also, the MixColumns matrix is the identity matrix:

    sage: sr.MixColumns
    [1 0 0 0]
    [0 1 0 0]
    [0 0 1 0]
    [0 0 0 1]

    Lin, however is not the identity matrix.

    sage: sr.Lin
    [          a^2 + 1                 1         a^3 + a^2           a^2 + 1]
    [                a                 a                 1 a^3 + a^2 + a + 1]
    [          a^3 + a               a^2               a^2                 1]
    [                1               a^3             a + 1             a + 1]


    M and Mstar are identical for SR(1,1,1,4)

    sage: sr.M
    [          a^2 + 1                 1         a^3 + a^2           a^2 + 1]
    [                a                 a                 1 a^3 + a^2 + a + 1]
    [          a^3 + a               a^2               a^2                 1]
    [                1               a^3             a + 1             a + 1]

    sage: sr.Mstar
    [          a^2 + 1                 1         a^3 + a^2           a^2 + 1]
    [                a                 a                 1 a^3 + a^2 + a + 1]
    [          a^3 + a               a^2               a^2                 1]
    [                1               a^3             a + 1             a + 1]


TESTS:
    sage: sr == loads(dumps(sr))
    True

REFERENCES:
   C. Cid , S. Murphy, and M.J.B. Robshaw; Small Scale Variants of the
   AES; in Proceedings of Fast Software Encryption 2005, LNCS 3557;
   Springer 2005; available at http://www.isg.rhul.ac.uk/~sean/smallAES-fse05.pdf

   C. Cid , S. Murphy, and M.J.B. Robshaw; Algebraic Aspects of the
   Advanced Encryption Standard; Springer 2006;
"""

from sage.rings.finite_field import GF
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import PolynomialRing

from sage.matrix.matrix import is_Matrix
from sage.matrix.constructor import Matrix, random_matrix
from sage.matrix.matrix_space import MatrixSpace

from sage.misc.misc import get_verbose, set_verbose
from sage.misc.flatten import flatten

from sage.modules.vector_modn_dense import Vector_modn_dense

from mpolynomialsystem import MPolynomialSystem, MPolynomialRoundSystem
from mpolynomialsystemgenerator import MPolynomialSystemGenerator

from sage.rings.polynomial.term_order import TermOrder

def SR(n=1,r=1,c=1,e=4, star=False, **kwargs):
    """
    Return a small scale variant of the AES polynomial system
    constructor subject to the following conditions:

    INPUT:
        n -- the number of rounds (default: 1)
        r -- the number of rows in the state array (default: 1)
        c -- the number of columns in the state array (default: 1)
        e -- the exponent of the finite extension field (default: 4)
        star -- determines if SR* or SR should be constructed (default: False)
        aes_mode -- as the SR key schedule specification differs slightly from the AES
                    key schedule this parameter controls which schedule to use (default: True)
        gf2 -- generate polynomial systems over GF(2) rather than over GF(2^n) (default: False)
        order -- a string to specify the term ordering of the variables
        postfix -- a string which is appended after the variable name (default: '')
        allow_zero_inversions -- a boolean to controll whether zero inversions raise
                                 an exception (default: False)

    ADDITIONAL INPUT FOR GF(2):
        correct_only -- only include correct inversion polynomials (default: False)
        biaffine_only -- only include bilinear and biaffine inversion polynomials (default: True)
    """
    if not kwargs.get("gf2",False):
        return SR_gf2n(n,r,c,e,star,**kwargs)
    else:
        return SR_gf2(n,r,c,e,star,**kwargs)

class SR_generic(MPolynomialSystemGenerator):
    """
    Small Scale Variants of the AES.

    EXAMPLES:
        sage: sr = mq.SR(1,1,1,4)
        sage: ShiftRows = sr.shift_rows_matrix()
        sage: MixColumns = sr.mix_columns_matrix()
        sage: Lin = sr.lin_matrix()
        sage: M = MixColumns * ShiftRows * Lin
        sage: print sr.hex_str_matrix(M)
         5 1 C 5
         2 2 1 F
         A 4 4 1
         1 8 3 3

    TESTS:
        sage: sr = mq.SR(1,2,1,4)
        sage: ShiftRows = sr.shift_rows_matrix()
        sage: MixColumns = sr.mix_columns_matrix()
        sage: Lin = sr.lin_matrix()
        sage: M = MixColumns * ShiftRows * Lin
        sage: print sr.hex_str_matrix(M)
         F 3 7 F A 2 B A
         A A 5 6 8 8 4 9
         7 8 8 2 D C C 3
         4 6 C C 5 E F F
         A 2 B A F 3 7 F
         8 8 4 9 A A 5 6
         D C C 3 7 8 8 2
         5 E F F 4 6 C C

        sage: sr = mq.SR(1,2,2,4)
        sage: ShiftRows = sr.shift_rows_matrix()
        sage: MixColumns = sr.mix_columns_matrix()
        sage: Lin = sr.lin_matrix()
        sage: M = MixColumns * ShiftRows * Lin
        sage: print sr.hex_str_matrix(M)
         F 3 7 F 0 0 0 0 0 0 0 0 A 2 B A
         A A 5 6 0 0 0 0 0 0 0 0 8 8 4 9
         7 8 8 2 0 0 0 0 0 0 0 0 D C C 3
         4 6 C C 0 0 0 0 0 0 0 0 5 E F F
         A 2 B A 0 0 0 0 0 0 0 0 F 3 7 F
         8 8 4 9 0 0 0 0 0 0 0 0 A A 5 6
         D C C 3 0 0 0 0 0 0 0 0 7 8 8 2
         5 E F F 0 0 0 0 0 0 0 0 4 6 C C
         0 0 0 0 A 2 B A F 3 7 F 0 0 0 0
         0 0 0 0 8 8 4 9 A A 5 6 0 0 0 0
         0 0 0 0 D C C 3 7 8 8 2 0 0 0 0
         0 0 0 0 5 E F F 4 6 C C 0 0 0 0
         0 0 0 0 F 3 7 F A 2 B A 0 0 0 0
         0 0 0 0 A A 5 6 8 8 4 9 0 0 0 0
         0 0 0 0 7 8 8 2 D C C 3 0 0 0 0
         0 0 0 0 4 6 C C 5 E F F 0 0 0 0

    """
    def __init__(self,n=1,r=1,c=1,e=4, star=False, **kwargs):
        """
        See help for SR.
        """
        if n-1 not in range(10):
            raise TypeError, "n must be between 1 and 10 (inclusive)"
        self._n = n

        if r not in (1,2,4):
            raise TypeError, "r must be in (1,2,4)"
        self._r = r

        if c not in (1,2,4):
            raise TypeError, "c must be in (1,2,4)"
        self._c = c

        if e not in (4,8):
            raise TypeError, "e must be either 4 or 8"
        self._e = e

        self._star = bool(star)

        self._base = self.base_ring()

        # to generate table, allow zero inversions
        self._allow_zero_inversions = True
        sub_byte_lookup = dict([(e,self.sub_byte(e)) for e in self._base])
        self._sub_byte_lookup = sub_byte_lookup

        self._postfix = kwargs.get("postfix", "")
        self._order = kwargs.get("order", "degrevlex")
        self._allow_zero_inversions= bool(kwargs.get("allow_zero_inversions", False))
        self._aes_mode = kwargs.get("aes_mode",True)
        self._gf2 = kwargs.get("gf2",False)

    def __getattr__(self, attr):
        if attr == "e":
            return self._e
        elif attr == "c":
            return self._c
        elif attr == "n":
            return self._n
        elif attr == "r":
            return self._r

        elif attr == "R":
            self.R = self.ring()
            return self.R
        elif attr == "k":
            self.k = self.base_ring()
            return self.k

        elif attr == "M":
            self.M = self.MixColumns * self.ShiftRows * self.Lin
            return self.M
        elif attr == "Mstar":
            self.Mstar = self.MixColumns * self.ShiftRows * self.Lin
            return self.Mstar
        elif attr == "ShiftRows":
            self.ShiftRows = self.shift_rows_matrix()
            return self.ShiftRows
        elif attr == "MixColumns":
            self.MixColumns = self.mix_columns_matrix()
            return self.MixColumns
        elif attr == "Lin":
            self.Lin = self.lin_matrix()
            return self.Lin

        raise AttributeError, "%s has no attribute %s"%(type(self),attr)

    def _repr_(self):
        if self._star:
            return "SR*(%d,%d,%d,%d)"%(self._n,self._r,self._c,self._e)
        else:
            return "SR(%d,%d,%d,%d)"%(self._n,self._r,self._c,self._e)

    def base_ring(self):
        """
        Return the base field of self as determined through self.e.

        EXAMPLE:
            sage: sr = mq.SR(10,2,2,4)
            sage: sr.base_ring().polynomial()
            a^4 + a + 1

            The Rijndael polynomial:

            sage: sr = mq.SR(10,4,4,8)
            sage: sr.base_ring().polynomial()
            a^8 + a^4 + a^3 + a + 1
        """
        try:
            return self._base
        except AttributeError:
            pass

        if self._e == 4:
            self._base = GF(2**4,'a',modulus=(1,1,0,0,1))
        elif self._e == 8:
            self._base = GF(2**8,'a',modulus=(1,1,0,1,1,0,0,0,1))

        return self._base

    def __cmp__(self, other):
        """
        The internal state dictionaries are compared.
        """
        return cmp(self.__dict__,other.__dict__)

    def sub_bytes(self, d):
        """
        Perform the non-linear transform on d.

        INPUT:
            d -- state array or something coercable to a state array
        """
        d = self.state_array(d)
        return Matrix(self.base_ring(), d.nrows(), d.ncols(), [self.sub_byte(b) for b in d.list()])

    def sub_byte(self, b):
        """
        Perform SubByte on a single byte/halfbyte b.

        A ZeroDivision exception is raised if an attempt is made to
        perform an inversion on the zero element. This can be disabled
        by passing allow_zero_inversion=True to the constructor. A zero inversion
        will result in an inconsisten equation system.

        INPUT:
            b -- an element in self.base_ring()

        EXAMPLE:

           The S-Box table for GF(2^4):

           sage: sr = mq.SR(1,1,1,4, allow_zero_inversions=True)
           sage: for e in sr.base_ring():
           ...    print '% 20s % 20s'%(e, sr.sub_byte(e))
                           0              a^2 + a
                           a              a^2 + 1
                         a^2                    a
                         a^3              a^3 + 1
                       a + 1                  a^2
                     a^2 + a          a^2 + a + 1
                   a^3 + a^2                a + 1
                 a^3 + a + 1            a^3 + a^2
                     a^2 + 1        a^3 + a^2 + a
                     a^3 + a    a^3 + a^2 + a + 1
                 a^2 + a + 1              a^3 + a
               a^3 + a^2 + a                    0
           a^3 + a^2 + a + 1                  a^3
               a^3 + a^2 + 1                    1
                     a^3 + 1        a^3 + a^2 + 1
                           1          a^3 + a + 1

        """
        if not b:
            if not self._allow_zero_inversions:
                raise ZeroDivisionError,  "A zero inversion occurred during an encryption or key schedule."
            else:
                return self.sbox_constant()
        try:
            return self._sub_byte_lookup[b]
        except AttributeError:
            pass

        e = self.e
        k = self.k
        a = k.gen()

        # inversion
        b = b ** ( 2**e - 2 )

        # GF(2) linear map
        if e == 4:
            if not hasattr(self,"_L"):
                self._L = Matrix(GF(2),4,4,[[1, 1, 1, 0],
                                            [0, 1, 1, 1],
                                            [1, 0, 1, 1],
                                            [1, 1, 0, 1]])

        elif e==8:
            if not hasattr(self,"_L"):
                self._L = Matrix(GF(2),8,8,[[1, 0, 0, 0, 1, 1, 1, 1],
                                            [1, 1, 0, 0, 0, 1, 1, 1],
                                            [1, 1, 1, 0, 0, 0, 1, 1],
                                            [1, 1, 1, 1, 0, 0, 0, 1],
                                            [1, 1, 1, 1, 1, 0, 0, 0],
                                            [0, 1, 1, 1, 1, 1, 0, 0],
                                            [0, 0, 1, 1, 1, 1, 1, 0],
                                            [0, 0, 0, 1, 1, 1, 1, 1]])

        b = k(self._L * b.vector())

        # constant addition
        if e == 4:
            b = b + k.fetch_int(6)
        elif e == 8:
            b = b + k.fetch_int(99)

        return b

    def sbox_constant(self):
        """
        Return the sbox constant which is added after $L(x^{-1})$ was
        performed. That is 0x63 if e == 8 or 0x6 if e == 4.

        EXAMPLE:
            sage: sr = mq.SR(10,1,1,8)
            sage: sr.sbox_constant()
            a^6 + a^5 + a + 1
        """
        k = self.k
        if self.e == 4:
            return k.fetch_int(6)
        elif self.e == 8:
            return k.fetch_int(99)
        else:
            raise TypeError, "sbox constant only defined for e in (4,8)"

    def shift_rows(self, d):
        """
        Perform the ShiftRows operation on d.

        INPUT:
            d -- state arrary or something coercable to an state array

        EXAMPLES:
            sage: sr = mq.SR(10,4,4,4)
            sage: E = sr.state_array() + 1; E
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: sr.shift_rows(E)
            [1 0 0 0]
            [1 0 0 0]
            [1 0 0 0]
            [1 0 0 0]
        """
        d = self.state_array(d)
        ret = []
        for i in range(d.nrows()):
            ret += list(d.row(i)[i%d.ncols():]) + list(d.row(i)[:i%d.ncols()])
        return Matrix(self.base_ring(), self._r, self._c, ret)

    def mix_columns(self, d):
        """
        Perform the MixColumns operation on d.

        INPUT:
            d -- state arrary or something coercable to an state array

        EXAMPLES:
            sage: sr = mq.SR(10,4,4,4)
            sage: E = sr.state_array() + 1; E
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

            sage: sr.mix_columns(E)
            [    a a + 1     1     1]
            [    1     a a + 1     1]
            [    1     1     a a + 1]
            [a + 1     1     1     a]
        """
        d = self.state_array(d)
        k = self.base_ring()
        a = k.gen()
        r = self._r
        if r == 1:
            M = Matrix(self.base_ring(),1,1, [[1]])
        elif r == 2:
            M = Matrix(self.base_ring(),2,2, [[a + 1, a],
                                              [a, a + 1]])

        elif r == 4:
            M = Matrix(self.base_ring(),4,4, [[a, a+1, 1, 1],
                                              [1, a, a+1, 1],
                                              [1, 1, a, a+1],
                                              [a+1, 1, 1, a]])
        ret =[]
        for column in d.columns():
            ret.append(M * column)
        # AES uses the column major ordering
        return Matrix(k, d.ncols(), d.nrows(), ret).transpose()


    def add_round_key(self, d, key):
        """
        Perform the AddRoundKey operation on d using key.

        INPUT:
            d -- state arrary or something coercable to an state array
            key -- state arrary or something coercable to an state array

        EXAMPLE:
            sage: sr = mq.SR(10,4,4,4)
            sage: D = sr.random_state_array()
            sage: K = sr.random_state_array()
            sage: sr.add_round_key(D,K) == K + D
            True
        """
        d = self.state_array(d)
        key = self.state_array(key)

        return d+key

    def state_array(self, d=None):
        """
        Convert the parameter to a state array.

        INPUT:
            d -- None, a matrix, a list, or a tuple

        EXAMPLES:
            sage: sr = mq.SR(2,2,2,4)
            sage: k = sr.base_ring()
            sage: e1 = [k.fetch_int(e) for e in range(2*2)]; e1
            [0, 1, a, a + 1]
            sage: e2 = sr.phi( Matrix(k,2*2,1,e1) )
            sage: sr.state_array(e1) # note the column major ordering
            [    0     a]
            [    1 a + 1]
            sage: sr.state_array(e2)
            [    0     a]
            [    1 a + 1]

            sage: sr.state_array()
            [0 0]
            [0 0]

        """
        r = self.r
        c = self.c
        e = self.e
        k = self.base_ring()

        if d is None:
            return Matrix(k, r, c)

        if is_Matrix(d):
            if d.nrows() == r*c*e:
                return Matrix(k, c, r, self.antiphi(d)).transpose()
            elif d.ncols() == c and d.nrows() == r and d.base_ring() == k:
                return d

        if isinstance(d, tuple([list,tuple])):
            return Matrix(k, c, r, d).transpose()

    def is_state_array(self, d):
        """
        Return True if d is a state array, i.e. has the correct dimensions and base field.
        """
        return is_Matrix(d) and \
               d.nrows() == self.r and \
               d.ncols() == self.c and \
               d.base_ring() == self.base_ring()

    def random_state_array(self):
        """
        Return a random element in MatrixSpace(self.base_ring(),self.r,self.c).
        """
        return random_matrix(self.base_ring(),self._r,self._c)

    def random_vector(self):
        """
        Return a random vector as it might appear in the algebraic
        expression of self.

        NOTE: Phi was already applied to the result.
        """
        return self.vector( self.random_state_array() )

    def random_element(self, type = "vector"):
        """
        Return a random element for self.

        INPUT:
            type -- either 'vector' or 'state array' (default: 'vector')

        EXAMPLE:
            sage: sr = mq.SR()
            sage: sr.random_element() # random
            [      a]
            [    a^2]
            [  a + 1]
            [a^2 + 1]
            sage: sr.random_element('state_array') # random
            [a^3 + a + 1]
        """
        if type == "vector":
            return self.random_vector()
        elif type == "state_array":
            return self.random_state_array()
        else:
            raise TypeError, "parameter type not understood"

    def key_schedule(self, kj, i):
        """
        Return $k_i$ for a given $i$ and $k_j$ with $j = i-1$.

        TESTS:
            sage: sr = mq.SR(10,4,4,8, star=True, allow_zero_inversions=True)
            sage: ki = sr.state_array()
            sage: for i in range(10):
            ...  ki = sr.key_schedule(ki,i+1)
            sage: print sr.hex_str_matrix(ki)
            B4 3E 23 6F
            EF 92 E9 8F
            5B E2 51 18
            CB 11 CF 8E
        """
        if i < 0:
            raise TypeError, "i must be >= i"

        if i == 0:
            return kj

        r = self.r
        c = self.c
        e = self.e
        F = self.base_ring()
        a = F.gen()
        SubByte = self.sub_byte

        rc = Matrix(F,r,c, ([a**(i-1)] * c) + [F(0)]*((r-1)*c) )
        ki = Matrix(F,r,c)

        if r == 1:
            s0 = SubByte(kj[0,c-1])

            if c > 1:
                for q in range(c):
                    ki[0,q] = s0 + sum([kj[0,t] for t in range(q+1) ])
            else:
                ki[0,0] = s0

        elif r == 2:
            s0 = SubByte(kj[1,c-1])
            s1 = SubByte(kj[0,c-1])

            if c > 1:
                for q in range(c):
                    ki[0, q] = s0 + sum([ kj[0,t] for t in range(q+1) ])
                    ki[1, q] = s1 + sum([ kj[1,t] for t in range(q+1) ])
            else:
                ki[0,0] = s0
                ki[1,0] = s1

        elif r == 4:

            if self._aes_mode:
                s0 = SubByte(kj[1,c-1])
                s1 = SubByte(kj[2,c-1])
                s2 = SubByte(kj[3,c-1])
                s3 = SubByte(kj[0,c-1])
            else:
                s0 = SubByte(kj[3,c-1])
                s1 = SubByte(kj[2,c-1])
                s2 = SubByte(kj[1,c-1])
                s3 = SubByte(kj[0,c-1])

            if c > 1:
                for q in range(c):
                    ki[0,q] = s0 + sum([ kj[0,t] for t in range(q+1) ])
                    ki[1,q] = s1 + sum([ kj[1,t] for t in range(q+1) ])
                    ki[2,q] = s2 + sum([ kj[2,t] for t in range(q+1) ])
                    ki[3,q] = s3 + sum([ kj[3,t] for t in range(q+1) ])

            else:
                ki[0,0] = s0
                ki[1,0] = s1
                ki[2,0] = s2
                ki[3,0] = s3
        ki += rc

        return ki

    def __call__(self, P, K):
        r"""
        Encrypts the plaintext $P$ using the key $K$.

        Both must be given as state arrays or coercable to state arrays.

        INPUTS:
            P -- plaintext as state array or something coercable to a state array
            K -- key as state array or something coercable to a state array

        TESTS:
            The official AES test vectors:

            sage: sr = mq.SR(10,4,4,8, star=True, allow_zero_inversions=True)
            sage: k = sr.base_ring()
            sage: plaintext = sr.state_array([k.fetch_int(e) for e in range(16)])
            sage: key = sr.state_array([k.fetch_int(e) for e in range(16)])
            sage: print sr.hex_str_matrix( sr(plaintext, key) )
            0A 41 F1 C6
            94 6E C3 53
            0B F0 94 EA
            B5 45 58 5A

            Brian Gladman's development vectors (dev_vec.txt):

            sage: sr = mq.SR(10,4,4,8,star=True,allow_zero_inversions=True, aes_mode=True)
            sage: k = sr.base_ring()
            sage: plain = '3243f6a8885a308d313198a2e0370734'
            sage: key = '2b7e151628aed2a6abf7158809cf4f3c'
            sage: set_verbose(2)
            sage: cipher = sr(plain,key)
            R[01].start   193DE3BEA0F4E22B9AC68D2AE9F84808
            R[01].s_box   D42711AEE0BF98F1B8B45DE51E415230
            R[01].s_row   D4BF5D30E0B452AEB84111F11E2798E5
            R[01].m_col   046681E5E0CB199A48F8D37A2806264C
            R[01].k_sch   A0FAFE1788542CB123A339392A6C7605
            R[02].start   A49C7FF2689F352B6B5BEA43026A5049
            R[02].s_box   49DED28945DB96F17F39871A7702533B
            R[02].s_row   49DB873B453953897F02D2F177DE961A
            R[02].m_col   584DCAF11B4B5AACDBE7CAA81B6BB0E5
            R[02].k_sch   F2C295F27A96B9435935807A7359F67F
            R[03].start   AA8F5F0361DDE3EF82D24AD26832469A
            R[03].s_box   AC73CF7BEFC111DF13B5D6B545235AB8
            R[03].s_row   ACC1D6B8EFB55A7B1323CFDF457311B5
            R[03].m_col   75EC0993200B633353C0CF7CBB25D0DC
            R[03].k_sch   3D80477D4716FE3E1E237E446D7A883B
            R[04].start   486C4EEE671D9D0D4DE3B138D65F58E7
            R[04].s_box   52502F2885A45ED7E311C807F6CF6A94
            R[04].s_row   52A4C89485116A28E3CF2FD7F6505E07
            R[04].m_col   0FD6DAA9603138BF6FC0106B5EB31301
            R[04].k_sch   EF44A541A8525B7FB671253BDB0BAD00
            R[05].start   E0927FE8C86363C0D9B1355085B8BE01
            R[05].s_box   E14FD29BE8FBFBBA35C89653976CAE7C
            R[05].s_row   E1FB967CE8C8AE9B356CD2BA974FFB53
            R[05].m_col   25D1A9ADBD11D168B63A338E4C4CC0B0
            R[05].k_sch   D4D1C6F87C839D87CAF2B8BC11F915BC
            R[06].start   F1006F55C1924CEF7CC88B325DB5D50C
            R[06].s_box   A163A8FC784F29DF10E83D234CD503FE
            R[06].s_row   A14F3DFE78E803FC10D5A8DF4C632923
            R[06].m_col   4B868D6D2C4A8980339DF4E837D218D8
            R[06].k_sch   6D88A37A110B3EFDDBF98641CA0093FD
            R[07].start   260E2E173D41B77DE86472A9FDD28B25
            R[07].s_box   F7AB31F02783A9FF9B4340D354B53D3F
            R[07].s_row   F783403F27433DF09BB531FF54ABA9D3
            R[07].m_col   1415B5BF461615EC274656D7342AD843
            R[07].k_sch   4E54F70E5F5FC9F384A64FB24EA6DC4F
            R[08].start   5A4142B11949DC1FA3E019657A8C040C
            R[08].s_box   BE832CC8D43B86C00AE1D44DDA64F2FE
            R[08].s_row   BE3BD4FED4E1F2C80A642CC0DA83864D
            R[08].m_col   00512FD1B1C889FF54766DCDFA1B99EA
            R[08].k_sch   EAD27321B58DBAD2312BF5607F8D292F
            R[09].start   EA835CF00445332D655D98AD8596B0C5
            R[09].s_box   87EC4A8CF26EC3D84D4C46959790E7A6
            R[09].s_row   876E46A6F24CE78C4D904AD897ECC395
            R[09].m_col   473794ED40D4E4A5A3703AA64C9F42BC
            R[09].k_sch   AC7766F319FADC2128D12941575C006E
            R[10].s_box   E9098972CB31075F3D327D94AF2E2CB5
            R[10].s_row   E9317DB5CB322C723D2E895FAF090794
            R[10].k_sch   D014F9A8C9EE2589E13F0CC8B6630CA6
            R[10].output  3925841D02DC09FBDC118597196A0B32
            sage: set_verbose(0)
        """
        r = self.r
        c = self.c
        n = self.n
        e = self.e
        F = self.base_ring()

        _type = self.state_array

        if isinstance(P,str):
            P = self.state_array([F.fetch_int(ZZ(P[i:i+2],16)) for i in range(0,len(P),2)])
        if isinstance(K,str):
            K = self.state_array([F.fetch_int(ZZ(K[i:i+2],16)) for i in range(0,len(K),2)])

        if self.is_state_array(P) and self.is_state_array(K):
            _type = self.state_array
        elif self.is_vector(P) and self.is_vector(K):
            _type = self.vector
        else:
            raise TypeError, "plaintext or key parameter not understood"

        P = self.state_array(P)
        K = self.state_array(K)

        AddRoundKey = self.add_round_key
        SubBytes = self.sub_bytes
        MixColumns = self.mix_columns
        ShiftRows = self.shift_rows
        KeyExpansion = self.key_schedule

        P = AddRoundKey(P,K)

        for r in range(self._n-1):
            if get_verbose() >= 2: print "R[%02d].start   %s"%(r+1,self.hex_str_vector(P))

            P = SubBytes(P)
            if get_verbose() >= 2: print "R[%02d].s_box   %s"%(r+1,self.hex_str_vector(P))

            P = ShiftRows(P)
            if get_verbose() >= 2: print "R[%02d].s_row   %s"%(r+1,self.hex_str_vector(P))

            P = MixColumns(P)
            if get_verbose() >= 2: print "R[%02d].m_col   %s"%(r+1,self.hex_str_vector(P))

            K = KeyExpansion(K,r+1)
            if get_verbose() >= 2: print "R[%02d].k_sch   %s"%(r+1,self.hex_str_vector(K))

            P = AddRoundKey(P, K)

        P = SubBytes(P)
        if get_verbose() >= 2: print "R[%02d].s_box   %s"%(self.n,self.hex_str_vector(P))

        P = ShiftRows(P)
        if get_verbose() >= 2: print "R[%02d].s_row   %s"%(self.n,self.hex_str_vector(P))

        if not self._star:
            P = MixColumns(P)
            if get_verbose() >= 2: print "R[%02d].m_col   %s"%(self.n,self.hex_str_vector(P))

        K = KeyExpansion(K,self._n)
        if get_verbose() >= 2: print "R[%02d].k_sch   %s"%(self.n,self.hex_str_vector(K))

        P = AddRoundKey(P, K)
        if get_verbose() >= 2: print "R[%02d].output  %s"%(self.n,self.hex_str_vector(P))

        return _type(P)

    def hex_str(self, M, type="matrix"):
        """
        Return a hex string for the provided AES state array/matrix.

        INPUT:
            M -- state array
            type -- controls what to return, either 'matrix' or 'vector'
                    (default: 'matrix')
        """
        if type == "matrix":
            return self._hex_str_matrix(M)
        elif type == "vector":
            return self._hex_str_vector(M)
        else:
            raise TypeError, "parameter type must either be 'matrix' or 'vector'"

    def hex_str_matrix(self, M):
        """
        Return a two-dimensional AES like representation of the matrix M.

        That is, show the finite field elements as hex strings.

        INPUT:
            M -- an AES state array
        """
        e = M.base_ring().degree()
        st = [""]
        for x in range(M.nrows()):
            for y in range(M.ncols()):
                if e == 8:
                    st.append("%02X"%(int(str(M[x,y].int_repr()))))
                else:
                    st.append("%X"%(int(str(M[x,y].int_repr()))))
            st.append("\n")
        return " ".join(st)

    def hex_str_vector(self, M):
        """
        Return a one dimensional AES like representation of the matrix M.

        That is, show the finite field elements as hex strings.

        INPUT:
            M -- an AES state array

        """
        e = M.base_ring().degree()
        st = [""]
        for y in range(M.ncols()):
            for x in range(M.nrows()):
                if e == 8:
                    st.append("%02X"%(int(str(M[x,y].int_repr()))))
                else:
                    st.append("%X"%(int(str(M[x,y].int_repr()))))
            #st.append("\n")
        return "".join(st)

    def _insert_matrix_into_matrix(self, dst, src, row, col):
        """
        Insert matrix src into matrix dst starting at row and col.

        INPUT:
            dst -- a matrix
            src -- a matrix
            row -- offset row
            col -- offset columns

        EXAMPLE:
            sage: sr = mq.SR(10,4,4,4)
            sage: a = sr.k.gen()
            sage: A = sr.state_array() + 1; A
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: B = Matrix(sr.base_ring(),2,2,[0,a,a+1,a^2]); B
            [    0     a]
            [a + 1   a^2]
            sage: sr._insert_matrix_into_matrix(A,B,1,1)
            [    1     0     0     0]
            [    0     0     a     0]
            [    0 a + 1   a^2     0]
            [    0     0     0     1]
        """
        for i in range(src.nrows()):
            for j in range(src.ncols()):
                dst[row+i,col+j] = src[i,j]
        return dst


    def varformatstr(self, name, n=None, rc=None, e=None):
        """
        Return a format string which is understood by print et al.

        If a numberical value is omitted the default value of self is
        used. The numerical values (n,rc,e) are used to determine the
        width of the respective fields in the format string.

        INPUT:
            name -- name of the variable
            n -- number of rounds (default: None)
            rc -- number of rows * number of cols (default: None)
            e -- exponent of base field (default: None)

        EXAMPLE:
            sage: sr = mq.SR(1,2,2,4)
            sage: sr.varformatstr('x')
            'x%01d%01d%01d'
            sage: sr.varformatstr('x', n=1000)
            'x%03d%03d%03d'
        """
        if n is None:
            n = self.n
        if rc is None:
            rc = self.r * self.c
        if e is None:
            e = self.e

        l = str(max([  len(str(rc-1)), len(str(n-1)), len(str(e-1)) ] ))
        if name != "k":
            pf = self._postfix
        else:
            pf = ""
        format_string = name + pf + "%0" + l + "d" + "%0" + l + "d" + "%0" + l + "d"
        return format_string

    def varstrs(self, name, round, rc = None, e = None):
        """
        Return a list of strings representing variables in self.

        INPUT:
            name -- variable name
            round -- number of round to create variable strings for
            rc -- number of rounds * number of columns in the state array
            e -- exponent of base field

        EXAMPLE:
            sr = mq.SR(10,1,2,4)
            sr._varstrs('x',2)
            ['x200', 'x201', 'x202', 'x203', 'x210', 'x211', 'x212', 'x213']

        """
        if rc is None:
            rc = self.r * self.c

        if e is None:
            e = self.e

        n = self._n

        format_string = self.varformatstr(name,n,rc,e)

        return [format_string%(round,rci,ei) for rci in range(rc) for ei in range(e)]

    def vars(self, name, round, rc=None, e=None):
        """
        Return a list ofvariables in self.

        INPUT:
            name -- variable name
            round -- number of round to create variable strings for
            rc -- number of rounds * number of columns in the state array
            e -- exponent of base field

        EXAMPLE:
            sr = mq.SR(10,1,2,4)
            sr._vars('x',2)
            [x200, x201, x202, x203, x210, x211, x212, x213]

        """
        gd = self.R.gens_dict()
        return [gd[e] for e in self.varstrs(name,round,rc,e)]

    def block_order(self):
        """
        Return a block order for self where each round is a block.
        """
        r = self.r
        c = self.c
        e = self.e
        n = self.n
        k = self.k

        T = None
        for _n in range(n):
            T = TermOrder('degrevlex', r*e + 3*r*c*e ) + T

        T += TermOrder('degrevlex',r*c*e);

        return T

    def ring(self, order=None):
        """
        Construct a ring as a base ring for the polynomial system.

        Variables are ordered in the reverse of their natural
        ordering, i.e. the reverse of as they appear.
        """
        r = self.r
        c = self.c
        e = self.e
        n = self.n
        if not self._gf2:
            k = self.base_ring()
        else:
            k = GF(2)

        if order is not None:
            self._order = order
        if self._order == 'block':
            self._order = self.block_order()

        names = []

        for _n in reversed(xrange(n)):
            names += self.varstrs("k",_n+1,r*c,e)
            names += self.varstrs("x",_n+1,r*c,e)
            names += self.varstrs("w",_n+1,r*c,e)
            names += self.varstrs("s",_n,r,e)

        names +=  self.varstrs("k",0,r*c,e)


        return PolynomialRing(k, 2*n*r*c*e + (n+1)*r*c*e + n*r*e, names, order=self._order)

    def round_polynomials(self, i, plaintext = None, ciphertext = None):
        r"""
        Return list of polynomials for a given round $i$.

        If $i == 0$ a plaintext must be provided, if $i == n$ a
        ciphertext must be provided.

        INPUT:
           i -- round number
           plaintext -- optional plaintext (mandatory in first round)
           ciphertext -- optional ciphertext (mandatory in last round)

        OUTPUT:
            MPolynomialRoundSystem

        EXAMPLE:
            sage: sr = mq.SR(1,1,1,4)
            sage: k = sr.base_ring()
            sage: p = [k.random_element() for _ in range(sr.r*sr.c)]
            sage: sr.round_polynomials(0,plaintext=p) # random
            [w100 + k000 + (a^2), w101 + k001 + (a + 1), w102 + k002 + (a^2 + 1), w103 + k003 + (a)]
        """
        r = self._r
        c = self._c
        e = self._e
        n = self._n
        R = self.R

        M = self.M

        _vars = self.vars

        if i == 0:
            w1 = Matrix(R, r*c*e, 1, _vars("w",1,r*c,e))
            k0 = Matrix(R, r*c*e, 1, _vars("k",0,r*c,e))
            if isinstance(plaintext,(tuple,list)) and len(plaintext) == r*c:
                plaintext = Matrix(R, r*c*e, 1, self.phi(plaintext))
            return MPolynomialRoundSystem(R, w1 + k0 + plaintext)

        elif i>0 and i<=n:
            xj = Matrix(R, r*c*e, 1, _vars("x",i,r*c,e))
            ki = Matrix(R, r*c*e, 1, _vars("k",i,r*c,e))
            rcon = Matrix(R, r*c*e, 1, self.phi([self.sbox_constant()]*r*c))

            if i < n:
                wj = Matrix(R, r*c*e, 1, _vars("w",i+1,r*c,e))
            if i == n:
                if isinstance(ciphertext,(tuple,list)) and len(ciphertext) == r*c:
                    ciphertext = Matrix(R, r*c*e, 1, self.phi(ciphertext))
                wj = ciphertext

            lin = (wj + ki + M * xj + rcon).list()


            wi = Matrix(R, r*c*e, 1, _vars("w",i,r*c,e))
            xi = Matrix(R, r*c*e, 1, _vars("x",i,r*c,e))
            sbox = []
            sbox += self.inversion_polynomials(xi,wi,r*c*e)
            sbox += self.field_polynomials("x",i)
            sbox += self.field_polynomials("w",i)
            return MPolynomialRoundSystem(R, lin + sbox)

    def key_schedule_polynomials(self, i):
        """
        Return polynomials for the $i$-th round of the key schedule.

        INPUT:
            i -- round (0 <= i <= n)
        """
        R = self.R
        r = self.r
        e = self.e
        c = self.c
        k = self.k
        a = k.gen()

        if i < 0:
            raise TypeError, "i must by >= 0"

        if i == 0:
            return MPolynomialRoundSystem(R, self.field_polynomials("k",i, r*c))
        else:
            L = self.lin_matrix(r)
            ki = Matrix(R, r*c*e, 1, self.vars("k", i  ,r*c,e))
            kj = Matrix(R, r*c*e, 1, self.vars("k", i-1,r*c,e))
            si = Matrix(R, r*e, 1, self.vars("s",i-1,r,e))

            rc = Matrix(R, r*e, 1, self.phi([a**(i-1)] + [k(0)]*(r-1)) )
            d  = Matrix(R, r*e, 1, self.phi([self.sbox_constant()]*r) )

            sbox = []

            sbox += self.field_polynomials("k",i)
            sbox += self.field_polynomials("s",i-1, r)

            if r == 1:
                sbox += self.inversion_polynomials(kj[(c - 1)*e:(c - 1)*e + e], si[0:e], e)
            if r == 2:
                sbox += self.inversion_polynomials( kj[(2*c -1)*e : (2*c -1)*e + e] , si[0:1*e], e )
                sbox += self.inversion_polynomials( kj[(2*c -2)*e : (2*c -2)*e + e] , si[e:2*e], e )
            if r == 4:
                if self._aes_mode:
                    sbox += self.inversion_polynomials( kj[(4*c-3)*e  : (4*c-3)*e + e] , si[0*e : 1*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-2)*e  : (4*c-2)*e + e] , si[1*e : 2*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-1)*e  : (4*c-1)*e + e] , si[2*e : 3*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-4)*e  : (4*c-4)*e + e] , si[3*e : 4*e] , e )
                else:
                    sbox += self.inversion_polynomials( kj[(4*c-1)*e  : (4*c-1)*e + e] , si[0*e : 1*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-2)*e  : (4*c-2)*e + e] , si[1*e : 2*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-3)*e  : (4*c-3)*e + e] , si[2*e : 3*e] , e )
                    sbox += self.inversion_polynomials( kj[(4*c-4)*e  : (4*c-4)*e + e] , si[3*e : 4*e] , e )

            si =  L * si + d + rc
            Sum = Matrix(R, r*e, 1)
            lin = []
            if c>1:
                for q in range(c):
                    t = range(r*e*(q) , r*e*(q+1) )
                    Sum += kj.matrix_from_rows(t)
                    lin += (ki.matrix_from_rows(t) + si + Sum).list()

            else:
                lin += (ki + si).list()
            return MPolynomialRoundSystem(R, lin + sbox )

    def polynomial_system(self, P=None, K=None):
        """
        Return a MPolynomialSystem for self for a given plaintext-key pair.

        If none are provided a random pair will be generated.

        INPUT:
            P -- vector, list, or tuple (default: None)
            K -- vector, list, or tuple (default: None)
        """
        plaintext = P
        key = K

        system = []
        n = self._n

        data = []

        for d in (plaintext, key):
            if d is None:
                data.append(self.random_element("vector"))
            elif isinstance(d, (tuple,list)):
                data.append( self.phi(self.state_array(d)) )
            elif self.is_state_array(d):
                data.append( self.phi(d) )
            elif self.is_vector(d):
                data.append( d )
            else:
                raise TypeError, "type %s of %s not understood"%(type(d),d)

        plaintext, key = data

        ciphertext = self(plaintext, key)

        for i in range(n+1):
            system.append( self.round_polynomials(i,plaintext, ciphertext) )
            system.append( self.key_schedule_polynomials(i) )

        return MPolynomialSystem(self.R, system), dict(zip(self.vars("k",0),key.list()))


class SR_gf2n(SR_generic):
    r"""
    Small Scale Variants of the AES polynomial system constructor over $GF(2^n)$.
    """
    def vector(self, d=None):
        """
        Constructs a vector suitable for the algebraic representation of SR, i.e. BES.

        INPUT:
            d -- values for vector, must be understood by self.phi (default:None)
        """
        r = self.r
        c = self.c
        e = self.e
        k = self.base_ring()

        if d is None:
            return Matrix(k, r*c*e, 1)
        elif d.ncols() == c and d.nrows() == r and d.base_ring() == k:
            return Matrix(k, r*c*e, 1, self.phi(d).transpose())

    def is_vector(self, d):
        """
        Return True if d can be used as a vector for self.
        """
        return is_Matrix(d) and \
               d.nrows() == self.r*self.c*self.e and \
               d.ncols() == 1 and \
               d.base_ring() == self.base_ring()

    def phi(self, l):
        r"""
        Projects state arrays to their algebraic representation.
        """
        ret = []
        if is_Matrix(l):
            for e in l.transpose().list():
                ret += [e**(2**i) for i in range(self.e)]
        else:
            for e in l:
                ret += [e**(2**i) for i in range(self.e)]
        if isinstance(l,list):
            return ret
        elif isinstance(l,tuple):
            return tuple(ret)
        elif is_Matrix(l):
            return Matrix(l.base_ring(),l.ncols(), l.nrows()*self.e, ret).transpose()
        else:
            raise TypeError

    def antiphi(self,l):
        """
        Inverse of self.phi.
        """
        if is_Matrix(l):
            ret = [e for e in l.transpose().list()[0:-1:self.e]]
        else:
            ret = [e for e in l[0:-1:self.e]]

        if isinstance(l,list):
            return ret
        elif isinstance(l,tuple):
            return tuple(ret)
        elif is_Matrix(l):
            return Matrix(self.base_ring(),l.ncols(), l.nrows()/self.e, ret).transpose()
        else:
            raise TypeError

    def shift_rows_matrix(self):
        """
        Return the ShiftRows matrix.

        EXAMPLE:
            sage: sr = mq.SR(1,2,2,4)
            sage: s = sr.random_state_array()
            sage: r1 = sr.shift_rows(s)
            sage: r2 = sr.state_array( sr.ShiftRows * sr.vector(s) )
            sage: r1 == r2
            True
        """
        e = self.e
        r = self.r
        c = self.c
        k = self.base_ring()
        bs = r*c*e
        shift_rows = Matrix(k,bs,bs)
        I = MatrixSpace(k,e,e)(1)
        for x in range(0, c):
            for y in range(0,r):
                _r = ((x*r)+y) * e
                _c = (((x*r)+((r+1)*y)) * e) % bs
                self._insert_matrix_into_matrix(shift_rows,I, _r, _c)

        return shift_rows

    def lin_matrix(self, length = None):
        """
        Return the Lin matrix.

        If no length is provided the standard state space size is
        used. The key schedule calls this method with an explicit
        length argument because only self.r S-Box applications are
        performed in the key schedule.

        INPUT:
            length -- length of state space. (default: None)
        """
        r = self.r
        c = self.c
        e = self.e
        k = self.k

        if length is None:
            length = r*c

        lin = Matrix(self.base_ring(), length*e, length*e)
        if e == 4:
            l = [ k.fetch_int(x) for x in  (5, 1, 12, 5) ]
            for k in range( 0, length ):
                for i in range(0,4):
                    for j in range(0,4):
                        lin[k*4+j,k*4+i]=  l[(i-j)%4] ** (2**j)
        elif e == 8:
            l = [ k.fetch_int(x) for x in  (5, 9, 249, 37, 244, 1, 181, 143) ]
            for k in range( 0, length ):
                for i in range(0,8):
                    for j in range(0,8):
                        lin[k*8+j,k*8+i]=  l[(i-j)%8] ** (2**j)

        return lin

    def mix_columns_matrix(self):
        """
        Return the MixColumns matrix.

        EXAMPLE:
            sage: sr = mq.SR(1,2,2,4)
            sage: s = sr.random_state_array()
            sage: r1 = sr.mix_columns(s)
            sage: r2 = sr.state_array(sr.MixColumns * sr.vector(s))
            sage: r1 == r2
            True
        """

        def D(b):
            D = Matrix(self.base_ring(), self._e, self._e)
            for i in range(self._e):
                D[i,i] = b**(2**i)
            return D

        r = self.r
        c = self.c
        e = self.e
        k = self.k
        a = k.gen()

        M = Matrix(k,r*e,r*e)

        if r == 1:
            self._insert_matrix_into_matrix(M,   D(1), 0, 0)

        elif r == 2:
            self._insert_matrix_into_matrix(M, D(a+1), 0, 0)
            self._insert_matrix_into_matrix(M, D(a+1), e, e)
            self._insert_matrix_into_matrix(M,   D(a), e, 0)
            self._insert_matrix_into_matrix(M,   D(a), 0, e)

        elif r == 4:
            self._insert_matrix_into_matrix(M,   D(a),   0,   0)
            self._insert_matrix_into_matrix(M,   D(a),   e,   e)
            self._insert_matrix_into_matrix(M,   D(a), 2*e, 2*e)
            self._insert_matrix_into_matrix(M,   D(a), 3*e, 3*e)

            self._insert_matrix_into_matrix(M, D(a+1),   0,   e)
            self._insert_matrix_into_matrix(M, D(a+1),   e, 2*e)
            self._insert_matrix_into_matrix(M, D(a+1), 2*e, 3*e)
            self._insert_matrix_into_matrix(M, D(a+1), 3*e,   0)

            self._insert_matrix_into_matrix(M,   D(1),   0, 2*e)
            self._insert_matrix_into_matrix(M,   D(1),   e, 3*e)
            self._insert_matrix_into_matrix(M,   D(1), 2*e,   0)
            self._insert_matrix_into_matrix(M,   D(1), 3*e, 1*e)

            self._insert_matrix_into_matrix(M,   D(1),   0, 3*e)
            self._insert_matrix_into_matrix(M,   D(1),   e,   0)
            self._insert_matrix_into_matrix(M,   D(1), 2*e, 1*e)
            self._insert_matrix_into_matrix(M,   D(1), 3*e, 2*e)

        mix_columns = Matrix(k,r*c*e,r*c*e)

        for i in range(c):
            self._insert_matrix_into_matrix(mix_columns, M, r*e*i, r*e*i)

        return mix_columns

    def inversion_polynomials(self, xi, wi, length):
        """
        Return polynomials to represent the inversion in the AES S-Box.

        INPUT:
            xi -- output variables
            wi -- input variables
            length -- length of both lists
        """
        return [xi[j,0]*wi[j,0] + 1 for j in range(length)]

    def field_polynomials(self, name, i, l=None):
        """
        Return list of conjugacy polynomials for a given round and name.

        INPUT:
            name -- variable name
            i -- round number
            l -- r*c (default:None)

        EXAMPLE:
            sage: sr = mq.SR(3,1,1,8)
            sage: sr.field_polynomials('x',2)
            [x200^2 + x201,
            x201^2 + x202,
            x202^2 + x203,
            x203^2 + x204,
            x204^2 + x205,
            x205^2 + x206,
            x206^2 + x207,
            x207^2 + x200]
        """
        r = self._r
        c = self._c
        e = self._e
        n = self._n

        if l is None:
            l = r*c

        fms = self.varformatstr(name,n,l,e)

        return [self.R( fms%(i,rci,ei) + "**2 + " + fms%(i, rci,(ei+1)%e) )  for rci in range(l)  for ei in range(e) ]



class SR_gf2(SR_generic):
    r"""
    Small Scale Variants of the AES polynomial system constructor over $GF(2)$.
    """
    def __init__(self,n=1,r=1,c=1,e=4, star=False, **kwargs):
        """
        See help for SR.

        """
        SR_generic.__init__(self,n,r,c,e,star,**kwargs)
        self._correct_only = kwargs.get("correct_only",False)
        self._biaffine_only = kwargs.get("biaffine_only",True)

    def vector(self, d=None):
        """
        Constructs a vector suitable for the algebraic representation of SR.

        INPUT:
            d -- values for vector(default:None)
        """
        r = self.r
        c = self.c
        e = self.e
        k = GF(2)

        if d is None:
            return Matrix(k, r*c*e, 1)
        elif is_Matrix(d) and d.ncols() == c and d.nrows() == r and d.base_ring() == self.k:
            l = flatten([self.phi(x) for x in d.transpose().list()], Vector_modn_dense)
            return Matrix(k, r*c*e, 1,l)
        elif isinstance(d,(list,tuple)):
            if len(d) == self.r*self.c:
                l = flatten([self.phi(x) for x in d], Vector_modn_dense)
                return Matrix(k, r*c*e, 1,l)
            elif len(d) == self.r*self.c*self.e:
                return Matrix(k, r*c*e, 1, d)
            else:
                raise TypeError
        else:
            raise TypeError

    def is_vector(self, d):
        """
        Return True if the given matrix satisfies the conditions for a vector
        as it appears in the algebraic expression of self.

        INPUT:
            d -- matrix
        """
        return is_Matrix(d) and \
               d.nrows() == self.r*self.c*self.e and \
               d.ncols() == 1 and \
               d.base_ring() == GF(2)

    def phi(self, l, diffusion_matrix=False):
        r"""
        Given a list/matrix of elements in $GF(2^n)$ return a matching
        list/matrix of elements in $GF(2)$.

        INPUT:
            l -- element to perform phi on.
            diffusion_matrix -- if True the given matrix $l$ is transformed
                                to a matrix which performs the same operation
                                over GF(2) as $l$ over $GF(2^n)$ (default: False).
        """
        ret = []
        r,c,e = self.r,self.c,self.e

        # handle diffusion layer matrices first
        if is_Matrix(l) and diffusion_matrix and \
           l.nrows() == r*c and l.ncols() == r*c and \
           l.base_ring() == self.k:
            B = Matrix(GF(2), r*c*e, r*c*e)
            for x in range(r*c):
                for y in range(r*c):
                    T = self._mul_matrix(l[x,y])
                    self._insert_matrix_into_matrix(B,T,x*e, y*e)
            return B

        # ground field elements
        if l in self.k:
            return list(reversed(l.vector()))

        # remaining matrices
        if is_Matrix(l):
            for x in l.transpose().list():
                ret += list(reversed(x.vector()))
        # or lists
        else:
            for x in l:
                ret += list(reversed(x.vector()))

        if isinstance(l,list): return ret
        elif isinstance(l,tuple): return tuple(ret)
        elif is_Matrix(l): return Matrix(GF(2), l.ncols(), l.nrows()*self.e, ret).transpose()
        else: raise TypeError

    def antiphi(self,l):
        """
        Inverse of self.phi.
        """
        e = self.e
        V = self.k.vector_space()

        if is_Matrix(l):
            l2 = l.transpose().list()
        else:
            l2 = l

        ret = []
        for i in range(0,len(l2),e):
            ret.append( self.k(V(list(reversed(l2[i:i+e])))) )

        if isinstance(l, list):
            return ret
        elif isinstance(l,tuple):
            return tuple(ret)
        elif is_Matrix(l):
            return Matrix(self.base_ring(), self.r *self.c, 1, ret)
        else:
            raise TypeError

    def shift_rows_matrix(self):
        """
        Return the ShiftRows matrix.

        EXAMPLE:
            sage: sr = mq.SR(1,2,2,4,gf2=True)
            sage: s = sr.random_state_array()
            sage: r1 = sr.shift_rows(s)
            sage: r2 = sr.state_array( sr.ShiftRows * sr.vector(s) )
            sage: r1 == r2
            True
        """
        r = self.r
        c = self.c
        k = self.k
        bs = r*c
        shift_rows = Matrix(k,r*c,r*c)
        for x in range(0, c):
            for y in range(0,r):
                _r = ((x*r)+y)
                _c = ((x*r)+((r+1)*y)) % bs
                shift_rows[_r,_c] = 1
        return self.phi(shift_rows, diffusion_matrix=True)

    def mix_columns_matrix(self):
        """
        Return the MixColumns matrix.

        EXAMPLE:
            sage: sr = mq.SR(1,2,2,4,gf2=True)
            sage: s = sr.random_state_array()
            sage: r1 = sr.mix_columns(s)
            sage: r2 = sr.state_array(sr.MixColumns * sr.vector(s))
            sage: r1 == r2
            True
        """
        r = self.r
        c = self.c
        k = self.k
        a = k.gen()



        if r == 1:
            M = Matrix(k,r,r,1)

        elif r == 2:
            M = Matrix(k,r,r,[a+1,a,a,a+1])

        elif r == 4:
            M = Matrix(k,r, [a,a+1,1,1,\
                             1,a,a+1,1,\
                             1,1,a,a+1,\
                             a+1,1,1,a])

        mix_columns = Matrix(k,r*c,r*c)

        for i in range(c):
            self._insert_matrix_into_matrix(mix_columns, M, r*i, r*i)

        return self.phi(mix_columns, diffusion_matrix=True)

    def lin_matrix(self, length=None):
        """
        Return the Lin matrix.

        If no length is provided the standard state space size is
        used. The key schedule calls this method with an explicit
        length argument because only self.r S-Box applications are
        performed in the key schedule.

        INPUT:
            length -- length of state space. (default: None)
        """
        r,c,e = self.r, self.c,self.e

        if length is None:
            length = r*c

        if e == 8:
            Z = Matrix(GF(2),8,8,[1,0,0,0,1,1,1,1,\
                                  1,1,0,0,0,1,1,1,\
                                  1,1,1,0,0,0,1,1,\
                                  1,1,1,1,0,0,0,1,\
                                  1,1,1,1,1,0,0,0,\
                                  0,1,1,1,1,1,0,0,\
                                  0,0,1,1,1,1,1,0,\
                                  0,0,0,1,1,1,1,1])
        else:
            Z = Matrix(GF(2),4,4,[1,1,1,0,\
                                  0,1,1,1,\
                                  1,0,1,1,\
                                  1,1,0,1])


        Z = Z.transpose() # account for endianess mismatch

        lin = Matrix(GF(2),length*e,length*e)

        for i in range(length):
            self._insert_matrix_into_matrix(lin,Z,i*e,i*e)
        return lin

    def _mul_matrix(self, x):
        """
        Given an element $x$ in self.base_ring() return a matrix which
        performs the same operation on a when interpreted over
        $GF(2)^e$ as $x$ over $GF(2^e)$.

        INPUT:
            x -- an element in self.base_ring()

        EXAMPLE:
            sage: sr = mq.SR(gf2=True)
            sage: a = sr.k.gen()
            sage: A = sr._mul_matrix(a^2+1)
            sage: sr.antiphi( A *  sr.vector([a+1]) )
            [a^3 + a^2 + a + 1]

            sage: (a^2 + 1)*(a+1)
            a^3 + a^2 + a + 1
        """
        a = self.k.gen()
        k = self.k
        e = self.e
        a = k.gen()

        columns = []
        for i in reversed(range(e)):
            columns.append( list(reversed((x * a**i).vector())) )
        return Matrix(GF(2), e, e, columns).transpose()

    def _square_matrix(self):
        """
        Return a matrix of dimension self.e x self.e which performs
        the squaring operation over GF(2^n) on vectors of length e.

        EXAMPLE:
            sage: sr = mq.SR(gf2=True)
            sage: a = sr.k.gen()
            sage: S = sr._square_matrix()
            sage: sr.antiphi( S *  sr.vector([a^3+1]) )
            [a^3 + a^2 + 1]

            sage: (a^3 + 1)^2
            a^3 + a^2 + 1

        """
        a = self.k.gen()
        e = self.e

        columns = []
        for i in reversed(range(e)):
            columns.append( list(reversed(((a**i)**2).vector())) )
        return Matrix(GF(2), e ,e, columns).transpose()

    def inversion_polynomials_single_sbox(self, x= None, w=None, biaffine_only=None, correct_only=None):
        """
        Generator for S-Box inversion polynomials of a single sbox.

        INPUT:
            x -- output variables (default: None)
            w -- input variables  (default: None)
            biaffine_only -- only include biaffine polynomials (default: object default)
            correct_only -- only include correct polynomials (default: object default)

        EXAMPLES:
            sage: sr = mq.SR(1,1,1,8,gf2=True)
            sage: len(sr.inversion_polynomials_single_sbox())
            24
            sage: len(sr.inversion_polynomials_single_sbox(correct_only=True))
            23
            sage: len(sr.inversion_polynomials_single_sbox(biaffine_only=False))
            40
            sage: len(sr.inversion_polynomials_single_sbox(biaffine_only=False, correct_only=True))
            39
        """
        e = self.e

        if biaffine_only is None:
            biaffine_only = self._biaffine_only
        if correct_only is None:
            correct_only = self._correct_only

        if x is None and w is None:
            # make sure it prints like in the book.
            names = ["w%d"%i for i in reversed(range(e))]+ ["x%d"%i for i in reversed(range(e))]
            P = PolynomialRing(GF(2),e*2, names, order='lex')
            x = Matrix(P,e,1,P.gens()[e:])
            w = Matrix(P,e,1,P.gens()[:e])
        else:
            if isinstance(x,(tuple,list)): P = x[0].parent()
            elif is_Matrix(x): P = x.base_ring()
            else: raise TypeError, "x not understood"

            if isinstance(x, (tuple,list)):
                x = Matrix(P,e,1, x)
            if isinstance(w, (tuple,list)):
                w = Matrix(P,e,1, w)

        T = self._mul_matrix(self.k.gen())
        o = Matrix(P,e,1,[0]*(e-1) + [1])

        columns = []
        for i in reversed(range(e)):
            columns.append((T**i * w).list())
        Cw = Matrix(P,e,e, columns).transpose()

        columns = []
        for i in reversed(range(e)):
            columns.append((T**i * x).list())
        Cx = Matrix(P,e,e, columns).transpose()

        S = self._square_matrix()

        l = []
        if correct_only:
            l.append( (Cw * x + o).list()[:-1] )
        else:
            l.append( (Cw * x + o).list() )
        l.append( (Cw * S *x  + x).list() )
        l.append( (Cx * S *w  + w).list() )
        if not biaffine_only:
            l.append( ((Cw * S**2 + Cx*S)*x).list() )
            l.append( ((Cx * S**2 + Cw*S)*w).list() )

        return sum(l,[])

    def inversion_polynomials(self, xi, wi, length):
        """
        Return polynomials to represent the inversion in the AES S-Box.

        INPUT:
            xi -- output variables
            wi -- input variables
            length -- length of both lists
        """
        if is_Matrix(xi):
            xi = xi.list()
        if is_Matrix(wi):
            wi = wi.list()

        e = self.e
        l = []
        for j in range(0,length, e):
            l += self.inversion_polynomials_single_sbox(xi[j:j+e], wi[j:j+e])
        return l


    def field_polynomials(self, name, i, l=None):
        """
        Return list of field polynomials for a given round -- given by its number $i$ -- and name.

        INPUT:
            name -- variable name
            i -- round number
            l -- length of variable list (default:None => r*c)

        EXAMPLE:
            sage: sr = mq.SR(3,1,1,8)
            sage: sr.field_polynomials('x',2)
            [x200^2 + x201,
            x201^2 + x202,
            x202^2 + x203,
            x203^2 + x204,
            x204^2 + x205,
            x205^2 + x206,
            x206^2 + x207,
            x207^2 + x200]
        """
        r = self._r
        c = self._c
        e = self._e
        n = self._n

        if l is None:
            l = r*c

        fms = self.varformatstr(name,n,l,e)
        return [self.R( fms%(i,rci,ei) + "**2 + " + fms%(i, rci,ei) )  for rci in range(l)  for ei in range(e) ]


def test_consistency(max_n=2, **kwargs):
    r"""
    Test all combinations of r,c,e and n in (1,2) for consistency of
    random encryptions and their polynomial systems. $GF(2)$ and $GF(2^e)$
    systems are tested. This test takes a while.

    INPUT:
        max_n -- maximal number of rounds to consider.
        kwargs -- are passed to the SR constructor

    TESTS:
        sage: from sage.crypto.mq.sr import test_consistency
        sage: test_consistency(2) # long time
        True
    """
    consistent = True
    for r in (1,2,4):
      for c in (1,2,4):
        for e in (4,8):
          for n in range(1,max_n+1):
            for gf2 in (True,False):
              zero_division = True
              while zero_division:
                sr = SR(n,r,c,e,gf2=gf2, **kwargs)
                try:
                  F, s = sr.polynomial_system()
                  F.subs(s)
                  consistent &= (F.groebner_basis('libsingular:slimgb')[0] != 1)
                  if not consistent:
                      print sr, " is not consistent"
                  zero_division = False

                except ZeroDivisionError:
                    pass
    return consistent
