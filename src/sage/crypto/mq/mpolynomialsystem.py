"""
Multivariate Polynomial Systems.

We call a finite set of multivariate polynomials an MPolynomialSystem.

Furthermore we assume that these multivariate polynomials have a
common solution if interpreted as equations where the left hand side
is the polynomial and the right hand side is equal to zero. Or in
other terms: The set of multivariate polynomials have common roots. In
many other computer algebra systems this class could be called Ideal
but -- strictly speaking -- an ideal is a very distinct object form its
generators and thus this is not an Ideal in \SAGE.

The main purpose of this class is to manipulate an MPolynomialSystem
to gather the common solution.

This idea is specialized to an MPolynomialSystem which consists of
several rounds. These kind of polynomial systems are often found in
symmetric algebraic cryptanalysis. The most prominent examples of these
kind of systems are: SR (AES), Flurry/Curry, and CTC(2).

AUTHOR: Martin Albrecht <malb@informatik.uni-bremen.de>

TEST:
    sage: P.<x,y> = PolynomialRing(QQ)
    sage: I = [[x^2 + y^2], [x^2 - y^2]]
    sage: F = mq.MPolynomialSystem(P,I)
    sage: loads(dumps(F)) == F
    True
"""

from sage.structure.sage_object import SageObject

from sage.rings.integer_ring import ZZ
from sage.rings.finite_field import FiniteField as GF

from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.polynomial.multi_polynomial import is_MPolynomial
from sage.rings.polynomial.polynomial_ring import PolynomialRing

from sage.matrix.matrix import is_Matrix
from sage.matrix.constructor import Matrix

from sage.interfaces.singular import singular
from sage.interfaces.magma import magma


def is_MPolynomialSystem(F):
    """
    Return True if F is an MPolynomialSystem

    EXAMPLE:
       sage: P.<x,y> = PolynomialRing(QQ)
       sage: I = [[x^2 + y^2], [x^2 - y^2]]
       sage: F = mq.MPolynomialSystem(P,I); F
       Polynomial System with 2 Polynomials in 2 Variables
       sage: mq.is_MPolynomialSystem(F)
       True
    """
    return isinstance(F,MPolynomialSystem_generic)

def is_MPolynomialRoundSystem(F):
    """
    Return True if F is an MPolynomialRoundSystem

    EXAMPLE:
       sage: P.<x,y> = PolynomialRing(QQ)
       sage: I = [[x^2 + y^2], [x^2 - y^2]]
       sage: F = mq.MPolynomialSystem(P,I); F
       Polynomial System with 2 Polynomials in 2 Variables
       sage: mq.is_MPolynomialRoundSystem(F.round(0))
       True
    """
    return isinstance(F,MPolynomialRoundSystem_generic)


def MPolynomialRoundSystem(R, gens):
    """
    Construct an object representing the equations of a single
    round (e.g. of a block cipher).

    INPUT:
        R -- base ring
        gens -- list (default: [])

    EXAMPLE:
        sage: P.<x,y,z> = MPolynomialRing(GF(2),3)
        sage: mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
        [x*y + 1, z + 1]
    """
    return MPolynomialRoundSystem_generic(R, gens)

def MPolynomialSystem(arg1, arg2=None):
    """
    Construct a new MPolynomialSystem.

    INPUT:
        arg1 -- a multivariate polynomial ring or an ideal
        arg2 -- an iterable object of rounds, preferable MPolynomialRoundSystem,
                or polynomials (default:None)

    EXAMPLES:
        sage: P.<a,b,c,d> = PolynomialRing(GF(127),4)
        sage: I = sage.rings.ideal.Katsura(P)

        If a list of MPolynomialRoundSystems is provided those form
        the rounds.

        sage: mq.MPolynomialSystem(I.ring(), [mq.MPolynomialRoundSystem(I.ring(),I.gens())])
        Polynomial System with 4 Polynomials in 4 Variables

        If an ideal is provided the generators are used.

        sage: mq.MPolynomialSystem(I)
        Polynomial System with 4 Polynomials in 4 Variables

        If a list of polynomials is provided the system has only one
        round.

        sage: mq.MPolynomialSystem(I.ring(), I.gens())
        Polynomial System with 4 Polynomials in 4 Variables
    """
    if is_MPolynomialRing(arg1):
        R = arg1
        rounds = arg2
    elif isinstance(arg1,MPolynomialIdeal):
        R = arg1.ring()
        rounds = arg1.gens()
    else:
        raise TypeError, "first parameter must be a MPolynomialRing"

    k = R.base_ring()

    if k.characteristic() != 2:
        return MPolynomialSystem_generic(R,rounds)
    elif k.degree() == 1:
        return MPolynomialSystem_gf2(R,rounds)
    elif k.degree() > 1:
        return MPolynomialSystem_gf2e(R,rounds)

class MPolynomialRoundSystem_generic(SageObject):
    """
    Represents a multivariate polynomial set e.g. of a single round of
    a block cipher.
    """
    def __init__(self, R, gens):
        """
        Construct an object representing the equations of a single
        round (e.g. of a block cipher).

        INPUT:
            R -- base ring
            gens -- list (default: [])

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(GF(2),3)
            sage: mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
            [x*y + 1, z + 1]
        """
        if is_MPolynomialRing(R):
            self._ring = R
        else:
            raise TypeError, "first parameter must be a MPolynomialRing"

        if is_Matrix(gens):
            self._gens = gens.list()
        else:
            self._gens = list(gens)

    def __cmp__(self, other):
        """
        Compare the ring and generators of self and other.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F == copy(F) # indirect doctest
            True
        """
        return cmp((self._ring, self._gens),(other._ring, other._gens))

    def ring(self):
        """
        Return the base ring.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R0 = F.round(0)
            sage: print R0.ring().repr_long()
            Polynomial Ring
             Base Ring : Finite Field of size 2
                  Size : 20 Variables
              Block  0 : Ordering : degrevlex
                         Names    : k100, k101, k102, k103, x100, x101, x102, x103, w100, w101, w102, w103, s000, s001, s002, s003
              Block  1 : Ordering : degrevlex
                         Names    : k000, k001, k002, k003
        """
        return self._ring

    def ngens(self):
        """
        Return number of polynomials in self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R0 = F.round(0)
            sage: R0.ngens()
            4
        """
        return len(self._gens)

    def gens(self):
        """
        Return list of polynomials in self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R1 = F.round(1)
            sage: l = R1.gens()
            sage: l[0]
            k000^2 + k000
        """
        return list(self)

    def variables(self):
        """
        Return unordered list of variables apprearing in polynomials
        in self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R1 = F.round(1)
            sage: sorted(R1.variables())
            [k003, k002, k001, k000]

        """
        V = set()
        for f in self._gens:
            for v in f.variables():
                V.add(v)
        return list(V)

    def monomials(self):
        """
        Return unordered list of monomials appearing in polynomials
        in self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R1 = F.round(1)
            sage: sorted(R1.monomials())
            [k003, k002, k001, k000, k003^2, k002^2, k001^2, k000^2]
        """
        M = set()
        for f in self._gens:
            for m in f.monomials():
                M.add(m)
        return list(M)

    def subs(self, *args, **kwargs):
        """
        Substitute variables for every polynomial in self. See
        MPolynomial.subs for calling convention.

        INPUT:
            args -- arguments to be passed to MPolynomial.subs
            kwargs -- keyword arguments to be passed to MPolynomial.subs

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R1 = F.round(1)
            sage: R1.subs(s) # the solution
            sage: R1
            [0, 0, 0, 0]
        """
        for i in range(len(self._gens)):
            self._gens[i] = self._gens[i].subs(*args,**kwargs)

    def _repr_(self):
        """
        Return string representation of self.

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(GF(2),3)
            sage: F = mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
            sage: str(F) # indirect doctest
            '[x*y + 1, z + 1]'
        """
        return "%s"%self._gens

    def __getitem__(self, i):
        """
        Return the i-th generator of self.

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(GF(2),3)
            sage: F = mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
            sage: F[0] # indirect doctest
            x*y + 1
        """
        return self._gens[i]

    def __add__(self, right):
        """
        Addition is the union of generators.

        INPUT:
            right -- MPolynomialSystem, list or tuple

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R1 = F.round(1)
            sage: len(R1 + F.round(2)) # indirect doctest
            28
            sage: len(R1 + list(F.round(2)))
            28

        """
        if is_MPolynomialRoundSystem(right) and self.ring() == right.ring():
            return MPolynomialRoundSystem(self.ring(), self._gens + right.gens())
        if isinstance(right, (list,tuple)) and self.ring() == right[0].parent():
            return MPolynomialRoundSystem(self.ring(), self._gens + right)
        else:
            raise ArithmeticError, "Cannot add MPolynomialRoundSystem and %s"%type(right)

    def __contains__(self, element):
        """
        Return True if element is in the list of generators for self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R1 = F.round(1)
            sage: P = F.ring()
            sage: f = P('k000^2 + k000')
            sage: f in R1
            True
            sage: f+1 in R1
            False
        """
        return (element in self._gens)

    def __len__(self):
        """
        Return self.ngens().

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(GF(2),3)
            sage: F = mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
            sage: len(F)
            2

        """
        return len(self._gens)

    def __iter__(self):
        """
        Iterate over the generators of self.

        EXAMPLE:
            sage: P.<x,y,z> = MPolynomialRing(GF(2),3)
            sage: F = mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
            sage: for f in F:
            ...     print f
            x*y + 1
            z + 1
        """
        return iter(self._gens)

    def _singular_(self):
        """
        Return SINGULAR ideal representation of self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R1 = F.round(1)
            sage: R1._singular_()
            k000^2+k000,
            k001^2+k001,
            k002^2+k002,
            k003^2+k003
        """
        return singular.ideal(self._gens)

    def _magma_(self):
        """
        Return MAGMA ideal representation of self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True)
            sage: F,s = sr.polynomial_system()
            sage: R1 = F.round(1)
            sage: R1._magma_() # optional, requires MAGMA
            Ideal of Polynomial ring of rank 20 over GF(2)
            Graded Reverse Lexicographical Order
            Variables: k100, k101, k102, k103, x100, x101, x102, x103, w100, w101, w102, w103, s000, s001, s002, s003, k000, k001, k002, k003
            Basis:
            [
            k000^2 + k000,
            k001^2 + k001,
            k002^2 + k002,
            k003^2 + k003
            ]
        """
        return magma.ideal(self._gens)

class MPolynomialSystem_generic(SageObject):
    """
    A system of multivariate polynomials. That is, a set of
    multivariate polynomials with at least one common root.
    """
    def __init__(self, R, rounds):
        """
        Construct a new MPolynomialSystem.

        INPUT:
            arg1 -- a multivariate polynomial ring or an ideal
            arg2 -- an iterable object of rounds, preferable MPolynomialRoundSystem,
                      or polynomials (default:None)

        EXAMPLES:
            sage: P.<a,b,c,d> = PolynomialRing(GF(127),4)
            sage: I = sage.rings.ideal.Katsura(P)

            If a list of MPolynomialRoundSystems is provided those
            form the rounds.

            sage: mq.MPolynomialSystem(I.ring(), [mq.MPolynomialRoundSystem(I.ring(),I.gens())])
            Polynomial System with 4 Polynomials in 4 Variables

            If an ideal is provided the generators are used.

            sage: mq.MPolynomialSystem(I)
            Polynomial System with 4 Polynomials in 4 Variables

            If a list of polynomials is provided the system has only
            one round.

            sage: mq.MPolynomialSystem(I.ring(), I.gens())
            Polynomial System with 4 Polynomials in 4 Variables
        """

        self._ring = R
        self._rounds = []

        # check for list of polynomials
        e = iter(rounds).next()
        if is_MPolynomial(e):
            rounds = [rounds]

        for b in rounds:
            if not is_MPolynomialRoundSystem(b):
                if isinstance(b, (tuple,list)):
                    self._rounds.append(MPolynomialRoundSystem(R, b))
                else:
                    self._rounds.append(MPolynomialRoundSystem(R, list(b)))
            elif b.ring() is R or b.ring() == R:
                self._rounds.append(b)
            else:
                raise TypeError, "parameter not supported"

    def __cmp__(self, other):
        """
        Compare the ring and rounds of self and other.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F == copy(F) # indirect doctest
            True
        """
        return cmp((self._ring, self._rounds),(other._ring, other._rounds))

    def ring(self):
        """
        Return base ring.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: print F.ring().repr_long()
            Polynomial Ring
             Base Ring : Finite Field of size 2
                  Size : 20 Variables
              Block  0 : Ordering : degrevlex
                         Names    : k100, k101, k102, k103, x100, x101, x102, x103, w100, w101, w102, w103, s000, s001, s002, s003
              Block  1 : Ordering : degrevlex
                         Names    : k000, k001, k002, k003

        """
        return self._ring

    def ngens(self):
        """
        Return number polynomials in self

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: F.ngens()
            56
        """
        return sum([e.ngens() for e in self._rounds])

    def gens(self):
        """
        Return list of polynomials in self

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: l = F.gens()
            sage: len(l), type(l)
            (40, <type 'list'>)
        """
        return list(self)

    def gen(self, ij):
        """
        Return an element of self.

        INPUT:
            ij -- tuple, slice, integer

        EXAMPLES:
            sage: P.<a,b,c,d> = PolynomialRing(GF(127),4)
            sage: F = mq.MPolynomialSystem(sage.rings.ideal.Katsura(P))

            $ij$-th polynomial overall

            sage: F[0] # indirect doctest
            a + 2*b + 2*c + 2*d - 1

            $i$-th to $j$-th polynomial overall

            sage: F[0:2]
            [a + 2*b + 2*c + 2*d - 1, a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a]

            $i$-th round, $j$-th polynomial

            sage: F[0,1]
            a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a
        """
        return self[ij]

    def nrounds(self):
        """
        Return number of rounds of self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.nrounds()
            4
        """
        return len(self._rounds)

    def rounds(self):
        """
        Return list of rounds of self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: l = F.rounds()
            sage: len(l)
            4
        """
        return list(self._rounds)

    def round(self, i):
        """
        Return $i$-th round of self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: R0 = F.round(1)
            sage: R0
            [k000^2 + k001, k001^2 + k002, k002^2 + k003, k003^2 + k000]
        """
        return self._rounds[i]

    def __iter__(self):
        """
        Iterate over the generators of self round by round.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: type(iter(F))
            <type 'generator'>
        """
        for b in self._rounds:
            for e in b:
                yield e

    def ideal(self):
        """
        Return SAGE ideal spanned by self.gens()

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: P = F.ring()
            sage: I = F.ideal()
            sage: I.elimination_ideal(P('s000*s001*s002*s003*w100*w101*w102*w103*x100*x101*x102*x103')) # random result
            Ideal (k002 + (a)*k003 + (a^3 + 1), k001 + (a^2 + 1)*k002
            + (a^3 + a + 1), k000 + (a^3 + a^2 + 1)*k003 + (a^3 + a^2
            + a), k103 + (a^2 + 1)*k000 + (a + 1)*k001 + (a^3 +
            a^2)*k002 + (a^2 + 1)*k003 + 1, k102 + (a^3 + a)*k103 +
            (a^2 + 1)*k001 + (a)*k002 + (a^2 + a + 1)*k003 + (a^3 + a
            + 1), k101 + (a^2 + 1)*k102 + (a^2 + a)*k103 + (a^2 +
            1)*k002 + (a), k100 + (a^2 + a)*k102 + (a^3 + a^2)*k103 +
            (a^3 + a^2)*k003 + (a^3 + a + 1), k003^2 + k000) of
            Multivariate Polynomial Ring in k100, k101, k102, k103,
            x100, x101, x102, x103, w100, w101, w102, w103, s000,
            s001, s002, s003, k000, k001, k002, k003 over Finite Field
            in a of size 2^4
        """
        return self._ring.ideal(self.gens())

    def groebner_basis(self, *args, **kwargs):
        """
        Compute and return a Groebner basis for self.

        INPUT:
            args -- list of arguments passed to MPolynomialIdeal.groebner_basis call
            kwargs -- dictionary of arguments passed to MPolynomialIdeal.groebner_basis call

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: gb = F.groebner_basis()
            sage: Ideal(gb).basis_is_groebner()
            True
        """
        return self.ideal().groebner_basis(*args, **kwargs)

    def monomials(self):
        """
        Return a list of monomials in self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: len(F.monomials())
            49
        """
        M = set()
        for r in self._rounds:
            for f in r._gens:
                for m in f.monomials():
                    M.add(m)
        return list(M)

    def nmonomials(self):
        """
        Return the number of monomials present in self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.nmonomials()
            49
        """
        return len(self.monomials())

    def variables(self):
        """
        Return all variables present in self. This list may or may not
        be equal to the generators of the ring of self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.variables()[:10]  # this output ordering depends on a hash and so might change if __hash__ changes
            [s000, k101, k100, s003, x102, x103, s002, w103, w102, x100]        # 32-bit
            [k101, k100, s003, x102, x103, s002, w103, w102, x100, x101]        # 64-bit
        """
        V = set()
        for r in self._rounds:
            for f in r._gens:
                for v in f.variables():
                    V.add(v)
        return list(V)


    def nvariables(self):
        """
        Return number of variables present in self.

        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.nvariables()
            20
        """
        return len(self.variables())

    def coefficient_matrix(self, sparse=True):
        """
        Return tuple (A,v) where A is the coefficent matrix of self
        and v the matching monomial vector. Monomials are order w.r.t.
        the term ordering of self.ring() in reverse order.

        INPUT:
            sparse -- construct a sparse matrix (default: True)

        EXAMPLE:
            sage: P.<a,b,c,d> = PolynomialRing(GF(127),4)
            sage: I = sage.rings.ideal.Katsura(P)
            sage: I.gens()
            (a + 2*b + 2*c + 2*d - 1, a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a, 2*a*b + 2*b*c
            + 2*c*d - b, b^2 + 2*a*c + 2*b*d - c)

            sage: F = mq.MPolynomialSystem(I)
            sage: A,v = F.coefficient_matrix()
            sage: A
            [  0   0   0   0   0   0   0   0   0   1   2   2   2 126]
            [  1   0   2   0   0   2   0   0   2 126   0   0   0   0]
            [  0   2   0   0   2   0   0   2   0   0 126   0   0   0]
            [  0   0   1   2   0   0   2   0   0   0   0 126   0   0]

            sage: v
            [a^2]
            [a*b]
            [b^2]
            [a*c]
            [b*c]
            [c^2]
            [b*d]
            [c*d]
            [d^2]
            [  a]
            [  b]
            [  c]
            [  d]
            [  1]

            sage: A*v
            [        a + 2*b + 2*c + 2*d - 1]
            [a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a]
            [      2*a*b + 2*b*c + 2*c*d - b]
            [        b^2 + 2*a*c + 2*b*d - c]

        """
        R = self.ring()

        m = sorted(self.monomials(),reverse=True)
        nm = len(m)
        f = self.gens()
        nf = len(f)

        #construct dictionary for fast lookups
        v = dict( zip( m , range(len(m)) ) )

        A = Matrix( R.base_ring(), nf, nm, sparse=sparse )

        for x in range( nf ):
            poly = f[x]
            for y in poly.monomials():
                A[ x , v[y] ] = poly.monomial_coefficient(y)

        return  A, Matrix(R,nm,1,m)

    def subs(self, *args, **kwargs):
        """
        Substitute variables for every polynomial in self. See
        MPolynomial.subs for calling convention.


        EXAMPLE:
            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system(); F
            Polynomial System with 40 Polynomials in 20 Variables
            sage: F.subs(s); F
            Polynomial System with 40 Polynomials in 16 Variables

        INPUT:
            args -- arguments to be passed to MPolynomial.subs
            kwargs -- keyword arguments to be passed to MPolynomial.subs
        """
        for r in self._rounds:
            r.subs(*args,**kwargs)

    def _singular_(self):
        """
        Return SINGULAR ideal representation of self.
        """
        return singular.ideal(list(self))

    def _magma_(self):
        """
        Return MAGMA ideal representation of self.
        """
        return magma.ideal(list(self))

    def _repr_(self):
        """
        Return a string representation of self.
        """
        return "Polynomial System with %d Polynomials in %d Variables"%(self.ngens(),self.nvariables())

    def __add__(self, right):
        """
        Add polynomial systems together, i.e. create a union of their
        polynomials.
        """
        if is_MPolynomialRoundSystem(right) and right.ring() == self.ring():
            return MPolynomialSystem(self.ring(),self.rounds() + right.rounds())
        elif is_MPolynomialRoundSystem(right) and right.ring() == self.ring():
            return MPolynomialSystem(self.ring(),self.rounds() + [right.gens()])
        else:
            raise TypeError, "right must be MPolynomialRing over same ring as self"

    def __getitem__(self, ij):
        """
        See self.gen().
        """
        if isinstance(ij, tuple):
            i,j = ij
            return self._rounds[i][j]
        elif isinstance(ij, slice):
            return sum(self._rounds,MPolynomialRoundSystem(self.ring(),[]))[ij]
        else:
            ij = int(ij)
            for r in self._rounds:
                if ij >= len(r):
                    ij = ij - len(r)
                else:
                    return r[ij]

    def __contains__(self, element):
        """
        Return True if element is in self or False else.
        """
        for r in self._rounds:
            if element in r:
                return True
        return False

    def __iter__(self):
        r"""
        Return an iterator for \code{self} where all polynomials in
        \code{self} are yielded in order as they appear in \code{self}.

        EXAMPLE:
            sage: P.<x0,x1,x2,x3> = PolynomialRing(GF(37))
            sage: I = sage.rings.ideal.Katsura(P)
            sage: F = mq.MPolynomialSystem(P,I.gens())
            sage: list(F)
            [x0 + 2*x1 + 2*x2 + 2*x3 - 1,
            x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 - x0,
            2*x0*x1 + 2*x1*x2 + 2*x2*x3 - x1,
            x1^2 + 2*x0*x2 + 2*x1*x3 - x2]
        """
        L = []
        for r in self._rounds:
            for f in r:
                yield f

class MPolynomialSystem_gf2(MPolynomialSystem_generic):
    """
    MPolynomialSystem over GF(2).
    """
    def cnf(self):
        """
        Return Canonical Normal Form (CNF) representation of self in a
        format MiniSAT et al. can understand.
        """
        raise NotImplemented

class MPolynomialSystem_gf2e(MPolynomialSystem_generic):
    r"""
    MPolynomialSystem over $GF(2^e)$.
    """

    def change_ring(self, k):
        """
        Project self onto $k$

        INPUT:
            k -- GF(2) (parameter only  for compatible syntax)

        EXAMPLE:
            sage: k.<a> = GF(2^2)
            sage: P.<x,y> = PolynomialRing(k,2)
            sage: a = P.base_ring().gen()
            sage: F = mq.MPolynomialSystem(P,[x*y + 1, a*x + 1])
            sage: F
            Polynomial System with 2 Polynomials in 2 Variables
            sage: F2 = F.change_ring(GF(2)); F2
            Polynomial System with 8 Polynomials in 4 Variables
            sage: F2.gens()
            [x1*y0 + x0*y1 + x1*y1,
            x0*y0 + x1*y1 + 1,
            x0 + x1,
            x1 + 1,
            x0^2 + x0,
            x1^2 + x1,
            y0^2 + y0,
            y1^2 + y1]

        NOTE: Based on SINGULAR implementation by Michael Brickenstein
        <brickenstein@googlemail.com>

        """
        R = self.ring()
        nvars = R.ngens()
        k = R.base_ring()

        helper = PolynomialRing(GF(2), nvars + 1, [str(k.gen())] + map(str, R.gens()), order='lex')
        myminpoly = helper(str(k.polynomial()))

        l = map(lambda x: helper(str(x)), self.gens())

        r = myminpoly.degree()

        intermediate_ring = PolynomialRing(GF(2), nvars*r+1, 'x', order='degrevlex')

        a = intermediate_ring.gen(0)

        # map e.g. x -> a^2*x_2 + a*x_1 + x_0, where x_0,..,x_2 represent
        # the bits of x
        map_ideal = [a]

        var_index=0
        for index in range(nvars):
           _sum=0
           for sum_index in range(r):
              var_index += 1
              _sum += a**sum_index * intermediate_ring.gen(var_index)
           map_ideal.append(_sum)

        myminpoly=myminpoly(*map_ideal)

        l = [f(*map_ideal).reduce([myminpoly]) for f in l]

        result = []
        # split e.g. a^2*x0+a*x1+x2 to x0,x1,x2
        for f in l:
            for i in reversed(range(r)):
               g = f.coefficient(a**i)
               result.append(g)
               f =  f - a**i * g

        # eliminate parameter, change order to lp
        new_var_names = [str(var)+"%d"%i for var in R.gens() for i in range(r)]

        result_ring = PolynomialRing(GF(2), nvars*r,new_var_names, order=R.term_order())

        map_ideal = (0,) + result_ring.gens()
        result = [f(*map_ideal) for f in result]
        result += [e**2 + e for e in result_ring.gens()]

        return MPolynomialSystem(result_ring,result)

