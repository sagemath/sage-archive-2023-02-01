"""
Multivariate Polynomial Systems.

We call a finite set of multivariate polynomials an ``MPolynomialSystem``.

In many other computer algebra systems (cf. Singular) this class would
be called ``Ideal`` but an ideal is a very distinct object from its
generators and thus this is not an ideal in Sage.

The idea of polynomial systems is specialized to systems which consist
of several "rounds" in Sage. These kind of polynomial systems arise
naturally in algebraic cryptanalysis of symmetric cryptographic
primitives. The most prominent examples of these systems are: the
small scale variants of the AES [CMR05]_ (cf. ``mq.SR``) and
Flurry/Curry [BPW06]_.

AUTHORS:

- Martin Albrecht (2007ff): initial version
- Martin Albrecht (2009): refactoring, clean-up, new functions

EXAMPLES:

As an example consider a small scale variant of the AES::

   sage: sr = mq.SR(2,1,2,4,gf2=True,polybori=True)
   sage: sr
   SR(2,1,2,4)

We can construct a polynomial system for a random plaintext-ciphertext
pair and study it::

   sage: F,s = sr.polynomial_system()
   sage: F
   Polynomial System with 112 Polynomials in 64 Variables

   sage: r2 = F.round(2); r2
   (w200 + k100 + x100 + x102 + x103,
    w201 + k101 + x100 + x101 + x103 + 1,
    w202 + k102 + x100 + x101 + x102 + 1,
    w203 + k103 + x101 + x102 + x103,
    w210 + k110 + x110 + x112 + x113,
    w211 + k111 + x110 + x111 + x113 + 1,
    w212 + k112 + x110 + x111 + x112 + 1,
    w213 + k113 + x111 + x112 + x113,
    w100*x100 + w100*x103 + w101*x102 + w102*x101 + w103*x100,
    w100*x100 + w100*x101 + w101*x100 + w101*x103 + w102*x102 + w103*x101,
    w100*x101 + w100*x102 + w101*x100 + w101*x101 + w102*x100 + w102*x103 + w103*x102,
    w100*x100 + w100*x101 + w100*x103 + w101*x101 + w102*x100 + w102*x102 + w103*x100 + x100,
    w100*x102 + w101*x100 + w101*x101 + w101*x103 + w102*x101 + w103*x100 + w103*x102 + x101,
    w100*x100 + w100*x101 + w100*x102 + w101*x102 + w102*x100 + w102*x101 + w102*x103 + w103*x101 + x102,
    w100*x101 + w101*x100 + w101*x102 + w102*x100 + w103*x101 + w103*x103 + x103,
    w100*x100 + w100*x102 + w100*x103 + w101*x100 + w101*x101 + w102*x102 + w103*x100 + w100,
    w100*x101 + w100*x103 + w101*x101 + w101*x102 + w102*x100 + w102*x103 + w103*x101 + w101,
    w100*x100 + w100*x102 + w101*x100 + w101*x102 + w101*x103 + w102*x100 + w102*x101 + w103*x102 + w102,
    w100*x101 + w100*x102 + w101*x100 + w101*x103 + w102*x101 + w103*x103 + w103,
    w100*x102 + w101*x101 + w102*x100 + w103*x103 + 1,
    w110*x110 + w110*x113 + w111*x112 + w112*x111 + w113*x110,
    w110*x110 + w110*x111 + w111*x110 + w111*x113 + w112*x112 + w113*x111,
    w110*x111 + w110*x112 + w111*x110 + w111*x111 + w112*x110 + w112*x113 + w113*x112,
    w110*x110 + w110*x111 + w110*x113 + w111*x111 + w112*x110 + w112*x112 + w113*x110 + x110,
    w110*x112 + w111*x110 + w111*x111 + w111*x113 + w112*x111 + w113*x110 + w113*x112 + x111,
    w110*x110 + w110*x111 + w110*x112 + w111*x112 + w112*x110 + w112*x111 + w112*x113 + w113*x111 + x112,
    w110*x111 + w111*x110 + w111*x112 + w112*x110 + w113*x111 + w113*x113 + x113,
    w110*x110 + w110*x112 + w110*x113 + w111*x110 + w111*x111 + w112*x112 + w113*x110 + w110,
    w110*x111 + w110*x113 + w111*x111 + w111*x112 + w112*x110 + w112*x113 + w113*x111 + w111,
    w110*x110 + w110*x112 + w111*x110 + w111*x112 + w111*x113 + w112*x110 + w112*x111 + w113*x112 + w112,
    w110*x111 + w110*x112 + w111*x110 + w111*x113 + w112*x111 + w113*x113 + w113,
    w110*x112 + w111*x111 + w112*x110 + w113*x113 + 1)

    sage: type(r2)
    <class 'sage.crypto.mq.mpolynomialsystem.MPolynomialRoundSystem_generic'>

As an example, we separate the system in independent subsystems or compute the coefficient matrix::

    sage: C = mq.MPolynomialSystem(r2).connected_components(); C
    [Polynomial System with 16 Polynomials in 16 Variables,
     Polynomial System with 16 Polynomials in 16 Variables]

    sage: C[0].groebner_basis()
    [w113*x111 + w113*w110 + w113*w111 + x111 + w110 + w111,
     w112*w111 + w113*x110 + w113*x112 + w113*w110 + w113*w111 + w113*w112 + x110 + w113 + 1,
     w112*w110 + w113*x112 + w113*x113 + w113*w110 + w113*w112 + x110 + x111 + x112 + w111 + 1,
     w112*x113 + w113*x110 + w113*x112 + w113*x113 + w113*w110 + w113*w111 + x111 + x112 + x113 + w110,
     w112*x112 + w113*x113 + x110 + x111 + w110 + w111 + 1,
     w112*x111 + w113*x112 + w113*x113 + x112 + w110 + w113,
     w112*x110 + w113*x112 + w113*x113 + w113*w110 + w113*w111 + x110 + x111 + w110 + w111 + w112 + 1,
     w111*w110 + w113*x110 + w113*x112 + w113*x113 + w113*w111 + w113*w112 + x110 + x111 + x112 + x113 + w111,
     w111*x113 + w113*x113 + w113*w110 + w113*w111 + x111 + x112 + w110 + w111 + w112 + 1,
     w111*x112 + w113*x110 + w113*x112 + w113*w110 + w113*w111 + x110 + x111 + x112 + w110 + w112,
     w111*x111 + w113*x110 + w113*x113 + w113*w110 + w113*w111 + x110 + x111 + x112 + x113,
     w111*x110 + w113*x110 + w113*x112 + x110 + x111 + w110 + w113 + 1,
     w110*x113 + w113*x110 + w113*x113 + x111 + x113 + w111 + w113 + 1,
     w110*x112 + w113*x110 + w113*x112 + w113*x113 + x112 + x113 + w110 + w111 + w112,
     w110*x111 + w113*x112 + w113*w110 + w113*w111 + x110 + x112 + x113 + w113,
     w110*x110 + w113*x110 + w113*w110 + w113*w111 + x110 + x113 + w111 + w112 + 1,
     x113*x112 + w113*x112 + w113*x113 + w113*w110 + x111 + w110 + w112 + 1,
     x113*x111 + w113*x110 + w113*x112 + w113*x113 + w113*w112 + x111 + x112 + x113 + w111 + w112,
     x113*x110 + w113*x110 + w113*x112 + w113*w110 + w113*w111 + w113*w112 + x110 + x111 + x113 + w110 + w111 + 1,
     x112*x111 + w113*w110 + x110 + x111 + x113 + w111 + w113,
     x112*x110 + w113*x112 + w113*x113 + w113*w110 + w113*w112 + x110 + x112 + 1,
     x111*x110 + w113*x110 + w113*x112 + w113*x113 + w113*w111 + w113*w112 + x111 + x113 + w110 + w111 + w112,
     w213 + k113 + x111 + x112 + x113,
     w212 + k112 + x110 + x111 + x112 + 1,
     w211 + k111 + x110 + x111 + x113 + 1,
     w210 + k110 + x110 + x112 + x113]

    sage: A,v = mq.MPolynomialSystem(r2).coefficient_matrix()
    sage: A.rank()
    32

Using these building blocks we can implement a simple XL algorithm
easily::

    sage: sr = mq.SR(1,1,1,4,gf2=True,polybori=True,order='lex')
    sage: F,s = sr.polynomial_system()

    sage: monomials = [a*b for a in F.variables() for b in F.variables() if a<b]
    sage: len(monomials)
    190
    sage: F2 = mq.MPolynomialSystem(map(mul, cartesian_product_iterator((monomials, F))))
    sage: A,v = F2.coefficient_matrix(sparse=False)
    sage: A.echelonize()
    sage: A
    6840 x 4474 dense matrix over Finite Field of size 2
    sage: A.rank()
    4056
    sage: A[4055]*v
    (k002*k003)

The last output tells us that ``k002`` and ``k003`` can't both be 1.

TEST::

    sage: P.<x,y> = PolynomialRing(QQ)
    sage: I = [[x^2 + y^2], [x^2 - y^2]]
    sage: F = mq.MPolynomialSystem(P,I)
    sage: loads(dumps(F)) == F
    True

REFERENCES:

.. [BPW06] J\. Buchmann, A\. Pychkine, R\.-P\. Weinmann *Block Ciphers
   Sensitive to Groebner Basis Attacks* in Topics in
   Cryptology –- CT RSA'06\; LNCS 3860\; pp. 313--331\; Springer Verlag 2006\;
   pre-print available at http://eprint.iacr.org/2005/200

"""

from sage.structure.sage_object import SageObject

from sage.rings.integer_ring import ZZ
from sage.rings.finite_field import FiniteField as GF

from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.polynomial.multi_polynomial import is_MPolynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.matrix.matrix import is_Matrix
from sage.matrix.constructor import Matrix

from sage.interfaces.singular import singular

def is_MPolynomialSystem(F):
    """
    Return ``True`` if ``F`` is an ``MPolynomialSystem``.

    INPUT::

    - ``F`` - anything

    EXAMPLE::

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
    Return ``True`` if ``F`` is an ``MPolynomialRoundSystem``.

    INPUT:

    - ``F`` - anything

    EXAMPLE::

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
    Construct an object representing the equations of a single round
    e.g. of a block cipher or however a "round" is defined.

    INPUT:

    -  ``R`` - the base ring
    -  ``gens`` - list (default: ``[]``)

    EXAMPLE::

        sage: P.<x,y,z> = PolynomialRing(GF(2),3)
        sage: mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
        (x*y + 1, z + 1)
    """
    return MPolynomialRoundSystem_generic(R, gens)

def MPolynomialSystem(arg1, arg2=None):
    """
    Construct a new polynomial system object.

    INPUT:

    - ``arg1`` - a multivariate polynomial ring or an ideal
    - ``arg2`` - an iterable object of rounds, preferable ``MPolynomialRoundSystem``,
      or polynomials (default:``None``)

    EXAMPLES::

        sage: P.<a,b,c,d> = PolynomialRing(GF(127),4)
        sage: I = sage.rings.ideal.Katsura(P)

    If a list of MPolynomialRoundSystems is provided, those form the
    rounds::

        sage: mq.MPolynomialSystem(I.ring(), [mq.MPolynomialRoundSystem(I.ring(),I.gens())])
        Polynomial System with 4 Polynomials in 4 Variables

    If an ideal is provided, the generators are used::

        sage: mq.MPolynomialSystem(I)
        Polynomial System with 4 Polynomials in 4 Variables

    If a list of polynomials is provided, the system has only one
    round::

        sage: mq.MPolynomialSystem(I.ring(), I.gens())
        Polynomial System with 4 Polynomials in 4 Variables
    """
    if is_MPolynomialRing(arg1):
        R = arg1
        rounds = arg2
    elif isinstance(arg1, MPolynomialRoundSystem_generic) and arg2 is None:
        R = arg1.ring()
        rounds = [list(arg1)]
    elif is_Matrix(arg1):
        R = arg1.base_ring()
        rounds = arg1.list()
    elif isinstance(arg1,MPolynomialIdeal) and arg2 is None:
        R = arg1.ring()
        rounds = arg1.gens()
    elif isinstance(arg1, (list,tuple)):
        rounds = arg1
        R = iter(rounds).next().parent()
        for f in rounds:
            if f.parent() is not R:
                raise TypeError("Generators must have same parent.")
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
    def __init__(self, R, gens):
        """
        See ``MPolynomialRoundSystem``.

        INPUT:

        -  ``R`` - base ring
        -  ``gens`` - list (default: [])

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(GF(2),3)
            sage: mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
            (x*y + 1, z + 1)
        """
        if is_MPolynomialRing(R):
            self._ring = R
        else:
            raise TypeError, "First parameter must be a MPolynomialRing."

        if is_Matrix(gens):
            self._gens = tuple(gens.list())
        else:
            self._gens = tuple(gens)

    def __copy__(self):
        """
        Return a copy of this round.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: r = F.round(0)
            sage: copy(r)
            (w100 + k000 + (a^3 + a + 1), w101 + k001 + (a^3 + 1),
            w102 + k002 + (a^3 + a^2 + 1), w103 + k003 + (a^3 + a^2 +
            a))
        """
        return MPolynomialRoundSystem_generic(self._ring, list(self._gens))

    def __cmp__(self, other):
        """
        Compare the ring and generators of ``self`` and ``other``.

        INPUT:

        - ``other`` - an ``MPolynomialRoungSystem``

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F == copy(F) # indirect doctest
            True
        """
        return cmp((self._ring, self._gens),(other._ring, other._gens))

    def ring(self):
        """
        Return the polynomial ring the members of this system live in.

        EXAMPLE::

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
        Return number of polynomials in the system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R0 = F.round(0)
            sage: R0.ngens()
            4
        """
        return len(self._gens)

    def gens(self):
        """
        Return list of polynomials in in the system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R1 = F.round(1)
            sage: l = R1.gens()
            sage: l[0]
            k000^2 + k000
        """
        return tuple(self)

    def variables(self):
        """
        Return an unordered list of variables appearing in polynomials
        in this system.

        EXAMPLE::

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
        Return an unordered list of monomials appearing in polynomials
        in this system.

        EXAMPLE::

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
        Substitute variables for every polynomial in this system and
        return a system. See ``MPolynomial.subs`` for calling
        convention.

        INPUT:

        - ``args`` - arguments to be passed to ``MPolynomial.subs``
        - ``kwargs`` - keyword arguments to be passed to ``MPolynomial.subs``

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: R1 = F.round(1)
            sage: R1 = R1.subs(s) # the solution
            sage: R1
            (0, 0, 0, 0)
        """
        return self.__class__(self._ring,[self._gens[i].subs(*args,**kwargs) for i in range(len(self._gens))])

    def _repr_(self):
        """
        Return string representation of this system.

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(GF(2),3)
            sage: F = mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
            sage: str(F) # indirect doctest
            '(x*y + 1, z + 1)'
        """
        return "%s"%(self._gens,)

    def __getitem__(self, i):
        """
        Return the i-th generator of this system.

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(GF(2),3)
            sage: F = mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
            sage: F[0] # indirect doctest
            x*y + 1
        """
        return self._gens[i]

    def __add__(self, right):
        """
        Addition is the union of generators.

        INPUT:

        -  ``right`` - MPolynomialSystem, list or tuple

        EXAMPLE::

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
            return MPolynomialRoundSystem(self.ring(), self._gens + tuple(right))
        else:
            raise ArithmeticError, "Cannot add MPolynomialRoundSystem and %s"%type(right)

    def __contains__(self, element):
        """
        Return True if element is in the list of generators for self.

        EXAMPLE::

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

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(GF(2),3)
            sage: F = mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
            sage: len(F)
            2
        """
        return len(self._gens)

    def __iter__(self):
        """
        Iterate over the generators of this system.

        EXAMPLE::

            sage: P.<x,y,z> = PolynomialRing(GF(2),3)
            sage: F = mq.MPolynomialRoundSystem(P,[x*y +1, z + 1])
            sage: for f in F:
            ...     print f
            x*y + 1
            z + 1
        """
        return iter(self._gens)

    def _singular_(self):
        """
        Return SINGULAR ideal representation of this system.

        EXAMPLE::

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

    def _magma_init_(self, magma):
        """
        Return Magma ideal representation of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True)
            sage: F,s = sr.polynomial_system()
            sage: R1 = F.round(1)
            sage: magma(R1)                               # implicit doctest; optional - magma
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
        P = magma(self.ring()).name()
        v = [x._magma_init_(magma) for x in self._gens]
        return 'ideal<%s|%s>'%(P, ','.join(v))

class MPolynomialSystem_generic(SageObject):
    def __init__(self, R, rounds):
        """
        Construct a new system of multivariate polynomials. That is, a set
        of multivariate polynomials with at least one common root.

        INPUT:


        -  ``arg1`` - a multivariate polynomial ring or an
           ideal

        -  ``arg2`` - an iterable object of rounds, preferably
           MPolynomialRoundSystem, or polynomials (default:None)


        EXAMPLES::

            sage: P.<a,b,c,d> = PolynomialRing(GF(127),4)
            sage: I = sage.rings.ideal.Katsura(P)

        If a list of MPolynomialRoundSystems is provided, those form the
        rounds.

        ::

            sage: mq.MPolynomialSystem(I.ring(), [mq.MPolynomialRoundSystem(I.ring(),I.gens())])
            Polynomial System with 4 Polynomials in 4 Variables

        If an ideal is provided, the generators are used.

        ::

            sage: mq.MPolynomialSystem(I)
            Polynomial System with 4 Polynomials in 4 Variables

        If a list of polynomials is provided, the system has only one
        round.

        ::

            sage: mq.MPolynomialSystem(I.ring(), I.gens())
            Polynomial System with 4 Polynomials in 4 Variables
        """

        self._ring = R
        self._rounds = []

        # check for list of polynomials
        try:
            e = iter(rounds).next()
        except StopIteration:
            self._rounds = list()
            return

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

    def __copy__(self):
        """
	Return a copy of this system.

        While this is not a deep copy, only mutable members of this
	system are copied.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: copy(F) # indirect doctest
            Polynomial System with 40 Polynomials in 20 Variables
        """
        return MPolynomialSystem_generic(self._ring, [r.__copy__() for r in self._rounds])

    def __cmp__(self, other):
        """
        Compare the ring and rounds of this system and ``other``.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F == copy(F) # indirect doctest
            True
        """
        return cmp((self._ring, self._rounds),(other._ring, other._rounds))

    def ring(self):
        """
        Return base ring.

        EXAMPLE::

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
        Return number polynomials in this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True,order='block')
            sage: F,s = sr.polynomial_system()
            sage: F.ngens()
            56
        """
        return sum([e.ngens() for e in self._rounds])

    def gens(self):
        """
        Return tuple of polynomials in this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: l = F.gens()
            sage: len(l), type(l)
            (40, <type 'tuple'>)
        """
        return tuple(self)

    def gen(self, ij):
        """
        Return an element in this system indexed by ``ij``.

        INPUT:

        -  ``ij`` - tuple, slice, integer


        EXAMPLES::

            sage: P.<a,b,c,d> = PolynomialRing(GF(127),4)
            sage: F = mq.MPolynomialSystem(sage.rings.ideal.Katsura(P))

        ``ij``-th polynomial overall::

            sage: F[0] # indirect doctest
            a + 2*b + 2*c + 2*d - 1

        ``i``-th to ``j``-th polynomial overall::

            sage: F[0:2]
            (a + 2*b + 2*c + 2*d - 1, a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a)

        ``i``-th round, ``j``-th polynomial::

            sage: F[0,1]
            a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a
        """
        return self[ij]

    def nrounds(self):
        """
        Return number of rounds of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.nrounds()
            4
        """
        return len(self._rounds)

    def rounds(self):
        """
        Return a tuple of rounds of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: l = F.rounds()
            sage: len(l)
            4
        """
        return tuple(self._rounds)

    def round(self, i):
        """
        Return ``i``-th round of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: R0 = F.round(1)
            sage: R0
            (k000^2 + k001, k001^2 + k002, k002^2 + k003, k003^2 + k000)
        """
        return self._rounds[i]

    def __iter__(self):
        """
        Iterate over the generators of self round by round.

        EXAMPLE::

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
        Return Sage ideal spanned by the elements of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: P = F.ring()
            sage: I = F.ideal()
            sage: I.elimination_ideal(P('s000*s001*s002*s003*w100*w101*w102*w103*x100*x101*x102*x103'))
            Ideal (k002 + (a^2)*k003 + 1,
                   k001 + (a^3)*k003 + (a + 1),
                   k000 + (a^3 + a^2 + a)*k003 + (a^3 + a^2 + a),
                   k103 + (a^3)*k003 + (a^2 + 1),
                   k102 + (a^3 + a^2 + a)*k003 + (a^3 + 1),
                   k101 + k003 + (a^2 + a),
                   k100 + (a^2)*k003 + (a^2 + a),
                   k003^2 + (a^3 + a^2 + a)*k003 + (a^3 + a^2 + a))
            of Multivariate Polynomial Ring in k100, k101, k102, k103,
            x100, x101, x102, x103, w100, w101, w102, w103, s000,
            s001, s002, s003, k000, k001, k002, k003 over Finite Field
            in a of size 2^4
        """
        return self._ring.ideal(self.gens())

    def groebner_basis(self, *args, **kwargs):
        """
        Compute and return a Groebner basis for the ideal spanned by
        the polynomials in this system (``self.gens()``).

        INPUT:

        - ``args`` - list of arguments passed to
           ``MPolynomialIdeal.groebner_basis`` call

        - ``kwargs`` - dictionary of arguments passed to
           ``MPolynomialIdeal.groebner_basis`` call

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: gb = F.groebner_basis()
            sage: Ideal(gb).basis_is_groebner()
            True
        """
        return self.ideal().groebner_basis(*args, **kwargs)

    def monomials(self):
        """
        Return an unordered tuple of monomials in this polynomial system.

        EXAMPLE::

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
        return tuple(M)

    def nmonomials(self):
        """
        Return the number of monomials present in this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.nmonomials()
            49
        """
        return len(self.monomials())

    def variables(self):
        """
        Return all variables present in this system. This tuple may or
        may not be equal to the generators of the ring of this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.variables()[:10]
            (k003, k002, k001, k000, s003, s002, s001, s000, w103, w102)
        """
        V = set()
        for r in self._rounds:
            for f in r._gens:
                for v in f.variables():
                    V.add(v)
        return tuple(sorted(V))

    def nvariables(self):
        """
        Return number of variables present in this system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: F.nvariables()
            20
        """
        return len(self.variables())

    def coefficient_matrix(self, sparse=True):
        """
        Return tuple ``(A,v)`` where ``A`` is the coefficient matrix
        of this system and ``v`` the matching monomial vector.

        Thus value of ``A[i,j]`` corresponds the coefficient of the
        monomial ``v[j]`` in the ``i``-th polynomial in this system.

        Monomials are order w.r.t. the term ordering of
        ``self.ring()`` in reverse order, i.e. such that the smallest
        entry comes last.

        INPUT:

        - ``sparse`` - construct a sparse matrix (default: ``True``)

        EXAMPLE::

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
        Substitute variables for every polynomial in this system and
        return a new system. See ``MPolynomial.subs`` for calling
        convention.

        INPUT:

        -  ``args`` - arguments to be passed to ``MPolynomial.subs``
        -  ``kwargs`` - keyword arguments to be passed to ``MPolynomial.subs``

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system(); F
            Polynomial System with 40 Polynomials in 20 Variables
            sage: F = F.subs(s); F
            Polynomial System with 40 Polynomials in 16 Variables

        """
        return MPolynomialSystem(self._ring, [r.subs(*args,**kwargs) for r in self._rounds])

    def _singular_(self):
        """
        Return Singular ideal representation of this system.

        EXAMPLE::

            sage: P.<a,b,c,d> = PolynomialRing(GF(127))
            sage: I = sage.rings.ideal.Katsura(P)
            sage: F = mq.MPolynomialSystem(I); F
            Polynomial System with 4 Polynomials in 4 Variables
            sage: F._singular_()
            a+2*b+2*c+2*d-1,
            a^2+2*b^2+2*c^2+2*d^2-a,
            2*a*b+2*b*c+2*c*d-b,
            b^2+2*a*c+2*b*d-c
        """
        return singular.ideal(list(self))

    def _magma_init_(self, magma):
        """
        Return Magma ideal representation of the ideal spanned by this
        system.

        EXAMPLE::

            sage: sr = mq.SR(allow_zero_inversions=True,gf2=True)
            sage: F,s = sr.polynomial_system()
            sage: magma(F)                          # implicit doctest; optional - magma
            Ideal of Polynomial ring of rank 20 over GF(2)
            Graded Reverse Lexicographical Order
            Variables: k100, k101, k102, k103, x100, x101, x102, x103, w100, w101, w102, w103, s000, s001, s002, s003, k000, k001, k002, k003
            Basis:
            [
            ...
            ]
        """
        P = magma(self.ring()).name()
        v = [x._magma_init_(magma) for x in list(self)]
        return 'ideal<%s|%s>'%(P, ','.join(v))

    def _repr_(self):
        """
        Return a string representation of this system.

        EXAMPLE::

            sage: P.<a,b,c,d> = PolynomialRing(GF(127))
            sage: I = sage.rings.ideal.Katsura(P)
            sage: F = mq.MPolynomialSystem(I); F # indirect doctest
            Polynomial System with 4 Polynomials in 4 Variables
        """
        return "Polynomial System with %d Polynomials in %d Variables"%(self.ngens(),self.nvariables())

    def __add__(self, right):
        """
        Add polynomial systems together, i.e. create a union of their
        polynomials.

        EXAMPLE::

            sage: P.<a,b,c,d> = PolynomialRing(GF(127))
            sage: I = sage.rings.ideal.Katsura(P)
            sage: F = mq.MPolynomialSystem(I)
            sage: F + [a^127 + a]
            Polynomial System with 5 Polynomials in 4 Variables

            sage: F + P.ideal([a^127 + a])
            Polynomial System with 5 Polynomials in 4 Variables

            sage: F + mq.MPolynomialSystem(P,[a^127 + a])
            Polynomial System with 5 Polynomials in 4 Variables
        """
        if is_MPolynomialSystem(right) and right.ring() == self.ring():
            return MPolynomialSystem(self.ring(),self.rounds() + right.rounds())
        elif is_MPolynomialRoundSystem(right) and right.ring() == self.ring():
            return MPolynomialSystem(self.ring(),self.rounds() + [right.gens()])
        elif isinstance(right,(tuple,list)) and all(map(lambda x: x.parent() == self.ring(), right)):
            return MPolynomialSystem(self.ring(),self.rounds() + (right,))
        elif isinstance(right,MPolynomialIdeal) and right.ring() == self.ring():
            return MPolynomialSystem(self.ring(),self.rounds() + (right.gens(),))
        else:
            raise TypeError, "right must be a system over same ring as self"

    def __getitem__(self, ij):
        """
        See ``self.gen()``.

        EXAMPLE::

            sage: P.<a,b,c,d> = PolynomialRing(GF(127),4)
                   sage: F = mq.MPolynomialSystem(sage.rings.ideal.Katsura(P))

        `ij`-th polynomial overall::

            sage: F[0] # indirect doctest
            a + 2*b + 2*c + 2*d - 1

        ``i``-th to ``j``-th polynomial overall::

            sage: F[0:2]
            (a + 2*b + 2*c + 2*d - 1, a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a)

        ``i``-th round, ``j``-th polynomial::

            sage: F[0,1]
            a^2 + 2*b^2 + 2*c^2 + 2*d^2 - a
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
        Return ``True`` if element is in ``self`` or ``False``
        otherwise. This method does not return an answer for the ideal
        spanned by the generators of this system but literately
        whether a polynomial is in the list of generators.

        EXAMPLE::

            sage: P.<x0,x1,x2,x3> = PolynomialRing(GF(37))
            sage: I = sage.rings.ideal.Katsura(P)
            sage: F = mq.MPolynomialSystem(I)
            sage: f = x0 + 2*x1 + 2*x2 + 2*x3 -1
            sage: f in F
            True
            sage: x0*f in F
            False
            sage: x0*f in F.ideal()
            True
        """
        for r in self._rounds:
            if element in r:
                return True
        return False

    def __iter__(self):
        """
        Return an iterator for this polynomial system where all
        polynomials are yielded in order as they appear.

        EXAMPLE::

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

    def connection_graph(self):
        """
        Return the graph which has the variables of this system as
        vertices and edges between two variables if they appear in the
        same polynomial.

        EXAMPLE::

            sage: B.<x,y,z> = BooleanPolynomialRing()
            sage: F = mq.MPolynomialSystem([x*y + y + 1, z + 1])
            sage: F.connection_graph()
            Graph on 3 vertices
        """
        V = sorted(self.variables())
        from sage.graphs.graph import Graph
        g = Graph()
        g.add_vertices(sorted(V))
        for f in self:
            v = f.variables()
            a,tail = v[0],v[1:]
            for b in tail:
                g.add_edge((a,b))
        return g

    def connected_components(self):
        """
        Split the polynomial system in systems which do not share any
        variables.

        EXAMPLE:

        As an example consider one round of AES, which naturally
        splits into four subsystems which are independent::

            sage: sr = mq.SR(2,4,4,8,gf2=True,polybori=True)
            sage: F,s = sr.polynomial_system()
            sage: Fz = mq.MPolynomialSystem(F.round(2))
            sage: Fz.connected_components()
            [Polynomial System with 128 Polynomials in 128 Variables,
             Polynomial System with 128 Polynomials in 128 Variables,
             Polynomial System with 128 Polynomials in 128 Variables,
             Polynomial System with 128 Polynomials in 128 Variables]
        """
        g = self.connection_graph()
        C = sorted(g.connected_components())

        P = [[] for _ in xrange(len(C))]
        for f in self:
            for i,c in enumerate(C):
                if len(set(f.variables()).difference(c)) == 0:
                    P[i].append(f)
                    break
        P = sorted([MPolynomialSystem(sorted(p)) for p in P])
        return P


class MPolynomialSystem_gf2(MPolynomialSystem_generic):
    """
    Polynomial Systems over `\mathbb{F}_2`.
    """
    #def cnf(self):
    #    """
    #    Return Canonical Normal Form (CNF) representation of self in a
    #    format MiniSAT et al. can understand.
    #    """
    #    raise NotImplemented
    def eliminate_linear_variables(self, maxlength=3, skip=lambda lm,tail: False):
        """
        Return a new system where "linear variables" are eliminated.

        In this function we call a variable "linear" if it appears as
        a leading term of a linear polynomial.

        INPUT:

        - ``maxlength`` - an optional upper bound on the number of
          monomials by which a variable is replaced.

        - ``skip`` - an optional callable to skip eliminations. It
          must accept two parameters and return either ``True`` or
          ``False``. The two parameters are the leading term and the
          tail of a polynomial (default: ``lambda lm,tail: False``).

        EXAMPLE::

            sage: B.<a,b,c,d> = BooleanPolynomialRing()
            sage: F = mq.MPolynomialSystem([c + d + b + 1, a + c + d, a*b + c, b*c*d + c])
            sage: F.eliminate_linear_variables().gens() # everything vanishes
            ()
            sage: F.eliminate_linear_variables(maxlength=2).gens()
            (b + c + d + 1, b*c + b*d + c, b*c*d + c)
            sage: F.eliminate_linear_variables(skip=lambda lm,tail: str(lm)=='a').gens()
            (a + c + d, a*c + a*d + a + c, c*d + c)

        .. note::

            This is called "massaging" in [CB07]_.

        REFERENCES:

        .. [CB07] Nicolas T\. Courtois, Gregory V\. Bard *Algebraic
           Cryptanalysis of the Data Encryption Standard*; Cryptography
           and Coding -- 11th IMA International Conference; 2007;
           available at http://eprint.iacr.org/2006/402

        """
        from polybori.ll import ll_encode
        from polybori.ll import ll_red_nf
        from sage.rings.polynomial.pbori import set_cring, BooleanPolynomialRing
        from sage.misc.misc import get_verbose

        R = self.ring()

        if not isinstance(R, BooleanPolynomialRing):
            raise NotImplementedError("Only BooleanPolynomialRing's are supported.")

        set_cring(R)

        F = self

        while True:
            linear = []
            higher = []

            for f in F:
                if f.degree() == 1 and len(f) <= maxlength + 1:
                    flm = f.lexLead()
                    if skip(flm, f-flm):
                        higher.append(f)
                        continue
                    lex_lead = map(lambda x: x.lexLead(), linear)
                    if not flm in lex_lead:
                        linear.append(f)
                    else:
                        higher.append(f)
                else:
                    higher.append(f)

            if linear == []:
                break

            assert len(set(linear)) == len(linear)

            rb = ll_encode(linear)

            F = []

            for f in linear:
                f = ll_red_nf(f, rb)
                if f:
                    F.append(f)

            for f in higher:
                f = ll_red_nf(f, rb)
                if f:
                    F.append(f)
            if get_verbose() > 0:
                print ".",
        if get_verbose() > 0:
            print
        return MPolynomialSystem(R, higher)

class MPolynomialSystem_gf2e(MPolynomialSystem_generic):
    """
    MPolynomialSystem over `\mathbb{F}_{2^e}`.
    """
    def change_ring(self, k):
        """
        Project self onto ``k`` using the Weil restriction of scalars.

        INPUT:

        -  ``k`` - GF(2) (parameter only for compatible syntax)

        EXAMPLE::

            sage: k.<a> = GF(2^2)
            sage: P.<x,y> = PolynomialRing(k,2)
            sage: a = P.base_ring().gen()
            sage: F = mq.MPolynomialSystem(P,[x*y + 1, a*x + 1])
            sage: F
            Polynomial System with 2 Polynomials in 2 Variables
            sage: F2 = F.change_ring(GF(2)); F2
            doctest... DeprecationWarning: The use of this function is deprecated please use the weil_restriction() function instead.
            Polynomial System with 8 Polynomials in 4 Variables
            sage: F2.gens()
            (x1*y0 + x0*y1 + x1*y1,
            x0*y0 + x1*y1 + 1,
            x0 + x1,
            x1 + 1,
            x0^2 + x0,
            x1^2 + x1,
            y0^2 + y0,
            y1^2 + y1)

        .. note::

            This function is deprecated use the :meth:`.weil_restriction()` function instead.
        """
        from sage.misc.misc import deprecation
        deprecation("The use of this function is deprecated please use the weil_restriction() function instead.")
        return self.weil_restriction()

    def weil_restriction(self):
        """
        Project this polynomial system to `\mathbb{F}_2`.

        That is, compute the Weil restriction of scalars for the
        variety corresponding to this polynomial system and express it
        as a polynomial system over `\mathbb{F}_2`.

        EXAMPLE::

            sage: k.<a> = GF(2^2)
            sage: P.<x,y> = PolynomialRing(k,2)
            sage: a = P.base_ring().gen()
            sage: F = mq.MPolynomialSystem(P,[x*y + 1, a*x + 1])
            sage: F
            Polynomial System with 2 Polynomials in 2 Variables
            sage: F2 = F.weil_restriction(); F2
            Polynomial System with 8 Polynomials in 4 Variables
            sage: F2.gens()
            (x1*y0 + x0*y1 + x1*y1,
            x0*y0 + x1*y1 + 1,
            x0 + x1,
            x1 + 1,
            x0^2 + x0,
            x1^2 + x1,
            y0^2 + y0,
            y1^2 + y1)

        Another bigger example for a small scale AES::

           sage: sr = mq.SR(1,1,1,4,gf2=False)
           sage: F,s = sr.polynomial_system(); F
           Polynomial System with 40 Polynomials in 20 Variables
           sage: F2 = F.weil_restriction(); F2
           Polynomial System with 240 Polynomials in 80 Variables
        """
        from sage.rings.ideal import FieldIdeal
        J = self.ideal().weil_restriction()
        J += FieldIdeal(J.ring())
        return MPolynomialSystem(J)
