"""
This module implements residue fields for various kinds of rings.

We can take the residue field of prime ideals in maximal order of number fields:

EXAMPLES:
    sage: K.<a> = NumberField(x^3-7)
    sage: P = K.ideal(29).factor()[0][0]
    sage: k = K.residue_field(P)
    sage: k
    Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)
    sage: k.order()
    841

AUTHORS:
    -- David Roe (2007-10-3): initial version
    -- William Stein (2007-12): bug fixes

TESTS:
    sage: K.<z> = CyclotomicField(7)
    sage: P = K.factor_integer(17)[0][0]
    sage: ff = K.residue_field(P)
    sage: a = ff(z)
    sage: parent(a*a)
    Residue field in zbar of Fractional ideal (17)
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"

from sage.rings.field import Field
from sage.rings.integer import Integer
from sage.categories.homset import Hom
from sage.categories.category_types import Fields, Rings
from sage.rings.all import ZZ, QQ, Integers
from sage.rings.number_field.number_field_ideal import is_NumberFieldIdeal
import weakref
from sage.rings.finite_field_givaro import FiniteField_givaro
from sage.rings.finite_field import FiniteField_prime_modn, GF
from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
from sage.structure.parent_base import ParentWithBase

residue_field_cache = {}

def ResidueField(p, names = None, check = True):
    """
    A function that returns the residue class field of a prime ideal p
    of the ring of integers of a number field.

    INPUT:
        p -- a prime integer or prime ideal of an order in a number
             field.
        names -- the variable name for the finite field created.
                 Defaults to the name of the number field variable but
                 with bar placed after it.
        check -- whether or not to check if p is prime.

    OUTPUT:
         -- The residue field at the prime p.

    EXAMPLES:
        sage: K.<a> = NumberField(x^3-7)
        sage: P = K.ideal(29).factor()[0][0]
        sage: ResidueField(P)
        Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)
        sage: k = K.residue_field(P); k
        Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)
        sage: k.order()
        841
    """
    if isinstance(names, tuple):
        if len(names) > 0:
            names = str(names[0])
        else:
            names = None
    key = (p, names)
    if residue_field_cache.has_key(key):
        k = residue_field_cache[key]()
        if k is not None:
            return k
    if PY_TYPE_CHECK(p, Integer):
        if check and not p.is_prime():
            raise ValueError, "p must be prime"
        if names is None:
            names = 'x'
        k = ResidueFiniteField_prime_modn(p, names)
    elif is_NumberFieldIdeal(p):
        if names is None:
            names = '%sbar'%(p.number_field().variable_name())
        if check and not p.is_prime():
            raise ValueError, "p must be prime"
        # Should generalize to allowing residue fields of relative extensions to be extensions of finite fields.
        characteristic = p.smallest_integer()
        K = p.number_field()
        OK = K.maximal_order() # should change to p.order once this works.
        f = OK.number_field().polynomial().change_ring(Integers(characteristic)) #polynomial() should change to absolute_polynomial()
        L = f.factor()
        for i in range(len(L)):
            g = L[i][0]
            a = K(g.change_ring(QQ))
            if a in p:
                break
        else:
            raise RuntimeError, "should have found a prime"
        if g.degree() == 1:
            k = ResidueFiniteField_prime_modn(p, names, im_gen = -g[0], intp = p.smallest_integer())
        else:
            q = characteristic**(g.degree())
            if q < Integer(2)**Integer(16):
                k = ResidueFiniteField_givaro(p, q, names, g, characteristic)
            else:
                k = ResidueFiniteField_ext_pari(p, q, names, g, characteristic)
    else: # Add support for primes in other rings later.
        raise TypeError, "p must be a prime in the integers or a number field"
    residue_field_cache[key] = weakref.ref(k)
    return k

class ResidueField_generic(Field):
    """
    The class representing a generic residue field.
    """
    def __init__(self, p, f, intp):
        """
        INPUT:
           p -- the prime (ideal) defining this residue field
           f -- the morphism from the order to self.
           intp -- the rational prime that p lives over.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3-17)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k = K.residue_field(P) # indirect doctest
        """
        self.p = p
        self.f = f
        if self.f is not None:
            ParentWithBase.__init__(self, GF(intp), coerce_from = [f])

    def __repr__(self):
        """
        Returns a string describing this residue field.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k = K.residue_field(P)
            sage: k
            Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)
        """
        return "Residue field %sof %s"%('in %s '%self.gen() if self.degree() > 1 else '', self.p)

    def lift(self, x):
        """
        Returns a lift of x to the Order, returning a "polynomial" in the generator with coefficients between 0 and p-1.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k =K.residue_field(P)
            sage: OK = K.maximal_order()
            sage: c = OK(a)
            sage: b = k(a)
            sage: k.lift(13*b + 5)
            13*a + 5
            sage: k.lift(12821*b+918)
            3*a + 19
        """
        if self.f is None:
            return x.lift()
        else:
            return self.f.lift(x)

    def __cmp__(self, x):
        """
        Compares two residue fields: they are equal iff the primes defining them are equal.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3-11)
            sage: F = K.ideal(37).factor(); F
            (Fractional ideal (37, a + 12)) * (Fractional ideal (-2*a + 5)) * (Fractional ideal (37, a + 9))
            sage: k =K.residue_field(F[0][0])
            sage: l =K.residue_field(F[1][0])
            sage: k == l
            False
        """
        if type(self) == type(x):
            try:
                return self.p.__cmp__(x.p)
            except AttributeError:
                return -1
        return cmp(type(self), type(x))

cdef class NFResidueFieldHomomorphism(ResidueFieldHomomorphism):
    """
    The class representing a homomorphism from the order of a number
    field to the residue field at a given prime.

    EXAMPLES:
        sage: K.<a> = NumberField(x^3-7)
        sage: P  = K.ideal(29).factor()[0][0]
        sage: k  = K.residue_field(P)
        sage: OK = K.maximal_order()
        sage: abar = k(OK.1); abar
        abar
        sage: (1+abar)^179
        24*abar + 12
        sage: k.coerce_map_from(OK)
        Ring morphism:
          From: Maximal Order in Number Field in a with defining polynomial x^3 - 7
          To:   Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)
    """
    def __init__(self, k, p, im_gen):
        """
        INPUT:
           k -- The residue field that is the codomain of this morphism.
           p -- The prime ideal defining this residue field
           im_gen -- The image of the generator of the number field.

        EXAMPLES:
        We create a residue field homomorphism:
            sage: K.<theta> = CyclotomicField(5)
            sage: P = K.factor_integer(7)[0][0]
            sage: P.residue_class_degree()
            4
            sage: kk.<a> = P.residue_field(); kk
            Residue field in a of Fractional ideal (7)
            sage: phi = kk.coerce_map_from(K.maximal_order()); phi
            Ring morphism:
              From: Maximal Order in Cyclotomic Field of order 5 and degree 4
              To:   Residue field in a of Fractional ideal (7)
            sage: type(phi)
            <type 'sage.rings.residue_field.NFResidueFieldHomomorphism'>

        """
        self.im_gen = im_gen
        self.p = p
        ResidueFieldHomomorphism.__init__(self,Hom(p.number_field().maximal_order(), k, Rings())) # should eventually change to p.order()

    cdef Element _call_c_impl(self, Element x):
        """
        Applies this morphism to an element

        EXAMPLES:
            sage: K.<a> = NumberField(x^3-x+8)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k =K.residue_field(P)
            sage: OK = K.maximal_order()
            sage: k.coerce_map_from(OK)(OK(a)^7)
            13*abar^2 + 7*abar + 21
        """
        y = x.polynomial().change_ring(self.codomain().base_ring())(self.im_gen)  #polynomial should change to absolute_polynomial?
        (<Element>y)._set_parent_c(self.codomain())
        return y

    def lift(self, x):
        """
        Returns a lift of x to the Order, returning a "polynomial" in
        the generator with coefficients between 0 and p-1.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k =K.residue_field(P)
            sage: OK = K.maximal_order()
            sage: f = k.coerce_map_from(OK)
            sage: c = OK(a)
            sage: b = k(a)
            sage: f.lift(13*b + 5)
            13*a + 5
            sage: f.lift(12821*b+918)
            3*a + 19
        """
        if self.domain() is ZZ:
            return x.lift()
        return self.domain()(x.polynomial().change_ring(self.domain().base_ring())(self.domain().ring_generators()[0]))  #polynomial should change to absolute_polynomial?


class ResidueFiniteField_prime_modn(ResidueField_generic, FiniteField_prime_modn):
    """
    The class representing residue fields of number fields that have prime order.

    EXAMPLES:
        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^3-7)
        sage: P = K.ideal(29).factor()[1][0]
        sage: k = ResidueField(P)
        sage: k
        Residue field of Fractional ideal (a^2 + 2*a + 2)
        sage: k.order()
        29
        sage: OK = K.maximal_order()
        sage: c = OK(a)
        sage: b = k(a)
        sage: k.f(c)
        16
        sage: k(4)
        4
        sage: k(c + 5)
        21
        sage: b + c
        3
    """
    def __init__(self, p, name, im_gen = None, intp = None):
        """
        INPUT:
           p -- An integral prime or a prime ideal of a number field.
           name -- the name of the generator of this extension
           im_gen -- the image of the generator of the number field in this finite field.
           intp -- the rational prime that p lies over.

        EXAMPLES:
            sage: k = ResidueField(17); k
            Residue field of 17
            sage: type(k)
            <class 'sage.rings.residue_field.ResidueFiniteField_prime_modn'>
        """
        self.p = p # Here because we have to create a NFResidueFieldHomomorphism before calling ResidueField_generic.__init__(self,...)
        if im_gen is None:
            FiniteField_prime_modn.__init__(self, p, name)
            ResidueField_generic.__init__(self, p, None, p)
        else:
            FiniteField_prime_modn.__init__(self, intp, name)
            self.f = NFResidueFieldHomomorphism(self, p, im_gen)
            ResidueField_generic.__init__(self, p, self.f, intp)

    def __call__(self, x):
        """
        INPUT:
           x -- something to cast in to self.

        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(29).factor()[1][0]
            sage: k = ResidueField(P)
            sage: k
            Residue field of Fractional ideal (a^2 + 2*a + 2)
            sage: OK = K.maximal_order()
            sage: c = OK(a)
            sage: b = k(a); b
            16
        """
        try:
            return self.coerce_map_from(self.f.domain())(self.f.domain()(x))
        except (AttributeError, TypeError):
            return FiniteField_prime_modn.__call__(self, x)

class ResidueFiniteField_ext_pari(ResidueField_generic, FiniteField_ext_pari):
    """
    The class representing residue fields of number fields that have non-prime order >= 2^16.

    EXAMPLES:
        sage: K.<a> = NumberField(x^3-7)
        sage: P = K.ideal(923478923).factor()[0][0]
        sage: k =K.residue_field(P)
        sage: k.degree()
        2
        sage: OK = K.maximal_order()
        sage: c = OK(a)
        sage: b = k(c)
        sage: b+c
        2*abar
        sage: b*c
        664346875*abar + 535606347
    """
    def __init__(self, p, q, name, g, intp):
        """
        EXAMPLES:
        We create an ext_pari residue field:
            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(923478923).factor()[0][0]
            sage: type(P.residue_field())
            <class 'sage.rings.residue_field.ResidueFiniteField_ext_pari'>
        """
        FiniteField_ext_pari.__init__(self, q, name, g)
        self.f = NFResidueFieldHomomorphism(self, p, GF(q, name = name, modulus = g).gen(0))
        ResidueField_generic.__init__(self, p, self.f, intp)

    def __call__(self, x):
        """
        Coerce x into self.

        EXAMPLES:
            sage: K.<aa> = NumberField(x^3 - 2)
            sage: P = K.factor_integer(10007)[0][0]
            sage: P.residue_class_degree()
            2
            sage: ff.<alpha> = P.residue_field(); ff
            Residue field in alpha of Fractional ideal (-12*aa^2 + 189*aa - 475)
            sage: type(ff)
            <class 'sage.rings.residue_field.ResidueFiniteField_ext_pari'>
            sage: ff(alpha^2 + 1)
            7521*alpha + 4131
            sage: ff(17/3)
            6677
        """
        try:
            return self.coerce_map_from(self.f.domain())(self.f.domain()(x))
        except (AttributeError, TypeError):
            return FiniteField_ext_pari.__call__(self, x)

class ResidueFiniteField_givaro(ResidueField_generic, FiniteField_givaro):
    """
    The class representing residue fields of number fields that have non-prime order < 2**16.

    EXAMPLES:
        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^3-7)
        sage: P = K.ideal(29).factor()[0][0]
        sage: k =K.residue_field(P)
        sage: k.degree()
        2
        sage: OK = K.maximal_order()
        sage: c = OK(a)
        sage: b = k(c)
        sage: b*c^2
        7
        sage: b*c
        13*abar + 5
    """
    def __init__(self, p, q, name, g, intp):
        """
        INPUT:
           p -- the prime ideal defining this residue field
           q -- the order of this residue field (a power of intp)
           name -- the name of the generator of this extension
           g -- the polynomial modulus for this extension
           intp -- the rational prime that p lies over.

        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: k =K.residue_field(P)
        """
        FiniteField_givaro.__init__(self, q, name, g)
        self.f = NFResidueFieldHomomorphism(self, p, GF(q, name = name, modulus = g).gen(0))
        ResidueField_generic.__init__(self, p, self.f, intp)

    def __call__(self, x):
        """
        INPUT:
          x -- Something to cast into self.

        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: k =K.residue_field(P)
            sage: k(77*a^7+4)
            2*abar + 4
        """
        try:
            return self.coerce_map_from(self.f.domain())(self.f.domain()(x))
        except (AttributeError, TypeError):
            return FiniteField_givaro.__call__(self, x)

