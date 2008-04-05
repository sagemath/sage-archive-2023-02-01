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

We reduce mod a prime for which the ring of integers is not
monogenic (i.e., 2 is an essential discriminant divisor):
    sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
    sage: F = K.factor(2); F
    (Fractional ideal (1/2*a^2 - 1/2*a + 1)) * (Fractional ideal (a^2 - 2*a + 3)) * (Fractional ideal (3/2*a^2 - 5/2*a + 4))
    sage: F[0][0].residue_field()
    Residue field of Fractional ideal (1/2*a^2 - 1/2*a + 1)
    sage: F[1][0].residue_field()
    Residue field of Fractional ideal (a^2 - 2*a + 3)
    sage: F[2][0].residue_field()
    Residue field of Fractional ideal (3/2*a^2 - 5/2*a + 4)

AUTHORS:
    -- David Roe (2007-10-3): initial version
    -- William Stein (2007-12): bug fixes

TESTS:
    sage: K.<z> = CyclotomicField(7)
    sage: P = K.factor(17)[0][0]
    sage: ff = K.residue_field(P)
    sage: a = ff(z)
    sage: parent(a*a)
    Residue field in zbar of Fractional ideal (17)

Reducing a curve modulo a prime:
    sage: K.<s> = NumberField(x^2+23)
    sage: OK = K.ring_of_integers()
    sage: E = EllipticCurve([0,0,0,K(1),K(5)])
    sage: pp = K.factor(13)[0][0]
    sage: Fpp = OK.residue_field(pp)
    sage: E.base_extend(Fpp)
    Elliptic Curve defined by y^2  = x^3 + x + 5 over Residue field of Fractional ideal (13, s - 4)
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
from sage.rings.finite_field import FiniteField as GF
from sage.rings.finite_field_givaro import FiniteField_givaro
from sage.rings.finite_field_prime_modn import FiniteField_prime_modn
from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
from sage.structure.parent_base import ParentWithBase

from sage.modules.free_module_element import FreeModuleElement

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

residue_field_cache = {}

def ResidueField(p, names = None, check = True):
    """
    A function that returns the residue class field of a prime ideal p
    of the ring of integers of a number field.

    INPUT:
        p -- a prime ideal of an order in a number field.
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

    The result is cached:
        sage: ResidueField(P) is ResidueField(P)
        True
        sage: k = K.residue_field(P); k
        Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)
        sage: k.order()
        841

    An example where the generator of the number field doesn't
    generate the residue class field.
        sage: K.<a> = NumberField(x^3-875)
        sage: P = K.ideal(5).factor()[0][0]; k = K.residue_field(P); k
        Residue field in abar of Fractional ideal (5, -2/25*a^2 - 1/5*a + 2)
        sage: k.polynomial()
        abar^2 + 3*abar + 4
        sage: k.0^3 - 875
        2

    An example where the residue class field is large but of degree 1:
        sage: K.<a> = NumberField(x^3-875); P = K.ideal(2007).factor()[0][0]; k = K.residue_field(P); k
        Residue field of Fractional ideal (-2/25*a^2 - 2/5*a - 3)
        sage: k(a)
        168
        sage: k(a)^3 - 875
        0

    In this example, 2 is an inessential discriminant divisor, so divides
    the index of ZZ[a] in the maximal order for all a.
        sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8); P = K.ideal(2).factor()[0][0]; P
        Fractional ideal (1/2*a^2 - 1/2*a + 1)
        sage: F = K.residue_field(P); F
        Residue field of Fractional ideal (1/2*a^2 - 1/2*a + 1)
        sage: F(a)
        0
        sage: B = K.maximal_order().basis(); B
        [1, 1/2*a^2 + 1/2*a, a^2]
        sage: F(B[1])
        1
        sage: F(B[2])
        0
        sage: F
        Residue field of Fractional ideal (1/2*a^2 - 1/2*a + 1)
        sage: F.degree()
        1
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
    if check:
        if not is_NumberFieldIdeal(p):
            raise TypeError, "p must be a prime ideal in the ring of integers of a number field."
        if not p.is_prime():
            raise ValueError, "p must be prime"

    if names is None:
        names = '%sbar'%(p.number_field().variable_name())
    # Should generalize to allowing residue fields of relative extensions to be extensions of finite fields.
    characteristic = p.smallest_integer()

    K = p.number_field()
    OK = K.maximal_order() # should change to p.order once this works.

    U, to_vs, to_order = p._p_quotient(characteristic)
    k = U.base_ring()
    R = PolynomialRing(k, names)
    n = p.residue_class_degree()
    gen_ok = False
    from sage.matrix.constructor import matrix
    try:
        x = K.gen()
        M = matrix(k, n+1, n, [to_vs(x**i).list() for i in range(n+1)])
        W = M.transpose().echelon_form()
        if M.rank() == n:
            PB = M.matrix_from_rows(range(n))
            gen_ok = True
            f = R((-W.column(n)).list() + [1])
    except (TypeError, ZeroDivisionError):
        pass
    if not gen_ok:
        bad = True
        for u in U: # using this iterator may not be optimal, we may get a long string of non-generators
            if u:
                x = to_order(u)
                M = matrix(k, n+1, n, [to_vs(x**i).list() for i in range(n+1)])
                W = M.transpose().echelon_form()
                if W.rank() == n:
                    f = R((-W.column(n)).list() + [1])
                    PB = M.matrix_from_rows(range(n))
                    bad = False
                    break
        assert not bad, "error -- didn't find a generator."
    if n == 1:
        k = ResidueFiniteField_prime_modn(p, names, im_gen = -f[0], intp = p.smallest_integer())
    else:
        q = characteristic**(f.degree())
        if q < Integer(2)**Integer(16):
            k = ResidueFiniteField_givaro(p, q, names, f, characteristic)
        else:
            k = ResidueFiniteField_ext_pari(p, q, names, f, characteristic)
    # end creating field.

    # The reduction map is just x |--> k(to_vs(x) * (PB**(-1)))
    # The lifting map is just x |--> to_order(x * PB)
    pi = ReductionMap(K, k, to_vs, PB**(-1))
    lift = LiftingMap(K, k, to_order, PB)
    k._structure = (pi, lift)

    residue_field_cache[key] = weakref.ref(k)
    return k

class ResidueField_generic(Field):
    """
    The class representing a generic residue field.

    EXAMPLES:
        sage: I = QQ[i].factor(2)[0][0]; I
        Fractional ideal (I + 1)
        sage: k = I.residue_field(); k
        Residue field of Fractional ideal (I + 1)
        sage: type(k)
        <class 'sage.rings.residue_field.ResidueFiniteField_prime_modn'>
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
        Returns a lift of x to the Order, returning a "polynomial" in the
        generator with coefficients between 0 and $p-1$.

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

    def reduction_map(self):
        """
        Return the partially defined reduction map from the number
        field to this residue class field.

        EXAMPLES:
            sage: I = QQ[2^(1/3)].factor(2)[0][0]; I
            Fractional ideal (-a)
            sage: k = I.residue_field(); k
            Residue field of Fractional ideal (-a)
            sage: pi = k.reduction_map(); pi
            Partially defined reduction map from Number Field in a with defining polynomial x^3 - 2 to Residue field of Fractional ideal (-a)
            sage: pi.domain()
            Number Field in a with defining polynomial x^3 - 2
            sage: pi.codomain()
            Residue field of Fractional ideal (-a)
        """
        return self._structure[0]

    def lift_map(self):
        """
        EXAMPLES:
            sage: I = QQ[3^(1/3)].factor(5)[1][0]; I
            Fractional ideal (-a + 2)
            sage: k = I.residue_field(); k
            Residue field of Fractional ideal (-a + 2)
            sage: f = k.lift_map(); f
            Lifting map from Residue field of Fractional ideal (-a + 2) to Number Field in a with defining polynomial x^3 - 3
            sage: f.domain()
            Residue field of Fractional ideal (-a + 2)
            sage: f.codomain()
            Number Field in a with defining polynomial x^3 - 3
            sage: f(k.0)
            1
        """
        return self._structure[1]

    def __cmp__(self, x):
        """
        Compares two residue fields: they are equal iff the primes
        defining them are equal.

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

class ReductionMap:
    """
    A reduction map from a (subset) of a number field to this residue
    class field.

    EXAMPLES:
        sage: I = QQ[sqrt(17)].factor(5)[0][0]; I
        Fractional ideal (5)
        sage: k = I.residue_field(); k
        Residue field in sqrt17bar of Fractional ideal (5)
        sage: R = k.reduction_map(); R
        Partially defined reduction map from Number Field in sqrt17 with defining polynomial x^2 - 17 to Residue field in sqrt17bar of Fractional ideal (5)
    """
    def __init__(self, K, F, to_vs, PBinv):
        """
        Create a reduction map.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: F = K.factor(2)[0][0].residue_field()
            sage: F.reduction_map()
            Partially defined reduction map from Number Field in a with defining polynomial x^3 + x^2 - 2*x + 8 to Residue field of Fractional ideal (1/2*a^2 - 1/2*a + 1)
        """
        self.__K = K
        self.__F = F   # finite field
        self.__to_vs = to_vs
        self.__PBinv = PBinv

    def domain(self):
        """
        Return the domain of this reduction map.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 32)
            sage: F = K.factor(2)[0][0].residue_field()
            sage: F.reduction_map().domain()
            Number Field in a with defining polynomial x^3 + x^2 - 2*x + 32
        """
        return self.__K

    def codomain(self):
        """
        Return the codomain of this reduction map.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3 + 128)
            sage: F = K.factor(2)[0][0].residue_field()
            sage: F.reduction_map().codomain()
            Residue field of Fractional ideal (-1/4*a)
        """
        return self.__F

    def __call__(self, x):
        """
        Apply this reduction map to an element that coerces into the number field.

        If x doesn't map because the denominator is not coprime to the
        prime ideal, then a ZeroDivisionError exception is raised.

        EXAMPLES:
            sage: K.<a> = NumberField(x^2 + 1)
            sage: F = K.factor(2)[0][0].residue_field()
            sage: r = F.reduction_map(); r
            Partially defined reduction map from Number Field in a with defining polynomial x^2 + 1 to Residue field of Fractional ideal (a + 1)
            sage: r(2 + a)
            1
            sage: r(a/2)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.
        """
        # The reduction map is just x |--> F(to_vs(x) * (PB**(-1)))
        x = self.__K(x)
        return self.__F(self.__to_vs(x) * self.__PBinv)

    def __repr__(self):
        """
        EXAMPLES:
            sage: K.<theta_5> = CyclotomicField(5)
            sage: F = K.factor(7)[0][0].residue_field()
            sage: F.reduction_map().__repr__()
            'Partially defined reduction map from Cyclotomic Field of order 5 and degree 4 to Residue field in theta_5bar of Fractional ideal (7)'
        """
        return "Partially defined reduction map from %s to %s"%(self.__K, self.__F)

class LiftingMap:
    """
    Lifting map from residue class field to number field.

    EXAMPLES:
        sage: K.<a> = NumberField(x^3 + 2)
        sage: F = K.factor(5)[0][0].residue_field()
        sage: F.degree()
        2
        sage: L = F.lift_map(); L
        Lifting map from Residue field in abar of Fractional ideal (a^2 + 2*a - 1) to Number Field in a with defining polynomial x^3 + 2
        sage: L(F.0^2)
        3*a + 1
        sage: L(3*a + 1) == F.0^2
        True
    """
    def __init__(self, K, F, to_order, PB):
        """
        Create a lifting map.

        EXAMPLES:
            sage: K.<theta_5> = CyclotomicField(5)
            sage: F = K.factor(7)[0][0].residue_field()
            sage: F.lift_map()
            Lifting map from Residue field in theta_5bar of Fractional ideal (7) to Cyclotomic Field of order 5 and degree 4
        """
        self.__K = K
        self.__F = F   # finite field
        self.__to_order = to_order
        self.__PB = PB

    def domain(self):
        """
        Return the domain of this lifting map.

        EXAMPLES:
            sage: K.<a> = NumberField(x^5 + 2)
            sage: F = K.factor(7)[0][0].residue_field()
            sage: L = F.lift_map(); L
            Lifting map from Residue field in abar of Fractional ideal (-2*a^4 + a^3 - 4*a^2 + 2*a - 1) to Number Field in a with defining polynomial x^5 + 2
            sage: L.domain()
            Residue field in abar of Fractional ideal (-2*a^4 + a^3 - 4*a^2 + 2*a - 1)
        """
        return self.__F

    def codomain(self):
        """
        Return the codomain of this lifting map.

        EXAMPLES:
            sage: K.<a> = CyclotomicField(7)
            sage: F = K.factor(5)[0][0].residue_field()
            sage: L = F.lift_map(); L
            Lifting map from Residue field in abar of Fractional ideal (5) to Cyclotomic Field of order 7 and degree 6
            sage: L.codomain()
            Cyclotomic Field of order 7 and degree 6
        """
        return self.__K

    def __call__(self, x):
        """
        Lift from this residue class field to the number field.

        EXAMPLES:
            sage: K.<a> = CyclotomicField(7)
            sage: F = K.factor(5)[0][0].residue_field()
            sage: L = F.lift_map(); L
            Lifting map from Residue field in abar of Fractional ideal (5) to Cyclotomic Field of order 7 and degree 6
            sage: L(F.0)
            a
            sage: F(a)
            abar
        """
        # The lifting map is just x |--> to_order(x * PB)
        x = self.__F(x)
        v = x.polynomial().padded_list(self.__F.degree())
        return self.__to_order(self.__PB.linear_combination_of_rows(v))

    def __repr__(self):
        """
        EXAMPLES:
            sage: K.<theta_12> = CyclotomicField(12)
            sage: F.<tmod> = K.factor(7)[0][0].residue_field()
            sage: F.lift_map().__repr__()
            'Lifting map from Residue field in tmod of Fractional ideal (-3*theta_12^2 + 1) to Cyclotomic Field of order 12 and degree 4'
        """
        return "Lifting map from %s to %s"%(self.__F, self.__K)

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
            sage: P = K.factor(7)[0][0]
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
        #y = x.polynomial().change_ring(self.codomain().base_ring())(self.im_gen)  #polynomial should change to absolute_polynomial?
        y = self.codomain()(x)
        (<Element>y)._set_parent_c(self.codomain())
        return y

    def lift(self, x):
        """
        Returns a lift of x to the Order, returning a "polynomial" in
        the generator with coefficients between 0 and p-1.

        EXAMPLES:
            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k = K.residue_field(P)
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
        else:
            return self.codomain()._structure[1](x)

        # return self.domain()(x.polynomial().change_ring(self.domain().base_ring())(self.domain().ring_generators()[0]))  #polynomial should change to absolute_polynomial?


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
           p -- A prime ideal of a number field.
           name -- the name of the generator of this extension
           im_gen -- the image of the generator of the number field in this finite field.
           intp -- the rational prime that p lies over.

        EXAMPLES:
            sage: K.<i> = QuadraticField(-1)
            sage: kk = ResidueField(K.factor(5)[0][0])
            sage: type(kk)
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
            return FiniteField_prime_modn.__call__(self, x)
        except TypeError:
            if isinstance(x, FreeModuleElement):
                return FiniteField_prime_modn.__call__(self, x[0])
            else:
                return self._structure[0](x)
        #try:
        #    return self.coerce_map_from(self.f.domain())(self.f.domain()(x))
        #except (AttributeError, TypeError):
        #    return FiniteField_prime_modn.__call__(self, x)

class ResidueFiniteField_ext_pari(ResidueField_generic, FiniteField_ext_pari):
    """
    The class representing residue fields of number fields that have non-prime order >= 2^16.

    EXAMPLES:
        sage: K.<a> = NumberField(x^3-7)
        sage: P = K.ideal(923478923).factor()[0][0]
        sage: k = K.residue_field(P)
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
            sage: P = K.factor(10007)[0][0]
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
            return FiniteField_ext_pari.__call__(self, x)
        except TypeError:
            return self._structure[0](x)
        #try:
        #    return self.coerce_map_from(self.f.domain())(self.f.domain()(x))
        #except (AttributeError, TypeError):
        #    return FiniteField_ext_pari.__call__(self, x)

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
            sage: k = K.residue_field(P)
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
            return FiniteField_givaro.__call__(self, x)
        except TypeError:
            try:
                return self._structure[0](x)
            except:
                raise TypeError



