r"""
Finite Fields of characteristic 2 and order > 2.

This implementation uses NTL's GF2E class to perform the arithmetic
and is the standard implementation for $GF(2^n)$ for $n >= 16$.

AUTHORS:
     -- Martin Albrecht <malb@informatik.uni-bremen.de> (2007-10)
"""


#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
include "../libs/ntl/decl.pxi"
include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"

from sage.structure.sage_object cimport SageObject

from sage.structure.parent  cimport Parent
from sage.structure.parent_base cimport ParentWithBase
from sage.structure.parent_gens cimport ParentWithGens
from sage.structure.element cimport Element, ModuleElement

from sage.rings.ring cimport Ring
from sage.structure.element cimport RingElement

from sage.rings.ring cimport FiniteField
from sage.structure.element cimport FiniteFieldElement

from sage.libs.pari.all import pari
from sage.libs.pari.gen import gen

from sage.interfaces.gap import is_GapElement

from sage.misc.randstate import current_randstate

from finite_field_ext_pari import FiniteField_ext_pari
from finite_field_element import FiniteField_ext_pariElement

from polynomial.polynomial_ring_constructor import PolynomialRing

cdef object is_IntegerMod
cdef object IntegerModRing_generic
cdef object Integer
cdef object Rational
cdef object is_Polynomial
cdef object ConwayPolynomials
cdef object conway_polynomial
cdef object MPolynomial
cdef object Polynomial
cdef object FreeModuleElement
cdef object GF
cdef object GF2, GF2_0, GF2_1

cdef void late_import():
    """
    Import a bunch of objects to the module name space late,
    i.e. after the module was loaded. This is needed to avoid circular
    imports.
    """
    global is_IntegerMod, \
           IntegerModRing_generic, \
           Integer, \
           Rational, \
           is_Polynomial, \
           ConwayPolynomials, \
           conway_polynomial, \
           MPolynomial, \
           Polynomial, \
           FreeModuleElement, \
           GF, \
           GF2, GF2_0, GF2_1

    if is_IntegerMod is not None:
        return

    import sage.rings.integer_mod
    is_IntegerMod = sage.rings.integer_mod.is_IntegerMod

    import sage.rings.integer_mod_ring
    IntegerModRing_generic = sage.rings.integer_mod_ring.IntegerModRing_generic

    import sage.rings.integer
    Integer = sage.rings.integer.Integer

    import sage.rings.rational
    Rational = sage.rings.rational.Rational

    import sage.rings.polynomial.polynomial_element
    is_Polynomial = sage.rings.polynomial.polynomial_element.is_Polynomial

    import sage.databases.conway
    ConwayPolynomials = sage.databases.conway.ConwayPolynomials

    import sage.rings.finite_field
    conway_polynomial = sage.rings.finite_field.conway_polynomial

    import sage.rings.polynomial.multi_polynomial_element
    MPolynomial = sage.rings.polynomial.multi_polynomial_element.MPolynomial

    import sage.rings.polynomial.polynomial_element
    Polynomial = sage.rings.polynomial.polynomial_element.Polynomial

    import sage.modules.free_module_element
    FreeModuleElement = sage.modules.free_module_element.FreeModuleElement

    import sage.rings.finite_field
    GF = sage.rings.finite_field.FiniteField
    GF2 = GF(2)
    GF2_0 = GF2(0)
    GF2_1 = GF2(1)

cdef extern from "arpa/inet.h":
    unsigned int htonl(unsigned int)

cdef little_endian():
    return htonl(1) != 1

cdef unsigned int switch_endianess(unsigned int i):
    cdef int j
    cdef unsigned int ret = 0
    for j from 0 <= j < sizeof(int):
        (<unsigned char*>&ret)[j] = (<unsigned char*>&i)[sizeof(int)-j-1]
    return ret

cdef class FiniteField_ntl_gf2e(FiniteField):
    def __init__(FiniteField_ntl_gf2e self, q, names="a",  modulus=None, repr="poly"):
        """
        Finite Field for characteristic 2 and order >= 2.

        INPUT:
            q     -- 2^n (must be 2 power)
            names  -- variable used for poly_repr (default: 'a')
            modulus -- you may provide a minimal polynomial to use for
                     reduction or 'random' to force a random
                     irreducible polynomial. (default: None, a conway
                     polynomial is used if found. Otherwise a random sparse
                     polynomial is used)
            repr  -- controls the way elements are printed to the user:
                     (default: 'poly')
                     'poly': polynomial representation

        OUTPUT:
            Finite field with characteristic 2 and cardinality 2^n.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: type(k)
            <type 'sage.rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e'>
            sage: k.<a> = GF(2^1024)
            sage: k.modulus()
            x^1024 + x^19 + x^6 + x + 1
            sage: set_random_seed(0)
            sage: k.<a> = GF(2^17, modulus='random')
            sage: k.modulus()
            x^17 + x^16 + x^15 + x^10 + x^8 + x^6 + x^4 + x^3 + x^2 + x + 1
            sage: k.modulus().is_irreducible()
            True
        """
        self._zero_element = self._new()
        GF2E_conv_long((<FiniteField_ntl_gf2eElement>self._zero_element).x,0)

        self._one_element = self._new()
        GF2E_conv_long((<FiniteField_ntl_gf2eElement>self._one_element).x,1)

    def __new__(FiniteField_ntl_gf2e self, q, names="a",  modulus=None, repr="poly"):
        cdef GF2X_c ntl_m
        cdef GF2_c c
        cdef GF2X_c ntl_tmp

        # we are calling late_import here because this constructor is
        # called at least once before any arithmetic is performed.
        late_import()

        q = Integer(q)
        if q < 2:
            raise ValueError, "q  must be a 2-power"
        F = q.factor()
        if len(F) > 1:
            raise ValueError, "q must be a 2-power"
        p = F[0][0]
        k = F[0][1]

        if p != 2:
            raise ValueError, "q must be a 2-power"

        ParentWithGens.__init__(self, GF(p), names, normalize=True)

        self._kwargs = {'repr':repr}
        self._is_conway = False

        if modulus is None or modulus == "random":
            if modulus == "random":
                current_randstate().set_seed_ntl(False)
                GF2X_BuildSparseIrred(ntl_tmp, k)
                GF2X_BuildRandomIrred(ntl_m, ntl_tmp)
            elif ConwayPolynomials().has_polynomial(p, k):
                modulus = conway_polynomial(p, k)
                modulus = modulus.list()
                self._is_conway = True
                for i from 0 <= i < len(modulus):
                    GF2_conv_long(c, int(modulus[i]))
                    GF2X_SetCoeff(ntl_m, i, c)
            else:
                GF2X_BuildSparseIrred(ntl_m, k)
            GF2EContext_construct_GF2X(&self.F, &ntl_m)
        else:
            if is_Polynomial(modulus):
                modulus = modulus.list()
            if PY_TYPE_CHECK(modulus, list) or PY_TYPE_CHECK(modulus, tuple):
                for i from 0 <= i < len(modulus):
                    GF2_conv_long(c, int(modulus[i]))
                    GF2X_SetCoeff(ntl_m, i, c)
                GF2EContext_construct_GF2X(&self.F, &ntl_m)
            else:
                raise TypeError, "Modulus parameter not understood"

    def __dealloc__(FiniteField_ntl_gf2e self):
        GF2EContext_destruct(&self.F)

    def __doctest_for_5340(self):
        r"""
        Every bug fix should have a doctest.  But #5340 only happens when
        a garbage collection happens between restoring the modulus and
        using it, so it can't be reliably doctested using any of the
        existing Cython functions in this module.  The sole purpose of
        this method is to doctest the fix for #5340.

        EXAMPLES:
            sage: k.<a> = GF(2^20)
            sage: k.__doctest_for_5340()
            [1 1 0 0 1 1 1 1 0 1 1 0 0 0 0 0 0 0 0 0 1]
            [1 1 0 0 1 1 1 1 0 1 1 0 0 0 0 0 0 0 0 0 1]
        """
        import gc

        # Do a garbage collection, so that this method is repeatable.
        gc.collect()

        # Create a new finite field.
        new_fld = GF(1<<30, 'a', modulus='random')
        # Create a garbage cycle.
        cycle = [new_fld]
        cycle.append(cycle)
        # Make the cycle inaccessible.
        new_fld = None
        cycle = None

        # Set the modulus.
        self.F.restore()
        # Print the current modulus.
        cdef GF2XModulus_c modulus = GF2E_modulus()
        cdef GF2X_c mod_poly = GF2XModulus_GF2X(modulus)
        print GF2X_to_PyString(&mod_poly)

        # do another garbage collection
        gc.collect()

        # and print the modulus again
        modulus = GF2E_modulus()
        mod_poly = GF2XModulus_GF2X(modulus)
        print GF2X_to_PyString(&mod_poly)

    cdef FiniteField_ntl_gf2eElement _new(FiniteField_ntl_gf2e self):
        """
        Return a new element in self. Use this method to construct
        'empty' elements.
        """
        cdef FiniteField_ntl_gf2eElement y
        self.F.restore()
        y = PY_NEW(FiniteField_ntl_gf2eElement)
        y._parent = <ParentWithBase>self
        return y

    def characteristic(FiniteField_ntl_gf2e self):
        """
        Return 2.

        EXAMPLE:
            sage: k.<a> = GF(2^16,modulus='random')
            sage: k.characteristic()
            2
        """
        return Integer(2)

    def order(self):
        """
        Return the cardinality of this field.

        EXAMPLE:
            sage: k.<a> = GF(2^64)
            sage: k.order()
            18446744073709551616
        """
        self.F.restore()
        return Integer(1) << GF2E_degree()

    def degree(FiniteField_ntl_gf2e self):
        r"""
        If \code{self.cardinality() == p^n} this method returns $n$.

        EXAMPLE:
            sage: k.<a> = GF(2^64)
            sage: k.degree()
            64
        """
        self.F.restore()
        return Integer(GF2E_degree())

    def is_atomic_repr(FiniteField_ntl_gf2e self):
        """
        Return whether elements of self are printed using an atomic
        representation.

        EXAMPLE:
            sage: k.<a> = GF(2^64)
            sage: k.is_atomic_repr()
            False
            sage: P.<x> = PolynomialRing(k)
            sage: (a+1)*x # indirect doctest
            (a + 1)*x
        """
        return False

    def __call__(FiniteField_ntl_gf2e self, e):
        """
        Coerces several data types to self.

        INPUT:
            e -- data to coerce

        EXAMPLE:
            sage: k.<a> = GF(2^20)
            sage: k(1)
            1
            sage: k(int(2))
            0

            sage: k('a+1')
            a + 1
            sage: k('b+1')
            Traceback (most recent call last):
            ...
            NameError: name 'b' is not defined

	    sage: R.<x>=GF(2)[]
	    sage: k(1+x+x^10+x^55)
	    a^19 + a^17 + a^16 + a^15 + a^12 + a^11 + a^8 + a^6 + a^4 + a^2 + 1

            sage: V = k.vector_space()
            sage: v = V.random_element(); v
            (1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1)
            sage: k(v)
            a^19 + a^15 + a^14 + a^13 + a^11 + a^10 + a^9 + a^6 + a^5 + a^4 + 1
            sage: vector(k(v)) == v
            True

            sage: k(pari('Mod(1,2)*a^20'))
            a^10 + a^9 + a^7 + a^6 + a^5 + a^4 + a + 1
        """
        cdef FiniteField_ntl_gf2eElement res = self._new()
        cdef FiniteField_ntl_gf2eElement x
        cdef FiniteField_ntl_gf2eElement g

        if PY_TYPE_CHECK(e, FiniteField_ntl_gf2eElement):
            if e.parent() is self:
                return e # this is safe because elements are immutable
            if e.parent() == self:
                res.x = (<FiniteField_ntl_gf2eElement>e).x
                return res
            if e.parent() is GF2 or e.parent() == GF2:
                GF2E_conv_long(res.x,int(e))
                return res

        elif PY_TYPE_CHECK(e, int) or \
             PY_TYPE_CHECK(e, Integer) or \
             PY_TYPE_CHECK(e, long) or is_IntegerMod(e):
            GF2E_conv_long(res.x,int(e))
            return res

        elif PY_TYPE_CHECK(e, float):
            GF2E_conv_long(res.x,int(e))
            return res

        elif PY_TYPE_CHECK(e, str):
            return self(eval(e.replace("^","**"),self.gens_dict()))

        elif PY_TYPE_CHECK(e, FreeModuleElement):
            if self.vector_space() != e.parent():
                raise TypeError, "e.parent must match self.vector_space"
            ret = self._zero_element
            for i in range(len(e)):
                ret = ret + self(int(e[i]))*self.gen()**i
            return ret

        elif PY_TYPE_CHECK(e, MPolynomial):
            if e.is_constant():
                return self(e.constant_coefficient())
            else:
                raise TypeError, "no coercion defined"

        elif PY_TYPE_CHECK(e, Polynomial):
            if e.is_constant():
                return self(e.constant_coefficient())
            else:
                return e(self.gen())

        elif PY_TYPE_CHECK(e, Rational):
            num = e.numer()
            den = e.denom()
            return self(num)/self(den)

        elif PY_TYPE_CHECK(e, gen):
            pass # handle this in next if clause

        elif PY_TYPE_CHECK(e,FiniteField_ext_pariElement):
            # reduce FiniteFieldElements to pari
            e = e._pari_()

        elif is_GapElement(e):
            from sage.interfaces.gap import gfq_gap_to_sage
            return gfq_gap_to_sage(e, self)
        else:
            raise TypeError, "unable to coerce"

        if PY_TYPE_CHECK(e, gen):
            e = e.lift().lift()
            try:
                GF2E_conv_long(res.x, int(e[0]))
            except TypeError:
                GF2E_conv_long(res.x, int(e))

            g = self._new()
            GF2E_from_str(&g.x, "[0 1]") # not the fastest
            x = self._new()
            GF2E_conv_long(x.x,1)

            for i from 0 < i <= e.poldegree():
                GF2E_mul(x.x, x.x, g.x)
                if e[i]:
                    GF2E_add(res.x, res.x, x.x )
            return res

        raise ValueError, "Cannot coerce element %s to self."%(e)

    cdef _coerce_c_impl(self, e):
        """
        Coercion accepts elements of self.parent(), ints, and prime subfield elements.
        """
        cdef FiniteField_ntl_gf2eElement res = self._new()

        if PY_TYPE_CHECK(e, int) \
               or PY_TYPE_CHECK(e, long) or PY_TYPE_CHECK(e, Integer):
            GF2E_conv_long(res.x,int(e))
            return res

        if PY_TYPE_CHECK(e, FiniteFieldElement) or \
               PY_TYPE_CHECK(e, FiniteField_ntl_gf2eElement) or is_IntegerMod(e):
            K = e.parent()
            if K is <object>self:
                return e
            if PY_TYPE_CHECK(K, IntegerModRing_generic) \
                   and K.characteristic() % self.characteristic() == 0:
                GF2E_conv_long(res.x,int(e))
                return res
            if K.characteristic() == self.characteristic():
                if K.degree() == 1:
                    GF2E_conv_long(res.x,int(e))
                    return res
                elif self.degree() % K.degree() == 0:
                    # This is where we *would* do coercion from one nontrivial finite field to another...
                    raise TypeError, 'no canonical coercion defined'

        raise TypeError, 'no canonical coercion defined'

    def gen(FiniteField_ntl_gf2e self, ignored=None):
        r"""
        Return a generator of self.

        EXAMPLE:
            sage: k.<a> = GF(2^19)
            sage: k.gen() == a
            True
            sage: a
            a
        """
        cdef FiniteField_ntl_gf2eElement x = self._new()
        GF2E_from_str(&x.x,"[0 1]")
        return x

    def prime_subfield(FiniteField_ntl_gf2e self):
        r"""
        Return the prime subfield $\FF_p$ of self if self is $\FF_{p^n}$.

        EXAMPLE:
            sage: F.<a> = GF(2^16)
            sage: F.prime_subfield()
            Finite Field of size 2
        """
        return GF2

    def fetch_int(FiniteField_ntl_gf2e self, number):
        r"""
        Given an integer $n$ return a finite field element in self
        which equals $n$ if self.gen() was set to
        self.characteristic().

        INPUT:
            number -- an integer

        EXAMPLES:
            sage: k.<a> = GF(2^48)
            sage: k.fetch_int(2^43 + 2^15 + 1)
            a^43 + a^15 + 1
            sage: k.fetch_int(33793)
            a^15 + a^10 + 1
            sage: 33793.digits(2) # little endian
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        """
        cdef FiniteField_ntl_gf2eElement a = self._new()
        cdef GF2X_c _a

        self.F.restore()

        cdef unsigned char *p
        cdef int i

        if number < 0 or number >= GF2E_degree():
            TypeError, "n must be between 0 and self.order()"

        if PY_TYPE_CHECK(number, int) or PY_TYPE_CHECK(number, long):
            from sage.misc.functional import log
            n = int(log(number,2))/8 + 1
        elif PY_TYPE_CHECK(number, Integer):
            n = int(number.nbits())/8 + 1
        else:
            raise TypeError, "number %s is not an integer"%number

        p = <unsigned char*>sage_malloc(n)
        for i from 0 <= i < n:
            p[i] = (number%256)
            number = number >> 8
        GF2XFromBytes(_a, p, n)
        GF2E_conv_GF2X(a.x, _a)
        sage_free(p)
        return a

    def polynomial(self, name = None):
        """
        Return the defining polynomial of this field as an element of
        self.polynomial_ring().

        This is the same as the characteristic polynomial of the
        generator of self.

        INPUT:
            name -- optional variable name

        EXAMPLE:
            sage: k.<a> = GF(2^20)
            sage: k.polynomial()
            a^20 + a^10 + a^9 + a^7 + a^6 + a^5 + a^4 + a + 1
            sage: k.polynomial('FOO')
            FOO^20 + FOO^10 + FOO^9 + FOO^7 + FOO^6 + FOO^5 + FOO^4 + FOO + 1
            sage: a^20
            a^10 + a^9 + a^7 + a^6 + a^5 + a^4 + a + 1

        """
        cdef FiniteField_ntl_gf2eElement P
        cdef GF2X_c _P
        cdef GF2_c c
        self.F.restore()
        cdef int i

        if self._polynomial is None or name is not None:
            P = -(self.gen()**(self.degree()))
            _P = GF2E_rep(P.x)

            ret = []
            for i from 0 <= i < GF2E_degree():
                c = GF2X_coeff(_P,i)
                if not GF2_IsZero(c):
                    ret.append(1)
                else:
                    ret.append(0)
            ret.append(1)

            R = self.polynomial_ring(name)
            if name is None:
                self._polynomial = R(ret)
                return self._polynomial
            else:
                return R(ret)
        else:
            return self._polynomial

    def _finite_field_ext_pari_(self):
        """
        Return a FiniteField_ext_pari isomorphic to self with the same
        defining polynomial.

        This method will vanish eventually because that implementation of
        finite fields will be deprecated.

        EXAMPLE:
            sage: k.<a> = GF(2^20)
            sage: kP = k._finite_field_ext_pari_()
            sage: kP
            Finite Field in a of size 2^20
            sage: type(kP)
            <class 'sage.rings.finite_field_ext_pari.FiniteField_ext_pari'>
        """
        f = self.polynomial()
        return FiniteField_ext_pari(self.order(), self.variable_name(), f)

    def _finite_field_ext_pari_modulus_as_str(self):
        """
        Return a PARI string representation of the modulus.

        This method might vanish eventually, because PARI is getting
        replaced by NTL as main implementation for finite extension
        fields.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: k.polynomial()
            a^16 + a^5 + a^3 + a^2 + 1
            sage: k._finite_field_ext_pari_modulus_as_str()
            'Mod(1, 2)*a^16 + Mod(1, 2)*a^5 + Mod(1, 2)*a^3 + Mod(1, 2)*a^2 + Mod(1, 2)'
        """
        return self._finite_field_ext_pari_().modulus()._pari_init_()

    def __richcmp__(left, right, int op):
        """
            sage: k1.<a> = GF(2^16)
            sage: k2.<a> = GF(2^17)
            sage: k1 == k2
            False
            sage: k3 = k1._finite_field_ext_pari_()
            sage: k1 == k3
            True
        """
        return (<Parent>left)._richcmp(right, op)

    def __hash__(FiniteField_ntl_gf2e self):
        """
            sage: k1.<a> = GF(2^16)
            sage: {k1:1} # indirect doctest
            {Finite Field in a of size 2^16: 1}
        """
        if self._hash is None:
            self._hash = hash((self.characteristic(),self.polynomial(),self.variable_name(),"ntl_gf2e"))
        return self._hash

    def _pari_modulus(self):
        """
        Return PARI object which is equivalent to the
        polynomial/modulus of self.

        EXAMPLE:
            sage: k1.<a> = GF(2^16)
            sage: k1._pari_modulus()
            Mod(1, 2)*a^16 + Mod(1, 2)*a^5 + Mod(1, 2)*a^3 + Mod(1, 2)*a^2 + Mod(1, 2)
        """
        f = pari(str(self.modulus()))
        return f.subst('x', 'a') * pari("Mod(1,%s)"%self.characteristic())

cdef class FiniteField_ntl_gf2eElement(FiniteFieldElement):
    """
    An element of an NTL:GF2E finite field.
    """

    def __init__(FiniteField_ntl_gf2eElement self, parent=None ):
        """
        Initializes an element in parent. It's much better to use
        parent(<value>) or any specialized method of parent
        (one,zero,gen) instead. Do not call this constructor directly,
        it doesn't make much sense.

        INPUT:
            parent -- base field

        OUTPUT:
            finite field element.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: from sage.rings.finite_field_ntl_gf2e import FiniteField_ntl_gf2eElement
            sage: FiniteField_ntl_gf2eElement(k)
            0
        """
        if parent is None:
            raise ValueError, "You must provide a parent to construct a finite field element"

    def __new__(FiniteField_ntl_gf2eElement self, parent=None ):
        if parent is None:
            return
        if PY_TYPE_CHECK(parent, FiniteField_ntl_gf2e):
            self._parent = parent
            (<FiniteField_ntl_gf2e>self._parent).F.restore()
            GF2E_construct(&self.x)

    def __dealloc__(FiniteField_ntl_gf2eElement self):
        GF2E_destruct(&self.x)

    cdef FiniteField_ntl_gf2eElement _new(FiniteField_ntl_gf2eElement self):
        cdef FiniteField_ntl_gf2eElement y
        (<FiniteField_ntl_gf2e>self._parent).F.restore()
        y = PY_NEW(FiniteField_ntl_gf2eElement)
        y._parent = self._parent
        return y

    def __repr__(FiniteField_ntl_gf2eElement self):
        """
        Polynomial representation of self.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: str(a^16) # indirect doctest
            'a^5 + a^3 + a^2 + 1'
            sage: k.<u> = GF(2^16)
            sage: u
            u
        """
        (<FiniteField_ntl_gf2e>self._parent).F.restore()
        cdef GF2X_c rep = GF2E_rep(self.x)
        cdef GF2_c c
        cdef int i

        if GF2E_IsZero(self.x):
            return "0"

        name = self._parent.variable_name()
        _repr = []

        c = GF2X_coeff(rep, 0)
        if not GF2_IsZero(c):
            _repr.append("1")

        c = GF2X_coeff(rep, 1)
        if not GF2_IsZero(c):
            _repr.append(name)

        for i from 1 < i <= GF2X_deg(rep):
            c = GF2X_coeff(rep, i)
            if not GF2_IsZero(c):
                _repr.append("%s^%d"%(name,i))

        return " + ".join(reversed(_repr))

    def parent(self):
        """
        Return parent finite field.

        EXAMPLE:
            sage: k.<a> = GF(2^20)
            sage: a.parent() is k
            True
        """
        return (<FiniteField_ntl_gf2e>self._parent)

    def __nonzero__(FiniteField_ntl_gf2eElement self):
        r"""
        Return True if \code{self != k(0)}.

        EXAMPLE:
            sage: k.<a> = GF(2^20)
            sage: bool(a) # indirect doctest
            True
            sage: bool(k(0))
            False
            sage: a.is_zero()
            False
        """
        (<FiniteField_ntl_gf2e>self._parent).F.restore()
        return not bool(GF2E_IsZero(self.x))

    def is_one(FiniteField_ntl_gf2eElement self):
        r"""
        Return True if \code{self == k(1)}.

        Return True if \code{self != k(0)}.

        EXAMPLE:
            sage: k.<a> = GF(2^20)
            sage: a.is_one() # indirect doctest
            False
            sage: k(1).is_one()
            True

        """
        (<FiniteField_ntl_gf2e>self._parent).F.restore()
        return bool(GF2E_equal(self.x,(<FiniteField_ntl_gf2eElement>(<FiniteField_ntl_gf2e>self._parent)._one_element).x))

    def is_unit(FiniteField_ntl_gf2eElement self):
        """
        Return True if self is nonzero, so it is a unit as an element
        of the finite field.

        EXAMPLE:
            sage: k.<a> = GF(2^17)
            sage: a.is_unit()
            True
            sage: k(0).is_unit()
            False
        """
        (<FiniteField_ntl_gf2e>self._parent).F.restore()
        if not GF2E_IsZero(self.x):
            return True
        else:
            return False

    def is_square(FiniteField_ntl_gf2eElement self):
        """
        Return True as every element in GF(2^n) is a square.

        EXAMPLE:
            sage: k.<a> = GF(2^18)
            sage: e = k.random_element()
            sage: e
            a^15 + a^14 + a^13 + a^11 + a^10 + a^9 + a^6 + a^5 + a^4 + 1
            sage: e.is_square()
            True
            sage: e.sqrt()
            a^16 + a^15 + a^14 + a^11 + a^9 + a^8 + a^7 + a^6 + a^4 + a^3 + 1
            sage: e.sqrt()^2 == e
            True
        """
        return True

    def sqrt(FiniteField_ntl_gf2eElement self, all=False, extend=False):
        """
        Return a square root of this finite field element in its
        parent, if there is one.  Otherwise, raise a ValueError.

        EXAMPLE:
            sage: k.<a> = GF(2^20)
            sage: a.is_square()
            True
            sage: a.sqrt()
            a^19 + a^15 + a^14 + a^12 + a^9 + a^7 + a^4 + a^3 + a + 1
            sage: a.sqrt()^2 == a
            True

        This failed before \#4899:
            sage: GF(2^16,'a')(1).sqrt()
            1

        """
        # this really should be handled special, its gf2 linear after
        # all
        if all:
            a = self.sqrt()
            return [a]
        cdef FiniteField_ntl_gf2eElement e
        if self.is_one():
            return self
        else:
            return self.nth_root(2)

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Add two elements.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: e = a^2 + a + 1
            sage: f = a^15 + a^2 + 1
            sage: e + f
            a^15 + a
        """
        cdef FiniteField_ntl_gf2eElement r = (<FiniteField_ntl_gf2eElement>self)._new()
        GF2E_add(r.x, (<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>right).x)
        return r

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        """
        Add two elements.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: e = a^2 + a + 1
            sage: f = a^15 + a^2 + 1
            sage: e + f
            a^15 + a
        """
        GF2E_add((<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>right).x)
        return self

    cpdef RingElement _mul_(self, RingElement right):
        """
        Multiply two elements.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: e = a^2 + a + 1
            sage: f = a^15 + a^2 + 1
            sage: e * f
            a^15 + a^6 + a^5 + a^3 + a^2
        """
        cdef FiniteField_ntl_gf2eElement r = (<FiniteField_ntl_gf2eElement>self)._new()
        GF2E_mul(r.x, (<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>right).x)
        return r

    cpdef RingElement _imul_(self, RingElement right):
        """
        Multiply two elements.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: e = a^2 * a + 1
            sage: f = a^15 * a^2 + 1
            sage: e * f
            a^9 + a^7 + a + 1
        """
        GF2E_mul((<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>right).x)
        return self

    cpdef RingElement _div_(self, RingElement right):
        """
        Divide two elements.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: e = a^2 + a + 1
            sage: f = a^15 + a^2 + 1
            sage: e / f
            a^11 + a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + 1
            sage: k(1) / k(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in finite field.
        """
        cdef FiniteField_ntl_gf2eElement r = (<FiniteField_ntl_gf2eElement>self)._new()
        if GF2E_IsZero((<FiniteField_ntl_gf2eElement>right).x):
            raise ZeroDivisionError, 'division by zero in finite field.'
        GF2E_div(r.x, (<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>right).x)
        return r

    cpdef RingElement _idiv_(self, RingElement right):
        """
        Divide two elements.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: e = a^2 / a + 1
            sage: f = a^15 / a^2 + 1
            sage: e / f
            a^15 + a^12 + a^10 + a^9 + a^6 + a^5 + a^3
        """
        GF2E_div((<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>right).x)
        return self

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Subtract two elements.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: e = a^2 - a + 1
            sage: f = a^15 - a^2 + 1
            sage: e - f
            a^15 + a
        """
        cdef FiniteField_ntl_gf2eElement r = (<FiniteField_ntl_gf2eElement>self)._new()
        GF2E_sub(r.x, (<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>right).x)
        return r

    cpdef ModuleElement _isub_(self, ModuleElement right):
        """
        Subtract two elements.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: e = a^2 - a + 1
            sage: f = a^15 - a^2 + 1
            sage: e - f
            a^15 + a
        """
        GF2E_sub((<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>right).x)
        return self

    def __neg__(FiniteField_ntl_gf2eElement self):
        """
        Return this element.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: -a
            a
        """
        cdef FiniteField_ntl_gf2eElement r = (<FiniteField_ntl_gf2eElement>self)._new()
        r.x = (<FiniteField_ntl_gf2eElement>self).x
        return r

    def __invert__(FiniteField_ntl_gf2eElement self):
        """
        Return the multiplicative inverse of an element.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: ~a
            a^15 + a^4 + a^2 + a
            sage: a * ~a
            1
        """
        cdef FiniteField_ntl_gf2eElement r = (<FiniteField_ntl_gf2eElement>self)._new()
        cdef FiniteField_ntl_gf2eElement o = (<FiniteField_ntl_gf2eElement>self)._parent._one_element
        GF2E_div(r.x, o.x, (<FiniteField_ntl_gf2eElement>self).x)
        return r

    def __pow__(FiniteField_ntl_gf2eElement self, exp, other):
        """
            sage: k.<a> = GF(2^63)
            sage: a^2
            a^2
            sage: a^64
            a^25 + a^24 + a^23 + a^18 + a^17 + a^16 + a^12 + a^10 + a^9 + a^5 + a^4 + a^3 + a^2 + a
            sage: a^(2^64)
            a^2
            sage: a^(2^128)
            a^4
        """
        cdef int exp_int
        cdef FiniteField_ntl_gf2eElement r = (<FiniteField_ntl_gf2eElement>self)._new()

        try:
            exp_int = exp
            if exp != exp_int:
                raise OverflowError
            GF2E_power(r.x, (<FiniteField_ntl_gf2eElement>self).x, exp_int)
            return r
        except OverflowError:
            # we could try to factor out the order first
            from sage.groups.generic import power
            return power(self,exp)

    def __richcmp__(left, right, int op):
        """
        Comparison of finite field elements.

        EXAMPLE:
            sage: k.<a> = GF(2^20)
            sage: e = k.random_element()
            sage: f = loads(dumps(e))
            sage: e is f
            False
            sage: e == f
            True
            sage: e != (e + 1)
            True

        NOTE: that in finite fields $<$ and $>$ don't make sense and
        that the result of these operators has no mathematical meaning
        and may vary across different finite field implementations.

        EXAMPLE:
        """
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Comparison of finite field elements.
        """
        (<FiniteField_ntl_gf2e>left._parent).F.restore()
        c = GF2E_equal((<FiniteField_ntl_gf2eElement>left).x, (<FiniteField_ntl_gf2eElement>right).x)
        if c == 1:
            return 0
        else:
            return 1

    def __int__(FiniteField_ntl_gf2eElement self):
        """
        Return the int representation of self.  When self is in the
        prime subfield, the integer returned is equal to self and not
        to \code{log_repr}.

        Elements of this field are represented as ints in as follows:
        for $e \in \FF_p[x]$ with $e = a_0 + a_1x + a_2x^2 + \cdots $, $e$ is
        represented as: $n= a_0 + a_1  p + a_2  p^2 + \cdots$.

        EXAMPLE:
            sage: k.<a> = GF(2^20)
            sage: int(a)
            2
            sage: int(a^2 + 1)
            5
            sage: k.<a> = GF(2^70)
            sage: int(a^65 + a^64 + 1)
            55340232221128654849L
        """
        cdef unsigned int i = 0
        ret = int(0)
        cdef unsigned long shift = 0
        cdef GF2X_c r = GF2E_rep(self.x)

        if GF2X_IsZero(r):
            return 0

        if little_endian():
            while GF2X_deg(r) >= sizeof(int)*8:
                BytesFromGF2X(<unsigned char *>&i, r, sizeof(int))
                ret += int(i) << shift
                shift += sizeof(int)*8
                GF2X_RightShift(r,r,(sizeof(int)*8))
            BytesFromGF2X(<unsigned char *>&i, r, sizeof(int))
            ret += int(i) << shift
        else:
            while GF2X_deg(r) >= sizeof(int)*8:
                BytesFromGF2X(<unsigned char *>&i, r, sizeof(int))
                ret += int(switch_endianess(i)) << shift
                shift += sizeof(int)*8
                GF2X_RightShift(r,r,(sizeof(int)*8))
            BytesFromGF2X(<unsigned char *>&i, r, sizeof(int))
            ret += int(switch_endianess(i)) << shift

        return int(ret)

    def polynomial(FiniteField_ntl_gf2eElement self, name=None):
        r"""
        Return self viewed as a polynomial over
        \code{self.parent().prime_subfield()}.

        INPUT:
            name -- (optional) variable name

        EXAMPLE:
            sage: k.<a> = GF(2^17)
            sage: e = a^15 + a^13 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^4 + a + 1
            sage: e.polynomial()
            a^15 + a^13 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^4 + a + 1

            sage: from sage.rings.polynomial.polynomial_element import is_Polynomial
            sage: is_Polynomial(e.polynomial())
            True

            sage: e.polynomial('x')
            x^15 + x^13 + x^11 + x^10 + x^9 + x^8 + x^7 + x^6 + x^4 + x + 1
        """
        cdef GF2X_c r = GF2E_rep(self.x)
        cdef int i

        C = []
        for i from 0 <= i <= GF2X_deg(r):
            C.append(GF2_conv_to_long(GF2X_coeff(r,i)))
        return self._parent.polynomial_ring(name)(C)

    def charpoly(self, var='x'):
        r"""
        Return the characteristic polynomial of self as a polynomial
        in var over the prime subfield.

        INPUT:
            var -- string (default: 'x')
        OUTPUT:
            polynomial

        EXAMPLES:
            sage: k.<a> = GF(2^8)
            sage: b = a^3 + a
            sage: b.minpoly()
            x^4 + x^3 + x^2 + x + 1
            sage: b.charpoly()
            x^8 + x^6 + x^4 + x^2 + 1
            sage: b.charpoly().factor()
            (x^4 + x^3 + x^2 + x + 1)^2
            sage: b.charpoly('Z')
            Z^8 + Z^6 + Z^4 + Z^2 + 1
        """
        f = self.minpoly(var)
        cdef int d = f.degree(), n = self.parent().degree()
        cdef int pow = n/d
        return f if pow == 1 else f**pow


    def minpoly(self, var='x'):
        r"""
        Return the minimal polynomial of self, which is the smallest
        degree polynomial $f \in \mathbf{F}_{2}[x]$ such that
        $f(self) = 0$.

        INPUT:
            var -- string (default: 'x')
        OUTPUT:
            polynomial

        EXAMPLES:
            sage: K.<a> = GF(2^100)
            sage: f = a.minpoly(); f
            x^100 + x^57 + x^56 + x^55 + x^52 + x^48 + x^47 + x^46 + x^45 + x^44 + x^43 + x^41 + x^37 + x^36 + x^35 + x^34 + x^31 + x^30 + x^27 + x^25 + x^24 + x^22 + x^20 + x^19 + x^16 + x^15 + x^11 + x^9 + x^8 + x^6 + x^5 + x^3 + 1
            sage: f(a)
            0
            sage: g = K.random_element()
            sage: g.minpoly()(g)
            0
        """
        cdef GF2X_c r = GF2X_IrredPolyMod(GF2E_rep(self.x), GF2E_modulus())
        cdef int i
        C = []
        for i from 0 <= i <= GF2X_deg(r):
            C.append(GF2_conv_to_long(GF2X_coeff(r,i)))
        return self._parent.polynomial_ring(var)(C)

    def trace(self):
        """
        Return the trace of self.

        EXAMPLES:
            sage: K.<a> = GF(2^25)
            sage: a.trace()
            0
            sage: a.charpoly()
            x^25 + x^8 + x^6 + x^2 + 1
            sage: parent(a.trace())
            Finite Field of size 2

            sage: b = a+1
            sage: b.trace()
            1
            sage: b.charpoly()[1]
            1
        """
        if GF2_IsOne(GF2E_trace(self.x)):
            return GF2_1
        else:
            return GF2_0

    def weight(self):
        """
        Returns the number of non-zero coefficients in the polynomial
        representation of self.

        EXAMPLES:
            sage: K.<a> = GF(2^21)
            sage: a.weight()
            1
            sage: (a^5+a^2+1).weight()
            3
            sage: b = 1/(a+1); b
            a^20 + a^19 + a^18 + a^17 + a^16 + a^15 + a^14 + a^13 + a^12 + a^11 + a^10 + a^9 + a^8 + a^7 + a^6 + a^4 + a^3 + a^2
            sage: b.weight()
            18
        """
        return GF2X_weight(GF2E_rep(self.x))

    def _finite_field_ext_pari_element(FiniteField_ntl_gf2eElement self, k=None):
        r"""
        Return an element of \var{k} supposed to match this
        element. No checks if \var{k} equals \code{self.parent()} are
        performed.

        INPUT:
            k -- (optional) FiniteField_ext_pari

        OUTPUT:
            equivalent of self in k

        EXAMPLE:
            sage: k.<a> = GF(2^17)
            sage: a._finite_field_ext_pari_element()
            a
        """
        if k is None:
            k = (<FiniteField_ntl_gf2e>self._parent)._finite_field_ext_pari_()
        elif not PY_TYPE_CHECK(k, FiniteField_ext_pari):
            raise TypeError, "k must be a pari finite field."
        cdef GF2X_c r = GF2E_rep(self.x)
        cdef int i

        g = k.gen()
        o = k(1)
        ret = k(0)
        for i from 0 <= i <= GF2X_deg(r):
            ret += GF2_conv_to_long(GF2X_coeff(r,i))*o
            o *= g
        return ret

    def _magma_init_(self, magma):
        r"""
        Return a string representation of self that \MAGMA can
        understand.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: a._magma_init_(magma)      # random; optional - magma
            '_sage_[...]'

        NOTE: This method calls \MAGMA to setup the parent.
        """
        km = magma(self.parent())
        vn_m = km.gen(1).name()
        vn_s = str(self.parent().polynomial_ring().gen())
        return str(self.polynomial()).replace(vn_s,vn_m)

    def __copy__(self):
        """
        Return a copy of this element.  Actually just returns self, since
        finite field elements are immutable.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: copy(a) is a
            True
        """
        return self

    def _pari_(self, var=None):
        r"""
        Return a \PARI representation of this element.

        EXAMPLE:
            sage: k.<a> = GF(2^17)
            sage: e = a^3 + a + 1
            sage: e._pari_()
            Mod(a^3 + a + 1, Mod(1, 2)*a^17 + Mod(1, 2)*a^3 + Mod(1, 2))

            sage: e._pari_('w')
            Mod(w^3 + w + 1, Mod(1, 2)*w^17 + Mod(1, 2)*w^3 + Mod(1, 2))
        """
        return pari(self._pari_init_(var))

    def _gap_init_(self):
        r"""
        Return a string that evaluates to the \GAP representation of
        this element.

        A \code{NotImplementedError} is raised if
        \code{self.parent().modulus()} is not a Conway polynomial, as
        the isomorphism of finite fields is not implemented yet.

        EXAMPLE:
            sage: k.<b> = GF(2^16)
            sage: b._gap_init_()
            'Z(65536)^1'
        """
        cdef FiniteField_ntl_gf2e F
        F = self._parent
        if not F._is_conway:
            raise NotImplementedError, "conversion of (NTL) finite field element to GAP not implemented except for fields defined by Conway polynomials."
        if F.order() > 65536:
            raise TypeError, "order (=%s) must be at most 65536."%F.order()
        if self == 0:
            return '0*Z(%s)'%F.order()
        assert F.degree() > 1
        g = F.multiplicative_generator()
        n = self.log(g)
        return 'Z(%s)^%s'%(F.order(), n)

    def __hash__(FiniteField_ntl_gf2eElement self):
        """
        Return the hash of this finite field element.  We hash the parent
        and the underlying integer representation of this element.

        EXAMPLE:
            sage: k.<a> = GF(2^18)
            sage: {a:1,a:0} # indirect doctest
            {a: 0}
        """
        return hash(int(self)) # todo, come up with a faster version

    def _vector_(FiniteField_ntl_gf2eElement self, reverse=False):
        r"""
        Return a vector in \code{self.parent().vector_space()}
        matching \code{self}. The most significant bit is to the
        right.

        INPUT:
            reverse -- reverse the order of the bits
                       from little endian to big endian.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: e = a^14 + a^13 + 1
            sage: vector(e) # little endian
            (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0)

            sage: e._vector_(reverse=True) # big endian
            (0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
        """
        #vector(foo) might pass in ZZ
        if PY_TYPE_CHECK(reverse, Parent):
            raise TypeError, "Base field is fixed to prime subfield."

        cdef GF2X_c r = GF2E_rep(self.x)
        cdef int i

        (<FiniteField_ntl_gf2e>self._parent).F.restore()

        C = []
        for i from 0 <= i < GF2E_degree():
            C.append(GF2_conv_to_long(GF2X_coeff(r,i)))
        if reverse:
            C = list(reversed(C))
        return self._parent.vector_space()(C)

    def __reduce__(FiniteField_ntl_gf2eElement self):
        """
        Used for supporting pickling of finite field elements.

        EXAMPLE:
            sage: k.<a> = GF(2^16)
            sage: loads(dumps(a)) == a
            True
        """
        return unpickleFiniteField_ntl_gf2eElement, (self._parent, str(self))

    def log(self, base):
        """
        Return $x$ such that $b^x = a$, where $x$ is $a$ and $b$
        is the base.

        INPUT:
            self -- finite field element
            b -- finite field element that generates the multiplicative group.

        OUTPUT:
            Integer $x$ such that $a^x = b$, if it exists.
            Raises a ValueError exception if no such $x$ exists.

        EXAMPLES:
            sage: F = GF(17)
            sage: F(3^11).log(F(3))
            11
            sage: F = GF(113)
            sage: F(3^19).log(F(3))
            19
            sage: F = GF(next_prime(10000))
            sage: F(23^997).log(F(23))
            997

            sage: F = FiniteField(2^10, 'a')
            sage: g = F.gen()
            sage: b = g; a = g^37
            sage: a.log(b)
            37
            sage: b^37; a
            a^8 + a^7 + a^4 + a + 1
            a^8 + a^7 + a^4 + a + 1

        AUTHOR: David Joyner and William Stein (2005-11)
        """
        from  sage.groups.generic import discrete_log

        b = self.parent()(base)
        return discrete_log(self, b)

def unpickleFiniteField_ntl_gf2eElement(parent, elem):
    """
    EXAMPLE:
        sage: k.<a> = GF(2^20)
        sage: e = k.random_element()
        sage: f = loads(dumps(e)) # indirect doctest
        sage: e == f
        True
    """
    return parent(elem)
