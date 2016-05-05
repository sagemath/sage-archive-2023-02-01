r"""
Finite Fields of characteristic 2.

This implementation uses NTL's GF2E class to perform the arithmetic
and is the standard implementation for ``GF(2^n)`` for ``n >= 16``.

AUTHORS:

- Martin Albrecht <malb@informatik.uni-bremen.de> (2007-10)
"""

#*****************************************************************************
#       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "cysignals/memory.pxi"
include "cysignals/signals.pxi"
include "sage/libs/ntl/decl.pxi"
from sage.libs.pari.paridecl cimport *

from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport Element, ModuleElement, RingElement

from sage.structure.parent cimport Parent

from sage.rings.ring cimport Ring

from sage.rings.finite_rings.finite_field_base cimport FiniteField

from sage.libs.pari.all import pari
from sage.libs.pari.gen cimport gen

from sage.interfaces.gap import is_GapElement

from sage.misc.randstate import current_randstate

from element_ext_pari import FiniteField_ext_pariElement
from element_pari_ffelt import FiniteFieldElement_pari_ffelt
from finite_field_ntl_gf2e import FiniteField_ntl_gf2e

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

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

cdef int late_import() except -1:
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
        return 0

    import sage.rings.finite_rings.integer_mod
    is_IntegerMod = sage.rings.finite_rings.integer_mod.is_IntegerMod

    import sage.rings.finite_rings.integer_mod_ring
    IntegerModRing_generic = sage.rings.finite_rings.integer_mod_ring.IntegerModRing_generic

    import sage.rings.integer
    Integer = sage.rings.integer.Integer

    import sage.rings.rational
    Rational = sage.rings.rational.Rational

    import sage.rings.polynomial.polynomial_element
    is_Polynomial = sage.rings.polynomial.polynomial_element.is_Polynomial

    import sage.databases.conway
    ConwayPolynomials = sage.databases.conway.ConwayPolynomials

    import sage.rings.finite_rings.conway_polynomials
    conway_polynomial = sage.rings.finite_rings.conway_polynomials.conway_polynomial

    import sage.rings.polynomial.multi_polynomial_element
    MPolynomial = sage.rings.polynomial.multi_polynomial_element.MPolynomial

    import sage.rings.polynomial.polynomial_element
    Polynomial = sage.rings.polynomial.polynomial_element.Polynomial

    import sage.modules.free_module_element
    FreeModuleElement = sage.modules.free_module_element.FreeModuleElement

    import sage.rings.finite_rings.finite_field_constructor
    GF = sage.rings.finite_rings.finite_field_constructor.FiniteField
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

cdef class Cache_ntl_gf2e(SageObject):
    """
    This class stores information for an NTL finite field in a Cython
    class so that elements can access it quickly.

    It's modeled on
    :class:`~sage.rings.finite_rings.integer_mod.NativeIntStruct`,
    but includes many functions that were previously included in
    the parent (see :trac:`12062`).
    """
    def __cinit__(self, parent, Py_ssize_t k, modulus):
        """
        Construction.

        TESTS::

            sage: from sage.rings.finite_rings.element_ntl_gf2e import Cache_ntl_gf2e
            sage: Cache_ntl_gf2e.__new__(Cache_ntl_gf2e, None, 2, [1,1,1])
            <type 'sage.rings.finite_rings.element_ntl_gf2e.Cache_ntl_gf2e'>
        """
        cdef GF2X_c ntl_m
        cdef GF2_c c
        cdef Py_ssize_t i

        for i in range(k + 1):
            GF2_conv_long(c, modulus[i])
            GF2X_SetCoeff(ntl_m, i, c)
        self.F = GF2EContext_c(ntl_m)

    def __init__(self, parent, Py_ssize_t k, modulus):
        """
        Initialization.

        TESTS::

            sage: k.<a> = GF(2^8, impl="ntl")
        """
        self._parent = <FiniteField?>parent
        self._zero_element = self._new()
        GF2E_conv_long((<FiniteField_ntl_gf2eElement>self._zero_element).x,0)
        self._one_element = self._new()
        GF2E_conv_long((<FiniteField_ntl_gf2eElement>self._one_element).x,1)
        if k > 1:
            self._gen = self._new()
            GF2E_from_str(&self._gen.x, "[0 1]")
        elif modulus[0]:
            self._gen = self._one_element
        else:
            self._gen = self._zero_element

        late_import()

    def _doctest_for_5340(self):
        r"""
        Every bug fix should have a doctest.  But :trac:`5340` only happens
        when a garbage collection happens between restoring the modulus and
        using it, so it can't be reliably doctested using any of the
        existing Cython functions in this module.  The sole purpose of
        this method is to doctest the fix for :trac:`5340`.

        EXAMPLES::

            sage: k.<a> = GF(2^20)
            sage: k._cache._doctest_for_5340()
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

    cdef FiniteField_ntl_gf2eElement _new(self):
        """
        Return a new element in ``self``. Use this method to construct
        'empty' elements.
        """
        cdef FiniteField_ntl_gf2eElement y
        self.F.restore()
        y = FiniteField_ntl_gf2eElement.__new__(FiniteField_ntl_gf2eElement)
        y._parent = self._parent
        y._cache = self
        return y

    def order(self):
        """
        Return the cardinality of the field.

        EXAMPLES::

            sage: k.<a> = GF(2^64)
            sage: k._cache.order()
            18446744073709551616
        """
        self.F.restore()
        return Integer(1) << GF2E_degree()

    def degree(self):
        r"""
        If the field has cardinality `2^n` this method returns `n`.

        EXAMPLES::

            sage: k.<a> = GF(2^64)
            sage: k._cache.degree()
            64
        """
        self.F.restore()
        return Integer(GF2E_degree())

    def import_data(self, e):
        """
        EXAMPLES::

            sage: k.<a> = GF(2^17)
            sage: V = k.vector_space()
            sage: v = [1,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0]
            sage: k._cache.import_data(v)
            a^13 + a^8 + a^5 + 1
            sage: k._cache.import_data(V(v))
            a^13 + a^8 + a^5 + 1

        TESTS:

        We check that :trac:`12584` is fixed::

            sage: k(2^63)
            0

        We can coerce from PARI finite field implementations::

            sage: K.<a> = GF(2^19, impl="ntl")
            sage: a^20
            a^6 + a^3 + a^2 + a
            sage: M.<c> = GF(2^19, impl="pari_ffelt")
            sage: K(c^20)
            a^6 + a^3 + a^2 + a
        """
        if isinstance(e, FiniteField_ntl_gf2eElement) and e.parent() is self._parent: return e
        cdef FiniteField_ntl_gf2eElement res = self._new()
        cdef FiniteField_ntl_gf2eElement x
        cdef FiniteField_ntl_gf2eElement g
        cdef Py_ssize_t i

        if is_IntegerMod(e):
            e = e.lift()
        if isinstance(e, int) or \
             isinstance(e, Integer) or \
             isinstance(e, long):
            GF2E_conv_long(res.x,int(e&1))
            return res

        elif isinstance(e, float):
            GF2E_conv_long(res.x,int(e))
            return res

        elif isinstance(e, str):
            return self._parent(eval(e.replace("^","**"),self._parent.gens_dict()))

        elif isinstance(e, FreeModuleElement):
            if self._parent.vector_space() != e.parent():
                raise TypeError, "e.parent must match self.vector_space"
            ztmp = Integer(e.list(),2)
            # Can't do the following since we can't cimport Integer because of circular imports.
            #for i from 0 <= i < len(e):
            #    if e[i]:
            #        mpz_setbit(ztmp.value, i)
            return self.fetch_int(ztmp)

        elif isinstance(e, (list, tuple)):
            if len(e) > self.degree():
                # could reduce here...
                raise ValueError, "list is too long"
            ztmp = Integer(e,2)
            return self.fetch_int(ztmp)

        elif isinstance(e, MPolynomial):
            if e.is_constant():
                return self._parent(e.constant_coefficient())
            else:
                raise TypeError, "no coercion defined"

        elif isinstance(e, Polynomial):
            if e.is_constant():
                return self._parent(e.constant_coefficient())
            else:
                return e(self._parent.gen())

        elif isinstance(e, Rational):
            num = e.numer()
            den = e.denom()
            if den % 2:
                if num % 2:
                    return self._one_element
                return self._zero_element
            raise ZeroDivisionError

        elif isinstance(e, gen):
            pass # handle this in next if clause

        elif isinstance(e, FiniteFieldElement_pari_ffelt) or \
             isinstance(e, FiniteField_ext_pariElement):
            # Reduce to pari
            e = e._pari_()

        elif is_GapElement(e):
            from sage.interfaces.gap import gfq_gap_to_sage
            return gfq_gap_to_sage(e, self._parent)
        else:
            raise TypeError("unable to coerce %r" % type(e))

        cdef GEN t
        if isinstance(e, gen):
            sig_on()
            t = (<gen>e).g
            if typ(t) == t_FFELT:
                t = FF_to_FpXQ(t)
            else:
                t = liftall_shallow(t)

            if typ(t) == t_INT:
                GF2E_conv_long(res.x, itos(t))
                sig_off()
            elif typ(t) == t_POL:
                g = self._gen
                x = self._new()
                GF2E_conv_long(x.x, 1)

                for i from 0 <= i <= degpol(t):
                    if gtolong(gel(t, i+2)):
                        GF2E_add(res.x, res.x, x.x)
                    GF2E_mul(x.x, x.x, g.x)
                sig_off()
            else:
                raise TypeError("bad PARI type %r" % e.type())

            return res

        raise ValueError, "Cannot coerce element %s to this field."%(e)

    cpdef FiniteField_ntl_gf2eElement fetch_int(self, number):
        """
        Given an integer less than `p^n` with base `2`
        representation `a_0 + a_1 \cdot 2 + \cdots + a_k 2^k`, this returns
        `a_0 + a_1 x + \cdots + a_k x^k`, where `x` is the
        generator of this finite field.

        INPUT:

        - ``number`` -- an integer, of size less than the cardinality

        EXAMPLES::

            sage: k.<a> = GF(2^48)
            sage: k._cache.fetch_int(2^33 + 2 + 1)
            a^33 + a + 1

        TESTS:

        We test that :trac:`17027` is fixed::

            sage: K.<a> = GF(2^16)
            sage: K._cache.fetch_int(0r)
            0
        """
        cdef FiniteField_ntl_gf2eElement a = self._new()
        cdef GF2X_c _a

        self.F.restore()

        cdef unsigned char *p
        cdef int i

        if number < 0 or number >= self.order():
            raise TypeError, "n must be between 0 and self.order()"

        if isinstance(number, int) or isinstance(number, long):
            if not number:
                n = 0
            else:
                from sage.misc.functional import log
                n = int(log(number,2))/8 + 1
        elif isinstance(number, Integer):
            n = int(number.nbits())/8 + 1
        else:
            raise TypeError, "number %s is not an integer"%number

        p = <unsigned char*>sig_malloc(n)
        for i from 0 <= i < n:
            p[i] = (number%256)
            number = number >> 8
        GF2XFromBytes(_a, p, n)
        GF2E_conv_GF2X(a.x, _a)
        sig_free(p)
        return a

    def polynomial(self):
        """
        Returns the list of 0's and 1's giving the defining polynomial of the
        field.

        EXAMPLES::

            sage: k.<a> = GF(2^20,modulus="minimal_weight")
            sage: k._cache.polynomial()
            [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        """
        cdef FiniteField_ntl_gf2eElement P
        cdef GF2X_c _P
        cdef GF2_c c
        self.F.restore()
        cdef int i

        P = -(self._gen**(self.degree()))
        _P = GF2E_rep(P.x)

        ret = []
        for i from 0 <= i < GF2E_degree():
            c = GF2X_coeff(_P,i)
            if not GF2_IsZero(c):
                ret.append(1)
            else:
                ret.append(0)
        ret.append(1)
        return ret

cdef class FiniteField_ntl_gf2eElement(FinitePolyExtElement):
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

        - ``parent`` -- base field

        OUTPUT:

        A finite field element.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: from sage.rings.finite_rings.element_ntl_gf2e import FiniteField_ntl_gf2eElement
            sage: FiniteField_ntl_gf2eElement(k)
            0
            sage: k.<a> = GF(2^20)
            sage: a.parent() is k
            True
        """
        if parent is None:
            raise ValueError, "You must provide a parent to construct a finite field element"

    def __cinit__(FiniteField_ntl_gf2eElement self, parent=None ):
        """
        Restores the cache and constructs the underlying NTL element.

        EXAMPLES::

            sage: k.<a> = GF(2^8, impl="ntl") # indirect doctest
        """
        if parent is None:
            return
        if isinstance(parent, FiniteField_ntl_gf2e):
            self._parent = parent
            (<Cache_ntl_gf2e>self._parent._cache).F.restore()

    cdef FiniteField_ntl_gf2eElement _new(FiniteField_ntl_gf2eElement self):
        cdef FiniteField_ntl_gf2eElement y
        (<Cache_ntl_gf2e>self._parent._cache).F.restore()
        y = FiniteField_ntl_gf2eElement.__new__(FiniteField_ntl_gf2eElement)
        y._parent = self._parent
        y._cache = self._cache
        return y

    def __repr__(FiniteField_ntl_gf2eElement self):
        """
        Polynomial representation of ``self``.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: str(a^16) # indirect doctest
            'a^5 + a^3 + a^2 + 1'
            sage: k.<u> = GF(2^16)
            sage: u
            u
        """
        (<Cache_ntl_gf2e>self._parent._cache).F.restore()
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

    def __nonzero__(FiniteField_ntl_gf2eElement self):
        r"""
        Return ``True`` if ``self != k(0)``.

        EXAMPLES::

            sage: k.<a> = GF(2^20)
            sage: bool(a) # indirect doctest
            True
            sage: bool(k(0))
            False
            sage: a.is_zero()
            False
        """
        (<Cache_ntl_gf2e>self._parent._cache).F.restore()
        return not bool(GF2E_IsZero(self.x))

    def is_one(FiniteField_ntl_gf2eElement self):
        r"""
        Return ``True`` if ``self == k(1)``.

        Equivalent to ``self != k(0)``.

        EXAMPLES::

            sage: k.<a> = GF(2^20)
            sage: a.is_one() # indirect doctest
            False
            sage: k(1).is_one()
            True

        """
        (<Cache_ntl_gf2e>self._parent._cache).F.restore()
        return self.x == self._cache._one_element.x

    def is_unit(FiniteField_ntl_gf2eElement self):
        """
        Return ``True`` if ``self`` is nonzero, so it is a unit as an element
        of the finite field.

        EXAMPLES::

            sage: k.<a> = GF(2^17)
            sage: a.is_unit()
            True
            sage: k(0).is_unit()
            False
        """
        (<Cache_ntl_gf2e>self._parent._cache).F.restore()
        if not GF2E_IsZero(self.x):
            return True
        else:
            return False

    def is_square(FiniteField_ntl_gf2eElement self):
        r"""
        Return ``True`` as every element in `\GF{2^n}` is a square.

        EXAMPLES::

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
        Return a square root of this finite field element in its parent.

        EXAMPLES::

            sage: k.<a> = GF(2^20)
            sage: a.is_square()
            True
            sage: a.sqrt()
            a^19 + a^15 + a^14 + a^12 + a^9 + a^7 + a^4 + a^3 + a + 1
            sage: a.sqrt()^2 == a
            True

        This failed before :trac:`4899`::

            sage: GF(2^16,'a')(1).sqrt()
            1

        """
        # this really should be handled special, its gf2 linear after
        # all
        a = self ** (self._cache.order() // 2)
        if all:
            return [a]
        else:
            return a

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Add two elements.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: e = a^2 + a + 1 # indirect doctest
            sage: f = a^15 + a^2 + 1
            sage: e + f
            a^15 + a
        """
        cdef FiniteField_ntl_gf2eElement r = (<FiniteField_ntl_gf2eElement>self)._new()
        GF2E_add(r.x, (<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>right).x)
        return r

    cpdef RingElement _mul_(self, RingElement right):
        """
        Multiply two elements.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: e = a^2 + a + 1 # indirect doctest
            sage: f = a^15 + a^2 + 1
            sage: e * f
            a^15 + a^6 + a^5 + a^3 + a^2
        """
        cdef FiniteField_ntl_gf2eElement r = (<FiniteField_ntl_gf2eElement>self)._new()
        GF2E_mul(r.x, (<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>right).x)
        return r

    cpdef RingElement _div_(self, RingElement right):
        """
        Divide two elements.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: e = a^2 + a + 1 # indirect doctest
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

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Subtract two elements.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: e = a^2 - a + 1 # indirect doctest
            sage: f = a^15 - a^2 + 1
            sage: e - f
            a^15 + a
        """
        cdef FiniteField_ntl_gf2eElement r = (<FiniteField_ntl_gf2eElement>self)._new()
        GF2E_sub(r.x, (<FiniteField_ntl_gf2eElement>self).x, (<FiniteField_ntl_gf2eElement>right).x)
        return r

    def __neg__(FiniteField_ntl_gf2eElement self):
        """
        Return this element.

        EXAMPLES::

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

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: ~a
            a^15 + a^4 + a^2 + a
            sage: a * ~a
            1
        """
        cdef FiniteField_ntl_gf2eElement r = (<FiniteField_ntl_gf2eElement>self)._new()
        cdef FiniteField_ntl_gf2eElement o = (<FiniteField_ntl_gf2eElement>self)._parent._cache._one_element
        GF2E_div(r.x, o.x, (<FiniteField_ntl_gf2eElement>self).x)
        return r

    def __pow__(FiniteField_ntl_gf2eElement self, exp, other):
        """
        EXAMPLES::

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

    cpdef int _cmp_(left, Element right) except -2:
        """
        Comparison of finite field elements.

        .. NOTE::

            Finite fields are unordered. However, we adopt the convention that
            an element ``e`` is bigger than element ``f`` if its polynomial
            representation is bigger.

        EXAMPLES::

            sage: k.<a> = GF(2^20)
            sage: e = k.random_element()
            sage: f = loads(dumps(e))
            sage: e is f
            False
            sage: e == f
            True
            sage: e != (e + 1)
            True

        ::

            sage: K.<a> = GF(2^100)
            sage: a < a^2
            True
            sage: a > a^2
            False
            sage: a+1 > a^2
            False
            sage: a+1 < a^2
            True
            sage: a+1 < a
            False
            sage: a+1 == a
            False
            sage: a == a
            True
        """
        (<Cache_ntl_gf2e>left._parent._cache).F.restore()
        cdef int c = (<FiniteField_ntl_gf2eElement>left).x == (<FiniteField_ntl_gf2eElement>right).x
        if c == 1:
            return 0
        else:
            r = cmp(GF2X_deg(GF2E_rep((<FiniteField_ntl_gf2eElement>left).x)), GF2X_deg(GF2E_rep((<FiniteField_ntl_gf2eElement>right).x)))
            if r:
                return r
            li = left.integer_representation()
            ri = right.integer_representation()
            return cmp(li,ri)

    def _integer_(FiniteField_ntl_gf2eElement self, Integer):
        """
        Convert ``self`` to an integer if it is in the prime subfield.

        EXAMPLES::

            sage: k.<b> = GF(2^16); k
            Finite Field in b of size 2^16
            sage: ZZ(k(1))
            1
            sage: ZZ(k(0))
            0
            sage: ZZ(b)
            Traceback (most recent call last):
            ...
            TypeError: not in prime subfield
            sage: GF(2)(b^(2^16-1)) # indirect doctest
            1
        """
        if self.is_zero(): return Integer(0)
        elif self.is_one(): return Integer(1)
        else:
            raise TypeError, "not in prime subfield"

    def __int__(FiniteField_ntl_gf2eElement self):
        """
        Return the int representation of ``self``.  When ``self`` is in the
        prime subfield, the integer returned is equal to ``self`` and
        otherwise raises an error.

        EXAMPLES::

            sage: k.<a> = GF(2^20)
            sage: int(k(0))
            0
            sage: int(k(1))
            1
            sage: int(a)
            Traceback (most recent call last):
            ...
            TypeError: Cannot coerce element to an integer.
            sage: int(a^2 + 1)
            Traceback (most recent call last):
            ...
            TypeError: Cannot coerce element to an integer.

        """
        if self == 0:
            return int(0)
        elif self == 1:
            return int(1)
        else:
            raise TypeError("Cannot coerce element to an integer.")

    def integer_representation(FiniteField_ntl_gf2eElement self):
        r"""
        Return the int representation of ``self``.  When ``self`` is in the
        prime subfield, the integer returned is equal to ``self`` and not
        to ``log_repr``.

        Elements of this field are represented as ints in as follows:
        for `e \in \GF{p}[x]` with `e = a_0 + a_1 x + a_2 x^2 + \cdots`,
        `e` is represented as: `n = a_0 + a_1  p + a_2  p^2 + \cdots`.

        EXAMPLES::

            sage: k.<a> = GF(2^20)
            sage: a.integer_representation()
            2
            sage: (a^2 + 1).integer_representation()
            5
            sage: k.<a> = GF(2^70)
            sage: (a^65 + a^64 + 1).integer_representation()
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
        Return ``self`` viewed as a polynomial over
        ``self.parent().prime_subfield()``.

        INPUT:

        - ``name`` -- (optional) variable name

        EXAMPLES::

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
        Return the characteristic polynomial of ``self`` as a polynomial
        in var over the prime subfield.

        INPUT:

        - ``var`` -- string (default: ``'x'``)

        OUTPUT:

        polynomial

        EXAMPLES::

            sage: k.<a> = GF(2^8, impl="ntl")
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
        Return the minimal polynomial of ``self``, which is the smallest
        degree polynomial `f \in \GF{2}[x]` such that ``f(self) == 0``.

        INPUT:

        - ``var`` -- string (default: ``'x'``)

        OUTPUT:

        polynomial

        EXAMPLES::

            sage: K.<a> = GF(2^100)
            sage: f = a.minpoly(); f
            x^100 + x^57 + x^56 + x^55 + x^52 + x^48 + x^47 + x^46 + x^45 + x^44 + x^43 + x^41 + x^37 + x^36 + x^35 + x^34 + x^31 + x^30 + x^27 + x^25 + x^24 + x^22 + x^20 + x^19 + x^16 + x^15 + x^11 + x^9 + x^8 + x^6 + x^5 + x^3 + 1
            sage: f(a)
            0
            sage: g = K.random_element()
            sage: g.minpoly()(g)
            0
        """
        (<Cache_ntl_gf2e>self._parent._cache).F.restore()
        cdef GF2X_c r = GF2X_IrredPolyMod(GF2E_rep(self.x), GF2E_modulus())
        cdef int i
        C = []
        for i from 0 <= i <= GF2X_deg(r):
            C.append(GF2_conv_to_long(GF2X_coeff(r,i)))
        return self._parent.polynomial_ring(var)(C)

    def trace(self):
        """
        Return the trace of ``self``.

        EXAMPLES::

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
        representation of ``self``.

        EXAMPLES::

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

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` that MAGMA can
        understand.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: a._magma_init_(magma)      # random; optional - magma
            '_sage_[...]'

        .. NOTE::

            This method calls MAGMA to setup the parent.
        """
        km = magma(self.parent())
        vn_m = km.gen(1).name()
        vn_s = str(self.parent().polynomial_ring().gen())
        return str(self.polynomial()).replace(vn_s,vn_m)

    def __copy__(self):
        """
        Return a copy of this element.  Actually just returns ``self``, since
        finite field elements are immutable.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: copy(a) is a
            True
        """
        return self

    def _gap_init_(self):
        r"""
        Return a string that evaluates to the GAP representation of
        this element.

        A ``NotImplementedError`` is raised if
        ``self.parent().modulus()`` is not a Conway polynomial, as
        the isomorphism of finite fields is not implemented yet.

        EXAMPLES::

            sage: k.<b> = GF(2^16)
            sage: b._gap_init_()
            'Z(65536)^1'
        """
        F = self._parent
        if not F.is_conway():
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

        EXAMPLES::

            sage: k.<a> = GF(2^18)
            sage: {a:1,a:0} # indirect doctest
            {a: 0}
        """
        return hash(self.integer_representation()) # todo, come up with a faster version

    def _vector_(FiniteField_ntl_gf2eElement self, reverse=False):
        r"""
        Return a vector in ``self.parent().vector_space()``
        matching ``self``. The most significant bit is to the
        right.

        INPUT:

        - ``reverse`` -- reverse the order of the bits from little endian to
          big endian.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: e = a^14 + a^13 + 1
            sage: vector(e) # little endian
            (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0)

            sage: e._vector_(reverse=True) # big endian
            (0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
        """
        #vector(foo) might pass in ZZ
        if isinstance(reverse, Parent):
            raise TypeError, "Base field is fixed to prime subfield."

        cdef GF2X_c r = GF2E_rep(self.x)
        cdef int i

        (<Cache_ntl_gf2e>self._parent._cache).F.restore()

        C = []
        for i from 0 <= i < GF2E_degree():
            C.append(GF2_conv_to_long(GF2X_coeff(r,i)))
        if reverse:
            C = list(reversed(C))
        return self._parent.vector_space()(C)

    def __reduce__(FiniteField_ntl_gf2eElement self):
        """
        Used for supporting pickling of finite field elements.

        EXAMPLES::

            sage: k.<a> = GF(2^16)
            sage: loads(dumps(a)) == a
            True
        """
        return unpickleFiniteField_ntl_gf2eElement, (self._parent, str(self))

    def log(self, base):
        """
        Return `x` such that `b^x = a`, where `x` is `a` and `b`
        is the base.

        INPUT:

        - ``base`` -- finite field element that generates the multiplicative
          group.

        OUTPUT:

        Integer `x` such that `a^x = b`, if it exists.
        Raises a ``ValueError`` exception if no such `x` exists.

        EXAMPLES::

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
    EXAMPLES::

        sage: k.<a> = GF(2^20)
        sage: e = k.random_element()
        sage: f = loads(dumps(e)) # indirect doctest
        sage: e == f
        True
    """
    return parent(elem)

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.rings.finite_field_ntl_gf2e', 'unpickleFiniteField_ntl_gf2eElement', unpickleFiniteField_ntl_gf2eElement)
