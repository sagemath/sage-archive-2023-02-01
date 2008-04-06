"""
Conversion routines from to SINGULAR's native data structures

AUTHOR: Martin Albrecht <malb@informatik.uni-bremen.de>

TODO: Figure out how to do the cdef public/extern stuff with C++
"""
###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

include "sage/libs/ntl/decl.pxi"
include "sage/ext/stdsage.pxi"

cdef extern from "limits.h":
    long INT_MAX
    long INT_MIN

from sage.rings.rational_field import RationalField
from sage.rings.finite_field_prime_modn import FiniteField_prime_modn
from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
from sage.libs.pari.all import pari

from sage.structure.parent_base cimport ParentWithBase
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular

cdef class Conversion:
    """
    A work around class to export the contained methods/functions
    """

    cdef public Rational si2sa_QQ(self, number *n, ring *_ring):
        """
        Converts a SINGULAR rational number to a SAGE rational number.

        INPUT:
            n -- number
            _ring -- singular ring, used to check type of n

        OUTPUT:
            SAGE rational number matching n
        """
        cdef number *nom
        cdef number *denom
        cdef mpq_t _z

        cdef mpz_t nom_z, denom_z

        cdef Rational z

        mpq_init(_z)

        ##  Immediate integers handles carry the tag 'SR_INT', i.e. the last bit is 1.
        ##  This distuingishes immediate integers from other handles which  point  to
        ##  structures aligned on 4 byte boundaries and therefor have last bit  zero.
        ##  (The second bit is reserved as tag to allow extensions of  this  scheme.)
        ##  Using immediates as pointers and dereferencing them gives address errors.
        nom = nlGetNom(n, _ring)
        mpz_init(nom_z)

        if (SR_HDL(nom) & SR_INT): mpz_set_si(nom_z, SR_TO_INT(nom))
        else: mpz_set(nom_z,&nom.z)

        mpq_set_num(_z,nom_z)
        n_Delete(&nom,_ring)
        mpz_clear(nom_z)

        denom = nlGetDenom(n, _ring)
        mpz_init(denom_z)

        if (SR_HDL(denom) & SR_INT): mpz_set_si(denom_z, SR_TO_INT(denom))
        else: mpz_set(denom_z,&denom.z)

        mpq_set_den(_z, denom_z)
        n_Delete(&denom,_ring)
        mpz_clear(denom_z)

        z = Rational()
        z.set_from_mpq(_z)
        mpq_clear(_z)
        return z

    cdef public FiniteField_givaroElement si2sa_GFqGivaro(self, number *n, ring *_ring, FiniteField_givaro base):
        """
        Convert a SINGULAR finite extension field element to a SAGE
        finite extension field element in a field with order
        $<2^{16}$.

        INPUT:
            n -- SINGULAR representation
            _ring -- SINGULAR ring
            base -- SAGE GF(q)

        OUTPUT:
            An Element in GF(q).

        EXAMPLE:
            sage: K.<a> = GF(5^3)
            sage: R.<x,y,z> = PolynomialRing(K)
            sage: K( (4*R(a)^2 + R(a))^3 )
            a^2
        """
        cdef napoly *z
        cdef int c, e
        cdef int a
        cdef int ret
        cdef int order

        if naIsZero(n):
            return base._zero_element
        elif naIsOne(n):
            return base._one_element
        z = (<lnumber*>n).z

        a = base.objectptr.sage_generator()
        ret = base.objectptr.zero
        order = base.objectptr.cardinality() - 1

        while z:
            c = base.objectptr.initi(c,<long>napGetCoeff(z))
            e = napGetExp(z,1)
            if e == 0:
                ret = base.objectptr.add(ret, c, ret)
            else:
                a = ( e * base.objectptr.sage_generator() ) % order
                ret = base.objectptr.axpy(ret, c, a, ret)
            z = napIter(z)
        return (<FiniteField_givaroElement>base._zero_element)._new_c(ret)

    cdef public FiniteField_ntl_gf2eElement si2sa_GFqNTLGF2E(self, number *n, ring *_ring, FiniteField_ntl_gf2e base):
        """
        Convert a SINGULAR finite extension field element to a SAGE
        finite extension field element implemented via NTL's GF2E

        INPUT:
            n -- SINGULAR representation
            _ring -- SINGULAR ring
            base -- SAGE GF(q)

        OUTPUT:
            An Element in GF(q).

        """
        cdef napoly *z
        cdef int c, e
        cdef FiniteField_ntl_gf2eElement a
        cdef FiniteField_ntl_gf2eElement ret

        if naIsZero(n):
            return base._zero_element
        elif naIsOne(n):
            return base._one_element
        z = (<lnumber*>n).z

        a = base.gen()
        ret = base._zero_element

        while z:
            c = <long>napGetCoeff(z)
            e = napGetExp(z,1)
            ret += c * a**e
            z = napIter(z)
        return ret

    cdef public object si2sa_GFqPari(self, number *n, ring *_ring, object base):
        """
        Convert a SINGULAR finite extension field element to a SAGE
        finite extension field element in a field.

        INPUT:
            n -- SINGULAR representation
            _ring -- SINGULAR ring
            base -- SAGE GF(q)

        OUTPUT:
            An Element in GF(q).
        """
        cdef napoly *z
        cdef int c, e
        cdef object a
        cdef object ret

        if naIsZero(n):
            return base._zero_element
        elif naIsOne(n):
            return base._one_element
        z = (<lnumber*>n).z

        a = pari("a")
        ret = pari(int(0)).Mod(int(_ring.ch))

        while z:
            c = <long>napGetCoeff(z)
            e = napGetExp(z,1)
            if e == 0:
                ret = ret + c
            elif c != 0:
                ret = ret  + c * a**e
            z = napIter(z)
        return base(ret)



    cdef public number *sa2si_QQ(self, Rational r, ring *_ring):
        """
        Convert a SAGE rational number to a SINGULAR rational number.

        INPUT:
            r -- SAGE notation
            _ring -- SINGULAR ring (ignored)
        """
        return nlInit2gmp( mpq_numref(r.value), mpq_denref(r.value) )

    cdef number *sa2si_GFqGivaro(self, int quo, ring *_ring):
        """
        Convert a SAGE GF(q) element to a SINGULAR GF(q) element.

        INPUT:
            quo -- int representing the finite field element
            _ring -- SINGULAR ring
        """
        cdef number *n1, *n2, *a, *coeff, *apow1, *apow2
        cdef int b

        rChangeCurrRing(_ring)
        b   = - _ring.ch;

        a = naPar(1)

        apow1 = naInit(1)
        n1 = naInit(0)

        while quo!=0:
            coeff = naInit(quo%b)

            if not naIsZero(coeff):
                n2 = naAdd( naMult(coeff, apow1),  n1)
                naDelete(&n1, _ring);
                n1= n2

            apow2 = naMult(apow1, a)
            naDelete(&apow1, _ring)
            apow1 = apow2

            quo = quo/b
            naDelete(&coeff, _ring)

        naDelete(&apow1, _ring)
        naDelete(&a, _ring)
        return n1

    cdef number *sa2si_GFqNTLGF2E(self, FiniteField_ntl_gf2eElement elem, ring *_ring):
        """
        """
        cdef int i
        cdef number *n1, *n2, *a, *coeff, *apow1, *apow2

        rChangeCurrRing(_ring)

        cdef GF2X_c rep = GF2E_rep(elem.x)

        if GF2X_deg(rep) >= 1:
            n1 = naInit(0)
            a = naPar(1)
            apow1 = naInit(1)

            for i from 0 <= i <= GF2X_deg(rep):
                coeff = naInit(GF2_conv_to_long(GF2X_coeff(rep,i)))

                if not naIsZero(coeff):
                    n2 = naAdd( naMult(coeff, apow1),  n1)
                    naDelete(&n1, _ring);
                    n1= n2

                apow2 = naMult(apow1, a)
                naDelete(&apow1, _ring)
                apow1 = apow2

                naDelete(&coeff, _ring)

            naDelete(&apow1, _ring)
            naDelete(&a, _ring)
        else:
            n1 = naInit(GF2_conv_to_long(GF2X_coeff(rep,i)))

        return n1

    cdef number *sa2si_GFqPari(self, object elem, ring *_ring):
        """
        """
        cdef int i
        cdef number *n1, *n2, *a, *coeff, *apow1, *apow2

        rChangeCurrRing(_ring)

        elem = elem._pari_().lift().lift()


        if len(elem) > 1:
            n1 = naInit(0)
            a = naPar(1)
            apow1 = naInit(1)

            for i from 0 <= i < len(elem):
                coeff = naInit(int(elem[i]))

                if not naIsZero(coeff):
                    n2 = naAdd( naMult(coeff, apow1),  n1)
                    naDelete(&n1, _ring);
                    n1= n2

                apow2 = naMult(apow1, a)
                naDelete(&apow1, _ring)
                apow1 = apow2

                naDelete(&coeff, _ring)

            naDelete(&apow1, _ring)
            naDelete(&a, _ring)
        else:
            n1 = naInit(int(elem))

        return n1

    cdef public number *sa2si_ZZ(self, Integer d, ring *_ring):
        """
        """
        cdef number *n
        if INT_MIN <= d <= INT_MAX:
            return nlInit(int(d))
        else:
            n = nlRInit(0)
            mpz_init_set(&n.z, d.value)
            return n

    cdef public object si2sa(self, number *n, ring *_ring, object base):
        if PY_TYPE_CHECK(base, FiniteField_prime_modn):
            return base(nInt(n))

        elif PY_TYPE_CHECK(base, RationalField):
            return self.si2sa_QQ(n,_ring)

        elif PY_TYPE_CHECK(base, FiniteField_givaro):
            return self.si2sa_GFqGivaro(n, _ring, base)

        elif PY_TYPE_CHECK(base, FiniteField_ext_pari):
            return self.si2sa_GFqPari(n, _ring, base)

        elif PY_TYPE_CHECK(base, FiniteField_ntl_gf2e):
            return self.si2sa_GFqNTLGF2E(n, _ring, base)

        else:
            raise ValueError, "cannot convert from SINGULAR number"

    cdef public number *sa2si(self, Element elem, ring * _ring):
        cdef int i
        if PY_TYPE_CHECK(elem._parent, FiniteField_prime_modn):
            return n_Init(int(elem),_ring)

        elif PY_TYPE_CHECK(elem._parent, RationalField):
            return self.sa2si_QQ(elem, _ring)

        elif PY_TYPE_CHECK(elem._parent, FiniteField_givaro):
            return self.sa2si_GFqGivaro( (<FiniteField_givaro>elem._parent).objectptr.convert(i, (<FiniteField_givaroElement>elem).element ), _ring )

        elif PY_TYPE_CHECK(elem._parent, FiniteField_ext_pari):
            return self.sa2si_GFqPari(elem, _ring)

        elif PY_TYPE_CHECK(elem._parent, FiniteField_ntl_gf2e):
            return self.sa2si_GFqNTLGF2E(elem, _ring)

        else:
            raise ValueError, "cannot convert to SINGULAR number"

    cdef public  MPolynomial_libsingular new_MP(self, MPolynomialRing_libsingular parent, poly *juice):
        """
        Construct MPolynomial_libsingular from parent and SINGULAR poly.
        """
        cdef MPolynomial_libsingular p
        p = PY_NEW(MPolynomial_libsingular)
        p._parent = <ParentWithBase>parent
        p._poly = juice
        p_Normalize(p._poly, parent._ring)
        return p

