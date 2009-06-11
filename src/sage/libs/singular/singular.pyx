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
from sage.rings.integer_ring cimport IntegerRing_class
from sage.rings.integer_mod_ring import IntegerModRing_generic
from sage.rings.finite_field_prime_modn import FiniteField_prime_modn
from sage.rings.finite_field_ext_pari import FiniteField_ext_pari
from sage.libs.pari.all import pari

from sage.structure.parent_base cimport ParentWithBase
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular

cdef Rational si2sa_QQ(number *n, ring *_ring):
    """
    TESTS::

        sage: P.<x,y,z> = QQ[]
        sage: P(1/3).lc()
        1/3
        sage: P(1).lc()
        1
        sage: P(0).lc()
        0
        sage: P(-1/3).lc()
        -1/3
        sage: type(P(3).lc())
        <type 'sage.rings.rational.Rational'>
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
    nlDelete(&nom,_ring)
    mpz_clear(nom_z)

    denom = nlGetDenom(n, _ring)
    mpz_init(denom_z)

    if (SR_HDL(denom) & SR_INT): mpz_set_si(denom_z, SR_TO_INT(denom))
    else: mpz_set(denom_z,&denom.z)

    mpq_set_den(_z, denom_z)
    nlDelete(&denom,_ring)
    mpz_clear(denom_z)

    z = Rational()
    z.set_from_mpq(_z)
    mpq_clear(_z)
    return z

cdef Integer si2sa_ZZ(number *n, ring *_ring):
    """
    TESTS::

        sage: P.<x,y,z> = ZZ[]
        sage: P(3).lc()
        3
        sage: P(0).lc()
        0
        sage: P(-3).lc()
        -3
        sage: P(-1234567890).lc()
        -1234567890
        sage: type(P(3).lc())
        <type 'sage.rings.integer.Integer'>
    """
    cdef Integer z
    z = Integer()
    z.set_from_mpz(<__mpz_struct*>n)
    return z

cdef FiniteField_givaroElement si2sa_GFqGivaro(number *n, ring *_ring, FiniteField_givaro base):
    """
    TESTS::

        sage: K.<a> = GF(5^3)
        sage: R.<x,y,z> = PolynomialRing(K)
        sage: K( (4*R(a)^2 + R(a))^3 )
        a^2
        sage: K(R(0))
        0
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

cdef FiniteField_ntl_gf2eElement si2sa_GFqNTLGF2E(number *n, ring *_ring, FiniteField_ntl_gf2e base):
    """
    TESTS::

        sage: K.<a> = GF(2^20)
        sage: P.<x,y,z> = K[]
        sage: f = a^21*x^2 + 1 # indirect doctest
        sage: f.lc()
        a^11 + a^10 + a^8 + a^7 + a^6 + a^5 + a^2 + a
        sage: type(f.lc())
        <type 'sage.rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2eElement'>
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

cdef object si2sa_GFqPari(number *n, ring *_ring, object base):
    """
    TESTS::

        sage: K.<a> = GF(3^16)
        sage: P.<x,y,z> = K[]
        sage: f = a^21*x^2 + 1 # indirect doctest
        sage: f.lc()
        a^12 + a^11 + a^9 + a^8 + a^7 + 2*a^6 + a^5
        sage: type(f.lc())
        <class 'sage.rings.finite_field_element.FiniteField_ext_pariElement'>
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

cdef object si2sa_NF(number *n, ring *_ring, object base):
    """
    TESTS::

        sage: K.<a> = NumberField(x^2 - 2)
        sage: P.<x,y,z> = K[]
        sage: f = a^21*x^2 + 1 # indirect doctest
        sage: f.lc()
        1024*a
        sage: type(f.lc())
        <type 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic'>
    """
    cdef napoly *z
    cdef number *c
    cdef int e
    cdef object a
    cdef object ret

    if naIsZero(n):
        return base._zero_element
    elif naIsOne(n):
        return base._one_element
    z = (<lnumber*>n).z

    a = base.gen()
    ret = base(0)

    while z:
        c = napGetCoeff(z)
        coeff = si2sa_QQ(c, _ring)
        e = napGetExp(z,1)
        if e == 0:
            ret = ret + coeff
        elif coeff != 0:
            ret = ret  + coeff * a**e
        z = napIter(z)
    return base(ret)

cdef inline object si2sa_ZZmod(number *n, ring *_ring, object base):
    """
    TESTS::

        sage: P.<x,y,z> = Integers(10)[]
        sage: P(3).lc()
        3
        sage: P(13).lc()
        3

        sage: P.<x,y,z> = Integers(16)[]
        sage: P(3).lc()
        3
        sage: P(19).lc()
        3

        sage: P.<x,y,z> = Integers(3**2)[]
        sage: P(3).lc()
        3
        sage: P(12).lc()
        3

        sage: P.<x,y,z> = Integers(2^32)[]
        sage: P(2^32-1).lc()
        4294967295

        sage: P(3).lc()
        3

        sage: P.<x,y,z> = Integers(17^20)[]
        sage: P(17^19 + 3).lc()
        239072435685151324847156

        sage: P(3)
        3
    """
    cdef Integer ret
    if _ring.ringtype == 1:
        return base(<long>n)
    else:
        ret = Integer()
        ret.set_from_mpz(<__mpz_struct*>n)
        return base(ret)

    return base(_ring.cf.nInt(n))

cdef number *sa2si_QQ(Rational r, ring *_ring):
    """
    TESTS::

        sage: P.<x,y,z> = QQ[]
        sage: P(0) + 1/2 - 2/4
        0
        sage: P(1/2) + 3/5 - 3/5
        1/2
        sage: P(2/3) + 1/4 - 1/4
        2/3
        sage: P(12345678901234567890/23) + 5/2 - 5/2
        12345678901234567890/23
    """
    return nlInit2gmp( mpq_numref(r.value), mpq_denref(r.value) )

cdef number *sa2si_GFqGivaro(int quo, ring *_ring):
    """
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

cdef number *sa2si_GFqNTLGF2E(FiniteField_ntl_gf2eElement elem, ring *_ring):
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
       n1 = naInit(GF2_conv_to_long(GF2X_coeff(rep,0)))

    return n1

cdef number *sa2si_GFqPari(object elem, ring *_ring):
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

cdef number *sa2si_NF(object elem, ring *_ring):
    """
    """
    cdef int i
    cdef number *n1, *n2, *a, *nlCoeff, *naCoeff, *apow1, *apow2

    rChangeCurrRing(_ring)

    elem = list(elem)

    if len(elem) > 1:
        n1 = naInit(0)
        a = naPar(1)
        apow1 = naInit(1)

        for i from 0 <= i < len(elem):
            nlCoeff = nlInit2gmp( mpq_numref((<Rational>elem[i]).value), mpq_denref((<Rational>elem[i]).value) )
            naCoeff = naMap00(nlCoeff)
            nlDelete(&nlCoeff, _ring)

            # faster would be to assign the coefficient directly
            n2 = naAdd( naMult(naCoeff, apow1),  n1)
            naDelete(&n1, _ring);
            naDelete(&naCoeff, _ring)
            n1 = n2

            apow2 = naMult(apow1, a)
            naDelete(&apow1, _ring)
            apow1 = apow2

        naDelete(&apow1, _ring)
        naDelete(&a, _ring)
    else:
        n1 = sa2si_QQ(elem[0], _ring)

    return n1

cdef number *sa2si_ZZ(Integer d, ring *_ring):
    """
    TESTS::

        sage: P.<x,y,z> = ZZ[]
        sage: P(0) + 1 - 1
        0
        sage: P(1) + 1 - 1
        1
        sage: P(2) + 1 - 1
        2
        sage: P(12345678901234567890) + 2 - 2
        12345678901234567890
    """
    cdef number *n = nrzInit(0)
    mpz_set(<__mpz_struct*>n, d.value)
    return <number*>n

cdef inline number *sa2si_ZZmod(IntegerMod_abstract d, ring *_ring):
    """
    TESTS::

        sage: P.<x,y,z> = Integers(10)[]
        sage: P(3)
        3
        sage: P(13)
        3

        sage: P.<x,y,z> = Integers(16)[]
        sage: P(3)
        3
        sage: P(19)
        3

        sage: P.<x,y,z> = Integers(3^2)[]
        sage: P(3)
        3
        sage: P(12)
        3

        sage: P.<x,y,z> = Integers(2^32)[]
        sage: P(2^32-1)
        -1

        sage: P(3)
        3

        sage: P.<x,y,z> = Integers(17^20)[]
        sage: P(17^19 + 3)
        239072435685151324847156

        sage: P(3)
        3

    """
    cdef long _d
    if _ring.ringtype == 1:
        _d = long(d)
        return nr2mMapZp(<number *>_d)
    else:
        return nrnMapGMP(<number *>(<Integer>Integer(d)).value)

cdef object si2sa(number *n, ring *_ring, object base):
    if PY_TYPE_CHECK(base, FiniteField_prime_modn):
        return base(nInt(n))

    elif PY_TYPE_CHECK(base, RationalField):
        return si2sa_QQ(n,_ring)

    elif PY_TYPE_CHECK(base, IntegerRing_class):
        return si2sa_ZZ(n,_ring)

    elif PY_TYPE_CHECK(base, FiniteField_givaro):
        return si2sa_GFqGivaro(n, _ring, base)

    elif PY_TYPE_CHECK(base, FiniteField_ext_pari):
        return si2sa_GFqPari(n, _ring, base)

    elif PY_TYPE_CHECK(base, FiniteField_ntl_gf2e):
        return si2sa_GFqNTLGF2E(n, _ring, base)

    elif PY_TYPE_CHECK(base, NumberField) and base.is_absolute():
        return si2sa_NF(n, _ring, base)

    elif PY_TYPE_CHECK(base, IntegerModRing_generic):
        if _ring.ringtype == 0:
            return base(nInt(n))
        return si2sa_ZZmod(n, _ring, base)

    else:
        raise ValueError, "cannot convert from SINGULAR number"

cdef number *sa2si(Element elem, ring * _ring):
    cdef int i
    if PY_TYPE_CHECK(elem._parent, FiniteField_prime_modn):
        return n_Init(int(elem),_ring)

    elif PY_TYPE_CHECK(elem._parent, RationalField):
        return sa2si_QQ(elem, _ring)

    elif PY_TYPE_CHECK(elem._parent, IntegerRing_class):
        return sa2si_ZZ(elem, _ring)

    elif PY_TYPE_CHECK(elem._parent, FiniteField_givaro):
        return sa2si_GFqGivaro( (<FiniteField_givaro>elem._parent).objectptr.convert(i, (<FiniteField_givaroElement>elem).element ), _ring )

    elif PY_TYPE_CHECK(elem._parent, FiniteField_ext_pari):
        return sa2si_GFqPari(elem, _ring)

    elif PY_TYPE_CHECK(elem._parent, FiniteField_ntl_gf2e):
        return sa2si_GFqNTLGF2E(elem, _ring)

    elif PY_TYPE_CHECK(elem._parent, NumberField) and elem._parent.is_absolute():
        return sa2si_NF(elem, _ring)
    elif PY_TYPE_CHECK(elem._parent, IntegerModRing_generic):
        if _ring.ringtype == 0:
            return n_Init(int(elem),_ring)
        return sa2si_ZZmod(elem, _ring)
    else:
        raise ValueError, "cannot convert to SINGULAR number"
