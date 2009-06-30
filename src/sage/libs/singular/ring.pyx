"""
Wrapper for Singular's Rings

AUTHOR:

- Martin Albrecht (2009-07): initial implementation
"""
#*****************************************************************************
#       Copyright (C) 2009 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../../ext/stdsage.pxi"

from sage.libs.gmp.types cimport __mpz_struct
from sage.libs.gmp.mpz cimport mpz_init_set_ui, mpz_init_set

from sage.libs.singular.decl cimport number, lnumber, napoly, ring, currRing
from sage.libs.singular.decl cimport rChangeCurrRing, rCopy0, rComplete, rDelete
from sage.libs.singular.decl cimport omAlloc0, omStrDup, omAlloc, omAlloc0Bin,  sip_sring_bin, rnumber_bin
from sage.libs.singular.decl cimport ringorder_dp, ringorder_Dp, ringorder_lp, ringorder_rp, ringorder_ds, ringorder_Ds, ringorder_ls, ringorder_C
from sage.libs.singular.decl cimport p_Copy

from sage.rings.integer cimport Integer
from sage.rings.integer_ring cimport IntegerRing_class
from sage.rings.integer_ring import ZZ
from sage.rings.integer_mod_ring import is_IntegerModRing
from sage.rings.number_field.number_field_base cimport NumberField
from sage.rings.rational_field import RationalField
from sage.rings.ring import FiniteField as FiniteField_generic

from sage.rings.polynomial.term_order import TermOrder
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular, MPolynomialRing_libsingular


from sage.misc.misc_c import is_64_bit


# mapping str --> SINGULAR representation
order_dict = {"dp":ringorder_dp,
              "Dp":ringorder_Dp,
              "lp":ringorder_lp,
              "rp":ringorder_rp,
              "ds":ringorder_ds,
              "Ds":ringorder_Ds,
              "ls":ringorder_ls,
              }

cdef ring *singular_ring_new(base_ring, n, names, term_order) except NULL:
    """
    Create a new Singular ring over the ``base_ring`` in ``n``
    variables with the names ``names`` and the term order
    ``term_order``.

    INPUT:

    - ``base_ring`` - a Sage ring

    - ``n`` - the number of variables (> 0)

    - ``names`` - a list of names of length ``n``

    - ``term_order`` - a term ordering

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: P
        Multivariate Polynomial Ring in x, y, z over Rational Field

        sage: P.term_order()
        Degree reverse lexicographic term order

        sage: P = PolynomialRing(GF(127),3,names='abc', order='lex')
        sage: P
        Multivariate Polynomial Ring in a, b, c over Finite Field of size 127

        sage: P.term_order()
        Lexicographic term order

        sage: z = QQ['z'].0
        sage: K.<s> = NumberField(z^2 - 2)
        sage: P.<x,y> = PolynomialRing(K, 2)

        sage: P.<x,y,z> = ZZ[]; P
        Multivariate Polynomial Ring in x, y, z over Integer Ring

        sage: P.<x,y,z> = Zmod(2^10)[]; P
        Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1024

        sage: P.<x,y,z> = Zmod(3^10)[]; P
        Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 59049

        sage: P.<x,y,z> = Zmod(2^100)[]; P
        Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 1267650600228229401496703205376

        sage: P.<x,y,z> = Zmod(2521352)[]; P
        Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 2521352

        sage: P.<x,y,z> = Zmod(25213521351515232)[]; P
        Multivariate Polynomial Ring in x, y, z over Ring of integers modulo 25213521351515232
    """
    cdef ring* _ring
    cdef char **_names
    cdef char *_name
    cdef int i
    cdef int nblcks
    cdef int offset
    cdef int characteristic
    cdef int ringtype = 0
    cdef MPolynomialRing_libsingular k
    cdef MPolynomial_libsingular minpoly
    cdef lnumber *nmp

    cdef __mpz_struct* ringflaga
    cdef unsigned long ringflagb

    is_extension = False

    n = int(n)
    if n<1:
        raise ArithmeticError("The number of variables must be at least 1.")

    order = TermOrder(term_order, n)

    _names = <char**>omAlloc0(sizeof(char*)*(len(names)))
    for i from 0 <= i < n:
        _name = names[i]
        _names[i] = omStrDup(_name)

    # from the SINGULAR source code documentation for the rInit function
    ##  characteristic --------------------------------------------------
    ##  input: 0 ch=0 : Q     parameter=NULL    ffChar=FALSE   float_len (done)
    ##         0    1 : Q(a,...)        *names         FALSE             (done)
    ##         0   -1 : R               NULL           FALSE  0
    ##         0   -1 : R               NULL           FALSE  prec. >6
    ##         0   -1 : C               *names         FALSE  prec. 0..?
    ##         p    p : Fp              NULL           FALSE             (done)
    ##         p   -p : Fp(a)           *names         FALSE             (done)
    ##         q    q : GF(q=p^n)       *names         TRUE              (todo)

    if base_ring.is_field() and base_ring.is_finite() and base_ring.is_prime_field():
        if base_ring.characteristic() <= 2147483629:
            characteristic = base_ring.characteristic()
        else:
            raise TypeError, "Characteristic p must be <= 2147483629."

    elif PY_TYPE_CHECK(base_ring, RationalField):
        characteristic = 0

    elif PY_TYPE_CHECK(base_ring, IntegerRing_class):
        ringflaga = <__mpz_struct*>omAlloc(sizeof(__mpz_struct))
        mpz_init_set_ui(ringflaga, 0)
        characteristic = 0
        ringtype = 4 # integer ring

    elif PY_TYPE_CHECK(base_ring, FiniteField_generic):
        if base_ring.characteristic() <= 2147483629:
            characteristic = -base_ring.characteristic() # note the negative characteristic
        else:
            raise TypeError, "characteristic must be <= 2147483629."
        # TODO: This is lazy, it should only call Singular stuff not MPolynomial stuff
        k = MPolynomialRing_libsingular(base_ring.prime_subfield(), 1, base_ring.variable_name(), 'lex')
        minpoly = base_ring.polynomial()(k.gen())
        is_extension = True

    elif PY_TYPE_CHECK(base_ring, NumberField) and base_ring.is_absolute():
        characteristic = 1
        k = MPolynomialRing_libsingular(RationalField(), 1, base_ring.variable_name(), 'lex')
        minpoly = base_ring.polynomial()(k.gen())
        is_extension = True

    elif is_IntegerModRing(base_ring):
        ch = base_ring.characteristic()
        if ch.is_power_of(2):
            exponent = ch.nbits() -1
            if is_64_bit:
                # it seems Singular uses ints somewhere
                # internally, cf. #6051 (Sage) and #138 (Singular)
                if exponent <= 30: ringtype = 1
                else: ringtype = 3
            else:
                if exponent <= 30: ringtype = 1
                else: ringtype = 3
            characteristic = exponent
            ringflaga = <__mpz_struct*>omAlloc(sizeof(__mpz_struct))
            mpz_init_set_ui(ringflaga, 2)
            ringflagb = exponent

        elif base_ring.characteristic().is_prime_power()  and ch < ZZ(2)**160:
            F = ch.factor()
            assert(len(F)==1)

            ringtype = 3
            ringflaga = <__mpz_struct*>omAlloc(sizeof(__mpz_struct))
            mpz_init_set(ringflaga, (<Integer>F[0][0]).value)
            ringflagb = F[0][1]
            characteristic = F[0][1]

        else:
            # normal modulus
            try:
                characteristic = ch
            except OverflowError:
                raise NotImplementedError("Characteristic %d too big."%ch)
            ringtype = 2
            ringflaga = <__mpz_struct*>omAlloc(sizeof(__mpz_struct))
            mpz_init_set_ui(ringflaga, characteristic)
            ringflagb = 1
    else:
        raise NotImplementedError("Base ring is not supported.")

    _ring = <ring*>omAlloc0Bin(sip_sring_bin)
    _ring.ch = characteristic
    _ring.ringtype = ringtype
    _ring.N = n
    _ring.names  = _names

    if is_extension:
        rChangeCurrRing(k._ring)
        _ring.algring = rCopy0(k._ring)
        rComplete(_ring.algring, 1)
        _ring.algring.pCompIndex = -1
        _ring.P = _ring.algring.N
        _ring.parameter = <char**>omAlloc0(sizeof(char*)*2)
        _ring.parameter[0] = omStrDup(_ring.algring.names[0])

        nmp = <lnumber*>omAlloc0Bin(rnumber_bin)
        nmp.z= <napoly*>p_Copy(minpoly._poly, _ring.algring) # fragile?
        nmp.s=2

        _ring.minpoly=<number*>nmp

    nblcks = len(order.blocks)
    offset = 0

    _ring.wvhdl  = <int **>omAlloc0((nblcks + 2) * sizeof(int*))
    _ring.order  = <int *>omAlloc0((nblcks + 2) * sizeof(int *))
    _ring.block0 = <int *>omAlloc0((nblcks + 2) * sizeof(int *))
    _ring.block1 = <int *>omAlloc0((nblcks + 2) * sizeof(int *))
    _ring.OrdSgn = 1

    for i from 0 <= i < nblcks:
        _ring.order[i] = order_dict.get(order[i].singular_str(), ringorder_dp)
        _ring.block0[i] = offset + 1
        if len(order[i]) == 0: # may be zero in some cases
            _ring.block1[i] = offset + n
        else:
            _ring.block1[i] = offset + len(order[i])
        offset = _ring.block1[i]

    # TODO: if we construct a free module don't hardcode! This
    # position determines whether we break ties at monomials first or
    # whether we break at indices first!
    _ring.order[nblcks] = ringorder_C

    if ringtype != 0:
        _ring.ringflaga = ringflaga
        _ring.ringflagb = ringflagb

    rComplete(_ring, 1)
    _ring.ShortOut = 0

    rChangeCurrRing(_ring)
    return _ring

cdef void singular_ring_delete(ring *doomed):
    """
    Carefully deallocate the ring, without changing "currRing" (since
    this method can be at unpredictable times due to garbage
    collection).

    TESTS:

    This example caused a segmentation fault with a previous version
    of this method::

        sage: import gc
        sage: from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
        sage: R1 = MPolynomialRing_libsingular(GF(5), 2, ('x', 'y'), TermOrder('degrevlex', 2))
        sage: R2 = MPolynomialRing_libsingular(GF(11), 2, ('x', 'y'), TermOrder('degrevlex', 2))
        sage: R3 = MPolynomialRing_libsingular(GF(13), 2, ('x', 'y'), TermOrder('degrevlex', 2))
        sage: _ = gc.collect()
        sage: foo = R1.gen(0)
        sage: del foo
        sage: del R1
        sage: _ = gc.collect()
        sage: del R2
        sage: _ = gc.collect()
        sage: del R3
        sage: _ = gc.collect()
    """
    cdef ring *oldRing = NULL
    if currRing != doomed:
        oldRing = currRing
        rChangeCurrRing(doomed)
        rDelete(doomed)
        rChangeCurrRing(oldRing)
    else:
        (&currRing)[0] = NULL
        rDelete(doomed)

