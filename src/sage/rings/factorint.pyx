r"""
Integer factorization functions

AUTHORS:

- Andre Apitzsch (2011-01-13): initial version

"""

#*****************************************************************************
#       Copyright (C) 2010 André Apitzsch <andre.apitzsch@st.ovgu.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../libs/pari/decl.pxi"
include "../ext/gmp.pxi"
include "../ext/stdsage.pxi"

from sage.rings.integer cimport Integer
from sage.structure.factorization_integer import IntegerFactorization

cdef extern from "limits.h":
    long LONG_MAX

def factor_cunningham(m, proof=None):
    r"""
    Return factorization of self obtained using trial division
    for all primes in the so called Cunningham table. This is
    efficient if self has some factors of type $b^n+1$ or $b^n-1$,
    with $b$ in $\{2,3,5,6,7,10,11,12\}$.

    You need to install an optional package to use this method,
    this can be done with the following command line:
    ``sage -i cunningham_tables-1.0``

    INPUT:

     - ``proof`` - bool (default: None) whether or not to
       prove primality of each factor, this is only for factors
       not in the Cunningham table.

    EXAMPLES::

        sage: (2^257-1)._factor_cunningham() # optional - cunningham
        535006138814359 * 1155685395246619182673033 * 374550598501810936581776630096313181393
        sage: ((3^101+1)*(2^60).next_prime())._factor_cunningham(proof=False) # optional - cunningham
        2^2 * 379963 * 1152921504606847009 * 1017291527198723292208309354658785077827527

    """
    from sage.databases import cunningham_tables
    cunningham_prime_factors = cunningham_tables.cunningham_prime_factors()
    if m.nbits() < 100 or len(cunningham_prime_factors) == 0:
        return m.factor(proof=proof)
    n = Integer(m)
    L = []
    for p in cunningham_prime_factors:
        if p > n:
            break
        if p.divides(n):
            v,n = n.val_unit(p)
            L.append( (p,v) )
    if n.is_one():
        return IntegerFactorization(L)
    else:
        return IntegerFactorization(L)*n.factor(proof=proof)

cpdef factor_trial_division(m, long limit=LONG_MAX):
    r"""
    Return partial factorization of self obtained using trial division
    for all primes up to limit, where limit must fit in a C signed long.

    INPUT:

     - ``limit`` - integer (default: LONG_MAX) that fits in a C signed long

    EXAMPLES::

        sage: from sage.rings.factorint import factor_trial_division
        sage: n = 920384092842390423848290348203948092384082349082
        sage: factor_trial_division(n, 1000)
        2 * 11 * 41835640583745019265831379463815822381094652231
        sage: factor_trial_division(n, 2000)
        2 * 11 * 1531 * 27325696005058797691594630609938486205809701
    """
    cdef Integer n = PY_NEW(Integer), unit = PY_NEW(Integer), p
    cdef unsigned long e

    n = Integer(m)
    if mpz_sgn(n.value) > 0:
        mpz_set_si(unit.value, 1)
    else:
        mpz_neg(n.value,n.value)
        mpz_set_si(unit.value, -1)

    F = []
    while mpz_cmpabs_ui(n.value, 1):
        p = n.trial_division(bound=limit)
        e = mpz_remove(n.value, n.value, p.value)
        F.append((p,e))

    return IntegerFactorization(F, unit=unit, unsafe=True,
                                   sort=False, simplify=False)

cpdef factor_using_pari(n, int_=False, debug_level=0, proof=None):
    r"""
    Factors this (positive) integer using PARI.

    This method returns a list of pairs, not a ``Factorization``
    object. The first element of each pair is the factor, of type
    ``Integer`` if ``int_`` is ``False`` or ``int`` otherwise,
    the second element is the positive exponent, of type ``int``.

    INPUT:

     - ``int_`` -- (default: ``False``), whether the factors are
       of type ``int`` instead of ``Integer``

     - ``debug_level`` -- (default: 0), debug level of the call
       to PARI

     - ``proof`` -- (default: ``None``), whether the factors are
       required to be proven prime;  if ``None``, the global default
       is used

    OUTPUT:

        -  a list of pairs

    EXAMPLES::

        sage: factor(-2**72 + 3, algorithm='pari')  # indirect doctest
        -1 * 83 * 131 * 294971519 * 1472414939
    """
    from sage.libs.pari.gen import pari

    if proof is None:
        from sage.structure.proof.proof import get_flag
        proof = get_flag(proof, "arithmetic")

    prev = pari.get_debug_level()

    if prev != debug_level:
        pari.set_debug_level(debug_level)

    F = pari(n).factor(proof=proof)
    B, e = F
    if int_:
        v = [(int(B[i]), int(e[i])) for i in xrange(len(B))]
    else:
        v = [(Integer(B[i]), int(e[i])) for i in xrange(len(B))]

    if prev != debug_level:
        pari.set_debug_level(prev)

    return v
