r"""
The Elliptic Curve Factorization Method, via libecm.

\sage includes GMP-ECM, which is a highly optimized implementation
of Lenstra's elliptic curve factorization method.  See
\url{http://ecm.gforge.inria.fr/}
for more about GMP-ECM.

AUTHOR:
    Robert L Miller (2008-01-21) -- library interface (clone of ecmfactor.c)

EXAMPLE:

    sage: import sage.interfaces.libecm
    sage: from sage.interfaces.libecm import ecmfactor
    sage: result = ecmfactor(999, 0.00)
    sage: result in [(True, 27), (True, 37), (True, 999)]
    True
    sage: result = ecmfactor(999, 0.00, verbose=True)
    Performing one curve with B1=0
    Found factor in step 1: ...
    sage: result in [(True, 27), (True, 37), (True, 999)]
    True

"""

include '../ext/cdefs.pxi'

from sage.rings.integer cimport Integer

cdef extern from "ecm.h":
    ctypedef struct __ecm_param_struct:
        pass
    ctypedef __ecm_param_struct ecm_params[1]
    int ecm_factor (mpz_t, mpz_t, double, ecm_params)
    int ECM_NO_FACTOR_FOUND

def ecmfactor(number, py_B1, verbose=False):
    cdef mpz_t n, f
    cdef int res
    cdef double B1
    cdef Integer sage_int_f, sage_int_number

    sage_int_f = Integer(0)
    sage_int_number = Integer(number)
    B1 = py_B1

    mpz_init(n)
    mpz_set(n, sage_int_number.value)
    mpz_init(f) # For potential factor

    if verbose:
        print "Performing one curve with B1=%1.0f"%B1

    res = ecm_factor(f, n, B1, NULL)

    if res > 0:
        mpz_set(sage_int_f.value, f)

    mpz_clear(f)
    mpz_clear(n)

    if res > 0:
        if verbose:
            print "Found factor in step %d: %d"%(res,sage_int_f)
        return (True, sage_int_f)
    elif res == ECM_NO_FACTOR_FOUND:
        if verbose:
            print "Found no factor."
        return (False, None)
    else:
        raise RuntimeError( "ECM lib error" )

