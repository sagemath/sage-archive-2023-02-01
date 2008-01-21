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
    sage: from sage.interfaces.libecm import py_ecm_factor
    sage: py_ecm_factor(999, 0.00)
    Performing one curve with B1=0
    Found factor in step 1:
    27



"""

include '../ext/cdefs.pxi'

from sage.rings.integer cimport Integer
import sys

#cdef extern from "gmp.h":
#    pass

cdef extern from "ecm.h":
    ctypedef struct __ecm_param_struct:
        pass
    ctypedef __ecm_param_struct ecm_params[1]
    int ecm_factor (mpz_t, mpz_t, double, ecm_params)
    int ECM_NO_FACTOR_FOUND

def py_ecm_factor(n, B1):
    ecmfactor(Integer(n), B1)

cdef void ecmfactor(Integer number, py_B1):
    cdef mpz_t n, f
    cdef int res
    cdef double B1
    cdef Integer sage_int_f

    sage_int_f = Integer(0)

    mpz_init(n)

    mpz_set(n, number.value)

    B1 = py_B1

    mpz_init(f) # For potential factor

    print "Performing one curve with B1=%1.0f"%B1

    res = ecm_factor(f, n, B1, NULL)

    if res > 0:
        print "Found factor in step %d: "%res
        mpz_set(sage_int_f.value, f)
        print sage_int_f
        #mpz_out_str(sys.stdout, 10, f)
        print "\n"
    elif res == ECM_NO_FACTOR_FOUND:
        print "Found no factor."
    else:
        print "ECM lib error"
    mpz_clear(f)
    mpz_clear(n)
