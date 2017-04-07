r"""
Cython wrapper for sirocco library

This is used to call the sirocco library directly from python.

The function contpath takes the following input:

- An integer, representing the degree of the polynomial

- A list of floating point numbers. Each four consecutive elements 
of this list represent the interval corresponding to a coefficient.
Coefficients are listed in increasing deglex order, and inside each
coefficients, the four numbers represent the lower real, upper real,
lower imaginary and real imaginary limits of the interval.

- A float representing the real part of the initial root approximation

- A float representing the imaginary part of the initial root.

The output is a list of tuples. Each tuple represents the x value (between 0 and 1)
and the real and imaginary parts of the y value of a vertex in the piecewise linear
approximation of the path tracked by the root.

AUTHORS:

- Miguel Marco (2016-07-19): initial version.
"""


include 'cysignals/signals.pxi'
include "sage/ext/stdsage.pxi" 

from sage.libs.mpfr cimport *
from sage.rings.real_mpfr cimport RealNumber
from sage.rings.real_mpfr import RealField


cdef extern from "stdlib.h":
    void free(void* ptr)


cdef extern from "sirocco.h":
    mpfr_t* homotopyPath_mp (int degree, mpfr_t *_coef, mpfr_t _y0R, mpfr_t _y0I, int prec)
    double * homotopyPath (int degree, double *_coef, double _y0R, double _y0I)



def devuelve(RealNumber a):
    cdef mpfr_prec_t ap = mpfr_get_prec(a.value)
    cdef RealNumber my_number = RealField(ap)()
    mpfr_set(my_number.value, a.value, MPFR_RNDN)
    return my_number


def contpath_mp(deg, values, RealNumber y0r, RealNumber y0i, int prec):
    cdef int cdeg = deg
    cdef mpfr_t* cvalues = <mpfr_t*> sage_malloc(sizeof(mpfr_t)*len(values))
    cdef mpfr_t* rop
    cdef int j
    for i in range(len(values)):
        j = <int>i
        sig_on()
        mpfr_init2(cvalues[j], prec)
        mpfr_set(cvalues[j], (<RealNumber>values[j]).value, MPFR_RNDN)
        sig_off()
    cdef mpfr_t y0R
    cdef mpfr_t y0I 
    sig_on()
    mpfr_init2(y0R, prec)
    mpfr_set(y0R, (<RealNumber>y0r).value, MPFR_RNDN)
    mpfr_init2(y0I, prec)
    mpfr_set(y0I, (<RealNumber>y0i).value, MPFR_RNDN)
    rop = homotopyPath_mp(cdeg, cvalues, y0R, y0I, prec)
    sig_off()
    for i in range(len(values)):
        j = <int>i
        sig_on()
        mpfr_clear(cvalues[j])
        sig_off()
    free(cvalues)
    if rop == NULL:
        raise ValueError("libsirocco could not guarantee one step")
    cdef int n = mpfr_get_si(rop[0], MPFR_RNDN)
    
    #cdef RealNumber res = RealField(53)()
    cdef double res
    l = []
    
    for i in range(1, 3*n+1):
        RN = RealField(prec)()
        sig_on()
        res = mpfr_set((<RealNumber>RN).value, rop[i], MPFR_RNDN)
        mpfr_clear(rop[i])
        sig_off()
        l.append(RN)
    free(rop)
    return [(l[3*i], l[3*i+1], l[3*i+2]) for i in range(len(l)/3)]
    



def contpath(deg,values,y0r,y0i):
    cdef double* rop
    cdef double* c_values = <double*> sage_malloc(sizeof(double)*len(values))
    cdef int clen = <int> len(values)
    for i,v in enumerate(values):
        c_values[i] = values[i]
        
        
    
    cdef double y0R = y0r
    cdef double y0I = y0i
    sig_on()
    rop = homotopyPath (int(deg), c_values, y0R, y0I)
    sig_off()
    if rop == NULL:
        raise ValueError("libsirocco could not guarantee one step")
    n=int(rop[0])
    l=[0 for i in range(n)]
    for i in range(n):
        l[i]=(rop[3*i+1],rop[3*i+2],rop[3*i+3])
    free(rop)
    free(c_values)
    return l
