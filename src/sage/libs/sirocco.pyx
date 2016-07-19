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
include 'cysignals/memory.pxi'


cdef extern from "stdlib.h":
    void free(void* ptr)


cdef extern from "sirocco.h":
    double* homotopyPath (int degree, double *_coef, double _y0R, double _y0I)



def contpath(deg,values,y0r,y0i):
    """
    Calls sirocco to perform a root continuation.
    
    EXAMPLES::
    
        sage: from sage.libs.sirocco import contpath
        sage: d = 2
        sage: l = [-2,-2,0,0, 0,0,0,0, 0,0,0,0, 1,1,0,0, 0,0,0,0, 1,1,0,0]
        sage: contpath(d, l, 1.414,0) # optional - sirocco
        [(0.0, 1.414, 0.0),
        (0.499924494298889, 1.3229041915422297, 0.0),
        (0.937443369298889, 1.0588672859983665, 0.0),
        (1.0, 1.0, 0.0)]

    """
    cdef double* rop
    cdef double* c_values = <double*> sig_malloc(sizeof(double)*len(values))
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