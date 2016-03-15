from sage.misc.superseded import deprecation
deprecation(20135, "the module sage.gsl.callback is deprecated")

from sage.libs.gsl.all cimport gsl_function, gsl_diff_central


foo = None
cdef double f(double x, void* p):
    return foo(x)

def diff(g, double x=1.0):
    cdef double result, abserr
    cdef gsl_function F

    global foo
    foo = g
    F.function = f
    F.params = NULL

    gsl_diff_central(&F, x, &result, &abserr)

    return result, abserr
