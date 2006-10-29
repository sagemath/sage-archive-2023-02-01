#import interactive_constructors_c
import sage.rings.all

SAGE_GLOBAL=None

def PolynomialRing(*args, **kwds):
    """
    Construct a polynomial ring and assign the variables of the polynomial
    ring to the global interactive interpreter.

    For the source code of this function, do
        sage.: import sage.rings.polynomial_ring_c
        sage.: sage.rings.polynomial_ring_c.PolynomialRing??
    -------------------------------------------------------------
    """
    kwds['scope'] = SAGE_GLOBAL
    return sage.rings.all.PolynomialRing(*args, **kwds)
    #return interactive_constructors_c.PolynomialRing(*args, **kwds)

PolynomialRing.__doc__ += sage.rings.all.PolynomialRing.__doc__
