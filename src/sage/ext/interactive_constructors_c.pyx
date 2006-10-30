# Interactive versions.

import sage.rings.all
SAGE_GLOBAL=None

def PolynomialRing(*args, **kwds):
    kwds['scope'] = True
    return sage.rings.all.PolynomialRing(*args, **kwds)



