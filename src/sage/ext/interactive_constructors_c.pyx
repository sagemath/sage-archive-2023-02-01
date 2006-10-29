# Interactive versions.

import sage.rings.all
def PolynomialRing(*args, **kwds):
    kwds['inject_variables'] = True
    return sage.rings.all.PolynomialRing(*args, **kwds)



