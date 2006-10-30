# Interactive versions.

import sage.rings.all

def inject(X):
    X.inject_variables()
    return X

def PolynomialRing(*args, **kwds):
    """
    Construct a polynomial ring and *assign* the variables of the
    polynomial ring to the global interactive interpreter.  This is a
    wrapper around the following function:
    <<<PolynomialRing>>>
    """
    R = sage.rings.all.PolynomialRing(*args, **kwds)
    return inject(R)





