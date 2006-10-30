import sage.rings.all

SAGE_GLOBAL=None

def inject(X):
    X.inject_variables(SAGE_GLOBAL)
    return X

def PolynomialRing(*args, **kwds):
    """
    Construct a polynomial ring and assign the variables of the polynomial
    ring to the global interactive interpreter.

    For the source code of the function that actualy creates the polynomial
    ring (without assigning variables), do:
        sage.: sage.rings.PolynomialRing??
    -------------------------------------------------------------
    """
    return inject(sage.rings.all.PolynomialRing(*args, **kwds))
PolynomialRing.__doc__ += sage.rings.all.PolynomialRing.__doc__



