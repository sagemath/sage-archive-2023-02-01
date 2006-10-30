# Interactive versions.

cdef inject(X, do):
    if do:
        X.inject_variables()
    return X

cdef do_inject(kwds):
    if kwds.has_key('inject') and kwds['inject'] == False:
        del kwds['inject']
        return False
    return True

def PolynomialRing(*args, **kwds):
    """
    Construct a polynomial ring and injects the variables of the
    polynomial ring to the global interactive interpreter.  Use
    inject=False to not inject the variables.  This is a wrapper
    around the following function: <<<PolynomialRing>>>

    MORE EXAMPLES:
    We illustrate creating a polynomial ring without injecting the variables
    into the interpreter.

        sage: PolynomialRing(QQ,'w')
        Univariate Polynomial Ring in w over Rational Field
        sage: parent(w)
        Univariate Polynomial Ring in w over Rational Field
        sage: PolynomialRing(GF(17), 'w', inject=False)
        Univariate Polynomial Ring in w over Finite Field of size 17
        sage: parent(w)
        Univariate Polynomial Ring in w over Rational Field
    """
    t = do_inject(kwds)
    import sage.rings.all
    R = sage.rings.all.PolynomialRing(*args, **kwds)
    return inject(R, t)

def FreeMonoid(*args, **kwds):
    """
    Construct a free monoid and inject the variables of the monoid
    into the global interactive interpreter.  Use inject=Fale to not
    inject the variables.  This is a wrapper around the following
    function: <<<FreeMonoid>>>

    MORE EXAMPLES:
    We illustrate creating a free monoid with and without injecting
    the variables into the interpreter.

        sage: FreeMonoid(4,'x')
        Free monoid on 4 generators (x0, x1, x2, x3)
        sage: x2
        x2
        sage: FreeMonoid(4,'y', inject=False)
        Free monoid on 4 generators (y0, y1, y2, y3)
        sage: y
        Traceback (most recent call last):
        ...
        NameError: name 'y' is not defined
    """
    t = do_inject(kwds)
    import sage.monoids.free_monoid
    R = sage.monoids.free_monoid.FreeMonoid(*args, **kwds)
    return inject(R, t)





###################### need to add a bunch more ############################
