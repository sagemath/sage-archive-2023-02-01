"""
Find generators for the ring of modular forms of given level.

AUTHORS:
    -- William Stein (20070824): first version
"""

from sage.rings.all import Integer, QQ, infinity

def span_of_series(v, prec=None, basis=False):
    """
    Return the free module spanned by the given list of power series
    or objects with a padded_list method.  If prec is not given, the
    precision used is the minimum of the precisions of the elements in
    the list (as determined by a prec method).

    INPUT:
        v    -- a list of power series
        prec -- optional; if given then the series do not have to be
                of finite precision, and will be considered to have
                precision prec.
        basis -- (default: False) if True the input is assumed to
                determine linearly independent vectors, and
                the resulting free module has that as basis.


    OUTPUT:
          A free module of rank prec over the base ring of the forms
          (actually, of the first form in the list).  If the list is
          empty, the free module is over QQ.

    EXAMPLES:
    An example involving modular forms:
        sage: v = ModularForms(11,2, prec=5).basis(); v
        [
        q - 2*q^2 - q^3 + 2*q^4 + O(q^5),
        1 + 12/5*q + 36/5*q^2 + 48/5*q^3 + 84/5*q^4 + O(q^5)
        ]
        sage: span_of_series(v)
        Vector space of degree 5 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0 12 12 12]
        [ 0  1 -2 -1  2]

    Next we make sure the vector give a basis:
        sage: span_of_series(v,basis=True)
        Vector space of degree 5 and dimension 2 over Rational Field
        User basis matrix:
        [   0    1   -2   -1    2]
        [   1 12/5 36/5 48/5 84/5]

    An example involving power series.
        sage: R.<x> = PowerSeriesRing(QQ, default_prec=5)
        sage: v = [1/(1-x), 1/(1+x), 2/(1+x), 2/(1-x)]; v
        [1 + x + x^2 + x^3 + x^4 + O(x^5),
         1 - x + x^2 - x^3 + x^4 + O(x^5),
         2 - 2*x + 2*x^2 - 2*x^3 + 2*x^4 + O(x^5),
         2 + 2*x + 2*x^2 + 2*x^3 + 2*x^4 + O(x^5)]
        sage: span_of_series(v)
        Vector space of degree 5 and dimension 2 over Rational Field
        Basis matrix:
        [1 0 1 0 1]
        [0 1 0 1 0]
        sage: span_of_series(v,10)
        Vector space of degree 10 and dimension 2 over Rational Field
        Basis matrix:
        [1 0 1 0 1 0 0 0 0 0]
        [0 1 0 1 0 0 0 0 0 0]

    An example involving polynomials.
        sage: x = polygen(QQ)
        sage: span_of_series([x^3, 2*x^2 + 17*x^3, x^2])
        Traceback (most recent call last):
        ...
        ValueError: please specify a precision
        sage: span_of_series([x^3, 2*x^2 + 17*x^3, x^2],5)
        Vector space of degree 5 and dimension 2 over Rational Field
        Basis matrix:
        [0 0 1 0 0]
        [0 0 0 1 0]
        sage: span_of_series([x^3, 2*x^2 + 17*x^3, x^2],3)
        Vector space of degree 3 and dimension 1 over Rational Field
        Basis matrix:
        [0 0 1]
    """
    if len(v) == 0:
        if not prec:
            prec = 0
        return (QQ**prec).zero_submodule()
    if prec:
        n = Integer(prec)
    else:
        n = min([g.prec() for g in v])
        if n == infinity:
            raise ValueError, "please specify a precision"

    K = v[0].parent().base_ring()
    V = K**n
    B = [V(g.padded_list(n)) for g in v]
    if basis:
        M = V.span_of_basis(B)
    else:
        M = V.span(B)
    return M
