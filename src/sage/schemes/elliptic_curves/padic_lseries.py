from sage.rings.integer_ring import ZZ

def find_alpha(E, p, prec=20):
    r"""
    Return a p-adic root $\alpha$ of the polynomial $x^2 - a_p x + p$
    with $\ord_p(\alpha) < 1$.
    """
    a_p = E.ap(p)

    R = ZZ['x']
    f = R([p, -a_p, 1])
    G = f.factor_padic(p, prec)
    for pr, e in G:
        a = -pr[0]
        if a.valuation() < 1:
            return a
    raise ValueError, "no root works"



