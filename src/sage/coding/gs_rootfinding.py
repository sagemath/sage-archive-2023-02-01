from sage.rings.infinity import infinity
from sage.functions.other import binomial

def _convert_Q_representation(Q):
    """Convert Q from either F[x,y], F[x][y] or F[x] list to F[x] list"""
    if isinstance(Q, list):
        Rx = Q[0].parent()
        if not hasattr(Rx,'gen'):
            raise ValueError("Q must be given as F[x][y], F[x,y] or as F[x] list.")
        return Q
    else:
        # Find out if Q is in F[x,y] or F[x][y]
        Qorig = Q
        Ryx = Q.parent()
        #TODO: Check Ryx is a polynomial ring over a field
        if len(Ryx.gens())==1:
            # Ok, Q is in F[x][y]
            pass
        elif len(Ryx.gens())==2:
            F = Ryx.base_ring()
            #(xs,ys) = Ryx.variable_names()
            #Rx.<x> = F[xs]
            #Ryx.<y> = Rx[ys]
            Rx = F['xs']
            Ryx = Rx['ys']
            x, y = Rx.gen(), Ryx.gen()
            #TODO: This conversion is wasteful since we will convert to F[x] list immediately after
            Q = Ryx(Q)
        else:
            raise ValueError("Q must be given as F[x][y], F[x,y] or as F[x] list.")
        # Then make sure Q is a list of F[x] elements
        return Q.list()

def _sanitise_rootfinding_input(Q, maxd, precision):
    """Local method for doing basic preparatory manipulations to rootfinding
    methods."""
    Qinp = Q
    Rx = Q[0].parent()
    F = Rx.base_ring()
    x = Rx.gen()

    if all(p.is_zero() for p in Q):
        if precision:
            return [(Rx.zero(), 0)]
        else:
            return ValueError("The zero polynomial has infinitely many roots.")
    if not maxd:
        if precision:
            maxd = precision-1
        else:
            # The maximal degree of a root is at most the degree of the first
            # non-zero coeff
            for p in Q:
                if not p.is_zero():
                    maxd = p.degree()
                    break
    if precision:
        for t in range(len(Q)):
            if Q[t].degree >= precision:
                Q[t] = Q[t].truncate(precision)
    return (Q, Qinp, F, Rx, x, maxd)

def _strip_x_pows(Q):
    r"""Take `Q \in F[x]^ell` and divide all elements by the largest possible
    power s of x such that all elements remain polynomials.
    Return (Q', s) where Q' is the new polynomial.
    """
    def lead_zeroes(p):
        if p.is_zero():
            return infinity
        i = 0
        while p[i].is_zero():
            i+=1
        return i
    strip = min([lead_zeroes(p) for p in Q])
    if strip == 0:
        return (Q, 0)
    if strip == infinity:
        return ([ Q[0].parent().zero() ], infinity)
    return ([ p.shift(-strip) for p in Q ] , strip)

def _roth_ruckenstein_i(Q, F, Rx, x, maxd, precision):
    r"""Internal function: The core of Roth-Ruckenstein's algorithm where all
    conversion, checks and parent-extraction, has been done."""
    sols = []
    g = [F.zero()] * (maxd+1)
    def roth_rec(Q, lam, k):
        """Recursion of the root finding:
        Q is the remaining poly, lam is the power of x whose coefficient we are
        to determine now, and k is the remaining precision to handle (if 'precision' is given)"""
        if precision and k <= 0:
            sols.append((Rx(g[:lam]), lam))
            return
        (T, strip) = _strip_x_pows(Q)
        #print "  "*lam, "lam:%s , Q:%s , T:%s , strip:%s, k:%s, \t\t g: %s" % (lam, Q, T, strip, k, g[:lam])
        if precision:
            k = k - strip
        Ty = Rx([ p[0] for p in T ]) # Ty = T(0,y)
        if Ty.is_zero() or (precision and k <= 0):
            if precision:
                sols.append((Rx(g[:lam]), lam))
            else:
                assert all(p.is_zero() for p in Q) , ("Q is not zero but Ty is?:\nQ = %s" % Q)
                sols.append(Rx(g[:lam]))
            return
        roots = Ty.roots(multiplicities=False)
        for gamma in roots:
            #print " "*lam, "Choosing root %s" % gamma
            g[lam] = gamma
            if lam<maxd:
                # Construct T(y=x*y + gamma)
                ell = len(T)-1
                yc = [ [ binomial(s, t)*x**t*gamma**(s-t) for t in range(s+1) ] for s in range(ell+1) ]#pd gamma, s-t
                Tg = []
                for t in range(ell+1):
                    Tg.append(sum( yc[s][t]*T[s] for s in range(t, ell+1)))
                roth_rec(Tg , lam+1, k)
            else:
                if precision:
                    sols.append((Rx(g[:lam+1]), lam+1))
                elif sum( Q[t] * gamma**t for t in range(len(Q)) ).is_zero():#pd gamma, t
                    sols.append(Rx(g[:lam+1]))
        return
    roth_rec(Q, 0, precision)
    return sols

def rootfind_roth_ruckenstein(Q, maxd=None, precision=None):
    r"""Use the Roth-Ruckenstein algorithm to find roots or roots
    modulo-up-to-some-precision of a `Q \in F[x][y]` where `F` is a field.

    Q should be given as a list of F[x] elements.

    If 'precision=None' then actual roots will be found, i.e. all `f \in \F[x]`
    such that `Q(f) = 0`. This will be returned as a list of `\F[x]` elements.

    If 'precision=k' for some integer `k`, then all `f \in \F[x]` such that
    `Q(f) \equiv 0 \mod x^k` will be returned. This set is infinite, and so it
    will be returned as a list of pairs in `\F[x] \times \mathbb{Z}_+`, where
    `(f, d)` denotes that `Q(f + x^d h) \equiv 0 \mod x^k` for any `h \in
    \F[[x]]`.

    If `maxd` is given, then find only `f` with `deg f \leq maxd`. In case
    `precision=k` setting `maxd` means to only find the roots up to precision
    `maxd`; otherwise, the precision will be `precision-1`.
    """
    (Q, Qinp, F, Rx, x, maxd) = _sanitise_rootfinding_input(Q, maxd, precision)
    return _roth_ruckenstein_i(Q, F, Rx, x, maxd, precision)

def rootfind_bivariate(Q, maxd=None, algorithm=None, **kwargs):
    """Find y-roots of the bivariate polynomial $Q \in F[x][y]$.
    If maxd is given, find only roots satisfying $\\deg f(x) <= maxd$. This has
    no effect when Alekhnovich's algorithm is used, however.

    Algorithm can be "roth_ruckenstein" or "alekhnovich". If not given, the
    former will be used for small polynomials while the latter will be used for
    bigger ones.

    TESTS:
        F = GF(11)
        Px.<x> = F[]
        P.<y> = Px[]
        Q = (y - x^2 + x + 1) * (y^2 - 1) * (y - x^3 + 4*x - 1)
        rootfind_bivariate(Q)
    ... [10, x^2 + 10*x + 10, x^3 + 7*x + 1]
    """
    Q = _convert_Q_representation(Q)
    if algorithm is None:
        xdeg = max( Qi.degree() for Qi in Q )
        if xdeg < THRESHOLD_ALEKHNOVICH:
            algorithm = "roth_ruckenstein"
        else:
            algorithm = "alekhnovich"
    if algorithm=="roth_ruckenstein":
        return rootfind_roth_ruckenstein(Q, maxd=maxd, **kwargs)
    elif algorithm=="alekhnovich":
        return rootfind_alekhnovich(Q, maxd=maxd, **kwargs)
    else:
        raise ValueError("No such root-finding algorithm: %s" % algorithm)

