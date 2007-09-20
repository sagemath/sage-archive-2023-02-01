
def _monomials(gens, R, n, i):
    # each power of the ith generator times all products
    # not involving the ith generator.
    if len(gens) == 1:
        b = gens[0]
        v = [R(1)]
        for _ in range(n-1):
            v.append(v[-1]*b)
        return v
    else:
        z = gens[i]
        w = list(gens)
        del w[i]
        v = monomials(w, n)
        k = len(v)
        for _ in range(n-1):
            for j in range(k):
                v.append(v[j]*z)
            z *= gens[i]
        return v

from sage.structure.sequence import Sequence

def monomials(v, n):
    """
    Given a list v of numbers and an integer n, return
    are monomials in the elements of v, where each variable
    in the monomial appears to degree less than n.

    EXAMPLES:
        sage: monomials([x], 3)
        [1, x, x^2]
        sage: R.<x,y,z> = QQ[]
        sage: monomials([x,y], 5)
        [1, y, y^2, y^3, y^4, x, x*y, x*y^2, x*y^3, x*y^4, x^2, x^2*y, x^2*y^2, x^2*y^3, x^2*y^4, x^3, x^3*y, x^3*y^2, x^3*y^3, x^3*y^4, x^4, x^4*y, x^4*y^2, x^4*y^3, x^4*y^4]
    """
    v = Sequence(v)
    R = v.universe()
    return _monomials(v, R, n, 0)
