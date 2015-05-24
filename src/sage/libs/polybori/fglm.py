if __name__ == "__main__":
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))


    def _test():
        import doctest
        doctest.testmod()

from sage.libs.polybori.PyPolyBoRi import BooleSet, Polynomial, BoolePolynomialVector, \
    FGLMStrategy, Monomial, Ring

from sage.libs.polybori.blocks import declare_ring


def _fglm(I, from_ring, to_ring):
    """Unchecked variant of fglm"""
    vec = BoolePolynomialVector(I)
    return FGLMStrategy(from_ring, to_ring, vec).main()


def fglm(I, from_ring, to_ring):
    """
    converts *reduced* Groebner Basis in from_ring to a GroebnerBasis in to_ring.
    It acts independend of the global ring, which is restored at the end of the
    computation,
    >>> from sage.libs.polybori import OrderCode, declare_ring
    >>> from sage.libs.polybori.fglm import fglm
    >>> dp_asc = OrderCode.dp_asc
    >>> r=declare_ring(['x','y','z'],dict())
    >>> old_ring = r
    >>> new_ring = old_ring.clone(ordering=dp_asc)
    >>> (x,y,z) = [old_ring.variable(i) for i in xrange(3)]
    >>> ideal=[x+z, y+z]# lp Groebner basis
    >>> list(fglm(ideal, old_ring, new_ring))
    [y + x, z + x]
    """
    for poly in I:
        if poly.ring().id() != from_ring.id():
            raise ValueError, "Ideal I must be from the first ring argument"
    return _fglm(I, from_ring, to_ring)


def vars_real_divisors(monomial, monomial_set):
    """
    returns all elements of of monomial_set, which result multiplied by a variable in monomial.
    >>> from sage.libs.polybori import OrderCode, Ring, BooleSet
    >>> from sage.libs.polybori.fglm import vars_real_divisors
    >>> dp_asc = OrderCode.dp_asc
    >>> r=Ring(1000)
    >>> x = r.variable
    >>> b=BooleSet([x(1)*x(2),x(2)])
    >>> vars_real_divisors(x(1)*x(2)*x(3),b)
    {{x(1),x(2)}}
    """
    return BooleSet(Polynomial(monomial_set.divisors_of(monomial)). \
        graded_part(monomial.deg() - 1))


def m_k_plus_one(completed_elements, variables):
    """ calculates $m_{k+1}$ from the FGLM algorithm as described in Wichmanns diploma thesis
    It would be nice to be able to efficiently extract the smallest term of a polynomial
    >>> from sage.libs.polybori import OrderCode, Ring, BooleSet, Monomial
    >>> from sage.libs.polybori.fglm import m_k_plus_one
    >>> dp_asc = OrderCode.dp_asc
    >>> r=Ring(1000)
    >>> x = r.variable
    >>> s=BooleSet([x(1)*x(2),x(1),x(2),Monomial(r),x(3)])
    >>> variables=BooleSet([x(1),x(2),x(3)])
    >>> m_k_plus_one(s,variables)
    x(2)*x(3)
    >>> r2 = r.clone(ordering=dp_asc)
    >>> m_k_plus_one(r2(s).set(),r2(variables).set())
    x(1)*x(3)
    """
    return sorted(completed_elements.cartesian_product(variables).diff(
        completed_elements))[0]


if __name__ == "__main__":
    _test()
