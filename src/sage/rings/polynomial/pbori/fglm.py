if __name__ == "__main__":
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))


    def _test():
        import doctest
        doctest.testmod()

from .PyPolyBoRi import BooleSet, Polynomial, BoolePolynomialVector, \
    FGLMStrategy, Monomial, Ring

from .blocks import declare_ring


def _fglm(I, from_ring, to_ring):
    r"""
    Unchecked variant of fglm
    """
    vec = BoolePolynomialVector(I)
    return FGLMStrategy(from_ring, to_ring, vec).main()


def fglm(I, from_ring, to_ring):
    r"""
    Converts *reduced* Groebner Basis in from_ring to a GroebnerBasis in to_ring.
    It acts independend of the global ring, which is restored at the end of the
    computation.
    
    TESTS::
    
        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import OrderCode
        sage: dp_asc = OrderCode.dp_asc
        sage: r=declare_ring(['x','y','z'],dict())
        sage: old_ring = r
        sage: new_ring = old_ring.clone(ordering=dp_asc)
        sage: (x,y,z) = [old_ring.variable(i) for i in xrange(3)]
        sage: ideal=[x+z, y+z]# lp Groebner basis
        sage: list(fglm(ideal, old_ring, new_ring))
        [y + x, z + x]
    """
    for poly in I:
        if poly.ring().id() != from_ring.id():
            raise ValueError("Ideal I must be from the first ring argument")
    return _fglm(I, from_ring, to_ring)


def vars_real_divisors(monomial, monomial_set):
    r"""
    Returns all elements of of monomial_set, which result multiplied by a variable in monomial.
    
    TESTS::
    
        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import OrderCode
        sage: dp_asc = OrderCode.dp_asc
        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import Ring
        sage: r=Ring(1000)
        sage: x = r.variable
        sage: b=BooleSet([x(1)*x(2),x(2)])
        sage: vars_real_divisors(x(1)*x(2)*x(3),b)
        {{x(1),x(2)}}
    """
    return BooleSet(Polynomial(monomial_set.divisors_of(monomial)). \
        graded_part(monomial.deg() - 1))


def m_k_plus_one(completed_elements, variables):
    r""" 
    Calculates $m_{k+1}$ from the FGLM algorithm as described in Wichmanns diploma thesis
    It would be nice to be able to efficiently extract the smallest term of a polynomial.
    
    TESTS::
    
        sage: from sage.rings.polynomial.pbori.PyPolyBoRi import OrderCode
        sage: dp_asc = OrderCode.dp_asc
        sage: r=Ring(1000)
        sage: x = r.variable
        sage: s=BooleSet([x(1)*x(2),x(1),x(2),Monomial(r),x(3)])
        sage: variables=BooleSet([x(1),x(2),x(3)])
        sage: m_k_plus_one(s,variables)
        x(2)*x(3)
        sage: r2 = r.clone(ordering=dp_asc)
        sage: m_k_plus_one(r2(s).set(),r2(variables).set())
        x(1)*x(3)
    """
    return sorted(completed_elements.cartesian_product(variables).diff(
        completed_elements))[0]


if __name__ == "__main__":
    _test()
