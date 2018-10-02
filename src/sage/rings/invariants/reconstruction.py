r"""
Reconstruction of Algebraic Forms

This module reconstructs algebraic forms from the values of their
invariants. Given a set of (classical) invariants, it returns a
form that attains this values as invariants (up to scaling).

AUTHORS:

- Jesper Noordsij (2018-06): initial version
"""

#*****************************************************************************
#     Copyright (C) 2018 Jesper Noordsij <jesper.noordsij@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ
from sage.rings.fraction_field import FractionField


def binary_form_from_invariants(degree, invariants):
    """
    Reconstruct a binary form from the values of its invariants.

    INPUT:

    - ``degree`` --  The degree of the binary form.

    - ``invariants`` --  The values of the invariants of the binary form.

    OUTPUT:

    A set of coefficients of a binary form, whose invariants are equal
    to the given ``invariants`` up to a scaling.

    EXAMPLES::

        sage: invariants = [1, 0, 0]
        sage: binary_form_from_invariants(5, invariants)
        (1, 0, 0, 0, 0, 1)

        sage: binary_form_from_invariants(6, invariants)
        Traceback (most recent call last):
        ...
        NotImplementedError: No reconstruction for binary forms of degree 6 implemented.
    """
    if degree == 2:
        return binary_quadratic_from_invariants(invariants)
    if degree == 3:
        return binary_cubic_from_invariants(invariants)
    if degree == 5:
        return binary_quintic_from_invariants(invariants)
    else:
        raise NotImplementedError('No reconstruction for binary forms of \
                                    degree %s implemented.' % degree)


######################################################################

def binary_quadratic_from_invariants(discriminant):
    """
    Reconstruct a binary quadratic from the value of its discriminant.

    INPUT:

    - ``discriminant`` --  The values of the the discriminant of the
      binary quadratic.

    OUTPUT:

    A set of coefficients of a binary quadratic, whose discriminant
    is equal to the given ``discriminant`` up to a scaling.

    EXAMPLES::

        sage: from sage.rings.invariants.reconstruction import binary_quadratic_from_invariants
        sage: binary_quadratic_from_invariants(1)
        (0, 0, 1)
    """
    if discriminant == 0:
        return (1, 0, 0) # x^2
    else:
        return (0, 0, 1) # x*y

def binary_cubic_from_invariants(discriminant):
    """
    Reconstruct a binary cubic from the value of its discriminant.

    INPUT:

    - ``discriminant`` --  The values of the the discriminant of the
      binary cubic.

    OUTPUT:

    A set of coefficients of a binary cubic, whose discriminant
    is equal to the given ``discriminant`` up to a scaling.

    EXAMPLES::

        sage: from sage.rings.invariants.reconstruction import binary_cubic_from_invariants
        sage: binary_cubic_from_invariants(1)
        (0, 1, 1, 0)
        sage: binary_cubic_from_invariants(0)
        Traceback (most recent call last):
        ...
        NotImplementedError: No distinction implemented for binary cubics with a double root.
    """
    if discriminant == 0:
        raise NotImplementedError('No distinction implemented for binary \
                                    cubics with a double root.')
    else:
        return (0, 1, 1, 0) # x * y * (x + y)

def binary_quintic_from_invariants(invariants, K=None, scaled=False, reduced=False):
    """
    Reconstruct a binary quintic from the values of its Clebsch invariants.

    INPUT:

    - ``invariants`` --  The values of the three or four invariants
      of the binary quintic.

    - ``K`` -- The field over which the quintic is defined.

    - ``scaled`` -- A boolean to determine wether the coefficients should
      be scaled so the result is independent of the scaling of the invariants.

    - ``reduced`` -- A boolean to determine wether the coefficients should
      reduced by dividing by their gcd.

    OUTPUT:

    A set of coefficients of a binary quintic, whose Clebsch invariants
    are equal to the given ``invariants`` up to a scaling.

    EXAMPLES::

        sage: from sage.rings.invariants.reconstruction import binary_quintic_from_invariants
        sage: R.<x0, x1> = QQ[]
        sage: p = 3*x1^5 + 6*x1^4*x0 + 3*x1^3*x0^2 + 4*x1^2*x0^3 - 5*x1*x0^4 + 4*x0^5
        sage: quintic = invariant_theory.binary_quintic(p, x0, x1)
        sage: invs = quintic.clebsch_invariants(as_tuple=True)
        sage: reconstructed = binary_quintic_from_invariants(invs)
        sage: reconstructed
        (9592267437341790539005557/244140625000000,
         2149296928207625556323004064707/610351562500000000,
         11149651890347700974453304786783/76293945312500000,
         122650775751894638395648891202734239/47683715820312500000,
         323996630945706528474286334593218447/11920928955078125000,
         1504506503644608395841632538558481466127/14901161193847656250000)

    The form obtained corresponds with the one found when setting the coordinates
    equal to the covariants alpha and beta::

        sage: alpha = quintic.alpha_covariant()
        sage: beta = quintic.beta_covariant()
        sage: g = matrix([[alpha(x0=1,x1=0),alpha(x0=0,x1=1)],[beta(x0=1,x1=0),beta(x0=0,x1=1)]])^-1
        sage: g = g*g.determinant()^-1
        sage: quintic.transformed(g).coeffs() == reconstructed
        True

    We can check that the invariants match by scaling the invariants
    B and C to see if they match::

        sage: newform = sum([ reconstructed[i]*x0^i*x1^(5-i) for i in range(6) ])
        sage: newquintic = invariant_theory.binary_quintic(newform, x0, x1)
        sage: scale = invs[0]/newquintic.A_invariant()
        sage: invs[1] == newquintic.B_invariant()*scale^2
        True
        sage: invs[2] == newquintic.C_invariant()*scale^3
        True

    If the invariant M vanishes, the coefficients are computed in a
    different way::

        sage: reconstructed = binary_quintic_from_invariants([3,1,2])
        sage: reconstructed
        (66741943359375/2097152,
         -125141143798828125/134217728,
         0,
         52793920040130615234375/34359738368,
         -19797720015048980712890625/1099511627776,
         -4454487003386020660400390625/17592186044416)
        sage: newform = sum([ reconstructed[i]*x0^i*x1^(5-i) for i in range(6) ])
        sage: newquintic = invariant_theory.binary_quintic(newform, x0, x1)
        sage: scale = 3/newquintic.A_invariant()
        sage: [3, newquintic.B_invariant()*scale^2, newquintic.C_invariant()*scale^3]
        [3, 1, 2]

    Several special cases::

        sage: quintic = invariant_theory.binary_quintic(x0^5 - x1^5, x0, x1)
        sage: invs = quintic.clebsch_invariants(as_tuple=True)
        sage: reconstructed = binary_quintic_from_invariants(invs)
        sage: reconstructed
        (1, 0, 0, 0, 0, 1)
        sage: quintic = invariant_theory.binary_quintic(x0*x1*(x0^3-x1^3), x0, x1)
        sage: invs = quintic.clebsch_invariants(as_tuple=True)
        sage: reconstructed = binary_quintic_from_invariants(invs)
        sage: reconstructed
        (0, 1, 0, 0, 1, 0)

    For fields of characteristic 2, 3 or 5, there is no reconstruction implemented::

        sage: binary_quintic_from_invariants([3,1,2], K=GF(5))
        Traceback (most recent call last):
        ...
        NotImplementedError: No reconstruction implemented for fields of characteristic 2, 3 or 5.
    """
    if reduced:
        if len(invariants) == 3:
            invariants = reduce_invariants(invariants, [1,2,3])
        else:
            invariants = reduce_invariants(invariants, [2,4,6,9])
    A, B, C = invariants[0:3]
    if K is None:
        K = FractionField(A.parent())
    if K.characteristic() in [2, 3, 5]:
        raise NotImplementedError('No reconstruction implemented for fields \
                                    of characteristic 2, 3 or 5.')
    M = 2*A*B - 3*C
    N = K(2)**-1 * (A*C-B**2)
    R2 = -K(2)**-1 * (A*N**2-2*B*M*N+C*M**2)
    scale = [1,1,1,1,1,1]
    from sage.functions.all import binomial, sqrt
    try:
        if R2.is_square():
            R = sqrt(R2)
        else:
            # if R2 is not a square, we scale the invariants in a suitable way
            # so that the 'new' R2 is a square
            # r = R2.squarefree_part() # slow!
            invariants = [R2*A, R2**2*B, R2**3*C]
            # we compute again with the new invariants, not reduced
            # (else the scaling is undone)
            coeffs = binary_quintic_from_invariants(invariants, K, scaled, \
                        reduced=False)
            if reduced:
                from sage.arith.misc import gcd
                return tuple([coeffs[i]/gcd(coeffs) for i in range(6)])
            else:
                return coeffs
    except (AttributeError, NotImplementedError):
        if len(invariants) > 3:
            R = invariants[3]
        else:
            raise ValueError('Value of R could not be determined.')
    if M == 0:
        if N == 0:
            if A == 0:
                raise NotImplementedError('No reconstruction implemented for \
                            quintics with a treefold linear factor.')
            else:
                if B == 0:
                    return (1,0,0,0,0,1) # x**5 + y**5
                else:
                    return (0,1,0,0,1,0) # x*y*(x**3+y**3)
        else:
            # case corresponding to using alpha and gamma as coordinates
            if A == 0:
                return (1,0,0,0,1,0) # x*(x**4+y**4)
            else:
                if scaled:
                    # subs = {y:(R/A**3)*y}
                    scale = [ A**-14*(R/A**3)**i for i in range(6) ]
                D = -N
                Delta = C
                a = [0]
                a.append((2*K(3)**-1*A**2-B)*N*B*K(2)**-1 - N**2*K(2)**-1)
                B0 = 2*K(3)**-1*A*R
                B1 = A*N*B*K(3)**-1
                C0 = 2*K(3)**-1*R
                C1 = B*N
    else:
        # case corresponding to using alpha and beta as coordinates
        if R == 0:
            if A == 0:
                return (1,0,10,0,-15,0) # x**5 + 10*x**3*y**2 - 15*x*y**4
            elif scaled:
                # subs = {x:x/A, y:y/sqrt(A)}
                scale = [ sqrt(A)**(i-18) for i in range(6) ]
        else:
            if A == 0:
                if B == 0:
                    return (1,0,0,1,0,0) # x**2*(x**3+y**3)
                elif scaled:
                    # subs = {y:(R/B**2)*y}
                    scale = [ R**-2*(R/B**2)**i for i in range(6) ]
            elif scaled:
                # subs = {y:(R/A**4)*y}
                scale = [ A**-9*(R/A**4)**i for i in range(6) ]
        D = -M
        Delta = A
        a = [0]
        a.append((2*K(3)**-1*A**2-B)*(N*A-M*B)*K(2)**-1 \
                    - M*(N*K(2)**-1-M*A*K(3)**-1))
        B0 = R
        B1 = K(2)**-1*(N*A-M*B)
        C0 = 0
        C1 = -M
    a[0] = (2*K(3)**-1*A**2-B)*R
    a.append(-D*B0 - K(2)**-1*Delta*a[0])
    a.append(-D*B1 - K(2)**-1*Delta*a[1])
    a.append(D**2*C0 + D*Delta*B0 + K(4)**-1*Delta**2*a[0])
    a.append(D**2*C1 + D*Delta*B1 + K(4)**-1*Delta**2*a[1])
    # D**(-5)*(a5*y**5 - 5*a4*x*y**4 + 10*a3*x**2*y**3 - 10*a2*x**3*y**2
    #           + 5*a1*x**4*y - a0*x**5)
    coeffs = tuple([K((-1)**i*binomial(5,i)*scale[5-i]*a[i]) for i in range(6)])
    if reduced:
        from sage.arith.misc import gcd
        return tuple([coeffs[i]/gcd(coeffs) for i in range(6)])
    else:
        return coeffs


######################################################################

def reduce_invariants(invariants, weights):
    """
    Reduce a list of invariants of given weights so that they have no common
    weighted divisor.

    INPUT:

    - ``invariants`` --  The values of the invariants.

    - ``weights`` --  The respective weights of the invariants.

    OUTPUT:

    A list of invariants that is equivalent to the input and
    has no common weighted divisor.

    EXAMPLES::

        sage: from sage.rings.invariants.reconstruction import reduce_invariants
        sage: invariants = [6/5, 12, 16]
        sage: weights = [1, 2, 3]
        sage: reduce_invariants(invariants, weights)
        [3, 75, 250]
    """
    factors = [dict(I.factor()) for I in invariants]
    scalar = ZZ(1)
    n = len(weights)
    from sage.arith.misc import gcd
    for prime in gcd(invariants).factor():
        p = prime[0]
        for D in factors:
            if not D.has_key(p):
                D[p] = 0
        scalar = scalar*p**min([factors[i][p]//weights[i] for i in range(n)])
    return [invariants[i]*scalar**-weights[i] for i in range(n)]

