r"""
Enumeration of rational points on affine schemes

Naive algorithms for enumerating rational points over `\QQ` or finite fields
over for general schemes.

.. WARNING::

    Incorrect results and infinite loops may occur if using a wrong function.

    (For instance using an affine function for a projective scheme or a finite
    field function for a scheme defined over an infinite field.)

EXAMPLES:

Affine, over `\QQ`::

    sage: from sage.schemes.affine.affine_rational_point import enum_affine_rational_field
    sage: A.<x,y,z> = AffineSpace(3, QQ)
    sage: S = A.subscheme([2*x-3*y])
    sage: enum_affine_rational_field(S, 2)
    [(0, 0, -2), (0, 0, -1), (0, 0, -1/2), (0, 0, 0),
     (0, 0, 1/2), (0, 0, 1), (0, 0, 2)]

Affine over a finite field::

    sage: from sage.schemes.affine.affine_rational_point import enum_affine_finite_field
    sage: A.<w,x,y,z> = AffineSpace(4, GF(2))
    sage: enum_affine_finite_field(A(GF(2)))
    [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0), (0, 0, 1, 1), (0, 1, 0, 0),
     (0, 1, 0, 1), (0, 1, 1, 0), (0, 1, 1, 1), (1, 0, 0, 0), (1, 0, 0, 1),
     (1, 0, 1, 0), (1, 0, 1, 1), (1, 1, 0, 0), (1, 1, 0, 1), (1, 1, 1, 0),
     (1, 1, 1, 1)]

AUTHORS:

- David R. Kohel <kohel@maths.usyd.edu.au>: original version.

- John Cremona and Charlie Turner <charlotteturner@gmail.com> (06-2010):
  improvements to clarity and documentation.

- Raghukul Raman <raghukul.raman01@gmail.com> (2018): Added sieve algorithm.

"""

#*****************************************************************************
#       Copyright (C) 2010 William Stein, David Kohel, John Cremona, Charlie Turner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from six.moves import range

from sage.rings.all import ZZ, QQ
from sage.misc.all import cartesian_product_iterator
from sage.schemes.generic.scheme import is_Scheme
from sage.schemes.projective.projective_rational_point import sieve as projective_sieve


def enum_affine_rational_field(X, B):
    """
    Enumerates affine rational points on scheme ``X`` up to bound ``B``.

    INPUT:

    - ``X`` -  a scheme or set of abstract rational points of a scheme.
    - ``B`` -  a positive integer bound.

    OUTPUT:

    - a list containing the affine points of ``X`` of height up to ``B``,
      sorted.

    EXAMPLES::

        sage: A.<x,y,z> = AffineSpace(3, QQ)
        sage: from sage.schemes.affine.affine_rational_point import enum_affine_rational_field
        sage: enum_affine_rational_field(A(QQ), 1)
        [(-1, -1, -1), (-1, -1, 0), (-1, -1, 1), (-1, 0, -1), (-1, 0, 0), (-1, 0, 1),
        (-1, 1, -1), (-1, 1, 0), (-1, 1, 1), (0, -1, -1), (0, -1, 0), (0, -1, 1),
        (0, 0, -1), (0, 0, 0), (0, 0, 1), (0, 1, -1), (0, 1, 0), (0, 1, 1), (1, -1, -1),
        (1, -1, 0), (1, -1, 1), (1, 0, -1), (1, 0, 0), (1, 0, 1), (1, 1, -1), (1, 1, 0),
        (1, 1, 1)]

    ::

        sage: A.<w,x,y,z> = AffineSpace(4, QQ)
        sage: S = A.subscheme([x^2-y*z+3, w^3+z+y^2])
        sage: enum_affine_rational_field(S(QQ), 2)
        []
        sage: enum_affine_rational_field(S(QQ), 3)
        [(-2, 0, -3, -1)]

    ::

        sage: A.<x,y> = AffineSpace(2, QQ)
        sage: C = Curve(x^2+y-x)
        sage: enum_affine_rational_field(C, 10)
        [(-2, -6), (-1, -2), (-2/3, -10/9), (-1/2, -3/4), (-1/3, -4/9),
        (0, 0), (1/3, 2/9), (1/2, 1/4), (2/3, 2/9), (1, 0),
        (4/3, -4/9), (3/2, -3/4), (5/3, -10/9), (2, -2), (3, -6)]

    AUTHORS:

    - David R. Kohel <kohel@maths.usyd.edu.au>: original version.

    - Charlie Turner (06-2010): small adjustments.

    - Raman Raghukul 2018: updated.
    """
    from sage.schemes.affine.affine_space import is_AffineSpace
    if(is_Scheme(X)):
        if (not is_AffineSpace(X.ambient_space())):
            raise TypeError("ambient space must be affine space over the rational field")
        X = X(X.base_ring())
    else:
        if (not is_AffineSpace(X.codomain().ambient_space())):
            raise TypeError("codomain must be affine space over the rational field")

    n = X.codomain().ambient_space().ngens()
    VR = X.value_ring()
    if VR is ZZ:
        R = [ 0 ] + [ s*k for k in range(1, B+1) for s in [1, -1] ]
        iters = [ iter(R) for _ in range(n) ]
    else:  # rational field
        iters = [ QQ.range_by_height(B + 1) for _ in range(n) ]
    pts = []
    P = [0] * n
    try:
        pts.append(X(P))
    except TypeError:
        pass

    for it in iters:
        next(it)
    i = 0
    while i < n:
        try:
            a = VR(next(iters[i]))
        except StopIteration:
            if VR is ZZ:
                iters[i] = iter(R)
            else:  # rational field
                iters[i] = QQ.range_by_height(B + 1)
            P[i] = next(iters[i]) # reset P[i] to 0 and increment
            i += 1
            continue
        P[i] = a
        try:
            pts.append(X(P))
        except TypeError:
            pass
        i = 0
    pts.sort()
    return pts


def enum_affine_number_field(X, **kwds):
    """
    Enumerates affine points on scheme ``X`` defined over a number field. Simply checks all of the
    points of absolute height up to ``B`` and adds those that are on the scheme to the list.

    This algorithm computes 2 lists: L containing elements x in `K` such that
    H_k(x) <= B, and a list L' containing elements x in `K` that, due to
    floating point issues,
    may be slightly larger then the bound. This can be controlled
    by lowering the tolerance.

    ALGORITHM:

    This is an implementation of the revised algorithm (Algorithm 4) in
    [Doyle-Krumm]_. Algorithm 5 is used for imaginary quadratic fields.


    INPUT:

    kwds:

    - ``bound`` - a real number

    - ``tolerance`` - a rational number in (0,1] used in doyle-krumm algorithm-4

    - ``precision`` - the precision to use for computing the elements of bounded height of number fields.

    OUTPUT:

     - a list containing the affine points of ``X`` of absolute height up to ``B``,
       sorted.

    EXAMPLES::

        sage: from sage.schemes.affine.affine_rational_point import enum_affine_number_field
        sage: u = QQ['u'].0
        sage: K = NumberField(u^2 + 2, 'v')
        sage: A.<x,y,z> = AffineSpace(K, 3)
        sage: X = A.subscheme([y^2 - x])
        sage: enum_affine_number_field(X(K), bound=2**0.5)
        [(0, 0, -1), (0, 0, -v), (0, 0, -1/2*v), (0, 0, 0), (0, 0, 1/2*v), (0, 0, v), (0, 0, 1),
        (1, -1, -1), (1, -1, -v), (1, -1, -1/2*v), (1, -1, 0), (1, -1, 1/2*v), (1, -1, v), (1, -1, 1),
        (1, 1, -1), (1, 1, -v), (1, 1, -1/2*v), (1, 1, 0), (1, 1, 1/2*v), (1, 1, v), (1, 1, 1)]

    ::

        sage: u = QQ['u'].0
        sage: K = NumberField(u^2 + 3, 'v')
        sage: A.<x,y> = AffineSpace(K, 2)
        sage: X=A.subscheme(x-y)
        sage: from sage.schemes.affine.affine_rational_point import enum_affine_number_field
        sage: enum_affine_number_field(X, bound=3**0.25)
        [(-1, -1), (-1/2*v - 1/2, -1/2*v - 1/2), (1/2*v - 1/2, 1/2*v - 1/2), (0, 0), (-1/2*v + 1/2, -1/2*v + 1/2),
        (1/2*v + 1/2, 1/2*v + 1/2), (1, 1)]
    """
    B = kwds.pop('bound')
    tol = kwds.pop('tolerance', 1e-2)
    prec = kwds.pop('precision', 53)
    from sage.schemes.affine.affine_space import is_AffineSpace
    if(is_Scheme(X)):
        if (not is_AffineSpace(X.ambient_space())):
            raise TypeError("ambient space must be affine space over a number field")
        X = X(X.base_ring())
    else:
        if (not is_AffineSpace(X.codomain().ambient_space())):
            raise TypeError("codomain must be affine space over a number field")

    R = X.codomain().ambient_space()

    pts = []
    for P in R.points_of_bounded_height(bound=B, tolerance=tol, precision=prec):
        try:
            pts.append(X(P))
        except TypeError:
            pass
    pts.sort()
    return pts


def enum_affine_finite_field(X):
    r"""
    Enumerates affine points on scheme ``X`` defined over a finite field.

    INPUT:

    - ``X`` -  a scheme defined over a finite field or a set of abstract
      rational points of such a scheme.

    OUTPUT:

    - a list containing the affine points of ``X`` over the finite field,
      sorted.

    EXAMPLES::

        sage: F = GF(7)
        sage: A.<w,x,y,z> = AffineSpace(4, F)
        sage: C = A.subscheme([w^2+x+4, y*z*x-6, z*y+w*x])
        sage: from sage.schemes.affine.affine_rational_point import enum_affine_finite_field
        sage: enum_affine_finite_field(C(F))
        []
        sage: C = A.subscheme([w^2+x+4, y*z*x-6])
        sage: enum_affine_finite_field(C(F))
        [(0, 3, 1, 2), (0, 3, 2, 1), (0, 3, 3, 3), (0, 3, 4, 4), (0, 3, 5, 6),
        (0, 3, 6, 5), (1, 2, 1, 3), (1, 2, 2, 5), (1, 2, 3, 1), (1, 2, 4, 6),
        (1, 2, 5, 2), (1, 2, 6, 4), (2, 6, 1, 1), (2, 6, 2, 4), (2, 6, 3, 5),
        (2, 6, 4, 2), (2, 6, 5, 3), (2, 6, 6, 6), (3, 1, 1, 6), (3, 1, 2, 3),
        (3, 1, 3, 2), (3, 1, 4, 5), (3, 1, 5, 4), (3, 1, 6, 1), (4, 1, 1, 6),
        (4, 1, 2, 3), (4, 1, 3, 2), (4, 1, 4, 5), (4, 1, 5, 4), (4, 1, 6, 1),
        (5, 6, 1, 1), (5, 6, 2, 4), (5, 6, 3, 5), (5, 6, 4, 2), (5, 6, 5, 3),
        (5, 6, 6, 6), (6, 2, 1, 3), (6, 2, 2, 5), (6, 2, 3, 1), (6, 2, 4, 6),
        (6, 2, 5, 2), (6, 2, 6, 4)]

    ::

        sage: A.<x,y,z> = AffineSpace(3, GF(3))
        sage: S = A.subscheme(x+y)
        sage: enum_affine_finite_field(S)
        [(0, 0, 0), (0, 0, 1), (0, 0, 2), (1, 2, 0), (1, 2, 1), (1, 2, 2),
        (2, 1, 0), (2, 1, 1), (2, 1, 2)]

    ALGORITHM:

    Checks all points in affine space to see if they lie on X.

    .. WARNING::

        If ``X`` is defined over an infinite field, this code will not finish!

    AUTHORS:

    - John Cremona and Charlie Turner (06-2010)
    """
    from sage.schemes.affine.affine_space import is_AffineSpace
    if(is_Scheme(X)):
        if (not is_AffineSpace(X.ambient_space())):
            raise TypeError("ambient space must be affine space over a finite field")
        X = X(X.base_ring())
    else:
        if (not is_AffineSpace(X.codomain().ambient_space())):
            raise TypeError("codomain must be affine space over a finite field")

    n = X.codomain().ambient_space().ngens()
    F = X.value_ring()
    pts = []
    for c in cartesian_product_iterator([F]*n):
        try:
            pts.append(X(c))
        except Exception:
            pass
    pts.sort()
    return pts


def sieve(X, bound):
    r"""
    Returns the list of all affine, rational points on scheme ``X`` of
    height up to ``bound``.

    This algorithm works correctly only if dimension of given
    scheme is positive.

    INPUT:

    - ``X`` - a scheme with ambient space defined over affine space

    - ``bound`` - a positive integer bound

    OUTPUT:

     - a list containing the affine rational points of ``X`` of height
    up to ``B``, sorted

    EXAMPLES::
        
        sage: from sage.schemes.affine.affine_rational_point import sieve
        sage: A.<x,y,z> = AffineSpace(3, QQ)
        sage: S = A.subscheme([x - 2*y - 3*z,x^2 - z^2])
        sage: sieve(S, 2)
        [(-2, 2, -2), (-1, -2, 1), (-1, 1, -1), (-1/2, -1, 1/2),
         (-1/2, 1/2, -1/2), (0, 0, 0), (1/2, -1/2, 1/2),
         (1/2, 1, -1/2), (1, -1, 1), (1, 2, -1), (2, -2, 2)]
    
    TESTS:

    Illustrates efficiency of algorithm::

        sage: A.<x,y> = AffineSpace(2, QQ)
        sage: C = Curve(x^2+y^2-x)
        sage: len(C.rational_points(20))
        20
        sage: from sage.schemes.affine.affine_rational_point import sieve
        sage: len(sieve(C, 20))
        20

    """
    # finds a projective embedding to use projective version of sieve
    pi = X.projective_embedding(0)
    P = pi.codomain()
    AA = P.affine_patch(0)
    
    proj_L = projective_sieve(P, bound)
    LL = set()
    for point in proj_L: # make them back into affine points
        pt = []
        denom = point[0]
        if denom == 0:
            continue

        for i in range(1,len(point)):
            pt.append(point[i] / denom)
        LL.add(AA(pt))

    return sorted(list(LL))