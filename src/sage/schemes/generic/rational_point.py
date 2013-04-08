r"""
Enumeration of rational points on affine and projective schemes

Naive algorithms for enumerating rational points over `\QQ` or finite fields
over for general schemes.

.. WARNING::

    Incorrect results and infinite loops may occur if using a wrong function.
    (For instance using an affine function for a projective scheme or a finite
    field function for a scheme defined over an infinite field.)

EXAMPLES:

Projective, over `\QQ`::

    sage: from sage.schemes.generic.rational_point import enum_projective_rational_field
    sage: P.<X,Y,Z> = ProjectiveSpace(2,QQ)
    sage: C = P.subscheme([X+Y-Z])
    sage: enum_projective_rational_field(C,3)
    [(-2 : 3 : 1), (-1 : 1 : 0), (-1 : 2 : 1), (-1/2 : 3/2 : 1),
     (0 : 1 : 1), (1/3 : 2/3 : 1), (1/2 : 1/2 : 1), (2/3 : 1/3 : 1),
     (1 : 0 : 1), (3/2 : -1/2 : 1), (2 : -1 : 1), (3 : -2 : 1)]

Affine, over `\QQ`::

    sage: from sage.schemes.generic.rational_point import enum_affine_rational_field
    sage: A.<x,y,z> = AffineSpace(3,QQ)
    sage: S = A.subscheme([2*x-3*y])
    sage: enum_affine_rational_field(S,2)
    [(0, 0, -2), (0, 0, -1), (0, 0, -1/2), (0, 0, 0),
     (0, 0, 1/2), (0, 0, 1), (0, 0, 2)]

Projective over a finite field::

    sage: from sage.schemes.generic.rational_point import enum_projective_finite_field
    sage: E = EllipticCurve('72').change_ring(GF(19))
    sage: enum_projective_finite_field(E)
    [(0 : 1 : 0), (1 : 0 : 1), (3 : 0 : 1), (4 : 9 : 1), (4 : 10 : 1),
     (6 : 6 : 1), (6 : 13 : 1), (7 : 6 : 1), (7 : 13 : 1), (9 : 4 : 1),
     (9 : 15 : 1), (12 : 8 : 1), (12 : 11 : 1), (13 : 8 : 1), (13 : 11 : 1),
     (14 : 3 : 1), (14 : 16 : 1), (15 : 0 : 1), (16 : 9 : 1), (16 : 10 : 1),
     (17 : 7 : 1), (17 : 12 : 1), (18 : 9 : 1), (18 : 10 : 1)]

Affine over a finite field::

    sage: from sage.schemes.generic.rational_point import enum_affine_finite_field
    sage: A.<w,x,y,z> = AffineSpace(4,GF(2))
    sage: enum_affine_finite_field(A(GF(2)))
    [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0), (0, 0, 1, 1), (0, 1, 0, 0),
     (0, 1, 0, 1), (0, 1, 1, 0), (0, 1, 1, 1), (1, 0, 0, 0), (1, 0, 0, 1),
     (1, 0, 1, 0), (1, 0, 1, 1), (1, 1, 0, 0), (1, 1, 0, 1), (1, 1, 1, 0),
     (1, 1, 1, 1)]

AUTHORS:

- David R. Kohel <kohel@maths.usyd.edu.au>: original version.

- John Cremona and Charlie Turner <charlotteturner@gmail.com> (06-2010):
  improvements to clarity and documentation.
"""


#*******************************************************************************
#  Copyright (C) 2010 William Stein, David Kohel, John Cremona, Charlie Turner
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************


from sage.rings.arith import gcd
from sage.rings.all import ZZ
from sage.misc.all import srange, cartesian_product_iterator
from sage.schemes.all import is_Scheme

def enum_projective_rational_field(X,B):
    r"""
    Enumerates projective, rational points on scheme ``X`` of height up to
    bound ``B``.

    INPUT:

    - ``X`` -  a scheme or set of abstract rational points of a scheme;
    - ``B`` -  a positive integer bound.

    OUTPUT:

    - a list containing the projective points of ``X`` of height up to ``B``,
      sorted.

    EXAMPLES::

        sage: P.<X,Y,Z> = ProjectiveSpace(2,QQ)
        sage: C = P.subscheme([X+Y-Z])
        sage: from sage.schemes.generic.rational_point import enum_projective_rational_field
        sage: enum_projective_rational_field(C(QQ),6)
        [(-5 : 6 : 1), (-4 : 5 : 1), (-3 : 4 : 1), (-2 : 3 : 1),
         (-3/2 : 5/2 : 1), (-1 : 1 : 0), (-1 : 2 : 1), (-2/3 : 5/3 : 1),
         (-1/2 : 3/2 : 1), (-1/3 : 4/3 : 1), (-1/4 : 5/4 : 1),
         (-1/5 : 6/5 : 1), (0 : 1 : 1), (1/6 : 5/6 : 1), (1/5 : 4/5 : 1),
         (1/4 : 3/4 : 1), (1/3 : 2/3 : 1), (2/5 : 3/5 : 1), (1/2 : 1/2 : 1),
         (3/5 : 2/5 : 1), (2/3 : 1/3 : 1), (3/4 : 1/4 : 1), (4/5 : 1/5 : 1),
         (5/6 : 1/6 : 1), (1 : 0 : 1), (6/5 : -1/5 : 1), (5/4 : -1/4 : 1),
         (4/3 : -1/3 : 1), (3/2 : -1/2 : 1), (5/3 : -2/3 : 1), (2 : -1 : 1),
         (5/2 : -3/2 : 1), (3 : -2 : 1), (4 : -3 : 1), (5 : -4 : 1),
         (6 : -5 : 1)]
        sage: enum_projective_rational_field(C,6) == enum_projective_rational_field(C(QQ),6)
        True

    ::

        sage: P3.<W,X,Y,Z> = ProjectiveSpace(3,QQ)
        sage: enum_projective_rational_field(P3,1)
        [(-1 : -1 : -1 : 1), (-1 : -1 : 0 : 1), (-1 : -1 : 1 : 0), (-1 : -1 : 1 : 1),
        (-1 : 0 : -1 : 1), (-1 : 0 : 0 : 1), (-1 : 0 : 1 : 0), (-1 : 0 : 1 : 1),
        (-1 : 1 : -1 : 1), (-1 : 1 : 0 : 0), (-1 : 1 : 0 : 1), (-1 : 1 : 1 : 0),
        (-1 : 1 : 1 : 1), (0 : -1 : -1 : 1), (0 : -1 : 0 : 1), (0 : -1 : 1 : 0),
        (0 : -1 : 1 : 1), (0 : 0 : -1 : 1), (0 : 0 : 0 : 1), (0 : 0 : 1 : 0),
        (0 : 0 : 1 : 1), (0 : 1 : -1 : 1), (0 : 1 : 0 : 0), (0 : 1 : 0 : 1),
        (0 : 1 : 1 : 0), (0 : 1 : 1 : 1), (1 : -1 : -1 : 1), (1 : -1 : 0 : 1),
        (1 : -1 : 1 : 0), (1 : -1 : 1 : 1), (1 : 0 : -1 : 1), (1 : 0 : 0 : 0),
        (1 : 0 : 0 : 1), (1 : 0 : 1 : 0), (1 : 0 : 1 : 1), (1 : 1 : -1 : 1),
        (1 : 1 : 0 : 0), (1 : 1 : 0 : 1), (1 : 1 : 1 : 0), (1 : 1 : 1 : 1)]

    ALGORITHM:

    We just check all possible projective points in correct dimension
    of projective space to see if they lie on ``X``.

    AUTHORS:

    - John Cremona and Charlie Turner (06-2010)
    """
    if is_Scheme(X):
        X = X(X.base_ring())
    n = X.codomain().ambient_space().ngens()
    zero = (0,) * n
    pts = []
    for c in cartesian_product_iterator([srange(-B,B+1) for _ in range(n)]):
        if gcd(c) == 1 and c > zero:
            try:
                pts.append(X(c))
            except TypeError:
                pass
    pts.sort()
    return pts

def enum_affine_rational_field(X,B):
    """
    Enumerates affine rational points on scheme ``X`` (defined over `\QQ`) up
    to bound ``B``.

    INPUT:

    - ``X`` -  a scheme or set of abstract rational points of a scheme;
    - ``B`` -  a positive integer bound.

    OUTPUT:

    - a list containing the affine points of ``X`` of height up to ``B``,
      sorted.

    EXAMPLES::

        sage: A.<x,y,z> = AffineSpace(3,QQ)
        sage: from sage.schemes.generic.rational_point import enum_affine_rational_field
        sage: enum_affine_rational_field(A(QQ),1)
        [(-1, -1, -1), (-1, -1, 0), (-1, -1, 1), (-1, 0, -1), (-1, 0, 0), (-1, 0, 1),
        (-1, 1, -1), (-1, 1, 0), (-1, 1, 1), (0, -1, -1), (0, -1, 0), (0, -1, 1),
        (0, 0, -1), (0, 0, 0), (0, 0, 1), (0, 1, -1), (0, 1, 0), (0, 1, 1), (1, -1, -1),
        (1, -1, 0), (1, -1, 1), (1, 0, -1), (1, 0, 0), (1, 0, 1), (1, 1, -1), (1, 1, 0),
        (1, 1, 1)]

    ::

        sage: A.<w,x,y,z> = AffineSpace(4,QQ)
        sage: S = A.subscheme([x^2-y*z+3,w^3+z+y^2])
        sage: enum_affine_rational_field(S(QQ),2)
        []
        sage: enum_affine_rational_field(S(QQ),3)
        [(-2, 0, -3, -1)]

    ::

        sage: A.<x,y> = AffineSpace(2,QQ)
        sage: C = Curve(x^2+y-x)
        sage: enum_affine_rational_field(C,10)
        [(-2, -6), (-1, -2), (0, 0), (1, 0), (2, -2), (3, -6)]


    AUTHORS:

    - David R. Kohel <kohel@maths.usyd.edu.au>: original version.

    - Charlie Turner (06-2010): small adjustments.
    """
    if is_Scheme(X):
        X = X(X.base_ring())
    n = X.codomain().ambient_space().ngens()
    if X.value_ring() is ZZ:
        Q = [ 1 ]
    else: # rational field
        Q = range(1, B + 1)
    R = [ 0 ] + [ s*k for k in range(1, B+1) for s in [1, -1] ]
    pts = []
    P = [0] * n
    m = ZZ(0)
    try:
        pts.append(X(P))
    except TypeError:
        pass
    iters = [ iter(R) for _ in range(n) ]
    for it in iters:
        it.next()
    i = 0
    while i < n:
        try:
            a = ZZ(iters[i].next())
        except StopIteration:
            iters[i] = iter(R) # reset
            P[i] = iters[i].next() # reset P[i] to 0 and increment
            i += 1
            continue
        m = m.gcd(a)
        P[i] = a
        for b in Q:
            if m.gcd(b) == 1:
                try:
                    pts.append(X([ num/b for num in P ]))
                except TypeError:
                    pass
        i = 0
        m = ZZ(0)
    pts.sort()
    return pts

def enum_projective_finite_field(X):
    """
    Enumerates projective points on scheme ``X`` defined over a finite field.

    INPUT:

    - ``X`` -  a scheme defined over a finite field or a set of abstract
      rational points of such a scheme.

    OUTPUT:

    - a list containing the projective points of ``X`` over the finite field,
      sorted.

    EXAMPLES::

        sage: F = GF(53)
        sage: P.<X,Y,Z> = ProjectiveSpace(2,F)
        sage: from sage.schemes.generic.rational_point import enum_projective_finite_field
        sage: len(enum_projective_finite_field(P(F)))
        2863
        sage: 53^2+53+1
        2863

    ::

        sage: F = GF(9,'a')
        sage: P.<X,Y,Z> = ProjectiveSpace(2,F)
        sage: C = Curve(X^3-Y^3+Z^2*Y)
        sage: enum_projective_finite_field(C(F))
        [(0 : 0 : 1), (0 : 1 : 1), (0 : 2 : 1), (1 : 1 : 0), (a + 1 : 2*a : 1),
        (a + 1 : 2*a + 1 : 1), (a + 1 : 2*a + 2 : 1), (2*a + 2 : a : 1),
        (2*a + 2 : a + 1 : 1), (2*a + 2 : a + 2 : 1)]

    ::

        sage: F = GF(5)
        sage: P2F.<X,Y,Z> = ProjectiveSpace(2,F)
        sage: enum_projective_finite_field(P2F)
        [(0 : 0 : 1), (0 : 1 : 0), (0 : 1 : 1), (0 : 2 : 1), (0 : 3 : 1), (0 : 4 : 1),
        (1 : 0 : 0), (1 : 0 : 1), (1 : 1 : 0), (1 : 1 : 1), (1 : 2 : 1), (1 : 3 : 1),
        (1 : 4 : 1), (2 : 0 : 1), (2 : 1 : 0), (2 : 1 : 1), (2 : 2 : 1), (2 : 3 : 1),
        (2 : 4 : 1), (3 : 0 : 1), (3 : 1 : 0), (3 : 1 : 1), (3 : 2 : 1), (3 : 3 : 1),
        (3 : 4 : 1), (4 : 0 : 1), (4 : 1 : 0), (4 : 1 : 1), (4 : 2 : 1), (4 : 3 : 1),
        (4 : 4 : 1)]

    ALGORITHM:

    Checks all points in projective space to see if they lie on X.

    .. WARNING::

        If ``X`` is defined over an infinite field, this code will not finish!

    AUTHORS:

    - John Cremona and Charlie Turner (06-2010).
    """
    if is_Scheme(X):
        X = X(X.base_ring())
    n = X.codomain().ambient_space().ngens()-1
    F = X.value_ring()
    pts = []
    for k in range(n+1):
        for c in cartesian_product_iterator([F for _ in range(k)]):
            try:
                pts.append(X(list(c)+[1]+[0]*(n-k)))
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
        sage: A.<w,x,y,z> = AffineSpace(4,F)
        sage: C = A.subscheme([w^2+x+4,y*z*x-6,z*y+w*x])
        sage: from sage.schemes.generic.rational_point import enum_affine_finite_field
        sage: enum_affine_finite_field(C(F))
        []
        sage: C = A.subscheme([w^2+x+4,y*z*x-6])
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

        sage: A.<x,y,z> = AffineSpace(3,GF(3))
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
    if is_Scheme(X):
        X = X(X.base_ring())
    n = X.codomain().ambient_space().ngens()
    F = X.value_ring()
    pts = []
    for c in cartesian_product_iterator([F]*n):
        try:
            pts.append(X(c))
        except StandardError:
            pass
    pts.sort()
    return pts
