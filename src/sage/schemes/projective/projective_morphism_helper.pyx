r"""
Morphisms on projective varieties (Cython helper)

This is the helper file providing functionality for projective_morphism.py.

AUTHORS:

- Dillon Rose (2014-01):  Speed enhancements

- Ben Hutz (2015-11): subscheme iteration

"""

#*****************************************************************************
#       Copyright (C) 2014 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.arith              import lcm
from sage.rings.finite_rings.constructor import GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.all                 import Set
from sage.misc.misc                import subsets

def _fast_possible_periods(self,return_points=False):
    r"""
    Returns the list of possible minimal periods of a periodic point
    over `\QQ` and (optionally) a point in each cycle.

    ALGORITHM:

    The list comes from B. Hutz. Good reduction of periodic points, Illinois Journal of
    Mathematics 53 (Winter 2009), no. 4, 1109-1126..

    INPUT:

    - ``return_points`` - Boolean (optional) - a value of True returns the points as well as the possible periods.

    OUTPUT:

    - a list of positive integers, or a list of pairs of projective points and periods if ``flag`` is 1.

    Examples::

            sage: from sage.schemes.projective.projective_morphism_helper import _fast_possible_periods
            sage: P.<x,y>=ProjectiveSpace(GF(23),1)
            sage: H=Hom(P,P)
            sage: f=H([x^2-2*y^2,y^2])
            sage: _fast_possible_periods(f,False)
            [1, 5, 11, 22, 110]

        ::

            sage: from sage.schemes.projective.projective_morphism_helper import _fast_possible_periods
            sage: P.<x,y> = ProjectiveSpace(GF(13),1)
            sage: H = End(P)
            sage: f = H([x^2-y^2,y^2])
            sage: sorted(_fast_possible_periods(f,True))
            [[(0 : 1), 2], [(1 : 0), 1], [(3 : 1), 3], [(3 : 1), 36]]

        ::

            sage: from sage.schemes.projective.projective_morphism_helper import _fast_possible_periods
            sage: PS.<x,y,z> = ProjectiveSpace(2,GF(7))
            sage: H = End(PS)
            sage: f = H([-360*x^3 + 760*x*z^2, y^3 - 604*y*z^2 + 240*z^3, 240*z^3])
            sage: _fast_possible_periods(f,False)
            [1, 2, 4, 6, 12, 14, 28, 42, 84]

        .. TODO::

        - more space efficient hash/pointtable
    """
    cdef int i, k
    cdef list pointslist

    if not self._is_prime_finite_field:
        raise TypeError("Must be prime field")
    from sage.schemes.projective.projective_space import is_ProjectiveSpace
    if is_ProjectiveSpace(self.domain()) == False or self.domain()!=self.codomain():
        raise NotImplementedError("Must be an endomorphism of projective space")

    PS = self.domain()
    p = PS.base_ring().order()
    N = PS.dimension_relative()

    point_table = [[0,0] for i in xrange(p**(N + 1))]
    index = 1
    periods = set()
    points_periods = []

    for P in _enum_points(p, N):

        hash_p = _hash(P, p)
        if point_table[hash_p][1] == 0:
            startindex = index
            while point_table[hash_p][1] == 0:
                point_table[hash_p][1] = index
                Q = self._fast_eval(P)
                Q = _normalize_coordinates(Q, p, N+1)
                hash_q = _hash(Q, p)
                point_table[hash_p][0] = hash_q
                P=Q
                hash_p=hash_q
                index+=1

            if point_table[hash_p][1] >= startindex:
                P_proj=PS(P)
                period=index-point_table[hash_p][1]
                periods.add(period)
                points_periods.append([P_proj,period])
                l=P_proj.multiplier(self,period,False)
                lorders=set()
                for poly,_ in l.charpoly().factor():
                    if poly.degree() == 1:
                        eig = -poly.constant_coefficient()
                        if not eig:
                            continue # exclude 0
                    else:
                        eig = GF(p ** poly.degree(), 't', modulus=poly).gen()
                    if eig:
                        lorders.add(eig.multiplicative_order())
                S = subsets(lorders)
                next(S)   # get rid of the empty set
                rvalues=set()
                for s in S:
                    rvalues.add(lcm(s))
                rvalues=list(rvalues)
                if N==1:
                    for k in xrange(len(rvalues)):
                        r=rvalues[k]
                        periods.add(period*r)
                        points_periods.append([P_proj,period*r])
                        if p == 2 or p == 3: #need e=1 for N=1, QQ
                            periods.add(period*r*p)
                            points_periods.append([P_proj,period*r*p])
                else:
                    for k in xrange(len(rvalues)):
                        r=rvalues[k]
                        periods.add(period*r)
                        periods.add(period*r*p)
                        points_periods.append([P_proj,period*r])
                        points_periods.append([P_proj,period*r*p])
                        if p==2:  #need e=3 for N>1, QQ
                            periods.add(period*r*4)
                            points_periods.append([P_proj,period*r*4])
                            periods.add(period*r*8)
                            points_periods.append([P_proj,period*r*8])

    if return_points==False:
        return sorted(periods)
    else:
        return(points_periods)

def _enum_points(int prime,int dimension):
    """
    Enumerate points in projective space over finite field with given prime and dimension.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_morphism_helper import _enum_points
        sage: list(_enum_points(3,2))
        [[1, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0], [0, 0, 1], [1, 0, 1], [2, 0, 1], [0, 1, 1], [1, 1, 1], [2, 1, 1], [0, 2, 1], [1, 2, 1], [2, 2, 1]]
    """
    cdef int current_range
    cdef int highest_range
    cdef int value

    current_range = 1
    highest_range = prime**dimension

    while current_range <= highest_range:
        for value in xrange(current_range, 2*current_range):
            yield _get_point_from_hash(value,prime,dimension)
        current_range = current_range*prime

def _hash(list Point,int prime):
    """
    Hash point given as list to unique number.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_morphism_helper import _hash
        sage: _hash([1, 2, 1], 3)
        16

    """
    cdef int hash_q
    cdef int coefficient

    Point.reverse()
    hash_q = 0

    for coefficient in Point:
        hash_q = hash_q * prime + coefficient

    Point.reverse()

    return hash_q

def _get_point_from_hash(int value,int prime,int dimension):
    """
    Hash unique number to point as a list.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_morphism_helper import _get_point_from_hash
        sage: _get_point_from_hash(16,3,2)
        [1, 2, 1]

    """
    cdef list P
    cdef int i
    P=[]

    for i in xrange(dimension + 1):
        P.append(value % prime)
        value = value / prime

    return P

def _mod_inv(int num, int prime):
    """
    Find the inverse of the number modulo the given prime.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_morphism_helper import _mod_inv
        sage: _mod_inv(2,7)
        4
    """
    cdef int a, b, q, t, x, y
    a = prime
    b = num
    x = 1
    y = 0
    while b != 0:
        t = b
        q = a/t
        b = a - q*t
        a = t
        t = x
        x = y - q*t
        y = t

    if y < 0:
        return y + prime
    else:
        return y

def _normalize_coordinates(list point, int prime, int len_points):
    """
    Normalize the coordinates of the point for the given prime.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_morphism_helper import _normalize_coordinates
        sage: _normalize_coordinates([1,5,1],3,3)
        [1, 2, 1]

    """
    cdef int last_coefficient, coefficient, mod_inverse

    for coefficient in xrange(len_points):
        point[coefficient] = (point[coefficient]+prime)%prime
        if point[coefficient] != 0:
            last_coefficient = point[coefficient]

    mod_inverse = _mod_inv(last_coefficient,prime)

    for coefficient in xrange(len_points):
        point[coefficient] = (point[coefficient]*mod_inverse)%prime

    return point

def _forward_image_subscheme(self, X):
    """
    Compute the forward image of a projective subscheme ``x`` by ``self``.

    The forward image is compute through elimination.
    In particular, let $X = V(h_1,\ldots, h_t)$ and define the ideal
    $I = (h_1,\ldots,h_t,y_0-f_0(\bar{x}), \ldots, y_n-f_n(\bar{x}))$.
    Then the elimination ideal $I_{n+1} = I \cap K[y_0,\ldots,y_n]$ is a homogeneous
    ideal and $self(X) = V(I_{n+1})$.

    INPUT:

    - ``X`` - a subscheme in the domain of ``self``

    OUTPUT:

    - a subscheme in the domain of ``self``

    EXAMPLES:

        sage: PS.<x,y,z,w> = ProjectiveSpace(CC, 3)
        sage: H = End(PS)
        sage: f = H([x^2 + y^2, y^2, z^2-y^2, w^2])
        sage: X = PS.subscheme([z-2*w])
        sage: f(X)
        Closed subscheme of Projective Space of dimension 3 over Complex Field
        with 53 bits of precision defined by:
          y + z + (-4.00000000000000)*w

    ::

        sage: R.<t> = PolynomialRing(QQ)
        sage: P.<x,y,z> = ProjectiveSpace(FractionField(R), 2)
        sage: H = End(P)
        sage: f = H([x^2 + 2*y*z, t^2*y^2, z^2])
        sage: f([t^2*y-z])
        Closed subscheme of Projective Space of dimension 2 over Fraction Field
        of Univariate Polynomial Ring in t over Rational Field defined by:
          y + (-1/t^2)*z

    ::

        sage: set_verbose(-1)
        sage: PS.<x,y,z> = ProjectiveSpace(Qp(3), 2)
        sage: H = End(PS)
        sage: f = H([x^2,2*y^2,z^2])
        sage: X = PS.subscheme([2*x-y,z])
        sage: f(X)
        Closed subscheme of Projective Space of dimension 2 over 3-adic Field
        with capped relative precision 20 defined by:
          z,
          x + (1 + 3^2 + 3^4 + 3^6 + 3^8 + 3^10 + 3^12 + 3^14 + 3^16 + 3^18 +
        O(3^20))*y
    """
    dom = self.domain()
    codom = self.codomain()
    CR_dom = dom.coordinate_ring()
    CR_codom = codom.coordinate_ring()
    n = CR_dom.ngens()
    m = CR_codom.ngens()
    #can't call eliminate if the base ring is polynomial so we do it ourselves
    #with a lex ordering
    R = PolynomialRing(self.base_ring(), n+m, 'y', order = 'lex')
    Rvars = R.gens()[0 : n]
    phi = CR_dom.hom(Rvars,R)
    zero = [0 for _ in range(n)]
    psi = R.hom(zero + list(CR_codom.gens()),CR_codom)
    #set up ideal
    L = R.ideal([phi(t) for t in X.defining_polynomials()] + [R.gen(n+i) - phi(self[i]) for i in range(m)])
    G = L.groebner_basis() #eliminate
    newL = []
    #get only the elimination ideal portion
    for i in range (len(G)-1,0,-1):
        v = G[i].variables()
        if all([Rvars[j] not in v for j in range(n)]):
            newL.append(psi(G[i]))
    return(codom.subscheme(newL))
