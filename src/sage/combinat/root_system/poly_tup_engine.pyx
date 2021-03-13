"""
Arithmetic Engine for polynomials as tuples
"""
# ****************************************************************************
#  Copyright (C) 2021 Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from functools import cmp_to_key
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular, MPolynomialRing_libsingular
from sage.rings.polynomial.polydict cimport ETuple
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.rational_field import QQ

from sage.arith.functions cimport LCM_list

#Pre-compute common values for speed
one = QQ.one()
degrevlex_sortkey = TermOrder().sortkey_degrevlex

###########
### API ###
###########

cpdef tuple poly_to_tup(MPolynomial_libsingular poly):
    """
    Convert a polynomial object into the internal representation as tuple of
    (ETuple exp, NumberFieldElement coeff) pairs

    EXAMPLES::

        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
        sage: R.<x,y> = PolynomialRing(QQ)
        sage: poly_to_tup(x**2 + 1)
        (((2, 0), 1), ((0, 0), 1))
        sage: poly_to_tup(x**2*y**4 - 4/5*x*y**2 + 1/3 * y)
        (((2, 4), 1), ((1, 2), -4/5), ((0, 1), 1/3))
    """
    return tuple(poly.dict().items())

cpdef MPolynomial_libsingular tup_to_poly(tuple eq_tup, MPolynomialRing_libsingular parent):
    """
    Return a polynomial object from its tuple representation. Inverse of poly_to_tup.
    poly_to_tup(tup_to_poly(eq_tup, ring)) == eq_tup && tup_to_poly(poly_to_tup(eq), eq.parent()) == eq

    NOTE:
      Assumes parent.ngens() == len(exp_tup) for exp_tup, c in eq_tup

    EXAMPLES::

        sage: from sage.combinat.root_system.poly_tup_engine import tup_to_poly
        sage: R.<x,y> = PolynomialRing(CyclotomicField(20))
        sage: poly_tup = (((2,0),1), ((0,0),1))
        sage: tup_to_poly(poly_tup, parent=R)
        x^2 + 1
        sage: poly = x**2*y**4 - 4/5*x*y**2 + 1/3 * y
        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
        sage: tup_to_poly(poly_to_tup(poly), parent=R) == poly
        True

    TESTS:

        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup, tup_to_poly
        sage: R.<x,y,z> = PolynomialRing(CyclotomicField(20))
        sage: r = R.random_element()
        sage: tup_to_poly(poly_to_tup(r), parent=R) == r
        True
        sage: R = PolynomialRing(QQ, 'fx', 100)
        sage: r = R.random_element()
        sage: tup_to_poly(poly_to_tup(r), parent=R) == r
        True
    """
    # Maybe the following is faster but we need to ensure all coefficients are
    # already in fmats._poly_ring.base_ring() so that implicit casting is avoided
    # (this is pretty slow)
    # return parent._element_constructor_({ exp : c for exp, c in tup_of_pairs }, check=False)
    return parent(dict(eq_tup))

######################
### "Change rings" ###
######################

def tup_to_univ_poly(tuple eq_tup, gen):
    """
    Given a tuple of pairs representing a univariate polynomial and a univariate
    polynomial ring generator, return a univariate polynomial object.

    Each pair in the tuple is assumed to be of the form (ETuple, coeff), where
    coeff is an element of gen.parent().base_ring().

    EXAMPLES::

      sage: from sage.combinat.root_system.poly_tup_engine import tup_to_univ_poly
      sage: from sage.rings.polynomial.polydict import ETuple
      sage: poly_tup = ((ETuple([0,3,0]),2), (ETuple([0,1,0]),-1), (ETuple([0,0,0]),-2/3))
      sage: b = QQ['b'].gen()
      sage: tup_to_univ_poly(poly_tup, b)
      2*b^3 - b - 2/3
      sage: poly_tup = ((ETuple([0, 0, 0]), -1/5),)
      sage: tup_to_univ_poly(poly_tup, b)
      -1/5
    """
    univ_tup = tuple((exp.nonzero_values()[0] if exp.nonzero_values() else 0, c) for exp, c in eq_tup)
    return sum(c * gen ** p for p, c in univ_tup)

cpdef tuple resize(tuple eq_tup, dict idx_map, int nvars):
    """
    Return a tuple representing a polynomial in a ring with len(sorted_vars) generators
    This method is used for creating polynomial objects with the "right number" of
    variables for computing Groebner bases of the partitioned equations graph
    and for adding constraints ensuring certain F-symbols are nonzero

    EXAMPLES::

      sage: from sage.combinat.root_system.poly_tup_engine import resize
      sage: from sage.rings.polynomial.polydict import ETuple
      sage: poly_tup = ((ETuple([0,3,0,2]),2), (ETuple([0,1,0,1]),-1), (ETuple([0,0,0,0]),-2/3))
      sage: idx_map = { 1 : 0, 3 : 1 }
      sage: resize(poly_tup,idx_map,2)
      (((3, 2), 2), ((1, 1), -1), ((0, 0), -2/3))

      sage: R = PolynomialRing(ZZ, 'fx', 20)
      sage: R.inject_variables()
      Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13, fx14, fx15, fx16, fx17, fx18, fx19
      sage: sparse_poly = fx0**2 * fx17 + fx3
      sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup, tup_to_poly
      sage: S.<x,y,z> = PolynomialRing(ZZ)
      sage: tup_to_poly(resize(poly_to_tup(sparse_poly),{0:0,3:1,17:2},3), parent=S)
      x^2*z + y
    """
    cdef ETuple new_e
    cdef list resized = list()
    for exp, c in eq_tup:
      new_e = ETuple({ idx_map[pos] : d for pos, d in exp.sparse_iter() }, nvars)
      resized.append((new_e,c))
    return tuple(resized)

###########################
### Convenience methods ###
###########################

cdef ETuple degrees(tuple poly_tup):
    """
    Return the maximal degree of each variable in the polynomial
    """
    #Deal with the empty tuple, representing the zero polynomial
    if not poly_tup: return ETuple()
    cdef ETuple max_degs, exp
    max_degs = poly_tup[0][0]
    for exp, c in poly_tup[1:]:
        max_degs = max_degs.emax(exp)
    return max_degs

cpdef ETuple get_variables_degrees(list eqns):
    """
    Find maximum degrees for each variable in equations

    EXAMPLES::

        sage: from sage.combinat.root_system.poly_tup_engine import get_variables_degrees
        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: polys = [x**2 + 1, x*y*z**2 - 4*x*y, x*z**3 - 4/3*y + 1]
        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
        sage: get_variables_degrees([poly_to_tup(p) for p in polys])
        (2, 1, 3)
    """
    if not eqns: return ETuple([])
    cdef tuple eq_tup
    cdef ETuple max_deg
    max_deg = degrees(eqns[0])
    for eq_tup in eqns[1:]:
        max_deg = max_deg.emax(degrees(eq_tup))
    return max_deg

cpdef list variables(tuple eq_tup):
    """
    Return indices of all variables appearing in eq_tup

    EXAMPLES::

      sage: from sage.combinat.root_system.poly_tup_engine import variables
      sage: from sage.rings.polynomial.polydict import ETuple
      sage: poly_tup = ((ETuple([0,3,0]),2), (ETuple([0,1,0]),-1), (ETuple([0,0,0]),-2/3))
      sage: variables(poly_tup)
      [1]
      sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
      sage: R.<x,y,z> = PolynomialRing(QQ)
      sage: variables(poly_to_tup(x*2*y + y**3 - 4/3*x))
      [0, 1]
      sage: variables(poly_to_tup(R(1/4)))
      []
    """
    return degrees(eq_tup).nonzero_positions()

cpdef constant_coeff(tuple eq_tup):
    """
    Return the constant coefficient of the polynomial represented by given tuple

    EXAMPLES::

        sage: from sage.combinat.root_system.poly_tup_engine import constant_coeff
        sage: from sage.rings.polynomial.polydict import ETuple
        sage: poly_tup = ((ETuple([0,3,0]),2), (ETuple([0,1,0]),-1), (ETuple([0,0,0]),-2/3))
        sage: constant_coeff(poly_tup)
        -2/3
        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
        sage: constant_coeff(poly_to_tup(x**5 + x*y*z - 9))
        -9
    """
    for exp, coeff in eq_tup:
      if exp.is_constant():
        return coeff
    return 0

cpdef tuple apply_coeff_map(tuple eq_tup, coeff_map):
    """
    Apply `coeff_map` to coefficients

    EXAMPLES::

        sage: from sage.combinat.root_system.poly_tup_engine import apply_coeff_map
        sage: sq = lambda x : x**2
        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup, tup_to_poly
        sage: tup_to_poly(apply_coeff_map(poly_to_tup(x + 2*y + 3*z), sq), parent=R)
        x + 4*y + 9*z
    """
    cdef list new_tup = list()
    for exp, coeff in eq_tup:
      new_tup.append((exp, coeff_map(coeff)))
    return tuple(new_tup)

cpdef bint tup_fixes_sq(tuple eq_tup):
    """
    Determine if given equation fixes the square of a variable. An equation
    fixes the sq of a variable if it is of the form `a*x^2 + c` for *nonzero* constants `a`, `c`

    EXAMPLES::

        sage: from sage.combinat.root_system.poly_tup_engine import tup_fixes_sq
        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: sq_fixer = 3*z**2 - 1/8
        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
        sage: tup_fixes_sq(poly_to_tup(sq_fixer))
        True
        sage: tup_fixes_sq(poly_to_tup(y**3 + 2))
        False
        sage: tup_fixes_sq(poly_to_tup(x*y + 2))
        False
        sage: tup_fixes_sq(poly_to_tup(x**2 + y**2 + 2))
        False
    """
    #Make this faster by combining two conditions into one... don't create temp variables
    return len(eq_tup) == 2 and len(variables(eq_tup)) == 1 and eq_tup[0][0].nonzero_values() == [2]

######################
### Simplification ###
######################

cpdef dict subs_squares(dict eq_dict, dict known_sq):
    """
    Substitutes for known squares in a given polynomial.
    The parameter known_sq is a dictionary of (int i, NumberFieldElement a) pairs such that x_i^2 - a == 0
    Returns a dictionary of (ETuple, coeff) pairs representing polynomial

    EXAMPLES::

        sage: from sage.combinat.root_system.poly_tup_engine import subs_squares
        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: poly = x**2 + y**3 + x*z**3
        sage: known_sq = { 0 : 2, 1 : -1, 2 : -1/2 }
        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
        sage: subs_squares(dict(poly_to_tup(poly)), known_sq)
        {(0, 0, 0): 2, (0, 1, 0): -1, (1, 0, 1): -1/2}
    """
    cdef dict subbed, new_e
    cdef ETuple exp, lm
    cdef int idx, power
    subbed = dict()
    for exp, coeff in eq_dict.items():
        new_e = dict()
        for idx, power in exp.sparse_iter():
            if idx in known_sq:
                coeff *= known_sq[idx] ** (power // 2)
                #New power is 1 if power is odd
                if power & True:
                  new_e[idx] = 1
            else:
                new_e[idx] = power
        exp = ETuple(new_e, len(exp))
        #If exponent tuple is already present in dictionary, coefficients are added
        if exp in subbed:
            subbed[exp] += coeff
        else:
            subbed[exp] = coeff
    return subbed

cdef dict remove_gcf(dict eq_dict, ETuple nonz):
    """
    Returns a dictionary of (ETuple, coeff) pairs describing the polynomial eq / GCF(eq)
    The input nonz is an ETuple indicating the positions of variables known to be nonzero.
    The entries of nonz are assumed to be some relatively large number, like 100
    """
    #Find common variables, filtered according to known nonzeros
    cdef ETuple common_powers, exp
    common_powers = nonz
    for exp, c in eq_dict.items():
        common_powers = common_powers.emin(exp)
    cdef dict ret = dict()
    for exp, c in eq_dict.items():
        ret[exp.esub(common_powers)] = c
    return ret

cdef tuple to_monic(dict eq_dict):
    """
    Return tuple of pairs (ETuple, coeff) describing the monic polynomial associated to eq_dict
    Here, the leading coefficient is chosen according to the degree reverse lexicographic ordering
    (default for multivariate polynomial rings)
    """
    if not eq_dict: return tuple()
    cdef list ord_monoms = sorted(eq_dict, key=degrevlex_sortkey)
    cdef ETuple lm = ord_monoms[-1]
    lc = eq_dict[lm]
    if not lc: return tuple()
    # cdef list ret = [(lm, one)]
    cdef list ret = [(lm, lc.parent().one())]
    inv_lc = lc.inverse_of_unit()
    cdef ETuple exp
    for exp in reversed(ord_monoms[:-1]):
        ret.append((exp, inv_lc * eq_dict[exp]))
    return tuple(ret)

cdef tuple reduce_poly_dict(dict eq_dict, ETuple nonz, dict known_sq):
    """
    Return a dictionary describing a monic polynomial with no known nonzero gcd and
    no known squares
    """
    if not eq_dict: return tuple()
    return to_monic(remove_gcf(subs_squares(eq_dict, known_sq), nonz))

# cdef int common_denom(tuple eq_tup):
#     #Compute the common denominator
#     cdef list denoms = list()
#     cdef int common_denom
#     cdef ETuple exp
#     for exp, c in eq_tup:
#         denoms.append(c.denominator())
#     return LCM_list(denoms)
#
# cdef tuple integralify(tuple eq_tup):
#     if not eq_tup: return tuple()
#     cdef list ret = list()
#     cdef int cd = common_denom(eq_tup)
#     for exp, c in eq_tup:
#         ret.append((exp, c * cd))
#     return tuple(ret)
#
# cdef tuple reduce_poly_dict(dict eq_dict, ETuple nonz, dict known_sq):
#     """
#     Return a dictionary describing a monic polynomial with no known nonzero gcd and
#     no known squares
#     """
#     if not eq_dict: return tuple()
#     return integralify(to_monic(remove_gcf(subs_squares(eq_dict, known_sq), nonz)))

####################
### Substitution ###
####################

cpdef dict compute_known_powers(ETuple max_deg, dict val_dict):
    """
    Pre-compute powers of known values for efficiency when preparing to substitute
    into a list of polynomials.

    INPUTS:

        max_deg is an ETuple indicating the maximal degree of each variable
        val_dict is a dictionary of (var_idx, poly_tup) key-value pairs. poly_tup
        is a tuple of (ETuple, coeff) pairs reperesenting a multivariate polynomial

    EXAMPLES::

        sage: from sage.combinat.root_system.poly_tup_engine import compute_known_powers
        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: polys = [x**3 + 1, x**2*y + z**3, y**2 - 3*y]
        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
        sage: known_val = { 0 : poly_to_tup(R(-1)), 2 : poly_to_tup(y**2) }
        sage: from sage.combinat.root_system.poly_tup_engine import get_variables_degrees
        sage: max_deg = get_variables_degrees([poly_to_tup(p) for p in polys])
        sage: compute_known_powers(max_deg, known_val)
        {0: [(((0, 0, 0), 1),),
        (((0, 0, 0), -1),),
        (((0, 0, 0), 1),),
        (((0, 0, 0), -1),)],
        2: [(((0, 0, 0), 1),),
        (((0, 2, 0), 1),),
        (((0, 4, 0), 1),),
        (((0, 6, 0), 1),)]}
    """
    if not max_deg: return dict()
    assert max(max_deg.nonzero_values(sort=False)) <= 100, "NotImplementedError: Cannot substitute for degree larger than 100"
    max_deg = max_deg.emin(ETuple({ idx : 100 for idx in val_dict }, len(max_deg)))
    cdef dict known_powers
    #Get polynomial unit as tuple to initialize list elements
    cdef tuple one_tup = ((max_deg._new(), one),)
    cdef int d, power, var_idx
    known_powers = { var_idx : [one_tup]*(d+1) for var_idx, d in max_deg.sparse_iter() }
    for var_idx, d in max_deg.sparse_iter():
        for power in range(d):
            known_powers[var_idx][power+1] = tup_mul(known_powers[var_idx][power],val_dict[var_idx])
    return known_powers

cpdef dict subs(tuple poly_tup, dict known_powers):
    """
    Substitute given variables into a polynomial tuple

    EXAMPLES::

        sage: from sage.combinat.root_system.poly_tup_engine import subs
        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: polys = [x**3 + 1, x**2*y + z**3, y**2 - 3*y]
        sage: from sage.combinat.root_system.poly_tup_engine import compute_known_powers, get_variables_degrees, poly_to_tup
        sage: known_val = { 0 : poly_to_tup(R(-1)), 2 : poly_to_tup(y**2) }
        sage: max_deg = get_variables_degrees([poly_to_tup(p) for p in polys])
        sage: poly_tup = poly_to_tup(polys[0])
        sage: subs(poly_tup, compute_known_powers(max_deg, known_val))
        {(0, 0, 0): 0}
        sage: poly_tup = poly_to_tup(polys[1])
        sage: subs(poly_tup, compute_known_powers(max_deg, known_val))
        {(0, 1, 0): 1, (0, 6, 0): 1}
    """
    cdef dict subbed = {}
    cdef ETuple exp, m, shifted_exp
    cdef int var_idx, power
    cdef tuple temp
    for exp, coeff in poly_tup:
        #Get polynomial unit as tuple
        # temp = ((exp._new(), one),)
        temp = ((exp._new(), coeff.parent().one()),)
        for var_idx, power in exp.sparse_iter():
            if var_idx in known_powers:
                exp = exp.eadd_p(-power,var_idx)
                temp = tup_mul(temp,known_powers[var_idx][power])
        for m, c in temp:
            shifted_exp = exp.eadd(m)
            if shifted_exp in subbed:
                subbed[shifted_exp] += coeff*c
            else:
                subbed[shifted_exp] = coeff*c
    return subbed

cdef tuple tup_mul(tuple p1, tuple p2):
    """
    Multiplication of two tuples... may have to make this faster
    """
    cdef dict prod = dict()
    cdef ETuple xi, yj, shifted_exp
    for xi, ai in p1:
        for yj, bj in p2:
            shifted_exp = xi.eadd(yj)
            if shifted_exp in prod:
                prod[shifted_exp] += ai*bj
            else:
                prod[shifted_exp] = ai*bj
    return tuple(prod.items())

###############
### Sorting ###
###############

#Implement richcmp comparator object that can be passed in as key to sorted method

cpdef int poly_tup_cmp(tuple tleft, tuple tright):
    """
    Determine which polynomial is larger with respect to the degrevlex ordering

    EXAMPLES::

        sage: from sage.combinat.root_system.poly_tup_engine import poly_tup_cmp
        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: p1 = x*y*z - x**2 + 3/2
        sage: p2 = x*y*z - x * y +1/2
        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
        sage: (p1 < p2) == (poly_tup_cmp(poly_to_tup(p1), poly_to_tup(p2)) < 0)
        True
        sage: R.<x,y,z> = PolynomialRing(CyclotomicField(20))
        sage: zeta20 = R.base_ring().gen()
        sage: p1 = zeta20**2 * x*z**2 - 2*zeta20
        sage: p2 = y**3 + 1/4
        sage: (p1 < p2) == (poly_tup_cmp(poly_to_tup(p1), poly_to_tup(p2)) < 0)
        True

    TESTS:

        sage: from sage.combinat.root_system.poly_tup_engine import poly_tup_cmp, poly_to_tup
        sage: R.<x,y,z> = PolynomialRing(CyclotomicField(20))
        sage: p1 = R.random_element()
        sage: p2 = R.random_element()
        sage: (p1 < p2) == (poly_tup_cmp(poly_to_tup(p1), poly_to_tup(p2)) < 0)
        True
        sage: (p1 > p2) == (poly_tup_cmp(poly_to_tup(p1), poly_to_tup(p2)) > 0)
        True
        sage: poly_tup_cmp(poly_to_tup(p1), poly_to_tup(p1)) == 0
        True
    """
    cdef int i, ret, sf, sg, val
    cdef ETuple f, g
    ret = 0
    for i in range(min(len(tleft),len(tright))):
        f, g = tleft[i][0], tright[i][0]
        if f == g:
            if tleft[i][1] != tright[i][1]:
                ret = -1 + 2*(tleft[i][1] > tright[i][1])
        else:
            sf, sg = 0, 0
            for val in f.nonzero_values(sort=False):
              sf += val
            for val in g.nonzero_values(sort=False):
              sg += val
            ret = -1 + 2*(sf > sg or ( sf == sg and f.reversed() < g.reversed() ))
        if ret != 0:
            return ret
    return len(tleft) - len(tright)
