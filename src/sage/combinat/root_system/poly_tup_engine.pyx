###################################################
### Arithmetic Engine for polynomials as tuples ###
###################################################
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
from sage.rings.rational_field import RationalField as QQ

###########
### API ###
###########

cpdef tuple poly_to_tup(MPolynomial_libsingular poly):
    """
    Convert a polynomial object into the internal representation as tuple of
    (ETuple exp, NumberFieldElement coeff) pairs
    """
    return tuple(poly.dict().items())

cpdef MPolynomial_libsingular tup_to_poly(tuple eq_tup, MPolynomialRing_libsingular parent):
    """
    Return a polynomial object from its tuple representation. Inverse of poly_to_tup.
    poly_to_tup(tup_to_poly(eq_tup, ring)) == eq_tup && tup_to_poly(poly_to_tup(eq), eq.parent()) == eq
    Assumes parent.ngens() == len(exp_tup) for exp_tup, c in eq_tup
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
    polynomial ring generator, return a univariate polynomial object
    """
    univ_tup = tuple((exp.nonzero_values()[0] if exp.nonzero_values() else 0, c) for exp, c in eq_tup)
    return sum(c * gen ** p for p, c in univ_tup)

cpdef tuple resize(tuple eq_tup, dict idx_map, int nvars):
    """
    Return a tuple representing a polynomial in a ring with len(sorted_vars) generators
    This method is used for creating polynomial objects with the "right number" of
    variables for computing Groebner bases of the partitioned equations graph
    and for adding constraints ensuring certain F-symbols are nonzero
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
    """
    cdef tuple eq_tup
    cdef ETuple max_deg
    max_deg = degrees(eqns[0])
    for eq_tup in eqns[1:]:
        max_deg = max_deg.emax(degrees(eq_tup))
    return max_deg

cpdef list variables(tuple eq_tup):
    """
    Return indices of all variables appearing in eq_tup
    """
    return degrees(eq_tup).nonzero_positions()

cpdef constant_coeff(tuple eq_tup):
    """
    Return the constant coefficient of the polynomial represented by given tuple
    """
    for exp, coeff in eq_tup:
      if exp.is_constant():
        return coeff
    return 0

cpdef tuple apply_coeff_map(tuple eq_tup, coeff_map):
    """
    Apply `coeff_map` to coefficients
    """
    cdef list new_tup = list()
    for exp, coeff in eq_tup:
      new_tup.append((exp, coeff_map(coeff)))
    return tuple(new_tup)

cpdef bint tup_fixes_sq(tuple eq_tup):
    """
    Determine if given equation fixes the square of a variable. An equation
    fixes the sq of a variable if it is of the form `a*x^2 + c` for *nonzero* constants `a`, `c`
    """
    return len(eq_tup) == 2 and len(variables(eq_tup)) == 1 and eq_tup[0][0].nonzero_values() == [2]

######################
### Simplification ###
######################

cpdef dict subs_squares(dict eq_dict, dict known_sq):
    """
    Substitutes for known squares in a given polynomial.
    The parameter known_sq is a dictionary of (int i, NumberFieldElement a) pairs such that x_i^2 - a == 0
    Returns a dictionary of (ETuple, coeff) pairs representing polynomial
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
    return { exp : a for exp, a in subbed.items() }

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
    return { exp.esub(common_powers) : c for exp, c in eq_dict.items() }

cdef tuple to_monic(dict eq_dict):
    """
    Return tuple of pairs (ETuple, coeff) describing the monic polynomial associated to eq_dict
    Here, the leading coefficient is chosen according to the degree reverse lexicographic ordering
    (default for multivariate polynomial rings)
    """
    if not eq_dict: return tuple()
    it = reversed(sorted(eq_dict, key=TermOrder().sortkey_degrevlex))
    lm = next(it)
    lc = eq_dict[lm]
    if not lc: return tuple()
    cdef ETuple exp
    cdef list ret = [(lm, QQ().one())]
    inv_lc = lc.inverse_of_unit()
    for exp in it:
        ret.append((exp, inv_lc * eq_dict[exp]))
    return tuple(ret)

cdef tuple reduce_poly_dict(dict eq_dict, ETuple nonz, dict known_sq):
    """
    Return a dictionary describing a monic polynomial with no known nonzero gcd and
    no known squares
    """
    if not eq_dict: return tuple()
    return to_monic(remove_gcf(subs_squares(eq_dict, known_sq), nonz))

####################
### Substitution ###
####################

#Pre-compute powers of known values for efficiency
cpdef dict compute_known_powers(ETuple max_deg, dict val_dict):
    """
    Pre-compute powers of known values for efficiency
    """
    assert max(max_deg.nonzero_values(sort=False)) <= 100, "NotImplementedError: Cannot substitute for degree larger than 100"
    max_deg = max_deg.emin(ETuple({ idx : 100 for idx in val_dict }, len(max_deg)))
    cdef dict known_powers
    #Get polynomial unit as tuple to initialize list elements
    cdef tuple one = tuple(((ETuple({},len(max_deg)),QQ().one()),))
    cdef int d, var_idx
    known_powers = { var_idx : [one]*(d+1) for var_idx, d in max_deg.sparse_iter() }
    for var_idx, d in max_deg.sparse_iter():
        for power in range(d):
            known_powers[var_idx][power+1] = tup_mul(known_powers[var_idx][power],val_dict[var_idx])
    return known_powers

cpdef dict subs(tuple poly_tup, dict known_powers):
    """
    Substitute given variables into a polynomial tuple
    """
    cdef dict subbed = {}
    cdef ETuple exp, m
    cdef int var_idx, power
    cdef tuple temp
    for exp, coeff in poly_tup:
        #Get polynomial unit as tuple
        temp = tuple(((ETuple({},len(exp)),QQ().one()),))
        for var_idx, power in exp.sparse_iter():
            if var_idx in known_powers:
                exp = exp.eadd_p(-power,var_idx)
                temp = tup_mul(temp,known_powers[var_idx][power])
        for m, c in temp:
            if exp.eadd(m) in subbed:
                subbed[exp.eadd(m)] += coeff*c
            else:
                subbed[exp.eadd(m)] = coeff*c
    return subbed

cdef tuple tup_mul(tuple p1, tuple p2):
    """
    Multiplication of two tuples... may have to make this faster
    """
    cdef dict prod = dict()
    cdef ETuple xi, yj
    for xi, ai in p1:
        for yj, bj in p2:
            if xi.eadd(yj) in prod:
                prod[xi.eadd(yj)] += ai*bj
            else:
                prod[xi.eadd(yj)] = ai*bj
    return tuple(prod.items())

###############
### Sorting ###
###############

#Implement richcmp comparator object that can be passed in as key to sorted method

cpdef int poly_tup_cmp(tuple tleft, tuple tright):
    """
    Determine which polynomial is larger with respect to the degrevlex ordering
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
