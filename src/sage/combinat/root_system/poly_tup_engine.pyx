from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular, MPolynomialRing_libsingular
from sage.rings.polynomial.polydict cimport ETuple, PolyDict
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.rational_field import RationalField as QQ

###NOTE: issues with importing NumberFieldElement_absolute and FMatrix
#from sage.rings.number_field.number_field_element cimport NumberFieldElement
#from sage.combinat.root_system.f_matrix import FMatrix

###################################################
### Arithmetic Engine for polynomials as tuples ###
###################################################

###########
### API ###
###########

#Convert a polynomial object into the internal representation as tuple of
#(ETuple exp, NumberFieldElement coeff) pairs
cpdef tuple poly_to_tup(MPolynomial_libsingular poly):
    return tuple(poly.dict().items())

#Return a polynomial object from its tuple representation. Inverse of poly_to_tup.
#poly_to_tup(tup_to_poly(eq_tup, ring)) == eq_tup && tup_to_poly(poly_to_tup(eq), eq.parent()) == eq
#Assumes parent.ngens() == len(exp_tup) for exp_tup, c in eq_tup
cpdef MPolynomial_libsingular tup_to_poly(tuple eq_tup, MPolynomialRing_libsingular parent):
    # Maybe the following is faster but we need to ensure all coefficients are
    # already in fmats._poly_ring.base_ring() so that implicit casting is avoided
    # (this is pretty slow)
    # return parent._element_constructor_({ exp : c for exp, c in tup_of_pairs }, check=False)
    return parent(dict(eq_tup))

######################
### "Change rings" ###
######################

#Given a tuple of pairs representing a univariate polynomial and a univariate
#polynomial ring generator, return a univariate polynomial object
def tup_to_univ_poly(eq_tup, gen):
    univ_tup = tuple((exp.nonzero_values()[0] if exp.nonzero_values() else 0, c) for exp, c in eq_tup)
    return sum(c * gen ** p for p, c in univ_tup)

#Return a tuple representing a polynomial in a ring with len(sorted_vars) generators
#This method is used for creating polynomial objects with the "right number" of
#variables for computing Groebner bases of the partitioned equations graph
#and for adding constraints ensuring certain F-symbols are nonzero
cpdef tuple resize(tuple eq_tup, dict idx_map, int nvars):
    cdef ETuple new_e
    cdef list resized = list()
    for exp, c in eq_tup:
      new_e = ETuple({ idx_map[pos] : d for pos, d in exp.sparse_iter() }, nvars)
      resized.append((new_e,c))
    return tuple(resized)

###########################
### Convenience methods ###
###########################

#Return the maximal degree of each variable in the polynomial
cdef ETuple degrees(tuple poly_tup):
    #Deal with the empty tuple, representing the zero polynomial
    if not poly_tup: return ETuple()
    cdef ETuple max_degs, exp
    max_degs = poly_tup[0][0]
    for exp, c in poly_tup[1:]:
        max_degs = max_degs.emax(exp)
    return max_degs

#Find maximum degrees for each variable in equations
cpdef ETuple get_variables_degrees(list eqns):
    cdef tuple eq_tup
    cdef ETuple max_deg
    max_deg = degrees(eqns[0])
    for eq_tup in eqns[1:]:
        max_deg = max_deg.emax(degrees(eq_tup))
    return max_deg

#Return indices of all variables appearing in eq_tup
cpdef list variables(tuple eq_tup):
    return degrees(eq_tup).nonzero_positions()

#Return the constant coefficient of the polynomial represented by given tuple
cpdef constant_coeff(tuple eq_tup):
    for exp, coeff in eq_tup:
      if exp.is_constant():
        return coeff
    return 0

cpdef tuple apply_coeff_map(tuple eq_tup, coeff_map):
    cdef list new_tup = list()
    for exp, coeff in eq_tup:
      new_tup.append((exp, coeff_map(coeff)))
    return tuple(new_tup)

#Determine if given equation fixes the square of a variable. An equation
#fixes the sq of a variable if it is of the form a*x^2 + c for NONZERO constants a, c
cpdef bint tup_fixes_sq(tuple eq_tup):
    return len(eq_tup) == 2 and len(variables(eq_tup)) == 1 and eq_tup[0][0].nonzero_values() == [2]

######################
### Simplification ###
######################

#Substitutes for known squares in a given polynomial.
#The parameter known_sq is a dictionary of (int i, NumberFieldElement a) pairs such that x_i^2 - a == 0
#Returns a dictionary of (ETuple, coeff) pairs representing polynomial
cpdef dict subs_squares(dict eq_dict, dict known_sq):
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

#Returns a dictionary of (ETuple, coeff) pairs describing the polynomial eq / GCF(eq)
#The input nonz is an ETuple indicating the positions of variables known to be nonzero.
#The entries of nonz are assumed to be some relatively large number, like 100
cdef dict remove_gcf(dict eq_dict, ETuple nonz):
    #Find common variables, filtered according to known nonzeros
    cdef ETuple common_powers, exp
    common_powers = nonz
    for exp, c in eq_dict.items():
        common_powers = common_powers.emin(exp)
    return { exp.esub(common_powers) : c for exp, c in eq_dict.items() }

#Return tuple of pairs (ETuple, coeff) describing the monic polynomial associated to eq_dict
#Here, the leading coefficient is chosen according to the degree reverse lexicographic ordering
#(default for multivariate polynomial rings)
cdef tuple to_monic(dict eq_dict):
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

# Return a dictionary describing a monic polynomial with no known nonzero gcd and
# no known squares
cdef tuple reduce_poly_dict(dict eq_dict, ETuple nonz, dict known_sq):
    if not eq_dict: return tuple()
    return to_monic(remove_gcf(subs_squares(eq_dict, known_sq), nonz))

####################
### Substitution ###
####################

#Pre-compute powers of known values for efficiency
cpdef dict compute_known_powers(ETuple max_deg, dict val_dict):
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

#Substitute given variables into a polynomial tuple
cpdef dict subs(tuple poly_tup, dict known_powers):
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

#Multiplication of two tuples... may have to make this faster
cdef tuple tup_mul(tuple p1, tuple p2):
    cdef dict prod = dict()
    cdef ETuple xi, yj
    for xi, ai in p1:
        for yj, bj in p2:
            if xi.eadd(yj) in prod:
                prod[xi.eadd(yj)] += ai*bj
            else:
                prod[xi.eadd(yj)] = ai*bj
    return tuple(prod.items())

#Substitute known values, known squares, and reduce!
cpdef update_reduce(tuple eq_tup, factory, mr_eng):
    cdef dict eq_dict = subs(eq_tup,factory._kp)
    cdef reduced
    if tup_fixes_sq(tuple(eq_dict.items())):
        reduced = to_monic(eq_dict)
    else:
        reduced = reduce_poly_dict(eq_dict,factory._nnz,factory._ks)
    mr_eng.worker_results.append(reduced)

###############
### Sorting ###
###############

#Determine which polynomial is larger with respect to the degrevlex ordering
cpdef int poly_tup_cmp(tuple tleft, tuple tright):
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

############################
### Fast FMatrix methods ###
############################
from itertools import product

#Given an FMatrix factory and a sextuple, return a pentagon equation as a polynomial object
cpdef tuple req_cy(factory, tuple sextuple, side="left"):
    a, b, c, d, e, g = sextuple
    #To add typing we need to ensure all fmats.fmat are of the same type?
    #Return fmats._poly_ring.zero() and fmats._poly_ring.one() instead of 0 and 1?
    lhs = factory.FR.r_matrix(a,c,e)*factory.fmat(a,c,b,d,e,g)*factory.FR.r_matrix(b,c,g)
    rhs = 0
    for f in factory.FR.basis():
      rhs += factory.fmat(c,a,b,d,e,f)*factory.FR.r_matrix(f,c,d)*factory.fmat(a,b,c,d,f,g)
    he = lhs - rhs
    if he:
        return reduce_poly_dict(he.dict(),factory._nnz,factory._ks)
    else:
        return tuple()

#Set up and reduce the hexagon equations corresponding to this worker
# cpdef get_reduced_hexagons(factory, int child_id, int n_proc):
#   cdef int i
#   cdef tuple sextuple
#   for i, sextuple in enumerate(product(factory.FR.basis(),repeat=6)):
#       if i % n_proc == child_id:
#           he = req_cy(factory,sextuple)
#           if he:
#               factory.temp_eqns.append(reduce_poly_dict(he.dict(),factory._nnz,factory._ks))

#Given an FMatrix factory and a nonuple, return a pentagon equation as a polynomial object

cpdef tuple feq_cy(factory, tuple nonuple, bint prune=False):
    a, b, c, d, e, f, g, k, l = nonuple
    cdef lhs = factory.fmat(f,c,d,e,g,l)*factory.fmat(a,b,l,e,f,k)
    if lhs == 0 and prune: # it is believed that if lhs=0, the equation carries no new information
      return factory._poly_ring.zero()
    cdef MPolynomial_libsingular rhs = factory._poly_ring.zero()
    for h in factory.FR.basis():
      rhs += factory.fmat(a,b,c,g,f,h)*factory.fmat(a,h,d,e,g,k)*factory.fmat(b,c,d,k,h,l)
    return reduce_poly_dict((lhs - rhs).dict(),factory._nnz,factory._ks)
  
#Set up and reduce the pentagon equations corresponding to this worker
# cpdef get_reduced_pentagons(factory, int child_id, int n_proc, bint prune=True):
#   cdef int i
#   cdef tuple nonuple
#   cdef MPolynomial_libsingular pe
#   for i, nonuple in enumerate(product(factory.FR.basis(),repeat=9)):
#       if i % n_proc == child_id:
#           pe = feq_cy(factory,nonuple,prune=prune)
#           if pe:
#               factory.temp_eqns.append(reduce_poly_dict(pe.dict(),factory._nnz,factory._ks))

####################
### Verification ###
####################

#Check the pentagon equation corresponding to the given nonuple
cpdef feq_verif(factory, tuple nonuple, float tol=5e-8):
    a, b, c, d, e, f, g, k, l = nonuple
    cdef float diff, lhs, rhs
    lhs = factory.fmat(f,c,d,e,g,l)*factory.fmat(a,b,l,e,f,k)
    rhs = 0.0
    for h in factory.FR.basis():
      rhs += factory.fmat(a,b,c,g,f,h)*factory.fmat(a,h,d,e,g,k)*factory.fmat(b,c,d,k,h,l)
    diff = lhs - rhs
    if diff > tol or diff < -tol:
      factory.temp_eqns.append(diff)

#Generate all the pentagon equations assigned to this process, and reduce them
from time import time
from os import getpid
cimport cython
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef pent_verify(tuple mp_params, factory):
    child_id, n_proc = mp_params
    cdef float t0
    cdef tuple nonuple
    cdef long i
    t0 = time()
    for i, nonuple in enumerate(product(factory.FR.basis(),repeat=9)):
        if i % n_proc == child_id:
          feq_verif(factory,nonuple)
        if i % 50000000 == 0 and i:
          print("{:5d}m equations checked in {:8.2f}... {} misses so far...".format(i // 1000000,time()-t0,len(factory.temp_eqns)))
    print("Ran through {} pentagons in {:8.2f}... Worker {} with {} reporting {} total misses...".format(i,time()-t0,child_id,getpid(),len(factory.temp_eqns)))

################################
### Well, this was a bust... ###
################################

#Given a tuple representation of an n-variate polynomial over a cyclotomic field,
#return a tuple representation of an (n+1)-variate polynomial over the rationals,
#with the cyclotomic field generator treated as an indeterminate and used as the
#last variable in the parent polynomial ring
cpdef tuple rationalify(factory, tuple eq_tup):
    F = factory.FR.field()
    cdef list rat_tup = list()
    cdef ETuple exp
    for exp, coeff in eq_tup:
      #In the future we should avoid casting by ensuring all coeffs created are
      #in fmats.FR.field (eg. cython reducers, FMatrix.fmat)
      for cyc_power, rat_coeff in F(coeff).polynomial().dict().items():
          nnz = { len(exp) : cyc_power }
          for idx, val in exp.sparse_iter():
              nnz[idx] = val
          new_e = ETuple(nnz, len(exp)+1)
          rat_tup.append((new_e, rat_coeff))
    return tuple(rat_tup)

#Given a tuple representation of an n-variate polynomial over the rationals,
#return a tuple repn of the polynomial with integer coefficients, obtained by
#multiplying all coefficients by the least common multiple of denominators
# cpdef tuple integralify(tuple rat_tup):
#     cdef int cdenom = LCM(coeff.denominator() for exp, coeff in rat_tup)
#     cdef list int_tup = list()
#     for exp, coeff in rat_tup:
#         int_tup.append((exp, int(cdenom * coeff)))
#     return tuple(int_tup)

#Left inverse of rationalify...
#cylcotomify(factory, rationalify(factory, eq_tup)) == eq_tup
cpdef tuple cyclotomify(factory, tuple rat_tup):
    gen = factory.FR.field().gen()
    cdef cyc_dict = dict()
    cdef ETuple exp
    for exp, coeff in rat_tup:
      cyc_pow = exp.get_exp(len(exp)-1)
      new_e = ETuple({ pos : val for pos, val in exp.sparse_iter() if pos != len(exp)-1 }, len(exp)-1)
      if new_e in cyc_dict:
          cyc_dict[new_e] += gen ** cyc_pow * coeff
      else:
          cyc_dict[new_e] = gen ** cyc_pow * coeff
    return tuple(cyc_dict.items())

# def clear_denom(eq):
#     for var in eq.variables():
#         if int(str(var)[1:]) >= fmats._poly_ring.ngens():
#             d = eq.degree(var)
#             eq *= var ** d
#     return eq
