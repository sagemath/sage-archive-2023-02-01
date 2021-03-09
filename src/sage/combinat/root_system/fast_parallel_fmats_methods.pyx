"""
Fast FMatrix methods
"""
# ****************************************************************************
#  Copyright (C) 2021 Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport cython
import ctypes
from itertools import product
from os import getpid
import sage
from sage.combinat.root_system.poly_tup_engine cimport *
from sage.rings.ideal import Ideal
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.qqbar import number_field_elements_from_algebraics

#Define a global temporary worker results repository
cdef list worker_results = list()

##############################
### Parallel code executor ###
##############################

def executor(params):
    """
    Execute a function defined in this module (sage.combinat.root_system.fast_parallel_fmats_methods)
    in a worker process, and supply the factory parameter by constructing a reference
    to the FMatrix object in the worker's memory adress space from its id
    ##  NOTE: When the parent process is forked, each worker gets a copy of
    every  global variable. The virtual memory address of object X in the parent
    process equals the VIRTUAL memory address of the copy of object X in each
    worker, so we may construct references to forked copies of X using an id
    obtained in the parent process

    TESTS:

        sage: from sage.combinat.root_system.fast_parallel_fmats_methods import executor
        sage: fmats = FMatrix(FusionRing("A1",3))
        sage: params = (('get_reduced_hexagons', id(fmats)), (0,1))
        sage: executor(params)
        sage: from sage.combinat.root_system.fast_parallel_fmats_methods import collect_eqns
        sage: len(collect_eqns(0)) == 63
        True
        sage: fmats = FMatrix(FusionRing("E8",2))
        sage: params = (('get_reduced_hexagons', id(fmats)), (0,1))
        sage: executor(params)
        sage: len(collect_eqns(0)) == 11
        True
    """
    (fn_name, fmats_id), args = params
    #Construct a reference to global FMatrix object in this worker's memory
    fmats_obj = ctypes.cast(fmats_id, ctypes.py_object).value
    mod = sage.combinat.root_system.fast_parallel_fmats_methods
    #Bind module method to FMatrix object in worker process, and call the method
    return getattr(mod,fn_name)(fmats_obj,args)


######################################
### Fast fusion coefficients cache ###
######################################

# from sage.misc.cachefunc import cached_function
# cdef dict _Nk_ij = dict()

# cpdef _Nk_ij(factory,proc):
#     cdef int coeff
#     for a,b,c in product(factory._FR.basis(),repeat=3):
#         try:
#             coeff = (a*b).monomial_coefficients(copy=False)[c.weight()]
#         except:
#             coeff = 0
#         _Nk_ij[a,b,c] = coeff

# cpdef int _Nk_ij(a,b,c):
#     try:
#         return (a*b).monomial_coefficients(copy=False)[c.weight()]
#     except KeyError:
#         return 0
#
# _Nk_ij = cached_function(_Nk_ij, name='_Nk_ij')


###############
### Mappers ###
###############

cdef _fmat(fvars, _Nk_ij, one, a, b, c, d, x, y):
      """
      Cython version of fmat class method. Using cdef for fastest dispatch
      """
      if _Nk_ij(a,b,x) == 0 or _Nk_ij(x,c,d) == 0 or _Nk_ij(b,c,y) == 0 or _Nk_ij(a,y,d) == 0:
              return 0
      #Some known F-symbols
      if a == one:
          if x == b and y == d:
              return 1
          else:
              return 0
      if b == one:
          if x == a and y == c:
              return 1
          else:
              return 0
      if c == one:
          if x == d and y == b:
              return 1
          else:
              return 0
      return fvars[a,b,c,d,x,y]

cdef req_cy(tuple basis, r_matrix, dict fvars, Nk_ij, one, tuple sextuple):
    """
    Given an FMatrix factory and a sextuple, return a hexagon equation as a polynomial object
    """
    a, b, c, d, e, g = sextuple
    #To add typing we need to ensure all fmats.fmat are of the same type?
    #Return fmats._poly_ring.zero() and fmats._poly_ring.one() instead of 0 and 1?
    lhs = r_matrix(a,c,e)*_fmat(fvars,Nk_ij,one,a,c,b,d,e,g)*r_matrix(b,c,g)
    rhs = 0
    for f in basis:
      rhs += _fmat(fvars,Nk_ij,one,c,a,b,d,e,f)*r_matrix(f,c,d)*_fmat(fvars,Nk_ij,one,a,b,c,d,f,g)
    return lhs-rhs

@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef get_reduced_hexagons(factory, tuple mp_params):
    """
    Set up and reduce the hexagon equations corresponding to this worker
    """
    #Set up multiprocessing parameters
    global worker_results
    cdef child_id, n_proc
    cdef unsigned long i
    child_id, n_proc = mp_params
    cdef tuple sextuple

    #Pre-compute common parameters for speed
    cdef tuple basis = tuple(factory._FR.basis())
    cdef dict fvars = factory._fvars
    r_matrix = factory._FR.r_matrix
    _Nk_ij = factory._FR.Nk_ij
    one = factory._FR.one()
    cdef ETuple _nnz = factory._nnz
    cdef dict _ks = factory._ks

    #Computation loop
    for i, sextuple in enumerate(product(basis,repeat=6)):
        if i % n_proc == child_id:
            he = req_cy(basis,r_matrix,fvars,_Nk_ij,one,sextuple)
            if he:
                worker_results.append(reduce_poly_dict(he.dict(),_nnz,_ks))

cdef MPolynomial_libsingular feq_cy(tuple basis, fvars, Nk_ij, one, zero, tuple nonuple, bint prune=False):
    """
    Given an FMatrix factory and a nonuple, return a pentagon equation as a polynomial object
    """
    a, b, c, d, e, f, g, k, l = nonuple
    cdef lhs = _fmat(fvars,Nk_ij,one,f,c,d,e,g,l)*_fmat(fvars,Nk_ij,one,a,b,l,e,f,k)
    if lhs == 0 and prune: # it is believed that if lhs=0, the equation carries no new information
        return zero
    cdef rhs = zero
    for h in basis:
      rhs += _fmat(fvars,Nk_ij,one,a,b,c,g,f,h)*_fmat(fvars,Nk_ij,one,a,h,d,e,g,k)*_fmat(fvars,Nk_ij,one,b,c,d,k,h,l)
    return lhs - rhs

@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef get_reduced_pentagons(factory, tuple mp_params, bint prune=True):
    """
    Set up and reduce the pentagon equations corresponding to this worker
    """
    #Set up multiprocessing parameters
    global worker_results
    cdef int child_id, n_proc
    child_id, n_proc = mp_params
    cdef unsigned long i
    cdef tuple nonuple, red
    cdef MPolynomial_libsingular pe

    #Pre-compute common parameters for speed
    cdef tuple basis = tuple(factory._FR.basis())
    cdef dict fvars = factory._fvars
    r_matrix = factory._FR.r_matrix
    _Nk_ij = factory._FR.Nk_ij
    one = factory._FR.one()
    cdef ETuple _nnz = factory._nnz
    cdef dict _ks = factory._ks
    cdef MPolynomial_libsingular zero = factory._poly_ring.zero()

    #Computation loop
    for i, nonuple in enumerate(product(basis,repeat=9)):
        if i % n_proc == child_id:
            pe = feq_cy(basis,fvars,_Nk_ij,one,zero,nonuple,prune=prune)
            if pe:
                red = reduce_poly_dict(pe.dict(),_nnz,_ks)
                worker_results.append(red)

cpdef update_reduce(factory, tuple eq_tup):
    """
    Substitute known values, known squares, and reduce!
    """
    global worker_results
    cdef dict eq_dict = subs(eq_tup,factory._kp)
    cdef tuple red = reduce_poly_dict(eq_dict,factory._nnz,factory._ks)
    worker_results.append(red)

cpdef compute_gb(factory, tuple args):
    """
    Compute Groebner basis for given equations iterable
    """
    global worker_results
    eqns, term_order = args
    #Define smaller poly ring in component vars
    sorted_vars = list()
    for eq_tup in eqns:
      for fx in variables(eq_tup):
        sorted_vars.append(fx)
    sorted_vars = sorted(set(sorted_vars))
    R = PolynomialRing(factory._FR.field(),len(sorted_vars),'a',order=term_order)

    #Zip tuples into R and compute Groebner basis
    cdef idx_map = { old : new for new, old in enumerate(sorted_vars) }
    nvars = len(sorted_vars)
    cdef list polys = list()
    for eq_tup in eqns:
        polys.append(tup_to_poly(resize(eq_tup,idx_map,nvars),parent=R))
    gb = Ideal(sorted(polys)).groebner_basis(algorithm="libsingular:slimgb")

    #Change back to fmats poly ring and append to temp_eqns
    cdef dict inv_idx_map = { v : k for k, v in idx_map.items() }
    cdef tuple t
    nvars = factory._poly_ring.ngens()
    for p in gb:
        t = resize(poly_to_tup(p),inv_idx_map,nvars)
        worker_results.append(t)

cpdef update_child_fmats(factory, tuple data_tup):
    """
    One-to-all communication used to update fvars after triangular elim step.
    """
    #factory object is assumed global before forking used to create the Pool object,
    #so each child has a global fmats variable. So it's enough to update that object
    factory._fvars, factory.solved, factory._ks, factory._var_degs = data_tup
    factory._nnz = factory.get_known_nonz()
    factory._kp = compute_known_powers(factory._var_degs,factory.get_known_vals())

def get_appropriate_number_field(factory,non_cyclotomic_roots):
    """
    If needed, find a NumberField containing the roots given roots
    """
    roots = [factory._FR.field().gen()]+[r[1] for r in non_cyclotomic_roots]
    return number_field_elements_from_algebraics(roots,minimal=True)

################
### Reducers ###
################

def collect_eqns(proc):
    """
    Helper function for returning processed results back to parent process.
    Trivial reducer: simply collects objects with the same key in the worker
    """
    #Discard the zero polynomial
    global worker_results
    reduced = set(worker_results)-set([tuple()])
    worker_results = list()
    return list(reduced)

# ####################
# ### Verification ###
# ####################
#
# cdef feq_verif(factory, tuple nonuple, float tol=5e-8):
#     """
#     Check the pentagon equation corresponding to the given nonuple
#     """
#     global worker_results
#     a, b, c, d, e, f, g, k, l = nonuple
#     cdef float diff, lhs, rhs
#     lhs = _fmat(factory,f,c,d,e,g,l)*_fmat(factory,a,b,l,e,f,k)
#     rhs = 0.0
#     for h in factory._FR.basis():
#       rhs += _fmat(factory,a,b,c,g,f,h)*_fmat(factory,a,h,d,e,g,k)*_fmat(factory,b,c,d,k,h,l)
#     diff = lhs - rhs
#     if diff > tol or diff < -tol:
#       worker_results.append(diff)
#
# @cython.wraparound(False)
# @cython.nonecheck(False)
# @cython.cdivision(True)
# cpdef pent_verify(factory, tuple mp_params):
#     """
#     Generate all the pentagon equations assigned to this process, and reduce them
#     """
#     child_id, n_proc = mp_params
#     cdef float t0
#     cdef tuple nonuple
#     cdef unsigned long i
#     for i, nonuple in enumerate(product(factory._FR.basis(),repeat=9)):
#         if i % n_proc == child_id:
#           feq_verif(factory,nonuple)
#         if i % 50000000 == 0 and i:
#           print("{:5d}m equations checked... {} potential misses so far...".format(i // 1000000,len(worker_results)))
#     print("Ran through {} pentagons... Worker {} with pid {} reporting {} potential misses...".format(i,child_id,getpid(),len(worker_results)))
