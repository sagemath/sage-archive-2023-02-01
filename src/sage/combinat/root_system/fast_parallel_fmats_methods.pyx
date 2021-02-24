############################
### Fast FMatrix methods ###
############################
from itertools import product
from os import getpid
from sage.combinat.root_system.poly_tup_engine cimport *
from sage.rings.ideal import Ideal
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
cimport cython

#Define a global temporary worker results repository
worker_results = list()

###############
### Mappers ###
###############

#Given an FMatrix factory and a sextuple, return a pentagon equation as a polynomial object
cdef req_cy(factory, tuple sextuple, side="left"):
    a, b, c, d, e, g = sextuple
    #To add typing we need to ensure all fmats.fmat are of the same type?
    #Return fmats._poly_ring.zero() and fmats._poly_ring.one() instead of 0 and 1?
    lhs = factory.FR.r_matrix(a,c,e)*factory.fmat(a,c,b,d,e,g)*factory.FR.r_matrix(b,c,g)
    rhs = 0
    for f in factory.FR.basis():
      rhs += factory.fmat(c,a,b,d,e,f)*factory.FR.r_matrix(f,c,d)*factory.fmat(a,b,c,d,f,g)
    return lhs-rhs

#Set up and reduce the hexagon equations corresponding to this worker
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef get_reduced_hexagons(factory, int child_id, int n_proc):
  global worker_results
  cdef int i
  cdef tuple sextuple
  for i, sextuple in enumerate(product(factory.FR.basis(),repeat=6)):
      if i % n_proc == child_id:
          he = req_cy(factory,sextuple)
          if he:
              worker_results.append(reduce_poly_dict(he.dict(),factory._nnz,factory._ks))

#Given an FMatrix factory and a nonuple, return a pentagon equation as a polynomial object
cdef MPolynomial_libsingular feq_cy(factory, tuple nonuple, bint prune=False):
    a, b, c, d, e, f, g, k, l = nonuple
    cdef lhs = factory.fmat(f,c,d,e,g,l)*factory.fmat(a,b,l,e,f,k)
    if lhs == 0 and prune: # it is believed that if lhs=0, the equation carries no new information
      return factory._poly_ring.zero()
    cdef rhs = factory._poly_ring.zero()
    for h in factory.FR.basis():
      rhs += factory.fmat(a,b,c,g,f,h)*factory.fmat(a,h,d,e,g,k)*factory.fmat(b,c,d,k,h,l)
    return lhs - rhs

#Set up and reduce the pentagon equations corresponding to this worker
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef get_reduced_pentagons(factory, int child_id, int n_proc, bint prune=True):
  cdef int i
  cdef tuple nonuple
  cdef MPolynomial_libsingular pe
  for i, nonuple in enumerate(product(factory.FR.basis(),repeat=9)):
      if i % n_proc == child_id:
          pe = feq_cy(factory,nonuple,prune=prune)
          if pe:
              worker_results.append(reduce_poly_dict(pe.dict(),factory._nnz,factory._ks))

#Substitute known values, known squares, and reduce!
cpdef update_reduce(factory, tuple eq_tup):
    global worker_results
    cdef dict eq_dict = subs(eq_tup,factory._kp)
    cdef reduced
    if tup_fixes_sq(tuple(eq_dict.items())):
        reduced = to_monic(eq_dict)
    else:
        reduced = reduce_poly_dict(eq_dict,factory._nnz,factory._ks)
    worker_results.append(reduced)

#Compute Groebner basis for given equations iterable
cpdef compute_gb(factory,eqns,term_order="degrevlex"):
    global worker_results
    #Define smaller poly ring in component vars
    sorted_vars = list()
    for eq_tup in eqns:
      for fx in variables(eq_tup):
        sorted_vars.append(fx)
    sorted_vars = sorted(set(sorted_vars))
    R = PolynomialRing(factory.FR.field(),len(sorted_vars),'a',order=term_order)

    #Zip tuples into R and compute Groebner basis
    cdef idx_map = dict()
    for new, old in enumerate(sorted_vars):
        idx_map[old] = new
    nvars = len(sorted_vars)
    cdef list polys = list()
    for eq_tup in eqns:
        polys.append(tup_to_poly(resize(eq_tup,idx_map,nvars),parent=R))
    gb = Ideal(sorted(polys)).groebner_basis(algorithm="libsingular:slimgb")

    #Change back to fmats poly ring and append to temp_eqns
    cdef dict inv_idx_map = dict()
    for k, v in idx_map.items():
        inv_idx_map[v] = k
    nvars = factory._poly_ring.ngens()
    for p in gb:
        worker_results.append(resize(poly_to_tup(p),inv_idx_map,nvars))

################
### Reducers ###
################

#Helper function for returning processed results back to parent process
#Trivial reducer: simply collects objects with the same key in the worker
def collect_eqns(proc):
    #Discard the zero polynomial
    global worker_results
    reduced = set(worker_results)-set([tuple()])
    worker_results = list()
    return list(reduced)

####################
### Verification ###
####################

#Check the pentagon equation corresponding to the given nonuple
cdef feq_verif(factory, tuple nonuple, float tol=5e-8):
    global worker_results
    a, b, c, d, e, f, g, k, l = nonuple
    cdef float diff, lhs, rhs
    lhs = factory.fmat(f,c,d,e,g,l)*factory.fmat(a,b,l,e,f,k)
    rhs = 0.0
    for h in factory.FR.basis():
      rhs += factory.fmat(a,b,c,g,f,h)*factory.fmat(a,h,d,e,g,k)*factory.fmat(b,c,d,k,h,l)
    diff = lhs - rhs
    if diff > tol or diff < -tol:
      worker_results.append(diff)

#Generate all the pentagon equations assigned to this process, and reduce them
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef pent_verify(tuple mp_params, factory):
    child_id, n_proc = mp_params
    cdef float t0
    cdef tuple nonuple
    cdef long i
    for i, nonuple in enumerate(product(factory.FR.basis(),repeat=9)):
        if i % n_proc == child_id:
          feq_verif(factory,nonuple)
        if i % 50000000 == 0 and i:
          print("{:5d}m equations checked... {} potential misses so far...".format(i // 1000000,len(worker_results)))
    print("Ran through {} pentagons... Worker {} with pid {} reporting {} potential misses...".format(i,child_id,getpid(),len(worker_results)))
