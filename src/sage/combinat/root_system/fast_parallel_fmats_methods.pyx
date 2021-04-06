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
from sage.combinat.root_system.poly_tup_engine cimport (
    compute_known_powers,
    get_variables_degrees, variables,
    poly_to_tup, _tup_to_poly,
    subs, subs_squares, reduce_poly_dict, resize,
    _flatten_coeffs, _unflatten_coeffs,
    has_appropriate_linear_term,
    resize
)
from sage.rings.number_field.number_field_element cimport NumberFieldElement_absolute
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular, MPolynomialRing_libsingular
from sage.rings.polynomial.polydict cimport ETuple

import ctypes
from itertools import product
from sage.rings.ideal import Ideal
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

##########################
### Fast class methods ###
##########################

cpdef _solve_for_linear_terms(factory,eqns=None):
    """
    Solve for a linear term occurring in a two-term equation, and for
    variables appearing in univariate single-term equations.

    EXAMPLES::

        sage: f = FMatrix(FusionRing("D3",1), inject_variables=True)
        creating variables fx1..fx27
        Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13, fx14, fx15, fx16, fx17, fx18, fx19, fx20, fx21, fx22, fx23, fx24, fx25, fx26
        sage: f.ideal_basis = [fx0**3, fx0 + fx3**4, fx2**2 - fx3, fx2 - fx3**2, fx4 - fx2]
        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
        sage: f.ideal_basis = [poly_to_tup(p) for p in f.ideal_basis]
        sage: from sage.combinat.root_system.poly_tup_engine import poly_tup_sortkey
        sage: f.ideal_basis.sort(key=poly_tup_sortkey)
        sage: f._fvars = {sextuple : poly_to_tup(rhs) for sextuple, rhs in f._fvars.items()}
        sage: from sage.combinat.root_system.fast_parallel_fmats_methods import _solve_for_linear_terms
        sage: _solve_for_linear_terms(f)
        True
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[0]])
        0
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[2]])
        fx4
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[3]])
        fx3
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[4]])
        fx4
    """
    if eqns is None:
            eqns = factory.ideal_basis

    linear_terms_exist = False
    for eq_tup in eqns:
        #Only unflatten relevant polynomials
        if len(eq_tup) > 2:
            continue
        eq_tup = _unflatten_coeffs(factory._field, eq_tup)

        if len(eq_tup) == 1:
            vars = variables(eq_tup)
            if len(vars) == 1 and vars[0] not in factory._solved:
                factory._fvars[factory._idx_to_sextuple[vars[0]]] = tuple()
                factory._solved.add(vars[0])
                linear_terms_exist = True
        if len(eq_tup) == 2:
            idx = has_appropriate_linear_term(eq_tup)
            if idx < 0: continue
            #The chosen term is guaranteed to be univariate in the largest variable
            max_var = eq_tup[idx][0].nonzero_positions()[0]
            if max_var not in factory._solved:
                rhs_exp = eq_tup[(idx+1) % 2][0]
                rhs_coeff = -eq_tup[(idx+1) % 2][1] / eq_tup[idx][1]
                factory._fvars[factory._idx_to_sextuple[max_var]] = ((rhs_exp,rhs_coeff),)
                factory._solved.add(max_var)
                linear_terms_exist = True
    return linear_terms_exist

cpdef _backward_subs(factory):
    """
    Perform backward substitution on ``self.ideal_basis``, traversing
    variables in reverse lexicographical order.

    EXAMPLES::

        sage: f = FMatrix(FusionRing("D3",1), inject_variables=True)
        creating variables fx1..fx27
        Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13, fx14, fx15, fx16, fx17, fx18, fx19, fx20, fx21, fx22, fx23, fx24, fx25, fx26
        sage: f.ideal_basis = [fx0**3, fx0 + fx3**4, fx2**2 - fx3, fx2 - fx3**2, fx4 - fx2]
        sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
        sage: f.ideal_basis = [poly_to_tup(p) for p in f.ideal_basis]
        sage: from sage.combinat.root_system.poly_tup_engine import poly_tup_sortkey
        sage: f.ideal_basis.sort(key=poly_tup_sortkey)
        sage: f._fvars = {sextuple : poly_to_tup(rhs) for sextuple, rhs in f._fvars.items()}
        sage: from sage.combinat.root_system.fast_parallel_fmats_methods import _solve_for_linear_terms
        sage: _solve_for_linear_terms(f)
        True
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[0]])
        0
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[2]])
        fx4
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[3]])
        fx3
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[4]])
        fx4
        sage: f.ideal_basis.append(poly_to_tup(fx4-1))
        sage: f.ideal_basis.sort(key=poly_tup_sortkey)
        sage: _solve_for_linear_terms(f)
        True
        sage: from sage.combinat.root_system.fast_parallel_fmats_methods import _backward_subs
        sage: _backward_subs(f)
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[0]])
        0
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[2]])
        1
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[3]])
        fx3
        sage: f._tup_to_fpoly(f._fvars[f._idx_to_sextuple[4]])
        1
    """
    one = factory._field.one()
    for var in reversed(factory._poly_ring.gens()):
        sextuple = factory._var_to_sextuple[var]
        rhs = factory._fvars[sextuple]
        d = {var_idx : factory._fvars[factory._idx_to_sextuple[var_idx]] for var_idx in variables(rhs) if var_idx in factory._solved}
        if d:
            kp = compute_known_powers(get_variables_degrees([rhs]), d, one)
            factory._fvars[sextuple] = tuple(subs_squares(subs(rhs,kp,one),factory._ks).items())

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

#Define a global temporary worker results repository
cdef list worker_results = list()

cdef _fmat(fvars, _Nk_ij, id_anyon, a, b, c, d, x, y):
      """
      Cython version of fmat class method. Using cdef for fastest dispatch
      """
      if _Nk_ij(a,b,x) == 0 or _Nk_ij(x,c,d) == 0 or _Nk_ij(b,c,y) == 0 or _Nk_ij(a,y,d) == 0:
              return 0
      #Some known F-symbols
      if a == id_anyon:
          if x == b and y == d:
              return 1
          else:
              return 0
      if b == id_anyon:
          if x == a and y == c:
              return 1
          else:
              return 0
      if c == id_anyon:
          if x == d and y == b:
              return 1
          else:
              return 0
      return fvars[a,b,c,d,x,y]

cdef req_cy(tuple basis, r_matrix, dict fvars, Nk_ij, id_anyon, tuple sextuple):
    """
    Given an FMatrix factory and a sextuple, return a hexagon equation as a polynomial object
    """
    a, b, c, d, e, g = sextuple
    #To add typing we need to ensure all fmats.fmat are of the same type?
    #Return fmats._poly_ring.zero() and fmats._poly_ring.one() instead of 0 and 1?
    lhs = r_matrix(a,c,e,base_coercion=False)*_fmat(fvars,Nk_ij,id_anyon,a,c,b,d,e,g)*r_matrix(b,c,g,base_coercion=False)
    rhs = 0
    for f in basis:
      rhs += _fmat(fvars,Nk_ij,id_anyon,c,a,b,d,e,f)*r_matrix(f,c,d,base_coercion=False)*_fmat(fvars,Nk_ij,id_anyon,a,b,c,d,f,g)
    return lhs-rhs

@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef get_reduced_hexagons(factory, tuple mp_params):
    """
    Set up and reduce the hexagon equations corresponding to this worker
    """
    #Set up multiprocessing parameters
    global worker_results
    cdef child_id, n_proc
    cdef unsigned long i
    child_id, n_proc = mp_params
    cdef tuple sextuple, red

    #Pre-compute common parameters for speed
    cdef tuple basis = tuple(factory._FR.basis())
    cdef dict fvars = factory._fvars
    r_matrix = factory._FR.r_matrix
    _Nk_ij = factory._FR.Nk_ij
    id_anyon = factory._FR.one()
    cdef NumberFieldElement_absolute one = factory._field.one()
    cdef ETuple _nnz = factory._nnz
    cdef dict _ks = factory._ks

    #Computation loop
    for i, sextuple in enumerate(product(basis,repeat=6)):
        if i % n_proc == child_id:
            he = req_cy(basis,r_matrix,fvars,_Nk_ij,id_anyon,sextuple)
            if he:
                red = reduce_poly_dict(he.dict(),_nnz,_ks,one)

                #Avoid pickling cyclotomic coefficients
                red = _flatten_coeffs(red)

                worker_results.append(red)

cdef MPolynomial_libsingular feq_cy(tuple basis, fvars, Nk_ij, id_anyon, zero, tuple nonuple, bint prune=False):
    """
    Given an FMatrix factory and a nonuple, return a pentagon equation as a polynomial object
    """
    a, b, c, d, e, f, g, k, l = nonuple
    cdef lhs = _fmat(fvars,Nk_ij,id_anyon,f,c,d,e,g,l)*_fmat(fvars,Nk_ij,id_anyon,a,b,l,e,f,k)
    if lhs == 0 and prune: # it is believed that if lhs=0, the equation carries no new information
        return zero
    cdef rhs = zero
    for h in basis:
      rhs += _fmat(fvars,Nk_ij,id_anyon,a,b,c,g,f,h)*_fmat(fvars,Nk_ij,id_anyon,a,h,d,e,g,k)*_fmat(fvars,Nk_ij,id_anyon,b,c,d,k,h,l)
    return lhs - rhs

@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef get_reduced_pentagons(factory, tuple mp_params):
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
    _Nk_ij = factory._FR.Nk_ij
    id_anyon = factory._FR.one()
    cdef NumberFieldElement_absolute one = factory._field.one()
    cdef ETuple _nnz = factory._nnz
    cdef dict _ks = factory._ks
    cdef MPolynomial_libsingular zero = factory._poly_ring.zero()

    #Computation loop
    it = product(basis,repeat=9)
    for i in range(len(basis)**9):
        nonuple = next(it)
        if i % n_proc == child_id:
            pe = feq_cy(basis,fvars,_Nk_ij,id_anyon,zero,nonuple,prune=True)
            if pe:
                red = reduce_poly_dict(pe.dict(),_nnz,_ks,one)

                #Avoid pickling cyclotomic coefficients
                red = _flatten_coeffs(red)

                worker_results.append(red)

cpdef update_reduce(factory, tuple eq_tup):
    """
    Substitute known values, known squares, and reduce!

    EXAMPLES::

        sage: f = FMatrix(FusionRing("B4",1), fusion_label="b", inject_variables=True)
        creating variables fx1..fx14
        Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
        sage: p = 2*fx0**2 - 3/7*fx0*fx1*fx3**3 + 5*fx0*fx4
        sage: from sage.combinat.root_system.poly_tup_engine import *
        sage: f.ideal_basis = [poly_to_tup(p)]
        sage: f._ks = {3: 2}
        sage: f._var_degs = get_variables_degrees(f.ideal_basis)
        sage: known_vals = {4: poly_to_tup(fx11**2)}
        sage: f._nnz
        (100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 0, 0, 0, 0)
        sage: f._kp = compute_known_powers(f._var_degs,known_vals,f._field.one())
        sage: from sage.combinat.root_system.fast_parallel_fmats_methods import update_reduce, collect_eqns
        sage: update_reduce(f, f.ideal_basis[0])
        sage: from sage.combinat.root_system.poly_tup_engine import _unflatten_coeffs
        sage: [f._tup_to_fpoly(_unflatten_coeffs(f._field,p)) for p in collect_eqns(0)]
        [fx1*fx3 - 35/6*fx11^2 - 7/3*fx0]
    """
    global worker_results
    cdef NumberFieldElement_absolute one = factory._field.one()

    #Construct cyclotomic field elts from list repn
    cdef tuple unflat = _unflatten_coeffs(factory._field, eq_tup)

    cdef dict eq_dict = subs(unflat,factory._kp,one)
    cdef tuple red = reduce_poly_dict(eq_dict,factory._nnz,factory._ks,one)

    #Avoid pickling cyclotomic coefficients
    red = _flatten_coeffs(red)

    worker_results.append(red)

cpdef compute_gb(factory, tuple args):
    """
    Compute the reduced Groebner basis for given equations iterable.

    EXAMPLES::

        sage: f = FMatrix(FusionRing("A2",2))
        sage: from multiprocessing import Pool, set_start_method
        sage: try:
        ....:     set_start_method('fork') # context can be set only once
        ....: except RuntimeError:
        ....:     pass
        sage: f.get_defining_equations('hexagons',worker_pool=Pool(),output=False)
        sage: partition = f._partition_eqns()
        Partitioned 261 equations into 57 components of size:
        [24, 12, 12, 12, 12, 12, 12, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1]
        sage: from sage.combinat.root_system.fast_parallel_fmats_methods import compute_gb, collect_eqns
        sage: args = (partition[(29, 40, 83, 111, 148, 154)], "degrevlex")
        sage: compute_gb(f, args)
        sage: from sage.combinat.root_system.poly_tup_engine import _unflatten_coeffs
        sage: [f._tup_to_fpoly(_unflatten_coeffs(f._field,t)) for t in collect_eqns(0)]
        [fx83*fx154 - fx111,
         fx29^2*fx154 + fx29,
         fx29*fx148 - fx40,
         fx29*fx154^2 + fx154,
         fx83*fx148 + fx29*fx154,
         fx40*fx83 - fx29,
         fx111*fx148 - fx154,
         fx40*fx154 + fx148,
         fx40*fx111 - fx29*fx154,
         fx29*fx111 + fx83]
    """
    global worker_results
    cdef list eqns, sorted_vars
    eqns, term_order = args
    #Define smaller poly ring in component vars
    sorted_vars = list()
    cdef tuple eq_tup
    cdef int fx
    for eq_tup in eqns:
      for fx in variables(eq_tup):
        sorted_vars.append(fx)
    sorted_vars = sorted(set(sorted_vars))
    cdef MPolynomialRing_libsingular R = PolynomialRing(factory._FR.field(),len(sorted_vars),'a',order=term_order)

    #Zip tuples into R and compute Groebner basis
    cdef idx_map = { old : new for new, old in enumerate(sorted_vars) }
    nvars = len(sorted_vars)
    F = factory.field()
    cdef list polys = list()
    for eq_tup in eqns:
        eq_tup = _unflatten_coeffs(F, eq_tup)
        polys.append(_tup_to_poly(resize(eq_tup,idx_map,nvars),parent=R))
    gb = Ideal(sorted(polys)).groebner_basis(algorithm="libsingular:slimgb")

    #Change back to fmats poly ring and append to temp_eqns
    cdef dict inv_idx_map = { v : k for k, v in idx_map.items() }
    cdef tuple t
    nvars = factory._poly_ring.ngens()
    for p in gb:
        t = resize(poly_to_tup(p),inv_idx_map,nvars)

        #Avoid pickling cyclotomic coefficients
        t = _flatten_coeffs(t)

        worker_results.append(t)

cpdef update_child_fmats(factory, tuple data_tup):
    """
    One-to-all communication used to update FMatrix object after each triangular
    elim step. We must update the algorithm's state values. These are:
    _fvars, _solved, _ks, _var_degs, _nnz, and _kp.

    TESTS::

        sage: f = FMatrix(FusionRing("A1",3))
        sage: f.get_orthogonality_constraints(output=False)
        sage: from multiprocessing import Pool, set_start_method
        sage: try:
        ....:     set_start_method('fork') # context can be set only once
        ....: except RuntimeError:
        ....:     pass
        sage: pool = Pool()
        sage: f.get_defining_equations('hexagons',worker_pool=pool,output=False)
        sage: f.ideal_basis = f._par_graph_gb(worker_pool=pool,verbose=False)
        sage: from sage.combinat.root_system.poly_tup_engine import poly_tup_sortkey
        sage: f.ideal_basis.sort(key=poly_tup_sortkey)
        sage: f.mp_thresh = 0
        sage: f._triangular_elim(worker_pool=pool)          # indirect doctest
        Elimination epoch completed... 10 eqns remain in ideal basis
        Elimination epoch completed... 0 eqns remain in ideal basis
        sage: f.ideal_basis
        []
    """
    #factory object is assumed global before forking used to create the Pool object,
    #so each child has a global fmats variable. So it's enough to update that object
    factory._fvars, factory._solved, factory._ks, factory._var_degs = data_tup
    factory._nnz = factory._get_known_nonz()
    factory._kp = compute_known_powers(factory._var_degs,factory._get_known_vals(),factory._field.one())

################
### Reducers ###
################

cpdef collect_eqns(proc):
    """
    Helper function for returning processed results back to parent process.

    Trivial reducer: simply collects objects with the same key in the worker.
    This method is only useful when called after :meth:`executor`, whose
    function argument appends output to the ``worker_results`` list.

    EXAMPLES::

        sage: from sage.combinat.root_system.fast_parallel_fmats_methods import executor
        sage: fmats = FMatrix(FusionRing("A1",3))
        sage: params = (('get_reduced_hexagons', id(fmats)), (0,1))
        sage: executor(params)
        sage: from sage.combinat.root_system.fast_parallel_fmats_methods import collect_eqns
        sage: len(collect_eqns(0)) == 63
        True
        sage: fmats = FMatrix(FusionRing("G2",2))
        sage: params = (('get_reduced_pentagons', id(fmats)), (0,1))
        sage: executor(params)
        sage: len(collect_eqns(0)) == 4911
        True
    """
    #Discard the zero polynomial
    global worker_results
    reduced = set(worker_results)-set([tuple()])
    worker_results = list()
    return list(reduced)

##############################
### Parallel code executor ###
##############################

#Hard-coded module __dict__-style attribute with visible cdef methods
cdef dict mappers = {
    "get_reduced_hexagons": get_reduced_hexagons,
    "get_reduced_pentagons": get_reduced_pentagons,
    "update_reduce": update_reduce,
    "compute_gb": compute_gb,
    "update_child_fmats": update_child_fmats,
    "pent_verify": pent_verify
    }

cpdef executor(params):
    r"""
    Execute a function defined in this module
    (``sage.combinat.root_system.fast_parallel_fmats_methods``) in a worker
    process, and supply the factory parameter by constructing a reference
    to the ``FMatrix`` object in the worker's memory adress space from
    its ``id``.

    INPUT:

    - ``params`` -- a tuple ``((fn_name, fmats_id), fn_args)`` where
      ``fn_name`` is the name of the function to be executed, ``fmats_id``
      is the ``id`` of the :class:`FMatrix` object, and ``fn_args`` is a
      tuple containing all arguments to be passed to the function ``fn_name``.

    .. NOTES::

        When the parent process is forked, each worker gets a copy of
        every  global variable. The virtual memory address of object `X` in
        the parent process equals the *virtual* memory address of the copy
        of object `X` in each worker, so we may construct references to
        forked copies of `X` using an ``id`` obtained in the parent process.

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
    #Bind module method to FMatrix object in worker process, and call the method
    return mappers[fn_name](fmats_obj,args)

####################
### Verification ###
####################

cdef feq_verif(factory, fvars, Nk_ij, id_anyon, tuple nonuple, float tol=5e-8):
    """
    Check the pentagon equation corresponding to the given nonuple
    """
    global worker_results
    a, b, c, d, e, f, g, k, l = nonuple
    cdef float diff, lhs, rhs

    lhs = _fmat(fvars,Nk_ij,id_anyon,f,c,d,e,g,l)*_fmat(fvars,Nk_ij,id_anyon,a,b,l,e,f,k)
    rhs = 0.0
    for h in factory._FR.basis():
      rhs += _fmat(fvars,Nk_ij,id_anyon,a,b,c,g,f,h)*_fmat(fvars,Nk_ij,id_anyon,a,h,d,e,g,k)*_fmat(fvars,Nk_ij,id_anyon,b,c,d,k,h,l)
    diff = lhs - rhs
    if diff > tol or diff < -tol:
      worker_results.append(diff)

@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef pent_verify(factory, tuple mp_params):
    """
    Generate all the pentagon equations assigned to this process, and reduce them
    """
    child_id, n_proc, verbose = mp_params
    cdef float t0
    cdef tuple nonuple
    cdef unsigned long long i

    #Pre-compute common parameters for speed
    Nk_ij = factory._FR.Nk_ij
    cdef dict fvars = factory._fvars
    id_anyon = factory._FR.one()
    for i, nonuple in enumerate(product(factory._FR.basis(),repeat=9)):
        if i % n_proc == child_id:
          feq_verif(factory,fvars,Nk_ij,id_anyon,nonuple)
        if i % 50000000 == 0 and i and verbose:
          print("{:5d}m equations checked... {} potential misses so far...".format(i // 1000000,len(worker_results)))
