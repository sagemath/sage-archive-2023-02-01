"""
Fast FusionRing methods for computing braid group representations
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
import sage
from sage.combinat.root_system.fast_parallel_fmats_methods cimport _fmat
from sage.misc.cachefunc import cached_function

#Define a global temporary worker results repository
worker_results = list()

##############################
### Parallel code executor ###
##############################

def executor(params):
    """
    Execute a function defined in this module (sage.combinat.root_system.fast_parallel_fusion_ring_braid_repn)
    in a worker process, and supply the factory parameter by constructing a reference
    to the FMatrix object in the worker's memory adress space from its id
    ##  NOTE: When the parent process is forked, each worker gets a copy of
    every  global variable. The virtual memory address of object X in the parent
    process equals the VIRTUAL memory address of the copy of object X in each
    worker, so we may construct references to forked copies of X
    """
    (fn_name, fr_id), args = params
    #Construct a reference to global FMatrix object in this worker's memory
    fusion_ring_obj = ctypes.cast(fr_id, ctypes.py_object).value
    mod = sage.combinat.root_system.fast_parallel_fusion_ring_braid_repn
    #Bind module method to FMatrix object in worker process, and call the method
    getattr(mod,fn_name)(fusion_ring_obj,args)

###############
### Mappers ###
###############

cpdef mid_sig_ij(fusion_ring,row,col,a,b):
    """
    Compute the (xi,yi), (xj,yj) entry of generator braiding the middle two strands
    in the tree b -> xi # yi -> (a # a) # (a # a), which results in a sum over j
    of trees b -> xj # yj -> (a # a) # (a # a)
    """
    xi, yi = row
    xj, yj = col
    entry = 0
    phi = fusion_ring.fmats.get_coerce_map_from_fr_cyclotomic_field()
    for c in fusion_ring.basis():
        for d in fusion_ring.basis():
            ##Warning: We assume F-matrices are orthogonal!!! (using transpose for inverse)
            f1 = _fmat(fusion_ring.fmats,a,a,yi,b,xi,c)
            f2 = _fmat(fusion_ring.fmats,a,a,a,c,d,yi)
            f3 = _fmat(fusion_ring.fmats,a,a,a,c,d,yj)
            f4 = _fmat(fusion_ring.fmats,a,a,yj,b,xj,c)
            r = fusion_ring.r_matrix(a,a,d)
            if not phi.is_identity():
                r = phi(r)
            entry += f1 * f2 * r * f3 * f4
    return entry

cpdef odd_one_out_ij(fusion_ring,xi,xj,a,b):
    """
    Compute the xi, xj entry of the braid generator on the right-most strands,
    corresponding to the tree b -> (xi # a) -> (a # a) # a, which results in a
    sum over j of trees b -> xj -> (a # a) # (a # a)
    """
    entry = 0
    phi = fusion_ring.fmats.get_coerce_map_from_fr_cyclotomic_field()
    for c in fusion_ring.basis():
        ##Warning: We assume F-matrices are orthogonal!!! (using transpose for inverse)
        f1 = _fmat(fusion_ring.fmats,a,a,a,b,xi,c)
        f2 = _fmat(fusion_ring.fmats,a,a,a,b,xj,c)
        r = fusion_ring.r_matrix(a,a,c)
        if not phi.is_identity():
            r = phi(r)
        entry += f1 * r * f2
    return entry

#Cache methods
mid_sig_ij = cached_function(mid_sig_ij, name='mid_sig_ij')
odd_one_out_ij = cached_function(odd_one_out_ij, name='odd_one_out_ij')

@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef sig_2k(fusion_ring, tuple args):
    """
    Compute entries of the 2k-th braid generator
    """
    child_id, n_proc, fn_args = args
    k, a, b, n_strands = fn_args
    cdef int ctr = -1
    global worker_results
    #Get computational basis
    cdef list comp_basis = fusion_ring.get_comp_basis(a,b,n_strands)
    cdef dict basis_dict = { elt : i for i, elt in enumerate(comp_basis) }
    cdef int dim = len(comp_basis)
    cdef set entries = set()
    cdef int i
    for i in range(dim):
        for f,e,q in product(fusion_ring.basis(),repeat=3):
            #Distribute work amongst processes
            ctr += 1
            if ctr % n_proc != child_id:
                continue

            #Compute appropriate possible nonzero row index
            nnz_pos = list(comp_basis[i])
            nnz_pos[k-1:k+1] = f,e
            #Handle the special case k = 1
            if k > 1:
                nnz_pos[n_strands//2+k-2] = q
            nnz_pos = tuple(nnz_pos)

            #Skip repeated entries when k = 1
            if nnz_pos in comp_basis and (basis_dict[nnz_pos],i) not in entries:
                m, l = comp_basis[i][:n_strands//2], comp_basis[i][n_strands//2:]
                #A few special cases
                top_left = m[0]
                if k >= 3:
                    top_left = l[k-3]
                root = b
                if k - 1 < len(l):
                    root = l[k-1]

                #Handle the special case k = 1
                if k == 1:
                    entry = mid_sig_ij(fusion_ring,m[:2],(f,e),a,root)
                    worker_results.append(((basis_dict[nnz_pos],i), entry))
                    continue

                entry = 0
                for p in fusion_ring.basis():
                    f1 = _fmat(fusion_ring.fmats,top_left,m[k-1],m[k],root,l[k-2],p)
                    f2 = _fmat(fusion_ring.fmats,top_left,f,e,root,q,p)
                    entry += f1 * mid_sig_ij(fusion_ring,(m[k-1],m[k]),(f,e),a,p) * f2
                worker_results.append(((basis_dict[nnz_pos],i), entry))
                entries.add((basis_dict[nnz_pos],i))

@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef odd_one_out(fusion_ring, tuple args):
    """
    Compute entries of the rightmost braid generator, in case we have an odd number
    of strands
    """
    global worker_results
    child_id, n_proc, fn_args = args
    a, b, n_strands = fn_args
    cdef ctr = -1
    #Get computational basis
    comp_basis = fusion_ring.get_comp_basis(a,b,n_strands)
    basis_dict = { elt : i for i, elt in enumerate(comp_basis) }
    dim = len(comp_basis)
    for i in range(dim):
        for f, q in product(fusion_ring.basis(),repeat=2):
            #Distribute work amongst processes
            ctr += 1
            if ctr % n_proc != child_id:
                continue

            #Compute appropriate possible nonzero row index
            nnz_pos = list(comp_basis[i])
            nnz_pos[n_strands//2-1] = f
            #Handle small special case
            if n_strands > 3:
                nnz_pos[-1] = q
            nnz_pos = tuple(nnz_pos)

            if nnz_pos in comp_basis:
                m, l = comp_basis[i][:n_strands//2], comp_basis[i][n_strands//2:]

                #Handle a couple of small special cases
                if n_strands == 3:
                    entry = odd_one_out_ij(fusion_ring,m[-1],f,a,b)
                    worker_results.append(((basis_dict[nnz_pos],i), entry))
                    continue
                top_left = m[0]
                if n_strands > 5:
                    top_left = l[-2]
                root = b

                #Compute relevant entry
                entry = 0
                for p in fusion_ring.basis():
                    f1 = _fmat(fusion_ring.fmats,top_left,m[-1],a,root,l[-1],p)
                    f2 = _fmat(fusion_ring.fmats,top_left,f,a,root,q,p)
                    entry += f1 * odd_one_out_ij(fusion_ring,m[-1],f,a,p) * f2
                worker_results.append(((basis_dict[nnz_pos],i), entry))

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
