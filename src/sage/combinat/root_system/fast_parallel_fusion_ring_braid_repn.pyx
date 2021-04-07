"""
Fast FusionRing methods for computing braid group representations
"""
# ****************************************************************************
#  Copyright (C) 2021 Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import ctypes
cimport cython
from sage.combinat.root_system.fast_parallel_fmats_methods cimport _fmat

from itertools import product
from sage.misc.cachefunc import cached_function
from sage.rings.qqbar import QQbar

#Define a global temporary worker results repository
worker_results = list()

###############
### Mappers ###
###############

cpdef mid_sig_ij(fusion_ring,row,col,a,b):
    """
    Compute the (xi,yi), (xj,yj) entry of generator braiding the middle two strands
    in the tree b -> xi # yi -> (a # a) # (a # a), which results in a sum over j
    of trees b -> xj # yj -> (a # a) # (a # a)

    .. WARNING::

        This method assumes F-matrices are orthogonal

    EXAMPLES::

        sage: from sage.combinat.root_system.fast_parallel_fusion_ring_braid_repn import mid_sig_ij
        sage: FR = FusionRing("A1",4)   # long time
        sage: FR.fusion_labels(['idd','one','two','three','four'],inject_variables=True)  # long time
        sage: one.weight()                                                                # long time
        (1/2, -1/2)
        sage: FR.get_computational_basis(one,two,4)                                       # long time
        [(two, two), (two, idd), (idd, two)]
        sage: FR.fmats.find_orthogonal_solution(verbose=False)                            # long time (~7s)
        sage: mid_sig_ij(FR, (two, two), (two, idd), one, two)                            # long time
        1/3*zeta48^10 - 2/3*zeta48^2

    This method works for all possible types of fields returned by
    ``self.fmats.field()``.

    TESTS::

        sage: FR = FusionRing("A1",3)                                                     # long time
        sage: FR.fusion_labels("a",inject_variables=True)                                 # long time
        sage: FR.fmats.find_orthogonal_solution(verbose=False)                            # long time
        sage: _, _, to_opt = FR.fmats.field().optimized_representation()                  # long time
        sage: a2**4                                                                       # long time
        2*a0 + 3*a2
        sage: FR.get_computational_basis(a2,a2,4)                                         # long time
        [(a2, a2), (a2, a0), (a0, a2)]
        sage: to_opt(mid_sig_ij(FR,(a2, a0),(a2, a2),a2,a2))                              # long time
        -2024728666370660589/10956322398441542221*a1^30 - 34142146914395596291/21912644796883084442*a1^28 - 21479437628091413631/21912644796883084442*a1^26 + 260131910217202103829/21912644796883084442*a1^24 + 69575612911670713183/10956322398441542221*a1^22 + 25621808994337724689/1992058617898462222*a1^20 - 1975139725303994650417/21912644796883084442*a1^18 - 1315664901396537703585/21912644796883084442*a1^16 - 2421451803369354765026/10956322398441542221*a1^14 - 5963323855935165859057/21912644796883084442*a1^12 - 4477124943233705460859/21912644796883084442*a1^10 - 2001454824483021618178/10956322398441542221*a1^8 - 2120319455379289595185/21912644796883084442*a1^6 - 15722612944437234961/755608441271830498*a1^4 - 39862668562651453480/10956322398441542221*a1^2 - 6967145776903524195/10956322398441542221
        sage: FR = FusionRing("G2",2)                             # long time
        sage: FR.fusion_labels("g",inject_variables=True)         # long time
        sage: FR.get_computational_basis(g1,g2,4)                 # long time
        [(g3, g2), (g3, g1), (g2, g3), (g2, g0), (g1, g3), (g1, g1), (g0, g2)]
        sage: FR.fmats.find_orthogonal_solution(verbose=False)    # long time (~11 s)
        sage: mid_sig_ij(FR,(g2, g3),(g1, g1),g1,g2)              # long time
        -0.4566723195695565? + 0.0805236512828312?*I
    """
    #Pre-compute common parameters for efficiency
    _fvars = fusion_ring.fmats._fvars
    _Nk_ij = fusion_ring.Nk_ij
    one = fusion_ring.one()

    xi, yi = row
    xj, yj = col
    entry = 0
    for c in fusion_ring.basis():
        for d in fusion_ring.basis():
            ##Warning: We assume F-matrices are orthogonal!!! (using transpose for inverse)
            f1 = _fmat(_fvars,_Nk_ij,one,a,a,yi,b,xi,c)
            f2 = _fmat(_fvars,_Nk_ij,one,a,a,a,c,d,yi)
            f3 = _fmat(_fvars,_Nk_ij,one,a,a,a,c,d,yj)
            f4 = _fmat(_fvars,_Nk_ij,one,a,a,yj,b,xj,c)
            r = fusion_ring.r_matrix(a,a,d)
            entry += f1 * f2 * r * f3 * f4
    return entry

cpdef odd_one_out_ij(fusion_ring,xi,xj,a,b):
    r"""
    Compute the `xi`, `xj` entry of the braid generator on the two right-most
    strands, corresponding to the tree b -> (xi # a) -> (a # a) # a, which
    results in a sum over j of trees b -> xj -> (a # a) # (a # a)

    .. WARNING::

        This method assumes F-matrices are orthogonal

    EXAMPLES::

        sage: from sage.combinat.root_system.fast_parallel_fusion_ring_braid_repn import odd_one_out_ij
        sage: FR = FusionRing("A1",4)                                                       # long time
        sage: FR.fusion_labels(["one","two","three","four","five"],inject_variables=True)   # long time
        sage: FR.get_computational_basis(two,two,5)                                         # long time
        [(three, three, one),
         (three, three, three),
         (three, one, three),
         (one, three, three),
         (one, one, one)]
        sage: FR.fmats.find_orthogonal_solution(verbose=False)                              # long time
        sage: odd_one_out_ij(FR,one,three,two,two)                                          # long time
        2/3*zeta48^12 - 1/3*zeta48^8 - 1/3*zeta48^4 - 1/3

    This method works for all possible types of fields returned by
    ``self.fmats.field()``.

    TESTS::

        sage: FR = FusionRing("A1",3)                                       # long time
        sage: FR.fusion_labels("a",inject_variables=True)                   # long time
        sage: FR.fmats.find_orthogonal_solution(verbose=False)              # long time
        sage: _, _, to_opt = FR.fmats.field().optimized_representation()    # long time
        sage: a2**3                                                         # long time
        a0 + 2*a2
        sage: FR.get_computational_basis(a2,a2,3)                           # long time
        [(a2,), (a0,)]
        sage: to_opt(odd_one_out_ij(FR,a0,a2,a2,a2))                        # long time
        6341990144855406911/21912644796883084442*a1^30 + 47313529044493641571/21912644796883084442*a1^28 - 6964289120109414595/10956322398441542221*a1^26 - 406719371329322780627/21912644796883084442*a1^24 + 87598732372849355687/10956322398441542221*a1^22 - 456724726845194775/19723352652460022*a1^20 + 3585892725441116840515/21912644796883084442*a1^18 - 645866255979227573282/10956322398441542221*a1^16 + 7958479159087829772639/21912644796883084442*a1^14 + 789748976956837633826/10956322398441542221*a1^12 + 3409710648897945752185/21912644796883084442*a1^10 + 903956381582048110980/10956322398441542221*a1^8 + 192973084151342020307/21912644796883084442*a1^6 - 9233312083438019435/755608441271830498*a1^4 + 667869266552877781/10956322398441542221*a1^2 + 17644302696056968099/21912644796883084442
        sage: FR = FusionRing("G2",2)                             # long time
        sage: FR.fusion_labels("g",inject_variables=True)         # long time
        sage: FR.fmats.find_orthogonal_solution(verbose=False)    # long time (~11 s)
        sage: odd_one_out_ij(FR,g1,g2,g1,g1)                      # long time
        -0.2636598866349343? + 0.4566723195695565?*I
    """
    #Pre-compute common parameters for efficiency
    _fvars = fusion_ring.fmats._fvars
    _Nk_ij = fusion_ring.Nk_ij
    one = fusion_ring.one()

    entry = 0
    for c in fusion_ring.basis():
        ##Warning: We assume F-matrices are orthogonal!!! (using transpose for inverse)
        f1 = _fmat(_fvars,_Nk_ij,one,a,a,a,b,xi,c)
        f2 = _fmat(_fvars,_Nk_ij,one,a,a,a,b,xj,c)
        r = fusion_ring.r_matrix(a,a,c)
        entry += f1 * r * f2
    return entry

#Cache methods
mid_sig_ij = cached_function(mid_sig_ij, name='mid_sig_ij')
odd_one_out_ij = cached_function(odd_one_out_ij, name='odd_one_out_ij')

@cython.nonecheck(False)
@cython.cdivision(True)
cdef sig_2k(fusion_ring, tuple args):
    """
    Compute entries of the 2k-th braid generator
    """
    #Pre-compute common parameters for efficiency
    _fvars = fusion_ring.fmats._fvars
    _Nk_ij = fusion_ring.Nk_ij
    one = fusion_ring.one()

    child_id, n_proc, fn_args = args
    k, a, b, n_strands = fn_args
    cdef int ctr = -1
    global worker_results
    #Get computational basis
    cdef list comp_basis = fusion_ring.get_computational_basis(a,b,n_strands)
    cdef dict basis_dict = { elt : i for i, elt in enumerate(comp_basis) }
    cdef int dim = len(comp_basis)
    cdef set coords = set()
    cdef int i
    #Avoid pickling cyclotomic field element objects
    must_flatten_coeff = fusion_ring.fvars_field() != QQbar
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
            if nnz_pos in comp_basis and (basis_dict[nnz_pos],i) not in coords:
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

                    #Avoid pickling cyclotomic field element objects
                    if must_flatten_coeff:
                        entry = entry.list()

                    worker_results.append(((basis_dict[nnz_pos],i), entry))
                    coords.add((basis_dict[nnz_pos],i))
                    continue

                entry = 0
                for p in fusion_ring.basis():
                    f1 = _fmat(_fvars,_Nk_ij,one,top_left,m[k-1],m[k],root,l[k-2],p)
                    f2 = _fmat(_fvars,_Nk_ij,one,top_left,f,e,root,q,p)
                    entry += f1 * mid_sig_ij(fusion_ring,(m[k-1],m[k]),(f,e),a,p) * f2

                #Avoid pickling cyclotomic field element objects
                if must_flatten_coeff:
                    entry = entry.list()

                worker_results.append(((basis_dict[nnz_pos],i), entry))

@cython.nonecheck(False)
@cython.cdivision(True)
cdef odd_one_out(fusion_ring, tuple args):
    """
    Compute entries of the rightmost braid generator, in case we have an odd number
    of strands
    """
    #Pre-compute common parameters for efficiency
    _fvars = fusion_ring.fmats._fvars
    _Nk_ij = fusion_ring.Nk_ij
    one = fusion_ring.one()

    global worker_results
    child_id, n_proc, fn_args = args
    a, b, n_strands = fn_args
    cdef ctr = -1
    #Get computational basis
    comp_basis = fusion_ring.get_computational_basis(a,b,n_strands)
    basis_dict = { elt : i for i, elt in enumerate(comp_basis) }
    dim = len(comp_basis)

    #Avoid pickling cyclotomic field element objects
    must_flatten_coeff = fusion_ring.fvars_field() != QQbar
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

                    #Avoid pickling cyclotomic field element objects
                    if must_flatten_coeff:
                        entry = entry.list()

                    worker_results.append(((basis_dict[nnz_pos],i), entry))
                    continue
                top_left = m[0]
                if n_strands > 5:
                    top_left = l[-2]
                root = b

                #Compute relevant entry
                entry = 0
                for p in fusion_ring.basis():
                    f1 = _fmat(_fvars,_Nk_ij,one,top_left,m[-1],a,root,l[-1],p)
                    f2 = _fmat(_fvars,_Nk_ij,one,top_left,f,a,root,q,p)
                    entry += f1 * odd_one_out_ij(fusion_ring,m[-1],f,a,p) * f2

                #Avoid pickling cyclotomic field element objects
                if must_flatten_coeff:
                    entry = entry.list()

                worker_results.append(((basis_dict[nnz_pos],i), entry))

################
### Reducers ###
################

cpdef collect_results(proc):
    """
    Helper function for returning processed results back to parent process.

    Trivial reducer: simply collects objects with the same key in the worker.
    This method is only useful when called after :meth:`executor`, whose
    function argument appends output to the ``worker_results`` list.

    EXAMPLES::

        sage: from sage.combinat.root_system.fast_parallel_fusion_ring_braid_repn import executor
        sage: FR = FusionRing("A1",4)                                                    # long time
        sage: FR.fusion_labels(['idd','one','two','three','four'],inject_variables=True) # long time
        sage: FR.fmats.find_orthogonal_solution(verbose=False)                           # long time
        sage: params = (('sig_2k',id(FR)),(0,1,(2,one,one,9)))                           # long time
        sage: executor(params)                                                           # long time
        sage: from sage.combinat.root_system.fast_parallel_fusion_ring_braid_repn import collect_results
        sage: len(collect_results(0)) == 171                                             # long time
        True
    """
    #Discard the zero polynomial
    global worker_results
    reduced = worker_results
    worker_results = list()
    return reduced

##############################
### Parallel code executor ###
##############################

#Hard-coded module __dict__-style attribute with visible cdef methods
cdef dict mappers = {
    "sig_2k": sig_2k,
    "odd_one_out": odd_one_out
}

cpdef executor(params):
    r"""
    Execute a function registered in this module's ``mappers``
    in a worker process, and supply the ``FusionRing`` parameter by
    constructing a reference to the FMatrix object in the worker's memory
    adress space from its ``id``.

    .. NOTE::

        When the parent process is forked, each worker gets a copy of
        every  global variable. The virtual memory address of object `X` in
        the parent process equals the *virtual* memory address of the copy of
        object `X` in each worker, so we may construct references to forked
        copies of `X`.

    TESTS::

        sage: from sage.combinat.root_system.fast_parallel_fusion_ring_braid_repn import executor
        sage: FR = FusionRing("A1",4)
        sage: FR.fusion_labels(['idd','one','two','three','four'],inject_variables=True)
        sage: FR.fmats.find_orthogonal_solution(verbose=False)    # long time
        sage: params = (('sig_2k',id(FR)),(0,1,(1,one,one,5)))    # long time
        sage: executor(params)                                    # long time
        sage: from sage.combinat.root_system.fast_parallel_fusion_ring_braid_repn import collect_results
        sage: len(collect_results(0)) == 13                       # long time
        True
        sage: from sage.combinat.root_system.fast_parallel_fusion_ring_braid_repn import executor, collect_results
        sage: FR = FusionRing("B2",2)                             # long time
        sage: FR.fusion_labels(['I0','Y1','X','Z','Xp','Y2'],inject_variables=True)         # long time
        sage: FR.fmats.find_orthogonal_solution(verbose=False)    # long time (~23 s)
        sage: params = (('odd_one_out',id(FR)),(0,1,(X,Xp,5)))    # long time
        sage: executor(params)                                    # long time
        sage: len(collect_results(0)) == 54                       # long time
        True
    """
    (fn_name, fr_id), args = params
    #Construct a reference to global FMatrix object in this worker's memory
    fusion_ring_obj = ctypes.cast(fr_id, ctypes.py_object).value
    #Bind module method to FMatrix object in worker process, and call the method
    mappers[fn_name](fusion_ring_obj,args)

######################################
### Pickling circumvention helpers ###
######################################

cpdef _unflatten_entries(fusion_ring, list entries):
    """
    Restore cyclotomic coefficient object from its tuple of rational
    coefficients representation.

    Used to circumvent pickling issue introduced by PARI settigs in trac
    ticket #30537

    EXAMPLES::

        sage: from sage.combinat.root_system.fast_parallel_fusion_ring_braid_repn import _unflatten_entries
        sage: fr = FusionRing("B2",2)
        sage: F = fr.field()
        sage: coeff = [F.random_element() for i in range(2)]
        sage: entries = [((0,0), coeff[0].list()), ((0,1), coeff[1].list())]
        sage: _unflatten_entries(fr, entries)
        sage: all(cyc_elt_obj == c for (coord, cyc_elt_obj), c in zip(entries, coeff))
        True
    """
    F = fusion_ring.fvars_field()
    fm = fusion_ring.fmats
    must_unflatten = F != QQbar
    if must_unflatten:
        for i, (coord, entry) in enumerate(entries):
            entries[i] = (coord, F(entry))
