
include '../ext/cdefs.pxi'
include '../ext/python_mem.pxi'
include '../ext/stdsage.pxi'
from sage.structure.element import is_Matrix
#from sage.graphs.graph_isom cimport OrbitPartition,\
#    _orbit_partition_from_list_perm, PartitionStack
from sage.misc.misc import cputime
from sage.rings.integer import Integer

cdef class BinaryCodeGraph:

    def __new__(self, arg1, arg2=None, arg3=None, arg4=None):
        cdef unsigned int i, j, k, z

        if arg2 is None and is_Matrix(arg1):

            if arg1.nrows() > 8*sizeof(int):
                raise NotImplementedError("Columns are stored as ints. This code is too big.")
            self.radix = 8*sizeof(int)
            self.ncols = arg1.ncols()
            self.nrows = arg1.nrows()
            self.columns = <int *> sage_malloc( arg1.ncols() * sizeof(int) )
            if not self.columns:
                raise MemoryError("Memory.")

            cols = arg1.columns()
            for i from 0 <= i < self.ncols:
                k = 0
                for j in cols[i].nonzero_positions():
                    k += (1 << j)
                self.columns[i] = k

            self.ham_wts = <int *> sage_malloc( 16 * sizeof(int) )
            if not self.ham_wts:
                sage_free(self.columns)
                raise MemoryError("Memory.")
            self.ham_wts[0]  = 0; self.ham_wts[1]  = 1; self.ham_wts[2]  = 1; self.ham_wts[3]  = 2
            self.ham_wts[4]  = 1; self.ham_wts[5]  = 2; self.ham_wts[6]  = 2; self.ham_wts[7]  = 3
            self.ham_wts[8]  = 1; self.ham_wts[9]  = 2; self.ham_wts[10] = 2; self.ham_wts[11] = 3
            self.ham_wts[12] = 2; self.ham_wts[13] = 3; self.ham_wts[14] = 3; self.ham_wts[15] = 4

            z = 1<<32
            self.test = <int *> sage_malloc( (z) * sizeof(int) )
            if not self.test:
                sage_free(self.columns)
                sage_free(self.ham_wts)
                raise MemoryError("Memory.")

            t = cputime()
            for i from 0 <= i < z:
                self.test[i] = 1
            t = cputime(t)
            print t

    def __dealloc__(self):
        sage_free(self.columns)
        sage_free(self.ham_wts)
        sage_free(test)

    cdef int has_edge(self, int word, int column):
        cdef int i, j, k, l
        i = 0
        k = word & self.columns[column]
        for j from 0 <= j < self.radix by 4:
            k = k >> 4
            i += self.ham_wts[15 & k]
        return i % 2

def search_tree( BinaryCodeGraph CG, lab=True, verbosity=0 ):

    cdef int i, j, # local variables
    cdef OrbitPartition Theta, OP
    cdef int index = 0, size = 1
    cdef int L = 100 # memory limit for storing values from fix and mcr:
                     # Phi and Omega store specific information about some
                     # of the automorphisms we already know about, and they
                     # are arrays of length L
    cdef int **Phi # stores the fixed point sets of each automorphism
    cdef int **Omega # stores the minimal elements of each cell of the
                     # orbit partition
    cdef int l = -1 # current index for storing values in Phi and Omega-
                    # we start at -1 so that when we increment first,
                    # the first place we write to is 0.
    cdef PartitionStack nu, zeta, rho
    cdef int k_rho # the number of partitions in rho
    cdef int k = 0 # the number of partitions in nu
    cdef int h = -1 # longest common ancestor of zeta and nu:
                    # zeta[h] == nu[h], zeta[h+1] != nu[h+1]
    cdef int hb     # longest common ancestor of rho and nu:
                    # rho[hb] == nu[hb], rho[hb+1] != nu[hb+1]
    cdef int hh = 1 # the height of the oldest ancestor of nu
                    # satisfying Lemma 2.25 in [1]
    cdef int ht # smallest such that all descendants of zeta[ht] are equivalent

    cdef mpz_t *Lambda_mpz, *zf_mpz, *zb_mpz # for tracking indicator values
    # zf and zb are indicator vectors remembering Lambda[k] for zeta and rho,
    # respectively
    cdef int hzf      # the max height for which Lambda and zf agree
    cdef int hzb = -1 # the max height for which Lambda and zb agree

    cdef int *gamma # for storing permutations
    cdef int *alpha # for storing pointers to cells of nu[k]:
                     # allocated to be length 4*n + 1 for scratch (see functions
                     # _sort_by_function and _refine_by_square_matrix)
    cdef int *v # list of vertices determining nu
    cdef int *e # 0 or 1, see states 12 and 17
    cdef int state # keeps track of place in algorithm
    cdef int tvc, tvh, n1 = 2**CG.nrows, n2 = CG.ncols, n = n1 + n2
    Pi = [range(n1), range(n1,n1+n2)]

    # trivial case
    if n == 0:
        return [], {}

    # allocate int pointers
    Phi = <int **> sage_malloc( L * sizeof(int *) )
    Omega = <int **> sage_malloc( L * sizeof(int *) )

    # allocate GMP int pointers
    Lambda_mpz = <mpz_t *> sage_malloc( (n+2) * sizeof(mpz_t) )
    zf_mpz = <mpz_t *> sage_malloc( (n+2) * sizeof(mpz_t) )
    zb_mpz = <mpz_t *> sage_malloc( (n+2) * sizeof(mpz_t) )

    # check for memory errors
    if not (Phi and Omega and Lambda_mpz and zf_mpz and zb_mpz):
        if Lambda_mpz: sage_free(Lambda_mpz)
        if zf_mpz: sage_free(zf_mpz)
        if zb_mpz: sage_free(zb_mpz)
        if Phi: sage_free(Phi)
        if Omega: sage_free(Omega)
        raise MemoryError("Error allocating memory.")

    # allocate int arrays
    gamma = <int *> sage_malloc( n * sizeof(int) )
    Phi[0] = <int *> sage_malloc( (L*n) * sizeof(int) )
    Omega[0] = <int *> sage_malloc( (L*n) * sizeof(int) )
    alpha = <int *> sage_malloc( (4*n + 1) * sizeof(int) )
    v = <int *> sage_malloc( n * sizeof(int) )
    e = <int *> sage_malloc( n * sizeof(int) )

    # check for memory errors
    if not (gamma and Phi[0] and Omega[0] and alpha and v and e):
        if gamma: sage_free(gamma)
        if Phi[0]: sage_free(Phi[0])
        if Omega[0]: sage_free(Omega[0])
        if alpha: sage_free(alpha)
        if v: sage_free(v)
        if e: sage_free(e)
        sage_free(Lambda_mpz)
        sage_free(zf_mpz)
        sage_free(zb_mpz)
        sage_free(Phi)
        sage_free(Omega)
        raise MemoryError("Error allocating memory.")

    # setup double index arrays
    for i from 0 < i < L:
        Phi[i] = Phi[0] + n*i
    for i from 0 < i < L:
        Omega[i] = Omega[0] + n*i

    # allocate GMP ints
    for i from 0 <= i < n+2:
        mpz_init(Lambda_mpz[i])
        mpz_init_set_si(zf_mpz[i], -1) # correspond to default values of
        mpz_init_set_si(zb_mpz[i], -1) # "infinity"
        # Note that there is a potential memory leak here - if a particular
        # mpz fails to allocate, this is not checked for

    # set up the rest of the variables
    nu = PartitionStack(Pi)
    Theta = OrbitPartition(n)
    output = []

##############################################################################
####### W, M, Pi, dig, dict, certify #########################################
##############################################################################

    # TODO: Implement W

    state = 1
    while state != -1:
        if verbosity > 0:
            print state

        if state == 1: # Entry point to algorithm
            # get alpha to point to cells of nu
            j = 1
            alpha[0] = 0
            for i from 0 < i < n:
                if nu.levels[i-1] == 0:
                    alpha[j] = i
                    j += 1
            alpha[j] = -1

            # "nu[0] := R(G, Pi, Pi)"
            nu._refine_by_square_matrix(k, alpha, n, M, _dig)

            if not _dig:
                if nu._sat_225(k, n): hh = k
            if nu._is_discrete(k): state = 18; continue

            # store the first smallest nontrivial cell in W[k], and set v[k]
            # equal to its minimum element
            v[k] = nu._first_smallest_nontrivial(W[k], k, n)
            mpz_set_ui(Lambda_mpz[k], 0)
            e[k] = 0 # see state 12, and 17
            state = 2

        elif state == 2: # Move down the search tree one level by refining nu
            k += 1

            # "nu[k] := nu[k-1] perp v[k-1]"
            nu._clear(k)
            alpha[0] = nu._split_vertex(v[k-1], k)
            alpha[1] = -1
            i = nu._refine_by_square_matrix(k, alpha, n, M, _dig)

            # add one, then multiply by the invariant
            mpz_add_ui(Lambda_mpz[k], Lambda_mpz[k-1], 1)
            mpz_mul_si(Lambda_mpz[k], Lambda_mpz[k], i)

            # only if this is the first time moving down the search tree:
            if h == -1: state = 5; continue

            # update hzf
            if hzf == k-1 and mpz_cmp(Lambda_mpz[k], zf_mpz[k]) == 0: hzf = k
            if not lab: state = 3; continue

            # "qzb := cmp(Lambda[k], zb[k])"
            if mpz_cmp_si(zb_mpz[k], -1) == 0: # if "zb[k] == oo"
                qzb = -1
            else:
                qzb = mpz_cmp( Lambda_mpz[k], zb_mpz[k] )
            # update hzb
            if hzb == k-1 and qzb == 0: hzb = k

            # if Lambda[k] > zb[k], then zb[k] := Lambda[k]
            # (zb keeps track of the indicator invariants corresponding to
            # rho, the closest canonical leaf so far seen- if Lambda is
            # bigger, then rho must be about to change
            if qzb > 0: mpz_set(zb_mpz[k], Lambda_mpz[k])
            state = 3

        elif state == 3: # attempt to rule out automorphisms while moving down
                         # the tree
            if hzf <= k or (lab and qzb >= 0): # changed hzb to hzf, == to <=
                state = 4
            else: state = 6
            # if k > hzf, then we know that nu currently does not look like
            # zeta, the first terminal node encountered. Then if we are not
            # looking for a canonical label, there is no reason to continue.
            # However, if we are looking for one, and qzb < 0, i.e.
            # Lambda[k] < zb[k], then the indicator is not maximal, and we
            # can't reach a canonical leaf.

        elif state == 4: # at this point we have -not- ruled out the presence
                         # of automorphisms
            if nu._is_discrete(k): state = 7; continue

            # store the first smallest nontrivial cell in W[k], and set v[k]
            # equal to its minimum element
            v[k] = nu._first_smallest_nontrivial(W[k], k, n)

            if _dig or not nu._sat_225(k, n): hh = k + 1
            e[k] = 0 # see state 12, and 17
            state = 2 # continue down the tree

        elif state == 5: # alternative to 3: since we have not yet gotten
                         # zeta, there are no automorphisms to rule out.
                         # instead we record Lambda to zf and zb
                         # (see state 3)
            mpz_set(zf_mpz[k], Lambda_mpz[k])
            mpz_set(zb_mpz[k], Lambda_mpz[k])
            state = 4

        elif state == 6: # at this stage, there is no reason to continue
                         # downward, and an automorphism has not been
                         # discovered
            j = k

            # return to the longest ancestor nu[i] of nu that could have a
            # descendant equivalent to zeta or could improve on rho.
            # All terminal nodes descending from nu[hh] are known to be
            # equivalent, so i < hh. Also, if i > hzb, none of the
            # descendants of nu[i] can improve rho, since the indicator is
            # off (Lambda(nu) < Lambda(rho)). If i >= ht, then no descendant
            # of nu[i] is equivalent to zeta (see [1, p67]).
            if ht-1 > hzb:
                if ht-1 < hh-1:
                    k = ht-1
                else:
                    k = hh-1
            else:
                if hzb < hh-1:
                    k = hzb
                else:
                    k = hh-1

            # TODO: investigate the following line
            if k == -1: k = 0 # not in BDM, broke at G = Graph({0:[], 1:[]}), Pi = [[0,1]], lab=False

            if j == hh: state = 13; continue
            # recall hh: the height of the oldest ancestor of zeta for which
            # Lemma 2.25 is satsified, which implies that all terminal nodes
            # descended from there are equivalent (or simply k if 2.25 does
            # not apply). If we are looking at such a node, then the partition
            # at nu[hh] can be used for later pruning, so we store its fixed
            # set and a set of representatives of its cells
            if l < L-1: l += 1
            for i from 0 <= i < n:
                Omega[l][i] = 0 # changed Lambda to Omega
                Phi[l][i] = 0
                if nu._is_min_cell_rep(i, hh):
                    Omega[l][i] = 1
                    if nu._is_fixed(i, hh):
                        Phi[l][i] = 1

            state = 12

        elif state == 7: # we have just arrived at a terminal node of the
                         # search tree T(G, Pi)
            # if this is the first terminal node, go directly to 18, to
            # process zeta
            if h == -1: state = 18; continue

            # hzf is the extremal height of ancestors of both nu and zeta,
            # so if k < hzf, nu is not equivalent to zeta, i.e. there is no
            # automorphism to discover.
            # TODO: investigate why, in practice, the same does not seem to be
            # true for hzf < k... BDM had !=, not <, and this broke at
            # G = Graph({0:[],1:[],2:[]}), Pi = [[0,1,2]]
            if k < hzf: state = 8; continue

            # get the permutation corresponding to this terminal node
            nu._get_permutation_from(zeta, gamma)

            if verbosity > 3:
                print 'automorphism discovered:'
                print [gamma[iii] for iii in range(n)]

            # if G^gamma == G, the permutation is an automorphism, goto 10
            if G_enum == _enumerate_graph_with_permutation(M, n, gamma):
                state = 10
            else:
                state = 8

        elif state == 8: # we have just ruled out the presence of automorphism
                         # and have not yet considered whether nu improves on
                         # rho
            # if we are not searching for a canonical label, there is nothing
            # to do here
            if (not lab) or (qzb < 0):
                state = 6; continue

            # if Lambda[k] > zb[k] or nu is shorter than rho, then we have
            # found an improvement for rho
            if (qzb > 0) or (k < k_rho):
                state = 9; continue

            # if G(nu) > G(rho), goto 9
            # if G(nu) < G(rho), goto 6
            # if G(nu) == G(rho), get the automorphism and goto 10
            m1 = nu._enumerate_graph_from_discrete(M, n)
            m2 = rho._enumerate_graph_from_discrete(M, n)

            if m1 > m2:
                state = 9; continue
            if m1 < m2:
                state = 6; continue

            rho._get_permutation_from(nu, gamma)
            if verbosity > 3:
                print 'automorphism discovered:'
                print [gamma[iii] for iii in range(n)]
            state = 10

        elif state == 9: # entering this state, nu is a best-so-far guess at
                         # the canonical label
            rho = PartitionStack(nu)
            k_rho = k

            qzb = 0
            hb = k
            hzb = k

            # set zb[k+1] = Infinity
            mpz_set_si(zb_mpz[k+1], -1)
            state = 6

        elif state == 10: # we have an automorphism to process
            # increment l
            if l < L - 1:
                l += 1

            # retrieve the orbit partition, and record the relevant
            # information
            # TODO: this step could be optimized. The variable OP is not
            # really necessary
            OP = orbit_partition_from_list_perm(gamma, n)
            for i from 0 <= i < n:
                Omega[l][i] = OP._is_min_cell_rep(i)
                Phi[l][i] = OP._is_fixed(i)

            # if each orbit of gamma is part of an orbit in Theta, then the
            # automorphism is already in the span of those we have seen
            if OP._is_finer_than(Theta, n):
                state = 11
                continue
            # otherwise, incorporate this into Theta
            Theta._vee_with(OP, n)

            # record the automorphism
            output.append([ Integer(gamma[i]) for i from 0 <= i < n ])

            # The variable tvc was set to be the minimum element of W[k]
            # the last time we were at state 13 and at a node descending to
            # zeta. If this is a minimal cell representative of Theta and
            # we are searching for a canonical label, goto state 11, i.e.
            # backtrack to the common ancestor of rho and nu, then goto state
            # 12, i.e. consider whether we still need to search downward from
            # there. TODO: explain why
            if Theta.elements[tvc] == -1 and lab: ## added "and lab"
                state = 11
                continue
            k = h
            state = 13

        elif state == 11: # if we are searching for a label, backtrack to the
                          # common ancestor of nu and rho
            k = hb
            state = 12

        elif state == 12: # we are looking at a branch we may have to continue
                          # to search downward on
            # e keeps track of the motion through the search tree. It is set to
            # 1 when you have just finished coming up the search tree, and are
            # at a node in the tree for which there may be more branches left
            # to explore. In this case, intersect W[k] with Omega[l], since
            # there may be an automorphism mapping one element of W[k] to
            # another, hence only one must be investigated.
            if e[k] == 1:
                for j from 0 <= j < n:
                    if W[k][j] and not Omega[l][j]:
                        W[k][j] = 0
            state = 13

        elif state == 13: # hub state
            if k == -1:
                # the algorithm has finished
                state = -1; continue
            if k > h:
                # if we are not at a node of zeta
                state = 17; continue
            if k == h:
                # if we are at a node of zeta, then state 14 can rule out
                # vertices to consider
                state = 14; continue

            # thus, it must be that k < h, and this means we are done
            # searching underneath zeta[k+1], so now, k is the new longest
            # ancestor of nu and zeta:
            h = k

            # set tvc and tvh to the minimum cell representative of W[k]
            # (see states 10 and 14)
            for i from 0 <= i < n:
                if W[k][i]:
                    tvc = i
                    break
            tvh = tvc
            state = 14

        elif state == 14: # iterate v[k] through W[k] until a minimum cell rep
                          # of Theta is found
            # The variable tvh was set to be the minimum element of W[k]
            # the last time we were at state 13 and at a node descending to
            # zeta. If this is in the same cell of Theta as v[k], increment
            # index (see Theorem 2.33 in [1])
            if Theta._find(v[k]) == Theta._find(tvh):
                index += 1

            # find the next v[k] in W[k]
            i = v[k] + 1
            while i < n and not W[k][i]:
                i += 1
            if i < n:
                v[k] = i
            else:
                # there is no new vertex to consider at this level
                v[k] = -1
                state = 16
                continue

            # if the new v[k] is not a minimum cell representative of Theta,
            # then we already considered that rep., and that subtree was
            # isomorphic to the one corresponding to v[k]
            if Theta.elements[v[k]] != -1: state = 14
            else:
                # otherwise, we do have a vertex to consider
                state = 15

        elif state == 15: # we have a new vertex, v[k], that we must split on
            # hh is smallest such that nu[hh] satisfies Lemma 2.25. If it is
            # larger than k+1, it must be modified, since we are changing that
            # part
            if k + 1 < hh:
                hh = k + 1
            # hzf is maximal such that indicators line up for nu and zeta
            if k < hzf:
                hzf = k
            if not lab or hb < k: # changed hzb to hb
                # in either case there is no need to update hb, which is the
                # length of the common ancestor of nu and rho
                state = 2; continue
            hb = k # changed hzb to hb
            qzb = 0
            state = 2

        elif state == 16: # backtrack one level in the search tree, recording
                          # information relevant to Theorem 2.33
            j = 0
            for i from 0 <= i < n:
                if W[k][i]: j += 1
            if j == index and ht == k+1: ht = k
            size = size*index
            index = 0
            k -= 1
            state = 13

        elif state == 17: # you have just finished coming up the search tree,
                          # and must now consider going back down.
            if e[k] == 0:
                # intersect W[k] with each Omega[i] such that {v_0,...,v_(k-1)}
                # is contained in Phi[i]
                for i from 0 <= i <= l:
                    # check if {v_0,...,v_(k-1)} is contained in Phi[i]
                    # i.e. fixed pointwise by the automorphisms so far seen
                    j = 0
                    while j < k and Phi[i][v[j]]:
                        j += 1
                    # if so, only check the minimal orbit representatives
                    if j == k:
                        for j from 0 <= j < n:
                            if W[k][j] and not Omega[i][j]:
                                W[k][j] = 0
            e[k] = 1 # see state 12

            # see if there is a relevant vertex to split on:
            i = v[k] + 1
            while i < n and not W[k][i]:
                i += 1
            if i < n:
                v[k] = i
                state = 15
                continue
            else:
                v[k] = -1

            # otherwise backtrack one level
            k -= 1
            state = 13

        elif state == 18: # The first time we encounter a terminal node, we
                          # come straight here to set up zeta. This is a one-
                          # time state.
            # initialize counters for zeta:
            h = k # zeta[h] == nu[h]
            ht = k # nodes descended from zeta[ht] are all equivalent
            hzf = k # max such that indicators for zeta and nu agree

            zeta = PartitionStack(nu)

            k -= 1
            if not lab: state = 13; continue

            rho = PartitionStack(nu)

            # initialize counters for rho:
            k_rho = k # number of partitions in rho
            hzb = k # max such that indicators for rho and nu agree - BDM had k+1
            hb = k # rho[hb] == nu[hb] - BDM had k+1

            qzb = 0 # Lambda[k] == zb[k], so...
            state = 13

    # free the GMP ints
    for i from 0 <= i < n+2:
        mpz_clear(Lambda_mpz[i])
        mpz_clear(zf_mpz[i])
        mpz_clear(zb_mpz[i])

    # free int arrays
    sage_free(gamma)
    sage_free(Phi[0])
    sage_free(Omega[0])
    sage_free(alpha)
    sage_free(v)
    sage_free(e)

    # free GMP int pointers
    sage_free(Lambda_mpz)
    sage_free(zf_mpz)
    sage_free(zb_mpz)

    # free int pointers
    sage_free(Phi)
    sage_free(Omega)

    return output, dd



















