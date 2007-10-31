
include '../ext/cdefs.pxi'
include '../ext/python_mem.pxi'
include '../ext/stdsage.pxi'
from sage.structure.element import is_Matrix
#from sage.graphs.graph_isom cimport OrbitPartition,\
#    _orbit_partition_from_list_perm
from sage.misc.misc import cputime
from math import log, ceil
from sage.rings.integer import Integer

cdef class BinaryCodeGraph:

    def __new__(self, arg1):
        cdef unsigned int i, k, j, z
        self.radix = 8*sizeof(int)
        self.ncols = arg1.ncols()
        self.nrows = arg1.nrows()
        self.nwords = 1 << self.nrows
        self.ptn_mask = ~0 << self.nrows

        if self.nrows >= self.radix or ceil(log(self.ncols,2)) >= self.radix:
            raise NotImplementedError("Columns and rows are stored as ints. This code is too big.")
        self.columns = <int *> sage_malloc( self.ncols * sizeof(int) )
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

    def __dealloc__(self):
        sage_free(self.columns)
        sage_free(self.ham_wts)

    cdef int has_edge_bip(self, int word, int column):
        cdef int i, j, k
        i = 0
        k = word & self.columns[column]
        for j from 0 <= j < self.radix by 4:
            i += self.ham_wts[15 & k]
            k = k >> 4
        return i % 2

    cdef int has_edge(self, int a, int b):
        cdef int i, j
        i = self.ptn_mask & a
        j = self.ptn_mask & b
        if (i==0)==(0==j):
            # both are on the same side of the bipartite graph
            return 0
        elif i:
            # first argument is a column
            return self.has_edge_bip(b, a - self.nwords)
        else:
            # first argument is a vector
            return self.has_edge_bip(a, b - self.nwords)

#    cdef int is_automorphism(self, int *wd_gamma, int *col_gamma):

        #????

        #cdef int i, j, k, l
        #for i from 0 <= i < self.nwords:
        #    for j from 0 <= j < self.ncols:
        #        k = CG.has_edge_bip(self.wd_ents[i], self.col_ents[j])
        #        l = CG.has_edge_bip(other.wd_ents[i], other.col_ents[j])
        #        if k ^ l:
        #            return k - l
        #return 0

cdef class PartitionStack:

    def __new__(self, arg1, arg2=None):
        cdef int k
        self.nwords = int(arg1)
        self.ncols = int(arg2)
        self.wd_ents = <int *> sage_malloc( self.nwords * sizeof(int) )
        self.wd_lvls = <int *> sage_malloc( self.nwords * sizeof(int) )
        self.col_ents = <int *> sage_malloc( self.ncols * sizeof(int) )
        self.col_lvls = <int *> sage_malloc( self.ncols * sizeof(int) )
        if not (self.wd_ents and self.wd_lvls and self.col_ents and self.col_lvls):
            if self.wd_ents: sage_free(self.wd_ents)
            if self.wd_lvls: sage_free(self.wd_lvls)
            if self.col_ents: sage_free(self.col_ents)
            if self.col_lvls: sage_free(self.col_lvls)
        for k from 0 <= k < self.nwords-1:
            self.wd_ents[k] = k
            self.wd_lvls[k] = self.nwords
        for k from 0 <= k < self.ncols-1:
            self.col_ents[k] = k
            self.col_lvls[k] = self.ncols
        self.wd_ents[self.nwords-1] = self.nwords-1
        self.wd_lvls[self.nwords-1] = -1
        self.col_ents[self.ncols-1] = self.ncols-1
        self.col_lvls[self.ncols-1] = -1

    def __dealloc__(self):
        sage_free(self.wd_ents)
        sage_free(self.wd_lvls)
        sage_free(self.col_ents)
        sage_free(self.col_lvls)

    def __repr__(self):
        cdef int i, j, k
        s = ''
        for i from 0 <= i < self.nwords:
            s += '({'
            for j from 0 <= j < self.nwords:
                s += str(self.wd_ents[j])
                if self.wd_lvls[j] <= i:
                    s += '},{'
                else:
                    s += ','
            s = s[:-2] + ')\n'
        for i from 0 <= i < self.ncols:
            s += '({'
            for j from 0 <= j < self.ncols:
                s += str(self.col_ents[j])
                if self.col_lvls[j] <= i:
                    s += '},{'
                else:
                    s += ','
            s = s[:-2] + ')\n'
        return s

    cdef int is_discrete(self, int k):
        cdef int i
        for i from 0 <= i < self.ncols:
            if self.col_lvls[i] > k:
                return 0
        return 1

    cdef int num_cells(self, int k):
        cdef int i, j = 0
        for i from 0 <= i < self.nwords:
            if self.wd_lvls[i] <= k:
                j += 1
        for i from 0 <= i < self.ncols:
            if self.col_lvls[i] <= k:
                j += 1
        return j

    cdef int sat_225(self, int k):
        cdef int i, n = self.nwords + self.ncols, in_cell = 0
        cdef int nontrivial_cells = 0, total_cells = self.num_cells(k)
        if n <= total_cells + 4:
            return 1
        for i from 0 <= i < self.nwords:
            if self.wd_lvls[i] <= k:
                if in_cell:
                    nontrivial_cells += 1
                in_cell = 0
            else:
                in_cell = 1
        in_cell = 0
        for i from 0 <= i < self.ncols:
            if self.col_lvls[i] <= k:
                if in_cell:
                    nontrivial_cells += 1
                in_cell = 0
            else:
                in_cell = 1
        if n == total_cells + nontrivial_cells:
            return 1
        if n == total_cells + nontrivial_cells + 1:
            return 1
        return 0

    cdef int is_min_cell_rep(self, int col_not_wd, int i, int k):
        if i == 0:
            return 1
        if col_not_wd:
            return self.col_lvls[i-1] <= k
        else:
            return self.wd_lvls[i-1] <= k

    cdef int is_fixed(self, int col_not_wd, int i, int k):
        """
        Assuming it is a min cell rep.
        """
        if col_not_wd:
            return self.col_lvls[i] <= k
        else:
            return self.wd_lvls[i] <= k

    ## TODO: first smallest nontrivial

    cdef void col_percolate(self, int start, int end):
        cdef int i, temp
        for i from end >= i > start:
            if self.col_ents[i] < self.col_ents[i-1]:
                temp = self.col_ents[i]
                self.col_ents[i] = self.col_ents[i-1]
                self.col_ents[i-1] = temp

    cdef void wd_percolate(self, int start, int end):
        cdef int i, temp
        for i from end >= i > start:
            if self.wd_ents[i] < self.wd_ents[i-1]:
                temp = self.wd_ents[i]
                self.wd_ents[i] = self.wd_ents[i-1]
                self.wd_ents[i-1] = temp

    cdef int split_vertex(self, int col_not_wd, int v, int k):
        cdef int i = 0, j
        if col_not_wd:
            while self.col_ents[i] != v: i += 1
            j = i
            while self.col_lvls[i] > k: i += 1
            if j == 0 or self.col_lvls[j-1] <= k:
                self.col_percolate(j+1, i)
            else:
                while j != 0 and self.col_lvls[j-1] > k:
                    self.col_ents[j] = self.col_ents[j-1]
                    j -= 1
                self.col_ents[j] = v
            self.col_lvls[j] = k
            return j
        else:
            while self.wd_ents[i] != v: i += 1
            j = i
            while self.wd_lvls[i] > k: i += 1
            if j == 0 or self.wd_lvls[j-1] <= k:
                self.wd_percolate(j+1, i)
            else:
                while j != 0 and self.wd_lvls[j-1] > k:
                    self.wd_ents[j] = self.wd_ents[j-1]
                    j -= 1
                self.wd_ents[j] = v
            self.wd_lvls[j] = k
            return j

    cdef int col_degree(self, BinaryCodeGraph CG, int col, int wd_ptr, int k):
        cdef int i = 0
        col = self.col_ents[col]
        while True:
            if CG.has_edge_bip(wd_ptr, col): i += 1
            if self.wd_lvls[wd_ptr] > k: wd_ptr += 1
            else: break
        return i

    cdef int wd_degree(self, BinaryCodeGraph CG, int wd, int col_ptr, int k):
        cdef int i = 0
        wd = self.wd_ents[wd]
        while True:
            if CG.has_edge_bip(wd, col_ptr): i += 1
            if self.col_lvls[col_ptr] > k: col_ptr += 1
            else: break
        return i

    cdef int sort_cols(self, int start, int *degrees, int k):
        cdef int i, j, max, max_location
        cdef int *counts = degrees + self.ncols, *output = degrees + 2*self.ncols

        for i from 0 <= i < self.ncols:
            counts[i] = 0
        i = 0
        while self.col_lvls[i+start] > k:
            counts[degrees[i]] += 1
            i += 1
        counts[degrees[i]] += 1

        # i+start is the right endpoint of the cell now
        max = counts[0]
        max_location = 0
        for j from 0 < j < self.ncols:
            if counts[j] > max:
                max = counts[j]
                max_location = j
            counts[j] += counts[j-1]

        for j from i >= j >= 0:
            counts[degrees[j]] -= 1
            output[counts[degrees[j]]] = self.col_ents[start+j]

        max_location = counts[max_location] + start

        for j from 0 <= j <= i:
            self.col_ents[start+j] = output[j]

        j = 1
        while j < self.ncols and counts[j] <= i:
            if counts[j] > 0:
                self.col_lvls[start + counts[j] - 1] = k
            self.col_percolate(start + counts[j-1], start + counts[j] - 1)
            j += 1

        return max_location

    cdef int sort_wds(self, int start, int *degrees, int k):
        cdef int i, j, max, max_location
        cdef int *counts = degrees + self.nwords, *output = degrees + 2*self.nwords

        for i from 0 <= i < self.nwords:
            counts[i] = 0
        i = 0
        while self.wd_lvls[i+start] > k:
            counts[degrees[i]] += 1
            i += 1
        counts[degrees[i]] += 1

        # i+start is the right endpoint of the cell now
        max = counts[0]
        max_location = 0
        for j from 0 < j < self.nwords:
            if counts[j] > max:
                max = counts[j]
                max_location = j
            counts[j] += counts[j-1]

        for j from i >= j >= 0:
            counts[degrees[j]] -= 1
            output[counts[degrees[j]]] = self.wd_ents[start+j]

        max_location = counts[max_location] + start

        for j from 0 <= j <= i:
            self.wd_ents[start+j] = output[j]

        j = 1
        while j < self.nwords and counts[j] <= i:
            if counts[j] > 0:
                self.wd_lvls[start + counts[j] - 1] = k
            self.wd_percolate(start + counts[j-1], start + counts[j] - 1)
            j += 1

        return max_location

    cdef int refine(self, int k, int *col_alpha, int *wd_alpha, BinaryCodeGraph CG):
        cdef int m = 0, j
        cdef int i, q, r, s, t
        cdef int invariant
        cdef int *col_degrees = col_alpha + self.ncols
        cdef int *wd_degrees = wd_alpha + self.nwords
        while not self.is_discrete(k) and col_alpha[m] != -1:
            invariant += 1
            j = 0
            while j < self.nwords:
                invariant += 8
                i = j; s = 0
                while True:
                    wd_degrees[i-j] = self.wd_degree(CG, i, col_alpha[m], k)
                    if s == 0 and wd_degrees[i-j] != wd_degrees[0]: s = 1
                    i += 1
                    if self.wd_lvls[i-1] <= k: break
                if s:
                    invariant += 8
                    t = self.sort_wds(j, wd_degrees, k)
                    invariant += t + wd_degrees[i-j-1]
                    q = m
                    while col_alpha[q] != -1:
                        if col_alpha[q] == j: col_alpha[q] = t
                        q += 1
                    r = j
                    while True:
                        if r == 0 or self.wd_lvls[r-1] == k:
                            if r != t:
                                col_alpha[q] = r
                                q += 1
                        r += 1
                        if r >= i: break
                    col_alpha[q] = -1
                    while self.wd_lvls[j] > k:
                        j += 1
                    j += 1
                    invariant += (i - j)
                else: j = i
            m += 1
        m = 0
        while not self.is_discrete(k) and wd_alpha[m] != -1:
            invariant += 1
            j = 0
            while j < self.ncols:
                invariant += 8
                i = j; s = 0
                while True:
                    col_degrees[i-j] = self.col_degree(CG, i, wd_alpha[m], k)
                    if s == 0 and col_degrees[i-j] != col_degrees[0]: s = 1
                    i += 1
                    if self.col_lvls[i-1] <= k: break
                if s:
                    invariant += 8
                    t = self.sort_cols(j, col_degrees, k)
                    invariant += t + col_degrees[i-j-1]
                    q = m
                    while wd_alpha[q] != -1:
                        if wd_alpha[q] == j: wd_alpha[q] = t
                        q += 1
                    r = j
                    while True:
                        if r == 0 or self.col_lvls[r-1] == k:
                            if r != t:
                                wd_alpha[q] = r
                                q += 1
                        r += 1
                        if r >= i: break
                    wd_alpha[q] = -1
                    while self.col_lvls[j] > k:
                        j += 1
                    j += 1
                    invariant += (i - j)
                else: j = i
            m += 1
        return invariant

    cdef void get_permutation(self, PartitionStack zeta, int *wd_gamma, int *col_gamma):
        cdef int i
        for i from 0 <= i < self.ncols:
            col_gamma[zeta.col_ents[i]] = self.col_ents[i]
        for i from 0 <= i < self.nwords:
            wd_gamma[zeta.wd_ents[i]] = self.wd_ents[i]

    cdef int cmp(self, PartitionStack other, BinaryCodeGraph CG):
        # if CG(self) > G(other): return 1
        # if CG(self) < G(other): return -1
        # else: return 0
        cdef int i, j, k, l
        for i from 0 <= i < CG.nwords:
            for j from 0 <= j < CG.ncols:
                k = CG.has_edge_bip(self.wd_ents[i], self.col_ents[j])
                l = CG.has_edge_bip(other.wd_ents[i], other.col_ents[j])
                if k ^ l:
                    return k - l
        return 0

