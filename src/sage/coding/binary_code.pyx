"""
Fast binary code routines.

Some computations with linear binary codes. Fix a basis for $GF(2)^n$.
A linear binary code is a linear subspace of $GF(2)^n$, together with
this choice of basis. A permutation $g \in S_n$ of the fixed basis
gives rise to a permutation of the vectors, or words, in $GF(2)^n$,
sending $(w_i)$ to $(w_{g(i)})$. The permutation automorphism group of
the code $C$ is the set of permutations of the basis that bijectively
map $C$ to itself. Note that if $g$ is such a permutation, then
$$g(a_i) + g(b_i) = (a_{g(i)} + b_{g(i)}) = g((a_i) + (b_i)).$$
Over other fields, it is also required that the map be linear, which
as per above boils down to scalar multiplication. However, over
$GF(2),$ the only scalars are 0 and 1, so the linearity condition has
trivial effect.

AUTHOR:
    Robert L Miller (Oct-Nov 2007)
        * Compiled code datastructure
        * union-find based orbit partition
        * optimized partition stack class

"""

#*****************************************************************************
#         Copyright (C) 2007 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include '../ext/cdefs.pxi'
include '../ext/python_mem.pxi'
include '../ext/stdsage.pxi'
from sage.structure.element import is_Matrix
from sage.misc.misc import cputime
from math import log, ceil
from sage.rings.integer import Integer

cdef class BinaryCode:
    """
    Minimal, but optimized, binary code object.

    EXAMPLES:
        sage: import sage.coding.binary_code
        sage: from sage.coding.binary_code import *
        sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
        sage: B = BinaryCode(M)
        sage: B
        Binary [8,4] linear code, generator matrix
        [11110000]
        [00111100]
        [00001111]
        [10101010]
        sage: B._is_one(7, 4)
        0
        sage: B._is_one(8, 4)
        1
        sage: B._is_automorphism([1,0,3,2,4,5,6,7], [1, 2, 4, 9])
        1

    """
    def __new__(self, arg1):
        """
        Create binary code from matrix. Input is assumed to be a matrix
        over $GF(2)$.

        EXAMPLE:
        sage: import sage.coding.binary_code
        sage: from sage.coding.binary_code import *
        sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
        sage: B = BinaryCode(M)
        sage: B
        Binary [8,4] linear code, generator matrix
        [11110000]
        [00111100]
        [00001111]
        [10101010]

        """
        cdef int i, j
        cdef unsigned int k
        self.radix = 8*sizeof(unsigned int)
        self.ncols = arg1.ncols()
        self.nrows = arg1.nrows()
        self.nwords = 1 << self.nrows

        if self.nrows > self.radix or ceil(log(self.ncols,2)) > self.radix:
            raise NotImplementedError("Columns and rows are stored as ints. This code is too big.")
        self.columns = <unsigned int *> sage_malloc( self.ncols * sizeof(unsigned int) )
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

    def __repr__(self):
        """
        String representation of self.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [8,4] linear code, generator matrix
            [11110000]
            [00111100]
            [00001111]
            [10101010]

        """
        cdef int i, j
        s = 'Binary [%d,%d] linear code, generator matrix\n'%(self.ncols, self.nrows)
        for i from 0 <= i < self.nrows:
            s += '['
            for j from 0 <= j < self.ncols:
                s += '%d'%self.is_one(1<<i,j)
            s += ']\n'
        return s

    def _is_one(self, word, col):
        """
        Returns the col-th letter of word, i.e. 0 or 1. Words are expressed
        as integers, which represent linear combinations of the rows of the
        generator matrix of the code.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [8,4] linear code, generator matrix
            [11110000]
            [00111100]
            [00001111]
            [10101010]
            sage: B._is_one(7, 4)
            0
            sage: B._is_one(8, 4)
            1
            sage: B._is_automorphism([1,0,3,2,4,5,6,7], [1, 2, 4, 9])
            1

        """
        return self.is_one(word, col)

    cdef int is_one(self, unsigned int word, int column):
        cdef int i, j
        cdef unsigned int k
        i = 0
        k = word & self.columns[column]
        for j from 0 <= j < self.radix by 4:
            i += self.ham_wts[15 & k]
            k = k >> 4
        return i % 2

    def _is_automorphism(self, col_gamma, basis_gamma):
        """
        Check whether a given permutation is an automorphism of the code.

        INPUT:
            col_gamma -- permutation sending i |--> col_gamma[i] acting
                on the columns.
            basis_gamma -- describes where the basis elements are mapped
                under gamma. basis_gamma[i] is where the i-th row is sent,
                as an integer expressing a linear combination of the rows
                of the generator matrix.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [8,4] linear code, generator matrix
            [11110000]
            [00111100]
            [00001111]
            [10101010]
            sage: B._is_automorphism([1,0,3,2,4,5,6,7], [1, 2, 4, 9])
            1

        """
        cdef int i
        cdef int *_col_gamma
        cdef unsigned int *_basis_gamma
        _basis_gamma = <unsigned int *> sage_malloc(self.nrows * sizeof(unsigned int))
        _col_gamma = <int *> sage_malloc(self.ncols * sizeof(int))
        if not (_col_gamma and _basis_gamma):
            if _basis_gamma: sage_free(_basis_gamma)
            if _col_gamma: sage_free(_col_gamma)
            raise MemoryError("Memory.")
        for i from 0 <= i < self.nrows:
            _basis_gamma[i] = basis_gamma[i]
        for i from 0 <= i < self.ncols:
            _col_gamma[i] = col_gamma[i]
        result = self.is_automorphism(_col_gamma, _basis_gamma)
        sage_free(_col_gamma)
        sage_free(_basis_gamma)
        return result

    cdef int is_automorphism(self, int *col_gamma, unsigned int *basis_gamma):
        cdef int i, j
        for i from 0 <= i < self.nrows:
            for j from 0 <= j < self.ncols:
                if self.is_one(1 << i, j) != self.is_one(basis_gamma[i], col_gamma[j]):
                    return 0
        return 1

cdef class OrbitPartition:
    """
    Structure which keeps track of which vertices are equivalent
    under the part of the automorphism group that has already been
    seen, during search. Essentially a disjoint-set data structure*,
    which also keeps track of the minimum element and size of each
    cell of the partition, and the size of the partition.

    * http://en.wikipedia.org/wiki/Disjoint-set_data_structure

    """
    def __new__(self, nrows, ncols):
        cdef unsigned int nwords, word
        cdef int col
        nwords = (1 << nrows)
        self.wd_parent =       <unsigned int *> sage_malloc( nwords * sizeof(unsigned int) )
        self.wd_rank =         <unsigned int *> sage_malloc( nwords * sizeof(unsigned int) )
        self.wd_min_cell_rep = <unsigned int *> sage_malloc( nwords * sizeof(unsigned int) )
        self.wd_size =         <unsigned int *> sage_malloc( nwords * sizeof(unsigned int) )
        self.col_parent =       <int *> sage_malloc( ncols * sizeof(int) )
        self.col_rank =         <int *> sage_malloc( ncols * sizeof(int) )
        self.col_min_cell_rep = <int *> sage_malloc( ncols * sizeof(int) )
        self.col_size =         <int *> sage_malloc( ncols * sizeof(int) )
        if not (self.wd_parent and self.wd_rank and self.wd_min_cell_rep and self.wd_size and self.col_parent and self.col_rank and self.col_min_cell_rep and self.col_size):
            if self.wd_parent: sage_free(self.wd_parent)
            if self.wd_rank: sage_free(self.wd_rank)
            if self.wd_min_cell_rep: sage_free(self.wd_min_cell_rep)
            if self.wd_size: sage_free(self.wd_size)
            if self.col_parent: sage_free(self.col_parent)
            if self.col_rank: sage_free(self.col_rank)
            if self.col_min_cell_rep: sage_free(self.col_min_cell_rep)
            if self.col_size: sage_free(self.col_size)
            raise MemoryError("Memory.")
        for word from 0 <= word < nwords:
            self.wd_parent[word] = word
            self.wd_rank[word] = 0
            self.wd_min_cell_rep[word] = word
            self.wd_size[word] = 1
        for col from 0 <= col < ncols:
            self.col_parent[col] = col
            self.col_rank[col] = 0
            self.col_min_cell_rep[col] = col
            self.col_size[col] = 1

    def __dealloc__(self):
        sage_free(self.wd_parent)
        sage_free(self.wd_rank)
        sage_free(self.wd_min_cell_rep)
        sage_free(self.wd_size)
        sage_free(self.col_parent)
        sage_free(self.col_rank)
        sage_free(self.col_min_cell_rep)
        sage_free(self.col_size)

    cdef unsigned int wd_find(self, unsigned int word):
        if self.wd_parent[word] == word:
            return word
        else:
            self.wd_parent[word] = self.wd_find(self.wd_parent[word])
            return self.wd_parent[word]

    cdef void wd_union(self, unsigned int x, unsigned int y):
        cdef unsigned int x_root, y_root
        x_root = self.wd_find(x)
        y_root = self.wd_find(y)
        if self.wd_rank[x_root] > self.wd_rank[y_root]:
            self.wd_parent[y_root] = x_root
            self.wd_min_cell_rep[y_root] = min(self.wd_min_cell_rep[x_root],self.wd_min_cell_rep[y_root])
            self.wd_size[y_root] += self.wd_size[x_root]
        elif self.wd_rank[x_root] < self.wd_rank[y_root]:
            self.wd_parent[x_root] = y_root
            self.wd_min_cell_rep[x_root] = min(self.wd_min_cell_rep[x_root],self.wd_min_cell_rep[y_root])
            self.wd_size[x_root] += self.wd_size[y_root]
        elif x_root != y_root:
            self.wd_parent[y_root] = x_root
            self.wd_min_cell_rep[y_root] = min(self.wd_min_cell_rep[x_root],self.wd_min_cell_rep[y_root])
            self.wd_size[y_root] += self.wd_size[x_root]
            self.wd_rank[x_root] += 1

    cdef int col_find(self, int col):
        if self.col_parent[col] == col:
            return col
        else:
            self.col_parent[col] = self.col_find(self.col_parent[col])
            return self.col_parent[col]

    cdef void col_union(self, int x, int y):
        cdef int x_root, y_root
        x_root = self.col_find(x)
        y_root = self.col_find(y)
        if self.col_rank[x_root] > self.col_rank[y_root]:
            self.col_parent[y_root] = x_root
            self.col_min_cell_rep[y_root] = min(self.col_min_cell_rep[x_root],self.col_min_cell_rep[y_root])
            self.col_size[y_root] += self.col_size[x_root]
        elif self.col_rank[x_root] < self.col_rank[y_root]:
            self.col_parent[x_root] = y_root
            self.col_min_cell_rep[x_root] = min(self.col_min_cell_rep[x_root],self.col_min_cell_rep[y_root])
            self.col_size[x_root] += self.col_size[y_root]
        elif x_root != y_root:
            self.col_parent[y_root] = x_root
            self.col_min_cell_rep[y_root] = min(self.col_min_cell_rep[x_root],self.col_min_cell_rep[y_root])
            self.col_size[y_root] += self.col_size[x_root]
            self.col_rank[x_root] += 1

    cdef int merge_perm(self, int *col_gamma, unsigned int *wd_gamma, int nrows, int ncols):
        """
        Returns 0 if nothing was done, otherwise returns 1.

        If gamma[a] = b, then after merge_perm, a and b will be in the same cell.

        """
        cdef unsigned int i, gamma_i_root
        cdef int j, gamma_j_root, return_value = 0
        for i from 0 <= i < (1 << nrows):
            if self.wd_parent[i] == i:
                gamma_i_root = self.wd_find(wd_gamma[i])
                if gamma_i_root != i:
                    return_value = 1
                    self.wd_union(i, gamma_i_root)
        for j from 0 <= j < ncols:
            if self.col_parent[j] == j:
                gamma_j_root = self.col_find(col_gamma[j])
                if gamma_j_root != j:
                    return_value = 1
                    self.col_union(j, gamma_j_root)
        return return_value

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
            raise MemoryError("Memory.")
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

    cdef int col_degree(self, BinaryCode CG, int col, int wd_ptr, int k):
        cdef int i = 0
        col = self.col_ents[col]
        while True:
            if CG.is_one(wd_ptr, col): i += 1
            if self.wd_lvls[wd_ptr] > k: wd_ptr += 1
            else: break
        return i

    cdef int wd_degree(self, BinaryCode CG, int wd, int col_ptr, int k):
        cdef int i = 0
        wd = self.wd_ents[wd]
        while True:
            if CG.is_one(wd, col_ptr): i += 1
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

    cdef int refine(self, int k, int *col_alpha, int *wd_alpha, BinaryCode CG):
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

    cdef int cmp(self, PartitionStack other, BinaryCode CG):
        # if CG(self) > G(other): return 1
        # if CG(self) < G(other): return -1
        # else: return 0
        cdef int i, j, k, l
        for i from 0 <= i < CG.nwords:
            for j from 0 <= j < CG.ncols:
                k = CG.is_one(self.wd_ents[i], self.col_ents[j])
                l = CG.is_one(other.wd_ents[i], other.col_ents[j])
                if k ^ l:
                    return k - l
        return 0

def classify(BinaryCode C, lab=True, verbosity=0):
    """
    """
    cdef int i, j # local variables
    cdef OrbitPartition Theta # keeps track of which vertices have been
                              # discovered to be equivalent
    cdef int index = 0, size = 1
    cdef int L = 100
    cdef int **Phi, **Omega
    cdef int l = -1
    cdef PartitionStack nu, zeta, rho
    cdef int k_rho
    cdef int k = 0
    cdef int h = -1
    cdef int hb
    cdef int hh = 1
    cdef int ht
    cdef mpz_t *Lambda_mpz, *zf_mpz, *zb_mpz
    cdef int hzf
    cdef int hzb = -1
    cdef unsigned int *basis_gamma
    cdef int *col_gamma
    cdef int *alpha
    cdef int *v
    cdef int *e
    cdef int state
    cdef int tvc, tvh, nwords = C.nwords, ncols = C.ncols, n = nwords + ncols, nrows = C.nrows

    # trivial case
    if ncols == 0:
        return [], {}
    elif nwords == 0 and ncols == 1:
        return [], {0:0}
    elif nwords == 0:
        output1 = []
        dd = {}
        for i from 0 <= i < ncols-1:
            dd[i] = i
            perm = range(ncols)
            perm[i] = i+1
            perm[i+1] = i
            output1.append(perm)
        dd[ncols-1] = ncols-1
        return output1, dd

    # allocate int pointers
    Phi = <int **> sage_malloc(L * sizeof(int *))
    Omega = <int **> sage_malloc(L * sizeof(int *))

    # allocate GMP int pointers
    Lambda_mpz = <mpz_t *> sage_malloc((ncols+2)*sizeof(mpz_t))
    zf_mpz     = <mpz_t *> sage_malloc((ncols+2)*sizeof(mpz_t))
    zb_mpz     = <mpz_t *> sage_malloc((ncols+2)*sizeof(mpz_t))

    # check for memory errors
    if not (Phi and Omega and Lambda_mpz and zf_mpz and zb_mpz):
        if Lambda_mpz: sage_free(Lambda_mpz)
        if zf_mpz: sage_free(zf_mpz)
        if zb_mpz: sage_free(zb_mpz)
        if Phi: sage_free(Phi)
        if Omega: sage_free(Omega)
        raise MemoryError("Error allocating memory.")

    # allocate int arrays
    basis_gamma = <unsigned int *> sage_malloc(nrows*sizeof(unsigned int))
    col_gamma = <int *> sage_malloc(ncols*sizeof(int))
    Phi[0] = <int *> sage_malloc(L*ncols*sizeof(int))
    Omega[0] = <int *> sage_malloc(L*ncols*sizeof(int))
    alpha = <int *> sage_malloc(4*ncols*sizeof(int))
    v = <int *> sage_malloc(ncols*sizeof(int))
    e = <int *> sage_malloc(ncols*sizeof(int))

    # check for memory errors
    if not (basis_gamma and col_gamma and Phi[0] and Omega[0] and alpha and v and e):
        if basis_gamma: sage_free(basis_gamma)
        if col_gamma: sage_free(col_gamma)
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