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
        * compiled code datastructure
        * union-find based orbit partition
        * optimized partition stack class

        * NICE-based partition refinement algorithm
        * canonical generation function

"""

#*******************************************************************************
#         Copyright (C) 2007 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

include '../ext/cdefs.pxi'
include '../ext/python_mem.pxi'
include '../ext/stdsage.pxi'
from sage.structure.element import is_Matrix
from sage.misc.misc import cputime
from math import log, ceil
from sage.rings.integer import Integer

## NOTE - Since most of the functions are used from within the module, cdef'd
## functions come without an underscore, and the def'd equivalents, which are
## essentially only for doctesting and debugging, have underscores.

cdef class BinaryCode:
    """
    Minimal, but optimized, binary code object.

    """
    def __new__(self, arg1):
        cdef unsigned int parity, word, combination
        cdef int nrows, i, j
        cdef unsigned int *self_words, *self_basis
        self.radix = 8*sizeof(unsigned int)
        self.ncols = arg1.ncols()
        self.nrows = arg1.nrows()
        self.nwords = 1 << self.nrows

        if self.nrows > self.radix or self.ncols > self.radix:
            raise NotImplementedError("Columns and rows are stored as ints. This code is too big.")

        self.words = <unsigned int *> sage_malloc( self.nwords * sizeof(unsigned int) )
        self.basis = <unsigned int *> sage_malloc( self.nrows * sizeof(unsigned int) )
        if not self.words or not self.basis:
            if self.words: sage_free(self.words)
            if self.basis: sage_free(self.basis)
            raise MemoryError("Memory.")

        nrows = self.nrows
        self_words = self.words
        self_basis = self.basis

        rows = arg1.rows()
        for i from 0 <= i < nrows:
            word = 0
            for j in rows[i].nonzero_positions():
                word += (1<<j)
            self_basis[i] = word

        word = 0
        parity = 0
        combination = 0
        while True:
            self_words[combination] = word
            parity ^= 1
            j = 0
            if not parity:
                while not combination & (1 << j): j += 1
                j += 1
            if j == nrows: break
            else:
                combination ^= (1 << j)
                word ^= self_basis[j]
        # NOTE: when implementing construction from parent
        # matrices, remember to grab this data from the parent

    def __dealloc__(self):
        sage_free(self.words)
        sage_free(self.basis)

    def delete(self):
        sage_free(self.words)
        sage_free(self.basis)

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
        if self.words[word] & (1 << column):
            return 1
        else:
            return 0

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
        cdef int i, j, self_nrows = self.nrows, self_ncols = self.ncols
        for i from 0 <= i < self_nrows:
            for j from 0 <= j < self_ncols:
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
        self.nwords = nwords
        self.ncols = ncols
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

    def delete(self):
        sage_free(self.wd_parent)
        sage_free(self.wd_rank)
        sage_free(self.wd_min_cell_rep)
        sage_free(self.wd_size)
        sage_free(self.col_parent)
        sage_free(self.col_rank)
        sage_free(self.col_min_cell_rep)
        sage_free(self.col_size)

    def __repr__(self):
        """
        Return a string representation of the orbit partition.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7

        """
        cdef unsigned int i
        cdef int j
        s = 'OrbitPartition on %d words and %d columns. Data:\nWords:\n'%(self.nwords, self.ncols)
        for i from 0 <= i < self.nwords:
            s += '%d,'%self.wd_parent[i]
        s = s[:-1] + '\nColumns:\n'
        for j from 0 <= j < self.ncols:
            s += '%d,'%self.col_parent[j]
        return s[:-1]

    def _wd_find(self, word):
        """
        Returns the root of word.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._wd_find(12)
            12L

        """
        return self.wd_find(word)

    cdef unsigned int wd_find(self, unsigned int word):
        if self.wd_parent[word] == word:
            return word
        else:
            self.wd_parent[word] = self.wd_find(self.wd_parent[word])
            return self.wd_parent[word]

    def _wd_union(self, x, y):
        """
        Join the cells containing x and y.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._wd_union(1,10)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,1,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._wd_find(10)
            1L

        """
        self.wd_union(x, y)

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

    def _col_find(self, col):
        """
        Returns the root of col.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._col_find(6)
            6

        """
        return self.col_find(col)

    cdef int col_find(self, int col):
        if self.col_parent[col] == col:
            return col
        else:
            self.col_parent[col] = self.col_find(self.col_parent[col])
            return self.col_parent[col]

    def _col_union(self, x, y):
        """
        Join the cells containing x and y.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._col_union(1,4)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,1,5,6,7
            sage: O._col_find(4)
            1

        """
        self.col_union(x, y)

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

    def _merge_perm(self, col_gamma, wd_gamma):
        """
        Merges the cells of self under the given permutation. If gamma[a] = b,
        then after merge_perm, a and b will be in the same cell. Returns 0 if
        nothing was done, otherwise returns 1.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: O = OrbitPartition(4, 8)
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
            Columns:
            0,1,2,3,4,5,6,7
            sage: O._merge_perm([1,0,3,2,4,5,6,7], [0,1,2,3,4,5,6,7,9,8,11,10,13,12,15,14])
            1
            sage: O
            OrbitPartition on 16 words and 8 columns. Data:
            Words:
            0,1,2,3,4,5,6,7,8,8,10,10,12,12,14,14
            Columns:
            0,0,2,2,4,5,6,7

        """
        cdef unsigned int i
        cdef int *_col_gamma
        cdef unsigned int *_wd_gamma
        _wd_gamma = <unsigned int *> sage_malloc(self.nwords * sizeof(unsigned int))
        _col_gamma = <int *> sage_malloc(self.ncols * sizeof(int))
        if not (_col_gamma and _wd_gamma):
            if _wd_gamma: sage_free(_wd_gamma)
            if _col_gamma: sage_free(_col_gamma)
            raise MemoryError("Memory.")
        for i from 0 <= i < self.nwords:
            _wd_gamma[i] = wd_gamma[i]
        for i from 0 <= i < self.ncols:
            _col_gamma[i] = col_gamma[i]
        result = self.merge_perm(_col_gamma, _wd_gamma)
        sage_free(_col_gamma)
        sage_free(_wd_gamma)
        return result

    cdef int merge_perm(self, int *col_gamma, unsigned int *wd_gamma):
        cdef unsigned int i, gamma_i_root
        cdef int j, gamma_j_root, return_value = 0
        cdef unsigned int *self_wd_parent = self.wd_parent
        cdef int *self_col_parent = self.col_parent
        for i from 0 <= i < self.nwords:
            if self_wd_parent[i] == i:
                gamma_i_root = self.wd_find(wd_gamma[i])
                if gamma_i_root != i:
                    return_value = 1
                    self.wd_union(i, gamma_i_root)
        for j from 0 <= j < self.ncols:
            if self_col_parent[j] == j:
                gamma_j_root = self.col_find(col_gamma[j])
                if gamma_j_root != j:
                    return_value = 1
                    self.col_union(j, gamma_j_root)
        return return_value

cdef class PartitionStack:
    """
    Partition stack structure for traversing the search tree during automorphism
    group computation.

    """
    def __new__(self, arg1, arg2=None):
        cdef int k
        self.nwords = int(arg1)
        self.ncols = int(arg2)
        self.radix      = 8*sizeof(unsigned int)

        # data
        self.wd_ents    = <unsigned int *> sage_malloc( self.nwords * sizeof(unsigned int) )
        self.wd_lvls    = <int *> sage_malloc( self.nwords * sizeof(int) )
        self.col_ents   = <int *> sage_malloc( self.ncols * sizeof(int) )
        self.col_lvls   = <int *> sage_malloc( self.ncols * sizeof(int) )

        # scratch space
        self.col_degs   = <unsigned int *> sage_malloc( self.ncols * sizeof(unsigned int) )
        self.col_counts = <int *> sage_malloc( self.nwords * sizeof(int) )
        self.col_output = <int *> sage_malloc( self.ncols * sizeof(int) )
        self.wd_degs    = <int *> sage_malloc( self.nwords * sizeof(int) )
        self.wd_counts  = <unsigned int *> sage_malloc( self.ncols * sizeof(unsigned int) )
        self.wd_output  = <unsigned int *> sage_malloc( self.nwords * sizeof(unsigned int) )

        if not (self.wd_ents and self.wd_lvls and self.col_ents and self.col_lvls \
            and self.col_degs and self.col_counts and self.col_output \
            and self.wd_degs and self.wd_counts and self.wd_output):
            if self.wd_ents: sage_free(self.wd_ents)
            if self.wd_lvls: sage_free(self.wd_lvls)
            if self.col_ents: sage_free(self.col_ents)
            if self.col_lvls: sage_free(self.col_lvls)
            if self.col_degs: sage_free(self.col_degs)
            if self.col_counts: sage_free(self.col_counts)
            if self.col_output: sage_free(self.col_output)
            if self.wd_degs: sage_free(self.wd_degs)
            if self.wd_counts: sage_free(self.wd_counts)
            if self.wd_output: sage_free(self.wd_output)
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
        sage_free(self.col_degs)
        sage_free(self.col_counts)
        sage_free(self.col_output)
        sage_free(self.wd_degs)
        sage_free(self.wd_counts)
        sage_free(self.wd_output)

    def delete(self):
        sage_free(self.wd_ents)
        sage_free(self.wd_lvls)
        sage_free(self.col_ents)
        sage_free(self.col_lvls)
        sage_free(self.col_degs)
        sage_free(self.col_counts)
        sage_free(self.col_output)
        sage_free(self.wd_degs)
        sage_free(self.wd_counts)
        sage_free(self.wd_output)

    def __repr__(self):
        """
        Return a string representation of self.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})

        """
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

    def _is_discrete(self, k):
        """
        Returns whether the partition at level k is discrete.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: [P._split_column(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3,4,5})
            ({0},{1,2,3,4,5})
            ({0},{1},{2,3,4,5})
            ({0},{1},{2},{3,4,5})
            ({0},{1},{2},{3},{4,5})
            ({0},{1},{2},{3},{4},{5})
            sage: P._is_discrete(4)
            0
            sage: P._is_discrete(5)
            1

        """
        return self.is_discrete(k)

    cdef int is_discrete(self, int k):
        cdef int i
        cdef int *self_col_lvls = self.col_lvls
        for i from 0 <= i < self.ncols:
            if self_col_lvls[i] > k:
                return 0
        return 1

    def _num_cells(self, k):
        """
        Returns the number of cells in the partition at level k.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: [P._split_column(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3,4,5})
            ({0},{1,2,3,4,5})
            ({0},{1},{2,3,4,5})
            ({0},{1},{2},{3,4,5})
            ({0},{1},{2},{3},{4,5})
            ({0},{1},{2},{3},{4},{5})
            sage: P._num_cells(3)
            5

        """
        return self.num_cells(k)

    cdef int num_cells(self, int k):
        cdef int i, j = 0
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_col_lvls = self.col_lvls
        for i from 0 <= i < self.nwords:
            if self_wd_lvls[i] <= k:
                j += 1
        for i from 0 <= i < self.ncols:
            if self_col_lvls[i] <= k:
                j += 1
        return j

    def _sat_225(self, k):
        """
        Returns whether the partition at level k satisfies the hypotheses of
        Lemma 2.25 in Brendan McKay's Practical Graph Isomorphism paper (see
        sage/graphs/graph_isom.pyx.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: [P._split_column(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P._sat_225(3)
            0
            sage: P._sat_225(4)
            1
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3,4,5})
            ({0},{1,2,3,4,5})
            ({0},{1},{2,3,4,5})
            ({0},{1},{2},{3,4,5})
            ({0},{1},{2},{3},{4,5})
            ({0},{1},{2},{3},{4},{5})

        """
        return self.sat_225(k)

    cdef int sat_225(self, int k):
        cdef int i, n = self.nwords + self.ncols, in_cell = 0
        cdef int nontrivial_cells = 0, total_cells = self.num_cells(k)
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_col_lvls = self.col_lvls
        if n <= total_cells + 4:
            return 1
        for i from 0 <= i < self.nwords:
            if self_wd_lvls[i] <= k:
                if in_cell:
                    nontrivial_cells += 1
                in_cell = 0
            else:
                in_cell = 1
        in_cell = 0
        for i from 0 <= i < self.ncols:
            if self_col_lvls[i] <= k:
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

    def _min_cell_reps(self, k):
        """
        Returns an integer whose bits represent which columns are minimal cell
        representatives.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: [P._split_column(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: a = P._min_cell_reps(2)
            sage: Integer(a).binary()
            '111'
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3,4,5})
            ({0},{1,2,3,4,5})
            ({0},{1},{2,3,4,5})
            ({0},{1},{2},{3,4,5})
            ({0},{1},{2},{3},{4,5})
            ({0},{1},{2},{3},{4},{5})

        """
        return self.min_cell_reps(k)

    cdef unsigned int min_cell_reps(self, int k):
        cdef int i
        cdef unsigned int reps = 1
        cdef int *self_col_lvls = self.col_lvls
        for i from 0 < i < self.ncols:
            if self_col_lvls[i-1] <= k:
                reps += (1 << i)
        return reps

    def _fixed_cols(self, mcrs, k):
        """
        Returns an integer whose bits represent which columns are fixed. For
        efficiency, mcrs is the output of min_cell_reps.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: [P._split_column(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: a = P._fixed_cols(7, 2)
            sage: Integer(a).binary()
            '11'
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3,4,5})
            ({0},{1,2,3,4,5})
            ({0},{1},{2,3,4,5})
            ({0},{1},{2},{3,4,5})
            ({0},{1},{2},{3},{4,5})
            ({0},{1},{2},{3},{4},{5})

        """
        return self.fixed_cols(mcrs, k)

    cdef unsigned int fixed_cols(self, unsigned int mcrs, int k):
        cdef int i
        cdef unsigned int fixed = 0
        cdef int *self_col_lvls = self.col_lvls
        for i from 0 <= i < self.ncols:
            if self_col_lvls[i] <= k:
                fixed += (1 << i)
        return fixed & mcrs

    def _first_smallest_nontrivial(self, k):
        """
        Returns an integer representing the first, smallest nontrivial cell.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: [P._split_column(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: a = P._first_smallest_nontrivial(2)
            sage: Integer(a).binary().zfill(32)
            '00000000000000000000000000111100'
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3,4,5})
            ({0},{1,2,3,4,5})
            ({0},{1},{2,3,4,5})
            ({0},{1},{2},{3,4,5})
            ({0},{1},{2},{3},{4,5})
            ({0},{1},{2},{3},{4},{5})

        """
        return self.first_smallest_nontrivial(k)

    cdef unsigned int first_smallest_nontrivial(self, int k):
        cdef unsigned int cell
        cdef int i = 0, j = 0, location = 0, ncols = self.ncols
        cdef int *self_col_lvls = self.col_lvls
        while True:
            if self_col_lvls[i] <= k:
                if i != j and ncols > i - j + 1:
                    ncols = i - j + 1
                    location = j
                j = i + 1
            if self_col_lvls[i] == -1: break
            i += 1
        # location now points to the beginning of the first, smallest,
        # nontrivial cell
        j = location
        while True:
            if self_col_lvls[j] <= k: break
            j += 1
        # j now points to the last element of the cell
        i = self.radix - j - 1                 # the cell is represented in binary, reading from the right:
        cell = (~0 << location) ^ (~0 << j+1)  # <-------            self.radix               ----->
        return cell                            # [0]*(radix-j-1) + [1]*(j-location+1) + [0]*location

    def _dangerous_dont_use_set_ents_lvls(self, col_ents, col_lvls, wd_ents, wd_lvls):
        """
        For debugging only.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            sage: P._dangerous_dont_use_set_ents_lvls([99]*6, [0,3,2,3,5,-1], [4,3,5,6], [3,2,1,-1])
            sage: P
            ({4,3,5,6})
            ({4,3,5},{6})
            ({4,3},{5},{6})
            ({4},{3},{5},{6})
            ({99},{99,99,99,99,99})
            ({99},{99,99,99,99,99})
            ({99},{99,99},{99,99,99})
            ({99},{99},{99},{99},{99,99})
            ({99},{99},{99},{99},{99,99})
            ({99},{99},{99},{99},{99},{99})

        """
        cdef unsigned int i
        for i from 0 <= i < len(col_ents):
            self.col_ents[i] = col_ents[i]
        for i from 0 <= i < len(col_lvls):
            self.col_lvls[i] = col_lvls[i]
        for i from 0 <= i < len(wd_ents):
            self.wd_ents[i] = wd_ents[i]
        for i from 0 <= i < len(wd_lvls):
            self.wd_lvls[i] = wd_lvls[i]

    def _col_percolate(self, start, end):
        """
        Do one round of bubble sort on ents.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: P._dangerous_dont_use_set_ents_lvls(range(5,-1,-1), [1,2,2,3,3,-1], range(3,-1,-1), [1,1,2,-1])
            sage: P
            ({3,2,1,0})
            ({3},{2},{1,0})
            ({3},{2},{1},{0})
            ({3},{2},{1},{0})
            ({5,4,3,2,1,0})
            ({5},{4,3,2,1,0})
            ({5},{4},{3},{2,1,0})
            ({5},{4},{3},{2},{1},{0})
            ({5},{4},{3},{2},{1},{0})
            ({5},{4},{3},{2},{1},{0})
            sage: P._wd_percolate(0,3)
            sage: P._col_percolate(0,5)
            sage: P
            ({0,3,2,1})
            ({0},{3},{2,1})
            ({0},{3},{2},{1})
            ({0},{3},{2},{1})
            ({0,5,4,3,2,1})
            ({0},{5,4,3,2,1})
            ({0},{5},{4},{3,2,1})
            ({0},{5},{4},{3},{2},{1})
            ({0},{5},{4},{3},{2},{1})
            ({0},{5},{4},{3},{2},{1})

        """
        self.col_percolate(start, end)

    cdef void col_percolate(self, int start, int end):
        cdef int i, temp
        cdef int *self_col_ents = self.col_ents
        for i from end >= i > start:
            if self_col_ents[i] < self_col_ents[i-1]:
                temp = self.col_ents[i]
                self_col_ents[i] = self_col_ents[i-1]
                self_col_ents[i-1] = temp

    def _wd_percolate(self, start, end):
        """
        Do one round of bubble sort on ents.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: P._dangerous_dont_use_set_ents_lvls(range(5,-1,-1), [1,2,2,3,3,-1], range(3,-1,-1), [1,1,2,-1])
            sage: P
            ({3,2,1,0})
            ({3},{2},{1,0})
            ({3},{2},{1},{0})
            ({3},{2},{1},{0})
            ({5,4,3,2,1,0})
            ({5},{4,3,2,1,0})
            ({5},{4},{3},{2,1,0})
            ({5},{4},{3},{2},{1},{0})
            ({5},{4},{3},{2},{1},{0})
            ({5},{4},{3},{2},{1},{0})
            sage: P._wd_percolate(0,3)
            sage: P._col_percolate(0,5)
            sage: P
            ({0,3,2,1})
            ({0},{3},{2,1})
            ({0},{3},{2},{1})
            ({0},{3},{2},{1})
            ({0,5,4,3,2,1})
            ({0},{5,4,3,2,1})
            ({0},{5},{4},{3,2,1})
            ({0},{5},{4},{3},{2},{1})
            ({0},{5},{4},{3},{2},{1})
            ({0},{5},{4},{3},{2},{1})

        """
        self.wd_percolate(start, end)

    cdef void wd_percolate(self, unsigned int start, unsigned int end):
        cdef unsigned int i, temp
        cdef unsigned int *self_wd_ents = self.wd_ents
        for i from end >= i > start:
            if self_wd_ents[i] < self_wd_ents[i-1]:
                temp = self.wd_ents[i]
                self_wd_ents[i] = self_wd_ents[i-1]
                self_wd_ents[i-1] = temp

    def _split_column(self, int v, int k):
        """
        Split column v out, placing it before the rest of the cell it was in.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: [P._split_column(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3,4,5})
            ({0},{1,2,3,4,5})
            ({0},{1},{2,3,4,5})
            ({0},{1},{2},{3,4,5})
            ({0},{1},{2},{3},{4,5})
            ({0},{1},{2},{3},{4},{5})
            sage: P = PartitionStack(4, 6)
            sage: P._split_column(0,1)
            0
            sage: P._split_column(2,2)
            1
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,2,1,3,4,5})
            ({0},{2,1,3,4,5})
            ({0},{2},{1,3,4,5})
            ({0},{2},{1,3,4,5})
            ({0},{2},{1,3,4,5})
            ({0},{2},{1,3,4,5})

        """
        return self.split_column(v, k)

    cdef int split_column(self, int v, int k):
        cdef int i = 0, j
        cdef int *self_col_ents = self.col_ents
        cdef int *self_col_lvls = self.col_lvls
        while self_col_ents[i] != v: i += 1
        j = i
        while self_col_lvls[i] > k: i += 1
        if j == 0 or self_col_lvls[j-1] <= k:
            self.col_percolate(j+1, i)
        else:
            while j != 0 and self_col_lvls[j-1] > k:
                self_col_ents[j] = self_col_ents[j-1]
                j -= 1
            self_col_ents[j] = v
        self_col_lvls[j] = k
        return j

    def _col_degree(self, C, col, wd_ptr, k):
        """
        Returns the number of words in the cell specified by wd_ptr that have a
        1 in the col-th column.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0],[0,0,1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [6,2] linear code, generator matrix
            [111100]
            [001111]
            sage: [P._split_column(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3,4,5})
            ({0},{1,2,3,4,5})
            ({0},{1},{2,3,4,5})
            ({0},{1},{2},{3,4,5})
            ({0},{1},{2},{3},{4,5})
            ({0},{1},{2},{3},{4},{5})
            sage: P._col_degree(B, 2, 0, 2)
            2L

        """
        return self.col_degree(C, col, wd_ptr, k)

    cdef unsigned int col_degree(self, BinaryCode CG, int col, unsigned int wd_ptr, int k):
        cdef unsigned int i = 0
        cdef int *self_wd_lvls = self.wd_lvls
        col = self.col_ents[col]
        while True:
            if CG.is_one(wd_ptr, col): i += 1
            if self_wd_lvls[wd_ptr] > k: wd_ptr += 1
            else: break
        return i

    def _wd_degree(self, C, wd, col_ptr, k):
        """
        Returns the number of columns in the cell specified by col_ptr that are
        1 in wd.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0],[0,0,1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [6,2] linear code, generator matrix
            [111100]
            [001111]
            sage: [P._split_column(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3,4,5})
            ({0},{1,2,3,4,5})
            ({0},{1},{2,3,4,5})
            ({0},{1},{2},{3,4,5})
            ({0},{1},{2},{3},{4,5})
            ({0},{1},{2},{3},{4},{5})
            sage: P._wd_degree(B, 1, 1, 1)
            3

        """
        return self.wd_degree(C, wd, col_ptr, k)

    cdef int wd_degree(self, BinaryCode CG, unsigned int wd, int col_ptr, int k):
        cdef int i = 0
        wd = self.wd_ents[wd]
        cdef int *self_col_lvls = self.col_lvls
        while True:
            if CG.is_one(wd, col_ptr): i += 1
            if self_col_lvls[col_ptr] > k: col_ptr += 1
            else: break
        return i

    def _sort_cols(self, start, degrees, k):
        """
        Essentially a counting sort, but on only one cell of the partition.

        INPUT:
            start -- location of the beginning of the cell
            k -- at what level of refinement the partition of interest lies
            degrees -- the counts to sort by

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 6)
            sage: [P._split_column(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P._sort_cols(1, [0,2,2,1,1], 1)
            2
            sage: P
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,2,3})
            ({0,1,4,5,2,3})
            ({0},{1},{4,5},{2,3})
            ({0},{1},{4,5},{2,3})
            ({0},{1},{4},{5},{2,3})
            ({0},{1},{4},{5},{2,3})
            ({0},{1},{4},{5},{2},{3})

        """
        cdef int i
        for i from 0 <= i < len(degrees):
            self.col_degs[i] = degrees[i]
        return self.sort_cols(start, k)

    cdef int sort_cols(self, int start, int k):
        cdef int i, j, max, max_location, self_ncols = self.ncols
        cdef int *self_col_counts = self.col_counts
        cdef int *self_col_lvls = self.col_lvls
        cdef unsigned int *self_col_degs = self.col_degs
        cdef int *self_col_ents = self.col_ents
        cdef int *self_col_output = self.col_output
        for i from 0 <= i < self_ncols:
            self_col_counts[i] = 0
        i = 0
        while self_col_lvls[i+start] > k:
            self_col_counts[self_col_degs[i]] += 1
            i += 1
        self_col_counts[self_col_degs[i]] += 1

        # i+start is the right endpoint of the cell now
        max = self_col_counts[0]
        max_location = 0
        for j from 0 < j < self_ncols:
            if self_col_counts[j] > max:
                max = self_col_counts[j]
                max_location = j
            self_col_counts[j] += self_col_counts[j-1]

        for j from i >= j >= 0:
            self_col_counts[self_col_degs[j]] -= 1
            self_col_output[self_col_counts[self_col_degs[j]]] = self_col_ents[start+j]

        max_location = self_col_counts[max_location] + start

        for j from 0 <= j <= i:
            self_col_ents[start+j] = self_col_output[j]

        j = 1
        while j < self_ncols and self_col_counts[j] <= i:
            if self_col_counts[j] > 0:
                self_col_lvls[start + self_col_counts[j] - 1] = k
            self.col_percolate(start + self_col_counts[j-1], start + self_col_counts[j] - 1)
            j += 1

        return max_location

    def _sort_wds(self, start, degrees, k):
        """
        Essentially a counting sort, but on only one cell of the partition.

        INPUT:
            start -- location of the beginning of the cell
            k -- at what level of refinement the partition of interest lies
            degrees -- the counts to sort by

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(8, 6)
            sage: P._sort_wds(0, [0,0,3,3,3,3,2,2], 1)
            4L
            sage: P
            ({0,1,6,7,2,3,4,5})
            ({0,1},{6,7},{2,3,4,5})
            ({0,1},{6,7},{2,3,4,5})
            ({0,1},{6,7},{2,3,4,5})
            ({0,1},{6,7},{2,3,4,5})
            ({0,1},{6,7},{2,3,4,5})
            ({0,1},{6,7},{2,3,4,5})
            ({0,1},{6,7},{2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})
            ({0,1,2,3,4,5})

        """
        cdef int i
        for i from 0 <= i < len(degrees):
            self.wd_degs[i] = degrees[i]
        return self.sort_wds(start, k)

    cdef unsigned int sort_wds(self, unsigned int start, int k):
        cdef unsigned int i, j, max, max_location, self_nwords = self.nwords
        cdef unsigned int *self_wd_counts = self.wd_counts
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_wd_degs = self.wd_degs
        cdef unsigned int *self_wd_ents = self.wd_ents
        cdef unsigned int *self_wd_output = self.wd_output

        for i from 0 <= i < self_nwords:
            self_wd_counts[i] = 0
        i = 0
        while self_wd_lvls[i+start] > k:
            self_wd_counts[self_wd_degs[i]] += 1
            i += 1
        self_wd_counts[self_wd_degs[i]] += 1

        # i+start is the right endpoint of the cell now
        max = self_wd_counts[0]
        max_location = 0
        for j from 0 < j < self_nwords:
            if self_wd_counts[j] > max:
                max = self_wd_counts[j]
                max_location = j
            self_wd_counts[j] += self_wd_counts[j-1]

        for j from i >= j >= 0:
            if j > i: break # cython bug with unsigned ints...
            self_wd_counts[self_wd_degs[j]] -= 1
            self_wd_output[self_wd_counts[self_wd_degs[j]]] = self_wd_ents[start+j]

        max_location = self_wd_counts[max_location] + start

        for j from 0 <= j <= i:
            self_wd_ents[start+j] = self_wd_output[j]

        j = 1
        while j < self_nwords and self_wd_counts[j] <= i:
            if self_wd_counts[j] > 0:
                self_wd_lvls[start + self_wd_counts[j] - 1] = k
            self.wd_percolate(start + self_wd_counts[j-1], start + self_wd_counts[j] - 1)
            j += 1

        return max_location

################################################################################
################################################################################
################################################################################

    def _refine(self, k, col_alpha, wd_alpha, Code):
        """
        Refines the partition at level k, using the list of cells alpha, and Code.

        EXAMPLE:

        """
        cdef unsigned int i
        cdef int j
        cdef unsigned int *_col_a = <unsigned int *> sage_malloc(4*self.ncols*sizeof(unsigned int))
        cdef int *_wd_a = <int *> sage_malloc(4*self.nwords*sizeof(int))
        if not _col_a or not _wd_a:
            if _col_a: sage_free(_col_a)
            if _wd_a: sage_free(_wd_a)
            raise MemoryError("Memory.")
    # TODO       for i from  TODO
        result = self.refine(k, _col_a, _wd_a, Code)
        sage_free(_col_a)
        sage_free(_wd_a)
        return result

    cdef unsigned int refine(self, int k, unsigned int *col_alpha, int *wd_alpha, BinaryCode CG):
        cdef int m = 0, j
        cdef int i, q, r, s, t
        cdef unsigned int t_w
        cdef unsigned int invariant
        cdef unsigned int *col_degrees = col_alpha + self.ncols
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
#FIX                 t_w = self.sort_wds(j, wd_degrees, k)
                    invariant += t_w + wd_degrees[i-j-1]
                    q = m
                    while col_alpha[q] != -1:
                        if col_alpha[q] == j: col_alpha[q] = t_w
                        q += 1
                    r = j
                    while True:
                        if r == 0 or self.wd_lvls[r-1] == k:
                            if r != t_w:
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
#FIX                 t = self.sort_cols(j, col_degrees, k)
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
        # if CG(self) > CG(other): return 1
        # if CG(self) < CG(other): return -1
        # else: return 0
        cdef int i, j, k, l
        for i from 0 <= i < CG.nwords:
            for j from 0 <= j < CG.ncols:
                k = CG.is_one(self.wd_ents[i], self.col_ents[j])
                l = CG.is_one(other.wd_ents[i], other.col_ents[j])
                if k != l:
                    return k - l
        return 0

cdef class BinaryCodeClassifier:

    def __new__(self):
        cdef int *self_ham_wts
        self.ham_wts = <int *> sage_malloc( 65536 * sizeof(int) )
        if not self.ham_wts:
            sage_free(self.ham_wts)
        self_ham_wts = self.ham_wts
        self_ham_wts[0]  = 0; self_ham_wts[1]  = 1; self_ham_wts[2]  = 1; self_ham_wts[3]  = 2
        self_ham_wts[4]  = 1; self_ham_wts[5]  = 2; self_ham_wts[6]  = 2; self_ham_wts[7]  = 3
        self_ham_wts[8]  = 1; self_ham_wts[9]  = 2; self_ham_wts[10] = 2; self_ham_wts[11] = 3
        self_ham_wts[12] = 2; self_ham_wts[13] = 3; self_ham_wts[14] = 3; self_ham_wts[15] = 4
        for i from 16 <= i < 256:
            self_ham_wts[i] = self_ham_wts[i & 15] + self_ham_wts[(i>>4) & 15]
        for i from 256 <= i < 65536:
            self_ham_wts[i] = self_ham_wts[i & 255] + self_ham_wts[(i>>8) & 255]
        # This may seem like overkill, but the worst case for storing the words
        # themselves is 65536- in this case, we are increasing memory usage by a
        # factor of 2.

    def __dealloc__(self):
        sage_free(self.ham_wts)



















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