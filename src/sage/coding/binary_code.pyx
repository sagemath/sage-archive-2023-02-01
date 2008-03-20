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
from math import log, floor
from sage.rings.integer import Integer

## NOTE - Since most of the functions are used from within the module, cdef'd
## functions come without an underscore, and the def'd equivalents, which are
## essentially only for doctesting and debugging, have underscores.

cdef int *hamming_weights():
    cdef int *ham_wts
    ham_wts = <int *> sage_malloc( 65536 * sizeof(int) )
    if not ham_wts:
        sage_free(ham_wts)
        raise MemoryError("Memory.")
    ham_wts[0] = 0
    ham_wts[1] = 1
    ham_wts[2] = 1
    ham_wts[3] = 2
    for i from 4 <= i < 16:
        ham_wts[i] = ham_wts[i & 3] + ham_wts[(i>>2) & 3]
    for i from 16 <= i < 256:
        ham_wts[i] = ham_wts[i & 15] + ham_wts[(i>>4) & 15]
    for i from 256 <= i < 65536:
        ham_wts[i] = ham_wts[i & 255] + ham_wts[(i>>8) & 255]
    return ham_wts

cdef class BinaryCode:
    """
    Minimal, but optimized, binary code object.

    EXAMPLE:
        sage: import sage.coding.binary_code
        sage: from sage.coding.binary_code import *
        sage: M = Matrix(GF(2), [[1,1,1,1]])
        sage: B = BinaryCode(M)     # create from matrix
        sage: C = BinaryCode(B, 60) # create using glue
        sage: D = BinaryCode(C, 240)
        sage: E = BinaryCode(D, 85)
        sage: B
        Binary [4,1] linear code, generator matrix
        [1111]
        sage: C
        Binary [6,2] linear code, generator matrix
        [111100]
        [001111]
        sage: D
        Binary [8,3] linear code, generator matrix
        [11110000]
        [00111100]
        [00001111]
        sage: E
        Binary [8,4] linear code, generator matrix
        [11110000]
        [00111100]
        [00001111]
        [10101010]

    """
    def __new__(self, arg1, arg2=None):
        cdef int nrows, i, j
        cdef int nwords, other_nwords, parity, word, combination, glue_word
        cdef BinaryCode other
        cdef int *self_words, *self_basis, *other_basis

        self.radix = sizeof(int) << 3

        if is_Matrix(arg1):
            self.ncols = arg1.ncols()
            self.nrows = arg1.nrows()
            nrows = self.nrows
            self.nwords = 1 << nrows
            nwords = self.nwords
        elif isinstance(arg1, BinaryCode):
            other = arg1
            self.nrows = other.nrows + 1
            glue_word = arg2
            self.ncols = max( other.ncols , floor(log(glue_word,2))+1 )
            other_nwords = other.nwords
            self.nwords = 2 * other_nwords
            nrows = self.nrows
            nwords = self.nwords
        else: raise NotImplementedError("!")

        if self.nrows >= self.radix or self.ncols > self.radix:
            raise NotImplementedError("Columns and rows are stored as ints. This code is too big.")

        self.words = <int *> sage_malloc( nwords * sizeof(int) )
        self.basis = <int *> sage_malloc( nrows * sizeof(int) )
        if not self.words or not self.basis:
            if self.words: sage_free(self.words)
            if self.basis: sage_free(self.basis)
            raise MemoryError("Memory.")
        self_words = self.words
        self_basis = self.basis

        if is_Matrix(arg1):
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

        else: # isinstance(arg1, BinaryCode)
            other_basis = other.basis
            for i from 0 <= i < nrows-1:
                self_basis[i] = other_basis[i]
            i = nrows - 1
            self_basis[i] = glue_word

            memcpy(self_words, other.words, other_nwords*(self.radix>>3))

            for combination from 0 <= combination < other_nwords:
                self_words[combination+other_nwords] = self_words[combination] ^ glue_word

    def __dealloc__(self):
        sage_free(self.words)
        sage_free(self.basis)

    def print_data(self):
        """
        Print all data for self.

        EXAMPLES:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: C = BinaryCode(B, 60)
            sage: D = BinaryCode(C, 240)
            sage: E = BinaryCode(D, 85)
            sage: B.print_data() # random - actually "print P.print_data()"
            ncols: 4
            nrows: 1
            nwords: 2
            radix: 32
            basis:
            1111
            words:
            0000
            1111
            sage: C.print_data() # random - actually "print P.print_data()"
            ncols: 6
            nrows: 2
            nwords: 4
            radix: 32
            basis:
            111100
            001111
            words:
            000000
            111100
            001111
            110011
            sage: D.print_data() # random - actually "print P.print_data()"
            ncols: 8
            nrows: 3
            nwords: 8
            radix: 32
            basis:
            11110000
            00111100
            00001111
            words:
            00000000
            11110000
            00111100
            11001100
            00001111
            11111111
            00110011
            11000011
            sage: E.print_data() # random - actually "print P.print_data()"
            ncols: 8
            nrows: 4
            nwords: 16
            radix: 32
            basis:
            11110000
            00111100
            00001111
            10101010
            words:
            00000000
            11110000
            00111100
            11001100
            00001111
            11111111
            00110011
            11000011
            10101010
            01011010
            10010110
            01100110
            10100101
            01010101
            10011001
            01101001
        """
        from sage.graphs.graph_fast import binary
        cdef int ui
        cdef int i
        s = ''
        s += "ncols:" + str(self.ncols)
        s += "\nnrows:" + str(self.nrows)
        s += "\nnwords:" + str(self.nwords)
        s += "\nradix:" + str(self.radix)
        s += "\nbasis:\n"
        for i from 0 <= i < self.nrows:
            b = list(binary(self.basis[i]).zfill(self.ncols))
            b.reverse()
            b.append('\n')
            s += ''.join(b)
        s += "\nwords:\n"
        for ui from 0 <= ui < self.nwords:
            b = list(binary(self.words[ui]).zfill(self.ncols))
            b.reverse()
            b.append('\n')
            s += ''.join(b)

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
            sage: B._is_automorphism([1,0,3,2,4,5,6,7], [0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 11, 10, 13, 12, 15, 14])
            1

        """
        return self.is_one(word, col) != 0

    cdef int is_one(self, int word, int column):
        return (self.words[word] & (1 << column)) >> column

    def _is_automorphism(self, col_gamma, word_gamma):
        """
        Check whether a given permutation is an automorphism of the code.

        INPUT:
            col_gamma -- permutation sending i |--> col_gamma[i] acting
                on the columns.
            word_gamma -- permutation sending i |--> word_gamma[i] acting
                on the words.

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
            sage: B._is_automorphism([1,0,3,2,4,5,6,7], [0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 11, 10, 13, 12, 15, 14])
            1

        """
        cdef int i
        cdef int *_col_gamma
        cdef int *_word_gamma
        _word_gamma = <int *> sage_malloc(self.nwords * sizeof(int))
        _col_gamma = <int *> sage_malloc(self.ncols * sizeof(int))
        if not (_col_gamma and _word_gamma):
            if _word_gamma: sage_free(_word_gamma)
            if _col_gamma: sage_free(_col_gamma)
            raise MemoryError("Memory.")
        for i from 0 <= i < self.nwords:
            _word_gamma[i] = word_gamma[i]
        for i from 0 <= i < self.ncols:
            _col_gamma[i] = col_gamma[i]
        result = self.is_automorphism(_col_gamma, _word_gamma)
        sage_free(_col_gamma)
        sage_free(_word_gamma)
        return result

    cdef int is_automorphism(self, int *col_gamma, int *word_gamma):
        # TODO: optimize? check only basis?
        cdef int i, j, self_nwords = self.nwords, self_ncols = self.ncols
        for i from 0 <= i < self_nwords:
            for j from 0 <= j < self_ncols:
                if self.is_one(i, j) != self.is_one(word_gamma[i], col_gamma[j]):
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
        cdef int col
        cdef int nwords, word
        nwords = (1 << nrows)
        self.nwords = nwords
        self.ncols = ncols
        self.wd_parent =        <int *> sage_malloc( nwords * sizeof(int) )
        self.wd_rank =          <int *> sage_malloc( nwords * sizeof(int) )
        self.wd_min_cell_rep =  <int *> sage_malloc( nwords * sizeof(int) )
        self.wd_size =          <int *> sage_malloc( nwords * sizeof(int) )
        self.col_parent =       <int *> sage_malloc( ncols * sizeof(int) )
        self.col_rank =         <int *> sage_malloc( ncols * sizeof(int) )
        self.col_min_cell_rep = <int *> sage_malloc( ncols * sizeof(int) )
        self.col_size =         <int *> sage_malloc( ncols * sizeof(int) )
        if not (self.wd_parent and self.wd_rank and self.wd_min_cell_rep and self.wd_size and self.col_parent and self.col_rank and self.col_min_cell_rep and self.col_size):
            if self.wd_parent: sage_free(self.wd_parent)
            if self.wd_rank: sage_free(self.wd_rank)
            if self.wd_min_cell_rep: sage_free(self.wd_min_cell_rep)
            if self.wd_size: sage_free(self.wd_size)
            if self.col_parent:       sage_free(self.col_parent)
            if self.col_rank:         sage_free(self.col_rank)
            if self.col_min_cell_rep: sage_free(self.col_min_cell_rep)
            if self.col_size:         sage_free(self.col_size)
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
        cdef int i
        cdef int j
        s = 'OrbitPartition on %d words and %d columns. Data:\n'%(self.nwords, self.ncols)
#        s += 'Parents::\n'
        s += 'Words:\n'
        for i from 0 <= i < self.nwords:
            s += '%d,'%self.wd_parent[i]
        s = s[:-1] + '\nColumns:\n'
        for j from 0 <= j < self.ncols:
            s += '%d,'%self.col_parent[j]
#        s = s[:-1] + '\n'
#        s += 'Min Cell Reps::\n'
#        s += 'Words:\n'
#        for i from 0 <= i < self.nwords:
#            s += '%d,'%self.wd_min_cell_rep[i]
#        s = s[:-1] + '\nColumns:\n'
#        for j from 0 <= j < self.ncols:
#            s += '%d,'%self.col_min_cell_rep[j]
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
            12

        """
        return self.wd_find(word)

    cdef int wd_find(self, int word):
#        print 'wd_find', word
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
            1

        """
        self.wd_union(x, y)

    cdef void wd_union(self, int x, int y):
#        print 'wd_union', x, y
        cdef int x_root, y_root
        x_root = self.wd_find(x)
        y_root = self.wd_find(y)
        if self.wd_rank[x_root] > self.wd_rank[y_root]:
            self.wd_parent[y_root] = x_root
            self.wd_min_cell_rep[x_root] = min(self.wd_min_cell_rep[x_root],self.wd_min_cell_rep[y_root])
            self.wd_size[x_root] += self.wd_size[y_root]
        elif self.wd_rank[x_root] < self.wd_rank[y_root]:
            self.wd_parent[x_root] = y_root
            self.wd_min_cell_rep[y_root] = min(self.wd_min_cell_rep[x_root],self.wd_min_cell_rep[y_root])
            self.wd_size[y_root] += self.wd_size[x_root]
        elif x_root != y_root:
            self.wd_parent[y_root] = x_root
            self.wd_min_cell_rep[x_root] = min(self.wd_min_cell_rep[x_root],self.wd_min_cell_rep[y_root])
            self.wd_size[x_root] += self.wd_size[y_root]
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
#        print 'col_find', col
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
#        print 'col_union', x, y
        cdef int x_root, y_root
        x_root = self.col_find(x)
        y_root = self.col_find(y)
        if self.col_rank[x_root] > self.col_rank[y_root]:
            self.col_parent[y_root] = x_root
            self.col_min_cell_rep[x_root] = min(self.col_min_cell_rep[x_root],self.col_min_cell_rep[y_root])
            self.col_size[x_root] += self.col_size[y_root]
        elif self.col_rank[x_root] < self.col_rank[y_root]:
            self.col_parent[x_root] = y_root
            self.col_min_cell_rep[y_root] = min(self.col_min_cell_rep[x_root],self.col_min_cell_rep[y_root])
            self.col_size[y_root] += self.col_size[x_root]
        elif x_root != y_root:
            self.col_parent[y_root] = x_root
            self.col_min_cell_rep[x_root] = min(self.col_min_cell_rep[x_root],self.col_min_cell_rep[y_root])
            self.col_size[x_root] += self.col_size[y_root]
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
        cdef int i
        cdef int *_col_gamma
        cdef int *_wd_gamma
        _wd_gamma = <int *> sage_malloc(self.nwords * sizeof(int))
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

    cdef int merge_perm(self, int *col_gamma, int *wd_gamma):
        cdef int i, gamma_i_root
        cdef int j, gamma_j_root, return_value = 0
        cdef int *self_wd_parent = self.wd_parent
        cdef int *self_col_parent = self.col_parent
#        print 'merge_perm'
#        print 'col_gamma:', [col_gamma[i] for i from 0 <= i < self.ncols]
#        print 'wd_gamma:', [wd_gamma[i] for i from 0 <= i < self.nwords]
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
        cdef int k, nwords, ncols, sizeof_int
        cdef PartitionStack other
        cdef int *wd_ents, *wd_lvls, *col_ents, *col_lvls
        cdef int *col_degs, *col_counts, *col_output, *wd_degs, *wd_counts, *wd_output
        sizeof_int = sizeof(int)

        try:
            self.nrows = int(arg1)
            self.nwords = 1 << self.nrows
            self.ncols = int(arg2)
        except:
            other = arg1
            self.nrows = other.nrows
            self.nwords = other.nwords
            self.ncols = other.ncols

        self.radix = sizeof_int << 3
        self.flag = (1 << (self.radix-1))

        # data
        self.wd_ents =    <int *> sage_malloc( self.nwords * sizeof_int )
        self.wd_lvls =    <int *> sage_malloc( self.nwords * sizeof_int )
        self.col_ents =   <int *> sage_malloc( self.ncols  * sizeof_int )
        self.col_lvls =   <int *> sage_malloc( self.ncols  * sizeof_int )

        # scratch space
        self.col_degs =   <int *> sage_malloc( self.ncols  * sizeof_int )
        self.col_counts = <int *> sage_malloc( self.nwords * sizeof_int )
        self.col_output = <int *> sage_malloc( self.ncols  * sizeof_int )
        self.wd_degs =    <int *> sage_malloc( self.nwords * sizeof_int )
        self.wd_counts =  <int *> sage_malloc( (self.ncols+1)  * sizeof_int )
        self.wd_output =  <int *> sage_malloc( self.nwords * sizeof_int )

        if not (self.wd_ents  and self.wd_lvls    and self.col_ents   and self.col_lvls  \
            and self.col_degs and self.col_counts and self.col_output \
            and self.wd_degs  and self.wd_counts  and self.wd_output):
            if self.wd_ents:         sage_free(self.wd_ents)
            if self.wd_lvls:         sage_free(self.wd_lvls)
            if self.col_ents:        sage_free(self.col_ents)
            if self.col_lvls:        sage_free(self.col_lvls)
            if self.col_degs:        sage_free(self.col_degs)
            if self.col_counts:      sage_free(self.col_counts)
            if self.col_output:      sage_free(self.col_output)
            if self.wd_degs:         sage_free(self.wd_degs)
            if self.wd_counts:       sage_free(self.wd_counts)
            if self.wd_output:       sage_free(self.wd_output)
            raise MemoryError("Memory.")

        nwords = self.nwords
        ncols = self.ncols

        if other:
            memcpy(self.wd_ents,  other.wd_ents, self.nwords * sizeof_int)
            memcpy(self.wd_lvls,  other.wd_lvls, self.nwords * sizeof_int)
            memcpy(self.col_ents, other.col_ents, self.ncols * sizeof_int)
            memcpy(self.col_lvls, other.col_lvls, self.ncols * sizeof_int)
        else:
            wd_ents = self.wd_ents
            wd_lvls = self.wd_lvls
            col_ents = self.col_ents
            col_lvls = self.col_lvls
            for k from 0 <= k < nwords-1:
                wd_ents[k] = k
                wd_lvls[k] = 2*ncols
            for k from 0 <= k < ncols-1:
                col_ents[k] = k
                col_lvls[k] = 2*ncols
            wd_ents[nwords-1] = nwords-1
            wd_lvls[nwords-1] = -1
            col_ents[ncols-1] = ncols-1
            col_lvls[ncols-1] = -1

        col_degs = self.col_degs
        col_counts = self.col_counts
        col_output = self.col_output
        wd_degs = self.wd_degs
        wd_counts = self.wd_counts
        wd_output = self.wd_output
        for k from 0 <= k < ncols:
            col_degs[k]=0
            col_output[k]=0
            wd_counts[k]=0
        wd_counts[ncols]=0
        for k from 0 <= k < nwords:
            col_counts[k]=0
            wd_degs[k]=0
            wd_output[k]=0

    def __dealloc__(self):
        if self.basis_locations: sage_free(self.basis_locations)
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

    def print_data(self):
        """
        Prints all data for self.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: print P.print_data()
            nwords:4
            nrows:2
            ncols:6
            radix:32
            wd_ents:
            0
            1
            2
            3
            wd_lvls:
            12
            12
            12
            -1
            col_ents:
            0
            1
            2
            3
            4
            5
            col_lvls:
            12
            12
            12
            12
            12
            -1
            col_degs:
            0
            0
            0
            0
            0
            0
            col_counts:
            0
            0
            0
            0
            col_output:
            0
            0
            0
            0
            0
            0
            wd_degs:
            0
            0
            0
            0
            wd_counts:
            0
            0
            0
            0
            0
            0
            0
            wd_output:
            0
            0
            0
            0

        """
        cdef int i, j
        s = ''
        s += "nwords:" + str(self.nwords) + '\n'
        s += "nrows:" + str(self.nrows) + '\n'
        s += "ncols:" + str(self.ncols) + '\n'
        s += "radix:" + str(self.radix) + '\n'
        s += "wd_ents:" + '\n'
        for i from 0 <= i < self.nwords:
            s += str(self.wd_ents[i]) + '\n'
        s += "wd_lvls:" + '\n'
        for i from 0 <= i < self.nwords:
            s += str(self.wd_lvls[i]) + '\n'
        s += "col_ents:" + '\n'
        for i from 0 <= i < self.ncols:
            s += str(self.col_ents[i]) + '\n'
        s += "col_lvls:" + '\n'
        for i from 0 <= i < self.ncols:
            s += str(self.col_lvls[i]) + '\n'
        s += "col_degs:" + '\n'
        for i from 0 <= i < self.ncols:
            s += str(self.col_degs[i]) + '\n'
        s += "col_counts:" + '\n'
        for i from 0 <= i < self.nwords:
            s += str(self.col_counts[i]) + '\n'
        s += "col_output:" + '\n'
        for i from 0 <= i < self.ncols:
            s += str(self.col_output[i]) + '\n'
        s += "wd_degs:" + '\n'
        for i from 0 <= i < self.nwords:
            s += str(self.wd_degs[i]) + '\n'
        s += "wd_counts:" + '\n'
        for i from 0 <= i < self.ncols + 1:
            s += str(self.wd_counts[i]) + '\n'
        s += "wd_output:" + '\n'
        for i from 0 <= i < self.nwords:
            s += str(self.wd_output[i]) + '\n'
        if self.basis_locations:
            s += "basis_locations:" + '\n'
            j = 1
            while (1 << j) < self.nwords:
                j += 1
            for i from 0 <= i < j:
                s += str(self.basis_locations[i]) + '\n'
        return s

    def __repr__(self):
        """
        Return a string representation of self.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})

        """
        cdef int i, j, k
        s = ''
        last = ''
        current = ''
        for k from 0 <= k < 2*self.ncols:
            current = self._repr_at_k(k)
            if current == last: break
            s += current
            last = current
        return s

    def _repr_at_k(self, k):
        s = '({'
        for j from 0 <= j < self.nwords:
            s += str(self.wd_ents[j])
            if self.wd_lvls[j] <= k:
                s += '},{'
            else:
                s += ','
        s = s[:-2] + ')  '
        s += '({'
        for j from 0 <= j < self.ncols:
            s += str(self.col_ents[j])
            if self.col_lvls[j] <= k:
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
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P._sort_wds(0, [0,2,3,1], 5)
            0
            sage: P
            ({0,3,1,2})  ({0,1,2,3,4,5})
            ({0,3,1,2})  ({0},{1,2,3,4,5})
            ({0,3,1,2})  ({0},{1},{2,3,4,5})
            ({0,3,1,2})  ({0},{1},{2},{3,4,5})
            ({0,3,1,2})  ({0},{1},{2},{3},{4,5})
            ({0},{3},{1},{2})  ({0},{1},{2},{3},{4},{5})
            sage: P._is_discrete(4)
            0
            sage: P._is_discrete(5)
            1

        """
        return self.is_discrete(k)

    cdef int is_discrete(self, int k):
        cdef int i, self_ncols = self.ncols, self_nwords = self.nwords
        cdef int *self_col_lvls = self.col_lvls
        cdef int *self_wd_lvls = self.wd_lvls
        for i from 0 <= i < self_ncols:
            if self_col_lvls[i] > k:
                return 0
        for i from 0 <= i < self_nwords:
            if self_wd_lvls[i] > k:
                return 0
        return 1

    def _num_cells(self, k):
        """
        Returns the number of cells in the partition at level k.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})
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
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P._sat_225(3)
            0
            sage: P._sat_225(4)
            1
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})

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

#    def _new_min_cell_reps(self, k): #TODO
#        """
#        Returns an integer whose bits represent which columns are minimal cell
#        representatives.
#
#        EXAMPLE:
#            sage: import sage.coding.binary_code
#            sage: from sage.coding.binary_code import *
#            sage: P = PartitionStack(2, 6)
#            sage: [P._split_column(i,i+1) for i in range(5)]
#            [0, 1, 2, 3, 4]
#            sage: a = P._min_cell_reps(2)
#            sage: Integer(a).binary()
#            '111'
#            sage: P
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3,4,5})
#            ({0},{1,2,3,4,5})
#            ({0},{1},{2,3,4,5})
#            ({0},{1},{2},{3,4,5})
#            ({0},{1},{2},{3},{4,5})
#            ({0},{1},{2},{3},{4},{5})
#
#        """
#        return self.min_cell_reps(k)
#
#    cdef int min_cell_reps(self, int k):
#        cdef int i
#        cdef int reps = 1
#        cdef int *self_col_lvls = self.col_lvls
#        for i from 0 < i < self.ncols:
#            if self_col_lvls[i-1] <= k:
#                reps += (1 << i)
#        return reps
#
    cdef void new_min_cell_reps(self, int k, int *Omega, int start):
        cdef int i, j
        cdef int *self_col_lvls = self.col_lvls, *self_wd_lvls = self.wd_lvls
        cdef int *self_col_ents = self.col_ents, *self_wd_ents = self.wd_ents
        cdef int reps = (1 << self_col_ents[0])
        cdef int radix = self.radix, nwords = self.nwords, ncols = self.ncols
        for i from 0 < i < ncols:
            if self_col_lvls[i-1] <= k:
                reps += (1 << self_col_ents[i])
        Omega[start] = reps
        reps = 1
        for i from 0 < i < min(radix, nwords):
            if self_wd_lvls[i-1] <= k:
                reps += (1 << self_wd_ents[i])
        Omega[start+1] = reps
        j = radix
        while j < nwords:
            reps = 0
            for i from 0 <= i < min(radix, nwords - j):
                if self_wd_lvls[j + i - 1] <= k:
                    reps += (1 << self_wd_ents[i])
            Omega[start+1+j] = reps
            j += radix

#    def _fixed_cols(self, mcrs, k): #TODO
#        """
#        Returns an integer whose bits represent which columns are fixed. For
#        efficiency, mcrs is the output of min_cell_reps.
#
#        EXAMPLE:
#            sage: import sage.coding.binary_code
#            sage: from sage.coding.binary_code import *
#            sage: P = PartitionStack(2, 6)
#            sage: [P._split_column(i,i+1) for i in range(5)]
#            [0, 1, 2, 3, 4]
#            sage: a = P._fixed_cols(7, 2)
#            sage: Integer(a).binary()
#            '11'
#            sage: P
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3,4,5})
#            ({0},{1,2,3,4,5})
#            ({0},{1},{2,3,4,5})
#            ({0},{1},{2},{3,4,5})
#            ({0},{1},{2},{3},{4,5})
#            ({0},{1},{2},{3},{4},{5})
#
#        """
#        return self.fixed_cols(mcrs, k)
#
#    cdef int fixed_cols(self, int mcrs, int k):
#        cdef int i
#        cdef int fixed = 0
#        cdef int *self_col_lvls = self.col_lvls
#        for i from 0 <= i < self.ncols:
#            if self_col_lvls[i] <= k:
#                fixed += (1 << i)
#        return fixed & mcrs
#
    cdef void fixed_vertices(self, int k, int *Phi, int *Omega, int start):
        cdef int i, j, ell
        cdef int fixed = 0, ncols = self.ncols, nwords = self.nwords
        cdef int *self_col_lvls = self.col_lvls, *self_wd_lvls = self.wd_lvls
        cdef int *self_col_ents = self.col_ents, *self_wd_ents = self.wd_ents
        for i from 0 <= i < ncols:
            if self_col_lvls[i] <= k:
                fixed += (1 << self_col_ents[i])
        Phi[start] = fixed & Omega[start]
        # zero out the rest of Phi
        ell = 1 + nwords/self.radix
        if nwords%self.radix:
            ell += 1
        for i from 0 < i < ell:
            Phi[start+i] = 0
        for i from 0 <= i < nwords:
            if self_wd_lvls[i] <= k:
                ell = self_wd_ents[i]
                j =   ell/self.radix
                ell = ell%self.radix
                if Omega[start+1+j]&(1 << ell):
                    Phi[start+1+j] ^= (1 << ell)

#    def _first_smallest_nontrivial(self, k): #TODO
#        """
#        Returns an integer representing the first, smallest nontrivial cell of columns.
#
#        EXAMPLE:
#            sage: import sage.coding.binary_code
#            sage: from sage.coding.binary_code import *
#            sage: P = PartitionStack(2, 6)
#            sage: [P._split_column(i,i+1) for i in range(5)]
#            [0, 1, 2, 3, 4]
#            sage: a = P._first_smallest_nontrivial(2)
#            sage: Integer(a).binary().zfill(32)
#            '00000000000000000000000000111100'
#            sage: P
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3,4,5})
#            ({0},{1,2,3,4,5})
#            ({0},{1},{2,3,4,5})
#            ({0},{1},{2},{3,4,5})
#            ({0},{1},{2},{3},{4,5})
#            ({0},{1},{2},{3},{4},{5})
#
#        """
#        return self.first_smallest_nontrivial(k)
#
#    cdef int first_smallest_nontrivial(self, int k):
#        cdef int cell
#        cdef int i = 0, j = 0, location = 0, ncols = self.ncols
#        cdef int *self_col_lvls = self.col_lvls
#        while True:
#            if self_col_lvls[i] <= k:
#                if i != j and ncols > i - j + 1:
#                    ncols = i - j + 1
#                    location = j
#                j = i + 1
#            if self_col_lvls[i] == -1: break
#            i += 1
#        # location now points to the beginning of the first, smallest,
#        # nontrivial cell
#        j = location
#        self.v = self.col_ents[j]
#        while True:
#            if self_col_lvls[j] <= k: break
#            j += 1
#        # j now points to the last element of the cell
##        print "fsnt:", location, j-location+1
#        i = self.radix - j - 1                 # the cell is represented in binary, reading from the right:
#        cell = (~0 << location) ^ (~0 << j+1)  # <-------            self.radix               ----->
#        return cell                            # [0]*(radix-j-1) + [1]*(j-location+1) + [0]*location
#
    cdef int new_first_smallest_nontrivial(self, int k, int *W, int start):
        cdef int ell
        cdef int i = 0, j = 0, location = 0, min = self.ncols, nwords = self.nwords
        cdef int min_is_col = 1, radix = self.radix
        cdef int *self_col_lvls = self.col_lvls, *self_wd_lvls = self.wd_lvls
        cdef int *self_col_ents = self.col_ents, *self_wd_ents = self.wd_ents
        while True:
            if self_col_lvls[i] <= k:
                if i != j and min > i - j + 1:
                    min = i - j + 1
                    location = j
                j = i + 1
            if self_col_lvls[i] == -1: break
            i += 1
        i = 0; j = 0
        while True:
            if self_wd_lvls[i] <= k:
                if i != j and min > i - j + 1:
                    min = i - j + 1
                    min_is_col = 0
                    location = j
                j = i + 1
            if self_wd_lvls[i] == -1: break
            i += 1
        # location now points to the beginning of the first, smallest,
        # nontrivial cell
        j = location
        #zero out this level of W:
        ell = 1 + nwords/radix
        if nwords%radix:
            ell += 1
        for i from 0 <= i < ell:
            W[start+i] = 0
        if min_is_col:
            while True:
                if self_col_lvls[j] <= k: break
                j += 1
            # j now points to the last element of the cell
            i = location
            while i <= j:
                W[start] ^= (1 << self_col_ents[i])
                i += 1
            return self_col_ents[location]
        else:
            while True:
                if self_wd_lvls[j] <= k: break
                j += 1
            # j now points to the last element of the cell
            i = location
            while i <= j:
                ell = self_wd_ents[i]
                W[start+1+ell/radix] ^= (1 << ell%radix)
                i += 1
            return self_wd_ents[location] ^ self.flag

    def _dangerous_dont_use_set_ents_lvls(self, col_ents, col_lvls, wd_ents, wd_lvls):
        """
        For debugging only.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            sage: P._dangerous_dont_use_set_ents_lvls([99]*6, [0,3,2,3,5,-1], [4,3,5,6], [3,2,1,-1])
            sage: P
            ({4,3,5,6})  ({99},{99,99,99,99,99})
            ({4,3,5},{6})  ({99},{99,99,99,99,99})
            ({4,3},{5},{6})  ({99},{99,99},{99,99,99})
            ({4},{3},{5},{6})  ({99},{99},{99},{99},{99,99})

        """
        cdef int i
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
            sage: P = PartitionStack(2, 6)
            sage: P._dangerous_dont_use_set_ents_lvls(range(5,-1,-1), [1,2,2,3,3,-1], range(3,-1,-1), [1,1,2,-1])
            sage: P
            ({3,2,1,0})  ({5,4,3,2,1,0})
            ({3},{2},{1,0})  ({5},{4,3,2,1,0})
            ({3},{2},{1},{0})  ({5},{4},{3},{2,1,0})
            ({3},{2},{1},{0})  ({5},{4},{3},{2},{1},{0})
            sage: P._wd_percolate(0,3)
            sage: P._col_percolate(0,5)
            sage: P
            ({0,3,2,1})  ({0,5,4,3,2,1})
            ({0},{3},{2,1})  ({0},{5,4,3,2,1})
            ({0},{3},{2},{1})  ({0},{5},{4},{3,2,1})
            ({0},{3},{2},{1})  ({0},{5},{4},{3},{2},{1})

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
            sage: P = PartitionStack(2, 6)
            sage: P._dangerous_dont_use_set_ents_lvls(range(5,-1,-1), [1,2,2,3,3,-1], range(3,-1,-1), [1,1,2,-1])
            sage: P
            ({3,2,1,0})  ({5,4,3,2,1,0})
            ({3},{2},{1,0})  ({5},{4,3,2,1,0})
            ({3},{2},{1},{0})  ({5},{4},{3},{2,1,0})
            ({3},{2},{1},{0})  ({5},{4},{3},{2},{1},{0})
            sage: P._wd_percolate(0,3)
            sage: P._col_percolate(0,5)
            sage: P
            ({0,3,2,1})  ({0,5,4,3,2,1})
            ({0},{3},{2,1})  ({0},{5,4,3,2,1})
            ({0},{3},{2},{1})  ({0},{5},{4},{3,2,1})
            ({0},{3},{2},{1})  ({0},{5},{4},{3},{2},{1})

        """
        self.wd_percolate(start, end)

    cdef void wd_percolate(self, int start, int end):
        cdef int i, temp
        cdef int *self_wd_ents = self.wd_ents
        for i from end >= i > start:
            if self_wd_ents[i] < self_wd_ents[i-1]:
                temp = self.wd_ents[i]
                self_wd_ents[i] = self_wd_ents[i-1]
                self_wd_ents[i-1] = temp

#    def _split_column(self, int v, int k): #TODO
#        """
#        Split column v out, placing it before the rest of the cell it was in.
#        Returns the location of the split column.
#
#        EXAMPLE:
#            sage: import sage.coding.binary_code
#            sage: from sage.coding.binary_code import *
#            sage: P = PartitionStack(2, 6)
#            sage: [P._split_column(i,i+1) for i in range(5)]
#            [0, 1, 2, 3, 4]
#            sage: P
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3,4,5})
#            ({0},{1,2,3,4,5})
#            ({0},{1},{2,3,4,5})
#            ({0},{1},{2},{3,4,5})
#            ({0},{1},{2},{3},{4,5})
#            ({0},{1},{2},{3},{4},{5})
#            sage: P = PartitionStack(2, 6)
#            sage: P._split_column(0,1)
#            0
#            sage: P._split_column(2,2)
#            1
#            sage: P
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,1,2,3})
#            ({0,2,1,3,4,5})
#            ({0},{2,1,3,4,5})
#            ({0},{2},{1,3,4,5})
#            ({0},{2},{1,3,4,5})
#            ({0},{2},{1,3,4,5})
#            ({0},{2},{1,3,4,5})
#
#        """
#        return self.split_column(v, k)
#
#    cdef int split_column(self, int v, int k):
#        cdef int i = 0, j
#        cdef int *self_col_ents = self.col_ents
#        cdef int *self_col_lvls = self.col_lvls
#        while self_col_ents[i] != v: i += 1
#        j = i
#        while self_col_lvls[i] > k: i += 1
#        if j == 0 or self_col_lvls[j-1] <= k:
#            self.col_percolate(j+1, i)
#        else:
#            while j != 0 and self_col_lvls[j-1] > k:
#                self_col_ents[j] = self_col_ents[j-1]
#                j -= 1
#            self_col_ents[j] = v
#        self_col_lvls[j] = k
#        return j
#

    def _split_vertex(self, v, k):
        """
        Split vertex v out, placing it before the rest of the cell it was in.
        Returns the location of the split vertex.

        NOTE:
            There is a convention regarding whether a vertex is a word or a
            column. See the 'flag' attribute of the PartitionStack object:
            If vertex&flag is not zero, it is a word.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})

        """
        return self.split_vertex(v, k)

    cdef int split_vertex(self, int v, int k):
        cdef int i = 0, j, flag = self.flag
        cdef int *ents
        cdef int *lvls
        if v & flag:
            ents = self.wd_ents
            lvls = self.wd_lvls
            v = v ^ flag
            while ents[i] != v: i += 1
            v = v ^ flag
        else:
            ents = self.col_ents
            lvls = self.col_lvls
            while ents[i] != v: i += 1
        j = i
        while lvls[i] > k: i += 1
        if j == 0 or lvls[j-1] <= k:
            if v & self.flag:
                self.wd_percolate(j+1, i)
            else:
                self.col_percolate(j+1, i)
        else:
            while j != 0 and lvls[j-1] > k:
                ents[j] = ents[j-1]
                j -= 1
            if v & flag:
                ents[j] = v ^ flag
            else:
                ents[j] = v
        lvls[j] = k
        if v & flag:
            return j ^ flag
        else:
            return j

    def _col_degree(self, C, col, wd_ptr, k):
        """
        Returns the number of words in the cell specified by wd_ptr that have a
        1 in the col-th column.

        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0],[0,0,1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [6,2] linear code, generator matrix
            [111100]
            [001111]
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})
            sage: P._col_degree(B, 2, 0, 2)
            2

        """
        return self.col_degree(C, col, wd_ptr, k)

    cdef int col_degree(self, BinaryCode CG, int col, int wd_ptr, int k):
        cdef int i = 0
        cdef int *self_wd_lvls = self.wd_lvls, *self_wd_ents = self.wd_ents
        while True:
            if CG.is_one(self_wd_ents[wd_ptr], col): i += 1
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
            sage: P = PartitionStack(2, 6)
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0],[0,0,1,1,1,1]])
            sage: B = BinaryCode(M)
            sage: B
            Binary [6,2] linear code, generator matrix
            [111100]
            [001111]
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})
            sage: P._wd_degree(B, 1, 1, 1)
            3

        """
        cdef int *ham_wts = hamming_weights()
        result = self.wd_degree(C, wd, col_ptr, k, ham_wts)
        sage_free(ham_wts)
        return result

    cdef int wd_degree(self, BinaryCode CG, int wd, int col_ptr, int k, int *ham_wts):

        cdef int *self_col_lvls = self.col_lvls, *self_col_ents = self.col_ents
        cdef int mask = (1 << self_col_ents[col_ptr])
        while self_col_lvls[col_ptr] > k:
            col_ptr += 1
            mask += (1 << self_col_ents[col_ptr])
        mask &= CG.words[wd]
        return ham_wts[mask & 65535] + ham_wts[(mask >> 16) & 65535]

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
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P._sort_cols(1, [0,2,2,1,1], 1)
            2
            sage: P
            ({0,1,2,3})  ({0,1,4,5,2,3})
            ({0,1,2,3})  ({0},{1},{4,5},{2,3})

        """
        cdef int i
        for i from 0 <= i < len(degrees):
            self.col_degs[i] = degrees[i]
        return self.sort_cols(start, k)

    cdef int sort_cols(self, int start, int k):
        cdef int i, j, max, max_location, self_ncols = self.ncols
        cdef int self_nwords = self.nwords, ii
        cdef int *self_col_counts = self.col_counts
        cdef int *self_col_lvls = self.col_lvls
        cdef int *self_col_degs = self.col_degs
        cdef int *self_col_ents = self.col_ents
        cdef int *self_col_output = self.col_output
        for ii from 0 <= ii < self_nwords:
            self_col_counts[ii] = 0
        i = 0
        while self_col_lvls[i+start] > k:
            self_col_counts[self_col_degs[i]] += 1
            i += 1
        self_col_counts[self_col_degs[i]] += 1

        # i+start is the right endpoint of the cell now
        max = self_col_counts[0]
        max_location = 0
        for ii from 0 < ii < self_nwords:
            if self_col_counts[ii] > max:
                max = self_col_counts[ii]
                max_location = ii
            self_col_counts[ii] += self_col_counts[ii-1]

        for j from i >= j >= 0:
            self_col_counts[self_col_degs[j]] -= 1
            self_col_output[self_col_counts[self_col_degs[j]]] = self_col_ents[start+j]

        max_location = self_col_counts[max_location] + start

        for j from 0 <= j <= i:
            self_col_ents[start+j] = self_col_output[j]

        ii = 1
        while ii < self_nwords and self_col_counts[ii] <= i:
            if self_col_counts[ii] > 0:
                self_col_lvls[start + self_col_counts[ii] - 1] = k
            self.col_percolate(start + self_col_counts[ii-1], start + self_col_counts[ii] - 1)
            ii += 1

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
            sage: P = PartitionStack(3, 6)
            sage: P._sort_wds(0, [0,0,3,3,3,3,2,2], 1)
            4
            sage: P
            ({0,1,6,7,2,3,4,5})  ({0,1,2,3,4,5})
            ({0,1},{6,7},{2,3,4,5})  ({0,1,2,3,4,5})

        """
        cdef int i
        for i from 0 <= i < len(degrees):
            self.wd_degs[i] = degrees[i]
        return self.sort_wds(start, k)

    cdef int sort_wds(self, int start, int k):
        cdef int i, j, max, max_location, self_nwords = self.nwords
        cdef int ii, self_ncols = self.ncols
        cdef int *self_wd_counts = self.wd_counts
        cdef int *self_wd_lvls = self.wd_lvls
        cdef int *self_wd_degs = self.wd_degs
        cdef int *self_wd_ents = self.wd_ents
        cdef int *self_wd_output = self.wd_output

        for ii from 0 <= ii < self_ncols+1:
            self_wd_counts[ii] = 0
        i = 0
        while self_wd_lvls[i+start] > k:
            self_wd_counts[self_wd_degs[i]] += 1
            i += 1
        self_wd_counts[self_wd_degs[i]] += 1

        # i+start is the right endpoint of the cell now
        max = self_wd_counts[0]
        max_location = 0
        for ii from 0 < ii < self_ncols+1:
            if self_wd_counts[ii] > max:
                max = self_wd_counts[ii]
                max_location = ii
            self_wd_counts[ii] += self_wd_counts[ii-1]

        for j from i >= j >= 0:
            if j > i: break # cython bug with ints...
            self_wd_counts[self_wd_degs[j]] -= 1
            self_wd_output[self_wd_counts[self_wd_degs[j]]] = self_wd_ents[start+j]

        max_location = self_wd_counts[max_location] + start

        for j from 0 <= j <= i:
            self_wd_ents[start+j] = self_wd_output[j]

        ii = 1
        while ii < self_ncols+1 and self_wd_counts[ii] <= i:
            if self_wd_counts[ii] > 0:
                self_wd_lvls[start + self_wd_counts[ii] - 1] = k
            self.wd_percolate(start + self_wd_counts[ii-1], start + self_wd_counts[ii] - 1)
            ii += 1

        return max_location

    def _refine(self, k, alpha, CG):
        """
        EXAMPLE:

            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: P = PartitionStack(4, 8)
            sage: P._refine(1, [[0,0],[1,0]], B)
            177
            sage: P._split_vertex(0, 2)
            0
            sage: P._refine(2, [[0,0]], B)
            289
            sage: P._split_vertex(1, 3)
            1
            sage: P._refine(3, [[0,1]], B)
            462
            sage: P._split_vertex(2, 4)
            2
            sage: P._refine(4, [[0,2]], B)
            1498
            sage: P._split_vertex(3, 5)
            3
            sage: P._refine(5, [[0,3]], B)
            641
            sage: P._split_vertex(4, 6)
            4
            sage: P._refine(6, [[0,4]], B)
            1218
            sage: P._is_discrete(5)
            0
            sage: P._is_discrete(6)
            1
            sage: P
            ({0,4,6,2,13,9,11,15,10,14,12,8,7,3,1,5})  ({0,1,2,3,4,7,6,5})
            ({0},{4,6,2,13,9,11,15,10,14,12,8,7,3,1},{5})  ({0,1,2,3,4,7,6,5})
            ({0},{4,6,2,13,9,11,15},{10,14,12,8,7,3,1},{5})  ({0},{1,2,3,4,7,6,5})
            ({0},{4,6,2},{13,9,11,15},{10,14,12,8},{7,3,1},{5})  ({0},{1},{2,3,4,7,6,5})
            ({0},{4},{6,2},{13,9},{11,15},{10,14},{12,8},{7,3},{1},{5})  ({0},{1},{2},{3,4,7,6,5})
            ({0},{4},{6,2},{13,9},{11,15},{10,14},{12,8},{7,3},{1},{5})  ({0},{1},{2},{3},{4,7,6,5})
            ({0},{4},{6},{2},{13},{9},{11},{15},{10},{14},{12},{8},{7},{3},{1},{5})  ({0},{1},{2},{3},{4},{7},{6},{5})

        """
        cdef int i, alpha_length = len(alpha)
        cdef int *_alpha = <int *> sage_malloc( (self.nwords + self.ncols) * sizeof(int) )
        cdef int *ham_wts = hamming_weights()
        if not _alpha:
            sage_free(_alpha)
            raise MemoryError("Memory.")
        for i from 0 <= i < alpha_length:
            if alpha[i][0]:
                _alpha[i] = alpha[i][1] ^ self.flag
            else:
                _alpha[i] = alpha[i][1]
        result = self.refine(k, _alpha, alpha_length, CG, ham_wts)
        sage_free(_alpha)
        sage_free(ham_wts)
        return result

    cdef int refine(self, int k, int *alpha, int alpha_length, BinaryCode CG, int *ham_wts):
        cdef int q, r, s, t, flag = self.flag, self_ncols = self.ncols
        cdef int t_w, self_nwords = self.nwords, invariant = 0, i, j, m = 0
        cdef int *self_wd_degs = self.wd_degs, *self_wd_lvls = self.wd_lvls, *self_wd_ents = self.wd_ents
        cdef int *self_col_degs = self.col_degs, *self_col_lvls = self.col_lvls, *self_col_ents = self.col_ents
        while not self.is_discrete(k) and m < alpha_length:
#            print "m:", m
#            print "alpha:", ','.join(['w'+str(alpha[i]^flag) if alpha[i]&flag else 'c'+str(alpha[i]) for i from 0 <= i < alpha_length])
#            print self
            invariant += 1
            j = 0
            if alpha[m] & flag:
#                print 'word'
                while j < self_ncols:
#                    print 'j', j
#                    print self
                    i = j; s = 0
                    invariant += 8
                    while True:
#                        print 'col_i', self_col_ents[i]
#                        print 'alpha[m]^flag', alpha[m]^flag
                        self_col_degs[i-j] = self.col_degree(CG, self_col_ents[i], alpha[m]^flag, k)
#                        print 'deg', self_col_degs[i-j]
                        if s == 0 and self_col_degs[i-j] != self_col_degs[0]: s = 1
                        i += 1
                        if self_col_lvls[i-1] <= k: break
                    if s:
#                        print 's'
                        invariant += 8
                        t = self.sort_cols(j, k)
                        invariant += t + self_col_degs[i-j-1]
                        q = m
                        while q < alpha_length:
                            if alpha[q] == j:
                                alpha[q] = t
                                break
                            q += 1
                        r = j
                        while True:
                            if r == j or self.col_lvls[r-1] == k:
                                if r != t:
                                    alpha[alpha_length] = r
                                    alpha_length += 1
                            r += 1
                            if r >= i: break
                        invariant += (i-j)
                    j = i
            else:
#                print 'col'
                while j < self.nwords:
#                    print 'j', j
#                    print self
                    i = j; s = 0
                    invariant += 64
                    while True:
#                        print 'i', i
                        self_wd_degs[i-j] = self.wd_degree(CG, self_wd_ents[i], alpha[m], k, ham_wts)
#                        print 'deg', self_wd_degs[i-j]
                        if s == 0 and self_wd_degs[i-j] != self_wd_degs[0]: s = 1
                        i += 1
                        if self_wd_lvls[i-1] <= k: break
                    if s:
                        invariant += 64
                        t_w = self.sort_wds(j, k)
                        invariant += t_w + self_wd_degs[i-j-1]
                        q = m
                        j ^= flag
                        while q < alpha_length:
                            if alpha[q] == j:
                                alpha[q] = t_w ^ flag
                                break
                            q += 1
                        j ^= flag
                        r = j
                        while True:
                            if r == j or self.wd_lvls[r-1] == k:
                                if r != t_w:
                                    alpha[alpha_length] = r^flag
                                    alpha_length += 1
                            r += 1
                            if r >= i: break
                        invariant += (i-j)
                    j = i
            m += 1
        if invariant != -1:
            return invariant
        else:
            return 0

    def _clear(self, k):
        """
        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(2, 6)
            sage: [P._split_vertex(i,i+1) for i in range(5)]
            [0, 1, 2, 3, 4]
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2,3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3,4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4,5})
            ({0,1,2,3})  ({0},{1},{2},{3},{4},{5})
            sage: P._clear(2)
            sage: P
            ({0,1,2,3})  ({0,1,2,3,4,5})
            ({0,1,2,3})  ({0},{1,2,3,4,5})

        """
        self.clear(k)

    cdef void clear(self, int k):
        cdef int i, j = 0, nwords = self.nwords, ncols = self.ncols
        cdef int *wd_lvls = self.wd_lvls, *col_lvls = self.col_lvls
        for i from 0 <= i < nwords:
            if wd_lvls[i] >= k:
                wd_lvls[i] += 1
            if wd_lvls[i] < k:
                self.wd_percolate(j, i)
                j = i + 1
        j = 0
        for i from 0 <= i < ncols:
            if col_lvls[i] >= k:
                col_lvls[i] += 1
            if col_lvls[i] < k:
                self.col_percolate(j, i)
                j = i + 1

    def _cmp(self, other, C):
        """
        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: P = PartitionStack(4, 8)
            sage: P._refine(0, [[0,0],[1,0]], B)
            177
            sage: P._split_vertex(0, 1)
            0
            sage: P._refine(1, [[0,0]], B)
            289
            sage: P._split_vertex(1, 2)
            1
            sage: P._refine(2, [[0,1]], B)
            462
            sage: P._split_vertex(2, 3)
            2
            sage: P._refine(3, [[0,2]], B)
            1498
            sage: P._split_vertex(4, 4)
            4
            sage: P._refine(4, [[0,4]], B)
            1218
            sage: P._is_discrete(4)
            1
            sage: Q = PartitionStack(P)
            sage: Q._clear(4)
            sage: Q._split_vertex(5, 4)
            4
            sage: Q._refine(4, [[0,4]], B)
            1219
            sage: Q._is_discrete(4)
            1
            sage: Q._cmp(P, B)
            0

        """
        return self.cmp(other, C)

    cdef int cmp(self, PartitionStack other, BinaryCode CG):
        cdef int *self_wd_ents = self.wd_ents
        cdef int *CG_words = CG.words
        cdef int i, j, l, m, span = 1, ncols = self.ncols, nwords = self.nwords
        for i from 0 <= i < nwords: # TODO: probably don't need to check i == 0 here!
            for j from 0 <= j < ncols:
                l = CG.is_one(self.wd_ents[i], self.col_ents[j])
                m = CG.is_one(other.wd_ents[i], other.col_ents[j])
                if l != m:
                    return l - m
        return 0

    def print_basis(self):
        """
        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 8)
            sage: P._dangerous_dont_use_set_ents_lvls(range(8), range(7)+[-1], [4,7,12,11,1,9,3,0,2,5,6,8,10,13,14,15], [0]*16)
            sage: P
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1,2,3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2,3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5},{6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5},{6},{7})
            sage: P._find_basis()
            sage: P.print_basis()
            basis_locations:
            4
            8
            0
            11

        """
        cdef int i, j
        if self.basis_locations:
            print "basis_locations:"
            j = 1
            while (1 << j) < self.nwords:
                j += 1
            for i from 0 <= i < j:
                print self.basis_locations[i]

    def _find_basis(self):
        """
        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: P = PartitionStack(4, 8)
            sage: P._dangerous_dont_use_set_ents_lvls(range(8), range(7)+[-1], [4,7,12,11,1,9,3,0,2,5,6,8,10,13,14,15], [0]*16)
            sage: P
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1,2,3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2,3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3,4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4,5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5,6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5},{6,7})
            ({4},{7},{12},{11},{1},{9},{3},{0},{2},{5},{6},{8},{10},{13},{14},{15})  ({0},{1},{2},{3},{4},{5},{6},{7})
            sage: P._find_basis()
            sage: P.print_basis()
            basis_locations:
            4
            8
            0
            11

        """
        cdef int i
        cdef int *ham_wts = hamming_weights()
        self.find_basis(ham_wts)
        sage_free(ham_wts)

    cdef void find_basis(self, int *ham_wts):
        cdef int i = 0, j, k, nwords = self.nwords, weight, basis_elts = 0, nrows = self.nrows
        cdef int *self_wd_ents = self.wd_ents
        if not self.basis_locations:
            self.basis_locations = <int *> sage_malloc( nrows * sizeof(int) )
        if not self.basis_locations:
            raise MemoryError("Memory.")
        while i < nwords:
            j = self_wd_ents[i]
            weight = ham_wts[j & 65535] + ham_wts[(j>>16) & 65535]
            if weight == 1:
                basis_elts += 1
                k = 0
                while not (1<<k) & j:
                    k += 1
                self.basis_locations[k] = i
                if basis_elts == nrows: break
            i += 1

    def _get_permutation(self, other):
        """
        EXAMPLE:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: M = Matrix(GF(2), [[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[1,0,1,0,1,0,1,0]])
            sage: B = BinaryCode(M)
            sage: P = PartitionStack(4, 8)
            sage: P._refine(0, [[0,0],[1,0]], B)
            177
            sage: P._split_vertex(0, 1)
            0
            sage: P._refine(1, [[0,0]], B)
            289
            sage: P._split_vertex(1, 2)
            1
            sage: P._refine(2, [[0,1]], B)
            462
            sage: P._split_vertex(2, 3)
            2
            sage: P._refine(3, [[0,2]], B)
            1498
            sage: P._split_vertex(4, 4)
            4
            sage: P._refine(4, [[0,4]], B)
            1218
            sage: P._is_discrete(4)
            1
            sage: Q = PartitionStack(P)
            sage: Q._clear(4)
            sage: Q._split_vertex(5, 4)
            4
            sage: Q._refine(4, [[0,4]], B)
            1219
            sage: Q._is_discrete(4)
            1
            sage: P._get_permutation(Q)
            ([0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 8, 9, 10, 11], [0, 1, 2, 3, 5, 4, 7, 6])

        """
        cdef int i
        cdef int *ham_wts = hamming_weights()
        cdef int *word_g = <int *> sage_malloc( self.nwords * sizeof(int) )
        cdef int *col_g = <int *> sage_malloc( self.ncols * sizeof(int) )
        if not (word_g and col_g):
            if word_g: sage_free(word_g)
            if col_g: sage_free(col_g)
            raise MemoryError("Memory.")
        self.get_permutation(other, word_g, col_g, ham_wts)
        sage_free(ham_wts)
        word_l = [word_g[i] for i from 0 <= i < self.nwords]
        col_l = [col_g[i] for i from 0 <= i < self.ncols]
        sage_free(word_g)
        sage_free(col_g)
        return word_l, col_l

    cdef void get_permutation(self, PartitionStack other, int *word_gamma, int *col_gamma, int *ham_wts):
        cdef int i
        cdef int *self_wd_ents = self.wd_ents, *other_wd_ents = other.wd_ents
        cdef int *self_col_ents = self.col_ents, *other_col_ents = other.col_ents
        # word_gamma[i] := image of the ith row as linear comb of rows
        for i from 0 <= i < self.nwords:
            word_gamma[other_wd_ents[i]] = self_wd_ents[i]
        for i from 0 <= i < self.ncols:
            col_gamma[other_col_ents[i]] = self_col_ents[i]

cdef class BinaryCodeClassifier:

    def __new__(self):
        self.radix = sizeof(int) << 3
        self.ham_wts = hamming_weights()
        self.L = 100 # memory limit for Phi and Omega- multiply by 8KB
        self.aut_gens_size = self.radix * 100

        self.w_gamma_size = 1 << (self.radix/2)
        self.alpha_size = self.w_gamma_size + self.radix
        self.Phi_size = self.w_gamma_size/self.radix + 1

        self.w_gamma =     <int *> sage_malloc( self.w_gamma_size              * sizeof(int) )
        self.alpha =       <int *> sage_malloc( self.alpha_size                * sizeof(int) )
        self.Phi =     <int *> sage_malloc( self.Phi_size * (self.L+1)     * sizeof(int) )
        self.Omega =   <int *> sage_malloc( self.Phi_size * self.L         * sizeof(int) )
        self.W =       <int *> sage_malloc( self.Phi_size * self.radix * 2 * sizeof(int) )

        self.aut_gp_gens = <int *> sage_malloc( self.aut_gens_size             * sizeof(int) )
        self.c_gamma =     <int *> sage_malloc( self.radix                     * sizeof(int) )
        self.labeling =    <int *> sage_malloc( self.radix * 2                 * sizeof(int) )
        self.Lambda1 =     <int *> sage_malloc( self.radix * 2                 * sizeof(int) )
        self.Lambda2 =     <int *> sage_malloc( self.radix * 2                 * sizeof(int) )
        self.Lambda3 =     <int *> sage_malloc( self.radix * 2                 * sizeof(int) )
        self.v =           <int *> sage_malloc( self.radix * 2                 * sizeof(int) )
        self.e =           <int *> sage_malloc( self.radix * 2                 * sizeof(int) )

        if not (self.Phi and self.Omega and self.W and self.Lambda1 and self.Lambda2 and self.Lambda3 \
            and self.w_gamma and self.c_gamma and self.alpha and self.v and self.e and self.aut_gp_gens \
            and self.labeling):
            if self.Phi:          sage_free(self.Phi)
            if self.Omega:        sage_free(self.Omega)
            if self.W:            sage_free(self.W)
            if self.Lambda1:      sage_free(self.Lambda1)
            if self.Lambda2:      sage_free(self.Lambda2)
            if self.Lambda3:      sage_free(self.Lambda3)
            if self.w_gamma:      sage_free(self.w_gamma)
            if self.c_gamma:      sage_free(self.c_gamma)
            if self.alpha:        sage_free(self.alpha)
            if self.v:            sage_free(self.v)
            if self.e:            sage_free(self.e)
            if self.aut_gp_gens:  sage_free(self.aut_gp_gens)
            if self.labeling:     sage_free(self.labeling)
            raise MemoryError("Memory.")

    def __dealloc__(self):
        if self.ham_wts: sage_free(self.ham_wts)
        if self.Phi:     sage_free(self.Phi)
        if self.Omega:   sage_free(self.Omega)
        if self.W:       sage_free(self.W)
        if self.Lambda1: sage_free(self.Lambda1)
        if self.Lambda2: sage_free(self.Lambda2)
        if self.Lambda3: sage_free(self.Lambda3)
        if self.c_gamma: sage_free(self.c_gamma)

        if self.w_gamma:   sage_free(self.w_gamma)
        if self.alpha:     sage_free(self.alpha)

        if self.v:           sage_free(self.v)
        if self.e:           sage_free(self.e)
        if self.aut_gp_gens: sage_free(self.aut_gp_gens)
        if self.labeling:    sage_free(self.labeling)

    cdef void record_automorphism(self, int *gamma, int ncols):
        cdef int i, j
        if self.aut_gp_index + ncols > self.aut_gens_size:
            self.aut_gens_size *= 2
            self.aut_gp_gens = <int *> sage_realloc( self.aut_gp_gens, self.aut_gens_size )
            if not self.aut_gp_gens:
                raise MemoryError("Memory.")
        j = self.aut_gp_index
        for i from 0 <= i < ncols:
            self.aut_gp_gens[i+j] = gamma[i]
        self.aut_gp_index += ncols

    def _aut_gp_and_can_label(self, CC, verbosity=0):
        """
        Compute the automorphism group and canonical label of the code CC.

        INPUT:
            CC - a BinaryCode object
            verbosity - a nonnegative integer

        OUTPUT:
            a tuple, (gens, labeling, size)
            gens -- a list of permutations (in list form) representing generators
                of the permutation automorphism group of the code CC.
            labeling -- a permutation representing the canonical labeling of the
                code. mostly for internal use; if the dimension of the code is k
                and the degree (number of columns) is n, then the first n entries
                describe the relabeling on the columns, and the next k describe
                where the basis is sent.
            size -- the order of the automorphism group.

        EXAMPLES:
            sage: import sage.coding.binary_code
            sage: from sage.coding.binary_code import *
            sage: BC = BinaryCodeClassifier()

            sage: M = Matrix(GF(2),[\
            ... [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],\
            ... [0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0],\
            ... [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],\
            ... [0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1],\
            ... [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]])
            sage: B = BinaryCode(M)
            sage: gens, labeling, size = BC._aut_gp_and_can_label(B)
            sage: S = SymmetricGroup(M.ncols())
            sage: L = [S([x+1 for x in g]) for g in gens]
            sage: PermutationGroup(L).order()
            322560
            sage: size
            322560

            sage: M = Matrix(GF(2),[\
            ... [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],\
            ... [0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0],\
            ... [0,0,0,0,0,1,0,1,0,0,0,1,1,1,1,1,1],\
            ... [0,0,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1]])
            sage: B = BinaryCode(M)
            sage: gens, labeling, size = BC._aut_gp_and_can_label(B)
            sage: S = SymmetricGroup(M.ncols())
            sage: L = [S([x+1 for x in g]) for g in gens]
            sage: PermutationGroup(L).order()
            2304
            sage: size
            2304

            sage: M=Matrix(GF(2),[\
            ... [1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0],\
            ... [0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0],\
            ... [0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0],\
            ... [0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0],\
            ... [0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0],\
            ... [0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0],\
            ... [0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0],\
            ... [0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1]])
            sage: B = BinaryCode(M)
            sage: gens, labeling, size = BC._aut_gp_and_can_label(B)
            sage: S = SymmetricGroup(M.ncols())
            sage: L = [S([x+1 for x in g]) for g in gens]
            sage: PermutationGroup(L).order()
            136
            sage: size
            136

            sage: M=Matrix(GF(2),[\
            ... [0,1,0,1,1,1,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1],\
            ... [1,0,1,1,1,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1,0],\
            ... [0,1,1,1,0,0,0,1,0,0,1,1,0,0,0,1,1,1,0,1,0,0],\
            ... [1,1,1,0,0,0,1,0,0,1,0,0,0,0,1,1,1,0,1,0,0,1],\
            ... [1,1,0,0,0,1,0,0,1,0,1,0,0,1,1,1,0,1,0,0,1,0],\
            ... [1,0,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,0,0,1,0,0],\
            ... [0,0,0,1,0,0,1,0,1,1,1,1,1,1,0,1,0,0,1,0,0,0],\
            ... [0,0,1,0,0,1,0,1,1,1,0,1,1,0,1,0,0,1,0,0,0,1],\
            ... [0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,0,0,1,1],\
            ... [1,0,0,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,0,1,1,1],\
            ... [0,0,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,0,1,1,1,0]])
            sage: B = BinaryCode(M)
            sage: gens, labeling, size = BC._aut_gp_and_can_label(B)
            sage: S = SymmetricGroup(M.ncols())
            sage: L = [S([x+1 for x in g]) for g in gens]
            sage: PermutationGroup(L).order()
            887040
            sage: size
            887040

            sage: B = BinaryCode(Matrix(GF(2),[[1,0,1],[0,1,1]]))
            sage: BC._aut_gp_and_can_label(B)
            ([[0, 2, 1], [1, 0, 2]], [0, 1, 2, 2, 1], 6)

            sage: B = BinaryCode(Matrix(GF(2),[[1,1,1,1]]))
            sage: BC._aut_gp_and_can_label(B)
            ([[0, 1, 3, 2], [0, 2, 1, 3], [1, 0, 2, 3]], [0, 1, 2, 3, 1], 24)

        """
        cdef int i, j
        cdef BinaryCode C = CC
        self.aut_gp_and_can_label(C, verbosity)
        i = 0
        py_aut_gp_gens = []
        while i < self.aut_gp_index:
            gen = [self.aut_gp_gens[i+j] for j from 0 <= j < C.ncols]
            py_aut_gp_gens.append(gen)
            i += C.ncols
        py_labeling = [self.labeling[i] for i from 0 <= i < C.ncols] \
                    + [self.labeling[i+C.ncols] for i from 0 <= i < C.nrows]
        aut_gp_size = self.aut_gp_size
        return py_aut_gp_gens, py_labeling, aut_gp_size

    cdef void aut_gp_and_can_label(self, BinaryCode C, int verbosity):

        # declare variables:
        cdef int i, j, ii, jj, iii, jjj, iiii # local variables

        cdef PartitionStack nu, zeta, rho # nu is the current position in the tree,
                                          # zeta the first terminal position,
                                          # and rho the best-so-far guess at canonical labeling position
        cdef int k = 0  # the number of partitions in nu
        cdef int k_rho  # the number of partitions in rho
        cdef int *v = self.v    # list of vertices determining nu
        cdef int h = -1 # longest common ancestor of zeta and nu: zeta[h] == nu[h], zeta[h+1] != nu[h+1]
                        # -1 indicates that zeta is not yet defined
        cdef int hb     # longest common ancestor of rho and nu:
                        # rho[hb] == nu[hb], rho[hb+1] != nu[hb+1]
        cdef int hh = 1 # the height of the oldest ancestor of nu satisfying Lemma 2.25 in [1]:
                        # if nu does not satisfy it at k, then hh = k
        cdef int ht # smallest such that all descendants of zeta[ht] are equivalent under
                    # the portion of the automorphism group so far discovered
        cdef int *alpha # for storing pointers to cells of nu[k]
        cdef int tvc    # tvc keeps track of which vertex is the first where nu and zeta differ-
                        # zeta was defined by splitting one vertex, and nu was defined by splitting tvc

        cdef OrbitPartition Theta # keeps track of which vertices have been discovered to be equivalent
        cdef int *Phi   = self.Phi      # Phi stores the fixed point sets of each automorphism
        cdef int *Omega = self.Omega    # Omega stores the minimal elements of each cell of the orbit partition
        cdef int l = -1                 # current index for storing values in Phi and Omega- we start at -1 so that when
                                        # we increment first, the first place we write to is 0.
        cdef int *W = self.W    # for each k, W[k] is a list (as int mask) of the vertices to be searched down from
                                # the current partition, at k. Phi and Omega are ultimately used to make the size of
                                # W as small as possible
        cdef int *e = self.e    # 0 or 1, whether or not we have used Omega and Phi to narrow down W[k] yet: see states 12 and 17

        cdef int index = 0, size = 1    # Define $\Gamma^{(-1)} := \text{Aut}(C)$, and
                                        # $\Gamma^{(i)} := \Gamma^{(-1)}_{v_0,...,v_i}$.
                                        # Then index = $|\Gamma^{(k-1)}|/|\Gamma^{(k)}|$ at (POINT A)
                                        # and size = $|\Gamma^{(k-1)}|$ at (POINT A) and (POINT B).

        cdef int *Lambda = self.Lambda1             # for tracking indicator values- zf and zb are
        cdef int *zf__Lambda_zeta = self.Lambda2    # indicator vectors remembering Lambda[k] for
        cdef int *zb__Lambda_rho = self.Lambda3     # zeta and rho, respectively
        cdef int qzb              # keeps track of Lambda[k] {>,<,=} zb[k]
        cdef int hzf__h_zeta      # the max height for which Lambda and zf agree
        cdef int hzb__h_rho = -1  # the max height for which Lambda and zb agree

        cdef int *word_gamma, *col_gamma = self.c_gamma # used for storing permutations
        cdef int nwords = C.nwords, ncols = C.ncols, nrows = C.nrows
        cdef int *ham_wts = self.ham_wts
        cdef int state  # keeps track of position in algorithm - see sage/graphs/graph_isom.pyx, search for "STATE DIAGRAM"

        self.aut_gp_index = 0

        if self.w_gamma_size < nwords:
            while self.w_gamma_size < nwords:
                self.w_gamma_size *= 2
            self.alpha_size = self.w_gamma_size + self.radix
            self.Phi_size = self.w_gamma_size/self.radix + 1
            self.w_gamma = <int *> sage_realloc(self.w_gamma,   self.w_gamma_size              * sizeof(int) )
            self.alpha =   <int *> sage_realloc(self.alpha,     self.alpha_size                * sizeof(int) )
            self.Phi =     <int *> sage_realloc(self.new_Phi,   self.Phi_size * self.L         * sizeof(int) )
            self.Omega =   <int *> sage_realloc(self.new_Omega, self.Phi_size * self.L         * sizeof(int) )
            self.W =       <int *> sage_realloc(self.new_W,     self.Phi_size * self.radix * 2 * sizeof(int) )
            if not (self.w_gamma and self.alpha and self.Phi and self.Omega and self.W):
                if self.w_gamma: sage_free(self.w_gamma)
                if self.alpha: sage_free(self.alpha)
                if self.Phi: sage_free(self.Phi)
                if self.Omega: sage_free(self.Omega)
                if self.W: sage_free(self.W)
                raise MemoryError("Memory.")
        word_gamma = self.w_gamma
        alpha = self.alpha # think of alpha as of length exactly nwords + ncols
        nu =    PartitionStack(nrows, ncols)
        Theta = OrbitPartition(nrows, ncols)

        # trivial case
        if ncols == 0 or nrows == 0:
            raise NotImplementedError("Must supply a nontrivial code.")

        state = 1
        while state != -1:
            if False:
                print '-----'
                print "k:", k
                print "h:", h
            if False:
                if k != -1:
                    if v[k]&nu.flag:
                        print "v[k]: word ", v[k]^nu.flag
                    else:
                        print "v[k]: col ", v[k]
                    if tvc&nu.flag:
                        print "tvc- wd", tvc^nu.flag
                    else:
                        print "tvc- col", tvc
                    if W[self.Phi_size * k]:
                        print "W[k]: cols", Integer(W[self.Phi_size * k]).binary()
                    else:
                        j = nwords/self.radix
                        if nwords%self.radix:
                            j += 1
                        L = ''
                        for i from 0 <= i < j:
                            if i == j - 1:
                                jj = nwords%self.radix
                                if jj == 0:
                                    jj = self.radix
                            else:
                                jj = self.radix
                            for ii from 0 <= ii < jj:
                                if W[self.Phi_size * k + 1 + i] & (1 << ii):
                                    L += '1'
                                else:
                                    L += '0'
                        print "W[k]: words", L#[Integer(W[self.Phi_size * k + 1 + i]).binary() for i from 0 <= i < j]
            if False:
                print 'nu'
                print nu
                if tvc&nu.flag:
                    print 'tvc is word', tvc^nu.flag
                else:
                    print 'tvc is col', tvc
                if v[k]&nu.flag:
                    print 'v[k] is word', v[k]^nu.flag
                else:
                    print 'v[k] is col', v[k]
            if False:
                if h != -1:
                    print 'zeta'
                    print zeta
                    print 'rho'
                    print rho
                print "hzf:", hzf__h_zeta
                print "aut_gp_index", self.aut_gp_index
                print 'hh', hh
                print 'ht', ht
                print 'hzf__h_zeta', hzf__h_zeta
                print 'qzb', qzb
            if False:
                print '-----'
                print "state:", state

            if state == 1: # Entry point: once only
                alpha[0] = 0
                alpha[1] = nu.flag
                nu.refine(k, alpha, 2, C, ham_wts)
                if nu.sat_225(k): hh = k
                if nu.is_discrete(k): state = 18; continue

                # store the first smallest nontrivial cell in W[k], and set v[k]
                # equal to its minimum element
                v[k] = nu.new_first_smallest_nontrivial(k, W, self.Phi_size * k)

                Lambda[k] = 0
                e[k] = 0
                state = 2

            elif state == 2: # Move down the search tree one level by refining nu:
                             # split out a vertex, and refine nu against it
                k += 1
                # TODO: Consider removing the following lines
                if k >= 2*self.radix: raise RuntimeError(\
                    "A counterexample to an assumption the author made while writing this software has been encountered.")
                nu.clear(k)

                alpha[0] = nu.split_vertex(v[k-1], k)
                Lambda[k] = nu.refine(k, alpha, 1, C, ham_wts) # store the invariant to Lambda[k]
                # only if this is the first time moving down the search tree:
                if h == -1: state = 5; continue

                # update hzf__h_zeta
                if hzf__h_zeta == k-1 and Lambda[k] == zf__Lambda_zeta[k]: hzf__h_zeta = k
                # update qzb
                if qzb == 0:
                    if zb__Lambda_rho[k] == -1 or Lambda[k] < zb__Lambda_rho[k]:
                        qzb = -1
                    elif Lambda[k] > zb__Lambda_rho[k]:
                        qzb = 1
                    else:
                        qzb = 0
                # update hzb
                if hzb__h_rho == k-1 and qzb == 0: hzb__h_rho = k
                # if Lambda[k] > zb[k], then zb[k] := Lambda[k]
                # (zb keeps track of the indicator invariants corresponding to
                # rho, the closest canonical leaf so far seen- if Lambda is
                # bigger, then rho must be about to change
                if qzb > 0: zb__Lambda_rho[k] = Lambda[k]
                state = 3

            elif state == 3: # attempt to rule out automorphisms while moving down the tree
                # if k > hzf, then we know that nu currently does not look like zeta, the first
                # terminal node encountered, thus there is no automorphism to discover. If qzb < 0,
                # i.e. Lambda[k] < zb[k], then the indicator is not maximal, and we can't reach a
                # canonical leaf. If neither of these is the case, then proceed to state 6.
                if hzf__h_zeta <= k or qzb >= 0: state = 4
                else: state = 6

            elif state == 4: # at this point we have -not- ruled out the presence of automorphisms
                if nu.is_discrete(k): state = 7; continue # we have a terminal node, so process it

                # otherwise, prepare to split out another column:
                # store the first smallest nontrivial cell in W[k], and set v[k]
                # equal to its minimum element
                v[k] = nu.new_first_smallest_nontrivial(k, W, self.Phi_size * k)
                if not nu.sat_225(k): hh = k + 1
                e[k] = 0 # see state 12 and 17
                state = 2 # continue down the tree

            elif state == 5: # same as state 3, but in the case where we haven't yet defined zeta
                             # i.e. this is our first time down the tree. Once we get to the bottom,
                             # we will have zeta = nu = rho, so we do:
                zf__Lambda_zeta[k] = Lambda[k]
                zb__Lambda_rho[k] = Lambda[k]
                state = 4

            elif state == 6: # at this stage, there is no reason to continue downward, so backtrack
                j = k
#                print 'current k', j
#                print 'ht', ht
#                print 'hzb__h_rho', hzb__h_rho
#                print 'hh', hh
                # return to the longest ancestor nu[i] of nu that could have a
                # descendant equivalent to zeta or could improve on rho.
                # All terminal nodes descending from nu[hh] are known to be
                # equivalent, so i < hh. Also, if i > hzb, none of the
                # descendants of nu[i] can improve rho, since the indicator is
                # off (Lambda(nu) < Lambda(rho)). If i >= ht, then no descendant
                # of nu[i] is equivalent to zeta (see [1, p67]).
                if ht-1 > hzb__h_rho:
                    if ht-1 < hh-1:
                        k = ht-1
                    else:
                        k = hh-1
                else:
                    if hzb__h_rho < hh-1:
                        k = hzb__h_rho
                    else:
                        k = hh-1
                # TODO: Consider removing the following lines
                if k >= 2*self.radix: raise RuntimeError(\
                    "A counterexample to an assumption the author made while writing this software has been encountered.")
                # TODO: is the following line necessary?
                if k == -1: k = 0

                if hb > k:# update hb since we are backtracking
                    hb = k
                # if j == hh, then all nodes lower than our current position are equivalent, so bail out
                if j == hh: state = 13; continue

                # recall hh: the height of the oldest ancestor of zeta for which Lemma 2.25 is
                # satsified, which implies that all terminal nodes descended from there are equivalent.
                # If we are looking at such a node, then the partition at nu[hh] can be used for later
                # pruning, so we store its fixed set and a set of representatives of its cells.
                if l < self.L-1: l += 1
                nu.new_min_cell_reps(hh, Omega, self.Phi_size*l)
                nu.fixed_vertices(hh, Phi, Omega, self.Phi_size*l)

                state = 12

            elif state == 7: # we have just arrived at a terminal node of the search tree T(G, Pi)
                # if this is the first terminal node, go directly to 18, to
                # process zeta
                if h == -1: state = 18; continue

                # hzf is the extremal height of ancestors of both nu and zeta, so if k < hzf, nu is not
                # equivalent to zeta, i.e. there is no automorphism to discover.
                if k < hzf__h_zeta: state = 8; continue

                nu.get_permutation(zeta, word_gamma, col_gamma, ham_wts)
#                print "gamma:", str([[word_gamma[i] for i from 0 <= i < nwords], [col_gamma[i] for i from 0 <= i < ncols]]).replace(' ','')
#                print Theta
                # if C^gamma == C, the permutation is an automorphism, goto 10
                if C.is_automorphism(col_gamma, word_gamma):
                    state = 10
                else:
                    state = 8

            elif state == 8: # we have just ruled out the presence of automorphism and have not yet
                             # considered whether nu improves on rho
                # if qzb < 0, then rho already has larger indicator tuple
                if qzb < 0: state = 6; continue

                # if Lambda[k] > zb[k] or nu is shorter than rho, then we have an improvement for rho
                if (qzb > 0) or (k < k_rho): state = 9; continue

                # now Lambda[k] == zb[k] and k == k_rho, so we appeal to an enumeration:
                j = nu.cmp(rho, C)
                # if C(nu) > C(rho), we have a new label, goto 9
                if j > 0: state = 9; continue

                # if C(nu) < C(rho), no new label, goto 6
                if j < 0: state = 6; continue

                # if C(nu) == C(rho), get the automorphism and goto 10
                rho.get_permutation(nu, word_gamma, col_gamma, ham_wts)
#                print "gamma:", str([[word_gamma[i] for i from 0 <= i < nwords], [col_gamma[i] for i from 0 <= i < ncols]]).replace(' ','')
#                print Theta
                state = 10

            elif state == 9: # nu is a better guess at the canonical label than rho
                rho = PartitionStack(nu)
                k_rho = k
                qzb = 0
                hb = k
                hzb__h_rho = k
                # set zb[k+1] = Infinity
                zb__Lambda_rho[k+1] = -1
                state = 6

            elif state == 10: # we have an automorphism to process
                # increment l
                if l < self.L-1: l += 1

                # store information about the automorphism to Omega and Phi
                ii = self.Phi_size*l
                Omega[ii] = ~(~0 << ncols)
                Phi[ii] = 0
                jj = 1 + nwords/self.radix
                if nwords%self.radix:
                    jj += 1
                for i from 0 < i < jj:
                    Omega[ii+i] = ~0
                    Phi[ii+i] = 0
                Omega[ii+jj-1] = ~(~0 << nwords%self.radix)
                # Omega stores the minimum cell representatives
                i = 0
                while i < ncols:
                    j = col_gamma[i]         # i is a minimum
                    while j != i:            # cell rep,
                        Omega[ii] ^= (1<<j)  # so cancel
                        j = col_gamma[j]     # cellmates
                    i += 1
                    while i < ncols and not Omega[ii]&(1<<i): # find minimal element
                        i += 1                                # of next cell
                i = 0
                jj = self.radix
                while i < nwords:
                    j = word_gamma[i]
                    while j != i:
                        Omega[ii+1+j/jj] ^= (1<<(j%jj))
                        j = word_gamma[j]
                    i += 1
                    while i < nwords and not Omega[ii+1+i/jj]&(1<<(i%jj)):
                        i += 1
                # Phi stores the columns fixed by the automorphism
                for i from 0 <= i < ncols:
                    if col_gamma[i] == i:
                        Phi[ii] ^= (1 << i)
                for i from 0 <= i < nwords:
                    if word_gamma[i] == i:
                        Phi[ii+1+i/jj] ^= (1<<(i%jj))

                # Now incorporate the automorphism into Theta
                j = Theta.merge_perm(col_gamma, word_gamma)

                # j stores whether anything happened or not- if not, then the automorphism we have
                # discovered is already in the subgroup spanned by the generators we have output
                if not j: state = 11; continue

                # otherwise, we have a new generator, so record it:
                self.record_automorphism(col_gamma, ncols)
                # The variable tvc was set to be the minimum element of W[k] the last time the
                # algorithm came up to meet zeta. At this point, we were considering the new
                # possibilities for descending away from zeta at this level.
                # if this is still a minimum cell representative of Theta, even in light
                # of this new automorphism, then the current branch off of zeta hasn't been
                # found equivalent to one already searched yet, so there may still be a
                # better canonical label downward.
                if tvc & nu.flag:
                    i = tvc^nu.flag
                    if Theta.wd_min_cell_rep[Theta.wd_find(i)] == i:
                        state = 11; continue
                else:
                    if Theta.col_min_cell_rep[Theta.col_find(tvc)] == tvc:
                        state = 11; continue

                # Otherwise, proceed to where zeta meets nu:
                k = h
                # TODO: Consider removing the following lines
                if k >= 2*self.radix: raise RuntimeError(\
                    "A counterexample to an assumption the author made while writing this software has been encountered.")
                state = 13

            elif state == 11: # We have just found a new automorphism, and deduced that there may
                # be a better canonical label below the current branch off of zeta. So go to where
                # nu meets rho
                k = hb
                # TODO: Consider removing the following lines
                if k >= 2*self.radix: raise RuntimeError(\
                    "A counterexample to an assumption the author made while writing this software has been encountered.")
                state = 12

            elif state == 12: # Coming here from either state 6 or 11, the algorithm has discovered
                              # some new information. 11 came from 10, where a new line in Omega and
                              # Phi was just recorded, and 6 stored information about implicit auto-
                              # morphisms in Omega and Phi
                if e[k] == 1:
                    # this means that the algorithm has come upward to this position (in state 17)
                    # before, so we have already intersected W[k] with the bulk of Omega and Phi, but
                    # we should still catch up with the latest ones
                    ii = self.Phi_size*l
                    jj = self.Phi_size*k
                    j = 1 + nwords/self.radix
                    if nwords%self.radix:
                        j += 1
                    W[jj] &= Omega[ii]
                    for i from 0 < i < j:
                        W[jj+i] &= Omega[ii+i]
                state = 13

            elif state == 13: # hub state
                if k == -1: state = -1; continue # exit point

                if k > h: state = 17; continue # we are still on the same principal branch from zeta

                if k == h: state = 14; continue # update the stabilizer index and check for new splits,
                                                # since we have returned to a partition of zeta
                # otherwise k < h, hence we have just backtracked up zeta, and are one level closer to done
                h = k
                tvc = 0
                jj = self.Phi_size*k
                if W[jj]:
#                    print 'W[jj]', W[jj]
#                    print tvc
                    while not (1 << tvc) & W[jj]:
                        tvc += 1
                else:
                    ii = 0
                    while not W[jj+1+ii]:
                        ii += 1
                    while not W[jj+1+ii] & (1 << tvc):
                        tvc += 1
                    tvc = (ii*self.radix + tvc) ^ nu.flag
                # now tvc points to the minimal cell representative of W[k]
                state = 14

            elif state == 14: # see if there are any more splits to make from this level of zeta (see state 17)
#                print Theta
                if v[k]&nu.flag == tvc&nu.flag:
                    if tvc&nu.flag:
 #                       print 'v[k] is word', v[k]^nu.flag
  #                      print 'tvc is word', tvc^nu.flag
                        if Theta.wd_find(v[k]^nu.flag) == Theta.wd_find(tvc^nu.flag):
                            index += 1
   #                         print 'index', index
                    else:
    #                    print 'v[k] is col', v[k]
     #                   print 'tvc is col', tvc
                        if Theta.col_find(v[k]) == Theta.col_find(tvc):
                            index += 1
      #                      print 'index', index
                            # keep tabs on how many elements are in the same cell of Theta as tvc
                # find the next split
                jj = self.Phi_size*k
                if v[k]&nu.flag:
                    ii = self.radix
                    i = (v[k]^nu.flag) + 1
                    while i < nwords and not (1 << i%ii) & W[jj+1+i/ii]:
                        i += 1
                    if i < nwords:
                        v[k] = i^nu.flag
                    else:
                        # there is no new split at this level
                        state = 16; continue
                    # new split column better be a minimal representative in Theta, or wasted effort
                    if Theta.wd_min_cell_rep[Theta.wd_find(i)] == i:
                        state = 15
                    else:
                        state = 14
                else:
                    i = v[k] + 1
                    while i < ncols and not (1 << i) & W[jj]:
                        i += 1
                    if i < ncols:
                        v[k] = i
                    else:
                        # there is no new split at this level
                        state = 16; continue
                    # new split column better be a minimal representative in Theta, or wasted effort
#                    print 'checking whether v[k] is a minimum cell rep of theta'
#                    print 'Theta.col_find(v[k]) = ', Theta.col_find(v[k])
#                    print 'Theta.col_min_cell_rep(^)', Theta.col_min_cell_rep[Theta.col_find(v[k])]
#                    print 'v[k]', v[k]
                    if Theta.col_min_cell_rep[Theta.col_find(v[k])] == v[k]:
                        state = 15
                    else:
                        state = 14

            elif state == 15: # split out the column v[k]
                # hh is smallest such that nu[hh] satisfies Lemma 2.25. If it is larger than k+1,
                # it must be modified, since we are changing that part
                if k + 1 < hh:
                    hh = k + 1
                # hzf is maximal such that indicators line up for nu and zeta
                if k < hzf__h_zeta:
                    hzf__h_zeta = k
                # hzb is longest such that nu and rho have the same indicators
                if hzb__h_rho >= k:
                    hzb__h_rho = k
                    qzb = 0
                state = 2

            elif state == 16: # backtrack up zeta, updating information about stabilizer vector
                jj = self.Phi_size*k
                if W[jj]:
                    i = W[jj]
                    j = ham_wts[i & 65535] + ham_wts[(i >> 16) & 65535]
                else:
                    i = 0; j = 0
                    ii = self.radix
                    while i*ii < nwords:
                        iii = W[jj+1+i]
                        j += ham_wts[iii & 65535] + ham_wts[(iii >> 16) & 65535]
                        i += 1
                if j == index and ht == k + 1: ht = k
#                print "POINT A, index =", index
                size = size*index
                # (POINT A)
                index = 0
                k -= 1
                if hb > k: # update hb since we are backtracking
                    hb = k
                state = 13

            elif state == 17: # see if there are any more splits to make from this level of nu (and not zeta)

                jjj = self.Phi_size*k
                if e[k] == 0: # now is the time to narrow down W[k] by Omega and Phi
                    # intersect W[k] with each Omega[i] such that v[0]...v[k-1] is in Phi[i]
                    jj = self.Phi_size*self.L
                    iii = nwords/self.radix
                    if nwords%self.radix:
                        iii += 1
                    for ii from 0 <= ii < iii:
                        Phi[jj+ii] = 0
                    for ii from 0 <= ii < k:
                        if v[ii]&nu.flag:
                            i = v[ii]^nu.flag
                            Phi[jj+1+i/self.radix] ^= (1 << i%self.radix)
                        else:
                            Phi[jj] ^= (1 << v[ii])
                    for i from 0 <= i <= l:
                        ii = self.Phi_size*i
                        iiii = 1
                        for j from 0 <= j < iii:
                            if Phi[ii + j] & Phi[jj + j] != Phi[jj + j]:
                                iiii = 0
                                break
                        if iiii:
                            for j from 0 <= j < iii:
                                W[jjj + j] &= Omega[ii + j]
                e[k] = 1

                # see if there is a vertex to split out
                if nu.flag&v[k]:
                    i = (v[k]^nu.flag) + 1
                    while i < nwords:
                        i += 1
                        if (1 << i%self.radix) & W[jjj+1+i/self.radix]: break
                    if i < nwords:
                        v[k] = i^nu.flag
                        state = 15; continue
                else:
                    i = v[k] + 1
                    while i < ncols:
                        i += 1
                        if (1 << i) & W[jjj]: break
                    if i < ncols:
                        v[k] = i
                        state = 15; continue

                k -= 1
                state = 13

            elif state == 18: # the first time nu becomes a discrete partition: set up zeta, our "identity" leaf
                # initialize counters for zeta:
                h = k # zeta[h] == nu[h]
                ht = k # nodes descended from zeta[ht] are all equivalent
                hzf__h_zeta = k # max such that indicators for zeta and nu agree
                zeta = PartitionStack(nu)
                zeta.find_basis(ham_wts)
                # (POINT B)
                k -= 1
                rho = PartitionStack(nu)
                # initialize counters for rho:
                k_rho = k+1 # number of partitions in rho
                hzb__h_rho = k # max such that indicators for rho and nu agree - BDM had k+1
                hb = k # rho[hb] == nu[hb] - BDM had k+1
                qzb = 0 # Lambda[k] == zb[k], so...
                state = 13

        # end big while loop
        rho.find_basis(ham_wts)
        for i from 0 <= i < ncols:
            self.labeling[rho.col_ents[i]] = i
        for i from 0 <= i < nrows:
            self.labeling[i+ncols] = rho.basis_locations[i]
        self.aut_gp_size = size




