"""
Partition backtrack functions for matrices

DOCTEST:
    sage: import sage.groups.perm_gps.partn_ref.refinement_matrices

REFERENCE:

    [1] McKay, Brendan D. Practical Graph Isomorphism. Congressus Numerantium,
        Vol. 30 (1981), pp. 45-87.

    [2] Leon, Jeffrey. Permutation Group Algorithms Based on Partitions, I:
        Theory and Algorithms. J. Symbolic Computation, Vol. 12 (1991), pp.
        533-583.

"""

#*****************************************************************************
#      Copyright (C) 2006 - 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include '../../../misc/bitset.pxi'
from sage.misc.misc import uniq
from sage.matrix.constructor import Matrix

cdef class MatrixStruct:

    def __new__(self, matrix):
        cdef int i, j
        cdef int *num_rows
        self.degree = matrix.ncols()
        self.nwords = matrix.nrows()
        cdef NonlinearBinaryCodeStruct S_temp
        self.matrix = matrix

        self.symbols = uniq(self.matrix.list())
        self.symbols.remove(0)
        self.nsymbols = len(self.symbols)

        self.symbol_structs = []
        num_rows = <int *> sage_malloc(self.nsymbols * sizeof(int))
        self.temp_col_ps = PS_new(self.degree, 1)
        if num_rows is NULL or self.temp_col_ps is NULL:
            if num_rows is not NULL:
                sage_free(num_rows)
            if self.temp_col_ps is not NULL:
                PS_dealloc(self.temp_col_ps)
            raise MemoryError

        for i from 0 <= i < self.nsymbols:
            num_rows[i] = 0
        for row in self.matrix.rows():
            row = uniq(row.list())
            if 0 in row: row.remove(0)
            for s in row:
                num_rows[self.symbols.index(s)] += 1
        for i from 0 <= i < self.nsymbols:
            S_temp = NonlinearBinaryCodeStruct( (self.degree, num_rows[i]) )
            self.symbol_structs.append(S_temp)

        for i from 0 <= i < self.nsymbols:
            num_rows[i] = 0
        for row in self.matrix.rows():
            row_list = row.list()
            row_list.reverse()
            for i in row.nonzero_positions():
                s = row[i]
                j = self.symbols.index(s)
                S_temp = <NonlinearBinaryCodeStruct>self.symbol_structs[j]
                bitset_set( &S_temp.words[num_rows[j]], i)
                if row_list.count(s) == 1 or row_list.index(s) == self.degree - i - 1:
                    num_rows[j] += 1
        sage_free(num_rows)
        self.output = NULL

    def __dealloc__(self):
        PS_dealloc(self.temp_col_ps)
        if self.output is not NULL:
            mpz_clear(self.output.order)
            sage_free(self.output.generators)
            sage_free(self.output.base)
            sage_free(self.output.relabeling)
            sage_free(self.output)

    def display(self):
        print self.matrix
        print
        cdef int i,j=0
        cdef NonlinearBinaryCodeStruct S_temp
        for S in self.symbol_structs:
            S_temp = <NonlinearBinaryCodeStruct>S
            for i from 0 <= i < S_temp.nwords:
                print bitset_string(&S_temp.words[i])
            print self.symbols[j]
            print
            j += 1

    def run(self, partition=None):
        """
        Perform the canonical labeling and automorphism group computation,
        storing results to self.

        INPUT:
        partition -- an optional list of lists partition of the columns.
            default is the unit partition.

        EXAMPLES:
            sage: from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct

            sage: M = MatrixStruct(matrix(GF(3),[[0,1,2],[0,2,1]]))
            sage: M.run()
            sage: M.automorphism_group()
            ([[0, 2, 1]], 2, [1])
            sage: M.canonical_relabeling()
            [0, 1, 2]

            sage: M = MatrixStruct(matrix(GF(3),[[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]]))
            sage: M.automorphism_group()[1] == 6
            True

            sage: M = MatrixStruct(matrix(GF(3),[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2]]))
            sage: M.automorphism_group()[1] == factorial(14)
            True

        """
        cdef int **part, i, j
        if partition is None:
            partition = [range(self.degree)]
        part = <int **> sage_malloc((len(partition)+1) * sizeof(int *))
        if part is NULL:
            raise MemoryError
        for i from 0 <= i < len(partition):
            part[i] = <int *> sage_malloc((len(partition[i])+1) * sizeof(int))
            if part[i] is NULL:
                for j from 0 <= j < i:
                    sage_free(part[j])
                sage_free(part)
                raise MemoryError
            for j from 0 <= j < len(partition[i]):
                part[i][j] = partition[i][j]
            part[i][len(partition[i])] = -1
        part[len(partition)] = NULL

        self.output = get_aut_gp_and_can_lab(self, part, self.degree, &all_matrix_children_are_equivalent, &refine_matrix, &compare_matrices, 1, 1, 1)

        for i from 0 <= i < len(partition):
            sage_free(part[i])
        sage_free(part)


    def automorphism_group(self):
        """
        Returns a list of generators of the automorphism group, along with its
        order and a base for which the list of generators is a strong generating
        set.

        EXAMPLE: (For more examples, see self.run())
            sage: from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct

            sage: M = MatrixStruct(matrix(GF(3),[[0,1,2],[0,2,1]]))
            sage: M.automorphism_group()
            ([[0, 2, 1]], 2, [1])

        """
        cdef int i, j
        cdef object generators, base
        cdef Integer order
        if self.output is NULL:
            self.run()
        generators = []
        for i from 0 <= i < self.output.num_gens:
            generators.append([self.output.generators[i*self.degree + j] for j from 0 <= j < self.degree])
        order = Integer()
        mpz_set(order.value, self.output.order)
        base = [self.output.base[i] for i from 0 <= i < self.output.base_size]
        return generators, order, base

    def canonical_relabeling(self):
        """
        Returns a canonical relabeling (in list permutation format).

        EXAMPLES: (For more examples, see self.run())
            sage: from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct

            sage: M = MatrixStruct(matrix(GF(3),[[0,1,2],[0,2,1]]))
            sage: M.canonical_relabeling()
            [0, 1, 2]

        """
        cdef int i
        if self.output is NULL:
            self.run()
        return [self.output.relabeling[i] for i from 0 <= i < self.degree]

cdef int refine_matrix(PartitionStack *PS, object S, int *cells_to_refine_by, int ctrb_len):
    cdef MatrixStruct M = <MatrixStruct> S
    cdef int i, temp_inv, invariant = 1
    cdef bint changed = 1
    while changed:
        PS_copy_from_to(PS, M.temp_col_ps)
        for BCS in M.symbol_structs:
            temp_inv = refine_by_bip_degree(PS, BCS, cells_to_refine_by, ctrb_len)
            invariant *= temp_inv + 1
        if (memcmp(PS.entries, M.temp_col_ps.entries, M.degree * sizeof(int)) == 0
         and memcmp(PS.levels, M.temp_col_ps.levels, M.degree * sizeof(int)) == 0):
            changed = 0
    return invariant

cdef int compare_matrices(int *gamma_1, int *gamma_2, object S):
    cdef MatrixStruct MS = <MatrixStruct> S
    M = MS.matrix
    cdef int i
    M1 = Matrix(M.base_ring(), M.nrows(), M.ncols(), sparse=M.is_sparse())
    M2 = Matrix(M.base_ring(), M.nrows(), M.ncols(), sparse=M.is_sparse())
    for i from 0 <= i < M.ncols():
        M1.set_column(i, M.column(gamma_1[i]))
        M2.set_column(i, M.column(gamma_2[i]))
    return cmp(sorted(M1.rows()), sorted(M2.rows()))

cdef bint all_matrix_children_are_equivalent(PartitionStack *PS, object S):
    cdef MatrixStruct M = <MatrixStruct> S
    cdef bint equiv = 1
    for BCS in M.symbol_structs:
        equiv &= all_binary_children_are_equivalent(PS, BCS)
    return equiv

def random_tests(t=10.0, nrows_max=50, ncols_max=50, nsymbols_max=20, perms_per_matrix=10, density_range=(.1,.9)):
    """
    Tests to make sure that C(gamma(M)) == C(M) for random permutations gamma
    and random matrices M.

    INPUT:
    t -- run tests for approximately this many seconds
    nrows_max -- test matrices with at most this many rows
    ncols_max -- test matrices with at most this many columns
    perms_per_matrix -- test each matrix with this many random permutations
    nsymbols_max -- maximum number of distinct symbols in the matrix

    DISCUSSION:

    Until t seconds have elapsed, this code generates a random matrix M on
    at most ncols_max columns and at most nrows_max rows. The density of entries
    in the basis is chosen randomly between 0 and 1.

    For each matrix M generated, we uniformly generate perms_per_matrix random
    permutations and verify that the canonical labels of M and the image of M
    under the generated permutation are equal.

    DOCTEST:
        sage: import sage.groups.perm_gps.partn_ref.refinement_matrices
        sage: sage.groups.perm_gps.partn_ref.refinement_matrices.random_tests()
        All passed: ... random tests on ... matrices.
        sage: sage.groups.perm_gps.partn_ref.refinement_matrices.random_tests(180.0, 100, 200, 40) # long time
        All passed: ... random tests on ... matrices.

    """
    from sage.misc.misc import walltime
    from sage.misc.prandom import random, randint
    from sage.combinat.permutation import Permutations
    from sage.matrix.constructor import random_matrix, matrix
    from sage.rings.finite_field import FiniteField as GF
    from sage.rings.arith import next_prime
    cdef int h, i, j, nrows, k, num_tests = 0, num_matrices = 0, passed = 1
    cdef MatrixStruct M, N
    t_0 = walltime()
    while walltime(t_0) < t:
        p = random()*(density_range[1]-density_range[0]) + density_range[0]
        nrows = randint(1, nrows_max)
        ncols = randint(1, ncols_max)
        nsymbols = next_prime(randint(1, nsymbols_max))
        S = Permutations(ncols)
        MM = random_matrix(GF(nsymbols), nrows, ncols, False, p)
        M = MatrixStruct( MM )
        M.run()

        for i from 0 <= i < perms_per_matrix:
            perm = [a-1 for a in list(S.random_element())]
            NN = matrix(GF(nsymbols), nrows, ncols)
            for j from 0 <= j < ncols:
                NN.set_column(perm[j], MM.column(j))
            N = MatrixStruct(NN)
            # now N is a random permutation of M
            N.run()

            M_relab = M.canonical_relabeling()
            N_relab = N.canonical_relabeling()

            M_C = matrix(GF(nsymbols), nrows, ncols)
            N_C = matrix(GF(nsymbols), nrows, ncols)

            for j from 0 <= j < ncols:
                M_C.set_column(M_relab[j], MM.column(j))
                N_C.set_column(N_relab[j], NN.column(j))

            M_C = matrix(GF(nsymbols), sorted(M_C.rows()))
            N_C = matrix(GF(nsymbols), sorted(N_C.rows()))

            if M_C != N_C:
                print "M:"
                print M
                print "perm:"
                print perm
                return

        num_tests += perms_per_matrix
        num_matrices += 2
    if passed:
        print "All passed: %d random tests on %d matrices."%(num_tests, num_matrices)




