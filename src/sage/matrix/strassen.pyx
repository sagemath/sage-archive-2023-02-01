"""
Generic Asymptotically Fast Strassen Algorithms

Sage implements asymptotically fast echelon form and matrix
multiplication algorithms.
"""

################################################################################
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL).
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
################################################################################


from matrix_window cimport MatrixWindow

include "../ext/interrupt.pxi"


def strassen_window_multiply(C, A,B, cutoff):
    """
    Multiplies the submatrices specified by A and B, places result in
    C. Assumes that A and B have compatible dimensions to be
    multiplied, and that C is the correct size to receive the product,
    and that they are all defined over the same ring.

    Uses strassen multiplication at high levels and then uses
    MatrixWindow methods at low levels. EXAMPLES: The following matrix
    dimensions are chosen especially to exercise the eight possible
    parity combinations that ocould ccur while subdividing the matrix
    in the strassen recursion. The base case in both cases will be a
    (4x5) matrix times a (5x6) matrix.

    ::

        sage: A = MatrixSpace(Integers(2^65), 64, 83).random_element()
        sage: B = MatrixSpace(Integers(2^65), 83, 101).random_element()
        sage: A._multiply_classical(B) == A._multiply_strassen(B, 3)
        True

    AUTHORS:

    - David Harvey
    """
    strassen_window_multiply_c(C, A, B, cutoff)

cdef strassen_window_multiply_c(MatrixWindow C, MatrixWindow A,
                                MatrixWindow B, Py_ssize_t cutoff):
    # todo -- I'm not sure how to interpret "cutoff". Should it be...
    # (a) the minimum side length of the matrices (currently implemented below)
    # (b) the maximum side length of the matrices
    # (c) the total number of entries being multiplied
    # (d) something else entirely?

    cdef Py_ssize_t A_nrows, A_ncols, B_ncols
    A_nrows = A._nrows
    A_ncols = A._ncols        # this should also be the number of rows of B
    B_ncols = B._ncols

    if (A_nrows <= cutoff) or (A_ncols <= cutoff) or (B_ncols <= cutoff):
        # note: this code is only reached if the TOP level is already beneath
        # the cutoff. In a typical large multiplication, the base case is
        # handled directly (see below).
        C.set_to_prod(A, B)
        return

    # Construct windows for the four quadrants of each matrix.
    # Note that if the side lengths are odd we're ignoring the
    # final row/column for the moment.

    cdef Py_ssize_t A_sub_nrows, A_sub_ncols, B_sub_ncols
    A_sub_nrows = A_nrows >> 1
    A_sub_ncols = A_ncols >> 1     # this is also like B_sub_nrows
    B_sub_ncols = B_ncols >> 1

    cdef MatrixWindow A00, A01, A10, A11, B00, B01, B10, B11
    A00 = A.matrix_window(0,           0,           A_sub_nrows, A_sub_ncols)
    A01 = A.matrix_window(0,           A_sub_ncols, A_sub_nrows, A_sub_ncols)
    A10 = A.matrix_window(A_sub_nrows, 0,           A_sub_nrows, A_sub_ncols)
    A11 = A.matrix_window(A_sub_nrows, A_sub_ncols, A_sub_nrows, A_sub_ncols)
    B00 = B.matrix_window(0,           0,           A_sub_ncols, B_sub_ncols)
    B01 = B.matrix_window(0,           B_sub_ncols, A_sub_ncols, B_sub_ncols)
    B10 = B.matrix_window(A_sub_ncols, 0,           A_sub_ncols, B_sub_ncols)
    B11 = B.matrix_window(A_sub_ncols, B_sub_ncols, A_sub_ncols, B_sub_ncols)

    # Allocate temp space.

    cdef MatrixWindow S0, S1, S2, S3, T0, T1 ,T2, T3, Q0, Q1, Q2
    cdef MatrixWindow tmp
    cdef Py_ssize_t tmp_cols, start_row
    tmp_cols = A_sub_ncols
    if (tmp_cols < B_sub_ncols):
        tmp_cols = B_sub_ncols # tmp_cols = max(A_sub_ncols, B_sub_ncols)
    tmp = A.new_empty_window(7*A_sub_nrows + 4*A_sub_ncols, tmp_cols)

    start_row = 0
    S0 = tmp.matrix_window(start_row, 0, A_sub_nrows, A_sub_ncols)
    start_row += A_sub_nrows
    S1 = tmp.matrix_window(start_row, 0, A_sub_nrows, A_sub_ncols)
    start_row += A_sub_nrows
    S2 = tmp.matrix_window(start_row, 0, A_sub_nrows, A_sub_ncols)
    start_row += A_sub_nrows
    S3 = tmp.matrix_window(start_row, 0, A_sub_nrows, A_sub_ncols)
    start_row += A_sub_nrows

    T0 = tmp.matrix_window(start_row, 0, A_sub_ncols, B_sub_ncols)
    start_row += A_sub_ncols
    T1 = tmp.matrix_window(start_row, 0, A_sub_ncols, B_sub_ncols)
    start_row += A_sub_ncols
    T2 = tmp.matrix_window(start_row, 0, A_sub_ncols, B_sub_ncols)
    start_row += A_sub_ncols
    T3 = tmp.matrix_window(start_row, 0, A_sub_ncols, B_sub_ncols)
    start_row += A_sub_ncols

    Q0 = tmp.matrix_window(start_row, 0, A_sub_nrows, B_sub_ncols)
    start_row += A_sub_nrows
    Q1 = tmp.matrix_window(start_row, 0, A_sub_nrows, B_sub_ncols)
    start_row += A_sub_nrows
    Q2 = tmp.matrix_window(start_row, 0, A_sub_nrows, B_sub_ncols)


    # Preparatory matrix additions/subtractions.

    # todo: we can probably save some memory in these
    # operations by reusing some of the buffers, if we interleave
    # these additions with the multiplications (below)

    # (I believe we can save on one S buffer and one T buffer)

    # S0 = A10 + A11,  T0 = B01 - B00
    # S1 = S0 - A00,   T1 = B11 - T0
    # S2 = A00 - A10,  T2 = B11 - B01
    # S3 = A01 - S1,   T3 = B10 - T1

    S0.set_to_sum(A10, A11)
    S1.set_to_diff(S0, A00)
    S2.set_to_diff(A00, A10)
    S3.set_to_diff(A01, S1)

    T0.set_to_diff(B01, B00)
    T1.set_to_diff(B11, T0)
    T2.set_to_diff(B11, B01)
    T3.set_to_diff(B10, T1)


    # The relations we need now are:

    # P0 = A00*B00
    # P1 = A01*B10
    # P2 =  S0*T0
    # P3 =  S1*T1
    # P4 =  S2*T2
    # P5 =  S3*B11
    # P6 = A11*T3

    # U0 = P0 + P1
    # U1 = P0 + P3
    # U2 = U1 + P4
    # U3 = U2 + P6
    # U4 = U2 + P2
    # U5 = U1 + P2
    # U6 = U5 + P5

    # We place the final answer into the matrix:
    # U0 U6
    # U3 U4

    cdef MatrixWindow U0, U6, U3, U4
    U0 = C.matrix_window(0,           0,           A_sub_nrows, B_sub_ncols)
    U6 = C.matrix_window(0,           B_sub_ncols, A_sub_nrows, B_sub_ncols)
    U3 = C.matrix_window(A_sub_nrows, 0,           A_sub_nrows, B_sub_ncols)
    U4 = C.matrix_window(A_sub_nrows, B_sub_ncols, A_sub_nrows, B_sub_ncols)

    if (A_sub_nrows <= cutoff) or (A_sub_ncols <= cutoff) or (B_sub_ncols <= cutoff):
        # This is the base case, so we use MatrixWindow methods directly.

        # This next chunk is arranged so that each output cell gets written
        # to exactly once. This is important because the output blocks might
        # be quite fragmented in memory, whereas our temporary buffers
        # (Q0, Q1, Q2) will be quite localised, so we can afford to do a bit
        # of arithmetic in them.

        Q0.set_to_prod(A00, B00)         # now Q0 holds P0
        Q1.set_to_prod(A01, B10)         # now Q1 holds P1
        U0.set_to_sum(Q0, Q1)            # now U0 is correct
        Q0.add_prod(S1, T1)              # now Q0 holds U1
        Q1.set_to_prod(S2, T2)           # now Q1 holds P4
        Q1.add(Q0)                       # now Q1 holds U2
        Q2.set_to_prod(A11, T3)          # now Q2 holds P6
        U3.set_to_sum(Q1, Q2)            # now U3 is correct
        Q2.set_to_prod(S0, T0)           # now Q2 holds P2
        U4.set_to_sum(Q2, Q1)            # now U4 is correct
        Q0.add(Q2)                       # now Q0 holds U5
        Q2.set_to_prod(S3, B11)          # now Q2 holds P5
        U6.set_to_sum(Q0, Q2)            # now U6 is correct

    else:
        # Recurse into sub-products.

        strassen_window_multiply_c(Q0, A00, B00, cutoff)   # now Q0 holds P0
        strassen_window_multiply_c(Q1, A01, B10, cutoff)   # now Q1 holds P1
        U0.set_to_sum(Q0, Q1)                              # now U0 is correct
        strassen_window_multiply_c(Q1, S1, T1, cutoff)     # now Q1 holds P3
        Q0.add(Q1)                                         # now Q0 holds U1
        strassen_window_multiply_c(Q1, S2, T2, cutoff)     # now Q1 holds P4
        Q1.add(Q0)                                         # now Q1 holds U2
        strassen_window_multiply_c(Q2, A11, T3, cutoff)    # now Q2 holds P6
        U3.set_to_sum(Q1, Q2)                              # now U3 is correct
        strassen_window_multiply_c(Q2, S0, T0, cutoff)     # now Q2 holds P2
        U4.set_to_sum(Q2, Q1)                              # now U4 is correct
        Q0.add(Q2)                                         # now Q0 holds U5
        strassen_window_multiply_c(Q2, S3, B11, cutoff)    # now Q2 holds P5
        U6.set_to_sum(Q0, Q2)                              # now U6 is correct


    # Now deal with the leftover row and/or column (if they exist).

    cdef MatrixWindow B_last_col, C_last_col, B_bulk, A_last_row, C_last_row, B_last_row, A_last_col, C_bulk

    if B_ncols & 1:
        B_last_col = B.matrix_window(0, B_ncols-1, A_ncols, 1)
        C_last_col = C.matrix_window(0, B_ncols-1, A_nrows, 1)
        C_last_col.set_to_prod(A, B_last_col)

    if A_nrows & 1:
        A_last_row = A.matrix_window(A_nrows-1, 0, 1, A_ncols)
        if B_ncols & 1:
            B_bulk = B.matrix_window(0, 0, A_ncols, B_ncols-1)
            C_last_row = C.matrix_window(A_nrows-1, 0, 1, B_ncols-1)
        else:
            B_bulk = B
            C_last_row = C.matrix_window(A_nrows-1, 0, 1, B_ncols)
        C_last_row.set_to_prod(A_last_row, B_bulk)

    if A_ncols & 1:
        A_last_col = A.matrix_window(0, A_ncols-1, A_sub_nrows << 1, 1)
        B_last_row = B.matrix_window(A_ncols-1, 0, 1, B_sub_ncols << 1)
        C_bulk = C.matrix_window(0, 0, A_sub_nrows << 1, B_sub_ncols << 1)
        C_bulk.add_prod(A_last_col, B_last_row)

cdef subtract_strassen_product(MatrixWindow result, MatrixWindow A, MatrixWindow B, Py_ssize_t cutoff):
    cdef MatrixWindow to_sub
    if (cutoff == -1 or result.ncols() <= cutoff or result.nrows() <= cutoff):
        result.subtract_prod(A, B)
    else:
        to_sub = A.new_empty_window(result.nrows(), result.ncols())
        strassen_window_multiply_c(to_sub, A, B, cutoff)
        result.subtract(to_sub)


def strassen_echelon(MatrixWindow A, cutoff):
    """
    Compute echelon form, in place. Internal function, call with
    M.echelonize(algorithm="strassen") Based on work of Robert Bradshaw
    and David Harvey at MSRI workshop in 2006.

    INPUT:


    -  ``A`` - matrix window

    -  ``cutoff`` - size at which algorithm reverts to
       naive gaussian elemination and multiplication must be at least 1.


    OUTPUT: The list of pivot columns

    EXAMPLE::

        sage: A = matrix(QQ, 7, [5, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 3, 1, 0, -1, 0, 0, -1, 0, 1, 2, -1, 1, 0, -1, 0, 1, 3, -1, 1, 0, 0, -2, 0, 2, 0, 1, 0, 0, -1, 0, 1, 0, 1])
        sage: B = A.copy(); B._echelon_strassen(1); B
        [ 1  0  0  0  0  0  0]
        [ 0  1  0 -1  0  1  0]
        [ 0  0  1  0  0  0  0]
        [ 0  0  0  0  1  0  0]
        [ 0  0  0  0  0  0  1]
        [ 0  0  0  0  0  0  0]
        [ 0  0  0  0  0  0  0]
        sage: C = A.copy(); C._echelon_strassen(2); C == B
        True
        sage: C = A.copy(); C._echelon_strassen(4); C == B
        True

    ::

        sage: n = 32; A = matrix(Integers(389),n,range(n^2))
        sage: B = A.copy(); B._echelon_in_place_classical()
        sage: C = A.copy(); C._echelon_strassen(2)
        sage: B == C
        True

    TESTS::

        sage: A = matrix(Integers(7), 4, 4, [1,2,0,3,0,0,1,0,0,1,0,0,0,0,0,1])
        sage: B = A.copy(); B._echelon_in_place_classical()
        sage: C = A.copy(); C._echelon_strassen(2)
        sage: B == C
        True

    ::

        sage: A = matrix(Integers(7), 4, 4, [1,0,5,0,2,0,3,6,5,1,2,6,4,6,1,1])
        sage: B = A.copy(); B._echelon_in_place_classical()
        sage: C = A.copy(); C._echelon_strassen(2)
        sage: B == C
        True

    AUTHORS:

    - Robert Bradshaw
    """
    if cutoff < 1:
        raise ValueError, "cutoff must be at least 1"
    _sig_on
    strassen_echelon_c(A, cutoff, A._matrix._strassen_default_cutoff(A._matrix))
    _sig_off

cdef strassen_echelon_c(MatrixWindow A, Py_ssize_t cutoff, Py_ssize_t mul_cutoff):
    # The following notation will be used in the comments below, which should be understood to give
    # the general idea of what's going on, as if there were no inconvenient non-pivot columns.
    # The original matrix is given by [ A B ]
    #                                 [ C D ]
    # For compactness, let A' denote the inverse of A
    # top_left, top_right, bottom_left, and bottom_right loosely correspond to A, B, C, and D respectively,
    # however, the "cut" between the top and bottom rows need not be the same.

    cdef Py_ssize_t nrows, ncols
    nrows = A.nrows()
    ncols = A.ncols()

    if (nrows <= cutoff or ncols <= cutoff):
        return A.echelon_in_place()

    cdef Py_ssize_t top_h, bottom_cut, bottom_h, bottom_start, top_cut
    cdef Py_ssize_t prev_pivot_count
    cdef Py_ssize_t split
    split = nrows / 2

    cdef MatrixWindow top, bottom, top_left, top_right, bottom_left, bottom_right, clear

    top = A.matrix_window(0, 0, split, ncols)
    bottom = A.matrix_window(split, 0, nrows-split, ncols)

    top_pivots = strassen_echelon_c(top, cutoff, mul_cutoff)
    # effectively "multiplied" top row by A^{-1}
    #     [  I  A'B ]
    #     [  C   D  ]

    top_pivot_intervals = int_range(top_pivots)
    top_h = len(top_pivots)

    if top_h == 0:
        #                                 [ 0 0 ]
        # the whole top is a zero matrix, [ C D ]. Run echelon on the bottom
        bottom_pivots = strassen_echelon_c(bottom, cutoff, mul_cutoff)
        #             [  0   0  ]
        # we now have [  I  C'D ], proceed to sorting

    else:
        bottom_cut = max(top_pivots) + 1
        bottom_left = bottom.matrix_window(0, 0, nrows-split, bottom_cut)

        if top_h == ncols:
            bottom.set_to_zero()
            # [ I ]
            # [ 0 ]
            # proceed to sorting

        else:
            if bottom_cut == top_h:
                clear = bottom_left
            else:
                clear = bottom_left.to_matrix().matrix_from_columns(top_pivots).matrix_window() # TODO: read only, can I do this faster? Also below
            # Subtract off C time top from the bottom_right
            if bottom_cut < ncols:
                bottom_right = bottom.matrix_window(0, bottom_cut, nrows-split, ncols-bottom_cut)
                subtract_strassen_product(bottom_right, clear, top.matrix_window(0, bottom_cut, top_h, ncols-bottom_cut), mul_cutoff);
            # [  I      A'B   ]
            # [  *   D - CA'B ]

            # Now subtract off C times the top from the bottom_left (pivots -> 0)
            if bottom_cut == top_h:
                bottom_left.set_to_zero()
                bottom_start = bottom_cut

            else:
                for cols in top_pivot_intervals:
                    bottom_left.matrix_window(0, cols[0], nrows-split, cols[1]).set_to_zero()
                non_pivots = int_range(0, bottom_cut) - top_pivot_intervals
                for cols in non_pivots:
                    if cols[0] == 0: continue
                    prev_pivot_count = len(top_pivot_intervals - int_range(cols[0]+cols[1], bottom_cut - cols[0]+cols[1]))
                    subtract_strassen_product(bottom_left.matrix_window(0, cols[0], nrows-split, cols[1]),
                                              clear.matrix_window(0, 0, nrows-split, prev_pivot_count),
                                              top.matrix_window(0, cols[0], prev_pivot_count, cols[1]),
                                              mul_cutoff)
                bottom_start = non_pivots._intervals[0][0]
            # [  I      A'B   ]
            # [  0   D - CA'B ]

            # Now recursively do echelon form on the bottom
            bottom_pivots_rel = strassen_echelon_c(bottom.matrix_window(0, bottom_start, nrows-split, ncols-bottom_start), cutoff, mul_cutoff)
            # [  I  A'B ]
            # [  0  I F ]
            bottom_pivots = []
            for pivot in bottom_pivots_rel:
                bottom_pivots.append(pivot + bottom_start)
            bottom_h = len(bottom_pivots)

            if bottom_h == 0:
                pass
                # [  I  A'B ]
                # [  0   0  ]
                # proceed to sorting

            else:
                #     [  I  A'B ]  =  [  I  E  G  ]
                # let [  0  I F ]  =  [  0  I  F  ]
                top_cut = max(max(bottom_pivots) + 1, bottom_cut)

                # Note: left with respect to leftmost non-zero column of bottom
                top_left = top.matrix_window(0, bottom_start, top_h, top_cut - bottom_start)

                if bottom_h + top_h < ncols:

                    if top_cut - bottom_start == bottom_h:
                        clear = top_left
                    else:
                        clear = top_left.to_matrix().matrix_from_columns(bottom_pivots_rel).matrix_window()

                # subtract off E times bottom from top right
                if top_cut < ncols:

                    top_right = top.matrix_window(0, top_cut, top_h, ncols - top_cut)
                    subtract_strassen_product(top_right, clear, bottom.matrix_window(0, top_cut, bottom_h, ncols - top_cut), mul_cutoff);

                # [  I  *  G - EF ]
                # [  0  I     F   ]

                # Now subtract of E times bottom from top left
                if top_cut - bottom_start == bottom_h:
                    top_left.set_to_zero()

                else:
                    bottom_pivot_intervals = int_range(bottom_pivots)
                    non_pivots = int_range(bottom_start, top_cut - bottom_start) - bottom_pivot_intervals - top_pivot_intervals
                    for cols in non_pivots:
                        if cols[0] == 0: continue
                        prev_pivot_count = len(bottom_pivot_intervals - int_range(cols[0]+cols[1], top_cut - cols[0]+cols[1]))
                        subtract_strassen_product(top.matrix_window(0, cols[0], top_h, cols[1]),
                                                  clear.matrix_window(0, 0, top_h, prev_pivot_count),
                                                  bottom.matrix_window(0, cols[0], prev_pivot_count, cols[1]),
                                                  mul_cutoff)
                    for cols in bottom_pivot_intervals:
                        top.matrix_window(0, cols[0], top_h, cols[1]).set_to_zero()

                # [  I  0  G - EF ]
                # [  0  I     F   ]
                # proceed to sorting

    # subrows already sorted...maybe I could do this more efficiently in cases with few pivot columns (e.g. merge sort)

    pivots = top_pivots
    pivots.extend(bottom_pivots)
    pivots.sort()

    cdef Py_ssize_t i, cur_row
    for cur_row from 0 <= cur_row < len(pivots):
        pivot = pivots[cur_row]
        for i from cur_row <= i < nrows:
            if not A.element_is_zero(i, pivot):
                break
        if i > cur_row and i < nrows:
            A.swap_rows(i, cur_row)

    return pivots



################################
# lots of room for optimization....
# eventually, should I just pass these around rather than lists of ints for pivots?
# would need new from_cols
class int_range:
    r"""
    Useful class for dealing with pivots in the strassen echelon, could
    have much more general application

    AUTHORS:

    - Robert Bradshaw
    """
    def __init__(self, indices=None, range=None):
        if indices is None:
            self._intervals = []
            return
        elif not range is None:
            self._intervals = [(int(indices), int(range))]
        else:
            self._intervals = []
            if len(indices) == 0:
                return
            indices.sort()
            start = None
            last = None
            for ix in indices:
                if last is None:
                    start = ix
                elif ix-last > 1:
                    self._intervals.append((start, last-start+1))
                    start = ix
                last = ix
            self._intervals.append((start, last-start+1))

    def __repr__(self):
        return str(self._intervals)

    def intervals(self):
        return self._intervals

    def to_list(self):
        all = []
        for iv in self._intervals:
            for i in range(iv[0], iv[0]+iv[1]):
                all.append(i)
        return all

    def __iter__(self):
        return self._intervals.__iter__()

    def __len__(self):
        len = 0
        for iv in self._intervals:
            len = len + iv[1]
        return len

    # Yes, these two could be a lot faster...
    # Basically, this class is for abstracting away what I was trying to do by hand in several places
    def __add__(self, right):
        all = self.to_list()
        for i in right.to_list():
            all.append(i)
        return int_range(all)

    def __sub__(self, right):
        all = self.to_list()
        for i in right.to_list():
            if i in all:
                all.remove(i)
        return int_range(all)

    def __mul__(self, right):
        intersection = []
        all = self.to_list()
        for i in right.to_list():
            if i in all:
                intersection.append(i)
        return int_range(intersection)




"""
Useful test code:

def test(n, m, R, c=2):
    A = matrix(R,n,m,range(n*m))
    B = A.copy(); B._echelon_in_place_classical()
    C = A.copy(); C._echelon_strassen(c)
    return B == C

E.g.,
for n in range(5):
    print n, test(2*n,n,Frac(QQ['x']),2)

"""

# This stuff gets tested extensively elsewhere, and the functions
# below aren't callable now without using Pyrex.


##     todo: doc cutoff parameter as soon as I work out what it really means

##     EXAMPLES:
##         The following matrix dimensions are chosen especially to exercise the
##         eight possible parity combinations that ocould ccur while subdividing
##         the matrix in the strassen recursion. The base case in both cases will
##         be a (4x5) matrix times a (5x6) matrix.

##         TODO -- the doctests below are currently not
##         tested/enabled/working -- enable them when linear algebra
##         restructing gets going.

##         sage: dim1 = 64; dim2 = 83; dim3 = 101
##         sage: R = MatrixSpace(QQ, dim1, dim2)
##         sage: S = MatrixSpace(QQ, dim2, dim3)
##         sage: T = MatrixSpace(QQ, dim1, dim3)


##         sage: A = R.random_element(range(-30, 30))
##         sage: B = S.random_element(range(-30, 30))
##         sage: C = T(0)
##         sage: D = T(0)

##         sage: A_window = A.matrix_window(0, 0, dim1, dim2)
##         sage: B_window = B.matrix_window(0, 0, dim2, dim3)
##         sage: C_window = C.matrix_window(0, 0, dim1, dim3)
##         sage: D_window = D.matrix_window(0, 0, dim1, dim3)

##         sage: from sage.matrix.strassen import strassen_window_multiply
##         sage: strassen_window_multiply(C_window, A_window, B_window, 2)   # use strassen method
##         sage: D_window.set_to_prod(A_window, B_window)             # use naive method
##         sage: C_window == D_window
##         True

##         sage: dim1 = 79; dim2 = 83; dim3 = 101
##         sage: R = MatrixSpace(QQ, dim1, dim2)
##         sage: S = MatrixSpace(QQ, dim2, dim3)
##         sage: T = MatrixSpace(QQ, dim1, dim3)

##         sage: A = R.random_element(range(30))
##         sage: B = S.random_element(range(30))
##         sage: C = T(0)
##         sage: D = T(0)

##         sage: A_window = A.matrix_window(0, 0, dim1, dim2)
##         sage: B_window = B.matrix_window(0, 0, dim2, dim3)
##         sage: C_window = C.matrix_window(0, 0, dim1, dim3)

##         sage: strassen_window_multiply(C, A, B, 2)   # use strassen method
##         sage: D.set_to_prod(A, B)                    # use naive method

##         sage: C == D
##         True

