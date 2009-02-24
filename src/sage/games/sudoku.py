"""
Sudoku Solver

Given a 9x9 Sudoku puzzle as an integer matrix, the program solves
it.
"""

def sudoku(A):
    """
    Solve the 9x9 Sudoku puzzle defined by the matrix `A`.

    EXAMPLE::

        sage: A = matrix(ZZ,9,[5,0,0, 0,8,0, 0,4,9, 0,0,0, 5,0,0, 0,3,0, 0,6,7, 3,0,0, 0,0,1,  1,5,0, 0,0,0, 0,0,0,  0,0,0, 2,0,8, 0,0,0,    0,0,0, 0,0,0, 0,1,8,     7,0,0, 0,0,4, 1,5,0,   0,3,0, 0,0,2, 0,0,0,  4,9,0, 0,5,0, 0,0,3])
        sage: A
        [5 0 0 0 8 0 0 4 9]
        [0 0 0 5 0 0 0 3 0]
        [0 6 7 3 0 0 0 0 1]
        [1 5 0 0 0 0 0 0 0]
        [0 0 0 2 0 8 0 0 0]
        [0 0 0 0 0 0 0 1 8]
        [7 0 0 0 0 4 1 5 0]
        [0 3 0 0 0 2 0 0 0]
        [4 9 0 0 5 0 0 0 3]
        sage: sudoku(A)
        [5 1 3 6 8 7 2 4 9]
        [8 4 9 5 2 1 6 3 7]
        [2 6 7 3 4 9 5 8 1]
        [1 5 8 4 6 3 9 7 2]
        [9 7 4 2 1 8 3 6 5]
        [3 2 6 7 9 5 4 1 8]
        [7 8 2 9 3 4 1 5 6]
        [6 3 5 1 7 2 8 9 4]
        [4 9 1 8 5 6 7 2 3]
    """
    i = 0
    j = 0
    while A[i,j] != 0:
        if j < 8:
            j += 1
        elif i < 8:
            j = 0
            i += 1
        else:
            break
    return solve_recursive(A, i, j)

copies = 0

R9 = range(9)
R10 = range(1,10)

def grid_has_k(A, i, j, k):
    r"""
    Return ``True`` precisely if the 3x3 submatrix that contains the
    `i,j` position of `A` already has a `k` in it.
    """
    ii = i//3
    jj = j//3
    for m in range(3*ii,3*ii+3):
        for n in range(3*jj,3*jj+3):
            if A[m,n] == k:
                return True
    return False

def solve_recursive(A, i, j):
    # determine excluded possibilities
    if A[i,j] != 0:
        return A
    v = []
    for k in R9:
       z = A[i,k]
       if z != 0:
          v.append(int(z))
       z = A[k,j]
       if z != 0:
          v.append(int(z))
    v = set(v)
    if len(v) == 9:
       return None  # failure

    # try each allowed possibility for the given ij position
    for k in R10:
       if not k in v:
          # We can make the move B[i,j] = k, only if the
          # 3x3 grid that contains the i,j position doesn't
          # have a k in it already.
          if grid_has_k(A, i,j,k):
              continue
          B = A.__copy__()
          B[i,j] = k
          ii = i; jj = j
          while B[ii,jj] != 0:
             if jj < 8:
                jj += 1
             elif ii < 8:
                jj = 0
                ii += 1
             else:
                ii = 8; jj = 8
                break
          C = solve_recursive(B, ii,jj)
          if not C is None:
             return C
    return None

