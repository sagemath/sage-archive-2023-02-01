r"""
Sudoku Solver

Given a 9x9 Sudoku puzzle as an integer matrix, this routine finds a single solution.
"""

def sudoku(A):
    r"""
    Solve the 9x9 Sudoku puzzle contained in the matrix `A`.

    INPUT:

    - `A` - a 9x9 matrix with integer entries from 0, 1..9. A `0` indicates an empty square

    OUTPUT:

    matrix - a 9x9 matrix over ZZ containing the first solution found

    ALGORITHM:

    A solution is found by examining blank cells in turn, determining which symbols
    are in use in the corresponding row, column and 3x3 sub-grid, and then making
    recursive calls exhausting all possibilities at that blank cell.  Basically this
    is a depth-first search to the first solution, with pruning accomplished according
    to only the basic requirements of a legitimate completed Sudoku.

    EXAMPLES:

    ::

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

    Now we perturb `A` slightly to make the `(3,5)` entry impossible to complete.

    ..link::

        sage: A[1,4], A[2,4] = 7, 9
        sage: A[3,6], A[3,7], A[3,8] = 3, 4, 6
        sage: A
        [5 0 0 0 8 0 0 4 9]
        [0 0 0 5 7 0 0 3 0]
        [0 6 7 3 9 0 0 0 1]
        [1 5 0 0 0 0 3 4 6]
        [0 0 0 2 0 8 0 0 0]
        [0 0 0 0 0 0 0 1 8]
        [7 0 0 0 0 4 1 5 0]
        [0 3 0 0 0 2 0 0 0]
        [4 9 0 0 5 0 0 0 3]
        sage: print sudoku(A)
        None
    """
    # locate first empty square to initiate recursion
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

R9 = range(9)
R10 = range(1,10)

def grid_has_k(A, i, j, k):
    r"""
    Checks for the presence of `k` in the 3x3 subgrid containing
    the location in row `i` and column `j`.

    INPUT:

    - `A` - a 9x9 matrix with entries from 0, 1..9

    - `i` - integer specifying row `i`

    - `j` - integer specifying column `j`

    - `k` - an integer from 1..9

    OUTPUT:

    boolean - True exactly when `k` is present in the 3x3 sub-grid
              that also has the entry in row `i` and column `j`.

    EXAMPLES:

    ::

        sage: from sage.games.sudoku import grid_has_k
        sage: B = matrix(ZZ, 9, 9, [ [0,0,0,0,1,0,9,0,0], [8,0,0,4,0,0,0,0,0], [2,0,0,0,0,0,0,0,0], [0,7,0,0,3,0,0,0,0], [0,0,0,0,0,0,2,0,4], [0,0,0,0,0,0,0,5,8], [0,6,0,0,0,0,1,3,0], [7,0,0,2,0,0,0,0,0], [0,0,0,8,0,0,0,0,0] ])
        sage: B
        [0 0 0 0 1 0 9 0 0]
        [8 0 0 4 0 0 0 0 0]
        [2 0 0 0 0 0 0 0 0]
        [0 7 0 0 3 0 0 0 0]
        [0 0 0 0 0 0 2 0 4]
        [0 0 0 0 0 0 0 5 8]
        [0 6 0 0 0 0 1 3 0]
        [7 0 0 2 0 0 0 0 0]
        [0 0 0 8 0 0 0 0 0]
        sage: grid_has_k(B, 3, 2, 7)
        True
        sage: grid_has_k(B, 3, 2, 1)
        False
    """
    ii = i//3
    jj = j//3
    for m in range(3*ii,3*ii+3):
        for n in range(3*jj,3*jj+3):
            if A[m,n] == k:
                return True
    return False

def solve_recursive(A, i, j):
    r"""
    Completes a Sudoku puzzle starting at (an empty) square in row `i` and column `j`

    INPUT:

    - `A` - a 9x9 matrix with integer entries from 0, 1..9

    - `i` - integer specifying row `i`

    - `j` - integer specifying column `j`

    OUTPUT:

    matrix - a 9x9 matrix over ZZ

           If square `(i,j)` is non-zero (not empty) then `A` is returned immediately

           If there is no way to complete the puzzle then None is returned

           Otherwise a recursive call is made, which will eventually return a (partial) solution or None

    NOTES:

    As a practical matter this should be called with `(i,j)` being the first
    empty square (in row-major order), though nothing prohibits starting at
    another square (empty or not), and possibly getting a partially complete
    square back as output.

    EXAMPLES:

    ::

    This puzzle has only 17 non-empty squares (hints) and has a unique solution.
    At this writing (2009/04/11) there are no known 16-hint Sudoku puzzles
    with a unique solution.  It is number 3000 in Gordon Royle's database of 17-hint
    uniquely-solvable puzzles [1].

    Here we just test finding a possible completion if
    we started midway through the bottom row.

        sage: from sage.games.sudoku import solve_recursive
        sage: B = matrix(ZZ, 9, 9, [ [0,0,0,0,1,0,9,0,0], [8,0,0,4,0,0,0,0,0], [2,0,0,0,0,0,0,0,0], [0,7,0,0,3,0,0,0,0], [0,0,0,0,0,0,2,0,4], [0,0,0,0,0,0,0,5,8], [0,6,0,0,0,0,1,3,0], [7,0,0,2,0,0,0,0,0], [0,0,0,8,0,0,0,0,0] ])
        sage: B
        [0 0 0 0 1 0 9 0 0]
        [8 0 0 4 0 0 0 0 0]
        [2 0 0 0 0 0 0 0 0]
        [0 7 0 0 3 0 0 0 0]
        [0 0 0 0 0 0 2 0 4]
        [0 0 0 0 0 0 0 5 8]
        [0 6 0 0 0 0 1 3 0]
        [7 0 0 2 0 0 0 0 0]
        [0 0 0 8 0 0 0 0 0]
        sage: solve_recursive(B, 8, 5)
        [0 0 0 0 1 0 9 0 0]
        [8 0 0 4 0 0 0 0 0]
        [2 0 0 0 0 0 0 0 0]
        [0 7 0 0 3 0 0 0 0]
        [0 0 0 0 0 0 2 0 4]
        [0 0 0 0 0 0 0 5 8]
        [0 6 0 0 0 0 1 3 0]
        [7 0 0 2 0 0 0 0 0]
        [0 0 0 8 0 1 4 2 5]

    Now we perturb `B` to make the (8,7) square impossible to complete.
    And we start at an empty square a little further back up the matrix.

    .. link::

        sage: B[0,7], B[1,7], B[2,7] = 2, 6, 7
        sage: B[8,0], B[8, 1] = 4, 9
        sage: B
        [0 0 0 0 1 0 9 2 0]
        [8 0 0 4 0 0 0 6 0]
        [2 0 0 0 0 0 0 7 0]
        [0 7 0 0 3 0 0 0 0]
        [0 0 0 0 0 0 2 0 4]
        [0 0 0 0 0 0 0 5 8]
        [0 6 0 0 0 0 1 3 0]
        [7 0 0 2 0 0 0 0 0]
        [4 9 0 8 0 0 0 0 0]
        sage: print solve_recursive(B, 7, 1)
        None

    REFERENCES:

    - [1] Gordon Royle, Minimum Sudoku, http://people.csse.uwa.edu.au/gordon/sudokumin.php, (2009/04/11)
    """
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
          # locate "next" empty square
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

