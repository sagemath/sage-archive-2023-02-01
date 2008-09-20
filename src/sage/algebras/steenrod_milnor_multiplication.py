r"""
Milnor multiplication for elements of the mod 2 Steenrod algebra

AUTHORS:
    - John H. Palmieri (2008-07-30: version 0.9)

See Milnor's paper [Mil] for proofs, etc.

To multiply Milnor basis elements $\text{Sq}(r_1, r_2, ...)$ and
$\text{Sq}(s_1, s_2,...)$, form all possible matrices $M$ with rows
and columns indexed starting at 0, with position (0,0) deleted (or
ignored), with $s_i$ equal to the sum of column $i$ for each $i$, and
with $r_j$ equal to the 'weighted' sum of row $j$.  The weights are as
follows: elements from column $i$ are multiplied by $2^i$.  For
example, to multiply $\text{Sq}(2)$ and $\text{Sq}(1,1)$, form the
matrices
\[
\begin{Vmatrix}
* & 1 & 1 \\
2 & 0 & 0
\end{Vmatrix}
\quad \text{and} \quad
\begin{Vmatrix}
* & 0 & 1 \\
0 & 1 & 0
\end{Vmatrix}
\]
(The $*$ is the ignored (0,0)-entry of the matrix.)  For each such
matrix $M$, compute a multinomial coefficient, mod 2: for each
diagonal $\{m_{ij}: i+j=n\}$, compute $(\sum m_{i,j}!) / (m_{0,n}!
m_{1,n-1}!  ... m_{n,0}!)$.  Multiply these together for all $n$.  (To
compute this mod 2, view the entries of the matrix as their base 2
expansions; then this coefficient is zero if and only if there is some
diagonal containing two numbers which have a summand in common in
their base 2 expansion.  For example, if 3 and 10 are in the same
diagonal, the coefficient is zero, because $3=1+2$ and $10=2+8$: they
both have a summand of 2.)

Now, for each matrix with multinomial coefficient 1, let $t_n$ be
the sum of the nth diagonal in the matrix; then
\[
\text{Sq}(r_1, r_2, ...) \text{Sq}(s_1, s_2, ...) = \sum \text{Sq}(t_1, t_2, ...)
\]

The function \code{milnor_multiplication} takes as input two tuples
of non-negative integers, $r$ and $s$, which represent
$\text{Sq}(r)=\text{Sq}(r_1, r_2, ...)$ and
$\text{Sq}(s)=\text{Sq}(s_1, s_2, ...)$; it returns as output a
dictionary whose keys are tuples $t=(t_1, t_2, ...)$ of non-negative
integers, and for each tuple the associated value is the coefficient
of $\text{Sq}(t)$ in the product formula.  Since we are working mod 2,
this coefficient is 1 (if it is zero, the the element is omitted from
the dictionary altogether).

EXAMPLES:
    sage: from sage.algebras.steenrod_milnor_multiplication import milnor_multiplication
    sage: milnor_multiplication((2,), (1,))
    {(0, 1): 1, (3,): 1}
    sage: milnor_multiplication((4,), (2,1))
    {(6, 1): 1, (0, 3): 1, (2, 0, 1): 1}
    sage: milnor_multiplication((2,4), (0,1))
    {(2, 5): 1, (2, 0, 0, 1): 1}

These examples correspond to the following product computations:
\begin{gather*}
\text{Sq}(2) \text{Sq}(1) = \text{Sq}(0,1) + \text{Sq}(3)
\text{Sq}(4) \text{Sq}(2,1) = \text{Sq}(6,1) + \text{Sq}(0,3) + \text{Sq}(2,0,1)
\text{Sq}(2,4) \text{Sq}(0,1) = \text{Sq}(2, 5) + \text{Sq}(2, 0, 0, 1)
\end{gather*}

REFERENCES:

    [Mil] J. W. Milnor, "The Steenrod algebra and its dual, Ann. of Math.
          (2) \textbf{67} (1958), 150--171.
"""

#*****************************************************************************
#       Copyright (C) 2008 John H. Palmieri <palmieri@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#*****************************************************************************

def milnor_multiplication(r,s):
    r"""
    Product of Milnor basis elements r and s.

    INPUT:
        r -- tuple of non-negative integers
        s -- tuple of non-negative integers

    OUTPUT:
        Dictionary of terms of the form (tuple: coeff), where 'tuple' is a
        tuple of non-negative integers and 'coeff' is 1.

    This computes Milnor matrices for the product of $\text{Sq}(r)$
    and $\text{Sq}(s)$, computes their multinomial coefficients, and
    for each matrix whose coefficient is 1, add $\text{Sq}(t)$ to the
    output, where $t$ is the tuple formed by the diagonals sums from
    the matrix.

    EXAMPLES:
        sage: from sage.algebras.steenrod_milnor_multiplication import milnor_multiplication
        sage: milnor_multiplication((2,), (1,))
        {(0, 1): 1, (3,): 1}
        sage: milnor_multiplication((4,), (2,1))
        {(6, 1): 1, (0, 3): 1, (2, 0, 1): 1}
        sage: milnor_multiplication((2,4), (0,1))
        {(2, 5): 1, (2, 0, 0, 1): 1}

    This uses the same algorithm Monks does in his Maple package.
    """
    result = {}
    rows = len(r) + 1
    cols = len(s) + 1
    diags = len(r) + len(s)
    # initialize matrix
    M = range(rows)
    for i in range(rows):
        M[i] = [0]*cols
    for j in range(1,cols):
        M[0][j] = s[j-1]
    for i in range(1,rows):
        M[i][0] = r[i-1]
        for j in range(1,cols):
            M[i][j] = 0
    found = True
    while found:
        # check diagonals
        n = 1
        okay = 1
        diagonal = [0]*diags
        while n <= diags and okay is not None:
            nth_diagonal = [M[i][n-i] for i in range(max(0,n-cols+1), min(1+n,rows))]
            okay = multinomial(nth_diagonal)
            diagonal[n-1] = okay
            n = n + 1
        if okay is not None:
            i = diags - 1
            while i >= 0 and diagonal[i] == 0:
                i = i - 1
            t = tuple(diagonal[:i+1])
            # reduce mod two:
            if result.has_key(t):
                del result[t]
            else:
                result[t] = 1
        # now look for new matrices:
        found = False
        i = 1
        while not found and i < rows:
            sum = M[i][0]
            j = 1
            while not found and j < cols:
                # check to see if column index j is small enough
                if sum >= 2**j:
                    # now check to see if there's anything above this entry
                    # to add to it
                    temp_col_sum = 0
                    for k in range(i):
                        temp_col_sum += M[k][j]
                    if temp_col_sum != 0:
                        found = True
                        for row in range(1,i):
                            M[row][0] = r[row-1]
                            for col in range(1,cols):
                                M[0][col] = M[0][col] + M[row][col]
                                M[row][col] = 0
                        for col in range(1,j):
                            M[0][col] = M[0][col] + M[i][col]
                            M[i][col] = 0
                        M[0][j] = M[0][j] - 1
                        M[i][j] = M[i][j] + 1
                        M[i][0] = sum - 2**j
                    else:
                        sum = sum + M[i][j] * 2**j
                else:
                        sum = sum + M[i][j] * 2**j
                j = j + 1
            i = i + 1
    return result


def multinomial(list):
    """
    Multinomial coefficient of list, mod 2.

    INPUT:
        list -- list of integers

    OUTPUT:
        None if the multinomial coefficient is 0, or sum of list if it is 1

    Given the input $[n_1, n_2, n_3, ...]$, this computes the
    multinomial coefficient $(n_1 + n_2 + n_3 + ...)! / (n_1! n_2!
    n_3! ...)$, mod 2.  The method is roughly this: expand each
    $n_i$ in binary.  If there is a 1 in the same digit for any $n_i$
    and $n_j$ with $i\neq j$, then the coefficient is 0; otherwise, it
    is 1.

    EXAMPLES:
        sage: from sage.algebras.steenrod_milnor_multiplication import multinomial
        sage: multinomial([1,2,4])
        7
        sage: multinomial([1,2,5])
        sage: multinomial([1,2,12,192,256])
        463

    This function does not compute any factorials, so the following are
    actually reasonable to do:
        sage: multinomial([1,65536])
        65537
        sage: multinomial([4,65535])
        sage: multinomial([32768,65536])
        98304
    """
    old_sum = list[0]
    okay = True
    i = 1
    while okay and i < len(list):
        j = 1
        while okay and j <= min(old_sum, list[i]):
            if j & old_sum == j:
                okay = (j & list[i] == 0)
            j = j << 1
        old_sum = old_sum + list[i]
        i = i + 1
    if okay:
        return old_sum
    else:
        return None
