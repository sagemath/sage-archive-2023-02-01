r"""
Milnor multiplication for elements of the odd primary Steenrod algebra

AUTHORS:
    - John H. Palmieri (2008-07-30: version 0.9)

See Milnor's paper [Mil] for proofs, etc.

Fix an odd prime $p$.  There are three steps to multiply Milnor basis
elements $Q_{f_1} Q_{f_2} ... \mathcal{P}(q_1, q_2, ...)$ and
$Q_{g_1} Q_{g_2} ... \mathcal{P}(s_1, s_2,...)$: first, use the formula
\[
\mathcal{P}(q_1, q_2, ...) Q_k = Q_k \mathcal{P}(q_1, q_2, ...)
  + Q_{k+1} \mathcal{P}(q_1 - p^k, q_2, ...)
  + Q_{k+2} \mathcal{P}(q_1, q_2 - p^k, ...)
  + ...
\]
Second, use the fact that the $Q_k$'s form an exterior algebra: $Q_k^2 =
0$ for all $k$, and if $i \neq j$, then $Q_i$ and $Q_j$ anticommute:
$Q_i Q_j = -Q_j Q_i$.  After these two steps, the product is of the form
\[
\sum Q_{e_1} Q_{e_2} ... \mathcal{P}(r_1, r_2, ...) \mathcal{P}(s_1, s_2, ...).
\]
Finally, use Milnor matrices to multiply the pairs of
$\mathcal{P}(...)$ terms: form all possible matrices $M$ with rows and
columns indexed starting at 0, with position (0,0) deleted (or
ignored), with $s_i$ equal to the sum of column $i$ for each $i$, and
with $r_j$ equal to the 'weighted' sum of row $j$.  The weights are as
follows: elements from column $i$ are multiplied by $p^i$.  For
example when $p=5$, to multiply $\mathcal{P}(5)$ and $\mathcal{P}(1,1)$,
form the matrices
\[
\begin{Vmatrix}
* & 1 & 1 \\
5 & 0 & 0
\end{Vmatrix}
\quad \text{and} \quad
\begin{Vmatrix}
* & 0 & 1 \\
0 & 1 & 0
\end{Vmatrix}
\]
(The $*$ is the ignored (0,0)-entry of the matrix.)  For each such
matrix $M$, compute a multinomial coefficient, mod $p$: for each
diagonal $\{m_{ij}: i+j=n\}$, compute $(\sum m_{i,j}!) / (m_{0,n}!
m_{1,n-1}!  ... m_{n,0}!)$.  Multiply these together for all $n$.

Now, for each matrix with nonzero multinomial coefficient $b_M$, let
$t_n$ be the sum of the $n$th diagonal in the matrix; then
\[
\mathcal{P}(r_1, r_2, ...) \mathcal{P}(s_1, s_2, ...)
     = \sum b_M \mathcal{P}(t_1, t_2, ...)
\]
For example when $p=5$, we have
\[
\mathcal{P}(5) \mathcal{P}(1,1) = \mathcal{P}(6,1) + 2 \mathcal{P}(0,2).
\]

The function \code{milnor_multiplication} takes as input two pairs of
tuples of non-negative integers, $(g,q)$ and $(f,s)$, which represent
$Q_{g_1} Q_{g_2} ... \mathcal{P}(q_1, q_2, ...)$ and
$Q_{f_1} Q_{f_2} ... \mathcal{P}(s_1, s_2, ...)$.  It returns as output a
dictionary whose keys are pairs of tuples $(e,t)$ of non-negative
integers, and for each tuple the associated value is the coefficient
in the product formula.

EXAMPLES:
    sage: from sage.algebras.steenrod_milnor_multiplication_odd import milnor_multiplication_odd
    sage: milnor_multiplication_odd(((0,2),(5,)), ((1,),(1,)), 5)
    {((0, 1, 2), (0, 1)): 4, ((0, 1, 2), (6,)): 4}
    sage: milnor_multiplication_odd(((0,2,4),()), ((1,3),()), 7)
    {((0, 1, 2, 3, 4), ()): 6}
    sage: milnor_multiplication_odd(((0,2,4),()), ((1,5),()), 7)
    {((0, 1, 2, 4, 5), ()): 1}
    sage: milnor_multiplication_odd(((),(6,)), ((),(2,)), 3)
    {((), (4, 1)): 1, ((), (8,)): 1, ((), (0, 2)): 1}

These examples correspond to the following product computations:
\begin{gather*}
p=5: \quad Q_0 Q_2 \mathcal{P}(5) Q_1 \mathcal{P}(1) = 4 Q_0 Q_1 Q_2 \mathcal{P}(0,1) + 4 Q_0 Q_1 Q_2 \mathcal{P}(6) \\
p=7: \quad (Q_0 Q_2 Q_4) (Q_1 Q_3) = 6 Q_0 Q_1 Q_2 Q_3 Q_4 \\
p=7: \quad (Q_0 Q_2 Q_4) (Q_1 Q_5) = Q_0 Q_1 Q_2 Q_3 Q_5 \\
p=3: \quad \mathcal{P}(6) \mathcal{P}(2) = \mathcal{P}(0,2) + \mathcal{P}(4,1) + \mathcal{P}(8)
\end{gather*}

REFERENCES:

    [Mil] J. W. Milnor, "The Steenrod algebra and its dual, Ann. of Math.
          (2) \textbf{67} (1958), 150--171.
"""

#*****************************************************************************
#       Copyright (C) 2008 John H. Palmieri <palmieri@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#*****************************************************************************

def milnor_multiplication_odd(m1,m2,p):
    r"""
    Product of Milnor basis elements defined by m1 and m2.

    INPUT:
        m1 -- pair of tuples (e,r), where e is an increasing tuple of
            non-negative integers and r is a tuple of non-negative integers
        m2 -- pair of tuples (f,s), same format as m1
        p -- odd prime number

    OUTPUT:
        Dictionary of terms of the form (tuple: coeff), where 'tuple' is a
        pair of tuples, as for r and s, and 'coeff' is an integer mod p.

    This computes the product of the Milnor basis elements
    $Q_e1 Q_e2 ... P(r_1, r_2, ...)$ and
    $Q_f1 Q_f2 ... P(s_1, s_2, ...)$.

    EXAMPLES:
        sage: from sage.algebras.steenrod_milnor_multiplication_odd import milnor_multiplication_odd
        sage: milnor_multiplication_odd(((0,2),(5,)), ((1,),(1,)), 5)
        {((0, 1, 2), (0, 1)): 4, ((0, 1, 2), (6,)): 4}
        sage: milnor_multiplication_odd(((0,2,4),()), ((1,3),()), 7)
        {((0, 1, 2, 3, 4), ()): 6}
        sage: milnor_multiplication_odd(((0,2,4),()), ((1,5),()), 7)
        {((0, 1, 2, 4, 5), ()): 1}
        sage: milnor_multiplication_odd(((),(6,)), ((),(2,)), 3)
        {((), (4, 1)): 1, ((), (8,)): 1, ((), (0, 2)): 1}

    This uses the same algorithm Monks does in his Maple package to
    iterate through the possible matrices.
    """
    from sage.rings.all import GF
    F = GF(p)
    (f,s) = m2
    # First compute Q_e0 Q_e1 ... P(r1, r2, ...) Q_f0 Q_f1 ...
    # Store results (as dictionary of pairs of tuples) in 'answer'.
    answer = {m1: F(1)}
    for k in f:
        old_answer = answer
        answer = {}
        for mono in old_answer:
            if k not in mono[0]:
                q_mono = set(mono[0])
                if len(q_mono) > 0:
                    ind = len(q_mono.intersection(range(k,1+max(q_mono))))
                else:
                    ind = 0
                coeff = (-1)**ind * old_answer[mono]
                lst = list(mono[0])
                if ind == 0:
                    lst.append(k)
                else:
                    lst.insert(-ind,k)
                q_mono = tuple(lst)
                p_mono = mono[1]
                answer[(q_mono, p_mono)] = F(coeff)
            for i in range(1,1+len(mono[1])):
                if (k+i not in mono[0]) and (p**k <= mono[1][i-1]):
                    q_mono = set(mono[0])
                    if len(q_mono) > 0:
                        ind = len(q_mono.intersection(range(k+i,1+max(q_mono))))
                    else:
                        ind = 0
                    coeff = (-1)**ind
                    lst = list(mono[0])
                    if ind == 0:
                        lst.append(k+i)
                    else:
                        lst.insert(-ind,k+i)
                    q_mono = tuple(lst)
                    p_mono = list(mono[1])
                    p_mono[i-1] = p_mono[i-1] - p**k
                    answer[(q_mono, tuple(p_mono))] = F(coeff)
    # Now for the Milnor matrices.  For each entry '(e,r): coeff' in answer,
    # multiply r with s.  Record coefficient for matrix and multiply by coeff.
    # Store in 'result'.
    if len(s) == 0:
        result = answer
    else:
        result = {}
        for (e, r) in answer:
            old_coeff = answer[(e,r)]
            # Milnor multiplication for r and s
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
                coeff = old_coeff
                diagonal = [0]*diags
                while n <= diags and coeff != 0:
                    nth_diagonal = [M[i][n-i] for i in range(max(0,n-cols+1), min(1+n,rows))]
                    coeff = coeff * multinomial_odd(nth_diagonal,p)
                    diagonal[n-1] = sum(nth_diagonal)
                    n = n + 1
                if coeff != 0:
                    i = diags - 1
                    while i >= 0 and diagonal[i] == 0:
                        i = i - 1
                    t = tuple(diagonal[:i+1])
                    if result.has_key((e,t)):
                        result[(e,t)] = F(coeff + result[t])
                    else:
                        result[(e,t)] = F(coeff)
                    # now look for new matrices:
                found = False
                i = 1
                while not found and i < rows:
                    temp_sum = M[i][0]
                    j = 1
                    while not found and j < cols:
                        # check to see if column index j is small enough
                        if temp_sum >= p**j:
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
                                M[i][0] = temp_sum - p**j
                            else:
                                temp_sum += M[i][j] * p**j
                        else:
                            temp_sum += M[i][j] * p**j
                        j = j + 1
                    i = i + 1
    return result


def multinomial_odd(list,p):
    """
    Multinomial coefficient of list, mod p.

    INPUT:
        list -- list of integers
        p -- a prime number

    OUTPUT:
        Associated multinomial coefficient, mod p

    Given the input $[n_1, n_2, n_3, ...]$, this computes the
    multinomial coefficient $(n_1 + n_2 + n_3 + ...)! / (n_1! n_2!
    n_3! ...)$, mod $p$.  The method is this: expand each $n_i$ in
    base $p$: $n_i = \sum_j p^j n_{ij}$.  Do the same for the sum of
    the $n_i$'s, which we call $m$: $m = \sum_j p^j m_j$.  Then the
    multinomial coefficient is congruent, mod $p$, to the product of
    the multinomial coefficients $m_j! / (n_{1j}! n_{2j}! ...)$.

    Furthermore, any multinomial coefficient $m! / (n_1! n_2! ...)$
    can be computed as a product of binomial coefficients: it equals
    \[
    \binom{n_1}{n_1} \binom{n_1 + n_2}{n_2} \binom{n_1 + n_2 + n_3}{n_3} ...
    \]
    This is convenient because Sage's binomial function returns
    integers, not rational numbers (as would be produced just by
    dividing factorials).

    EXAMPLES:
        sage: from sage.algebras.steenrod_milnor_multiplication_odd import multinomial_odd
        sage: multinomial_odd([1,2,4], 2)
        1
        sage: multinomial_odd([1,2,4], 7)
        0
        sage: multinomial_odd([1,2,4], 11)
        6
        sage: multinomial_odd([1,2,4], 101)
        4
        sage: multinomial_odd([1,2,4], 107)
        105
    """
    from sage.rings.arith import factorial
    from sage.rings.all import GF
    from sage.rings.arith import binomial
    from sage.algebras.steenrod_algebra_element import base_p_expansion
    n = sum(list)
    answer = 1
    F = GF(p)
    n_expansion = base_p_expansion(n,p)
    list_expansion = [base_p_expansion(k,p) for k in list]
    index = 0
    while answer != 0 and index < len(n_expansion):
        multi = F(1)
        partial_sum = 0
        for exp in list_expansion:
            if index < len(exp):
                partial_sum = partial_sum + exp[index]
                multi = F(multi * binomial(partial_sum, exp[index]))
        answer = F(answer * multi)
        index += 1
    return answer
