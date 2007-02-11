r"""
Combinatorial Functions.

AUTHORS:
        -- David Joyner (2006-07), initial implementation.
        -- William Stein (2006-07), editing of docs and code; many optimizations,
                      refinements, and bug fixes in corner cases
        -- DJ (2006-09): bug fix for combinations, added permutations_iterator,
                      combinations_iterator from Python Cookbook, edited docs.

This module implements some combinatorial functions, as listed
below. For a more detailed description, see the relevant docstrings.

Sequences:
\begin{itemize}
\item Bell numbers, \code{bell_number}

\item Bernoulli numbers, \code{bernoulli_number} (though PARI's bernoulli is
  better)

\item Catalan numbers, \code{catalan_number} (not to be confused with the
  Catalan constant)

\item Eulerian/Euler numbers, \code{euler_number} (Maxima)

\item Fibonacci numbers, \code{fibonacci} (PARI) and \code{fibonacci_number} (GAP)
  The PARI version is better.

\item Lucas numbers, \code{lucas_number1}, \code{lucas_number2}.

\item Stirling numbers, \code{stirling_number1}, \code{stirling_number2}.
\end{itemize}

Set-theoretic constructions:
\begin{itemize}
\item Combinations of a multiset, \code{combinations}, \code{combinations_iterator},
and \code{number_of_combinations}. These are unordered selections without
repetition of k objects from a multiset S.

\item Arrangements of a multiset, \code{arrangements} and \code{number_of_arrangements}
  These are ordered selections without repetition of k objects from a
  multiset S.

\item Derangements of a multiset, \code{derangements} and \code{number_of_derangements}.

\item Tuples of a multiset, \code{tuples} and \code{number_of_tuples}.
  An ordered tuple of length k of set S is a ordered selection with
  repetitions of S and is represented by a sorted list of length k
  containing elements from S.

\item Unordered tuples of a set, \code{unordered_tuple} and \code{number_of_unordered_tuples}.
  An unordered tuple of length k of set S is a unordered selection with
  repetitions of S and is represented by a sorted list of length k
  containing elements from S.

\item Permutations of a multiset, \code{permutations}, \code{permutations_iterator},
\code{number_of_permutations}. A permutation is a list that contains exactly the same elements but
possibly in different order.
\end{itemize}

Partitions:
\begin{itemize}
\item Partitions of a set, \code{partitions_set}, \code{number_of_partitions_set}.
  An unordered partition of set S is a set of pairwise disjoint
  nonempty sets with union S and is represented by a sorted list of
  such sets.

\item Partitions of an integer, \code{partitions_list}, \code{number_of_partitions_list}.
  An unordered partition of n is an unordered sum
  $n = p_1+p_2 +\ldots+ p_k$ of positive integers and is represented by
  the list $p = [p_1,p_2,\ldots,p_k]$, in nonincreasing order, i.e.,
  $p1\geq p_2 ...\geq p_k$.

\item Ordered partitions of an integer, \code{ordered_partitions},
  \code{number_of_ordered_partitions}.
  An ordered partition of n is an ordered sum $n = p_1+p_2 +\ldots+ p_k$
  of positive integers and is represented by
  the list $p = [p_1,p_2,\ldots,p_k]$, in nonincreasing order, i.e.,
  $p1\geq p_2 ...\geq p_k$.

\item Restricted partitions of an integer, \code{partitions_restricted},
  \code{number_of_partitions_restricted}.
  An unordered restricted partition of n is an unordered sum
  $n = p_1+p_2 +\ldots+ p_k$ of positive integers $p_i$ belonging to a
  given set $S$, and is represented by the list $p = [p_1,p_2,\ldots,p_k]$,
  in nonincreasing order, i.e., $p1\geq p_2 ...\geq p_k$.

\item \code{partitions_greatest}
   implements a special type of restricted partition.

\item \code{partitions_greatest_eq} is another type of restricted partition.

\item Tuples of partitions, \code{partition_tuples},
  \code{number_of_partition_tuples}.
  A $k$-tuple of partitions is represented by a list of all $k$-tuples
  of partitions which together form a partition of $n$.

\item Powers of a partition, \code{partition_power(pi, k)}.
  The power of a partition corresponds to the $k$-th power of a
  permutation with cycle structure $\pi$.

\item Sign of a partition, \code{partition_sign( pi ) }
  This means the sign of a permutation with cycle structure given by the
  partition pi.

\item Associated partition, \code{partition_associated( pi )}
  The ``associated'' (also called ``conjugate'' in the literature)
  partition of the partition pi which is obtained by transposing the
  corresponding Ferrers diagram.

\item Ferrers diagram, \code{ferrers_diagram}.
  Analogous to the Young diagram of an irredicible representation
  of $S_n$.
  \end{itemize}

Related functions:

\begin{itemize}
\item Bernoulli polynomials, \code{bernoulli_polynomial}
\end{itemize}

Implemented in other modules (listed for completeness):

The \module{sage.rings.arith} module contains
the following combinatorial functions:
\begin{itemize}
 \item binomial
     the binomial coefficient (wrapped from PARI)
 \item factorial (wrapped from PARI)
 \item partition (from the Python Cookbook)
     Generator of the list of all the partitions of the integer $n$.
 \item \code{number_of_partitions} (wrapped from PARI)
     the *number* of partitions:
 \item \code{falling_factorial}
     Definition: for integer $a \ge 0$ we have $x(x-1) \cdots (x-a+1)$.
     In all other cases we use the GAMMA-function:
     $\frac {\Gamma(x+1)} {\Gamma(x-a+1)}$.
 \item \code{rising_factorial}
     Definition: for integer $a \ge 0$ we have $x(x+1) \cdots (x+a-1)$.
     In all other cases we use the GAMMA-function:
     $\frac {\Gamma(x+a)} {\Gamma(x)}$.
 \item gaussian_binomial
     the gaussian binomial
     $$
        \binom{n}{k}_q = \frac{(1-q^m)(1-q^{m-1})\cdots (1-q^{m-r+1})}
                             {(1-q)(1-q^2)\cdots (1-q^r)}.
     $$
\end{itemize}

The \module{sage.groups.perm_gps.permgroup_elements} contains the following
combinatorial functions:
\begin{itemize}
\item matrix method of PermutationGroupElement yielding the permutation
     matrix of the group element.
\end{itemize}

\begin{verbatim}
TODO:
   GUAVA commands:
    * MOLS returns a list of n Mutually Orthogonal Latin Squares (MOLS).
    * HadamardMat returns a Hadamard matrix of order n.
    * VandermondeMat
    * GrayMat returns a list of all different vectors of length n over
      the field F, using Gray ordering.
   Not in GAP:
    * Rencontres numbers
      http://en.wikipedia.org/wiki/Rencontres_number
\end{verbatim}

REFERENCES:
      \url{http://en.wikipedia.org/wiki/Twelvefold_way} (general reference)

"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner <wdjoyner@gmail.com>,
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.interfaces.all import gap, maxima
from sage.rings.all import QQ, RR, ZZ
from sage.rings.arith import binomial
from sage.misc.sage_eval import sage_eval
from sage.libs.all import pari

######### combinatorial sequences

def bell_number(n):
    r"""
    Returns the n-th Bell number (the number of ways to partition a
    set of n elements into pairwise disjoint nonempty subsets).

    If $n \leq 0$, returns $1$.

    Wraps GAP's Bell.

    EXAMPLES:
        sage: bell_number(10)
        115975
        sage: bell_number(2)
        2
        sage: bell_number(-10)
        1
        sage: bell_number(1)
        1
        sage: bell_number(1/3)
        Traceback (most recent call last):
        ...
        TypeError: no coercion of this rational to integer
    """
    ans=gap.eval("Bell(%s)"%ZZ(n))
    return eval(ans)

## def bernoulli_number(n):
##     r"""
##     Returns the n-th Bernoulli number $B_n$; $B_n/n!$ is the
##     coefficient of $x^n$ in the power series of $x/(e^x-1)$.
##     Wraps GAP's Bernoulli.
##     EXAMPLES:
##         sage: bernoulli_number(50)
##         495057205241079648212477525/66
##     """
##     ans=gap.eval("Bernoulli(%s)"%n)
##     return QQ(ans)   ## use QQ, not eval, here

def catalan_number(n):
    r"""
    Returns the n-th Catalan number

      Catalan numbers: The $n$-th Catalan number is given directly in terms of
      binomial coefficients by
      \[
      C_n = \frac{1}{n+1}{2n\choose n} = \frac{(2n)!}{(n+1)!\,n!}
            \qquad\mbox{ for }n\ge 0.
      \]

      Consider the set $S = \{ 1, ..., n \}$. A {\it noncrossing partition} of $S$
      is a partition in which no two blocks "cross" each other, i.e., if a
      and b belong to one block and x and y to another, they are not arranged
      in the order axby. $C_n$ is the number of noncrossing partitions of the set $S$.
      There are many other interpretations (see REFERENCES).

    When $n=-1$, this function raises a ZeroDivisionError; for other $n<0$ it
    returns $0$.

    EXAMPLES:
        sage: [catalan_number(i) for i in range(7)]
        [1, 1, 2, 5, 14, 42, 132]
        sage: maxima.eval("-(1/2)*taylor (sqrt (1-4*x^2), x, 0, 15)")
        '-1/2 + x^2 + x^4 + 2*x^6 + 5*x^8 + 14*x^10 + 42*x^12 + 132*x^14'
        sage: [catalan_number(i) for i in range(-7,7) if i != -1]
        [0, 0, 0, 0, 0, 0, 1, 1, 2, 5, 14, 42, 132]
        sage: catalan_number(-1)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Rational division by zero

    REFERENCES:
    \begin{itemize}
      \item \url{http://en.wikipedia.org/wiki/Catalan_number}
      \item \url{http://www-history.mcs.st-andrews.ac.uk/~history/Miscellaneous/CatalanNumbers/catalan.html}
    \end{itemize}

    """
    n = ZZ(n)
    return binomial(2*n,n)/(n+1)

def euler_number(n):
    """
    Returns the n-th Euler number

    IMPLEMENTATION: Wraps Maxima's euler.

    EXAMPLES:
        sage: [euler_number(i) for i in range(10)]
        [1, 0, -1, 0, 5, 0, -61, 0, 1385, 0]
        sage: maxima.eval("taylor (2/(exp(x)+exp(-x)), x, 0, 10)")
        '1 - x^2/2 + 5*x^4/24 - 61*x^6/720 + 277*x^8/8064 - 50521*x^10/3628800'
        sage: [euler_number(i)/factorial(i) for i in range(11)]
        [1, 0, -1/2, 0, 5/24, 0, -61/720, 0, 277/8064, 0, -50521/3628800]
        sage: euler_number(-1)
        Traceback (most recent call last):
        ...
        ValueError: n (=-1) must be a nonnegative integer

    REFERENCES:
        http://en.wikipedia.org/wiki/Euler_number
    """
    n = ZZ(n)
    if n < 0:
        raise ValueError, "n (=%s) must be a nonnegative integer"%n
    return ZZ(maxima.eval("euler(%s)"%n))

def fibonacci(n, algorithm="pari"):
    """
    Returns then n-th Fibonacci number. The Fibonacci sequence $F_n$
    is defined by the initial conditions $F_1=F_2=1$ and the
    recurrence relation $F_{n+2} = F_{n+1} + F_n$. For negative $n$ we
    define $F_n = (-1)^{n+1}F_{-n}$, which is consistent with the
    recurrence relation.

    For negative $n$, define $F_{n} = -F_{|n|}$.

    INPUT:
        algorithm -- string:
                     "pari" -- (default) -- use the PARI C library's fibo function.
                     "gap" -- use GAP's Fibonacci function

    NOTE: PARI is tens to hundreds of times faster than GAP here;
          moreover, PARI works for every large input whereas GAP
          doesn't.

    EXAMPLES:
        sage: fibonacci(10)
        55
        sage: fibonacci(10, algorithm='gap')
        55

        sage: fibonacci(-100)
        -354224848179261915075
        sage: fibonacci(100)
        354224848179261915075

        sage: fibonacci(0)
        0
        sage: fibonacci(1/2)
        Traceback (most recent call last):
        ...
        TypeError: no coercion of this rational to integer
    """
    n = ZZ(n)
    if algorithm == 'pari':
        return ZZ(pari(n).fibonacci())
    elif algorithm == 'gap':
        return ZZ(gap.eval("Fibonacci(%s)"%n))
    else:
        raise ValueError, "no algorithm %s"%algorithm

def lucas_number1(n,P,Q):
    """
    Returns then n-th Lucas number ``of the first kind'' (this is not
    standard terminology). The Lucas sequence $L^{(1)}_n$ is defined
    by the initial conditions $L^{(1)}_1=0$, $L^{(1)}_2=1$ and the recurrence
    relation $L^{(1)}_{n+2} = P*L^{(1)}_{n+1} - Q*L^{(1)}_n$.

    Wraps GAP's Lucas(...)[1].

    P=1, Q=-1 fives the Fibonacci sequence.

    INPUT:
        n -- integer
        P, Q -- integer or rational numbers

    OUTPUT:
        integer or rational number

    EXAMPLES:
        sage: lucas_number1(5,1,-1)
        5
        sage: lucas_number1(6,1,-1)
        8
        sage: lucas_number1(7,1,-1)
        13
        sage: lucas_number1(7,1,-2)
        43

        sage: lucas_number1(5,2,3/5)
        229/25
        sage: lucas_number1(5,2,1.5)
        Traceback (most recent call last):
        ...
        TypeError: no implicit coercion of element to the rational numbers

    There was a conjecture that the sequence $L_n$ defined by
    $L_{n+2} = L_{n+1} + L_n$, $L_1=1$, $L_2=3$, has the property
    that $n$ prime implies that $L_n$ is prime.

        sage: lucas = lambda n:(5/2)*lucas_number1(n,1,-1)+(1/2)*lucas_number2(n,1,-1)
        sage: [[lucas(n),is_prime(lucas(n)),n+1,is_prime(n+1)] for n in range(15)]
        [[1, False, 1, False],
         [3, True, 2, True],
         [4, False, 3, True],
         [7, True, 4, False],
         [11, True, 5, True],
         [18, False, 6, False],
         [29, True, 7, True],
         [47, True, 8, False],
         [76, False, 9, False],
         [123, False, 10, False],
         [199, True, 11, True],
         [322, False, 12, False],
         [521, True, 13, True],
         [843, False, 14, False],
         [1364, False, 15, False]]

    Can you use SAGE to find a counterexample to the conjecture?
    """
    ans=gap.eval("Lucas(%s,%s,%s)[1]"%(QQ._coerce_(P),QQ._coerce_(Q),ZZ(n)))
    return sage_eval(ans)

def lucas_number2(n,P,Q):
    r"""
    Returns then n-th Lucas number ``of the second kind'' (this is not
    standard terminology). The Lucas sequence $L^{(2)}_n$ is defined
    by the initial conditions $L^{(2)}_1=2$, $L^{(2)}_2=P$ and the recurrence
    relation $L^{(2)}_{n+2} = P*L^{(2)}_{n+1} - Q*L^{(2)}_n$.

    Wraps GAP's Lucas(...)[2].

    INPUT:
        n -- integer
        P, Q -- integer or rational numbers

    OUTPUT:
        integer or rational number

    EXAMPLES:
        sage: [lucas_number2(i,1,-1) for i in range(10)]
        [2, 1, 3, 4, 7, 11, 18, 29, 47, 76]
        sage: [fibonacci(i-1)+fibonacci(i+1) for i in range(10)]
        [2, 1, 3, 4, 7, 11, 18, 29, 47, 76]

        sage: n = lucas_number2(5,2,3); n
        2
        sage: type(n)
        <type 'sage.rings.integer.Integer'>
        sage: n = lucas_number2(5,2,-3/9); n
        418/9
        sage: type(n)
        <type 'sage.rings.rational.Rational'>

    The case P=1, Q=-1 is the Lucas sequence in Brualdi's
    {\bf Introductory Combinatorics}, 4th ed., Prentice-Hall, 2004:
        sage: [lucas_number2(n,1,-1) for n in range(10)]
        [2, 1, 3, 4, 7, 11, 18, 29, 47, 76]
    """
    ans=gap.eval("Lucas(%s,%s,%s)[2]"%(QQ._coerce_(P),QQ._coerce_(Q),ZZ(n)))
    return sage_eval(ans)

def stirling_number1(n,k):
    """
    Returns the n-th Stilling number $S_1(n,k)$ of the first kind (the number of
    permutations of n points with k cycles).
    Wraps GAP's Stirling1.

    EXAMPLES:
        sage: stirling_number1(3,2)
        3
        sage: stirling_number1(5,2)
        50
        sage: 9*stirling_number1(9,5)+stirling_number1(9,4)
        269325
        sage: stirling_number1(10,5)
        269325

    Indeed, $S_1(n,k) = S_1(n-1,k-1) + (n-1)S_1(n-1,k)$.
    """
    return ZZ(gap.eval("Stirling1(%s,%s)"%(ZZ(n),ZZ(k))))

def stirling_number2(n,k):
    """
    Returns the n-th Stirling number $S_2(n,k)$ of the second kind (the
    number of ways to partition a set of n elements into k pairwise
    disjoint nonempty subsets). (The n-th Bell number is the sum of
    the $S_2(n,k)$'s, $k=0,...,n$.)
    Wraps GAP's Stirling2.

    EXAMPLES:
    Stirling numbers satisfy $S_2(n,k) = S_2(n-1,k-1) + kS_2(n-1,k)$:

        sage: 5*stirling_number2(9,5) + stirling_number2(9,4)
        42525
        sage: stirling_number2(10,5)
        42525

        sage: n = stirling_number2(20,11); n
        1900842429486
        sage: type(n)
        <type 'sage.rings.integer.Integer'>
    """
    return ZZ(gap.eval("Stirling2(%s,%s)"%(ZZ(n),ZZ(k))))

#### combinatorial sets/lists

def combinations(mset,k):
    r"""
    A {\it combination} of a multiset (a list of objects which may
    contain the same object several times) mset is an unordered
    selection without repetitions and is represented by a sorted
    sublist of mset.  Returns the set of all combinations of the
    multiset mset with k elements.

    WARNING: Wraps GAP's Combinations.  Hence mset must be a list of
    objects that have string representations that can be interpreted by
    the GAP intepreter.  If mset consists of at all complicated SAGE
    objects, this function does *not* do what you expect.  A proper
    function should be written! (TODO!)

    EXAMPLES:
        sage: mset = [1,1,2,3,4,4,5]
        sage: combinations(mset,2)
        [[1, 1],
         [1, 2],
         [1, 3],
         [1, 4],
         [1, 5],
         [2, 3],
         [2, 4],
         [2, 5],
         [3, 4],
         [3, 5],
         [4, 4],
         [4, 5]]
         sage: mset = ["d","a","v","i","d"]
         sage: combinations(mset,3)
         ['add', 'adi', 'adv', 'aiv', 'ddi', 'ddv', 'div']

    NOTE: For large lists, this raises a string error.
    """
    ans=gap.eval("Combinations(%s,%s)"%(mset,ZZ(k))).replace("\n","")
    return eval(ans)

def combinations_iterator(mset,n=None):
    """
    Posted by Raymond Hettinger, 2006/03/23, to the Python Cookbook:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/474124

    Much faster than combinations.

    EXAMPLES:
        sage: X = combinations_iterator([1,2,3,4,5],3)
        sage: [x for x in X]
        [[1, 2, 3],
         [1, 2, 4],
         [1, 2, 5],
         [1, 3, 4],
         [1, 3, 5],
         [1, 4, 5],
         [2, 3, 4],
         [2, 3, 5],
         [2, 4, 5],
         [3, 4, 5]]
    """
    items = mset
    if n is None:
        n = len(items)
    for i in range(len(items)):
        v = items[i:i+1]
        if n == 1:
            yield v
        else:
            rest = items[i+1:]
            for c in combinations_iterator(rest, n-1):
                yield v + c


def number_of_combinations(mset,k):
    """
    Returns the size of combinations(mset,k).
    IMPLEMENTATION: Wraps GAP's NrCombinations.


    NOTE: mset must be a list of integers or strings (i.e., this is very restricted).

    EXAMPLES:
        sage: mset = [1,1,2,3,4,4,5]
        sage: number_of_combinations(mset,2)
        12
    """
    return ZZ(gap.eval("NrCombinations(%s,%s)"%(mset,ZZ(k))))

def arrangements(mset,k):
    r"""
    An arrangement of mset is an ordered selection without repetitions
    and is represented by a list that contains only elements from
    mset, but maybe in a different order.

    \code{arrangements} returns the set of arrangements of the
    multiset mset that contain k elements.

    IMPLEMENTATION: Wraps GAP's Arrangements.

    WARNING: Wraps GAP -- hence mset must be a list of objects that
    have string representations that can be interpreted by the GAP
    intepreter.  If mset consists of at all complicated SAGE objects,
    this function does *not* do what you expect.  A proper function
    should be written! (TODO!)

    EXAMPLES:
        sage: mset = [1,1,2,3,4,4,5]
        sage: arrangements(mset,2)
        [[1, 1],
         [1, 2],
         [1, 3],
         [1, 4],
         [1, 5],
         [2, 1],
         [2, 3],
         [2, 4],
         [2, 5],
         [3, 1],
         [3, 2],
         [3, 4],
         [3, 5],
         [4, 1],
         [4, 2],
         [4, 3],
         [4, 4],
         [4, 5],
         [5, 1],
         [5, 2],
         [5, 3],
         [5, 4]]
         sage: arrangements( ["c","a","t"], 2 )
         ['ac', 'at', 'ca', 'ct', 'ta', 'tc']
         sage: arrangements( ["c","a","t"], 3 )
         ['act', 'atc', 'cat', 'cta', 'tac', 'tca']
    """
    ans=gap.eval("Arrangements(%s,%s)"%(mset,k))
    return eval(ans)

def number_of_arrangements(mset,k):
    """
    Returns the size of arrangements(mset,k).
    Wraps GAP's NrArrangements.

    EXAMPLES:
        sage: mset = [1,1,2,3,4,4,5]
        sage: number_of_arrangements(mset,2)
        22
    """
    return ZZ(gap.eval("NrArrangements(%s,%s)"%(mset,ZZ(k))))

def derangements(mset):
    """
    A derangement is a fixed point free permutation of list and is
    represented by a list that contains exactly the same elements as
    mset, but possibly in different order.  Derangements returns the
    set of all derangements of a multiset.

    Wraps GAP's Derangements.

    WARNING: Wraps GAP -- hence mset must be a list of objects that
    have string representations that can be interpreted by the GAP
    intepreter.  If mset consists of at all complicated SAGE objects,
    this function does *not* do what you expect.  A proper function
    should be written! (TODO!)

    EXAMPLES:
        sage: mset = [1,2,3,4]
        sage: derangements(mset)
        [[2, 1, 4, 3],
         [2, 3, 4, 1],
         [2, 4, 1, 3],
         [3, 1, 4, 2],
         [3, 4, 1, 2],
         [3, 4, 2, 1],
         [4, 1, 2, 3],
         [4, 3, 1, 2],
         [4, 3, 2, 1]]
         sage: derangements(["c","a","t"])
         ['atc', 'tca']

    """
    ans=gap.eval("Derangements(%s)"%mset)
    return eval(ans)

def number_of_derangements(mset):
    """
    Returns the size of derangements(mset).
    Wraps GAP's NrDerangements.

    EXAMPLES:
        sage: mset = [1,2,3,4]
        sage: number_of_derangements(mset)
        9
    """
    ans=gap.eval("NrDerangements(%s)"%mset)
    return ZZ(ans)

def tuples(S,k):
    """
    An ordered tuple of length k of set is an ordered selection with
    repetition and is represented by a list of length k containing
    elements of set.
    tuples returns the set of all ordered tuples of length k of the set.

    EXAMPLES:
        sage: S = [1,2]
        sage: tuples(S,3)
	[[1, 1, 1], [2, 1, 1], [1, 2, 1], [2, 2, 1], [1, 1, 2], [2, 1, 2], [1, 2, 2], [2, 2, 2]]
        sage: mset = ["s","t","e","i","n"]
        sage: tuples(mset,2)
	[['s', 's'], ['t', 's'], ['e', 's'], ['i', 's'], ['n', 's'], ['s', 't'], ['t', 't'],
	 ['e', 't'], ['i', 't'], ['n', 't'], ['s', 'e'], ['t', 'e'], ['e', 'e'], ['i', 'e'],
         ['n', 'e'], ['s', 'i'], ['t', 'i'], ['e', 'i'], ['i', 'i'], ['n', 'i'], ['s', 'n'],
	 ['t', 'n'], ['e', 'n'], ['i', 'n'], ['n', 'n']]

    The Set(...) comparisons are necessary because finite fields are not
    enumerated in a standard order.
	sage: K.<a> = GF(4, 'a')
	sage: mset = [x for x in K if x!=0]
	sage: tuples(mset,2)
        [[a, a], [a + 1, a], [1, a], [a, a + 1], [a + 1, a + 1], [1, a + 1], [a, 1], [a + 1, 1], [1, 1]]

    AUTHOR: Jon Hanke (2006-08?)
    """
    import copy
    if k<=0:
        return [[]]
    if k==1:
        return [[x] for x in S]
    ans = []
    for s in S:
        for x in tuples(S,k-1):
            y = copy.copy(x)
            y.append(s)
            ans.append(y)
    return ans
    ## code wrapping GAP's Tuples:
    #ans=gap.eval("Tuples(%s,%s)"%(S,k))
    #return eval(ans)


def number_of_tuples(S,k):
    """
    Returns the size of tuples(S,k).
    Wraps GAP's NrTuples.

    EXAMPLES:
        sage: S = [1,2,3,4,5]
        sage: number_of_tuples(S,2)
        25
        sage: S = [1,1,2,3,4,5]
        sage: number_of_tuples(S,2)
        25

    """
    ans=gap.eval("NrTuples(%s,%s)"%(S,ZZ(k)))
    return ZZ(ans)

def unordered_tuples(S,k):
    """
    An unordered tuple of length k of set is a unordered selection
    with repetitions of set and is represented by a sorted list of
    length k containing elements from set.

    unordered_tuples returns the set of all unordered tuples of length k
    of the set.
    Wraps GAP's UnorderedTuples.

    WARNING: Wraps GAP -- hence mset must be a list of objects that
    have string representations that can be interpreted by the GAP
    intepreter.  If mset consists of at all complicated SAGE objects,
    this function does *not* do what you expect.  A proper function
    should be written! (TODO!)

    EXAMPLES:
        sage: S = [1,2]
        sage: unordered_tuples(S,3)
        [[1, 1, 1], [1, 1, 2], [1, 2, 2], [2, 2, 2]]
        sage: unordered_tuples(["a","b","c"],2)
        ['aa', 'ab', 'ac', 'bb', 'bc', 'cc']

    """
    ans=gap.eval("UnorderedTuples(%s,%s)"%(S,ZZ(k)))
    return eval(ans)

def number_of_unordered_tuples(S,k):
    """
    Returns the size of unordered_tuples(S,k).
    Wraps GAP's NrUnorderedTuples.

    EXAMPLES:
        sage: S = [1,2,3,4,5]
        sage: number_of_unordered_tuples(S,2)
        15
    """
    ans=gap.eval("NrUnorderedTuples(%s,%s)"%(S,ZZ(k)))
    return ZZ(ans)

def permutations(mset):
    """
    A {\it permutation} is represented by a list that contains exactly the same
    elements as mset, but possibly in different order. If mset is a
    proper set there are $|mset| !$ such permutations. Otherwise if the
    first elements appears $k_1$ times, the second element appears $k_2$ times
    and so on, the number of permutations is $|mset|! / (k_1! k_2! \ldots)$,
    which is sometimes called a {\it multinomial coefficient}.

    permutations returns the set of all permutations of a multiset.
    Wraps GAP's PermutationsList.

    WARNING: Wraps GAP -- hence mset must be a list of objects that
    have string representations that can be interpreted by the GAP
    intepreter.  If mset consists of at all complicated SAGE objects,
    this function does *not* do what you expect.  A proper function
    should be written! (TODO!)

    EXAMPLES:
        sage: mset = [1,1,2,2,2]
        sage: permutations(mset)
        [[1, 1, 2, 2, 2],
         [1, 2, 1, 2, 2],
         [1, 2, 2, 1, 2],
         [1, 2, 2, 2, 1],
         [2, 1, 1, 2, 2],
         [2, 1, 2, 1, 2],
         [2, 1, 2, 2, 1],
         [2, 2, 1, 1, 2],
         [2, 2, 1, 2, 1],
         [2, 2, 2, 1, 1]]

    """
    ans=gap.eval("PermutationsList(%s)"%mset)
    return eval(ans)

def permutations_iterator(mset,n=None):
    """
    Posted by Raymond Hettinger, 2006/03/23, to the Python Cookbook:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/474124

    EXAMPLES:
        sage: X = permutations_iterator(range(3),2)
        sage: [x for x in X]
        [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]]
    """
    items = mset
    if n is None:
        n = len(items)
    for i in range(len(items)):
        v = items[i:i+1]
        if n == 1:
            yield v
        else:
            rest = items[:i] + items[i+1:]
            for p in permutations_iterator(rest, n-1):
                yield v + p

def number_of_permutations(mset):
    """
    Returns the size of permutations(mset).
    Wraps GAP's NrPermutationsList.

    EXAMPLES:
        sage: mset = [1,1,2,2,2]
        sage: number_of_permutations(mset)
        10

    """
    ans=gap.eval("NrPermutationsList(%s)"%mset)
    return ZZ(ans)

#### partitions

def partitions_set(S,k=None):
    r"""
    An {\it unordered partition} of a set $S$ is a set of pairwise disjoint
    nonempty subsets with union $S$ and is represented by a sorted
    list of such subsets.

    partitions_set returns the set of all unordered partitions of the
    list $S$ of increasing positive integers into k pairwise disjoint
    nonempty sets. If k is omitted then all partitions are returned.

    The Bell number $B_n$, named in honor of Eric Temple Bell, is
    the number of different partitions of a set with n elements.

    WARNING: Wraps GAP -- hence S must be a list of objects that have
    string representations that can be interpreted by the GAP
    intepreter.  If mset consists of at all complicated SAGE objects,
    this function does *not* do what you expect.  A proper function
    should be written! (TODO!)

    Wraps GAP's PartitionsSet.

    EXAMPLES:
        sage: S = [1,2,3,4]
        sage: partitions_set(S,2)
        [[[1], [2, 3, 4]],
         [[1, 2], [3, 4]],
         [[1, 2, 3], [4]],
         [[1, 2, 4], [3]],
         [[1, 3], [2, 4]],
         [[1, 3, 4], [2]],
         [[1, 4], [2, 3]]]

    REFERENCES:
       http://en.wikipedia.org/wiki/Partition_of_a_set
    """
    if k==None:
        ans=gap.eval("PartitionsSet(%s)"%S)
    else:
        ans=gap.eval("PartitionsSet(%s,%s)"%(S,k))
    return eval(ans)

def number_of_partitions_set(S,k):
    r"""
    Returns the size of \code{partitions_set(S,k)}.  Wraps GAP's
    NrPartitionsSet.

    The Stirling number of the second kind is the number of partitions
    of a set of size n into k blocks.

    EXAMPLES:
        sage: mset = [1,2,3,4]
        sage: number_of_partitions_set(mset,2)
        7
        sage: stirling_number2(4,2)
        7

    REFERENCES
        http://en.wikipedia.org/wiki/Partition_of_a_set

    """
    if k==None:
        ans=gap.eval("NrPartitionsSet(%s)"%S)
    else:
        ans=gap.eval("NrPartitionsSet(%s,%s)"%(S,ZZ(k)))
    return ZZ(ans)

def partitions_list(n,k=None):
    r"""
    An {\it unordered partition of $n$} is an unordered sum
    $n = p_1+p_2 +\ldots+ p_k$ of positive integers and is represented by
    the list $p = [p_1,p_2,\ldots,p_k]$, in nonincreasing order, i.e.,
    $p1\geq p_2 ...\geq p_k$.

    \code{partitions_list(n,k)} returns the list of all (unordered)
    partitions of the positive integer n into sums with k summands. If
    k is omitted then all partitions are returned.

    Do not call partitions_list with an n much larger than 40, in
    which case there are 37338 partitions, since the list will simply
    become too large.

    Wraps GAP's Partitions.

    The function \code{partitions} (a wrapper for the corresponding
    PARI function) returns not a list but rather a generator for a
    list. It is also a function of only one argument.

    EXAMPLES:
        sage: partitions_list(10,2)
        [[5, 5], [6, 4], [7, 3], [8, 2], [9, 1]]
        sage: partitions_list(5)
        [[1, 1, 1, 1, 1], [2, 1, 1, 1], [2, 2, 1], [3, 1, 1], [3, 2], [4, 1], [5]]

    However, partitions(5) returns ``<generator object at ...>''.
    """
    if k==None:
        ans=gap.eval("Partitions(%s)"%(n))
    else:
        ans=gap.eval("Partitions(%s,%s)"%(n,k))
    return eval(ans.replace('\n',''))

def number_of_partitions_list(n,k=None):
    r"""
    Returns the size of partitions_list(n,k).

    Wraps GAP's NrPartitions.

    It is possible to associate with every partition of the integer n
    a conjugacy class of permutations in the symmetric group on n
    points and vice versa.  Therefore p(n) = NrPartitions(n) is the
    number of conjugacy classes of the symmetric group on n points.

    \code{number_of_partitions(n)} is also available in PARI, however
    the speed seems the same until $n$ is in the thousands (in which
    case PARI is faster).

    EXAMPLES:
        sage: number_of_partitions_list(10,2)
        5
        sage: number_of_partitions_list(10)
        42

    A generating function for p(n) is given by the reciprocal of Euler's function:
    \[
    \sum_{n=0}^\infty p(n)x^n = \prod_{k=1}^\infty \left(\frac {1}{1-x^k} \right).
    \]
    SAGE verifies that the first several coefficients do instead agree:

        sage: q = PowerSeriesRing(QQ, 'q', default_prec=9).gen()
        sage: prod([(1-q^k)^(-1) for k in range(1,9)])  ## partial product of
        1 + q + 2*q^2 + 3*q^3 + 5*q^4 + 7*q^5 + 11*q^6 + 15*q^7 + 22*q^8 + O(q^9)
        sage: [number_of_partitions_list(k) for k in range(2,10)]
        [2, 3, 5, 7, 11, 15, 22, 30]

    REFERENCES:
        http://en.wikipedia.org/wiki/Partition_%28number_theory%29

    """
    if k==None:
        ans=gap.eval("NrPartitions(%s)"%(ZZ(n)))
    else:
        ans=gap.eval("NrPartitions(%s,%s)"%(ZZ(n),ZZ(k)))
    return ZZ(ans)

def ferrers_diagram(pi):
    """
    Return the Ferrers diagram of pi.

    INPUT:
        pi -- a partition, given as a list of integers.

    EXAMPLES:
        sage: print ferrers_diagram([5,5,2,1])
        *****
        *****
        **
        *
        sage: pi = partitions_list(10)[30] ## [6,1,1,1,1]
        sage: print ferrers_diagram(pi)
        ******
        *
        *
        *
        *
        sage: pi = partitions_list(10)[33] ## [6, 3, 1]
        sage: print ferrers_diagram(pi)
        ******
        ***
        *
    """
    return '\n'.join(['*'*p for p in pi])


def ordered_partitions(n,k=None):
    r"""
    An {\it ordered partition of $n$} is an ordered sum
    $$
       n = p_1+p_2 + \cdots + p_k
    $$
    of positive integers and is represented by the list $p = [p_1,p_2,\cdots ,p_k]$.
    If $k$ is omitted then all ordered partitions are returned.

    \code{ordered_partitions(n,k)} returns the set of all (ordered)
    partitions of the positive integer n into sums with k summands.

    Do not call \code{ordered_partitions} with an n much larger than
    15, since the list will simply become too large.

    Wraps GAP's OrderedPartitions.

    The number of ordered partitions $T_n$ of $\{ 1, 2, ..., n \}$ has the
    generating function is
    \[
    \sum_n {T_n \over n!} x^n = {1 \over 2-e^x}.
    \]

    EXAMPLES:
        sage: ordered_partitions(10,2)
        [[1, 9], [2, 8], [3, 7], [4, 6], [5, 5], [6, 4], [7, 3], [8, 2], [9, 1]]

        sage: ordered_partitions(4)
        [[1, 1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 3], [2, 1, 1], [2, 2], [3, 1], [4]]

    REFERENCES:
        http://en.wikipedia.org/wiki/Ordered_partition_of_a_set

    """
    if k==None:
        ans=gap.eval("OrderedPartitions(%s)"%(ZZ(n)))
    else:
        ans=gap.eval("OrderedPartitions(%s,%s)"%(ZZ(n),ZZ(k)))
    return eval(ans.replace('\n',''))

def number_of_ordered_partitions(n,k=None):
    """
    Returns the size of ordered_partitions(n,k).
    Wraps GAP's NrOrderedPartitions.

    It is possible to associate with every partition of the integer n a conjugacy
    class of permutations in the symmetric group on n points and vice versa.
    Therefore p(n) = NrPartitions(n) is the number of conjugacy classes of the
    symmetric group on n points.


    EXAMPLES:
        sage: number_of_ordered_partitions(10,2)
        9
        sage: number_of_ordered_partitions(15)
        16384
    """
    if k==None:
        ans=gap.eval("NrOrderedPartitions(%s)"%(n))
    else:
        ans=gap.eval("NrOrderedPartitions(%s,%s)"%(n,k))
    return ZZ(ans)

def partitions_greatest(n,k):
    """
    Returns the set of all (unordered) ``restricted'' partitions of the integer n having
    parts less than or equal to the integer k.

    Wraps GAP's PartitionsGreatestLE.

    EXAMPLES:
        sage: partitions_greatest(10,2)
        [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         [2, 1, 1, 1, 1, 1, 1, 1, 1],
         [2, 2, 1, 1, 1, 1, 1, 1],
         [2, 2, 2, 1, 1, 1, 1],
         [2, 2, 2, 2, 1, 1],
         [2, 2, 2, 2, 2]]
    """
    return eval(gap.eval("PartitionsGreatestLE(%s,%s)"%(ZZ(n),ZZ(k))))

def partitions_greatest_eq(n,k):
    """
    Returns the set of all (unordered) ``restricted'' partitions of the
    integer n having at least one part equal to the integer k.

    Wraps GAP's PartitionsGreatestEQ.

    EXAMPLES:
        sage: partitions_greatest_eq(10,2)
        [[2, 1, 1, 1, 1, 1, 1, 1, 1],
         [2, 2, 1, 1, 1, 1, 1, 1],
         [2, 2, 2, 1, 1, 1, 1],
         [2, 2, 2, 2, 1, 1],
         [2, 2, 2, 2, 2]]

    """
    ans = gap.eval("PartitionsGreatestEQ(%s,%s)"%(n,k))
    return eval(ans)

def partitions_restricted(n,S,k=None):
    r"""
    A {\it restricted partition} is, like an ordinary partition, an
    unordered sum $n = p_1+p_2+\ldots+p_k$ of positive integers and is
    represented by the list $p = [p_1,p_2,\ldots,p_k]$, in nonincreasing
    order. The difference is that here the $p_i$ must be elements
    from the set $S$, while for ordinary partitions they may be
    elements from $[1..n]$.

    Returns the set of all restricted partitions of the positive integer
    n into sums with k summands with the summands of the partition coming
    from the set $S$. If k is not given all restricted partitions for all
    k are returned.

    Wraps GAP's RestrictedPartitions.

    EXAMPLES:
        sage: partitions_restricted(8,[1,3,5,7])
        [[1, 1, 1, 1, 1, 1, 1, 1],
         [3, 1, 1, 1, 1, 1],
         [3, 3, 1, 1],
         [5, 1, 1, 1],
         [5, 3],
         [7, 1]]
        sage: partitions_restricted(8,[1,3,5,7],2)
        [[5, 3], [7, 1]]
    """
    if k==None:
        ans=gap.eval("RestrictedPartitions(%s,%s)"%(n,S))
    else:
        ans=gap.eval("RestrictedPartitions(%s,%s,%s)"%(n,S,k))
    return eval(ans)

def number_of_partitions_restricted(n,S,k=None):
    """
    Returns the size of partitions_restricted(n,S,k).
    Wraps GAP's NrRestrictedPartitions.

    EXAMPLES:
        sage: number_of_partitions_restricted(8,[1,3,5,7])
        6
        sage: number_of_partitions_restricted(8,[1,3,5,7],2)
        2

    """
    if k==None:
        ans=gap.eval("NrRestrictedPartitions(%s,%s)"%(ZZ(n),S))
    else:
        ans=gap.eval("NrRestrictedPartitions(%s,%s,%s)"%(ZZ(n),S,ZZ(k)))
    return ZZ(ans)

def partitions_tuples(n,k):
    """
    partition_tuples( n, k ) returns the list of all k-tuples of partitions
    which together form a partition of n.

    k-tuples of partitions describe the classes and the characters of
    wreath products of groups with k conjugacy classes with the symmetric
    group $S_n$.

    Wraps GAP's PartitionTuples.

    EXAMPLES:
        sage: partitions_tuples(3,2)
        [[[1, 1, 1], []],
         [[1, 1], [1]],
         [[1], [1, 1]],
         [[], [1, 1, 1]],
         [[2, 1], []],
         [[1], [2]],
         [[2], [1]],
         [[], [2, 1]],
         [[3], []],
         [[], [3]]]
    """
    ans=gap.eval("PartitionTuples(%s,%s)"%(ZZ(n),ZZ(k)))
    return eval(ans)

def number_of_partitions_tuples(n,k):
    r"""
    number_of_partition_tuples( n, k ) returns the number of partition_tuples(n,k).

    Wraps GAP's NrPartitionTuples.

    EXAMPLES:
        sage: number_of_partitions_tuples(3,2)
        10
        sage: number_of_partitions_tuples(8,2)
        185

    Now we compare that with the result of the following GAP
    computation:
 \begin{verbatim}
        gap> S8:=Group((1,2,3,4,5,6,7,8),(1,2));
        Group([ (1,2,3,4,5,6,7,8), (1,2) ])
        gap> C2:=Group((1,2));
        Group([ (1,2) ])
        gap> W:=WreathProduct(C2,S8);
        <permutation group of size 10321920 with 10 generators>
        gap> Size(W);
        10321920     ## = 2^8*Factorial(8), which is good:-)
        gap> Size(ConjugacyClasses(W));
        185
\end{verbatim}
    """
    ans=gap.eval("NrPartitionTuples(%s,%s)"%(ZZ(n),ZZ(k)))
    return ZZ(ans)

def partition_power(pi,k):
    """
    partition_power( pi, k ) returns the partition corresponding to the
    $k$-th power of a permutation with cycle structure pi
    (thus describes the powermap of symmetric groups).

    Wraps GAP's PowerPartition.

    EXAMPLES:
        sage: partition_power([5,3],1)
        [5, 3]
        sage: partition_power([5,3],2)
        [5, 3]
        sage: partition_power([5,3],3)
        [5, 1, 1, 1]
        sage: partition_power([5,3],4)
        [5, 3]

     Now let us compare this to the power map on $S_8$:

        sage: G = SymmetricGroup(8)
        sage: g = G([(1,2,3,4,5),(6,7,8)])
        sage: g
        (1,2,3,4,5)(6,7,8)
        sage: g^2
        (1,3,5,2,4)(6,8,7)
        sage: g^3
        (1,4,2,5,3)
        sage: g^4
        (1,5,4,3,2)(6,7,8)

    """
    ans=gap.eval("PowerPartition(%s,%s)"%(pi,ZZ(k)))
    return eval(ans)

def partition_sign(pi):
    r"""
    partition_sign( pi ) returns the sign of a permutation with cycle structure
    given by the partition pi.

    This function corresponds to a homomorphism from the symmetric group
    $S_n$ into the cyclic group of order 2, whose kernel is exactly the
    alternating group $A_n$. Partitions of sign $1$ are called {\it even partitions}
    while partitions of sign $-1$ are called {\it odd}.

    Wraps GAP's SignPartition.

    EXAMPLES:
        sage: partition_sign([5,3])
        1
        sage: partition_sign([5,2])
        -1

    {\it Zolotarev's lemma} states that the Legendre symbol
    $ \left(\frac{a}{p}\right)$ for an integer $a \pmod p$ ($p$ a prime number),
    can be computed as sign(p_a), where sign denotes the sign of a permutation
    and p_a the permutation of the residue classes $\pmod p$ induced by
    modular multiplication by $a$, provided $p$ does not divide $a$.

    We verify this in some examples.

        sage: F = GF(11)
        sage: a = F.multiplicative_generator();a
        2
        sage: plist = [int(a*F(x)) for x in range(1,11)]; plist
        [2, 4, 6, 8, 10, 1, 3, 5, 7, 9]

    This corresponds ot the permutation (1, 2, 4, 8, 5, 10, 9, 7, 3, 6)
    (acting the set $\{1,2,...,10\}$) and to the partition [10].

        sage: p = PermutationGroupElement('(1, 2, 4, 8, 5, 10, 9, 7, 3, 6)')
        sage: p.sign()
        -1
        sage: partition_sign([10])
        -1
        sage: kronecker_symbol(11,2)
        -1

    Now replace $2$ by $3$:

        sage: plist = [int(F(3*x)) for x in range(1,11)]; plist
        [3, 6, 9, 1, 4, 7, 10, 2, 5, 8]
        sage: range(1,11)
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        sage: p = PermutationGroupElement('(3,4,8,7,9)')
        sage: p.sign()
        1
        sage: kronecker_symbol(3,11)
        1
        sage: partition_sign([5,1,1,1,1,1])
        1

    In both cases, Zolotarev holds.

    REFERENCES:
        http://en.wikipedia.org/wiki/Zolotarev's_lemma
    """
    ans=gap.eval("SignPartition(%s)"%(pi))
    return sage_eval(ans)

def partition_associated(pi):
    """
    partition_associated( pi ) returns the ``associated'' (also called
    ``conjugate'' in the literature) partition of the partition pi which is
    obtained by transposing the corresponding Ferrers diagram.

    Wraps GAP's AssociatedPartition.

    EXAMPLES:
        sage: partition_associated([2,2])
        [2, 2]
        sage: partition_associated([6,3,1])
        [3, 2, 2, 1, 1, 1]
        sage: print ferrers_diagram([6,3,1])
        ******
        ***
        *
        sage: print ferrers_diagram([3,2,2,1,1,1])
        ***
        **
        **
        *
        *
        *
    """
    ans=gap.eval("AssociatedPartition(%s)"%(pi))
    return eval(ans)

## related functions

def bernoulli_polynomial(x,n):
    r"""
    The generating function for the Bernoulli polynomials is
    \[
     \frac{t e^{xt}}{e^t-1}= \sum_{n=0}^\infty B_n(x) \frac{t^n}{n!}.
    \]

    One has $B_n(x) = - n\zeta(1 - n,x)$, where $\zeta(s,x)$ is the
    Hurwitz zeta function.  Thus, in a certain sense, the Hurwitz zeta
    generalizes the Bernoulli polynomials to non-integer values of n.

    EXAMPLES:
        sage: y = QQ['y'].0
        sage: bernoulli_polynomial(y,5)
        y^5 - 5/2*y^4 + 5/3*y^3 - 1/6*y

    REFERENCES:
        http://en.wikipedia.org/wiki/Bernoulli_polynomials
    """
    return sage_eval(maxima.eval("bernpoly(x,%s)"%n), {'x':x})


#def hurwitz_zeta(s,x,N):
#    """
#    Returns the value of the $\zeta(s,x)$ to $N$ decimals, where s and x are real.
#
#    The Hurwitz zeta function is one of the many zeta functions. It defined as
#    \[
#    \zeta(s,x) = \sum_{k=0}^\infty (k+x)^{-s}.
#    \]
#    When $x = 1$, this coincides with Riemann's zeta function. The Dirichlet L-functions
#    may be expressed as a linear combination of Hurwitz zeta functions.
#
#    EXAMPLES:
#        sage: hurwitz_zeta(3,1/2,6)
#        '8.41439b0'
#        sage: hurwitz_zeta(1.1,1/2,6)
#        'Warning:Floattobigfloatconversionof12.1040625294406841.21041b1'
#
#    "b0" can be ignored, but "b1" means that the decimal point was shifted
#    1 to the left (so the answer is not 1.21041b1 but 12.10406).
#
#    REFERENCES:
#        http://en.wikipedia.org/wiki/Hurwitz_zeta_function
#
#    """
#    maxima.eval('load ("bffac")')
#    s = maxima.eval("bfhzeta (%s,%s,%s)"%(s,x,N))
#    return s  ## returns an odd string
