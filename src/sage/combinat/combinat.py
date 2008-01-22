r"""
Combinatorial Functions.

AUTHORS:
        -- David Joyner (2006-07), initial implementation.
        -- William Stein (2006-07), editing of docs and code; many optimizations,
                      refinements, and bug fixes in corner cases
        -- DJ (2006-09): bug fix for combinations, added permutations_iterator,
                      combinations_iterator from Python Cookbook, edited docs.
        -- DJ (2007-11): changed permutations, added hadamard_matrix

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
#                     2007 Mike Hansen <mhansen@gmail.com>,
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

import os

from sage.interfaces.all import gap, maxima
from sage.rings.all import QQ, ZZ
from sage.rings.arith import binomial
from sage.misc.sage_eval import sage_eval
from sage.libs.all import pari
from sage.rings.arith import factorial
from random import randint
from sage.misc.misc import prod
from sage.structure.sage_object import SageObject
import __builtin__
from sage.algebras.algebra import Algebra
from sage.algebras.algebra_element import AlgebraElement
import sage.structure.parent_base
import partitions as partitions_ext

######### combinatorial sequences

def hadamard_matrix(n):
    """
    Returns an n x n Hadamard matrix of order $n$, if possible.

    If the construction of this matrix is not implemented in GUAVA or there is
    no such matrix, raises a NotImplementedError.

    EXAMPLES:
        sage: hadamard_matrix(4)
        [ 1  1  1  1]
        [ 1 -1  1 -1]
        [ 1  1 -1 -1]
        [ 1 -1 -1  1]
        sage: hadamard_matrix(6)
        Traceback (most recent call last):
        ...
        NotImplementedError: Hadamard matrix of order 6 does not exist or is not implemented yet.
    """
    try:
        ans = gap("HadamardMat(%s)"%ZZ(n))
        return ans._matrix_(ZZ)
    except:
        raise NotImplementedError, "Hadamard matrix of order %s does not exist or is not implemented yet."%n

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
    return ZZ(eval(ans))

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
        '-1/2+x^2+x^4+2*x^6+5*x^8+14*x^10+42*x^12+132*x^14'
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
        '1-x^2/2+5*x^4/24-61*x^6/720+277*x^8/8064-50521*x^10/3628800'
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

def mod_stirling(q,n,k):
    """
    """
    if k>n or k<0 or n<0:
        raise ValueError, "n (= %s) and k (= %s) must greater than or equal to 0, and n must be greater than or equal to k"

    if k == 0:
        return 1
    elif k == 1:
        return (n**2+(2*q+1)*n)/2
    elif k == n:
        return prod( [ q+i for i in range(1, n+1) ] )
    else:
        return mod_stirling(q,n-1,k)+(q+n)*mod_stirling(q, n-1, k-1)




class CombinatorialObject(SageObject):
    def __init__(self, l):
        self.list = l
        self._hash = None

    def __str__(self):
        return str(self.list)

    def __repr__(self):
        return self.list.__repr__()

    def __eq__(self, other):
        if isinstance(other, CombinatorialObject):
            return self.list.__eq__(other.list)
        else:
            return self.list.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, CombinatorialObject):
            return self.list.__lt__(other.list)
        else:
            return self.list.__lt__(other)

    def __le__(self, other):
        if isinstance(other, CombinatorialObject):
            return self.list.__le__(other.list)
        else:
            return self.list.__le__(other)

    def __gt__(self, other):
        if isinstance(other, CombinatorialObject):
            return self.list.__gt__(other.list)
        else:
            return self.list.__gt__(other)

    def __ge__(self, other):
        if isinstance(other, CombinatorialObject):
            return self.list.__ge__(other.list)
        else:
            return self.list.__ge__(other)

    def __ne__(self, other):
        if isinstance(other, CombinatorialObject):
            return self.list.__ne__(other.list)
        else:
            return self.list.__ne__(other)

    def __add__(self, other):
        return self.list + other

    def __hash__(self):
        if self._hash is None:
            self._hash = str(self.list).__hash__()
        return self._hash

    #def __cmp__(self, other):
    #    return self.list.__cmp__(other)

    def __len__(self):
        return self.list.__len__()

    def __getitem__(self, key):
        return self.list.__getitem__(key)

    def __iter__(self):
        return self.list.__iter__()

    def __contains__(self, item):
        return self.list.__contains__(item)


    def index(self, key):
        return self.list.index(key)


class CombinatorialClass(SageObject):
    def __len__(self):
        """
        Returns the number of elements in the combinatorial class.

        EXAMPLES:
            sage: len(Partitions(5))
            7
        """
        return self.count()

    def __getitem__(self, i):
        """
        Returns the combinatorial object of rank i.

        EXAMPLES:
            sage: p5 = Partitions(5)
            sage: p5[0]
            [5]
            sage: p5[6]
            [1, 1, 1, 1, 1]
            sage: p5[7]
            Traceback (most recent call last):
            ...
            ValueError: the value must be between 0 and 6 inclusive
        """
        return self.unrank(i)

    def __str__(self):
        """
        Returns a string representation of self.

        EXAMPLES:
            sage: str(Partitions(5))
            'Partitions of the integer 5'
        """
        return self.__repr__()

    def __repr__(self):
        """
        EXAMPLES:
            sage: repr(Partitions(5))
            'Partitions of the integer 5'
        """
        if hasattr(self, '_name') and self._name:
            return self._name
        else:
            return "Combinatorial Class -- REDEFINE ME!"

    def __contains__(self, x):
        """
        Tests whether or not the combinatorial class contains the
        object x.  This raises a NotImplementedError as a default
        since _all_ subclasses of CombinatorialClass should
        override this.

        Note that we could replace this with a default implementation
        that just iterates through the elements of the combinatorial
        class and checks for equality.  However, since we use __contains__
        for type checking, this operation should be cheap and should be
        implemented manually for each combinatorial class.
        """
        raise NotImplementedError

    def __iter__(self):
        """
        Allows the combinatorial class to be treated as an iterator.

        EXAMPLES:
            sage: p5 = Partitions(5)
            sage: [i for i in p5]
            [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
        """
        return self.iterator()

    def __cmp__(self, x):
        """
        Compares two different combinatorial classes.  For now, the comparison
        is done just on their repr's.

        EXAMPLES:
            sage: p5 = Partitions(5)
            sage: p6 = Partitions(6)
            sage: repr(p5) == repr(p6)
            False
            sage: p5 == p6
            False
        """
        return cmp(repr(self), repr(x))

    def __count_from_iterator(self):
        """
        Default implmentation of count which just goes through the iterator
        of the combinatorial class to count the number of objects.
        """
        c = 0
        for x in self.iterator():
            c += 1
        return c
    count = __count_from_iterator

    def __call__(self, x):
        """
        Returns x as an element of the combinatorial class's object class.

        EXAMPLES:
            sage: p5 = Partitions(5)
            sage: a = [2,2,1]
            sage: type(a)
            <type 'list'>
            sage: a = p5(a)
            sage: type(a)
            <class 'sage.combinat.partition.Partition_class'>
            sage: p5([2,1])
            Traceback (most recent call last):
            ...
            ValueError: [2, 1] not in Partitions of the integer 5
        """
        if x in self:
            return self.object_class(x)
        else:
            raise ValueError, "%s not in %s"%(x, self)

    def __list_from_iterator(self):
        """
        The default implementation of list which builds the list from
        the iterator.
        """
        return [x for x in self.iterator()]

    #Set list to the default implementation
    list  = __list_from_iterator

    #Set the default object class to be CombinatorialObject
    object_class = CombinatorialObject

    def __iterator_from_next(self):
        """
        An iterator to use when .first() and .next() are provided.
        """
        f = self.first()
        yield f
        while True:
            try:
                f = self.next(f)
            except:
                break

            if f == None:
                break
            else:
                yield f

    def __iterator_from_previous(self):
        """
        An iterator to use when .last() and .previous() are provided.
        """
        l = self.last()
        yield l
        while True:
            try:
                l = self.previous(l)
            except:
                break

            if l == None:
                break
            else:
                yield l

    def __iterator_from_unrank(self):
        """
        An iterator to use when .unrank() is provided.
        """
        r = 0
        u = self.unrank(r)
        yield f
        while True:
            r += 1
            try:
                u = self.unrank(l)
            except:
                break

            if u == None:
                break
            else:
                yield u

    def __iterator_from_list(self):
        """
        An iterator to use when .list() is provided()
        """
        for x in self.list():
            yield x

    def iterator(self):
        """
        Default implementation of iterator.
        """
        #Check to see if .first() and .next() are overridden in the subclass
        if ( self.first != self.__first_from_iterator and
             self.next  != self.__next_from_iterator ):
            return self.__iterator_from_next()
        #Check to see if .last() and .previous() are overridden in the subclass
        elif ( self.last != self.__last_from_iterator and
               self.previous != self.__previous_from_iterator):
            return self.__iterator_from_previous()
        #Check to see if .unrank() is overridden in the subclass
        elif self.unrank != self.__unrank_from_iterator:
            return self.__iterator_from_unrank()
        #Finally, check to see if .list() is overridden in the subclass
        elif self.list != self.__list_from_iterator:
            return self.__iterator_from_list()
        else:
            raise NotImplementedError, "iterator called but not implemented"


    def __list_from_unrank_and_count(self):
        return [ self.unrank(i) for i in range(self.count()) ]

    def __unrank_from_iterator(self, r):
        """
        Default implementation of unrank which goes through the iterator.
        """
        counter = 0
        for u in self.iterator():
            if counter == r:
                return u
            counter += 1
        raise ValueError, "the value must be between %s and %s inclusive"%(0,counter-1)

    #Set the default implementation of unrank
    unrank = __unrank_from_iterator


    def __random_from_unrank(self):
        """
        Default implementation of random which uses unrank.
        """
        c = self.count()
        r = randint(0, c-1)
        if hasattr(self, 'object_class'):
            return self.object_class(self.unrank(r))
        else:
            return self.unrank(r)

    #Set the default implementation of random
    random = __random_from_unrank


    def __rank_from_iterator(self, obj):
        r = 0
        for i in self.iterator():
            if i == obj:
                return r
            r += 1

    rank =  __rank_from_iterator

    def __first_from_iterator(self):
        for i in self.iterator():
            return i

    first = __first_from_iterator

    def __last_from_iterator(self):
        for i in self.iterator():
            pass
        return i

    last = __last_from_iterator

    def __next_from_iterator(self, obj):
        if hasattr(obj, 'next'):
            res = obj.next()
            if res:
                return res
            else:
                return None
        found = False
        for i in self.iterator():
            if found:
                return i
            if i == obj:
                found = True
        return None

    next = __next_from_iterator

    def __previous_from_iterator(self, obj):
        if hasattr(obj, 'previous'):
            res = obj.previous()
            if res:
                return res
            else:
                return None
        prev = None
        for i in self.iterator():
            if i == obj:
                break
            prev = i
        return prev

    previous = __previous_from_iterator

    def filter(self, f, name=None):
        """
        Returns the combinatorial subclass of f which consists of
        the elements x of self such that f(x) is True.

        EXAMPLES:
            sage: P = Permutations(3).filter(lambda x: x.avoids([1,2]))
            sage: P.list()
            [[3, 2, 1]]
        """
        return FilteredCombinatorialClass(self, f, name=name)

    def union(self, right_cc, name=None):
        """
        Returns the combinatorial class representing the union
        of self and right_cc.

        EXAMPLES:
            sage: P = Permutations(2).union(Permutations(1))
            sage: P.list()
            [[1, 2], [2, 1], [1]]
        """
        return UnionCombinatorialClass(self, right_cc, name=name)

class FilteredCombinatorialClass(CombinatorialClass):
    def __init__(self, combinatorial_class, f, name=None):
        """
        TESTS:
            sage: P = Permutations(3).filter(lambda x: x.avoids([1,2]))
        """
        self.f = f
        self.combinatorial_class = combinatorial_class
        self._name = name

    def __repr__(self):
        if self._name:
            return self._name
        else:
            return "Filtered sublass of " + repr(self.combinatorial_class)

    def __contains__(self, x):
        return x in self.combinatorial_class and self.f(x)

    def count(self):
        c = 0
        for x in self.iterator():
            c += 1
        return c

    def list(self):
        res = []
        for x in self.combinatorial_class.iterator():
            if self.f(x):
                res.append(x)
        return res

    def iterator(self):
        for x in self.combinatorial_class.iterator():
            if self.f(x):
                yield x

    def __unrank_from_iterator(self, r):
        """
        Default implementation of unrank which goes through the iterator.
        """
        counter = 0
        for u in self.iterator():
            if counter == r:
                return u
            counter += 1
        raise ValueError, "the value must be between %s and %s inclusive"%(0,counter-1)

    #Set the default implementation of unrank
    unrank = __unrank_from_iterator


    def __random_from_unrank(self):
        """
        Default implementation of random which uses unrank.
        """
        c = self.count()
        r = randint(0, c-1)
        if hasattr(self, 'object_class'):
            return self.object_class(self.unrank(r))
        else:
            return self.unrank(r)

    #Set the default implementation of random
    random = __random_from_unrank


    def __rank_from_iterator(self, obj):
        r = 0
        for i in self.iterator():
            if i == obj:
                return r
            r += 1

    rank =  __rank_from_iterator

    def __first_from_iterator(self):
        for i in self.iterator():
            return i

    first = __first_from_iterator

    def __last_from_iterator(self):
        for i in self.iterator():
            pass
        return i

    last = __last_from_iterator

    def __next_from_iterator(self, obj):
        if hasattr(obj, 'next'):
            res = obj.next()
            if res:
                return res
            else:
                return None
        found = False
        for i in self.iterator():
            if found:
                return i
            if i == obj:
                found = True
        return None

    next = __next_from_iterator

    def __previous_from_iterator(self, obj):
        if hasattr(obj, 'previous'):
            res = obj.previous()
            if res:
                return res
            else:
                return None
        prev = None
        for i in self.iterator():
            if i == obj:
                break
            prev = i
        return prev

    previous = __previous_from_iterator


class UnionCombinatorialClass(CombinatorialClass):
    def __init__(self, left_cc, right_cc, name=None):
        """
        TESTS:
            sage: P = Permutations(3).union(Permutations(2))
            sage: P == loads(dumps(P))
            True
        """
        self.left_cc = left_cc
        self.right_cc = right_cc
        self._name = name

    def __repr__(self):
        """
        TESTS:
            sage: print repr(Permutations(3).union(Permutations(2)))
            Union combinatorial class of
                Standard permutations of 3
            and
                Standard permutations of 2

        """
        if self._name:
            return self._name
        else:
            return "Union combinatorial class of \n    %s\nand\n    %s"%(self.left_cc, self.right_cc)

    def __contains__(self, x):
        return x in self.left_cc or x in self.right_cc

    def count(self):
        return self.left_cc.count() + self.right_cc.count()

    def list(self):
        res = []
        for x in self.left_cc.iterator():
            res.append(x)
        for x in self.right_cc.iterator():
            res.append(x)
        return res

    def iterator(self):
        for x in self.left_cc.iterator():
            yield x
        for x in self.right_cc.iterator():
            yield x

    def __unrank_from_iterator(self, r):
        """
        Default implementation of unrank which goes through the iterator.
        """
        counter = 0
        for u in self.iterator():
            if counter == r:
                return u
            counter += 1
        raise ValueError, "the value must be between %s and %s inclusive"%(0,counter-1)

    #Set the default implementation of unrank
    unrank = __unrank_from_iterator

    def __random_from_unrank(self):
        """
        Default implementation of random which uses unrank.
        """
        c = self.count()
        r = randint(0, c-1)
        if hasattr(self, 'object_class'):
            return self.object_class(self.unrank(r))
        else:
            return self.unrank(r)

    #Set the default implementation of random
    random = __random_from_unrank

    def __rank_from_iterator(self, obj):
        r = 0
        for i in self.iterator():
            if i == obj:
                return r
            r += 1

    rank =  __rank_from_iterator

    def __first_from_iterator(self):
        for i in self.iterator():
            return i

    first = __first_from_iterator

    def __last_from_iterator(self):
        for i in self.iterator():
            pass
        return i

    last = __last_from_iterator

    def __next_from_iterator(self, obj):
        if hasattr(obj, 'next'):
            res = obj.next()
            if res:
                return res
            else:
                return None
        found = False
        for i in self.iterator():
            if found:
                return i
            if i == obj:
                found = True
        return None

    next = __next_from_iterator

    def __previous_from_iterator(self, obj):
        if hasattr(obj, 'previous'):
            res = obj.previous()
            if res:
                return res
            else:
                return None
        prev = None
        for i in self.iterator():
            if i == obj:
                break
            prev = i
        return prev

    previous = __previous_from_iterator


def hurwitz_zeta(s,x,N):
    """
    Returns the value of the $\zeta(s,x)$ to $N$ decimals, where s and x are real.

    The Hurwitz zeta function is one of the many zeta functions. It defined as
    \[
    \zeta(s,x) = \sum_{k=0}^\infty (k+x)^{-s}.
    \]
    When $x = 1$, this coincides with Riemann's zeta function. The Dirichlet L-functions
    may be expressed as a linear combination of Hurwitz zeta functions.

    EXAMPLES:
        sage: hurwitz_zeta(3,1/2,6)
        8.41439000000000
        sage: hurwitz_zeta(1.1,1/2,6)
        12.1041000000000
        sage: hurwitz_zeta(1.1,1/2,50)
        12.103813495683744469025853545548130581952676591199

    REFERENCES:
        http://en.wikipedia.org/wiki/Hurwitz_zeta_function

    """
    maxima.eval('load ("bffac")')
    s = maxima.eval("bfhzeta (%s,%s,%s)"%(s,x,N))

    #Handle the case where there is a 'b' in the string
    #'1.2000b0' means 1.2000 and
    #'1.2000b1' means 12.000
    i = s.rfind('b')
    if i == -1:
        return sage_eval(s)
    else:
        if s[i+1:] == '0':
            return sage_eval(s[:i])
        else:
            return sage_eval(s[:i])*10**sage_eval(s[i+1:])

    return s  ## returns an odd string


#####################################################
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
    Calls a function written by Mike Hansen, not GAP.

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
        sage: MS = MatrixSpace(GF(2),2,2)
        sage: A = MS([1,0,1,1])
        sage: permutations(A.rows())
        [[(1, 0), (1, 1)], [(1, 1), (1, 0)]]

    """
    from sage.combinat.permutation import Permutations
    ans = Permutations(mset)
    return ans.list()

def permutations_iterator(mset,n=None):
    """
    Posted by Raymond Hettinger, 2006/03/23, to the Python Cookbook:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/474124

    Note-- This function considers repeated elements as different entries,
    so for example:
        sage: from sage.combinat.combinat import permutations, permutations_iterator
        sage: mset = [1,2,2]
        sage: permutations(mset)
        [[1, 2, 2], [2, 1, 2], [2, 2, 1]]
        sage: for p in permutations_iterator(mset): print p
        [1, 2, 2]
        [1, 2, 2]
        [2, 1, 2]
        [2, 2, 1]
        [2, 1, 2]
        [2, 2, 1]


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

    AUTHOR: Robert L. Miller

    EXAMPLES:
        sage: mset = [1,1,2,2,2]
        sage: number_of_permutations(mset)
        10

    """
    from sage.rings.arith import factorial, prod
    m = len(mset)
    n = []
    seen = []
    for element in mset:
        try:
            n[seen.index(element)] += 1
        except:
            n.append(1)
            seen.append(element)
    return factorial(m)/prod([factorial(k) for k in n])

def cyclic_permutations(mset):
    """
    Returns a list of all cyclic permutations of mset. Treats mset as a list,
    not a set, i.e. entries with the same value are distinct.

    AUTHOR: Emily Kirkman

    EXAMPLES:
        sage: from sage.combinat.combinat import cyclic_permutations, cyclic_permutations_iterator
        sage: cyclic_permutations(range(4))
        [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3], [0, 2, 3, 1], [0, 3, 1, 2], [0, 3, 2, 1]]
        sage: for cycle in cyclic_permutations(['a', 'b', 'c']):
        ...       print cycle
        ['a', 'b', 'c']
        ['a', 'c', 'b']

    Note that lists with repeats are not handled intuitively:
        sage: cyclic_permutations([1,1,1])
        [[1, 1, 1], [1, 1, 1]]

    """
    return list(cyclic_permutations_iterator(mset))

def cyclic_permutations_iterator(mset):
    """
    Iterates over all cyclic permutations of mset in cycle notation. Treats
    mset as a list, not a set, i.e. entries with the same value are distinct.

    AUTHOR: Emily Kirkman

    EXAMPLES:
        sage: from sage.combinat.combinat import cyclic_permutations, cyclic_permutations_iterator
        sage: cyclic_permutations(range(4))
        [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3], [0, 2, 3, 1], [0, 3, 1, 2], [0, 3, 2, 1]]
        sage: for cycle in cyclic_permutations(['a', 'b', 'c']):
        ...       print cycle
        ['a', 'b', 'c']
        ['a', 'c', 'b']

    Note that lists with repeats are not handled intuitively:
        sage: cyclic_permutations([1,1,1])
        [[1, 1, 1], [1, 1, 1]]

    """
    if len(mset) > 2:
        for perm in permutations_iterator(mset[1:]):
            yield [mset[0]] + perm
    else:
        yield mset

#### partitions

def partitions_set(S,k=None, use_file=True):
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

    WARNING: This function is inefficient.  The runtime is dominated
    by parsing the output from GAP.

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
    if k is None:
        ans=gap("PartitionsSet(%s)"%S).str(use_file=use_file)
    else:
        ans=gap("PartitionsSet(%s,%s)"%(S,k)).str(use_file=use_file)
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
    if k is None:
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

    INPUT:
        n -- a positive integer

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
    n = ZZ(n)
    if n <= 0:
        raise ValueError, "n (=%s) must be a positive integer"%n
    if k is None:
        ans=gap.eval("Partitions(%s)"%(n))
    else:
        ans=gap.eval("Partitions(%s,%s)"%(n,k))
    return eval(ans.replace('\n',''))

def number_of_partitions(n,k=None, algorithm='default'):
    r"""
    Returns the size of partitions_list(n,k).

    INPUT:
        n -- an integer
        k -- (default: None); if specified, instead returns the
             cardinality of the set of all (unordered) partitions of
             the positive integer n into sums with k summands.
        algorithm -- (default: 'default')
            'default' -- If k is not None, then use Gap (very slow).
                         If k is None, use Jon Bober's highly
                         optimized implementation (this is the fastest
                         code in the world for this problem).
            'bober' -- use Jonathon Bober's implementation
            'gap' -- use GAP (VERY *slow*)
            'pari' -- use PARI.  Speed seems the same as GAP until $n$ is
                      in the thousands, in which case PARI is faster. *But*
                      PARI has a bug, e.g., on 64-bit Linux PARI-2.3.2
                      outputs numbpart(147007)%1000 as 536, but it
                      should be 533!.  So do not use this option.

    IMPLEMENTATION: Wraps GAP's NrPartitions or PARI's numbpart function.

    Use the function \code{partitions(n)} to return a generator over
    all partitions of $n$.

    It is possible to associate with every partition of the integer n
    a conjugacy class of permutations in the symmetric group on n
    points and vice versa.  Therefore p(n) = NrPartitions(n) is the
    number of conjugacy classes of the symmetric group on n points.

    EXAMPLES:
        sage: v = list(partitions(5)); v
        [(1, 1, 1, 1, 1), (1, 1, 1, 2), (1, 2, 2), (1, 1, 3), (2, 3), (1, 4), (5,)]
        sage: len(v)
        7
        sage: number_of_partitions(5, algorithm='gap')
        7
        sage: number_of_partitions(5, algorithm='pari')
        7
        sage: number_of_partitions(5, algorithm='bober')
        7

    The input must be a nonnegative integer or a ValueError is raised.
        sage: number_of_partitions(-5)
        Traceback (most recent call last):
        ...
        ValueError: n (=-5) must be a nonnegative integer

        sage: number_of_partitions(10,2)
        5
        sage: number_of_partitions(10)
        42
        sage: number_of_partitions(3)
        3
        sage: number_of_partitions(10)
        42
        sage: number_of_partitions(3, algorithm='pari')
        3
        sage: number_of_partitions(10, algorithm='pari')
        42
        sage: number_of_partitions(40)
        37338
        sage: number_of_partitions(100)
        190569292
        sage: number_of_partitions(100000)
        27493510569775696512677516320986352688173429315980054758203125984302147328114964173055050741660736621590157844774296248940493063070200461792764493033510116079342457190155718943509725312466108452006369558934464248716828789832182345009262853831404597021307130674510624419227311238999702284408609370935531629697851569569892196108480158600569421098519

    A generating function for p(n) is given by the reciprocal of
    Euler's function:

    \[
    \sum_{n=0}^\infty p(n)x^n = \prod_{k=1}^\infty \left(\frac {1}{1-x^k} \right).
    \]

    We use SAGE to verify that the first several coefficients do
    instead agree:

        sage: q = PowerSeriesRing(QQ, 'q', default_prec=9).gen()
        sage: prod([(1-q^k)^(-1) for k in range(1,9)])  ## partial product of
        1 + q + 2*q^2 + 3*q^3 + 5*q^4 + 7*q^5 + 11*q^6 + 15*q^7 + 22*q^8 + O(q^9)
        sage: [number_of_partitions(k) for k in range(2,10)]
        [2, 3, 5, 7, 11, 15, 22, 30]

    REFERENCES:
        http://en.wikipedia.org/wiki/Partition_%28number_theory%29

    TESTS:
        sage: n = 500 + randint(0,500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1500 + randint(0,1500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 100000000 + randint(0,100000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000000 + randint(0,1000000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0      # takes a long time
        True

    Another consistency test for n up to 500:
        sage: len([n for n in [1..500] if number_of_partitions(n) != number_of_partitions(n,algorithm='pari')])
        0
    """
    n = ZZ(n)
    if n < 0:
        raise ValueError, "n (=%s) must be a nonnegative integer"%n
    elif n == 0:
        return ZZ(1)

    if algorithm == 'default':
        if k is None:
            algorithm = 'bober'
        else:
            algorithm = 'gap'

    if algorithm == 'gap':
        if k is None:
            ans=gap.eval("NrPartitions(%s)"%(ZZ(n)))
        else:
            ans=gap.eval("NrPartitions(%s,%s)"%(ZZ(n),ZZ(k)))
        return ZZ(ans)

    if k is not None:
        raise ValueError, "only the GAP algorithm works if k is specified."

    if algorithm == 'bober':
        return partitions_ext.number_of_partitions(n)

    elif algorithm == 'pari':
        return ZZ(pari(ZZ(n)).numbpart())

    raise ValueError, "unknown algorithm '%s'"%algorithm

def partitions(n):
    r"""
    Generator of all the partitions of the integer $n$.

    INPUT:
        n -- int

    To compute the number of partitions of $n$ use
    \code{number_of_partitions(n)}.

    EXAMPLES:
        sage: partitions(3)          # random location
        <generator object at 0xab3b3eac>
        sage: list(partitions(3))
        [(1, 1, 1), (1, 2), (3,)]


    AUTHOR: Adapted from David Eppstein, Jan Van lent, George Yoshida;
    Python Cookbook 2, Recipe 19.16.
    """
    n == ZZ(n)
    # base case of the recursion: zero is the sum of the empty tuple
    if n == 0:
        yield ( )
        return
    # modify the partitions of n-1 to form the partitions of n
    for p in partitions(n-1):
        yield (1,) + p
        if p and (len(p) < 2 or p[1] > p[0]):
            yield (p[0] + 1,) + p[1:]

def cyclic_permutations_of_partition(partition):
    """
    Returns all combinations of cyclic permutations of each cell of the
    partition.

    AUTHOR: Robert L. Miller

    EXAMPLES:
        sage: from sage.combinat.combinat import cyclic_permutations_of_partition
        sage: cyclic_permutations_of_partition([[1,2,3,4],[5,6,7]])
        [[[1, 2, 3, 4], [5, 6, 7]],
         [[1, 2, 4, 3], [5, 6, 7]],
         [[1, 3, 2, 4], [5, 6, 7]],
         [[1, 3, 4, 2], [5, 6, 7]],
         [[1, 4, 2, 3], [5, 6, 7]],
         [[1, 4, 3, 2], [5, 6, 7]],
         [[1, 2, 3, 4], [5, 7, 6]],
         [[1, 2, 4, 3], [5, 7, 6]],
         [[1, 3, 2, 4], [5, 7, 6]],
         [[1, 3, 4, 2], [5, 7, 6]],
         [[1, 4, 2, 3], [5, 7, 6]],
         [[1, 4, 3, 2], [5, 7, 6]]]

    Note that repeated elements are not considered equal:
        sage: cyclic_permutations_of_partition([[1,2,3],[4,4,4]])
        [[[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]],
         [[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]]]

    """
    return list(cyclic_permutations_of_partition_iterator(partition))

def cyclic_permutations_of_partition_iterator(partition):
    """
    Iterates over all combinations of cyclic permutations of each cell of the
    partition.

    AUTHOR: Robert L. Miller

    EXAMPLES:
        sage: from sage.combinat.combinat import cyclic_permutations_of_partition
        sage: cyclic_permutations_of_partition([[1,2,3,4],[5,6,7]])
        [[[1, 2, 3, 4], [5, 6, 7]],
         [[1, 2, 4, 3], [5, 6, 7]],
         [[1, 3, 2, 4], [5, 6, 7]],
         [[1, 3, 4, 2], [5, 6, 7]],
         [[1, 4, 2, 3], [5, 6, 7]],
         [[1, 4, 3, 2], [5, 6, 7]],
         [[1, 2, 3, 4], [5, 7, 6]],
         [[1, 2, 4, 3], [5, 7, 6]],
         [[1, 3, 2, 4], [5, 7, 6]],
         [[1, 3, 4, 2], [5, 7, 6]],
         [[1, 4, 2, 3], [5, 7, 6]],
         [[1, 4, 3, 2], [5, 7, 6]]]

    Note that repeated elements are not considered equal:
        sage: cyclic_permutations_of_partition([[1,2,3],[4,4,4]])
        [[[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]],
         [[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]]]

    """
    if len(partition) == 1:
        for i in cyclic_permutations_iterator(partition[0]):
            yield [i]
    else:
        for right in cyclic_permutations_of_partition_iterator(partition[1:]):
            for perm in cyclic_permutations_iterator(partition[0]):
                yield [perm] + right

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
    if k is None:
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
    if k is None:
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
    if k is None:
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
    if k is None:
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


def fibonacci_sequence(start, stop=None, algorithm=None):
    r"""
    Returns an iterator over the Fibonacci sequence, for all fibonacci numbers
    $f_n$ from \code{n = start} up to (but not including) \code{n = stop}

    INPUT:
        start -- starting value
        stop -- stopping value
        algorithm -- default (None) -- passed on to fibonacci function (or
                     not passed on if None, i.e., use the default).


    EXAMPLES:
        sage: fibs = [i for i in fibonacci_sequence(10, 20)]
        sage: fibs
        [55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181]

        sage: sum([i for i in fibonacci_sequence(100, 110)])
        69919376923075308730013

    SEE ALSO: fibonacci_xrange

    AUTHOR:
        Bobby Moretti
    """
    from sage.rings.integer_ring import ZZ
    if stop is None:
        stop = ZZ(start)
        start = ZZ(0)
    else:
        start = ZZ(start)
        stop = ZZ(stop)

    if algorithm:
        for n in xrange(start, stop):
            yield fibonacci(n, algorithm=algorithm)
    else:
        for n in xrange(start, stop):
            yield fibonacci(n)

def fibonacci_xrange(start, stop=None, algorithm='pari'):
    r"""
    Returns an iterator over all of the Fibonacci numbers in the given range,
    including \code{f_n = start} up to, but not including, \code{f_n = stop}.

    EXAMPLES:
        sage: fibs_in_some_range =  [i for i in fibonacci_xrange(10^7, 10^8)]
        sage: len(fibs_in_some_range)
        4
        sage: fibs_in_some_range
        [14930352, 24157817, 39088169, 63245986]

        sage: fibs = [i for i in fibonacci_xrange(10, 100)]
        sage: fibs
        [13, 21, 34, 55, 89]

        sage: list(fibonacci_xrange(13, 34))
        [13, 21]

    A solution to the second Project Euler problem:
        sage: sum([i for i in fibonacci_xrange(10^6) if is_even(i)])
        1089154

    SEE ALSO: fibonacci_sequence

    AUTHOR:
        Bobby Moretti
    """
    from sage.rings.integer_ring import ZZ
    if stop is None:
        stop = ZZ(start)
        start = ZZ(0)
    else:
        start = ZZ(start)
        stop = ZZ(stop)

    # iterate until we've gotten high enough
    fn = 0
    n = 0
    while fn < start:
        n += 1
        fn = fibonacci(n)

    while True:
        fn = fibonacci(n)
        n += 1
        if fn < stop:
            yield fn
        else:
            return


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
