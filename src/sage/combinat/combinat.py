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

from sage.interfaces.all import gap, maxima
from sage.rings.all import QQ, RR, ZZ
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
        return str(self.list).__hash__()

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
        return self.count()

    def __getitem__(self, i):
        return self.unrank(i)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "Combinatorial Class"

    def __contains__(self, x):
        return False

    def __iter__(self):
        return self.iterator()

    def __cmp__(self, x):
        return cmp(repr(self), repr(x))

    def __count_from_iterator(self):
        c = 0
        for x in self.iterator():
            c += 1
        return c
    count = __count_from_iterator

    def __call__(self, x):
        if x in self:
            return self.object_class(x)
        else:
            raise ValueError, "%s not in %s"%(x, self)

    def __list_from_iterator(self):
        res = []
        for x in self.iterator():
            res.append(x)
        return res

    list  = __list_from_iterator

    object_class = CombinatorialObject

    def __iterator_from_next(self):
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
        for x in self.list():
            yield x

    def iterator(self):
        if ( self.first != self.__first_from_iterator and
             self.next  != self.__next_from_iterator ):
            return self.__iterator_from_next()
        elif ( self.last != self.__last_from_iterator and
               self.previous != self.__previous_from_iterator):
            return self.__iterator_from_previous()
        elif self.unrank != self.__unrank_from_iterator:
            return self.__iterator_from_unrank()
        elif self.list != self.__list_from_iterator:
            return self.__iterator_from_list()
        else:
            raise NotImplementedError, "iterator called but not implemented"


    def __list_from_unrank_and_count(self):
        return [ self.unrank(i) for i in range(self.count()) ]

    def __unrank_from_iterator(self, r):
        counter = 0
        for u in self.iterator():
            if counter == r:
                return u
            counter += 1
        raise ValueError, "the value must be between %s and %s inclusive"%(0,counter-1)

    def __unrank_from_iterator(self, r):
        counter = 0
        for u in self.iterator():
            if counter == r:
                return u
            counter += 1

    unrank = __unrank_from_iterator

    def __random_from_unrank(self):
        c = self.count()
        r = randint(0, c-1)
        if hasattr(self, 'object_class'):
            return self.object_class(self.unrank(r))
        else:
            return self.unrank(r)

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

class ExampleCombinatorialClass_size(CombinatorialClass):
    def __init__(self, size):
        self.size = size

    def __contains__(self, x):
        if not isinstance(x, __builtin__.list):
            return False
        if len(x) != self.size:
            return False
        for i in x:
            if not ( i == 0 or i == 1 ):
                return False
        return True

    def size(self):
        return self.size

    def first(self):
        return [0]*int(self.size)

    def next(self, x):
        try:
            pos = x.index(0)
        except:
            return None

        return [0]*int(pos) + [1] + x[pos+1:]

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


