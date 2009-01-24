r"""
Bounds for Parameters of Codes

AUTHORS:
        -- David Joyner (2006-07), initial implementation.
        -- William Stein (2006-07), minor editing of docs and code
                                    (fixed bug in elias_bound_asymp)
        -- D. Joyner (2006-07), fixed dimension_upper_bound to return an integer,
                      added example to elias_bound_asymp.

This module provided some upper and lower bounds for the parameters of codes.

Let $ F$ be a finite field (we denote the finite field with $q$ elements
$GF(q)$ by $\FF_q$). A subset $ C$ of $ V=F^n$ is called
a {\it code} of length $ n$. A subspace of $ V$ (with the standard basis)
is called a {\it linear code} of length $ n$. If its
dimension is denoted $k$ then we typically store a basis
of $C$ as a $k\times n$ matrix (the rows are the basis vectors).
If $ F=\FF_2$ then $ C$ is called a {\it binary code}.
If $ F$ has $ q$ elements then $ C$ is called a
{\it $ q$-ary code}. The elements of a code $ C$ are called
{\it codewords}. The {\it information rate} of $ C$ is

\[
R={\frac{\log_q\vert C\vert}{n}},
\]
where $ \vert C\vert$ denotes the number of elements of $ C$.
If $ {\bf v}=(v_1,v_2,...,v_n)$, $ {\bf w}=(w_1,w_2,...,w_n)$
are vectors in $ V=F^n$ then we define

\[
d({\bf v},{\bf w}) =\vert\{i\ \vert\ 1\leq i\leq n,\ v_i\not= w_i\}\vert
\]
to be the {\it Hamming distance} between $ {\bf v}$ and $ {\bf w}$.
The function $ d:V\times V\rightarrow \mathbf{N}$ is called the
{\it Hamming metric}. The {\it weight} of a vector (in the
Hamming metric) is $ d({\bf v},{\bf 0})$. The {\it minimum distance}
of a linear code is the smallest non-zero weight of a codeword
in $C$. The {\it relatively minimum distance} is denoted

\[
\delta = d/n.
\]
A linear code with length $n$, dimension $k$,
and minimum distance $d$ is called an {\it $[n,k,d]_q$-code} and
$n,k,d$ are called its {\it parameters}.
A (not necessarily linear) code $C$ with length $n$, size $M=|C|$,
and minimum distance $d$ is called an {\it $(n,M,d)_q$-code} (using
parantheses instead of square brackets).
Of course, $k=\log_q(M)$ for linear codes.

What is the ``best'' code of a given length?
Let $ F$ be a finite field with $q$ elements. Let $A_q(n,d)$ denote
the largest $M$ such that there exists a $(n,M,d)$ code in $ F^n$.
Let $B_q(n,d)$ (also denoted $A^{lin}_q(n,d)$) denote
the largest $k$ such that there exists a $[n,k,d]$ code in $ F^n$.
(Of course, $A_q(n,d)\geq B_q(n,d)$.)
Determining $A_q(n,d)$ and $B_q(n,d)$ is one of the main problems in the theory of
error-correcting codes.

These quantities related to solving a generalization of the childhood
game of ``20 questions''.

GAME: Player 1 secretly chooses a number from $1$ to $M$ ($M$ is large but fixed).
Player 2 asks a series of ``yes/no questions'' in an attempt to determine that number.
Player 1 may lie at most $e$ times ($e\geq 0$ is fixed). What is the minimum
number of ``yes/no questions'' Player 2 must ask to (always) be able to correctly
determine the number Player 1 chose?

If feedback is not allowed (the only situation considered here), call this
minimum number $g(M,e)$.

Lemma: For fixed $e$ and $M$, $g(M,e)$ is the smallest $n$ such that
$A_2(n,2e+1)\geq M$.

Thus, solving the solving a generalization of the game of ``20 questions''
is equivalent to determining $A_2(n,d)$! Using SAGE, you can determine the best
known estimates for this number in 2 ways:

(1) Indirectly, using minimum_distance_lower_bound(n,k,F) and
    minimum_distance_upper_bound(n,k,F) (both of which which connect to the
    internet using Steven Sivek's linear_code_bound(q,n,k))
(2) codesize_upper_bound(n,d,q), dimension_upper_bound(n,d,q), which use
    GUAVA's UpperBound( n, d, q )

This module implements:
\begin{itemize}
\item codesize_upper_bound(n,d,q), for the best known (as of May, 2006)
  upper bound A(n,d) for the size of a code of length n, minimum distance d
  over a field of size q.
\item dimension_upper_bound(n,d,q), an upper bound $B(n,d)=B_q(n,d)$ for the
  dimension of a linear code of length n, minimum distance d over a field of size q.
\item gilbert_lower_bound(n,q,d), a lower bound for number of elements in
  the largest code of min distance d in $\FF_q^n$.
\item gv_info_rate(n,delta,q), $log_q(GLB)/n$, where GLB is the Gilbert lower bound
  and delta = d/n.
\item gv_bound_asymp(delta,q), asymptotic analog of Gilbert lower bound.
\item plotkin_upper_bound(n,q,d)
\item plotkin_bound_asymp(delta,q), asymptotic analog of Plotkin bound.
\item griesmer_upper_bound(n,q,d)
\item elias_upper_bound(n,q,d)
\item elias_bound_asymp(delta,q), asymptotic analog of Elias bound.
\item hamming_upper_bound(n,q,d)
\item hamming_bound_asymp(delta,q), asymptotic analog of Hamming bound.
\item singleton_upper_bound(n,q,d)
\item singleton_bound_asymp(delta,q), asymptotic analog of Singleton bound.
\item mrrw1_bound_asymp(delta,q), ``first'' asymptotic McEliese-Rumsey-Rodemich-Welsh
  bound for the information rate.
\end{itemize}

PROBLEM:
    In this module we shall typically either
    (a) seek bounds on k, given n, d, q,
    (b) seek bounds on R, delta, q (assuming n is ``infinity'').


TODO: * Johnson bounds for binary codes.

      * mrrw2_bound_asymp(delta,q), ``second'' asymptotic
        McEliese-Rumsey-Rodemich-Welsh bound for the information rate.

REFERENCES:
    C. Huffman, V. Pless, {\bf Fundamentals of error-correcting codes},
    Cambridge Univ. Press, 2003.

"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner <wdj@usna.edu>
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.interfaces.all import gap
from sage.rings.all import QQ, RR, ZZ, RDF
from sage.rings.arith import factorial
from sage.misc.functional import log, sqrt

def codesize_upper_bound(n,d,q):
    r"""
    Returns the best known upper bound $A(n,d)=A_q(n,d)$ for the size of a code of length n,
    minimum distance d over a field of size q. The function first checks for trivial cases
    (like d=1 or n=d), and if the value is in the built-in table. Then it calculates the
    minimum value of the upper bound using the methods of Singleton, Hamming, Johnson, Plotkin
    and Elias. If the code is binary, $A(n, 2\ell-1) = A(n+1,2\ell)$, so the function takes
    the minimum of the values obtained from all methods for the parameters
    $(n, 2\ell-1)$ and $(n+1, 2\ell)$.

    Wraps GUAVA's UpperBound( n, d, q ).

    EXAMPLES:
        sage: codesize_upper_bound(10,3,2)
        85
    """
    return int(gap.eval("UpperBound(%s,%s,%s)"%( n, d, q )))

def dimension_upper_bound(n,d,q):
    r"""
    Returns an upper bound $B(n,d) = B_q(n,d)$ for the dimension of a linear code of length n,
    minimum distance d over a field of size q.

    EXAMPLES:
        sage: dimension_upper_bound(10,3,2)
        6

    """
    q = ZZ(q)
    return int(log(codesize_upper_bound(n,d,q),q))

def volume_hamming(n,q,r):
    r"""
    Returns number of elements in a Hamming ball of radius r in $\FF_q^n$.

    EXAMPLES:
        sage: volume_hamming(10,2,3)
        176

    """
    ans=sum([factorial(n)/(factorial(i)*factorial(n-i))*(q-1)**i for i in range(r+1)])
    return ans

def gilbert_lower_bound(n,q,d):
    r"""
    Returns lower bound for number of elements in
    the largest code of minimum distance d in $\FF_q^n$.

    EXAMPLES:
        sage: gilbert_lower_bound(10,2,3)
        128/7

    """
    ans=q**n/volume_hamming(n,q,d-1)
    return ans

def plotkin_upper_bound(n,q,d):
    r"""
    Returns Plotkin upper bound for number of elements in
    the largest code of minimum distance d in $\FF_q^n$. Wraps GAP's
    UpperBoundPlotkin.

    EXAMPLES:
        sage: plotkin_upper_bound(10,2,3)
        192

    """

    ans=gap.eval("UpperBoundPlotkin(%s,%s,%s)"%(n,d,q))
    return QQ(ans)

def griesmer_upper_bound(n,q,d):
    r"""
    Returns the Griesmer upper bound for number of elements in
    the largest code of minimum distance d in $\FF_q^n$. Wraps GAP's
    UpperBoundGriesmer.

    EXAMPLES:
        sage: griesmer_upper_bound(10,2,3)
        128

    """

    ans=gap.eval("UpperBoundGriesmer(%s,%s,%s)"%(n,d,q))
    return QQ(ans)

def elias_upper_bound(n,q,d):
    r"""
    Returns the Elias upper bound for number of elements in
    the largest code of minimum distance d in $\FF_q^n$. Wraps GAP's
    UpperBoundElias.

    EXAMPLES:
        sage: elias_upper_bound(10,2,3)
        232

    """

    ans=gap.eval("UpperBoundElias(%s,%s,%s)"%(n,d,q))
    return QQ(ans)

def hamming_upper_bound(n,q,d):
    r"""
    Returns the Hamming upper bound for number of elements in
    the largest code of minimum distance d in $\FF_q^n$. Wraps GAP's
    UpperBoundHamming.

    The Hamming bound (also known as the sphere packing bound) returns an
    upper bound on the size of a code of length n, minimum distance d, over
    a field of size q. The Hamming bound is obtained by dividing the contents
    of the entire space $\FF_q^n$ by the contents of a ball with radius
    floor((d-1)/2). As all these balls are disjoint, they can never contain more
    than the whole vector space.

    \[ M \leq {q^n \over V(n,e)}, \]

    where M is the maxmimum number of codewords and $V(n,e)$ is equal to the
    contents of a ball of radius e. This bound is useful for small values of d.
    Codes for which equality holds are called {\it perfect}.


    EXAMPLES:
        sage: hamming_upper_bound(10,2,3)
        93

    """

    ans=gap.eval("UpperBoundHamming(%s,%s,%s)"%(n,d,q))
    return QQ(ans)

def singleton_upper_bound(n,q,d):
    r"""
    Returns the Singleton upper bound for number of elements in
    the largest code of minimum distance d in $\FF_q^n$. Wraps GAP's
    UpperBoundSingleton.

    This bound is based on the shortening of codes. By shortening an
    $(n, M, d)$ code d-1 times, an $(n-d+1,M,1)$ code results, with
    $M \leq q^n-d+1$. Thus

    \[ M \leq q^{n-d+1}. \]

    Codes that meet this bound are called {\it maximum distance separable} (MDS).

    EXAMPLES:
        sage: singleton_upper_bound(10,2,3)
        256

    """

    ans=gap.eval("UpperBoundSingleton(%s,%s,%s)"%(n,d,q))
    return QQ(ans)

def gv_info_rate(n,delta,q):
    """
    GV lower bound for information rate of a q-ary code of length n
    minimum distance delta*n

    EXAMPLES:
        sage: gv_info_rate(100,1/4,3)
        0.367049926083

    """
    q = ZZ(q)
    ans=log(gilbert_lower_bound(n,q,int(n*delta)),q)/n
    return ans

def entropy(x,q):
    """
    Computes the entropy on the q-ary symmetric channel.

    """
    q = ZZ(q)
    H = x*log(q-1,q)-x*log(x,q)-(1-x)*log(1-x,q)
    return H

def gv_bound_asymp(delta,q):
    """
    Computes the asymptotic GV bound for the information rate, R.

    EXAMPLES:
        sage: f = lambda x: gv_bound_asymp(x,2)
        sage: plot(f,0,1)
    """
    return (1-entropy(delta,q))


def hamming_bound_asymp(delta,q):
    """
    Computes the asymptotic Hamming bound for the information rate.

    EXAMPLES:
        sage: f = lambda x: hamming_bound_asymp(x,2)
        sage: plot(f,0,1)
    """
    return (1-entropy(delta/2,q))

def singleton_bound_asymp(delta,q):
    """
    Computes the asymptotic Singleton bound for the information rate.

    EXAMPLES:
        sage: f = lambda x: singleton_bound_asymp(x,2)
        sage: plot(f,0,1)
    """
    return (1-delta)

def plotkin_bound_asymp(delta,q):
    """
    Computes the asymptotic Plotkin bound for the information rate,
    provided 0 < delta < 1-1/q.

    EXAMPLES:
        sage: plotkin_bound_asymp(1/4,2)
        1/2
    """
    r = 1-1/q
    return (1-delta/r)

def elias_bound_asymp(delta,q):
    """
    Computes the asymptotic Elias bound for the information rate,
    provided 0 < delta < 1-1/q.

    EXAMPLES:
        sage: elias_bound_asymp(1/4,2)
        0.39912396330...
    """
    r = 1-1/q
    return RDF((1-entropy(r-sqrt(r*(r-delta)), q)))

def mrrw1_bound_asymp(delta,q):
    """
    Computes the ``first asymptotic McEliese-Rumsey-Rodemich-Welsh
    bound for the information rate, provided 0 < delta < 1-1/q.

    EXAMPLES:
        sage: mrrw1_bound_asymp(1/4,2)
        0.354578902665

    """
    return RDF(entropy((q-1-delta*(q-2)-2*sqrt((q-1)*delta*(1-delta)))/q,q))





