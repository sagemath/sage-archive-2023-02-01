"""
Linear feedback shift register (LFSR) sequence commands


Stream ciphers have been used for
a long time as a source of pseudo-random number generators.

S. Golomb [G] gives a list of three
statistical properties a
sequence of numbers ${\bf a}=\{a_n\}_{n=1}^\infty$,
$a_n\in \{0,1\}$, should display to be considered
``random''. Define the {\bf autocorrelation} of
${\bf a}$ to be
\[
C(k)=C(k,{\bf a})=\lim_{N\rightarrow \infty}
{1\over N}\sum_{n=1}^N (-1)^{a_n+a_{n+k}}.
\]
In the case where ${\bf a}$ is periodic with
period $P$ then this reduces to
\[
C(k)={1\over P}\sum_{n=1}^P (-1)^{a_n+a_{n+k}}.
\]
Assume ${\bf a}$ is periodic with period $P$.
\begin{itemize}
\item[] {\bf balance}: $|\sum_{n=1}^P(-1)^{a_n}|\leq 1$.
\item[] {\bf low autocorrelation}:
\[
C(k)=
\left\{
\begin{array}{cc}
1,& k=0,\\
\epsilon, & k\not= 0.
\end{array}
\right.
\]
(For sequences satisfying these first two properties,
it is known that $\epsilon=-1/P$ must hold.)
\item[] {\bf proportional runs property}:
In each period, half the runs have length $1$,
one-fourth have length $2$, etc. Moveover, there
are as many runs of $1$'s as there are of
$0$'s.
\end{itemize}

A {\bf general feedback shift register} is a map
$f:{\bf F}_q^d\rightarrow {\bf F}_q^d$
of the form
\[
\begin{array}{c}
f(x_0,...,x_{n-1})=(x_1,x_2,...,x_n),\\
x_n=C(x_0,...,x_{n-1}),
\end{array}
\]
where $C:{\bf F}_q^d\rightarrow {\bf F}_q$ is a given
function. When $C$ is of the form
\[
C(x_0,...,x_{n-1})=a_0x_0+...+a_{n-1}x_{n-1},
\]
for some given constants $a_i\in {\bf F}_q$, the
map is called a {\bf linear feedback shift register
(LFSR)}.

{\bf Example of a LFSR} Let
\[
f(x)=a_{{0}}+a_{{1}}x+...+a_{{n}}{x}^n+...,
\]
\[
g(x)=b_{{0}}+b_{{1}}x+...+b_{{n}}{x}^n+...,
\]
be given polynomials in ${\bf F}_2[x]$ and let
\[
h(x)={f(x)\over g(x)}=c_0+c_1x+...+c_nx^n+... \ .
\]
We can compute a recursion formula which allows us to rapidly compute
the coefficients of $h(x)$ (take $f(x)=1$):
\[
c_{n}=\sum_{i=1}^n {{-b_i\over b_0}c_{n-i}}.
\]

The coefficients of $h(x)$ can, under certain conditions on
$f(x)$ and $g(x)$, be considered ``random'' from certain statistical
points of view.

{\bf Example}:
For instance, if
\[
f(x)=1,\ \ \ \ g(x)=x^4+x+1,
\]
then
\[
h(x)=1+x+x^2+x^3+x^5+x^7+x^8+...\ .
\]
The coefficients of $h$ are
\[
\begin{array}{c}
1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, \\
1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, ...\ .
\end{array}
\]
The sequence of $0,1$'s is periodic with period
$P=2^4-1=15$ and satisfies Golomb's three
randomness conditions.
However, this sequence of period 15 can be ``cracked''
(i.e., a procedure to reproduce $g(x)$)
by knowing only 8 terms! This is the function of the
Berlekamp-Massey algorithm [M], implemented as
\code{berlekamp_massey.py}.

[G] Solomon Golomb, {\bf Shift register sequences},
Aegean Park Press, Laguna Hills, Ca,
1967

[M] James L. Massey, ``Shift-Register Synthesis and BCH Decoding.''
IEEE Trans. on Information Theory, vol. 15(1), pp. 122-127, Jan 1969.

AUTHOR:
    -- Timothy Brock

Created 11-24-2005 by wdj. Last updated 12-02-2005.
"""

###########################################################################
#  Copyright (C) 2006 Timothy Brock
#  and William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
###########################################################################

import copy

from sage.structure.all import Sequence
from sage.rings.all import is_FiniteField

def lfsr_sequence(key, fill, n):
    r"""
    This function creates an lfsr sequence.

    INPUT:
        key -- a list of finite field elements, [c_0,c_1,...,c_k].
        fill -- the list of the initial terms of the lfsr sequence,
                [x_0,x_1,...,x_k].
        n -- number of terms of the sequence that the
             function returns.

    OUTPUT:
        The lfsr sequence defined by $x_{n+1} = c_kx_n+...+c_0x_{n-k}$,
        for $n \leq k$.

    EXAMPLES:
        sage: F = GF(2); l = F(1); o = F(0)
        sage: F = GF(2); S = LaurentSeriesRing(F,'x'); x = S.gen()
        sage: fill = [l,l,o,l]; key = [1,o,o,l]; n = 20
        sage: L = lfsr_sequence(key,fill,20); L
        [1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0]
        sage: g = berlekamp_massey(L); g
        x^4 + x^3 + 1
        sage: (1)/(g.reverse()+O(x^20))
        1 + x + x^2 + x^3 + x^5 + x^7 + x^8 + x^11 + x^15 + x^16 + x^17 + x^18 + O(x^20)
        sage: (1+x^2)/(g.reverse()+O(x^20))
        1 + x + x^4 + x^8 + x^9 + x^10 + x^11 + x^13 + x^15 + x^16 + x^19 + O(x^20)
        sage: (1+x^2+x^3)/(g.reverse()+O(x^20))
        1 + x + x^3 + x^5 + x^6 + x^9 + x^13 + x^14 + x^15 + x^16 + x^18 + O(x^20)
        sage: fill = [l,l,o,l]; key = [l,o,o,o]; n = 20
        sage: L = lfsr_sequence(key,fill,20); L
        [1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1]
        sage: g = berlekamp_massey(L); g
        x^4 + 1
        sage: (1+x)/(g.reverse()+O(x^20))
        1 + x + x^4 + x^5 + x^8 + x^9 + x^12 + x^13 + x^16 + x^17 + O(x^20)
        sage: (1+x+x^3)/(g.reverse()+O(x^20))
        1 + x + x^3 + x^4 + x^5 + x^7 + x^8 + x^9 + x^11 + x^12 + x^13 + x^15 + x^16 + x^17 + x^19 + O(x^20)

    AUTHOR:
        -- Timothy Brock (11-2005, with code modified from Python Cookbook,
            \url{http://aspn.activestate.com/ASPN/Python/Cookbook/}
    """
    if not isinstance(key, list):
        raise TypeError, "key must be a list"
    key = Sequence(key)
    F = key.universe()
    if not is_FiniteField(F):
        raise TypeError, "universe of sequence must be a finite field"

    s = fill
    k = len(fill)
    L = []
    for i in range(n):
        s0 = copy.copy(s)
        L.append(s[0])
        s = s[1:k]
        s.append(sum([key[i]*s0[i] for i in range(k)]))
    return L

