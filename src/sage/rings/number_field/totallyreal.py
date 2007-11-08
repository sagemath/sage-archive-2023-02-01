"""
Enumeration of Totally Real Fields

AUTHORS:
    -- Craig Citro and John Voight (2007-11-04):
        * Additional doctests and type checking.
    -- John Voight (2007-10-27):
        * Separated DSage component.
    -- John Voight (2007-10-17):
        * Added pari functions to avoid recomputations.
    -- John Voight (2007-10-09):
        * Added DSage module.
    -- John Voight (2007-09-19):
        * Various optimization tweaks.
    -- John Voight (2007-09-01):
        * Initial version.
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein and John Voight
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

import math, sys, bisect

from sage.rings.polynomial.polynomial_ring import PolynomialRing
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.libs.pari.gen import pari

ZZx = PolynomialRing(IntegerRing(), 'x')

from sage.rings.number_field.totallyreal_data import tr_data, int_has_small_square_divisor

def odlyzko_bound_totallyreal(n):
    r"""
    This function returns the unconditional Odlyzko bound for the root
    discriminant of a totally real number field of degree n.

    NOTE:
    The bounds for n > 50 are not necessarily optimal.

    INPUT:
    n -- integer, the degree

    OUTPUT:
    a lower bound on the root discriminant

    EXAMPLES:
    sage: [sage.rings.number_field.totallyreal.odlyzko_bound_totallyreal(n) for n in range(1,5)]
    [1.0, 2.2229999999999999, 3.6099999999999999, 5.0670000000000002]

    AUTHORS:
    - John Voight (2007-09-03)

    NOTES:
    The values are calculated by Martinet [M].

        [M] Jacques Martinet, Petits discriminants des corps de nombres, Journ. Arithm. 1980,
            Cambridge Univ. Press, 1982, 151--193.
    """

    if n <= 10:
        dB = [1.,2.223,3.610,5.067,6.523,7.941,9.301,10.596,11.823,12.985][n-1]
    elif n <= 20:
        dB = [14.0831,15.1217,16.1047,17.0359,17.9192,18.7580,19.5555,20.3148,21.0386,21.7294][n-11]
    elif n <= 30:
        dB = [22.3896,23.0212,23.6261,24.2061,24.7628,25.2976,25.2976,26.3071,26.3071,27.2440][n-21]
    elif n <= 40:
        dB = [27.2440,28.1165,28.1165,28.9315,28.9315,29.6948,29.6948,30.4117,30.4117,31.0865][n-31]
    elif n <= 50:
        dB = [31.0865,31.7232,31.7232,32.3252,32.3252,32.8954,32.8954,33.4365,33.4365,33.9508][n-41]
    else:
        dB = 33.9508
    return dB

def enumerate_totallyreal_fields(n, B, a = [], verbose=0, return_seqs=False, phc=False):
    r"""
    This function enumerates (primitive) totally real fields of
    degree $n>1$ with discriminant $d \leq B$; optionally one can
    specify the first few coefficients, where the sequence $a$
    corresponds to a polynomial by
        $$ a[d]*x^n + ... + a[0]*x^(n-d) $$
    if length(a) = d+1, so in particular always a[d] = 1.
    If verbose == 1 (or 2), then print to the screen (really) verbosely; if
    verbose is a string, then print verbosely to the file specified by verbose.
    If return_seqs, then return the polynomials as sequences (for easier
    exporting to a file).

    NOTE:
    This is guaranteed to give all primitive such fields, and
    seems in practice to give many imprimitive ones.

    INPUT:
    n -- integer, the degree
    B -- integer, the discriminant bound
    a -- list (default: []), the coefficient list to begin with
    verbose -- boolean or string (default: False)
    phc -- boolean or integer (default: False)

    OUTPUT:
    the list of fields with entries [d,f], where
      d is the discriminant and f is a defining polynomial,
    sorted by discriminant.

    EXAMPLES:
    In this first simple example, we compute the totally real quadratic
    fields of discriminant <= 50.

    sage: enumerate_totallyreal_fields(2,50)
    [[5, x^2 - x - 1],
     [8, x^2 - 2],
     [12, x^2 - 3],
     [13, x^2 - x - 3],
     [17, x^2 - x - 4],
     [21, x^2 - x - 5],
     [24, x^2 - 6],
     [28, x^2 - 7],
     [29, x^2 - x - 7],
     [33, x^2 - x - 8],
     [37, x^2 - x - 9],
     [40, x^2 - 10],
     [41, x^2 - x - 10],
     [44, x^2 - 11],
     [5, x^2 - 3*x + 1]]
    sage: [ d for d in range(5,50) if (is_squarefree(d) and d%4 == 1) or (d%4 == 0 and is_squarefree(d/4)) ]
    [5, 8, 12, 13, 17, 20, 21, 24, 28, 29, 33, 37, 40, 41, 44]

    Next, we compute all totally real quintic fields of discriminant <= 10^5.

    sage: enumerate_totallyreal_fields(5,10^5)
    [[14641, x^5 - x^4 - 4*x^3 + 3*x^2 + 3*x - 1],
     [24217, x^5 - 5*x^3 - x^2 + 3*x + 1],
     [36497, x^5 - 2*x^4 - 3*x^3 + 5*x^2 + x - 1],
     [38569, x^5 - 5*x^3 + 4*x - 1],
     [65657, x^5 - x^4 - 5*x^3 + 2*x^2 + 5*x + 1],
     [70601, x^5 - x^4 - 5*x^3 + 2*x^2 + 3*x - 1],
     [81509, x^5 - x^4 - 5*x^3 + 3*x^2 + 5*x - 2],
     [81589, x^5 - 6*x^3 + 8*x - 1],
     [89417, x^5 - 6*x^3 - x^2 + 8*x + 3]]

    We see that there are 9 such fields (up to isomorphism!).

    NOTES:
    This function uses Hunter's algorithm [C, Section 9.3] and
    modifications due to Takeuchi [T] and the author (not yet published).

    We enumerate polynomials
        f(x) = x^n + a[n-1]*x^(n-1) + ... + a[0].
    Hunter's theorem gives bounds on a[n-1] and a[n-2];
    then given a[n-1] and a[n-2], one can recursively compute bounds on
    a[n-3], ..., a[0] using the fact that the polynomial is totally real
    by looking at the zeros of successive derivatives and applying
    Rolle's theorem!  See [T] for more details.

        REFERENCES:
            [C] Henri Cohen, Advanced topics in computational number
                theory, Graduate Texts in Mathematics, vol. 193,
                Springer-Verlag, New York, 2000.
            [T] Kisao Takeuchi, Totally real algebraic number fields of
                degree 9 with small discriminant, Saitama Math. J.
                17 (1999), 63--85 (2000).

    AUTHORS:
    - John Voight (2007-09-03)
    """

    if not isinstance(n, Integer):
        try:
            n = Integer(n)
        except:
            raise TypeError, "cannot coerce n (= %s) to an integer"%n
    if (n < 1):
        raise ValueError, "n must be at least 1."

    # Initialize
    T = tr_data(n,B,a)
    S = []
    dB_odlyzko = odlyzko_bound_totallyreal(n)
    dB = math.ceil(40000*dB_odlyzko**n)
    counts = [0,0,0,0]

    # Trivial case
    if n == 1:
        if return_seqs:
            return [[0,0,0,0],[[1,[-1,1]]]]
        else:
            return [[1,pari('x-1')]]

    if verbose:
        saveout = sys.stdout
        if type(verbose) == str:
            fsock = open(verbose, 'w')
            sys.stdout = fsock
        # Else, print to screen
    f_out = [0]*n + [1]
    if verbose == 2:
        T.incr(f_out,verbose,phc=phc)
    else:
        T.incr(f_out,phc=phc)

    while f_out[n] <> 0:
        if verbose:
            print "==>", f_out,

        nf = pari(str(f_out)).Polrev()
        d = nf.poldisc()
        counts[0] += 1
        if d > 0 and nf.polsturm_full() == n:
            da = int_has_small_square_divisor(Integer(d))
            if d > dB or d <= B*da:
                counts[1] += 1
                if nf.polisirreducible():
                    counts[2] += 1
                    [zk,d] = nf.nfbasis_d()

                    if d <= B:
                        if verbose:
                            print "has discriminant", d,

                        # Find a minimal lattice element
                        counts[3] += 1
                        ng = pari([nf,zk]).polredabs()

                        # Check if K is contained in the list.
                        found = False
                        ind = bisect.bisect_left(S, [d,ng])
                        while ind < len(S) and S[ind][0] == d:
                            if S[ind][1] == ng:
                                if verbose:
                                    print "but is not new"
                                found = True
                                break
                            ind += 1
                        if not found:
                            if verbose:
                                print "and is new!"
                            S.insert(ind, [d,ng])

                    else:
                        if verbose:
                            print "has discriminant", abs(d), "> B"
                else:
                    if verbose:
                        print "is not irreducible"
            else:
                if verbose:
                    print "has discriminant", abs(d), "with no large enough square divisor"
        else:
            if verbose:
                if d == 0:
                    print "is not squarefree"
                else:
                    print "is not totally real"

        if verbose == 2:
            T.incr(f_out,verbose=verbose,phc=phc)
        else:
            T.incr(f_out,phc=phc)

    # In the application of Smyth's theorem above, we exclude finitely
    # many possibilities which we must now throw back in.
    if n == 2 and B >= 5:
        S = [[5,pari('x^2-3*x+1')]] + S
    elif n == 3 and B >= 49:
        S = [[49,pari('x^3-5*x^2+6*x-1')]] + S
    # The polynomials with n = 4 define imprimitive number fields.

    # Now check for isomorphic fields
    weed_fields(S)

    # Output.
    if verbose:
        print "="*80
        print "Polynomials tested:", counts[0]
        print "Polynomials with sssd poldisc:", counts[1]
        print "Irreducible polynomials:", counts[2]
        print "Polynomials with nfdisc <= B:", counts[3]
        for i in range(len(S)):
            print S[i]
        if type(verbose) == str:
            fsock.close()
        sys.stdout = saveout

    if return_seqs:
        return [counts,[[s[0],s[1].reverse().Vec()] for s in S]]
    else:
        return S

def weed_fields(S):
    r"""
    Function used internally by the enumerate_totallyreal_fields()
    routine.  (Weeds the fields listed by [discriminant, polynomial]
    for isomorphism classes.)

    EXAMPLES:
        sage: ls = [[5,pari('x^2-3*x+1')],[5,pari('x^2-5')]]
        sage: sage.rings.number_field.totallyreal.weed_fields(ls); ls
        [[5, x^2 - 3*x + 1]]
    """
    i = 0
    if len(S) == 0:
       return
    n = len(S[0][1])-1
    while i < len(S)-1:
       j = i+1
       while j < len(S) and S[i][0] == S[j][0]:
           if S[i][1].nfisisom(S[j][1]):
               # Keep the one with a smallest T_2
               T_2i = S[i][1][n-1]**2 - 2*S[i][1][n-2]
               T_2j = S[j][1][n-1]**2 - 2*S[j][1][n-2]
               if T_2i <= T_2j:
                   S.pop(j)
               else:
                   s = S.pop(i)
                   S.insert(i, S.pop(j))
           else:
               j += 1
       i += 1

def timestr(m):
    r"""
    Converts seconds to a human-readable time string.

    INPUT:
        m -- integer, number of seconds

    OUTPUT:
        The time in days, hours, etc.

    EXAMPLES:
        sage: sage.rings.number_field.totallyreal.timestr(3765)
        '1h 2m 45.0s'
    """

    n = math.floor(m)
    p = m-n
    outstr = ''
    if m >= 60*60*24:
        t = n//(60*60*24)
        outstr += str(t)[:len(str(t))-2] + 'd '
        n -= t*(60*60*24)
    if m >= 60*60:
        t = n//(60*60)
        outstr += str(t)[:len(str(t))-2] + 'h '
        n -= t*(60*60)
    if m >= 60:
        t = n//60
        outstr += str(t)[:len(str(t))-2] + 'm '
        n -= t*60
    n += p
    outstr += '%.1f'%n + 's'

    return outstr

def __selberg_zograf_bound(n, g):
    r"""
    Returns an upper bound on the possible root discriminant of a
    totally real field of degree n which gives rise to an arithmetic
    Fuchsian group of genus g.  The bound is:
       (16/3*(g+1))^(2/(3*n))*(2*pi)^(4/3).

    INPUT:
        n -- integer, the degree
        g -- integer, the genus

    OUTPUT:
        the upper bound.

    AUTHORS:
        - John Voight (2007-09-19)

    EXAMPLES:
        sage: sage.rings.number_field.totallyreal.__selberg_zograf_bound(8,7)
        15.851871776151311
    """
    return ((16./3)*(g+1))**(2./(3*n))*(2*3.1415926535897931)**(4./3)
