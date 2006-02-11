"""nodoctest
Hyperelliptic curves over a finite field

TOTALLY BROKEN!!

AUTHOR: 2005-11-13, David Joyner <wdj@usna.edu>
AUTHOR: 2005-11-13, William Stein <wstein@ucsd.edu>
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#                     2005 David Joyner <wdj@usna.edu>
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

import hyperell_generic
from sage.schemes.plane_curves.projective_curve import ProjectiveCurve_generic
from sage.rings.all import (factor, FractionField, MPolynomial, MPolynomialRing, degree_lowest_rational_function)
from sage.misc.all import add, mul


class HyperellipticCurve_finite_field(hyperell_generic.HyperellipticCurve_generic,
                                      ProjectiveCurve_generic):
    def riemann_roch_space(self, D, E):
        r"""
        Return a maximal subbasis $C$ of $L(E)$ no non-constant
        element of which is contained in $L(D)$.

        Using \code{self.riemann_roch_space(D,E)} and
        \code{self.riemann_roch_space([0,...,0],D)}, one can compute
        bases $B$ of $L(D)$ and $C$ of $L(E)$ such that $B\subset C$.

        This is an implementation of the pivot method, described in
        Joyner-Ksir, {\em Bases for Riemann-Roch spaces of hyperelliptic
        curves: Effective case}.

        Let $X$ denote a hyperelliptic curve over an algebraically
        closed field $F$ and let $D<E$ denote effective divisors on
        $X$.  See the docstring of the function \code{rr_space1_3} for
        more details of the algorithm.


        INPUT:
            self -- a hyperelliptic curve defined by a polynomial
                    equation f = y^2 - h(x) = 0
                    over a prime finite field F = GF(p)
                    having n F-rational points (see the SAGE function
                    places_on_curve)
           D,E   -- n-tuples of non-negative integers $(d1, ..., dn)$,
                    $(e1, ..., en)$ that represent
                    effective divisors $Div(D) = d1*P1+...+dn*Pn$,
                    $Div(E) = e1*P1+...+en*Pn$, where
                    $di\leq ei$ and $X(F) = \{P1,...,Pn\}$.
                    **The ordering is that dictated by places_on_curve.**

        OUTPUT:
           list --  A sub-basis of L(Div(E)) no element of which is
                    contained in $L(D)$.  Moreover, every element of
                    this basis is of the form a(x)+y*b(x), for some
                    rational functions a(x), b(x).

        EXAMPLES:
            sage: F = GF(5)
            sage: R = MPolynomialRing(F,2, names = ['x','y'])
            sage: x,y = R.gens()
            sage: f = y**2 - x**9 - x
            sage: C = HyperellipticCurve(f)
            sage: C.places_on_curve()
                  [(-2, -1, 1), (-2, 1, 1), (0, 0, 1), (0, 1, 0), (2, -2, 1), (2, 2, 1)]
            sage: D = [0,0,0,0,0,0]
            sage: E = [0,4,2,0,0,0]
            sage: b = C.riemann_roch_space(D,E)
            sage: len(b)
                  2
            sage: C.divisor_of_function(b[0])
                  [[-4, (-2, -1, 1)], [-2, (0, 0, 1)]]
            sage: C.divisor_of_function(b[1])
                  [[-4, (-2, -1, 1)], [-1, (0, 0, 1)]]
            sage: D = [0,4,2,0,0,0]
            sage: E = [0,6,3,0,0,0]
            sage: b = C.riemann_roch_space(D,E)
            sage: len(b)
                  3
            sage: C.divisor_of_function(b[0])
                  [[-6, (-2, -1, 1)], [-3, (0, 0, 1)]]
            sage: C.divisor_of_function(b[1])
                  [[-6, (-2, -1, 1)], [-2, (0, 0, 1)]]
            sage: C.divisor_of_function(b[2])
                  [[-5, (-2, -1, 1)], [-2, (0, 0, 1)]]
        """
        Ecopy = list(E)
        f = self.polynomial()
        R = f.parent()
        if f.coefficient(R.gen(1)) != 0:
            raise NotImplementedError, "The defining equation (=%s) must be of the form y^2 - h(x) = 0"%f
        F = self.base_ring()
        p = F.characteristic()
        x,y = R.gens()
        deg = self.polynomial().degree(x)
        genus = self.geometric_genus()
        dim = sum(Ecopy)+1-genus
        suppDE = []
        suppDEidx = []
        pts = self.places_on_curve()
        for i in range(len(Ecopy)):
            if E[i]>D[i]:
                for j in range(Ecopy[i]-D[i]):
                    suppDE.append([pts[i],i])
                suppDEidx.append(i)
        basis = []
        while (dim>1 and len(suppDE)>0):
            #print Ecopy        ###### comment this out
            b = rr_space1_3(f,D,Ecopy)
            basis.append(b[0])
            PP = suppDE.pop()
            P = PP[0]
            idxP = PP[1]
            Ecopy[idxP] = Ecopy[idxP]-1
            dim = sum(Ecopy)+1-genus
        return basis



def rr_space1_3(f,D,E):
    r"""
    Compute and return a Riemann-Roch space.

    The hyperelliptic curve $X$ is defined by $f(x,y) = y^2-h(x) = 0$,
    h(x) is a polynomial over a finite field F.  $D \leq E$ are
    effective divisors supported on F-rational (degree 1) places on
    $X$. The divisors $D$ and $E$ are represented by

    Let $E = \sum_i e_i P_i, D = \sum_i d_i P_i$
        $L = \{ i \ |\ e_i > d_i \} = \{ i_0 < i_1 < ... < i_k \}$

    \begin{enumerate}
    \item Pick $j = i_0$, pivot pt = $P_j^*$,
       seed fcn = $f_0(x,y) = (y+y_j)/\prod_i (x-x_i)$,
       pivot functions = $f_{j,a}(x,y) = (x - x_j)^{-a}$

    \item Compute local coords x = x(t), y = y(t) about pivot point

    \item Compute $c_1,...,c_{e_j}$ such that
          $f^1_0(x,y)= f_0(x,y)-c_{e_j}f_{j,e_j}(x,y)-...-c_1f_{j,1}(x,y)$
       has not pole at the pivot pt $P_j^*$

    \item Replace j by next element of L, replace $f_0(x,y)$ by
       $f^1_0(x,y)$, and redefine the pivot fcns.
       Repeat steps 2, 3.

    \item Repeat step 4 until elements of L are exhausted and all poles at
       $P_{i_0}^*, ..., P_{i_k}^*$ are  removed. Add the resulting fcn
       $f^k_0(x,y)$ to the basis.

    \item Replace E by E-P, where P in supp(E)-supp(D) is arbitrary.

    \item Repeat 1-6 until E = D.
    \end{enumerate}

    Steps 1--3 are implemented in this function. The remaining steps
    in the riemann_roch_space function above.

    EXAMPLES:

        sage: F = GF(5)
        sage: R = MPolynomialRing(F,2, names = ['x','y'])
        sage: x,y = R.gens()
        sage: f = y**2 - x**9 - x
        sage: C = HyperellipticCurve(f)
        sage: C.places_on_curve()
              [(-2, -1, 1), (-2, 1, 1), (0, 0, 1), (0, 1, 0), (2, -2, 1), (2, 2, 1)]
        sage: D = [0,0,0,0,0,0]
        sage: E = [0,4,0,0,0,0]
        sage: rr_space1_3(f,D,E)
        ((1 + y + x^2 + 2*x^3)/(1 + 2*x + 4*x^2 + 3*x^3 + x^4),
        [(1 + y)/(1 + 2*x + 4*x^2 + 3*x^3 + x^4),
        3/(1 + 2*x + 4*x^2 + 3*x^3 + x^4),
        4/(4 + 4*x + x^2),
        2/(2 + x)])

    """
    import constructor
    F = f.base_ring()
    C = constructor.HyperellipticCurve(f)
    pts = C.places_on_curve()
    numpts = len(pts)
    R = f.parent()
    x,y = R.gens()
    genus = C.geometric_genus()
    dim = sum(E)+1-genus
    suppE = []
    suppEidx = []
    for i in range(numpts):
        if E[i]>0:
            suppE.append(pts[i])
            suppEidx.append(i)
    ptindx = min(suppEidx)
    if pts[ptindx][2] != F(1):
        print "Point choosen: ",pts[ptindx]
        raise NotImplementedError
    pt0 = [pts[ptindx][0],pts[ptindx][1]]
    pivot_pt = [pts[ptindx][0],-pts[ptindx][1]]
    lcs = C.local_coordinates(pt0,5)
    yt = lcs[1]
    xt = lcs[0]
    f0 = (y+pt0[1])/mul([(x-pts[i][0])**E[i] for i in range(numpts)])
    f1 = [1/(x-pts[ptindx][0])**j for j in range(1,E[ptindx]+1)]
    f1.reverse()
    R0 = MPolynomialRing(F,3,names = [str(x),str(y),"t"])
    vars0 = R0.gens()
    #print R0,vars0
    t = vars0[2]
    g = f0
    glist = [f0]
    coeff0 = degree_lowest_rational_function(g(xt,yt),t)[1]
    for j in range(1,E[ptindx]+1):
        #print "j,f1,range: ",j,"  ",f1,"      ",range(1,E[ptindx])
        h = f1[j-1]
        ldg = degree_lowest_rational_function(g(xt,yt),t)
        d0 = ldg[0]
        ldh = degree_lowest_rational_function(h(xt,yt),t)
        d1 = ldh[0]
        #print "j,d0,d1,h= ",j,"  ",d0,"  ", d1,"  ",h
        if d1==0:
            return g
        if d1==d0:
            coeff1 = ldh[1]
            fcn_matching_pole = h*ldg[1]/coeff1
            g = g-fcn_matching_pole
            glist.append(-fcn_matching_pole)
    return g,glist

