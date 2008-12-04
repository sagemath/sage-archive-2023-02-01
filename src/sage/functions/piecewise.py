r"""
Piecewise-defined Functions.

\sage implements a very simple class of piecewise-defined functions.  Functions
may be any type of symbolic expression.  Infinite intervals are not supported.
The endpoints of each interval must line up.

Implemented methods:
    latex outout
    __call__
    __eq__
    plot (using matplotlib)
    integral
    convolution
    derivative
    critical points
    tangent_line
    trapezoid
    trapezoid_integral_approximation
    riemann_sum   (right and left end points)
    riemann_sum_integral_approximation (right and left end points)
    fourier series
    fourier series
       value
       coefficients (also sine series and cosine series, and Cesaro sum)
       partial sum (in string format)
       plot of partial sum
       plot_fourier_series_partial_sum_cesaro
       plot_fourier_series_partial_sum_hann
       plot_fourier_series_partial_sum_filter
    laplace transform
       latex output option
    domain
    range
    list
    __add__  - addition (of functions)
    __mul__  - multiplication (of functions, or fcn*scalar - ie, *right* multiplication by QQ)
    __rmul__ - *left* multiplication by QQ
    extend_by_zero_to
    unextend

TODO:
   [] Implement (a) max/min location and values,
   [] (b) multiplication by a scalar.
   [] (c) Extend the implementation of the trick to pass \sage's pi back
      and forth with Maxima's %pi to other constants (e, for example)
      [[Passing the constants to maxima is already implemented; maybe
       need passing them back?]]
   [] Need: parent object -- ring of piecewise functions
   [] This class should derive from an element-type class, and should
      define _add_, _mul_, etc.   That will automatically take care
      of left multiplication and proper coercion.
      The coercion mentioned below for scalar mult on right is
      bad, since it only allows ints and rationals.  The right
      way is to use an element class and only define _mul_, and
      have a parent, so anything gets coerced properly.

AUTHOR: David Joyner (2006-04) -- initial version
        DJ (2006-09) -- added __eq__, extend_by_zero_to, unextend, convolution,
                        trapezoid, trapezoid_integral_approximation, riemann_sum,
                        riemann_sum_integral_approximation, tangent_line
                        fixed bugs in __mul__, __add__
	DJ (2007-03) -- adding Hann filter for FS, added general FS filter methods
	                for computing and plotting, added options to plotting of FS
	                (eg, specifying rgb values arenow allowed). Fixed bug in
                        documentation reported by Pablo De Napoli.
	DJ (2007-09) - bug fixes due to behaviour of SymbolicArithmetic
	DJ (2008-04) - fixed docstring bugs reported by J Morrow; added support for
                       laplace transform of functions with infinite support.
        DJ (2008-07) - fixed a left multiplication bug reported by C. Boncelet
                       (by defining __rmul__ = __mul__).

TESTS:
    sage: R.<x> = QQ[]
    sage: f = Piecewise([[(0,1),1*x^0]])
    sage: 2*f
    Piecewise defined function with 1 parts, [[(0, 1), 2]]

"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdjoyner@gmail.com>
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

import sage.plot.plot
import sage.interfaces.all
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField
from sage.rings.real_mpfr import RealField_constructor as RealField
from sage.misc.sage_eval import sage_eval
from sage.rings.all import QQ, RR, Integer, Rational
from sage.calculus.functional import derivative

from sage.calculus.calculus import SR, var, maxima

def meval(x):
    from sage.calculus.calculus import symbolic_expression_from_maxima_element
    return symbolic_expression_from_maxima_element(maxima(x))

def piecewise(list_of_pairs):
    """
    Returns a piecewise function from a list of (interval, function)
    pairs.

    \code{list_of_pairs} is a list of pairs (I, fcn), where fcn is
    a SAGE function (such as a polynomial over RR, or functions
    using the lambda notation), and I is an interval such as I = (1,3).
    Two consecutive intervals must share a common endpoint.

    We assume that these definitions are consistent (ie, no checking is
    done).

    EXAMPLES:
        sage: f1 = lambda x:-1
        sage: f2 = lambda x:2
        sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
        sage: f(1)
        -1
        sage: f(3)
        2
    """
    return PiecewisePolynomial(list_of_pairs)

Piecewise = piecewise

class PiecewisePolynomial:
    """
    Returns a piecewise function from a list of (interval, function)
    pairs.

    EXAMPLES:
        sage: f1 = lambda x:-1
        sage: f2 = lambda x:2
        sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
        sage: f(1)
        -1
        sage: f(3)
        2
    """
    def __init__(self, list_of_pairs):
        r"""
        \code{list_of_pairs} is a list of pairs (I, fcn), where fcn is
        a SAGE function (such as a polynomial over RR, or functions
        using the lambda notation), and I is an interval such as I = (1,3).
        Two consecutive intervals must share a common endpoint.

        We assume that these definitions are consistent (ie, no checking is
        done).
        """
        self._length = len(list_of_pairs)
        self._intervals = [x[0] for x in list_of_pairs]
        self._functions = [x[1] for x in list_of_pairs]
        self._list = list_of_pairs

    def list(self):
        return self._list

    def length(self):
        return self._length

    def __repr__(self):
        return 'Piecewise defined function with %s parts, %s'%(
            self.length(),self.list())

    def _latex_(self):
	r"""
	EXAMPLES:
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: latex(f)
            \begin{cases}
            x \ {\mapsto}\ 1 &\text{on $(0, 1)$}\cr
            x \ {\mapsto}\ 1 - x &\text{on $(1, 2)$}
            \end{cases}

            sage: f(x) = sin(x*pi/2)
            sage: g(x) = 1-(x-1)^2
            sage: h(x) = -x
            sage: P = Piecewise([[(0,1), f], [(1,3),g], [(3,5), h]])
            sage: latex(P)
            \begin{cases}
            x \ {\mapsto}\ \sin \left( \frac{{\pi x}}{2} \right) &\text{on $(0, 1)$}\cr
            x \ {\mapsto}\ 1 - {\left( x - 1 \right)}^{2}  &\text{on $(1, 3)$}\cr
            x \ {\mapsto}\ -x &\text{on  $(3, 5)$}
            \end{cases}
	"""
        intervals = self.intervals()
        funcs = [p[1] for p in self.list()]
        tex = ['\\begin{cases}\n']
        for i in range(len(funcs)):
            f = funcs[i]
            # left endpoint
            a = intervals[i][0]
            # right endpoint
            b = intervals[i][1]
            tex.append(r'%s &\text{on $(%s, %s)$}' % (f._latex_(), a, b))
            tex.append('\\cr\n')
        # remove the last TeX newline
        tex = tex[:-1]
        tex.append('\n\\end{cases}')
        return ''.join(tex)

    def intervals(self):
        """
	A piecewise non-polynomial example.

        EXAMPLES:
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f.intervals()
            [(0, 1), (1, 2), (2, 3), (3, 10)]
        """
        return self._intervals

    def domain(self):
        """
        Returns the domain of the function.

        EXAMPLES:
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f.domain()
            (0, 10)
        """
        endpoints = sum(self.intervals(), ())
        return (min(endpoints), max(endpoints))

    def functions(self):
        """
        Returns the list of functions (the "pieces").

        EXAMPLES:
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f.functions()
            [x |--> 1, x |--> 1 - x, x |--> e^x, x |--> sin(2*x)]
        """
        return self._functions

    def extend_by_zero_to(self,xmin=-1000,xmax=1000):
        """
        This function simply returns the piecewise defined
        function which is extended by 0 so it is defined on all of (xmin,xmax).
        This is needed to add two piecewise functions in a reasonable way.

        EXAMPLES:

        """
        fcns = self.functions()
        f0 = fcns[0]
        R1 = f0.parent()
        xx = R1.gen()
        endpts = self.end_points()
        a = min(endpts)
        b = max(endpts)
        if xmin<a and xmax>b:
            F = Piecewise([[(xmin,a),0*xx]]+[[p[0],p[1]] for p in self.list()]+[[(b,xmax),0*xx]])
            return F
        if xmin>=a and xmax>b:
            F = Piecewise([[p[0],p[1]] for p in self.list()]+[[(b,xmax),0*xx]])
            return F
        if xmin<a and xmax<=b:
            F = Piecewise([[(xmin,a),0*xx]]+[[p[0],p[1]] for p in self.list()])
            return F
        return self

    def unextend(self):
        """
        This removes any parts in the front or back of the function which is
        zero (the inverse to extend_by_zero_to).

        EXAMPLES:
            sage: x = PolynomialRing(QQ,'x').gen()
            sage: f = Piecewise([[(-3,-1),1+2+x],[(-1,1),1-x^2]])
            sage: e = f.extend_by_zero_to(-10,10); e
            Piecewise defined function with 4 parts, [[(-10, -3), 0], [(-3, -1), x + 3], [(-1, 1), -x^2 + 1], [(1, 10), 0]]
            sage: d = e.unextend(); d
            Piecewise defined function with 2 parts, [[(-3, -1), x + 3], [(-1, 1), -x^2 + 1]]
            sage: d==f
            True

        """
        fcns = self.functions()
        intvls = self.intervals()
        n = len(fcns)
        if fcns[0]==0 and fcns[n-1]!=0:
            L = self.list()
            F = Piecewise([[p[0],p[1]] for p in L[1:]])
            return F
        if fcns[0]!=0 and fcns[n-1]==0:
            L = self.list()
            F = Piecewise([[p[0],p[1]] for p in L[:-1]])
            return F
        if fcns[0]==0 and fcns[n-1]==0:
            L = self.list()
            F = Piecewise([[p[0],p[1]] for p in L[1:-1]])
            return F
        return self

    def riemann_sum_integral_approximation(self,N,mode=None):
        """
        Returns the piecewise line function defined by the
        Riemann sums in numerical integration based on a subdivision
        into N subintervals.
    	Set mode="midpoint" for the height of the rectangles to be
        determined by the midpoint of the subinterval;
        set mode="right" for the height of the rectangles to be
        determined by the right-hand endpoint of the subinterval;
        the default is mode="left" (the height of the rectangles to be
        determined by the leftt-hand endpoint of the subinterval).

        EXAMPLES:
            sage: f1 = x^2                   ## example 1
            sage: f2 = 5-x^2
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.riemann_sum_integral_approximation(6)
            17/6
            sage: f.riemann_sum_integral_approximation(6,mode="right")
            19/6
            sage: f.riemann_sum_integral_approximation(6,mode="midpoint")
            3
            sage: f.integral()
            3
        """
        x = PolynomialRing(QQ,'x').gen()
        b = max(self.end_points())
        a = min(self.end_points())
        rsum = 0
        if mode==None:
	    for i in range(N):
                x0 = a+i*(b-a)/N
                x1 = a+(i+1)*(b-a)/N
                f0 = self(x0)
                rsum = rsum + (x1-x0)*f0
            return rsum
	if mode=="right":
	    for i in range(N):
                x0 = a+i*(b-a)/N
                x1 = a+(i+1)*(b-a)/N
                f0 = self(x1)
                rsum = rsum + (x1-x0)*f0
            return rsum
	if mode=="midpoint":
	    for i in range(N):
                x0 = a+i*(b-a)/N
                x1 = a+(i+1)*(b-a)/N
                f0 = self((x0+x1)/2)
                rsum = rsum + (x1-x0)*f0
            return rsum

    def riemann_sum(self,N,mode=None):
        """
        Returns the piecewise line function defined by the
        Riemann sums in numerical integration based on a subdivision
        into N subintervals.
	Set mode="midpoint" for the height of the rectangles to be
	determined by the midpoint of the subinterval;
        set mode="right" for the height of the rectangles to be
	determined by the right-hand endpoint of the subinterval;
	the default is mode="left" (the height of the rectangles to be
	determined by the leftt-hand endpoint of the subinterval).

        EXAMPLES:
            sage: f1(x) = x^2
            sage: f2(x) = 5-x^2
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.riemann_sum(6,mode="midpoint")
            Piecewise defined function with 6 parts, [[(0, 1/3), 1/36], [(1/3, 2/3), 1/4], [(2/3, 1), 25/36], [(1, 4/3), 131/36], [(4/3, 5/3), 11/4], [(5/3, 2), 59/36]]
            sage: f = Piecewise([[(-1,1),1-x^2]])
            sage: rsf = f.riemann_sum(7)
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: Q = rsf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: L = add([line([[pf[0][0],0],[pf[0][0],pf[1](x=pf[0][0])]],rgbcolor=(0.7,0.6,0.6)) for pf in rsf.list()])
            sage: ## To view this, type show(P+Q+L).
            sage: f = Piecewise([[(-1,1),1/2+x-x^3]]) ## example 3
            sage: rsf = f.riemann_sum(8)
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: Q = rsf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: L = add([line([[pf[0][0],0],[pf[0][0],pf[1](x=pf[0][0])]],rgbcolor=(0.7,0.6,0.6)) for pf in rsf.list()])
            sage: ## To view this, type show(P+Q+L).
        """
        x = PolynomialRing(QQ,'x').gen()
        b = max(self.end_points())
        a = min(self.end_points())
        rlist=[]
	if mode==None:
            for i in range(N):
                x0 = a+i*(b-a)/N
                x1 = a+(i+1)*(b-a)/N
                f0 = self(x0)
                rlist.append([(x0,x1),f0*x**0])
            return Piecewise(rlist)
	if mode=="right":
            for i in range(N):
                x0 = a+i*(b-a)/N
                x1 = a+(i+1)*(b-a)/N
                f0 = self(x1)
                rlist.append([(x0,x1),f0*x**0])
            return Piecewise(rlist)
	if mode=="midpoint":
            for i in range(N):
                x0 = a+i*(b-a)/N
                x1 = a+(i+1)*(b-a)/N
                f0 = self((x0+x1)/2)
                rlist.append([(x0,x1),f0*x**0])
            return Piecewise(rlist)


    def trapezoid(self,N):
        """
        Returns the piecewise line function defined by the
        trapezoid rule for numerical integration based on a subdivision
        into N subintervals.

        EXAMPLES:
            sage: f1 = lambda x:x^2                   ## example 1
            sage: f2 = lambda x:5-x^2
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.trapezoid(4)
            Piecewise defined function with 4 parts, [[(0, 1/2), 1/2*x], [(1/2, 1), 9/2*x - 2], [(1, 3/2), 1/2*x + 2], [(3/2, 2), -7/2*x + 8]]
            sage: x = PolynomialRing(QQ,'x').gen()        ## example 2
            sage: f = Piecewise([[(-1,1),1-x^2]])
            sage: tf = f.trapezoid(4)
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: Q = tf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: L = add([line([[pf[0][0],0],[pf[0][0],pf[1](pf[0][0])]],rgbcolor=(0.7,0.6,0.6)) for pf in tf.list()])
            sage: ## To view this, type show(P+Q+L).
            sage: f = Piecewise([[(-1,1),1/2+x-x^3]]) ## example 3
            sage: tf = f.trapezoid(6)
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: Q = tf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: L = add([line([[pf[0][0],0],[pf[0][0],pf[1](pf[0][0])]],rgbcolor=(0.7,0.6,0.6)) for pf in tf.list()])
            sage: ## To view this, type show(P+Q+L).
        """
        x = PolynomialRing(QQ,'x').gen()
        b = max(self.end_points())
        a = min(self.end_points())
        traplist=[]
        for i in range(N):
            x0 = a+i*(b-a)/N
            x1 = a+(i+1)*(b-a)/N
            f0 = self(x0)
            f1 = self(x1)
            traplist.append([(x0,x1),f0+(f1-f0)*(x1-x0)**(-1)*(x-x0)])
        return Piecewise(traplist)

    def trapezoid_integral_approximation(self,N):
        """
        Returns the approximation given by the
        trapezoid rule for numerical integration based on a subdivision
        into N subintervals.

        EXAMPLES:
            sage: f1(x) = x^2                      ## example 1
            sage: f2(x) = 1-(1-x)^2
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: tf = f.trapezoid(6)
            sage: Q = tf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: ta = f.trapezoid_integral_approximation(6)
            sage: t = text('trapezoid approximation = %s'%ta, (1.5, 0.25))
            sage: a = f.integral()
            sage: tt = text('area under curve = %s'%a, (1.5, -0.5))
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])  ## example 2
            sage: tf = f.trapezoid(4)
            sage: ta = f.trapezoid_integral_approximation(4)
            sage: Q = tf.plot(rgbcolor=(0.7,0.6,0.6), plot_points=40)
            sage: t = text('trapezoid approximation = %s'%ta, (1.5, 0.25))
            sage: a = f.integral()
            sage: tt = text('area under curve = %s'%a, (1.5, -0.5))
            sage: ## To view this, type show(P+Q+L).
        """
        x = PolynomialRing(QQ, 'x').gen()
        b = max(self.end_points())
        a = min(self.end_points())
        trapapprox = 0
        for i in range(N):
            x0 = a+i*(b-a)/N
            x1 = a+(i+1)*(b-a)/N
            f0 = self(x0)
            f1 = self(x1)
            trapapprox = trapapprox + ((f1+f0)/2)*(x1-x0)
        return trapapprox

    def critical_points(self):
        """
        Return the critical points of this piecewise function.

        WARNINGS: Uses maxima, which prints the warning to use results
        with caution. Only works for piecewise functions whose parts
        are polynomials with real critical not occurring on the
        interval endpoints.

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: f1 = x^0
            sage: f2 = 10*x - x^2
            sage: f3 = 3*x^4 - 156*x^3 + 3036*x^2 - 26208*x
            sage: f = Piecewise([[(0,3),f1],[(3,10),f2],[(10,20),f3]])
            sage: expected = [5, 12, 13, 14]
            sage: all(abs(e-a) < 0.001 for e,a in zip(expected, f.critical_points()))
            True
        """
        maxima = sage.interfaces.all.maxima
        x = PolynomialRing(QQ,'x').gen()
        fcns = self.functions()
        N = len(fcns)
        crit_pts = []
        I = self.intervals()
        for i in range(N):
            maxima.eval("eqn:diff(%s,x)=0"%fcns[i])
            ans = maxima.eval("allroots(eqn)")
            while True:
                start = ans.find('x=')
                if start == -1:
                    break
                ans = ans[start+2:]
                end = ans.find(',')
                if end == -1:
                    end = ans.find(']')
                r = float(ans[:end])
                if I[i][0] < r < I[i][1]:
                    crit_pts.append(r)
                ans = ans[end+1:]
        return crit_pts

    def base_ring(self):
        """
        Returns the base-ring (ie, QQ[x]) - useful when this
        class is extended.
        """
        return (self.functions()[0]).base_ring()

    def end_points(self):
        n = self.length()
        endpts = [self.intervals()[0][0]]
        for i in range(n):
            endpts.append(self.intervals()[i][1])
        return endpts

    def __call__(self,x0):
        """
        Evaluates self at x0. Returns the average value of the jump if x0 is
        an interior endpoint of one of the intervals of self and the
        usual value otherwise.

        EXAMPLES:
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f(0.5)
            1
            sage: f(2.5)
            12.1824939607
            sage: f(1)
            1/2
        """
        #x0 = QQ(x0) ## does not allow for evaluation at pi
        n = self.length()
        endpts = self.end_points()
        for i in range(1,n):
            if x0 == endpts[i]:
                return (self.functions()[i-1](x0) + self.functions()[i](x0))/2
        if x0 == endpts[0]:
            return self.functions()[0](x0)
        if x0 == endpts[n]:
            return self.functions()[n-1](x0)
        for i in range(n):
            if endpts[i] < x0 < endpts[i+1]:
                return self.functions()[i](x0)
        raise ValueError,"Value not defined outside of domain."

    def which_function(self,x0):
        """
        Returns the function piece used to evaluate self at x0.

        EXAMPLES:
            sage: f1(z) = z
            sage: f2(x) = 1-x
            sage: f3(y) = exp(y)
            sage: f4(t) = sin(2*t)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f.which_function(3/2)
            x |--> 1 - x
        """
        n = self.length()
        endpts = self.end_points()
        for i in range(1,n):
            if x0 == endpts[i]:
                return self.functions()[i-1]
        if x0 == endpts[0]:
            return self.functions()[0]
        if x0 == endpts[n]:
            return self.functions()[n-1]
        for i in range(n):
            if endpts[i] < x0 < endpts[i+1]:
                return self.functions()[i]
        raise ValueError,"Function not defined outside of domain."

    def integral(self, x=None):
        r"""
        Returns the definite integral (as computed by maxima)
        $\sum_I \int_I self|_I$, as I runs over the intervals
        belonging to self.

        EXAMPLES:
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.integral()
            1/2
            sage: f1(x) = -1
            sage: f2(x) = 2
            sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.integral()
            pi/2
        """
        #maxima = sage.interfaces.all.maxima
        #x = PolynomialRing(QQ,'x').gen()
        #ints = [maxima('%s'%p[1](x)).integral('x', p[0][0], p[0][1]) \
        #         for p in self.list()]
        #RETURN MEVAl(repr(sum(ints)))
        funcs = self.functions()
        invs = self.intervals()
        n = len(funcs)
        return sum([funcs[i].integral(x,invs[i][0],invs[i][1]) for i in range(n)])

    def convolution(self,other):
        """
        Returns the convolution function, $f*g(t)=\int_{-\infty}^\infty f(u)g(t-u)du$,
        for compactly supported $f,g$.

        EXAMPLES:
            sage: x = PolynomialRing(QQ,'x').gen()
            sage: f = Piecewise([[(0,1),1*x^0]])  ## example 0
            sage: g = f.convolution(f)
            sage: h = f.convolution(g)
            sage: P = f.plot(); Q = g.plot(rgbcolor=(1,1,0)); R = h.plot(rgbcolor=(0,1,1));
            sage: # Type show(P+Q+R) to view
            sage: f = Piecewise([[(0,1),1*x^0],[(1,2),2*x^0],[(2,3),1*x^0]])  ## example 1
            sage: g = f.convolution(f)
            sage: h = f.convolution(g)
            sage: P = f.plot(); Q = g.plot(rgbcolor=(1,1,0)); R = h.plot(rgbcolor=(0,1,1));
            sage: # Type show(P+Q+R) to view
            sage: f = Piecewise([[(-1,1),1]])                             ## example 2
            sage: g = Piecewise([[(0,3),x]])
            sage: f.convolution(g)
            Piecewise defined function with 3 parts, [[(-1, 1), 0], [(1, 2), -3/2*x], [(2, 4), -3/2*x]]
            sage: g = Piecewise([[(0,3),1*x^0],[(3,4),2*x^0]])
            sage: f.convolution(g)
            Piecewise defined function with 5 parts, [[(-1, 1), x + 1], [(1, 2), 3], [(2, 3), x], [(3, 4), -x + 8], [(4, 5), -2*x + 10]]

        """
        maxima = sage.interfaces.all.maxima
        f = self
        g = other
        M = min(min(f.end_points()),min(g.end_points()))
        N = max(max(f.end_points()),max(g.end_points()))
        R2 = PolynomialRing(QQ,2,names=["tt","uu"])
        tt,uu = R2.gens()
        conv = 0
        f0 = f.functions()[0]
        g0 = g.functions()[0]
        R1 = f0.parent()
        xx = R1.gen()
        var = repr(xx)
        if len(f.intervals())==1 and len(g.intervals())==1:
            f = f.unextend()
            g = g.unextend()
            a1 = f.intervals()[0][0]
            a2 = f.intervals()[0][1]
            b1 = g.intervals()[0][0]
            b2 = g.intervals()[0][1]
            i1 = repr(f0).replace(var,repr(uu))
            i2 = repr(g0).replace(var,"("+repr(tt-uu)+")")
            cmd1 = "integrate((%s)*(%s),%s,%s,%s)"%(i1,i2, uu, a1,    tt-b1)    ## if a1+b1 < tt < a2+b1
            cmd2 = "integrate((%s)*(%s),%s,%s,%s)"%(i1,i2, uu, tt-b2, tt-b1)    ## if a1+b2 < tt < a2+b1
            cmd3 = "integrate((%s)*(%s),%s,%s,%s)"%(i1,i2, uu, tt-b2, a2)       ## if a1+b2 < tt < a2+b2
            cmd4 = "integrate((%s)*(%s),%s,%s,%s)"%(i1,i2, uu, a1, a2)          ## if a2+b1 < tt < a1+b2
            conv1 = maxima.eval(cmd1)
            conv2 = maxima.eval(cmd2)
            conv3 = maxima.eval(cmd3)
            conv4 = maxima.eval(cmd4)
            # this is a very, very, very ugly hack
            x = PolynomialRing(QQ,'x').gen()
            fg1 = sage_eval(conv1.replace("tt",var), {'x':x}) ## should be = R2(conv1)
            fg2 = sage_eval(conv2.replace("tt",var), {'x':x}) ## should be = R2(conv2)
            fg3 = sage_eval(conv3.replace("tt",var), {'x':x}) ## should be = R2(conv3)
            fg4 = sage_eval(conv4.replace("tt",var), {'x':x}) ## should be = R2(conv4)
            if a1-b1<a2-b2:
                if a2+b1!=a1+b2:
                    h = Piecewise([[(a1+b1,a1+b2),fg1],[(a1+b2,a2+b1),fg4],[(a2+b1,a2+b2),fg3]])
                else:
                    h = Piecewise([[(a1+b1,a1+b2),fg1],[(a1+b2,a2+b2),fg3]])
            else:
                if a1+b2!=a2+b1:
                    h = Piecewise([[(a1+b1,a2+b1),fg1],[(a2+b1,a1+b2),fg2],[(a1+b2,a2+b2),fg3]])
                else:
                    h = Piecewise([[(a1+b1,a2+b1),fg1],[(a2+b1,a2+b2),fg3]])
            return h

        if len(f.intervals())>1 or len(g.intervals())>1:
            z = Piecewise([[(-3*abs(N-M),3*abs(N-M)),0*xx**0]])
            ff = f.functions()
            gg = g.functions()
            intvlsf = f.intervals()
            intvlsg = g.intervals()
            for i in range(len(ff)):
                for j in range(len(gg)):
                    f0 = Piecewise([[intvlsf[i],ff[i]]])
                    g0 = Piecewise([[intvlsg[j],gg[j]]])
                    h = g0.convolution(f0)
                    z = z + h
            return z.unextend()

    def derivative(self):
        r"""
        Returns the derivative (as computed by maxima)
        Piecewise(I,$(d/dx)(self|_I)$), as I runs over the intervals
        belonging to self. self must be piecewise polynomial.

        EXAMPLES:
            sage: f1 = lambda x:1
            sage: f2 = lambda x:1-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.derivative()
            Piecewise defined function with 2 parts, [[(0, 1), 0], [(1, 2), -1]]
	    sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.derivative()
            Piecewise defined function with 2 parts, [[(0, pi/2), 0], [(pi/2, pi), 0]]

            sage: f = Piecewise([[(0,1), x * 2]])
            sage: f.derivative()
            Piecewise defined function with 1 parts, [[(0, 1), 2]]
        """
        x = var('x')
        dlist = [[(p[0][0], p[0][1]), derivative(p[1](x), x)] for p in self.list()]
        return Piecewise(dlist)

    def tangent_line(self,pt):
        """
        Computes the linear function defining the tangent line of
        the piecewise function self.

        EXAMPLES:
            sage: f1 = lambda x:x^2
            sage: f2 = lambda x:5-x^3+x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: tf = f.tangent_line(0.9) ## tangent line at x=0.9
            sage: P = f.plot(rgbcolor=(0.7,0.1,0.5), plot_points=40)
            sage: Q = tf.plot(rgbcolor=(0.7,0.2,0.2), plot_points=40)

        Type show(P+Q) to view the graph of the function and the tangent line.
        """
        pt = QQ(pt)
        R = PolynomialRing(QQ,'x')
        x = R.gen()
        der = self.derivative()
        tanline = (x-pt)*der(pt)+self(pt)
        dlist = [[(p[0][0], p[0][1]), tanline] for p in self.list()]
        return Piecewise(dlist)

    def plot(self, *args, **kwds):
        """
        Returns the plot of self.

        Keyword arguments are passed onto the plot command for each
        piece of the function.  E.g., the plot_points keyword affects
        each segment of the plot.

        EXAMPLES:
            sage: f1 = lambda x:1
	    sage: f2 = lambda x:1-x
            sage: f3 = lambda x:exp(x)
            sage: f4 = lambda x:sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: P = f.plot(rgbcolor=(0.7,0.1,0), plot_points=40)

	Remember: to view this, type show(P) or P.save("<path>/myplot.png") and
        then open it in a graphics viewer such as GIMP.
        """
        plot = sage.plot.plot.plot
        return sum([plot(p[1], p[0][0], p[0][1], *args, **kwds ) for p in self.list()])

    def fourier_series_cosine_coefficient(self,n,L):
        r"""
        Returns the n-th Fourier series coefficient of $\cos(n\pi x/L)$, $a_n$.

        INPUT:
            self -- the function f(x), defined over -L < x < L
            n    -- an integer n>=0
            L    -- (the period)/2

        OUTPUT:
            $a_n = \frac{1}{L}\int_{-L}^L f(x)\cos(n\pi x/L)dx$

        EXAMPLES:
            sage: f = lambda x:x^2
            sage: f = Piecewise([[(-1,1),f]])
            sage: f.fourier_series_cosine_coefficient(2,1)
            1/pi^2
	    sage: f = lambda x:x^2
            sage: f = Piecewise([[(-pi,pi),f]])
            sage: float(f.fourier_series_cosine_coefficient(2,pi))
            1.0
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_cosine_coefficient(5,pi)
            -3/(5*pi)
        """
        maxima = sage.interfaces.all.maxima
        x = PolynomialRing(QQ,'x').gen()
        ints = []
        for p in self.list():
            fcn = '(%s)*cos('%p[1](x) + 'pi*x*%s/%s)/%s'%(n,L,L)
            fcn = fcn.replace("pi","%"+"pi")
	    a = repr(p[0][0]).replace("pi","%"+"pi")
	    b = repr(p[0][1]).replace("pi","%"+"pi")
	    cmd = "integrate("+fcn+", x, %s, %s )"%(a, b)
	    int = maxima(cmd).trigsimp()
            ints.append(int)
        ans = sum(ints)
        return meval(repr(ans))

    def fourier_series_sine_coefficient(self,n,L):
        r"""
        Returns the n-th Fourier series coefficient of $\sin(n\pi x/L)$, $b_n$.

        INPUT:
            self -- the function f(x), defined over -L < x < L
            n    -- an integer n>0
            L    -- (the period)/2

        OUTPUT:
            $b_n = \frac{1}{L}\int_{-L}^L f(x)\sin(n\pi x/L)dx$

        EXAMPLES:
            sage: f = lambda x:x^2
            sage: f = Piecewise([[(-1,1),f]])
            sage: f.fourier_series_sine_coefficient(2,1)  # L=1, n=2
            0
        """
	maxima = sage.interfaces.all.maxima
        x = PolynomialRing(QQ,'x').gen()
        ints = []
        for p in self.list():
            fcn = '(%s)*sin('%p[1](x) + 'pi*x*%s/%s)/%s'%(n,L,L)
            fcn = fcn.replace("pi","%"+"pi")
	    a = repr(p[0][0]).replace("pi","%"+"pi")
	    b = repr(p[0][1]).replace("pi","%"+"pi")
	    cmd = "integrate("+fcn+", x, %s, %s )"%(a, b)
	    int = maxima(cmd).trigsimp()
            ints.append(int)
        ans = sum(ints)
        return meval(repr(ans))

    def fourier_series_partial_sum(self,N,L):
        r"""
        Returns the partial sum
        \[
        f(x) \sim \frac{a_0}{2} +
                   \sum_{n=1}^N [a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],
        \]
        as a string.

        EXAMPLE:
            sage: f = lambda x:x^2
            sage: f = Piecewise([[(-1,1),f]])
            sage: f.fourier_series_partial_sum(3,1)
            cos(2*pi*x)/pi^2 - 4*cos(pi*x)/pi^2 + 1/3
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_partial_sum(3,pi)
            -3*sin(2*x)/pi + 3*sin(x)/pi - 3*cos(x)/pi - 1/4
        """
        a0 = self.fourier_series_cosine_coefficient(0,L)
        A = [repr(self.fourier_series_cosine_coefficient(n,L))+"*cos(%s*pi*x/%s)"%(n,L) for n in range(1,N)]
        B = [repr(self.fourier_series_sine_coefficient(n,L))+"*sin(%s*pi*x/%s)"%(n,L) for n in range(1,N)]
        FS =  ["("+A[i] +" + " + B[i]+")" for i in range(0,N-1)]
        sumFS = repr(a0/2)+" + "
        for s in FS:
            sumFS = sumFS+s+ " + "
        return meval(sumFS[:-3])

    def fourier_series_partial_sum_cesaro(self,N,L):
        r"""
        Returns the Cesaro partial sum
        \[
        f(x) \sim \frac{a_0}{2} +
                   \sum_{n=1}^N (1-n/N)*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],
        \]
        as a string. This is a "smoother" partial sum - the Gibbs phenomenon is mollified.

        EXAMPLE:
            sage: f = lambda x:x^2
            sage: f = Piecewise([[(-1,1),f]])
            sage: f.fourier_series_partial_sum_cesaro(3,1)
            cos(2*pi*x)/(3*pi^2) - 8*cos(pi*x)/(3*pi^2) + 1/3
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_partial_sum_cesaro(3,pi)
            -sin(2*x)/pi + 2*sin(x)/pi - 2*cos(x)/pi - 1/4

        """
        a0 = self.fourier_series_cosine_coefficient(0,L)
        A = [repr((1-n/N)*self.fourier_series_cosine_coefficient(n,L))+"*cos(%s*pi*x/%s)"%(n,L) for n in range(1,N)]
        B = [repr((1-n/N)*self.fourier_series_sine_coefficient(n,L))+"*sin(%s*pi*x/%s)"%(n,L) for n in range(1,N)]
        FS =  ["("+A[i] +" + " + B[i]+")" for i in range(0,N-1)]
        sumFS = repr(a0/2)+" + "
        for s in FS:
            sumFS = sumFS+s+ " + "
        return meval(sumFS[:-3])

    def fourier_series_partial_sum_hann(self,N,L):
        r"""
        Returns the Hann-filtered partial sum (named after von Hann, not Hamming)
        \[
        f(x) \sim \frac{a_0}{2} +
                   \sum_{n=1}^N H_N(n)*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],
        \]
        as a string, where $H_N(x) = (1+\cos(\pi x/N))/2$. This is a "smoother" partial sum - the Gibbs phenomenon is mollified.

        EXAMPLE:
            sage: f = lambda x:x^2
            sage: f = Piecewise([[(-1,1),f]])
            sage: f.fourier_series_partial_sum_hann(3,1)
            0.5*(cos(2*pi/3) + 1)*cos(2*pi*x)/pi^2 - 2.0*(cos(pi/3) + 1)*cos(pi*x)/pi^2 + 1/3
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_partial_sum_hann(3,pi)
            -1.5*(cos(2*pi/3) + 1)*sin(2*x)/pi + 1.5*(cos(pi/3) + 1)*sin(x)/pi - 1.5*(cos(pi/3) + 1)*cos(x)/pi - 1/4

        """
        a0 = self.fourier_series_cosine_coefficient(0,L)
        A = ["(1+cos(pi*%s/%s))*"%(n,N)+repr((0.5)*self.fourier_series_cosine_coefficient(n,L))+"*cos(%s*pi*x/%s)"%(n,L) for n in range(1,N)]
        B = ["(1+cos(pi*%s/%s))*"%(n,N)+repr((0.5)*self.fourier_series_sine_coefficient(n,L))+"*sin(%s*pi*x/%s)"%(n,L) for n in range(1,N)]
        FS =  ["("+A[i] +" + " + B[i]+")" for i in range(0,N-1)]
        sumFS = repr(a0/2)+" + "
        for s in FS:
            sumFS = sumFS+s+ " + "
        return meval(sumFS[:-3])

    def fourier_series_partial_sum_filtered(self,N,L,F):
        r"""
        Returns the "filtered" partial sum
        \[
        f(x) \sim \frac{a_0}{2} +
                   \sum_{n=1}^N F_n*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],
        \]
        as a string, where $F = [F_1,F_2, ..., F_{N}]$ is a list of length $N$ consisting of real numbers.
	This can be used to plot FS solutions to the heat and wave PDEs.

        EXAMPLE:
            sage: f = lambda x:x^2
            sage: f = Piecewise([[(-1,1),f]])
            sage: f.fourier_series_partial_sum_filtered(3,1,[1,1,1])
            cos(2*pi*x)/pi^2 - 4*cos(pi*x)/pi^2 + 1/3
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_partial_sum_filtered(3,pi,[1,1,1])
            -3*sin(2*x)/pi + 3*sin(x)/pi - 3*cos(x)/pi - 1/4
        """
        a0 = self.fourier_series_cosine_coefficient(0,L)
        A = [repr((F[n])*self.fourier_series_cosine_coefficient(n,L))+"*cos(%s*pi*x/%s)"%(n,L) for n in range(1,N)]
        B = [repr((F[n])*self.fourier_series_sine_coefficient(n,L))+"*sin(%s*pi*x/%s)"%(n,L) for n in range(1,N)]
        FS =  ["("+A[i] +" + " + B[i]+")" for i in range(0,N-1)]
        sumFS = repr(a0/2)+" + "
        for s in FS:
            sumFS = sumFS+s+ " + "
        return meval(sumFS[:-3])

    def plot_fourier_series_partial_sum(self,N,L,xmin,xmax, **kwds):
        r"""
        Plots the partial sum
        \[
        f(x) \sim \frac{a_0}{2} +
                   \sum_{n=1}^N [a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],
        \]
        over xmin < x < xmin.

        EXAMPLE:
            sage: f1 = lambda x:-2
            sage: f2 = lambda x:1
            sage: f3 = lambda x:-1
            sage: f4 = lambda x:2
            sage: f = Piecewise([[(-pi,-pi/2),f1],[(-pi/2,0),f2],[(0,pi/2),f3],[(pi/2,pi),f4]])
            sage: P = f.plot_fourier_series_partial_sum(3,pi,-5,5)    # long time
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: P = f.plot_fourier_series_partial_sum(15,pi,-5,5)   # long time

	Remember, to view this type show(P) or P.save("<path>/myplot.png") and then
        open it in a graphics viewer such as GIMP.
        """
        line = sage.plot.plot.line
        pts = []
        h = QQ(1)/QQ(10)
        n = int((xmax - xmin)/h) + 1
        Pi = 3.14159265
        ff = self.fourier_series_partial_sum(N,L)
        for i in range(n):
            pi = 3.14159265
            xi = xmin + i*h
	    yi = ff(xi)
            pts.append([xi,yi])
        return line(pts, **kwds)

    def plot_fourier_series_partial_sum_cesaro(self,N,L,xmin,xmax, **kwds):
        r"""
        Plots the partial sum
        \[
        f(x) \sim \frac{a_0}{2} +
                   \sum_{n=1}^N (1-n/N)*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],
        \]
        over xmin < x < xmin. This is a "smoother" partial sum - the Gibbs phenomenon is mollified.

        EXAMPLE:
            sage: f1 = lambda x:-2
            sage: f2 = lambda x:1
            sage: f3 = lambda x:-1
            sage: f4 = lambda x:2
            sage: f = Piecewise([[(-pi,-pi/2),f1],[(-pi/2,0),f2],[(0,pi/2),f3],[(pi/2,pi),f4]])
            sage: P = f.plot_fourier_series_partial_sum_cesaro(3,pi,-5,5)    # long time
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: P = f.plot_fourier_series_partial_sum_cesaro(15,pi,-5,5)   # long time

	Remember, to view this type show(P) or P.save("<path>/myplot.png") and then
        open it in a graphics viewer such as GIMP.
        """
        line = sage.plot.plot.line
        pts = []
        h = QQ(1)/QQ(10)
        n = int((xmax - xmin)/h) + 1
        Pi = 3.14159265
        ff = self.fourier_series_partial_sum_cesaro(N,L)
        for i in range(n):
            pi = 3.14159265
            xi = xmin + i*h
	    yi = ff(xi)
            pts.append([xi,yi])
        return line(pts, **kwds)

    def plot_fourier_series_partial_sum_hann(self,N,L,xmin,xmax, **kwds):
        r"""
        Plots the partial sum
        \[
        f(x) \sim \frac{a_0}{2} +
                   \sum_{n=1}^N H_N(n)*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],
        \]
        over xmin < x < xmin, where H_N(x) = (0.5)+(0.5)*cos(x*pi/N) is the N-th Hann filter.

        EXAMPLE:
            sage: f1 = lambda x:-2
            sage: f2 = lambda x:1
            sage: f3 = lambda x:-1
            sage: f4 = lambda x:2
            sage: f = Piecewise([[(-pi,-pi/2),f1],[(-pi/2,0),f2],[(0,pi/2),f3],[(pi/2,pi),f4]])
            sage: P = f.plot_fourier_series_partial_sum_hann(3,pi,-5,5)    # long time
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(-pi,pi/2),f1],[(pi/2,pi),f2]])
            sage: P = f.plot_fourier_series_partial_sum_hann(15,pi,-5,5)   # long time

	Remember, to view this type show(P) or P.save("<path>/myplot.png") and then
        open it in a graphics viewer such as GIMP.
        """
        line = sage.plot.plot.line
        pts = []
        h = QQ(1)/QQ(10)
        n = int((xmax - xmin)/h) + 1
        Pi = 3.14159265
        ff = self.fourier_series_partial_sum_hann(N,L)
        for i in range(n):
            pi = 3.14159265
            xi = xmin + i*h
	    yi = ff(xi)
            pts.append([xi,yi])
        return line(pts, **kwds)

    def plot_fourier_series_partial_sum_filtered(self,N,L,F,xmin,xmax, **kwds):
        r"""
        Plots the partial sum
        \[
        f(x) \sim \frac{a_0}{2} +
                   \sum_{n=1}^N F_n*[a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],
        \]
        over xmin < x < xmin, where $F = [F_1,F_2, ..., F_{N}]$ is a list of length
	$N$ consisting of real numbers. This can be used to plot FS solutions to the heat
	and wave PDEs.

        EXAMPLE:
            sage: f1 = lambda x:-2
            sage: f2 = lambda x:1
            sage: f3 = lambda x:-1
            sage: f4 = lambda x:2
            sage: f = Piecewise([[(-pi,-pi/2),f1],[(-pi/2,0),f2],[(0,pi/2),f3],[(pi/2,pi),f4]])
            sage: P = f.plot_fourier_series_partial_sum_filtered(3,pi,[1]*3,-5,5)    # long time
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(-pi,-pi/2),f1],[(-pi/2,0),f2],[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: P = f.plot_fourier_series_partial_sum_filtered(15,pi,[1]*15,-5,5)   # long time

	Remember, to view this type show(P) or P.save("<path>/myplot.png") and then
        open it in a graphics viewer such as GIMP.
        """
        line = sage.plot.plot.line
        ff = self.fourier_series_partial_sum_filtered(N,L,F)
	pts = []
        h = QQ(1)/QQ(10)
        n = int((xmax - xmin)/h) + 1
        Pi = 3.14159265
        for i in range(n):
            pi = 3.14159265
            xi = xmin + i*h
	    yi = ff(xi)
            pts.append([xi,yi])
        return line(pts, **kwds)

    def fourier_series_value(self,x,L):
        r"""
        Returns the value of the Fourier series coefficient of self at $x$,

        \[
        f(x) \sim \frac{a_0}{2} +
                   \sum_{n=1}^\infty [a_n\cos(\frac{n\pi x}{L}) + b_n\sin(\frac{n\pi x}{L})],
        \ \ \ -L<x<L.
        \]
        This method applies to piecewise non-polynomial functions as well.

        INPUT:
            self -- the function f(x), defined over -L < x < L
            x    -- a real number
            L    -- (the period)/2

        OUTPUT:
            $(f^*(x+)+f^*(x-)/2$, where $f^*$ denotes the function $f$
            extended to $\R$ with period $2L$ (Dirichlet's Theorem for
	    Fourier series).

        EXAMPLES:
            sage: f1 = lambda x:1
            sage: f2 = lambda x:1-x
            sage: f3 = lambda x:exp(x)
            sage: f4 = lambda x:sin(2*x)
            sage: f = Piecewise([[(-10,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f.fourier_series_value(101,10)
            1/2
            sage: f.fourier_series_value(100,10)
            1
            sage: f.fourier_series_value(10,10)
            sin(20)
            sage: f.fourier_series_value(20,10)
            1
            sage: f.fourier_series_value(30,10)
            sin(20)
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
	    sage: f = Piecewise([[(-pi,0),lambda x:0],[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_value(-1,pi)
            0
            sage: f.fourier_series_value(20,pi)
            -1
            sage: f.fourier_series_value(pi/2,pi)
            1/2
        """
        xnew = x - int(RR(x/(2*L)))*2*L
        endpts = self.end_points()
        n = self.length()
        if xnew == endpts[0] or xnew == endpts[n-1]:
            ave = (self.list()[0][1](endpts[0]) + self.list()[n-1][1](endpts[n-1]))/2
        return self(xnew)

    def cosine_series_coefficient(self,n,L):
        r"""
        Returns the n-th cosine series coefficient of $\cos(n\pi x/L)$, $a_n$.

        INPUT:
            self -- the function f(x), defined over 0 < x < L (no checking is done
	                                                       to insure this)
            n    -- an integer n>=0
            L    -- (the period)/2

        OUTPUT:
            $a_n = \frac{2}{L}\int_{-L}^L f(x)\cos(n\pi x/L)dx$ such that
        \[
        f(x) \sim \frac{a_0}{2} +
                   \sum_{n=1}^\infty a_n\cos(\frac{n\pi x}{L}),\ \ 0<x<L.
        \]

        EXAMPLES:
            sage: f = lambda x:x
            sage: f = Piecewise([[(0,1),f]])
            sage: f.cosine_series_coefficient(2,1)
            0
            sage: f.cosine_series_coefficient(3,1)
            -4/(9*pi^2)
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.cosine_series_coefficient(2,pi)
            0
            sage: f.cosine_series_coefficient(3,pi)
            2/pi
            sage: f.cosine_series_coefficient(111,pi)
            2/(37*pi)

        """
	maxima = sage.interfaces.all.maxima
        x = PolynomialRing(QQ,'x').gen()
        ints = []
        for p in self.list():
            fcn = '2*(%s)*cos('%p[1](x) + 'pi*x*%s/%s)/%s'%(n,L,L)
            fcn = fcn.replace("pi","%"+"pi")
	    a = repr(p[0][0]).replace("pi","%"+"pi")
	    b = repr(p[0][1]).replace("pi","%"+"pi")
	    cmd = "integrate("+fcn+", x, %s, %s )"%(a, b)
	    I = maxima(cmd).trigsimp()
            ints.append(I)
        ans = sum(ints)
        return meval(repr(ans))


    def sine_series_coefficient(self,n,L):
        r"""
        Returns the n-th sine series coefficient of $\sin(n\pi x/L)$, $b_n$.

        INPUT:
            self -- the function f(x), defined over 0 < x < L (no checking is done
	                                                       to insure this)
            n    -- an integer n>0
            L    -- (the period)/2

        OUTPUT:
            $b_n = \frac{2}{L}\int_{-L}^L f(x)\sin(n\pi x/L)dx$ such that
        \[
        f(x) \sim \sum_{n=1}^\infty b_n\sin(\frac{n\pi x}{L}),\ \ 0<x<L.
        \]

        EXAMPLES:
            sage: f = lambda x:1
            sage: f = Piecewise([[(0,1),f]])
            sage: f.sine_series_coefficient(2,1)
            0
            sage: f.sine_series_coefficient(3,1)
            4/(3*pi)
        """
	maxima = sage.interfaces.all.maxima
        x = PolynomialRing(QQ,'x').gen()
        ints = []
        for p in self.list():
            fcn = '2*(%s)*sin('%p[1](x) + 'pi*x*%s/%s)/%s'%(n,L,L)
            fcn = fcn.replace("pi","%"+"pi")
	    a = repr(p[0][0]).replace("pi","%"+"pi")
	    b = repr(p[0][1]).replace("pi","%"+"pi")
	    cmd = "integrate("+fcn+", x, %s, %s )"%(a, b)
	    I = maxima(cmd).trigsimp()
            ints.append(I)
        ans = sum(ints)
        return meval(repr(ans))

    def laplace(self, x, s):
        r"""
        Returns the Laplace transform of self with respect to the
        variable var.

        INPUT:
            x -- variable of self
            s -- variable of Laplace transform.

        We assume that a piecewise function is 0 outside of its domain
        and that the left-most endpoint of the domain is 0.

        EXAMPLES:
            sage: x, s, w = var('x, s, w')
            sage: f = Piecewise([[(0,1),1],[(1,2), 1-x]])
            sage: f.laplace(x, s)
            -e^(-s)/s - e^(-s)/s^2 + (s + 1)*e^(-(2*s))/s^2 + 1/s
            sage: f.laplace(x, w)
            -e^(-w)/w - e^(-w)/w^2 + (w + 1)*e^(-(2*w))/w^2 + 1/w

            sage: y, t = var('y, t')
            sage: f = Piecewise([[(1,2), 1-y]])
            sage: f.laplace(y, t)
            (t + 1)*e^(-(2*t))/t^2 - e^(-t)/t^2

            sage: s = var('s')
            sage: t = var('t')
            sage: f1(t) = -t
            sage: f2(t) = 2
            sage: f = Piecewise([[(0,1),f1],[(1,infinity),f2]])
            sage: f.laplace(t,s)
            (s + 1)*e^(-s)/s^2 + 2*e^(-s)/s - 1/s^2
        """
        from sage.calculus.equations import assume
        x = var(x)
        s = var(s)
        assume(s>0)
        ints = []
        for p in self.list():
            g = SR(p[1])
            fcn = maxima('(%s)*exp(-%s*%s)'%(g._maxima_init_(), s, x))
            ints.append(fcn.integral(x, p[0][0], p[0][1]))
        ans = ""
        for i in range(len(ints)-1):
            ans = ans+repr(ints[i]) + " + "
        ans = ans+repr(ints[len(ints)-1])
        return meval(ans)

    def __add__(self,other):
	"""
	Returns the piecewise defined function which is the sum of
	self and other. Does not require both domains be the same.

	EXAMPLES:
	    sage: x = PolynomialRing(QQ,'x').gen()
	    sage: f1 = x^0
            sage: f2 = 1-x
            sage: f3 = 2*x
            sage: f4 = 10-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
	    sage: g1 = x-2
            sage: g2 = x-5
            sage: g = Piecewise([[(0,5),g1],[(5,10),g2]])
	    sage: h = f+g
	    sage: h
            Piecewise defined function with 5 parts, [[(0, 1), x - 1], [(1, 2), -1], [(2, 3), 3*x - 2], [(3, 5), 8], [(5, 10), 5]]

        Note that in this case the functions must be defined using polynomial
	expressions *not* using the lambda notation.
	"""
	self_endpts = self.end_points()
        a1 = min(self_endpts); a2 = max(self_endpts)
	other_endpts = other.end_points()
        b1 = min(other_endpts); b2 = max(other_endpts)
        c1 = min(a1,b1); c2 = max(a2,b2)
        F = self.extend_by_zero_to(c1,c2)
        G = other.extend_by_zero_to(c1,c2)
	f = F.functions()
	g = G.functions()
	endpts = list(set(F.end_points()).union(set(G.end_points())))
        endpts.sort()
	N = len(list(endpts))
	fcn = []
	for j in range(N-1):
	    x0 = endpts[j+1]
	    fcn.append([(endpts[j],endpts[j+1]),F.which_function(x0)+G.which_function(x0)])
	return Piecewise(fcn)

    def __mul__(self,other):
	r"""
	Returns the piecewise defined function which is the product of
	one piecewise function (self) with another one (other).

	EXAMPLES:
	    sage: x = PolynomialRing(QQ,'x').gen()
	    sage: f1 = x^0
            sage: f2 = 1-x
            sage: f3 = 2*x
            sage: f4 = 10-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
	    sage: g1 = x-2
            sage: g2 = x-5
            sage: g = Piecewise([[(0,5),g1],[(5,10),g2]])
	    sage: h = f*g
	    sage: h
	    Piecewise defined function with 5 parts, [[(0, 1), x - 2], [(1, 2), -x^2 + 3*x - 2], [(2, 3), 2*x^2 - 4*x], [(3, 5), -x^2 + 12*x - 20], [(5, 10), -x^2 + 15*x - 50]]
            sage: g*(11/2)
            Piecewise defined function with 2 parts, [[(0, 5), 11/2*x - 11], [(5, 10), 11/2*x - 55/2]]

        Note that in this method the functions must be defined using polynomial
	expressions *not* using the lambda notation.
	"""
        R = PolynomialRing(QQ,'x')
        fcn = []
        if isinstance(other,Rational) or isinstance(other,Integer):    ## needed for scalar multiplication
            endpts = self.end_points()
            N = len(list(endpts))
            for j in range(N-1):
	        x0 = endpts[j+1]
	        fcn.append([(endpts[j],endpts[j+1]),R(other)*self.which_function(x0)])
            return Piecewise(fcn)
	self_endpts = self.end_points()     ## we assume these start
	other_endpts = other.end_points()   ## and end at the same point
	f = self.functions()
	g = other.functions()
	endpts = list(set(self_endpts).union(set(other_endpts)))
	N = len(list(endpts))
	for j in range(N-1):
	    x0 = endpts[j+1]
	    fcn.append([(endpts[j],endpts[j+1]),self.which_function(x0)*other.which_function(x0)])
	return Piecewise(fcn)

    __rmul__ = __mul__

    def __eq__(self,other):
        """
        Implements Boolean == operator.
        """
        return self.list()==other.list()
