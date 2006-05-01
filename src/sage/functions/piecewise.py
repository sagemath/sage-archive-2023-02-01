"""
    Implements a very simple class of piecewise-defined functions.
    Functions must be piecewise polynomial, though some methods
    apply more generally. Only compactly supported functions are
    currently implemented. Moreover, the coefficients should be
    rational and the support should be "connected". The intervals of
    polynomial support can be in terms of rationals and $\pi$, or in terms
    of floats.

    Implemented methods:
	latex outout
	__call__
	plotting
	fourier series
	fourier series
	   value
	   coefficients (also sine series and cosine series)
	   partial sum (in string format)
	   plot of partial sum
	laplace transform
	   latex output option
	domain
	range
	list
	addition (of functions)
	multiplication (of functions, or fcn*scalar - ie, *right* multiplication by QQ)
	critical points

    TODO:
        Implement (a) functions defined on infinite intervals,
        (b) max/min location and values,
	(c) left multiplication by a scalar.
        (d) Extend the implementation of the trick to pass SAGE's pi back
        and forth with MAXIMA's %pi to other constants (e, for example)
    (For more general non-polynomial piecewise functions, it appears
    a new class of functions (for example, "ElementaryFunctionRing") is
    needed. This a preliminary "todo".)

    AUTHOR: David Joyner (2006-04)

"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdj@usna.edu>
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
from sage.rings.polynomial_ring import PolynomialRing
from sage.rings.rational_field import RationalField
from sage.rings.real_field import RealField
from sage.misc.sage_eval import *
from sage.rings.rational import *
from sage.rings.integer import *
QQ = RationalField()
RR = RealField()

class PiecewisePolynomial:

    def __init__(self, list_of_pairs):
        r"""
        \code{list_of_pairs} is a list of pairs (fcn,I), where fcn is
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

    def latex(self):
	"""
	EXAMPLES:
            sage: f1 = lambda x:1
            sage: f2 = lambda x:1-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.latex()
            '\\begin{array}{ll} \\left\\{ 1,& 0 < x < 1 ,\\right. \\end{array}'

	"""
        x = PolynomialRing(RationalField()).gen()
        intvls = self.intervals()
        fcn_list = [p[1] for p in self.list()]
        tex = ["\\begin{array}{ll} \left\\{"]
        for i in range(len(fcn_list)):
	    f = fcn_list[i]
	    a = intvls[i][0]
	    b = intvls[i][1]
            tex.append(str(f(x)))
	    tex.append(",& %s < x < %s ,\\"%(a,b))
        tex = tex[:-2]
	tex.append("\right\. \end{array}")
	ltex = ""
        for i in range(len(tex)-1):
            ltex = ltex + tex[i]
        ltex = ltex + str(tex[len(tex)-1]).replace("%","")
        ltex = ltex.replace("\\\right\\.","\\right.")
        ltex = ltex.replace("\\left\\{","\\left\{ ")
        return ltex

    def intervals(self):
        """
	A piecewise non-polynomial example.

        EXAMPLES:
            sage: f1 = lambda x:1
            sage: f2 = lambda x:1-x
            sage: f3 = lambda x:exp(x)
            sage: f4 = lambda x:sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f.intervals()
            [(0, 1), (1, 2), (2, 3), (3, 10)]
        """
        return self._intervals

    def domain(self):
        return (min(self.intervals()),max(self.intervals()))

    def functions(self):
        """
        Returns the list of functions (the "pieces").
        """
        return self._functions

    def critical_points(self):
        """
        Function to return the critical points. Uses maxima, which prints the
        warning to use results with caution. Only works for piecewise functions
        whose parts are polynomials with real critical not occurring on the
        interval endpoints.

        EXAMPLES:
            sage: x = PolynomialRing(RationalField()).gen()
            sage: f1 = x^0
            sage: f2 = 1-x
            sage: f3 = 2*x
            sage: f4 = 10*x-x^2
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f.critical_points()
            [5.0]
        """
        maxima = sage.interfaces.all.maxima
        x = PolynomialRing(RationalField()).gen()
        fcns = self.functions()
        N = len(fcns)
        crit_pts = []
        for i in range(N):
            maxima.eval("eqn:diff(%s,x)=0"%fcns[i])
            ans = maxima.eval("allroots(eqn)")
            if "[x =" in ans:
                i1 = ans.index("[x =")
                i2 = ans.index("]")
                r = eval(ans[i1+4:i2])
                if self.intervals()[i][0] < r < self.intervals()[i][1]:
                    crit_pts.append(r)
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
            sage: f1 = lambda x:1
            sage: f2 = lambda x:1-x
            sage: f3 = lambda x:exp(x)
            sage: f4 = lambda x:sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f(0.5)
            1
            sage: f(2.5)
            12.182493960703473
            sage: f(1)
            1/2

        """
        n = self.length()
        endpts = self.end_points()
        for i in range(1,n):
            if x0 == endpts[i]:
                return (self.functions()[i-1](x0)+self.functions()[i](x0))/2
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
        Returns the function piece used to evealuate self at x0.

        EXAMPLES:
	    sage: x = PolynomialRing(RationalField()).gen()
            sage: f1 = lambda z:1
            sage: f2 = 1-x
            sage: f3 = lambda y:exp(y)
            sage: f4 = lambda t:sin(2*t)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f.which_function(3/2)
	    -x + 1

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

    def integral(self):
        """
        Returns the definite integral (as computed by maxima)
        $\sum_I \int_I self|_I$, as I runs over the intervals
        belonging to self.

        EXAMPLES:
            sage: f1 = lambda x:1
            sage: f2 = lambda x:1-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.integral()
            1/2
	    sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.integral()
            (pi/2)
        """
        maxima = sage.interfaces.all.maxima
        x = PolynomialRing(RationalField()).gen()
        ints = [maxima('%s'%p[1](x)).integral('x', p[0][0], p[0][1]) \
                 for p in self.list()]
        return sage_eval(str(sum(ints)).replace("%",""))

    def plot(self):
        """
        Returns the plot of self.

        EXAMPLES:
            sage: f1 = lambda x:1
	    sage: f2 = lambda x:1-x
            sage: f3 = lambda x:exp(x)
            sage: f4 = lambda x:sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: P = f.plot()

	Remember: to view this, type P.save("<path>/myplot.png") and then open it in a graphics
	viewer such as gimp.
        """
        plot = sage.plot.plot.plot
        plots = [plot(p[1], p[0][0], p[0][1], rgbcolor=(1,0,0) ) \
                 for p in self.list()]
        n = len(plots)
	P = plots[0]
	for i in range(1,n):
	    P = P + plots[i]
	return P

    def fourier_series_cosine_coefficient(self,n,L):
        """
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
            (1/(pi^2))
	    sage: f = lambda x:x^2
            sage: f = Piecewise([[(-pi,pi),f]])
            sage: f.fourier_series_cosine_coefficient(2,pi)
            1
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_cosine_coefficient(5,pi)
            (-3/(5*pi))

        """
        maxima = sage.interfaces.all.maxima
        x = PolynomialRing(RationalField()).gen()
        ints = []
        for p in self.list():
            fcn = '(%s)*cos('%p[1](x) + 'pi*x*%s/%s)/%s'%(n,L,L)
            fcn = fcn.replace("pi","%"+"pi")
	    a = str(p[0][0]).replace("pi","%"+"pi")
	    b = str(p[0][1]).replace("pi","%"+"pi")
	    cmd = "integrate("+fcn+", x, %s, %s )"%(a, b)
	    int = maxima(cmd)
            ints.append(int)
        ans = sum(ints)
        return sage_eval(str(ans).replace("%",""))

    def fourier_series_sine_coefficient(self,n,L):
        """
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
        x = PolynomialRing(RationalField()).gen()
        ints = []
        for p in self.list():
            fcn = '(%s)*sin('%p[1](x) + 'pi*x*%s/%s)/%s'%(n,L,L)
            fcn = fcn.replace("pi","%"+"pi")
	    a = str(p[0][0]).replace("pi","%"+"pi")
	    b = str(p[0][1]).replace("pi","%"+"pi")
	    cmd = "integrate("+fcn+", x, %s, %s )"%(a, b)
	    int = maxima(cmd)
            ints.append(int)
        ans = sum(ints)
        return sage_eval(str(ans).replace("%",""))

    def fourier_series_partial_sum(self,N,L):
        """
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
            '1/3 + ((-4/(pi^2))*cos(1*pi*x/1) + 0*sin(1*pi*x/1)) + ((1/(pi^2))*cos(2*pi*x/1) + 0*sin(2*pi*x/1))'
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.fourier_series_partial_sum(3,pi)
            '1/4 + ((-3/pi)*cos(1*pi*x/pi) + (1/pi)*sin(1*pi*x/pi)) + (0*cos(2*pi*x/pi) + (-3/pi)*sin(2*pi*x/pi))'

        """
        a0 = self.fourier_series_cosine_coefficient(0,L)
        A = [str(self.fourier_series_cosine_coefficient(n,L))+"*cos(%s*pi*x/%s)"%(n,L) for n in range(1,N)]
        B = [str(self.fourier_series_sine_coefficient(n,L))+"*sin(%s*pi*x/%s)"%(n,L) for n in range(1,N)]
        FS =  ["("+A[i] +" + " + B[i]+")" for i in range(0,N-1)]
        sumFS = str(a0/2)+" + "
        for s in FS:
            sumFS = sumFS+s+ " + "
        return sumFS[:-3]

    def plot_fourier_series_partial_sum(self,N,L,xmin,xmax):
        """
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
            sage: P = f.plot_fourier_series_partial_sum(3,pi,-5,5)
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: P = f.plot_fourier_series_partial_sum(15,pi,-5,5)

	Remember, to view this type P.save("<path>/myplot.png") and then open it in a graphics
	viewer such as gimp.

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
            yi = ff.replace("pi",str(RR(pi)))
            yi = sage_eval(yi.replace("x",str(xi)))
            pts.append([xi,yi])
        return line(pts)

    def fourier_series_value(self,x,L):
        """
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
            0.91294525072762767
            sage: f.fourier_series_value(20,10)
            1
            sage: f.fourier_series_value(30,10)
            0.91294525072762767
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
        """
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
            (-4/(9*(pi^2)))
            sage: f1 = lambda x:-1
            sage: f2 = lambda x:2
            sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
            sage: f.cosine_series_coefficient(2,pi)
            0
            sage: f.cosine_series_coefficient(3,pi)
            (2/pi)
            sage: f.cosine_series_coefficient(111,pi)
            (2/(37*pi))

        """
	maxima = sage.interfaces.all.maxima
        x = PolynomialRing(RationalField()).gen()
        ints = []
        for p in self.list():
            fcn = '2*(%s)*cos('%p[1](x) + 'pi*x*%s/%s)/%s'%(n,L,L)
            fcn = fcn.replace("pi","%"+"pi")
	    a = str(p[0][0]).replace("pi","%"+"pi")
	    b = str(p[0][1]).replace("pi","%"+"pi")
	    cmd = "integrate("+fcn+", x, %s, %s )"%(a, b)
	    int = maxima(cmd)
            ints.append(int)
        ans = sum(ints)
        return sage_eval(str(ans).replace("%",""))

    def sine_series_coefficient(self,n,L):
        """
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
            (4/(3*pi))

        """
	maxima = sage.interfaces.all.maxima
        x = PolynomialRing(RationalField()).gen()
        ints = []
        for p in self.list():
            fcn = '2*(%s)*sin('%p[1](x) + 'pi*x*%s/%s)/%s'%(n,L,L)
            fcn = fcn.replace("pi","%"+"pi")
	    a = str(p[0][0]).replace("pi","%"+"pi")
	    b = str(p[0][1]).replace("pi","%"+"pi")
	    cmd = "integrate("+fcn+", x, %s, %s )"%(a, b)
	    int = maxima(cmd)
            ints.append(int)
        ans = sum(ints)
        return sage_eval(str(ans).replace("%",""))

    def laplace_transform(self,var = "s",latex_output=0):
        """
        Returns the laplace transform of self, as a function of var.
        We assume that a piecewise function is 0 outside of its domain
        and that the left-most endpoint of the domain is 0.

        EXAMPLES:
            sage: f1 = lambda x:1
            sage: f2 = lambda x:1-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.laplace_transform()
            '1/s - e^-s/s + (s + 1)*e^-(2*s)/s^2 - e^-s/s^2'
            sage: f.laplace_transform("w",latex_output=1)
            ' - \\frac{e^{ - w}}{w} - \\frac{e^{ - w}}{w^2} + \\frac{\\left(w + 1\\right)e^{ - 2w}}{w^2} + \\frac{1}{w}'
            sage: f.laplace_transform("w",True)
            ' - \\frac{e^{ - w}}{w} - \\frac{e^{ - w}}{w^2} + \\frac{\\left(w + 1\\right)e^{ - 2w}}{w^2} + \\frac{1}{w}'
            sage: f.laplace_transform("w")
            '1/w - e^-w/w + (w + 1)*e^-(2*w)/w^2 - e^-w/w^2'

        """
        maxima = sage.interfaces.all.maxima
        x = PolynomialRing(RationalField()).gen()
        ints = []
        for p in self.list():
            fcn = '(%s)*exp(-%s*x)'%(p[1](x),var)
            ints.append(maxima(fcn).integral('x', p[0][0], p[0][1]))
        ans = ""
        ans_latex = ""
        for i in range(len(ints)-1):
            ans = ans+str(ints[i]).replace("%","")+" + "
            ans_latex = ans_latex+str(ints[i])+" + "
        ans = ans+str(ints[len(ints)-1]).replace("%","")
        ans_latex = ans_latex + str(ints[len(ints)-1])

        if latex_output == 0:
            return ans
        if latex_output == 1:
            ans0 = maxima.eval("tex("+ans_latex+")")
            ans0 = ans0.replace("$$","")
            ans0 = ans0.replace("false","")
            return ans0

    def __add__(self,other):
	"""
	Returns the piecewise defined function which is the sum of
	self and other.

	EXAMPLES:
	    sage: x = PolynomialRing(RationalField()).gen()
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
            Piecewise defined function with 5 parts, [[[0, 1], x - 1], [[1, 2], -1], [[2, 3], 3*x - 2], [[3, 5], 8], [[5, 10], 5]]

        Note that in this case the functions must be defined using polynomial
	expressions *not* using the lambda notation.
	"""
	self_endpts = self.end_points()     ## we assume these start
	other_endpts = other.end_points()   ## and end at the same point
	f = self.functions()
	g = other.functions()
	endpts = list(set(self_endpts).union(set(other_endpts)))
	N = len(list(endpts))
	fcn = []
	for j in range(N-1):
	    x0 = endpts[j+1]
	    fcn.append([[endpts[j],endpts[j+1]],self.which_function(x0)+other.which_function(x0)])
	return Piecewise(fcn)

    def __mul__(self,other):
	"""
	Returns the piecewise defined function which is the product of
	one piecewise function (self) with another one (other).

	EXAMPLES:
	    sage: x = PolynomialRing(RationalField()).gen()
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
	    Piecewise defined function with 5 parts, [[[0, 1], x - 2], [[1, 2], -x^2 + 3*x - 2], [[2, 3], 2*x^2 - 4*x], [[3, 5], -x^2 + 12*x - 20], [[5, 10], -x^2 + 15*x - 50]]
            sage: g*(11/2)
            Piecewise defined function with 2 parts, [[[0, 5], 11/2*x - 11], [[5, 10], 11/2*x - 55/2]]

        Note that in this method the functions must be defined using polynomial
	expressions *not* using the lambda notation.
	"""
        R = PolynomialRing(RationalField())
        fcn = []
        if isinstance(self,Rational) or isinstance(self,Integer):    ## needed for scalar multiplication
            endpts = other.end_points()
            N = len(list(endpts))
            for j in range(N-1):
	        x0 = endpts[j+1]
	        fcn.append([[endpts[j],endpts[j+1]],R(self)*other.which_function(x0)])
            return Piecewise(fcn)
        if isinstance(other,Rational) or isinstance(other,Integer):    ## needed for scalar multiplication
            endpts = self.end_points()
            N = len(list(endpts))
            for j in range(N-1):
	        x0 = endpts[j+1]
	        fcn.append([[endpts[j],endpts[j+1]],R(other)*self.which_function(x0)])
            return Piecewise(fcn)
	self_endpts = self.end_points()     ## we assume these start
	other_endpts = other.end_points()   ## and end at the same point
	f = self.functions()
	g = other.functions()
	endpts = list(set(self_endpts).union(set(other_endpts)))
	N = len(list(endpts))
	for j in range(N-1):
	    x0 = endpts[j+1]
	    fcn.append([[endpts[j],endpts[j+1]],self.which_function(x0)*other.which_function(x0)])
	return Piecewise(fcn)

Piecewise = PiecewisePolynomial  ## added so that functions/all.py does not need to be changed