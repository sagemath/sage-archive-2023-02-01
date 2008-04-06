r"""nodoctest
Elementary Functions

!!!TEMPORARILY TOTALLY BROKEN!!!

Call a univariate function an "elementary function" if it can be
written as a sum of functions of the form "polynomial times and
exponential times a sine or a cosine", ie
$$
       f(x) = (a_nx^n + ... + a_0)\exp(ax)\cos(bx)
$$
or
$$
       f(x) = (a_nx^n + ... + a_0)\exp(ax)\sin(cx).
$$

These arise naturally in solving constant coefficient ordinary
differential equations using the method of undetermined coefficients.

The set E of elementary functions is an algebra over RR.  If D is
differentiation and A = RR[D] is the polynomial ring in D over RR, let
us define a smooth function f to be *finite* if the vector space A(f)
is finite dimensional.

Theorem: E is the algebra of all finite functions.

\begin{verbatim}
Methods implemented:
    * addition
    * multiplication (right scalar multiplication only)
    * positive integer powers
    * integration (definite -, returns a RR; indefinite -, usually
      returns a *string* but in simple cases returns an instance
      of ElementaryFunction)
    * differentiation (returns an instance of ElementaryFunction)
    * laplace transform (returns a *string*)
    * __call__ (ie, evaluation) for ElementaryFunction and
      for ElementaryFunctionRing
    * _coercise_ for ElementaryFunctionRing
    * latex
    * simplify (combines terms using the distributive law)
    * print methods
    * parent method for ElementaryFunctions
    * unary - (ie, multiply by -1)
    * desolve - solves a constant coefficient DE of the form
      Phi(D)y(x) = e(x), where Phi(D) is a
      constant coefficient polynomial in D and
      e(x) is an elementary function (using the method of
      Laplace transforms). Note: returns a *string*.

TODO:
    * extend "integrate" so that it always returns a class instance
      instead of a string
    * extend "desolve" so that it always returns a class instance
      instead of a string
    * extend piecewise piecewise functions to allow ElementaryFunctions.
\end{verbatim}

AUTHOR: David Joyner (2006-06)
        -- (2006-09) - minor bug fixes
        -- (2006-11) -- completely rewrite

     REFERENCE:
        * Abramowitz and Stegun: Handbook of Mathematical Functions,
          http://www.math.sfu.ca/~cbm/aands/
"""

################################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdj@usna.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
################################################################################

import copy
import sage.plot.plot
import sage.interfaces.all
from sage.rings.polynomial.polynomial_ring import PolynomialRing
from sage.rings.rational_field import RationalField
from sage.rings.real_mpfr import RealField_constructor as RealField
from sage.misc.sage_eval import sage_eval
from sage.rings.all import QQ, RR, RDF, ZZ
import sage.rings.commutative_ring as commutative_ring
import sage.rings.ring as ring
from constants import e as E
from functions import sin,cos,expo

from sage.structure.parent_base import ParentWithBase
from sage.structure.element import CommutativeRingElement, generic_power


############################################################
##   the algebra generators of the elementary functions
############################################################

def cosine(a,var):
    """
    Cosine -- one of the algebra generators of the space of elementary
    functions in the given variable.

    EXAMPLES:
        sage: R = ElementaryFunctionRing(QQ,"t")
        sage: t = R.polygen()
        sage: f = cosine(2,t); f
        Elementary function (1)cos(2*t)
        sage: f.diff()
        Elementary function (-2)sin(2*t)
        sage: f.int([])
        Elementary function (1/2)sin(2*t)
        sage: f.int([0,1])
        0.45464871341284085

    """
    if is_Polynomial(var):
        var = ElementaryFunctionRing(QQ,var).var()
    R = ElementaryFunctionRing(QQ,var)
    return R.cosine(a)

def sine(a,var):
    """
    An algebra generator of the space of "elementary functions".

    EXAMPLES:
        sage: R = ElementaryFunctionRing(QQ,"t")
        sage: t = R.polygen()
        sage: f = sine(2,t); f
        Elementary function (1)sin(2*t)
        sage: f.diff()
        Elementary function (2)cos(2*t)
        sage: f.int([])
        Elementary function (-1/2)cos(2*t)
        sage: f.int([0,1])
        0.70807341827357118
        sage: f^2
        Elementary function  1/2  + (-1/2)cos(4*t)

    """
    if is_Polynomial(var):
        var = ElementaryFunctionRing(QQ,var).var()
    R = ElementaryFunctionRing(QQ,var)
    return R.sine(a)

def exponential(a,var):
    """
    An algebra generator of the space of "elementary functions".

    EXAMPLES:
        sage: R = ElementaryFunctionRing(QQ,"t"); R
        ElementaryFunctionRing over Rational Field in t
        sage: t = R.polygen(); t
        t
        sage: f = exponential(2,t); f
        Elementary function (1)exp(2*t)
        sage: f.diff()
        Elementary function (2)exp(2*t)
        sage: f.int([])
        Elementary function (1/2)exp(2*t)
        sage: f.latex()
        '(1)e^{2t}\\cos(0t)'
        sage: f(1)       ## somewhat randomish output
        7.3890560989306495
        sage: f.laplace_transform("s")
        '1/(s - 2)'
        sage: f^2
        Elementary function (1)exp(4*t)

    """
    if is_Polynomial(var):
        var = ElementaryFunctionRing(QQ,var).var()
    R = ElementaryFunctionRing(QQ,var)
    return R.exponential(a)


############################################################



class ElementaryFunctionRing(commutative_ring.CommutativeRing):
    """
    Implements ring of elementary functions.

    EXAMPLES:
       sage: ElementaryFunctionRing(RR,"x")
       ElementaryFunctionRing over Real Field with 53 bits of precision in x
       sage: R = ElementaryFunctionRing(QQ,"t"); R
       ElementaryFunctionRing over Rational Field in t
       sage: R = ElementaryFunctionRing(QQ,"t"); R
       ElementaryFunctionRing over Rational Field in t
       sage: t = R.polygen(); t
       t
       sage: g = cosine(2,t); g
       Elementary function (1)cos(2*t)
       sage: x = g.polygen()
       sage: ElementaryFunction(R,[(-2*x^0,3,5,0),(x,-1,0,2)])
       Elementary function (-2)exp(3*t)cos(5*t) + (t)exp(-t)sin(2*t)

    """
    def __init__(self, base_ring, varname='x'):
        ParentWithBase.__init__(self, base_ring)
        self._varname = varname

    def var(self):
        return self._varname

    def _repr_(self):
        return "ElementaryFunctionRing over %s in %s"%(self.base_ring(),self.var())

    def polygen(self): ## added for consistency of notation with PolynomialRing
        return PolynomialRing(self.base_ring(),str(self._varname)).gen()

    def __call__(self,z):
        """
        Helps with coersion.
        """
        x = self.polygen()
        F = self.base_ring()
        if z is sin:
            return ElementaryFunction(self, [(x**0,0,0,1)])
        if z is cos:
            return ElementaryFunction(self, [(x**0,0,1,0)])
        if z is expo:
            return ElementaryFunction(self, [(x**0,1,0,0)])
        try:
            if z in RR:
                return ElementaryFunction(self, [(z*x**0,0,0,0)])
            P = PolynomialRing(F,str(x))
            return ElementaryFunction(self, [(P(z),0,0,0)])
        except TypeError:
            raise TypeError,"Not coercible."
        raise TypeError,"Not coercible."

    def _coerce_impl(self,x):
        return self._coerce_try(x, [RR])

    def poly(self,p):
        """
        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: p = x^2 + 1
            sage: R = ElementaryFunctionRing(QQ,"x"); R
            ElementaryFunctionRing over Rational Field in x
            sage: R.poly(p)
            Elementary function  x^2 + 1
            sage: R.cosine(2)
            Elementary function (1)cos(2*x)
            sage: R.cosine(2)*R.poly(p)
            Elementary function (x^2 + 1)cos(2*x)
            sage: R.cosine(2)*R.poly(p)*R.exponential(3)
            Elementary function (x^2 + 1)exp(3*x)cos(2*x)

        """
        return ElementaryFunction(self, [(p,0,0,0)])

    def exponential(self,a):
        """
        EXAMPLES:
            sage: R = ElementaryFunctionRing(QQ,"t"); R
            ElementaryFunctionRing over Rational Field in t
            sage: R.exponential(2)
            Elementary function (1)exp(2*t)

        """
        var = self.polygen()
        base_ring = self.base_ring()
        x = PolynomialRing(base_ring,str(var)).gen()
        return ElementaryFunction(self, [(x**0,a,0,0)])

    def sine(self,a):
        """
        EXAMPLES:
            sage: R = ElementaryFunctionRing(QQ,"t"); R
            ElementaryFunctionRing over Rational Field in t
            sage: R.sine(2)
            Elementary function (1)sin(2*t)
        """
        var = self.polygen()
        base_ring = self.base_ring()
        x = PolynomialRing(base_ring,str(var)).gen()
        return ElementaryFunction(self, [(x**0,0,0,a)])

    def cosine(self,a):
        """
        EXAMPLES:
            sage: R = ElementaryFunctionRing(QQ,"t"); R
            ElementaryFunctionRing over Rational Field in t
            sage: R.cosine(2)
            Elementary function (1)cos(2*t)

        """
        var = self.polygen()
        base_ring = self.base_ring()
        x = PolynomialRing(base_ring,str(var)).gen()
        return ElementaryFunction(self, [(x**0,0,a,0)])

# The default global elementary function ring.
EFR = ElementaryFunctionRing(RDF, 'x')

def ElementaryFunction(list_of_fcns):
   """
   INPUT:
        list_of_fcns --
        \code{list_of_fcns} is a list of 4-tuples (p,a,b,c),
        representating p*exp(a*x)*cos(b*x) (if c=0) or
        p*exp(a*x)*sin(b*x) (if b=0). Here a,b,c are
        in RR and p is in RR[x].
        We assume that b=0 or c=0.

   OUTPUT:
        A double precision elementary function.
   """
   return ElementaryFunction_class(EFR, list_of_fcns)


class ElementaryFunction_class(CommutativeRingElement):
    def __init__(self, parent, list_of_fcns):
        r"""
        \code{list_of_fcns} is a list of 4-tuples (p,a,b,c),
        representating p*exp(a*x)*cos(b*x) (if c=0) or
        p*exp(a*x)*sin(b*x) (if b=0). Here a,b,c are
        in RR and p is in RR[x].

        We assume that b=0 or c=0.
        """
        if not isinstance(list_of_fcns, list):
           raise TypeError, "list_of_fcns must be a list of 4-tuples"
        self._length = len(list_of_fcns)
        self._list = list_of_fcns
        self._set_parent(parent)
        for x in list_of_fcns:
           if not (isinstance(x, tuple) and len(x) == 4):
              raise TypeError, "list_of_fcns must be a list of 4-tuples"
        if is_Polynomial(list_of_fcns[0][0]):
            self._var = list_of_fcns[0][0].parent().gen()
        else:
            F = list_of_fcns[0][0].base_ring()
            var = self.variable()                        #### todo: wrong -- FIX THIS
            R = PolynomialRing(F,str(var))
            self._var = R.gen()

    def plot(self,*args, **kwds):
        """
        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: g = ElementaryFunction(self, [(x^2-x+1,3,5,0),(x^3-2,-1,0,2)])
            sage: P = plot(g,-1,1)

        Now type show(P) to view this.
        """
        return sage.plot.plot.plot(self, *args, **kwds)

    def list(self):
        return self._list

    def variable(self):
        """

        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: g = ElementaryFunction(self, [(x^2-x+1,3,5,0),(x^3-2,-1,0,2)])
            sage: g.variable()
            x
        """
        return self._var

    def polygen(self):
        var = self.variable()
        F = self.variable().base_ring()
        return PolynomialRing(F,str(var)).gen()

    def parent(self):
        var = self.variable()
        F = self.variable().base_ring()
        R = ElementaryFunctionRing(F,str(var))
        return R

    def base_ring(self):
        return self.parent().base_ring()

    def _neg_(self):
        """
        EXAMPLES:
            sage: R = ElementaryFunctionRing(QQ,"t")
            sage: t = R.polygen()
            sage: g = ElementaryFunction(self, [((4*t^3+t),3,5,0),(t,-1,0,2)])
            sage: -g
            Elementary function (-4*t^3 - t)exp(3*t)cos(5*t) + (-t)exp(-t)sin(2*t)
        """
        return self*(-1)

    def exponential(self,a):
        """
        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: f = ElementaryFunction([(x^0,0,0,0)])
            sage: f.exponential(2)
            Elementary function (1)exp(2*x)
        """
        var = self.variable()
        return ElementaryFunction([(var**0,a,0,0)])

    def cosine(self,a):
        """
        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: f = ElementaryFunction([(x^0,0,0,0)])
            sage: print f.cosine(2)
            (1)*cos(2*x)
        """
        var = self.variable()
        return ElementaryFunction([(var**0,0,a,0)])

    def sine(self,a):
        """
        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: f = ElementaryFunction([(x^0,0,0,0)])
            sage: f.sine(2)
            Elementary function (1)sin(2*x)
        """
        var = self.variable()
        return ElementaryFunction([(var**0,0,0,a)])

    def length(self):
        return self._length

    def _repr_(self):
        """

        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: g = ElementaryFunction([(x^2-x+1,3,5,0),(x^3-2,-1,0,2)])
            sage: g
            Elementary function (x^2 - x + 1)exp(3*x)cos(5*x) + (x^3 - 2)exp(-x)sin(2*x)
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: f = ElementaryFunction([(x^0,0,0,0)])
            sage: f
            Elementary function  1
            sage: print f
            1
            sage: R = ElementaryFunctionRing(QQ,"t"); R
            ElementaryFunctionRing over Rational Field in t
            sage: t = R.polygen(); t
            t
            sage: f = cosine(2,t); f
            Elementary function (1)cos(2*t)
            sage: f.int([])
            Elementary function (1/2)sin(2*t)

        """
        var = str(self.variable())
        fcn = ""
        nonzero = "false"  ## silly trick
        for f in self.list():
            p = f[0] # the polynomial
            a = f[1] # the exponential coeff
            b = f[2] # the cosine coeff
            c = f[3] # the sine coeff
            if p==0*var:
                pass
            else:
                nonzero = "true"
                if c==0 and b!=0 and a==1:
                    fcn = fcn + " + " + "(" + str(p) + ")" + "exp(" + var + ")cos(" + str(b) + "*" + var+ ")"
                elif c==0 and b!=0 and a==-1:
                    fcn = fcn + " + " + "(" + str(p) + ")" + "exp(-" + var + ")cos(" + str(b) + "*" + var+ ")"
                elif c==0 and b==0 and a==-1:
                    fcn = fcn + " + " + "(" + str(p) + ")" + "exp(-" + var + ")"
                elif c==0 and b==0 and a==1:
                    fcn = fcn + " + " + "(" + str(p) + ")" + "exp(" + var + ")"
                elif c==0 and b==0 and a==0:
                    fcn = fcn + " + " + " " + str(p) + " "
                elif c!=0 and b==0 and a==0:
                    fcn = fcn + " + " + "(" + str(p) + ")sin(" + str(c) + "*" + var+ ")"
                elif c==0 and b!=0 and a==0:
                    fcn = fcn + " + " + "(" + str(p) + ")cos(" + str(b) + "*" + var+ ")"
                elif c==0 and b==0:
                    fcn = fcn + " + " + "(" + str(p) + ")" + "exp(" + str(a) + "*" + var + ")"
                elif c==0 and a!=-1 and a!=1:
                    fcn = fcn + " + " + "(" + str(p) + ")" + "exp(" + str(a) + "*" + var + ")cos(" + str(b) + "*" + var+ ")"
                elif b==0 and a==-1:
                    fcn = fcn + " + " + "(" + str(p) + ")" + "exp(-" + var + ")sin(" + str(c) + "*" + var+ ")"
                elif b==0 and a==1:
                    fcn = fcn + " + " + "(" + str(p) + ")" + "exp(" + var + ")sin(" + str(c) + "*" + var+ ")"
                else:
                    fcn = fcn + " + " + "(" + str(p) + ")" + "exp(" + str(a) + "*" + var + ")sin(" + str(c) + "*" + var+ ")"
        if nonzero == "true":
            return 'Elementary function ' + fcn[3:]
        else:
            return 'Elementary function 0'

    # TODO: Why is this here?  Is it the same as _repr_ above?!!?  -- William Stein
    def __str__(self):
        """

        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: g = ElementaryFunction([(x^2-x+1,3,5,0),(x^3-2,-1,0,2)])
            sage: print g
            (x^2 - x + 1)*exp(3*x)*cos(5*x) + (x^3 - 2)*exp(-x)*sin(2*x)
        """
        var = str(self.variable())
        fcn = ""
        nonzero = "false"  ## silly trick
        for f in self.list():
            p = f[0] # the polynomial
            a = f[1] # the exponential coeff
            b = f[2] # the cosine coeff
            c = f[3] # the sine coeff
            # cases a,b,c = 0, 1, -1
            if p==0*var:
                pass
            else:
                nonzero = "true"
                if c==0 and b!=0 and a==1:
                    fcn = fcn + " + " + "(" + str(p) + ")*" + "exp(" + var + ")*cos(" + str(b) + "*" + var+ ")"
                elif c==0 and b!=0 and a==-1:
                    fcn = fcn + " + " + "(" + str(p) + ")*" + "exp(-" + var + ")*cos(" + str(b) + "*" + var+ ")"
                elif c==0 and b==0 and a==-1:
                    fcn = fcn + " + " + "(" + str(p) + ")*" + "exp(-" + var + ")"
                elif c==0 and b==0 and a==1:
                    fcn = fcn + " + " + "(" + str(p) + ")*" + "exp(" + var + ")"
                elif c==0 and b==0 and a==0:
                    fcn = fcn + " + " + " " + str(p) + " "
                elif c!=0 and b==0 and a==0:
                    fcn = fcn + " + " + "(" + str(p) + ")*sin(" + str(c) + "*" + var+ ")"
                elif c==0 and b!=0 and a==0:
                    fcn = fcn + " + " + "(" + str(p) + ")*cos(" + str(b) + "*" + var+ ")"
                elif c==0 and b==0:
                    fcn = fcn + " + " + "(" + str(p) + ")*" + "exp(" + str(a) + "*" + var + ")"
                elif c==0 and a!=-1 and a!=1:
                    fcn = fcn + " + " + "(" + str(p) + ")*" + "exp(" + str(a) + "*" + var + ")*cos(" + str(b) + "*" + var+ ")"
                elif b==0 and a==-1:
                    fcn = fcn + " + " + "(" + str(p) + ")*" + "exp(-" + var + ")*sin(" + str(c) + "*" + var+ ")"
                elif b==0 and a==1:
                    fcn = fcn + " + " + "(" + str(p) + ")*" + "exp(" + var + ")*sin(" + str(c) + "*" + var+ ")"
                else:
                    fcn = fcn + " + " + "(" + str(p) + ")*" + "exp(" + str(a) + "*" + var + ")*sin(" + str(c) + "*" + var+ ")"
        if nonzero == "true":
            return fcn[3:]
        else:
            return '0'


    def latex(self):
	"""
	EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: g = ElementaryFunction([(x^2-x+1,3,5,0),(x^3-2,-1,0,2)])
            sage: g.variable()
            x
            sage: g.latex()
            '(x^2 - x + 1)e^{3x}\\cos(5x) + (x^3 - 2)e^{-1x}\\sin(2x)'

	"""
        var = self.variable()
        fcn = ""
        nonzero = "false"  ## silly trick
        for f in self.list():
            p = f[0] # the polynomial
            a = f[1] # the exponential coeff
            b = f[2] # the cosine coeff
            c = f[3] # the sine coeff
            if p==0 and a==0 and b==0 and c==0:
                return '0'
            if p==0*var:
                pass
            else:
                nonzero = "true"
                if c==0:
                    fcn = fcn + " + " + "(" + str(p) + ")" + "e^{"+str(a) + var + "}\cos(" + str(b) + var+ ")"
                else:
                    fcn = fcn + " + " + "(" + str(p) + ")" + "e^{"+str(a) + var + "}\sin(" + str(c) + var+ ")"
        if nonzero == "true":
            return fcn[3:]
        else:
            return '0'

    def __call__(self,x0):
        """
        Evaluates self at x0.

        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: g = ElementaryFunction([(x^2-x+1,3,5,0),(x^3-2,-1,0,2)])
            sage: g(1)       ## somewhat randomish output
            5.3629954705944760
            sage: exp(3)*cos(5)-exp(-1)*sin(2)
            5.3629954705944769

        """
        val = RR(0)
        x = RR(x0)
        for f in self.list():
            p = f[0] # the polynomial
            a = f[1] # the exponential coeff
            b = f[2] # the cosine coeff
            c = f[3] # the sine coeff
            if c!=0:
                val = val + p(x)*expo(a*x)*cos(b*x)*sin(c*x)
            else:
                val = val + p(x)*expo(a*x)*cos(b*x)
        return val

    def exact_value(self,a):
        """
        Returns symbolic value (as string) if a is not in RR

        """
        if type(axs)==real_mpfr.RealNumber:
            return self(a)
        maxima = sage.interfaces.all.maxima
        fcn = self.__str__()
        var = self.variable()
        maxima.eval("expr: " + fcn)
        a = maxima.eval("ev(expr," + var +"=" + str(x0)+")")
        return a.replace("%","")

    def integrate(self,intvl=[0,1]):
        r"""
        Returns the integral (as computed by maxima).

        WARNING:
           Maxima asks the sign of the variable of integration to
           evaluate the expression. We assume that this variable is > 0.

           The indefinite integration option is raised by
           using the "intvl=[]" option. This returns a *string*.
           Definite integrals are the default.

        TODO:
           Make indefinite integration return an instance of ElementaryFunction.
           This is not easy using Maxima's "args" command, since the
           ordering poly*exp*trig is not preserved.

        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: g = ElementaryFunction([(-2*x^0,3,5,0),(x,-1,0,2)]) # note: x^0 is required
            sage: g
            Elementary function (-2)exp(3*x)cos(5*x) + (x)exp(-x)sin(2*x)
            sage: g.int([])
            '-e^(3*x)*(5*sin(5*x) + 3*cos(5*x))/17 - e^-x*((5*x - 3)*sin(2*x) + (10*x + 4)*cos(2*x))/25'
            sage: g = ElementaryFunction([((4*x^3+x),3,5,0),(x,-1,0,2)])
            sage: g
            Elementary function (4*x^3 + x)exp(3*x)cos(5*x) + (x)exp(-x)sin(2*x)
            sage: g.int([])
            '2*((24565*x^3 - 13005*x^2 + 255*x + 720)*e^(3*x)*sin(5*x) + (14739*x^3 + 6936*x^2 - 5049*x + 483)*e^(3*x)*cos(5*x))/83521 + ((85*x - 15)*e^(3*x)*sin(5*x) + (51*x + 8)*e^(3*x)*cos(5*x))/578 - e^-x*((5*x - 3)*sin(2*x) + (10*x + 4)*cos(2*x))/25'
            sage: g.int()
            -5.0045254279249374
            sage: g.int([0,1])   ## the interval [0,1] is the default.
            -5.0045254279249374
            sage: R = ElementaryFunctionRing(QQ,"t"); R
            ElementaryFunctionRing over Rational Field in t
            sage: t = R.polygen(); t
            t
            sage: f = cosine(2,t); f
            Elementary function (1)cos(2*t)
            sage: f.int([])
            Elementary function (1/2)sin(2*t)
        """
        maxima = sage.interfaces.all.maxima
        fcn = self.__str__()
        var = self.variable()
        if self.length()==1 and intvl==[]:
            f = self.list()[0]
            p = var**0*f[0]
            a = f[1]
            b = f[2]
            c = f[3]
            if p.derivative()==0 and c==0: # ie, if f = const*exp*cos
                intf = [(p*a/(a**2 + b**2),a,b,0)]
                intf.append((p*b/(a**2 + b**2),a,0,b))
                return ElementaryFunction(intf).simplify()
            elif p.derivative()==0 and b==0: # ie, if f = const*exp*sin
                intf = [(-p*c/(a**2 + c**2),a,c,0)]
                intf.append((p*a/(a**2 + c**2),a,0,c))
                return ElementaryFunction(intf).simplify()
            else:
                pass
        maxima.eval("expr: " + fcn)
        if intvl==[]:
            maxima.eval("assume("+ var + " > 0)")   # maxima requires knowing the sign of var
            a = maxima.eval("integrate(expr," + var + ")")
            return a.replace("%","")
        else:
            a = maxima.eval("integrate(expr," + var + "," + str(intvl[0]) + "," + str(intvl[1]) + ")")
            return RR(sage_eval(a.replace("%","")))

    int = integrate

    def differentiate(self,verbose=None):
        r"""
        Returns the definite integral (as computed by maxima)
        $\sum_I \int_I self|_I$, as I runs over the intervals
        belonging to self.

        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: g = ElementaryFunction([(-2*x^0,3,5,0),(x,-1,0,2)]) # note: x^0 is required
            sage: g.variable()
            x
            sage: g
            Elementary function (-2)exp(3*x)cos(5*x) + (x)exp(-x)sin(2*x)
            sage: g.diff()
            Elementary function (-6)exp(3*x)cos(5*x) + (10)exp(3*x)sin(5*x) + (-x + 1)exp(-x)sin(2*x) + (2*x)exp(-x)cos(2*x)
            sage: g.diff(verbose=True)
            <BLANKLINE>
            (Elementary function (-6)exp(3*x)cos(5*x) + (10)exp(3*x)sin(5*x) + (-x + 1)exp(-x)sin(2*x) + (2*x)exp(-x)cos(2*x),
            '10*e^(3*x)*sin(5*x) - 6*e^(3*x)*cos(5*x) - x*e^-x*sin(2*x) + e^-x*sin(2*x) + 2*x*e^-x*cos(2*x)')

        """
        maxima = sage.interfaces.all.maxima
        f_list = self.list()
        df_list = []
        for f in f_list:
            p = f[0]
            a = f[1]
            b = f[2]
            c = f[3]
            if c==0:
                df_list.append((p.derivative()+a*p,a,b,0))
                df_list.append((-b*p,a,0,b))
            else:
                df_list.append((p.derivative()+a*p,a,0,c))
                df_list.append((c*p,a,c,0))
        df = ElementaryFunction(df_list)
        fcn = self.__str__()
        var = str(self.variable())
        maxima.eval("expr: " + fcn)
        a = maxima.eval("diff(expr," + var +")")
        if verbose!=None:
            return df,a.replace("%","")
        else:
            return df

    diff = differentiate
    derivative = differentiate

    def laplace_transform(self,var = "s",latex_output=0):
        r"""
        Returns the laplace transform of self, as a function of var.

        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: g = ElementaryFunction([(-2*x^0,3,5,0)])
            sage: print g
            (-2)*exp(3*x)*cos(5*x)
            sage: g.laplace_transform()
            '-2*(s - 3)/(s^2 - 6*s + 34)'

        """
        maxima = sage.interfaces.all.maxima
        fcn = self.__str__()
        var = str(self.variable())
        maxima.eval("expr: " + fcn)
        a = maxima.eval("laplace(expr," + var +", s)")
        return a.replace("%","")

    def _add_(self,other):
	"""
	Returns the elementary function which is the sum of
	self and other.

	EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: f = ElementaryFunction([(x^2+1,1,1,0)])
            sage: g = ElementaryFunction([(-2*x^0,3,5,0)])
            sage: print f
            (x^2 + 1)*exp(x)*cos(1*x)
            sage: print g
            (-2)*exp(3*x)*cos(5*x)
            sage: h = f+g
            sage: print h
            (x^2 + 1)*exp(x)*cos(1*x) + (-2)*exp(3*x)*cos(5*x)

	"""
        # case 1 : self is a poly, other is elementary:
        # ...
        # case 2 : self is a elementary, other is poly:
        # ...
        # case 3 : self, other are elementary:
        return ElementaryFunction(self.list()+other.list())     ## gotta love Python:-)

    def _mul_(self,other):
	r"""
	Returns the elementary function which is the product of self
        with other. Right multiplication by a constant or a polynomial is
        possible and the program does its own coersion. (Left
        multiplication by a polynomial doesn't work poly*elem calls poly._mul_(elem),
        which is a non-existent method of Polynomial.)

	EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: f = ElementaryFunction([(x^2+1,1,1,0)])
            sage: g = ElementaryFunction([(-2*x^0,3,5,0)])
            sage: print f
            (x^2 + 1)*exp(x)*cos(1*x)
            sage: print g
            (-2)*exp(3*x)*cos(5*x)
            sage: h = f*g
            sage: print h
            (-x^2 - 1)*exp(4*x)*cos(6*x) + (-x^2 - 1)*exp(4*x)*cos(4*x)
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: p = x^2 + 1
            sage: R = ElementaryFunctionRing(QQ,"x"); R
            ElementaryFunctionRing over Rational Field in x
            sage: f = R.cosine(2)
            sage: f
            Elementary function (1)cos(2*x)
            sage: f*p
            Elementary function (x^2 + 1)cos(2*x)
            sage: R = ElementaryFunctionRing(QQ,"t"); R
            ElementaryFunctionRing over Rational Field in t
            sage: t = R.polygen()
            sage: is_Polynomial(t)
            True
            sage: g = cosine(2,t); g
            Elementary function (1)cos(2*t)
            sage: g*t
            Elementary function (t)cos(2*t)
            sage: g*2
            Elementary function (2)cos(2*t)

	"""
        self0 = self
        other0 = other
        F = self.base_ring()
        # case 1 : self is a poly or const, other is elementary (does not get called)
        if self in F:
            R = other.parent()
            var = R.polygen()
            self0 = R.poly(self*var**0)
        if is_Polynomial(self):
            R = other.parent() # = ElementaryFunctionRing(...)
            self0 = R.poly(self)
        # case 2 : self is a elementary, other is poly or const (this icase works)
        if other in F:
            R = self.parent()
            var = R.polygen()
            other0 = R.poly(other*var**0)
        if is_Polynomial(other):
            R = self.parent() # = ElementaryFunctionRing(...)
            other0 = R.poly(other)
        # case 3 : self, other are elementary:
        h = []
        f = self0.list()
        g = other0.list()
        for r in f:
            for s in g:
                if r[3]==0 and s[3]==0: ## no sine terms, use add law of cos
                    h.append((r[0]*s[0]/2,r[1] + s[1],r[2] + s[2],0))
                    h.append((r[0]*s[0]/2,r[1] + s[1],r[2] - s[2],0))
                elif r[2]==0 and s[3]==0: ## sine*cosine
                    h.append(((r[0]*s[0])/2,r[1] + s[1],0,r[3] + s[2]))
                    h.append(((r[0]*s[0])/2,r[1] + s[1],0,r[3] - s[2]))
                elif r[3]==0 and s[2]==0: ## cosine*sine
                    h.append(((r[0]*s[0])/2,r[1] + s[1],0,r[2] + s[3]))
                    h.append((-(r[0]*s[0])/2,r[1] + s[1],0,r[2] - s[3]))
                elif r[2]==0 and s[2]==0: ## no cosine terms, use add law of sin
                    h.append(((r[0]*s[0])/2,r[1] + s[1],r[3] - s[3],0))
                    h.append((-(r[0]*s[0])/2,r[1] + s[1],r[3] + s[3],0))
        return ElementaryFunction(h).simplify()


    def _sub_(self,other):
	"""
	Returns the elementary function which is the difference of
	self and other.
        """
        return (self + other*(-1)).simplify()

    def __pow__(self,n):
        """
        Implements powers.

        EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: f = ElementaryFunction([(x^0,1,0,0)])
            sage: f
            Elementary function (1)exp(x)
            sage: print f
            (1)*exp(x)
            sage: f^2
            Elementary function (1)exp(2*x)
            sage: R = ElementaryFunctionRing(QQ,"t"); R
            ElementaryFunctionRing over Rational Field in t
            sage: t = R.polygen(); t
            t
            sage: f = cosine(2,t); f
            Elementary function (1)cos(2*t)
            sage: f^2
            Elementary function (1/2)cos(4*t) +  1/2

        """
        if n < 0:
            raise TypeError,"Negative powers are not allowed"
        return generic_power(self, n)

    def simplify(self):
	r"""
	Returns the simplified list if possible.
        Check if there are two terms of the form
        p1*exp(ax)cos(bx), p2*exp(ax)cos(bx)
        or there are two terms of the form
        p1*exp(ax)sin(bx), p2*exp(ax)sin(bx).
        If so, it combines them.

	EXAMPLES:
            sage: x = PolynomialRing(QQ,"x").gen()
            sage: f = ElementaryFunction([(x^2+1,1,1,0)])
            sage: g = ElementaryFunction([(-2*x^0,1,1,0)])
            sage: print f
            (x^2 + 1)*exp(x)*cos(1*x)
            sage: print g
            (-2)*exp(x)*cos(1*x)
            sage: print f+g
            (x^2 + 1)*exp(x)*cos(1*x) + (-2)*exp(x)*cos(1*x)
            sage: print (f+g).simplify()
            (x^2 - 1)*exp(x)*cos(1*x)
            sage: R = ElementaryFunctionRing(QQ,"t"); R
            ElementaryFunctionRing over Rational Field in t
            sage: t = R.polygen()
            sage: g1 = R.cosine(2)
            sage: g2 = R.cosine(2)*(-1)
            sage: (g1+g2).simplify()
            Elementary function 0
            sage: (g1+g2).simplify().latex()
            '0'
            sage: print (g1+g2).simplify()
            0
            sage: print (g1+g2)
            (1)*cos(2*t) + (-1)*cos(2*t)
            sage: (g1+g2).latex()
            '(1)e^{0t}\\cos(2t) + (-1)e^{0t}\\cos(2t)'
            sage: g1+g2
            Elementary function (1)cos(2*t) + (-1)cos(2*t)

	"""
        fcn = self.list()
        h = copy.deepcopy(fcn)
        n = self.length()
        I = range(n)
        for i in I:
            f = list(fcn[i])
            if f[3]<0: # sin(-a*x) occurs
                f[0] = -f[0]
                f[3] = -f[3]
            if f[2]<0: # cos(-a*x) occurs
                f[2] = -f[2]
            combined = "no"   ## silly hack - must be a cleaner way.
            for j in range(i,n):
                g = list(fcn[j])
                if g[3]<0: # sin(-a*x) occurs
                    g[0] = -g[0]
                    g[3] = -g[3]
                if g[2]<0: # cos(-a*x) occurs
                    g[2] = -g[2]
                if j in I and i != j:
                    if (f[1] == g[1] and f[2] == g[2] and f[3] == g[3]):
                        h[i] = (f[0]+g[0],f[1],f[2],f[3])
                        f = h[i]
                        combined = "yes"
                        I.remove(j)
            if combined == "no":
                h[i] = f
        return ElementaryFunction([h[i] for i in I])

    def desolve(self,de_poly,soln="soln",ics=None):
        """
        Solves an ODE using laplace transforms.
        INPUT: -- de_poly is a poly in D over Q
               -- forcing_fcn is an elmt of ElemFcnRing (in t say)
               -- soln only need a string name for the solution
                  to the DE (default "soln")

        TODO:
           Implement ics -- a list of numbers representing initial conditions,
                            with symbols allowed which are represented by strings
                            (eg, f(0)=1, f'(0)=2 is ics = [0,1,2])

        The example below shows how to solve x'' - x = sin(2*t)

        EXAMPLES:
            sage: DR = PolynomialRing(QQ,"D")
            sage: D = DR.gen(); Phi = D^2 - 1
            sage: R = ElementaryFunctionRing(QQ,"t")
            sage: t = R.polygen()
            sage: g = ElementaryFunction([(1*t^0,0,0,2)])
            sage: g.desolve(Phi,"x")
            "x(t) = e^t*(5*(?at('diff(x(t),t,1),t = 0)) + 5*x(0) + 2)/10 - e^-t*(5*(?at('diff(x(t),t,1),t = 0)) - 5*x(0) + 2)/10 - sin(2*t)/5"

            sage: Phi = D^4
            sage: g = ElementaryFunction([(1*t^0,0,0,2)])
            sage: g.desolve(Phi,"x")
            "x(t) = t^3*(?at('diff(x(t),t,3),t = 0))/6 + t^2*(?at('diff(x(t),t,2),t = 0))/2 + t*(8*(?at('diff(x(t),t,1),t = 0)) - 1)/12 + t*(?at('diff(x(t),t,1),t = 0))/3 + sin(2*t)/16 + t^3/12 - t/24 + x(0)"
            sage: Phi = D^4 - 1
            sage: g.desolve(Phi,"x",[0,1,1,1,1])
            'x(t) = sin(2*t)/15 - sin(t)/3 + 11*e^t/10 - e^-t/10'

        This last command gives the solution to

         x'''' - x = 1 + sin(2t), x(0) = x'(0) =x''(0) =x'''(0) = 1.

        WARNINGS: (1) The ics will set the values of soln(0) and soln'(0) in
        Maxima, so subsequent ODEs involving these variables will have
        these initial conditions automatically imposed (unless of course
        you change the different string name for the solution).
        (2) If Maxima cannot find the solution, you will see "ilt"
        somewhere in the returned string (meaning it cannot compute the
        inverse Laplace transforms of that expression).
        """
        coeffs = de_poly.list()
        deg = len(coeffs)
        var = self.polygen()
        de = ""
        for i in range(len(coeffs)):
            de = de + "%s*diff(%s(%s),%s,%s)+"%(coeffs[i],soln,var,var,i)
        cmd = "de:" + de[:-1] + "=" + str(self) + ";"
        #print cmd
        maxima(cmd)
        vars = [var,soln]
        if ics!=None:
            d = len(ics)
            for i in range(0,d-1):
                ic = "atvalue(diff("+vars[1]+"("+vars[0]+"),"+str(vars[0])+","+str(i)+"),"+str(vars[0])+"="+str(ics[0])+","+str(ics[1+i])+")"
                maxima(ic)
                #print i,ic
        cmd = "desolve(de," + vars[1] + "(" + vars[0] + "));"
        #print cmd
        ans = str(maxima(cmd))
        #print ans
        ic = ["(at(\'diff(%s(%s),%s,%s),%s = 0)"%(soln,var,var,j,var) for j in range(deg)]
        #print ic1
        ans = ans.replace("%","")
        #print ic1, str(ans)
        for j in range(deg):  ## make output nicer looking
            if ic[j] in str(ans):
                if j==0:
                    ans = ans.replace(ic[j],soln + "(0)")
                elif j==1:
                    ans = ans.replace(ic[j],soln + "\'(0)")
                elif j==2:
                    ans = ans.replace(ic[j],soln + "\'\'(0)")
                else:
                    ans = ans.replace(ic[j],soln + "^(%s)(0)"%j)
        #print ans
        return ans
