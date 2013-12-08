r"""
Orthogonal Polynomials

This module wraps some of the orthogonal/special functions in the
Maxima package "orthopoly". This package was written by Barton
Willis of the University of Nebraska at Kearney. It is released
under the terms of the General Public License (GPL). Send
Maxima-related bug reports and comments on this module to
willisb@unk.edu. In your report, please include Maxima and specfun
version information.


-  The Chebyshev polynomial of the first kind arises as a solution
   to the differential equation

   .. math::

         (1-x^2)\,y'' - x\,y' + n^2\,y = 0


   and those of the second kind as a solution to

   .. math::

         (1-x^2)\,y'' - 3x\,y' + n(n+2)\,y = 0.


   The Chebyshev polynomials of the first kind are defined by the
   recurrence relation

   .. math::

     T_0(x) = 1 \, T_1(x) = x \, T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x). \,


   The Chebyshev polynomials of the second kind are defined by the
   recurrence relation

   .. math::

     U_0(x) = 1 \, U_1(x) = 2x \, U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x). \,



   For integers `m,n`, they satisfy the orthogonality
   relations

   .. math::

     \int_{-1}^1 T_n(x)T_m(x)\,\frac{dx}{\sqrt{1-x^2}} =\left\{ \begin{matrix} 0 &: n\ne m~~~~~\\ \pi &: n=m=0\\ \pi/2 &: n=m\ne 0 \end{matrix} \right.


   and


   .. math::

     \int_{-1}^1 U_n(x)U_m(x)\sqrt{1-x^2}\,dx =\frac{\pi}{2}\delta_{m,n}.



   They are named after Pafnuty Chebyshev (alternative
   transliterations: Tchebyshef or Tschebyscheff).

-  The Hermite polynomials are defined either by

   .. math::

     H_n(x)=(-1)^n e^{x^2/2}\frac{d^n}{dx^n}e^{-x^2/2}


   (the "probabilists' Hermite polynomials"), or by


   .. math::

     H_n(x)=(-1)^n e^{x^2}\frac{d^n}{dx^n}e^{-x^2}


   (the "physicists' Hermite polynomials"). Sage (via Maxima)
   implements the latter flavor. These satisfy the orthogonality
   relation

   .. math::

     \int_{-\infty}^\infty H_n(x)H_m(x)\,e^{-x^2}\,dx ={n!2^n}{\sqrt{\pi}}\delta_{nm}



   They are named in honor of Charles Hermite.

-  Each *Legendre polynomial* `P_n(x)` is an `n`-th degree polynomial.
   It may be expressed using Rodrigues' formula:

   .. math::

      P_n(x) = (2^n n!)^{-1} {\frac{d^n}{dx^n} } \left[ (x^2 -1)^n \right].

   These are solutions to Legendre's differential equation:

   .. math::

      {\frac{d}{dx}} \left[ (1-x^2) {\frac{d}{dx}} P(x) \right] + n(n+1)P(x) = 0.

   and satisfy the orthogonality relation

   .. math::

      \int_{-1}^{1} P_m(x) P_n(x)\,dx = {\frac{2}{2n + 1}} \delta_{mn}

   The *Legendre function of the second kind* `Q_n(x)` is another
   (linearly independent) solution to the Legendre differential equation.
   It is not an "orthogonal polynomial" however.

   The associated Legendre functions of the first kind
   `P_\ell^m(x)` can be given in terms of the "usual"
   Legendre polynomials by

   .. math::

     \begin{array}{ll} P_\ell^m(x)    &=  (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}P_\ell(x) \\ &=  \frac{(-1)^m}{2^\ell \ell!} (1-x^2)^{m/2}\frac{d^{\ell+m}}{dx^{\ell+m}}(x^2-1)^\ell. \end{array}


   Assuming `0 \le m \le \ell`, they satisfy the orthogonality
   relation:

   .. math::

      \int_{-1}^{1} P_k ^{(m)} P_\ell ^{(m)} dx  = \frac{2 (\ell+m)!}{(2\ell+1)(\ell-m)!}\ \delta _{k,\ell},


   where `\delta _{k,\ell}` is the Kronecker delta.

   The associated Legendre functions of the second kind
   `Q_\ell^m(x)` can be given in terms of the "usual"
   Legendre polynomials by


   .. math::

     Q_\ell^m(x)   =  (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}Q_\ell(x).



   They are named after Adrien-Marie Legendre.

-  Laguerre polynomials may be defined by the Rodrigues formula

   .. math::

      L_n(x)=\frac{e^x}{n!}\frac{d^n}{dx^n}\left(e^{-x} x^n\right).


   They are solutions of Laguerre's equation:


   .. math::

      x\,y'' + (1 - x)\,y' + n\,y = 0\,

   and satisfy the orthogonality relation


   .. math::

      \int_0^\infty L_m(x) L_n(x) e^{-x}\,dx = \delta_{mn}.



   The generalized Laguerre polynomials may be defined by the
   Rodrigues formula:


   .. math::

       L_n^{(\alpha)}(x)   = {\frac{x^{-\alpha} e^x}{n!}}{\frac{d^n}{dx^n}} \left(e^{-x} x^{n+\alpha}\right) .


   (These are also sometimes called the associated Laguerre
   polynomials.) The simple Laguerre polynomials are recovered from
   the generalized polynomials by setting `\alpha =0`.

   They are named after Edmond Laguerre.

-  Jacobi polynomials are a class of orthogonal polynomials. They
   are obtained from hypergeometric series in cases where the series
   is in fact finite:

   .. math::

     P_n^{(\alpha,\beta)}(z) =\frac{(\alpha+1)_n}{n!} \,_2F_1\left(-n,1+\alpha+\beta+n;\alpha+1;\frac{1-z}{2}\right) ,


   where `()_n` is Pochhammer's symbol (for the rising
   factorial), (Abramowitz and Stegun p561.) and thus have the
   explicit expression


   .. math::

     P_n^{(\alpha,\beta)} (z) = \frac{\Gamma (\alpha+n+1)}{n!\Gamma (\alpha+\beta+n+1)} \sum_{m=0}^n {n\choose m} \frac{\Gamma (\alpha + \beta + n + m + 1)}{\Gamma (\alpha + m + 1)} \left(\frac{z-1}{2}\right)^m .



   They are named after Carl Jacobi.

-  Ultraspherical or Gegenbauer polynomials are given in terms of
   the Jacobi polynomials `P_n^{(\alpha,\beta)}(x)` with
   `\alpha=\beta=a-1/2` by


   .. math::

     C_n^{(a)}(x)= \frac{\Gamma(a+1/2)}{\Gamma(2a)}\frac{\Gamma(n+2a)}{\Gamma(n+a+1/2)} P_n^{(a-1/2,a-1/2)}(x).


   They satisfy the orthogonality relation

   .. math::

     \int_{-1}^1(1-x^2)^{a-1/2}C_m^{(a)}(x)C_n^{(a)}(x)\, dx =\delta_{mn}2^{1-2a}\pi \frac{\Gamma(n+2a)}{(n+a)\Gamma^2(a)\Gamma(n+1)} ,


   for `a>-1/2`. They are obtained from hypergeometric series
   in cases where the series is in fact finite:


   .. math::

     C_n^{(a)}(z) =\frac{(2a)^{\underline{n}}}{n!} \,_2F_1\left(-n,2a+n;a+\frac{1}{2};\frac{1-z}{2}\right)


   where `\underline{n}` is the falling factorial. (See
   Abramowitz and Stegun p561)

   They are named for Leopold Gegenbauer (1849-1903).


For completeness, the Pochhammer symbol, introduced by Leo August
Pochhammer, `(x)_n`, is used in the theory of special
functions to represent the "rising factorial" or "upper factorial"

.. math::

         (x)_n=x(x+1)(x+2)\cdots(x+n-1)=\frac{(x+n-1)!}{(x-1)!}.


On the other hand, the "falling factorial" or "lower factorial" is

.. math::

     x^{\underline{n}}=\frac{x!}{(x-n)!} ,


in the notation of Ronald L. Graham, Donald E. Knuth and Oren
Patashnik in their book Concrete Mathematics.

.. note::

   The first call of any of these will usually cost a bit extra
   (it loads "specfun", but I'm not sure if that is the real reason).
   The next call is usually faster but not always.

.. TODO::

    Implement associated Legendre polynomials and Zernike
    polynomials. (Neither is in Maxima.)
    :wikipedia:`Associated_Legendre_polynomials`
    :wikipedia:`Zernike_polynomials`

REFERENCES:

.. [ASHandbook] Abramowitz and Stegun: Handbook of Mathematical Functions,
    http://www.math.sfu.ca/ cbm/aands/

.. :wikipedia:`Chebyshev_polynomials`

.. :wikipedia:`Legendre_polynomials`

.. :wikipedia:`Hermite_polynomials`

.. http://mathworld.wolfram.com/GegenbauerPolynomial.html

.. :wikipedia:`Jacobi_polynomials`

.. :wikipedia:`Laguerre_polynomia`

.. :wikipedia:`Associated_Legendre_polynomials`

.. [EffCheby] Wolfram Koepf: Effcient Computation of Chebyshev Polynomials 
    in Computer Algebra
    Computer Algebra Systems: A Practical Guide. 
    John Wiley, Chichester (1999): 79-99.

AUTHORS:

- David Joyner (2006-06)
- Stefan Reiterer (2010-)
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdj@usna.edu>
#                     2010 Stefan Reiterer <maldun.finsterschreck@gmail.com>
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

import warnings

from sage.misc.sage_eval import sage_eval
from sage.rings.all import ZZ, RR, CC, RIF, CIF
from sage.calculus.calculus import maxima


from sage.symbolic.function import BuiltinFunction, GinacFunction, is_inexact
from sage.symbolic.expression import is_Expression
import sage.functions.special
from sage.functions.special import MaximaFunction, meval
from sage.functions.other import floor, gamma, factorial, abs, binomial
from sage.functions.other import sqrt, conjugate
from sage.functions.trig import sin, cos
from sage.functions.log import ln
import sage.symbolic.expression as expression
from sage.structure.parent import Parent
from sage.structure.coerce import parent

_done = False
def _init():
    """
    Internal function which checks if Maxima has loaded the
    "orthopoly" package.  All functions using this in this
    file should call this function first.

    TEST:

    The global starts ``False``::

        sage: sage.functions.orthogonal_polys._done
        False

    Then after using one of these functions, it changes::

        sage: from sage.functions.orthogonal_polys import legendre_P 
        sage: legendre_P(2,x)
        3/2*(x - 1)^2 + 3*x - 2
        sage: sage.functions.orthogonal_polys._done
        True


    Note that because here we use a Pynac variable ``x``,
    the representation of the function is different from
    its actual doctest, where a polynomial indeterminate
    ``x`` is used.
    """
    global _done
    if _done:
        return
    maxima.eval('load("orthopoly");')
    # TODO -- make it possible to use the intervals returned
    # instead of just discarding this info!
    maxima.eval('orthopoly_returns_intervals:false;')
    _done = True



class OrthogonalPolynomial(BuiltinFunction):
    """
    Base class for orthogonal polynomials.

    This class is an abstract base class for all orthogonal polynomials since
    they share similar properties. The evaluation as a polynomial
    is either done via maxima, or with pynac.

    Convention: The first argument is always the order of the polynomial, 
    the last one is always the value `x` where the polynomial is evaluated.
    """
    def __init__(self, name, nargs=2, latex_name=None, conversions={}):
        """
        :class:`OrthogonalPolynomial` class needs the same input parameter as
        it's parent class.
        
        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import OrthogonalPolynomial
            sage: new = OrthogonalPolynomial('testo_P')
            sage: new
            testo_P
        """
        try:
            self._maxima_name = conversions['maxima'] 
        except KeyError:
            self._maxima_name = None

        super(OrthogonalPolynomial,self).__init__(name=name, nargs=nargs, 
                                 latex_name=latex_name, conversions=conversions)

    def _maxima_init_evaled_(self, *args):
        r"""
        Return a string which represents this function evaluated at
        ``*args`` in Maxima.

        In fact these are thought to be the old wrappers for the orthogonal
        polynomials. They are used when the other evaluation methods fail,
        or are not fast enough. It appears that direct computation
        with pynac is in most cases faster than maxima. Maxima comes into
        play when all other methods fail.

        EXAMPLES::

            sage: chebyshev_T(3,x)
            4*x^3 - 3*x
        """
        return None

    def _apply_formula_(self, *args):
        """
        Method which uses the three term recursion of the polynomial, 
        or explicit formulas instead of maxima to evaluate the polynomial 
        efficiently, if the `x` argument is not a symbolic expression. 

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import OrthogonalPolynomial
            sage: new = OrthogonalPolynomial('testo_P')
            sage: new._apply_formula_(1,2.0)
            Traceback (most recent call last):
            ...
            NotImplementedError: no recursive calculation of values implemented
        """
        raise NotImplementedError("no recursive calculation of values implemented")

    def _eval_special_values_(self,*args):
        """
        Evaluate the polynomial explicitly for special values.

        EXAMPLES::
            
            sage: var('n')
            n
            sage: chebyshev_T(n,-1)
            (-1)^n
        """
        raise ValueError("no special values known")

    def _eval_(self, *args):
        """
        The _eval_ method decides which evaluation suits best
        for the given input, and returns a proper value.
        
        EXAMPLES::

            sage: chebyshev_T(5,x)
            16*x^5 - 20*x^3 + 5*x  
            sage: var('n')
            n
            sage: chebyshev_T(n,-1)
            (-1)^n
            sage: chebyshev_T(-7,x)
            64*x^7 - 112*x^5 + 56*x^3 - 7*x
            sage: chebyshev_T(3/2,x)
            chebyshev_T(3/2, x)
            sage: x = PolynomialRing(QQ, 'x').gen()
            sage: chebyshev_T(2,x)
            2*x^2 - 1
            sage: chebyshev_U(2,x)
            4*x^2 - 1
            sage: parent(chebyshev_T(4, RIF(5)))
            Real Interval Field with 53 bits of precision
            sage: RR2 = RealField(5)
            sage: chebyshev_T(100000,RR2(2))
            8.9e57180
            sage: chebyshev_T(5,Qp(3)(2))
            2 + 3^2 + 3^3 + 3^4 + 3^5 + O(3^20)
            sage: chebyshev_T(100001/2, 2)
            doctest:500: RuntimeWarning: Warning: mpmath returns NoConvergence exception! Use other method instead.
            chebyshev_T(100001/2, 2)
            sage: chebyshev_U._eval_(1.5, Mod(8,9)) is None
            True
        """
        if not is_Expression(args[0]):
        # If x is no expression and is inexact or n is not an integer -> make numerical evaluation
            if (not is_Expression(args[-1])) and (is_inexact(args[-1]) or not args[0] in ZZ):
                try:
                    import sage.libs.mpmath.all as mpmath
                    return self._evalf_(*args)
                except AttributeError:
                    pass
                except mpmath.NoConvergence:
                    warnings.warn("Warning: mpmath returns NoConvergence exception! Use other method instead.",
                                  RuntimeWarning)
                except ValueError:
                    pass

            # n is not an integer and x is an expression -> return symbolic expression.
            if not args[0] in ZZ:
                if is_Expression(args[-1]):
                    return None

        # Check for known identities 
        try:
            return self._eval_special_values_(*args)
        except ValueError:
            pass

        #if negative indices are not specified 
        #in _eval_special_values only return symbolic
        #value
        if args[0] < 0 and args[0] in ZZ:
                return None

        if args[0] in ZZ:
            try: 
                return self._apply_formula_(*args)
            except NotImplementedError:
                pass

        if self._maxima_name is None:
            return None

        if args[0] in ZZ: # use maxima as last resort
            return self._old_maxima_(*args)
        else:
            return None
        
    def __call__(self,*args,**kwds):
        """ 
        This overides the call method from SageObject to avoid problems with coercions,
        since the _eval_ method is able to handle more data types than symbolic functions
        would normally allow.
        Thus we have the distinction between algebraic objects (if n is an integer),
        and else as symbolic function.
        
        EXAMPLES:: 
        
            sage: K.<a> = NumberField(x^3-x-1) 
            sage: chebyshev_T(5, a) 
            16*a^2 + a - 4 
            sage: chebyshev_T(5,MatrixSpace(ZZ, 2)([1, 2, -4, 7]))
            [-40799  44162]
            [-88324  91687]
            sage: R.<x> = QQ[]
            sage: parent(chebyshev_T(5, x))
            Univariate Polynomial Ring in x over Rational Field
        """ 
        if 'hold' not in kwds:
            kwds['hold'] = False
        if 'coerce' not in kwds:
            kwds['coerce']=True
        
        if args[0] in ZZ and kwds['hold'] is False: #check if n is in ZZ->consider polynomial as algebraic structure
            return self._eval_(*args) # Let eval methode decide which is best
        else: # Consider OrthogonalPolynomial as symbol
            return super(OrthogonalPolynomial,self).__call__(*args,**kwds)

    def _old_maxima_(self,*args):
        """
        Method which holds the old maxima wrappers as last alternative.
        It returns None per default, and it only needs to be implemented,
        if it is necessary.

        EXAMPLES::
        
            sage: chebyshev_T._old_maxima_(-7,x) is None
            True
        """
        None

class Func_chebyshev_T(OrthogonalPolynomial): 
    """
    Chebyshev polynomials of the first kind.

    REFERENCE:

    - [ASHandbook]_ 22.5.31 page 778 and 6.1.22 page 256.

    EXAMPLES::

       sage: chebyshev_T(3,x)
       4*x^3 - 3*x
       sage: chebyshev_T(5,x)
       16*x^5 - 20*x^3 + 5*x
       sage: var('k')
       k
       sage: test = chebyshev_T(k,x)
       sage: test
       chebyshev_T(k, x)
    """
    def __init__(self):
        """
        Init method for the chebyshev polynomials of the first kind.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import Func_chebyshev_T
            sage: chebyshev_T2 = Func_chebyshev_T()
            sage: chebyshev_T2(1,x)
            x
        """
        super(Func_chebyshev_T,self).__init__("chebyshev_T", nargs=2,
                                      conversions=dict(maxima='chebyshev_t',
                                                       mathematica='ChebyshevT'))
    
    def _eval_special_values_(self,*args):
        """
        Values known for special values of x.
        For details see [ASHandbook]_ 22.4 (p. 777)

        EXAMPLES:

            sage: var('n')
            n
            sage: chebyshev_T(n,1)
            1
            sage: chebyshev_T(n,-1)
            (-1)^n
            sage: chebyshev_T(-7, x) - chebyshev_T(7,x)
            0
            sage: chebyshev_T._eval_special_values_(3/2,x)
            Traceback (most recent call last):
            ...
            ValueError: No special values for non integral indices!
            sage: chebyshev_T._eval_special_values_(n, 0.1)
            Traceback (most recent call last):
            ...
            ValueError: Value not found!
            sage: chebyshev_T._eval_special_values_(26, Mod(9,9))
            Traceback (most recent call last):
            ...
            ValueError: Value not found!
        """
        if (not is_Expression(args[0])) and (not args[0] in ZZ):
            raise ValueError("No special values for non integral indices!")
        
        if args[-1] == 1:
            return args[-1]
        
        if args[-1] == -1:
            return args[-1]**args[0]

        if (args[-1] == 0 and args[-1] in CC):
            return (1+(-1)**args[0])*(-1)**(args[0]/2)/2

        if args[0] < 0 and args[0] in ZZ:
            return self._eval_(-args[0],args[-1])

        raise ValueError("Value not found!")
    
    def _evalf_(self, *args,**kwds):
        """
        Evaluates :class:`chebyshev_T` numerically with mpmath.
        If the index is an integer we use the recursive formula since
        it is faster.

        EXAMPLES::

            sage: chebyshev_T(10,3).n(75)
            2.261953700000000000000e7
            sage: chebyshev_T(10,I).n()
            -3363.00000000000
            sage: chebyshev_T(5,0.3).n()
            0.998880000000000
            sage: chebyshev_T(1/2, 0)
            0.707106781186548
            sage: chebyshev_T._evalf_(1.5, Mod(8,9))
            Traceback (most recent call last):
            ...
            ValueError: No compatible type!

        """
        if args[0] in ZZ and args[0] >= 0:
            return self._cheb_recur_(*args)[0]
        
        try:
            real_parent = kwds['parent']
        except KeyError:
            real_parent = parent(args[-1])

        x_set = False
        if hasattr(real_parent,"precision"): # Check if we have a data type with precision 
            x = args[-1]
            step_parent = real_parent
            x_set = True
        else:
            if args[-1] in RR:
                x = RR(args[-1])
                step_parent = RR
                x_set = True
            elif args[-1] in CC:
                x = CC(args[-1])
                step_parent = CC
                x_set = True

        if not x_set:
            raise ValueError("No compatible type!")
                
        from sage.libs.mpmath.all import call as mpcall
        from sage.libs.mpmath.all import chebyt as mpchebyt

        return mpcall(mpchebyt,args[0],x,parent=step_parent)

    def _maxima_init_evaled_(self, *args):
        """
        Evaluate the Chebyshev polynomial ``self`` with maxima.

        EXAMPLES::

            sage: chebyshev_T._maxima_init_evaled_(1,x)
            'x'
            sage: var('n')
            n
            sage: maxima(chebyshev_T(n,x))
            chebyshev_t(n,x)

        """
        n = args[0]
        x = args[1]
        return maxima.eval('chebyshev_t({0},{1})'.format(n,x))

        
    def _apply_formula_(self,*args):
        """
        Applies explicit formulas for :class:`chebyshev_T`. 
        This is much faster for numerical evaluation than maxima!
        See [ASHandbook]_ 227 (p. 782) for details for the recurions.
        See also [EffCheby]_ for fast evaluation techniques.

        EXAMPLES::

            sage: chebyshev_T._apply_formula_(2,0.1) == chebyshev_T._evalf_(2,0.1)
            True
            sage: chebyshev_T(51,x)
            2*(2*(2*(2*(2*(2*x^2 - 1)^2 - 1)*(2*(2*x^2 - 1)*x - x) - x)*(2*(2*(2*x^2 - 1)*x - x)^2 - 1) - x)^2 - 1)*(2*(2*(2*(2*(2*x^2 - 1)^2 - 1)*(2*(2*x^2 - 1)*x - x) - x)*(2*(2*(2*x^2 - 1)*x - x)^2 - 1) - x)*(2*(2*(2*(2*x^2 - 1)*x - x)^2 - 1)^2 - 1) - x) - x
            sage: chebyshev_T._apply_formula_(10,x)
            512*x^10 - 1280*x^8 + 1120*x^6 - 400*x^4 + 50*x^2 - 1

        """
        k = args[0]
        x = args[1]

        if k == 0:
            return 1
        if k == 1:
            return x
        
        help1 = 1
        help2 = x
        if is_Expression(x) and k <= 25:
            # Recursion gives more compact representations for large k
            help3 = 0
            for j in xrange(0,floor(k/2)+1):
                f = factorial(k-j-1) / factorial(j) / factorial(k-2*j)
                help3 = help3 + (-1)**j * (2*x)**(k-2*j) * f
            help3 = help3 * k / 2
            return help3
        else:
            return self._cheb_recur_(k,x)[0]

    def _cheb_recur_(self,n, x, both=False): 
        """ 
        Generalized recursion formula for Chebyshev polynomials. 
        Implementation suggested by Frederik Johansson. 
        returns (T(n,x), T(n-1,x)), or (T(n,x), _) if both=False 
        
        EXAMPLES:: 
        
            sage: chebyshev_T._cheb_recur_(5,x) 
            (2*(2*(2*x^2 - 1)*x - x)*(2*x^2 - 1) - x, False) 
        """ 
        if n == 0: 
            return 1, x 
        if n == 1: 
            return x, 1 
        a, b = self._cheb_recur_((n+1)//2, x, both or n % 2) 
        if n % 2 == 0: 
            return 2*a**2 - 1, both and 2*a*b - x 
        else: 
            return 2*a*b - x, both and 2*b**2 - 1 


    def _eval_numpy_(self, *args):
        """
        Evaluate ``self`` using numpy.

        EXAMPLES::

            sage: import numpy
            sage: z = numpy.array([1,2])
            sage: z2 = numpy.array([[1,2],[1,2]])
            sage: z3 = numpy.array([1,2,3.])
            sage: chebyshev_T(1,z)
            array([1, 2])
            sage: chebyshev_T(1,z2)
            array([[1, 2],
                   [1, 2]])
            sage: chebyshev_T(1,z3)
            array([ 1.,  2.,  3.])
            sage: chebyshev_T(z,0.1)
            array([ 0.1 , -0.98])
        """
        from scipy.special import eval_chebyt
        return eval_chebyt(args[0],args[-1])

    def _derivative_(self, *args, **kwds):
        """
        Return the derivative of :class:`chebyshev_T` in form of the Chebyshev
        polynomial of the second kind :class:`chebyshev_U`.

        EXAMPLES::

            sage: var('k')
            k
            sage: derivative(chebyshev_T(k,x),x)
            k*chebyshev_U(k - 1, x)
            sage: derivative(chebyshev_T(3,x),x)
            12*x^2 - 3
            sage: derivative(chebyshev_T(k,x),k)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative w.r.t. to the index is not supported yet
        """
        diff_param = kwds['diff_param']
        if diff_param == 0: 
            raise NotImplementedError("derivative w.r.t. to the index is not supported yet")

        return args[0]*chebyshev_U(args[0]-1,args[1])
    
chebyshev_T = Func_chebyshev_T()

class Func_chebyshev_U(OrthogonalPolynomial):
    """
    Class for the Chebyshev polynomial of the second kind.

    REFERENCE:

    - [ASHandbook]_ 22.8.3 page 783 and 6.1.22 page 256.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: chebyshev_U(2,x)
        4*x^2 - 1
        sage: chebyshev_U(3,x)
        8*x^3 - 4*x
    """
    def __init__(self):
        """
        Init method for the chebyshev polynomials of the second kind.

        EXAMPLES::
        
            sage: from sage.functions.orthogonal_polys import Func_chebyshev_U
            sage: chebyshev_U2 = Func_chebyshev_U()
            sage: chebyshev_U2(1,x)
            2*x
        """
        OrthogonalPolynomial.__init__(self, "chebyshev_U", nargs=2,
                                      conversions=dict(maxima='chebyshev_u',
                                                       mathematica='ChebyshevU'))
        
    def _apply_formula_(self,*args):
        """
        Applies explicit formulas for :class:`chebyshev_U`.
        This is much faster for numerical evaluation than maxima.
        See [ASHandbook]_ 227 (p. 782) for details on the recurions.
        See also [EffCheby]_ for the recursion formulas.

        EXAMPLES::

            sage: chebyshev_U._apply_formula_(2,0.1) == chebyshev_U._evalf_(2,0.1)
            True
        """
        k = args[0]
        x = args[1]
           
        if k == 0:
            return 1
        if k == 1:
            return 2*x

        help1 = 1
        help2 = 2*x
        if is_Expression(x) and k <= 25:
            # Recursion gives more compact representations for large k
            help3 = 0
            for j in xrange(0,floor(k/2)+1):
                f = factorial(k-j) / factorial(j) / factorial(k-2*j) # Change to a binomial?
                help3 = help3 + (-1)**j * (2*x)**(k-2*j) * f
            return help3
            
        else:
            return self._cheb_recur_(k,x)[0]

 
    def _cheb_recur_(self,n, x, both=False): 
        """ 
        Generalized recursion formula for Chebyshev polynomials.
        Implementation suggested by Frederik Johansson. 
        returns (U(n,x), U(n-1,x)), or (U(n,x), _) if both=False 
        
        EXAMPLES:: 
        
            sage: chebyshev_U._cheb_recur_(3,x) 
            (4*(2*x^2 - 1)*x, False) 
            sage: chebyshev_U._cheb_recur_(5,x)[0]
            -2*((2*x + 1)*(2*x - 1)*x - 4*(2*x^2 - 1)*x)*(2*x + 1)*(2*x - 1)
            sage: abs(pari('polchebyshev(5, 2, 0.1)') - chebyshev_U(5,0.1)) < 1e-10
            True
        """ 
        
        if n == 0: 
            return 1, both and 2*x 
        if n == 1: 
            return 2*x, both and 4*x**2-1 
            
        a, b = self._cheb_recur_((n-1)//2, x, True) 
        if n % 2 == 0: 
            return (b+a)*(b-a), both and 2*b*(x*b-a) 
        else: 
            return 2*a*(b-x*a), both and (b+a)*(b-a) 

    def _maxima_init_evaled_(self, *args):
        """
        Uses maxima to evaluate ``self``.

        EXAMPLES::

            sage: maxima(chebyshev_U(5,x))
            32*x^5-32*x^3+6*x
            sage: var('n')
            n
            sage: maxima(chebyshev_U(n,x))
            chebyshev_u(n,x)
            sage: maxima(chebyshev_U(2,x))
            4*x^2-1
        """
        n = args[0]
        x = args[1]
        return maxima.eval('chebyshev_u({0},{1})'.format(n,x))

    def _evalf_(self, *args,**kwds):
        """
        Evaluate :class:`chebyshev_U` numerically with mpmath.
        If index is an integer use recursive formula since it is faster,
        for chebyshev polynomials.

        EXAMPLES::

            sage: chebyshev_U(5,-4+3.*I)
            98280.0000000000 - 11310.0000000000*I
            sage: chebyshev_U(10,3).n(75)
            4.661117900000000000000e7
            sage: chebyshev_U._evalf_(1.5, Mod(8,9))
            Traceback (most recent call last):
            ...
            ValueError: No compatible type!
        """
        if args[0] in ZZ and args[0] >= 0:
            return self._cheb_recur_(*args)[0]
        try:
            real_parent = kwds['parent']
        except KeyError:
            real_parent = parent(args[-1])

        x_set = False
        if hasattr(real_parent,"precision"): # Check if we have a data type with precision 
            x = args[-1]
            step_parent = real_parent
            x_set = True
        else:
            if args[-1] in RR:
                x = RR(args[-1])
                step_parent = RR
                x_set = True
            elif args[-1] in CC:
                x = CC(args[-1])
                step_parent = CC
                x_set = True

        if not x_set:
            raise ValueError("No compatible type!")

        from sage.libs.mpmath.all import call as mpcall
        from sage.libs.mpmath.all import chebyu as mpchebyu

        return mpcall(mpchebyu,args[0],args[-1],parent = step_parent)

    def _eval_special_values_(self,*args):
        """
        Special values that known. [ASHandbook]_ 22.4 (p.777).

        EXAMPLES::
        
            sage: var('n')
            n
            sage: chebyshev_U(n,1)
            n + 1
            sage: chebyshev_U(n,-1)
            (-1)^n*(n + 1)
            sage: chebyshev_U._eval_special_values_(26, Mod(0,9))
            Traceback (most recent call last):
            ...
            ValueError: Value not found!
            sage: parent(chebyshev_U(3, Mod(8,9)))
            Ring of integers modulo 9
            sage: parent(chebyshev_U(3, Mod(1,9)))
            Ring of integers modulo 9
            sage: chebyshev_U(n, 0)
            1/2*(-1)^(1/2*n)*((-1)^n + 1)
            sage: chebyshev_U(-3,x) + chebyshev_U(1,x)
            0
            sage: chebyshev_U._eval_special_values_(1.5, Mod(8,9))
            Traceback (most recent call last):
            ...
            ValueError: No special values for non integral indices!
            sage: chebyshev_U(-1,Mod(5,8))
            0
            sage: parent(chebyshev_U(-1,Mod(5,8)))
            Ring of integers modulo 8
        """
        if (not is_Expression(args[0])) and (not args[0] in ZZ):
            raise ValueError("No special values for non integral indices!")

        if args[0] == -1:
          return args[-1]*0
            
        if args[-1] == 1:
            return args[-1]*(args[0]+1)
        
        if args[-1] == -1:
            return args[-1]**args[0]*(args[0]+1)

        if (args[-1] == 0 and args[-1] in CC):
            return (1+(-1)**args[0])*(-1)**(args[0]/2)/2

        if args[0] < 0 and args[0] in ZZ:
            return -self._eval_(-args[0]-2,args[-1])

        raise ValueError("Value not found!")

    def _eval_numpy_(self, *args):
        """
        Evaluate ``self`` using numpy.

        EXAMPLES::

            sage: import numpy
            sage: z = numpy.array([1,2])
            sage: z2 = numpy.array([[1,2],[1,2]])
            sage: z3 = numpy.array([1,2,3.])
            sage: chebyshev_U(1,z)
            array([2, 4])
            sage: chebyshev_U(1,z2)
            array([[2, 4],
                   [2, 4]])
            sage: chebyshev_U(1,z3)
            array([ 2.,  4.,  6.])
            sage: chebyshev_U(z,0.1)
            array([ 0.2 , -0.96])
        """
        from scipy.special import eval_chebyu
        return eval_chebyu(args[0],args[1])


    def _derivative_(self, *args, **kwds):
        """
        Return the derivative of :class:`chebyshev_U` in form of the Chebyshev
        polynomials of the first and second kind.
        
        EXAMPLES::

            sage: var('k')
            k
            sage: derivative(chebyshev_U(k,x),x)
            ((k + 1)*chebyshev_T(k + 1, x) - x*chebyshev_U(k, x))/(x^2 - 1)
            sage: derivative(chebyshev_U(3,x),x)
            24*x^2 - 4
            sage: derivative(chebyshev_U(k,x),k)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative w.r.t. to the index is not supported yet
        """
        diff_param = kwds['diff_param']
        if diff_param == 0: 
            raise NotImplementedError("derivative w.r.t. to the index is not supported yet")

        return ((args[0]+1)*chebyshev_T(args[0]+1,args[1])-args[1]*
                chebyshev_U(args[0],args[1]))/(args[1]**2-1)

chebyshev_U = Func_chebyshev_U()


def gen_laguerre(n,a,x):
    """
    Returns the generalized Laguerre polynomial for integers `n > -1`.
    Typically, `a = 1/2` or `a = -1/2`.

    REFERENCES:

    - Table on page 789 in [ASHandbook]_.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: gen_laguerre(2,1,x)
        1/2*x^2 - 3*x + 3
        sage: gen_laguerre(2,1/2,x)
        1/2*x^2 - 5/2*x + 15/8
        sage: gen_laguerre(2,-1/2,x)
        1/2*x^2 - 3/2*x + 3/8
        sage: gen_laguerre(2,0,x)
        1/2*x^2 - 2*x + 1
        sage: gen_laguerre(3,0,x)
        -1/6*x^3 + 3/2*x^2 - 3*x + 1
    """
    from sage.functions.all import sqrt
    _init()
    return sage_eval(maxima.eval('gen_laguerre(%s,%s,x)'%(ZZ(n),a)), locals={'x':x})

def gen_legendre_P(n,m,x):
    r"""
    Returns the generalized (or associated) Legendre function of the
    first kind for integers `n > -1, m > -1`.

    The awkward code for when m is odd and 1 results from the fact that
    Maxima is happy with, for example, `(1 - t^2)^3/2`, but
    Sage is not. For these cases the function is computed from the
    (m-1)-case using one of the recursions satisfied by the Legendre
    functions.

    REFERENCE:

    - Gradshteyn and Ryzhik 8.706 page 1000.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: gen_legendre_P(2, 0, t)
        3/2*t^2 - 1/2
        sage: gen_legendre_P(2, 0, t) == legendre_P(2, t)
        True
        sage: gen_legendre_P(3, 1, t)
        -3/2*(5*t^2 - 1)*sqrt(-t^2 + 1)
        sage: gen_legendre_P(4, 3, t)
        105*(t^2 - 1)*sqrt(-t^2 + 1)*t
        sage: gen_legendre_P(7, 3, I).expand()
        -16695*sqrt(2)
        sage: gen_legendre_P(4, 1, 2.5)
        -583.562373654533*I
    """
    from sage.functions.all import sqrt
    _init()
    if m.mod(2).is_zero() or m.is_one():
        return sage_eval(maxima.eval('assoc_legendre_p(%s,%s,x)'%(ZZ(n),ZZ(m))), locals={'x':x})
    else:
        return sqrt(1-x**2)*(((n-m+1)*x*gen_legendre_P(n,m-1,x)-(n+m-1)*gen_legendre_P(n-1,m-1,x))/(1-x**2))

def gen_legendre_Q(n,m,x):
    """
    Returns the generalized (or associated) Legendre function of the
    second kind for integers `n>-1`, `m>-1`.

    Maxima restricts m = n. Hence the cases m n are computed using the
    same recursion used for gen_legendre_P(n,m,x) when m is odd and
    1.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: gen_legendre_Q(2,0,t)
        3/4*t^2*log(-(t + 1)/(t - 1)) - 3/2*t - 1/4*log(-(t + 1)/(t - 1))
        sage: gen_legendre_Q(2,0,t) - legendre_Q(2, t)
        0
        sage: gen_legendre_Q(3,1,0.5)
        2.49185259170895
        sage: gen_legendre_Q(0, 1, x)
        -1/sqrt(-x^2 + 1)
        sage: gen_legendre_Q(2, 4, x).factor()
        48*x/((x + 1)^2*(x - 1)^2)
    """
    from sage.functions.all import sqrt
    if m <= n:
        _init()
        return sage_eval(maxima.eval('assoc_legendre_q(%s,%s,x)'%(ZZ(n),ZZ(m))), locals={'x':x})
    if m == n + 1 or n == 0:
        if m.mod(2).is_zero():
            denom = (1 - x**2)**(m/2)
        else:
            denom = sqrt(1 - x**2)*(1 - x**2)**((m-1)/2)
        if m == n + 1:
            return (-1)**m*(m-1).factorial()*2**n/denom
        else:
            return (-1)**m*(m-1).factorial()*((x+1)**m - (x-1)**m)/(2*denom)
    else:
        return ((n-m+1)*x*gen_legendre_Q(n,m-1,x)-(n+m-1)*gen_legendre_Q(n-1,m-1,x))/sqrt(1-x**2)

def hermite(n,x):
    """
    Returns the Hermite polynomial for integers `n > -1`.

    REFERENCE:

    - [ASHandbook]_ 22.5.40 and 22.5.41, page 779.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: hermite(2,x)
        4*x^2 - 2
        sage: hermite(3,x)
        8*x^3 - 12*x
        sage: hermite(3,2)
        40
        sage: S.<y> = PolynomialRing(RR)
        sage: hermite(3,y)
        8.00000000000000*y^3 - 12.0000000000000*y
        sage: R.<x,y> = QQ[]
        sage: hermite(3,y^2)
        8*y^6 - 12*y^2
        sage: w = var('w')
        sage: hermite(3,2*w)
        8*(8*w^2 - 3)*w
    """
    _init()
    return sage_eval(maxima.eval('hermite(%s,x)'%ZZ(n)), locals={'x':x})

def jacobi_P(n,a,b,x):
    r"""
    Returns the Jacobi polynomial `P_n^{(a,b)}(x)` for
    integers `n > -1` and a and b symbolic or `a > -1`
    and `b > -1`. The Jacobi polynomials are actually defined
    for all a and b. However, the Jacobi polynomial weight
    `(1-x)^a(1+x)^b` isn't integrable for `a \leq -1`
    or `b \leq -1`.

    REFERENCE:

    - Table on page 789 in [ASHandbook]_.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: jacobi_P(2,0,0,x)
        3/2*x^2 - 1/2
        sage: jacobi_P(2,1,2,1.2)        # random output of low order bits
        5.009999999999998
    """
    _init()
    return sage_eval(maxima.eval('jacobi_p(%s,%s,%s,x)'%(ZZ(n),a,b)), locals={'x':x})

def laguerre(n,x):
    """
    Return the Laguerre polynomial for integers `n > -1`.

    REFERENCE:

    - [ASHandbook]_ 22.5.16, page 778 and page 789.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: laguerre(2,x)
        1/2*x^2 - 2*x + 1
        sage: laguerre(3,x)
        -1/6*x^3 + 3/2*x^2 - 3*x + 1
        sage: laguerre(2,2)
        -1
    """
    _init()
    return sage_eval(maxima.eval('laguerre(%s,x)'%ZZ(n)), locals={'x':x})

def legendre_P(n,x):
    """
    Returns the Legendre polynomial of the first kind for integers
    `n > -1`.

    REFERENCE:

    - [ASHandbook]_ 22.5.35 page 779.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: legendre_P(2,t)
        3/2*t^2 - 1/2
        sage: legendre_P(3, 1.1)
        1.67750000000000
        sage: legendre_P(3, 1 + I)
        7/2*I - 13/2
        sage: legendre_P(3, MatrixSpace(ZZ, 2)([1, 2, -4, 7]))
        [-179  242]
        [-484  547]
        sage: legendre_P(3, GF(11)(5))
        8
    """
    _init()
    return sage_eval(maxima.eval('legendre_p(%s,x)'%ZZ(n)), locals={'x':x})

def legendre_Q(n,x):
    """
    Returns the Legendre function of the second kind for integers
    `n > -1`.

    Computed using Maxima.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: legendre_Q(2, t)
        3/4*t^2*log(-(t + 1)/(t - 1)) - 3/2*t - 1/4*log(-(t + 1)/(t - 1))
        sage: legendre_Q(3, 0.5)
        -0.198654771479482
        sage: legendre_Q(4, 2)
        443/16*I*pi + 443/16*log(3) - 365/12
        sage: legendre_Q(4, 2.0)
        0.00116107583162324 + 86.9828465962674*I
    """
    _init()
    return sage_eval(maxima.eval('legendre_q(%s,x)'%ZZ(n)), locals={'x':x})

def ultraspherical(n,a,x):
    """
    Returns the ultraspherical (or Gegenbauer) polynomial for integers
    `n > -1`.

    Computed using Maxima.

    REFERENCE:

    - [ASHandbook]_ 22.5.27

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: ultraspherical(2,3/2,x)
        15/2*x^2 - 3/2
        sage: ultraspherical(2,1/2,x)
        3/2*x^2 - 1/2
        sage: ultraspherical(1,1,x)
        2*x
        sage: t = PolynomialRing(RationalField(),"t").gen()
        sage: gegenbauer(3,2,t)
        32*t^3 - 12*t
    """
    _init()
    return sage_eval(maxima.eval('ultraspherical(%s,%s,x)'%(ZZ(n),a)), locals={'x':x})

gegenbauer = ultraspherical
