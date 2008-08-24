###############################################################################
#   SAGE: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

"""
EXAMPLES:
We mix Singular variables with symbolic variables:
    sage: R.<u,v> = QQ[]
    sage: var('a,b,c', ns=1)
    (a, b, c)
    sage: expand((u + v + a + b + c)^2)
    2*a*b + 2*a*c + 2*b*c + a^2 + b^2 + c^2 + (2*u + 2*v)*a + (2*u + 2*v)*b + (2*u + 2*v)*c + u^2 + 2*u*v + v^2
"""

include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"

import ring
import sage.rings.integer

from sage.structure.element cimport ModuleElement, RingElement, Element


from sage.rings.rational import Rational  # Used for sqrt.

cdef class Expression(CommutativeRingElement):
    def __dealloc__(self):
        """
        Delete memory occupied by this expression.
        """
        GEx_destruct(&self._gobj)

    def _repr_(self):
        """
        Return string representation of this symbolic expression.

        EXAMPLES:
            sage: var("x y", ns=1)
            (x, y)
            sage: (x+y)._repr_()
            'x + y'
        """
        return GEx_to_str(&self._gobj)

    def __float__(self):
        """
        Return float conversion of self, assuming self is constant.
        Otherwise, raise a TypeError.

        OUTPUT:
            float -- double precision evaluation of self

        EXAMPLES:
            sage: x = var('x', ns=1); SR = x.parent()
            sage: float(SR(12))
            12.0
            sage: float(SR(2/3))
            0.66666666666666663
            sage: float(sqrt(SR(2)))
            1.4142135623730951
            sage: float(x^2 + 1)
            Traceback (most recent call last):
            ...
            TypeError: float() argument must be a string or a number
            sage: float(SR(RIF(2)))
            Traceback (most recent call last):
            ...
            TypeError: float() argument must be a string or a number
        """
        cdef bint success
        cdef double ans = GEx_to_double(self._gobj, &success)
        if not success:
            raise TypeError, "float() argument must be a string or a number"
        return ans

    def __hash__(self):
        """
        Return hash of this expression.

        EXAMPLES:
            sage: x, y = var("x y", ns=1); S = x.parent()
            sage: hash(x)
            2013265920

        The hash of an object in Python or its coerced version into
        the symbolic ring is the same.
            sage: hash(S(19/23))
            4
            sage: hash(19/23)
            4
            sage: hash(x+y)   # random -- the hash for some expression is unfortunately random
            1631713410
            sage: d = {x+y: 5}
            sage: d
            {x + y: 5}
        """
        return self._gobj.gethash()

    def __richcmp__(left, right, int op):
        """
        Create a formal symbolic inequality or equality.

        EXAMPLES:
            sage: var('x, y', ns=1)
            (x, y)
            sage: x + 2/3 < y^2
            x + 2/3 < y^2
            sage: x^3 -y <= y + x
            x^3 - y <= x + y
            sage: x^3 -y == y + x
            x^3 - y == x + y
            sage: x^3 - y^10 >= y + x^10
            x^3 - y^10 >= x^10 + y
            sage: x^2 > x
            x^2 > x
        """
        cdef Expression l, r
        # We coerce left and right to be Expressions, but
        # we do not know which is already an Expression.
        # We know at least one is.
        try:
            l = left
        except TypeError:
            # if l = left failed, we can definitely use
            # r=right to coerce in l.
            r = right
            l = r.coerce_in(left)
        else:
            # if l = left succeeded, we can certainly
            # use l to coerce in right.
            r = l.coerce_in(right)

        cdef GEx e
        if op == Py_LT:
            e = g_lt(l._gobj, r._gobj)
        elif op == Py_EQ:
            e = g_eq(l._gobj, r._gobj)
        elif op == Py_GT:
            e = g_gt(l._gobj, r._gobj)
        elif op == Py_LE:
            e = g_le(l._gobj, r._gobj)
        elif op == Py_NE:
            e = g_ne(l._gobj, r._gobj)
        elif op == Py_GE:
            e = g_ge(l._gobj, r._gobj)
        return new_Expression_from_GEx(e)

    def __nonzero__(self):
        """
        Return True if self is nonzero.

        EXAMPLES:
            sage: sage.symbolic.ring.NSR(0).__nonzero__()
            False
            sage: sage.symbolic.ring.NSR(1).__nonzero__()
            True
        """
        # TODO: Problem -- if self is a symbolic equality then
        # this is_zero isn't the right thing at all:
        #  sage: bool(x == x+1)
        #  True  # BAD
        # Solution is to probably look something up in ginac manual.
        return not self._gobj.is_zero()

    cdef Expression coerce_in(self, z):
        """
        Quickly coerce z to be an Expression.
        """
        cdef Expression w
        try:
            w = z
            return w
        except TypeError:
            return self._parent._coerce_c(z)

    cdef ModuleElement _add_c_impl(left, ModuleElement right):
        """
        Add left and right.

        EXAMPLES;
            sage.: var("x y", ns=1)
            (x, y)
            sage.: x + y + y + x
            2*x+2*y
        """
        return new_Expression_from_GEx(gadd(left._gobj, (<Expression>right)._gobj))

    cdef ModuleElement _sub_c_impl(left, ModuleElement right):
        """
            sage.: var("x y", ns=1)
            (x, y)
            sage.: x - x
            x-y
        """
        return new_Expression_from_GEx(gsub(left._gobj, (<Expression>right)._gobj))

    cdef RingElement _mul_c_impl(left, RingElement right):
        """
        Multiply left and right.

        EXAMPLES:
            sage: var("x y", ns=1)
            (x, y)
            sage: x*y*y
            x*y^2
        """
        return new_Expression_from_GEx(gmul(left._gobj, (<Expression>right)._gobj))

    cdef RingElement _div_c_impl(left, RingElement right):
        """
        Divide left and right.

            sage: var("x y", ns=1)
            (x, y)
            sage: x/y/y
            x*y^(-2)
        """
        return new_Expression_from_GEx(gdiv(left._gobj, (<Expression>right)._gobj))

    cdef int _cmp_c_impl(left, Element right) except -2:
        # TODO: this is never called, maybe, since I define
        # richcmp above to make formal symbolic expressions?
        return left._gobj.compare((<Expression>right)._gobj)

    def __cmp__(self, right):
        """
        Compare self and right, returning -1, 0, or 1, depending on if
        self < right, self == right, or self > right, respectively.

        Use this instead of the operators <=, <, etc. to compare symbolic
        expressions when you do not want to get a formal inequality back.

        IMPORTANT: Both self and right *must* have the same type, or
        this function won't be called.

        EXAMPLES:
            sage: x,y = var('x,y', ns=1); S = x.parent()
            sage: x.__cmp__(y)
            -1
            sage: x < y
            x < y
            sage: cmp(x,y)
            -1
            sage: cmp(S(0.5), S(0.7))
            -1
            sage: S(0.5) < S(0.7)
            0.500000000000000 < 0.700000000000000

        This is confusing because 0.7 is not a symbolic expression:
            sage: cmp(S(0.5), 0.7)
            0
            sage: cmp(sin(S(2)), sin(S(1)))
            1
            sage: float(sin(S(2)))
            0.90929742682568171
            sage: float(sin(S(1)))
            0.8414709848078965
        """
        cdef Expression r = self.coerce_in(right)
        return self._gobj.compare(r._gobj)

    def __pow__(Expression self, exp, ignored):
        """
        Return self raised to the power of exp.

        INPUT:
            self -- symbolic expression
            exp -- something that coerces to a symbolic expressions

        OUTPUT:
            symbolic expression

        EXAMPLES:
            sage: var('x,y',ns=1); S=x.parent()
            (x, y)
            sage: x.__pow__(y)
            x^y
            sage: x^(3/5)
            x^(3/5)
            sage: x^sin(x)^cos(y)
            x^(sin(x)^cos(y))
        """
        cdef Expression nexp = self.coerce_in(exp)
        return new_Expression_from_GEx(g_pow(self._gobj, nexp._gobj))

    def diff(self, symb, deg=1):
        """
        Return the deg-th (partial) derivative of self with respect to symb.

        EXAMPLES:
            sage: var("x y", ns=1)
            (x, y)
            sage: b = (x+y)^5
            sage: b.diff(x, 2)
            20*(x + y)^3

            sage: from sage.symbolic.function import function as myfunc
            sage: foo = myfunc('foo',2)
            sage: foo(x^2,x^2).diff(x)
            2*D[0](foo)(x^2,x^2)*x + 2*D[1](foo)(x^2,x^2)*x
        """
        if not isinstance(deg, (int, long, sage.rings.integer.Integer)) \
                or deg < 1:
            raise TypeError, "argument deg should be an integer >1."
        if not isinstance(symb, Expression):
            try:
                symb = ring.NSR(symb)
            except TypeError, err:
                raise TypeError, "argument symb must be a symbol"
        if not is_a_symbol((<Expression>symb)._gobj):
            raise TypeError, "argument symb must be a symbol"
        _sig_on
        cdef GEx x = self._gobj.diff(ex_to_symbol((<Expression>symb)._gobj), deg)
        _sig_off
        return new_Expression_from_GEx(x)

    def expand(Expression self):
        """
        Return expanded form of his expression, obtained by multiplying out
        all products.

        OUTPUT:
            symbolic expression

        EXAMPLES:
            sage: var('x,y',ns=1)
            (x, y)
            sage: ((x + (2/3)*y)^3).expand()
            4/3*x*y^2 + 2*x^2*y + x^3 + 8/27*y^3
            sage: expand( (x*sin(x) - cos(y)/x)^2 )
            -2*sin(x)*cos(y) + sin(x)^2*x^2 + cos(y)^2*x^(-2)
            sage: f = (x-y)*(x+y); f
            (x - y)*(x + y)
            sage: f.expand()
            x^2 - y^2
        """
        _sig_on
        cdef GEx x = self._gobj.expand(0)
        _sig_off
        return new_Expression_from_GEx(x)

    def collect(Expression self, s):
        """
        INPUT:
            s -- a symbol

        OUTPUT:
            expression

        EXAMPLES:
            sage: var('x,y,z',ns=1)
            (x, y, z)
            sage: f = 4*x*y + x*z + 20*y^2 + 21*y*z + 4*z^2 + x^2*y^2*z^2
            sage: f.collect(x)
            (4*y + z)*x + 21*y*z + x^2*y^2*z^2 + 20*y^2 + 4*z^2
            sage: f.collect(y)
            (x^2*z^2 + 20)*y^2 + (4*x + 21*z)*y + x*z + 4*z^2
            sage: f.collect(z)
            (x^2*y^2 + 4)*z^2 + (x + 21*y)*z + 4*x*y + 20*y^2
        """
        cdef Expression s0 = self.coerce_in(s)
        _sig_on
        cdef GEx x = self._gobj.collect(s0._gobj, False)
        _sig_off
        return new_Expression_from_GEx(x)

    def __abs__(self):
        """
        Return the absolute value of this expression.

        EXAMPLES:
            sage: var('x, y', ns=1); S = parent(x)
            (x, y)

        The absolute value of a symbolic expression
            sage: abs(x^2+y^2)
            abs(x^2 + y^2)

        The absolute value of a number in the symbolic ring:
            sage: abs(S(-5))
            5
            sage: type(abs(S(-5)))
            <type 'sage.symbolic.expression.Expression'>
        """
        return new_Expression_from_GEx(g_abs(self._gobj))

    def step(self):
        """
        Return the value of the Heaviside step function, which is 0 for
        negative x, 1/2 for 0, and 1 for positive x.

        EXAMPLES:
            sage: x = var('x',ns=1); SR = x.parent()
            sage: SR(1.5).step()
            1
            sage: SR(0).step()
            1/2
            sage: SR(-1/2).step()
            0
            sage: SR(float(-1)).step()
            0
        """
        return new_Expression_from_GEx(g_step(self._gobj))

    def csgn(self):
        """
        Return the sign of self, which is -1 if self < 0, 0 if self ==
        0, and 1 if self > 0, or unevaluated when self is a nonconstant
        symbolic expression.

        It can be somewhat arbitrary when self is not real.

        EXAMPLES:
            sage: x = var('x', ns=1); SR = x.parent()
            sage: SR(-2).csgn()
            -1
            sage: SR(0.0).csgn()
            0
            sage: SR(10).csgn()
            1
            sage: x.csgn()
            csgn(x)
            sage: SR(CDF.0).csgn()
            1
        """
        return new_Expression_from_GEx(g_csgn(self._gobj))

    def conjugate(self):
        """
        Return the complex conjugate of self.

        EXAMPLES:
            sage: x = var('x', ns=1); SR = x.parent()
            sage: SR(CDF.0).conjugate()
            -I
            sage: x.conjugate()
            conjugate(x)
            sage: SR(RDF(1.5)).conjugate()
            1.5
            sage: SR(float(1.5)).conjugate()
            1.5
            sage: I = SR(CDF.0)
            sage: I.conjugate()
            -I
            sage: ( 1+I  + (2-3*I)*x).conjugate()
            (2.0+3.0*I)*conjugate(x) + 1.0-I
        """
        return new_Expression_from_GEx(g_conjugate(self._gobj))

    def real_part(self):
        """
        Return the real part of self.

        EXAMPLES:
            sage: x = var('x', ns=1); SR = x.parent()
            sage: x.real_part()
            real_part(x)
            sage: SR(CDF(2,3)).real_part()
            2.0
            sage: SR(CC(2,3)).real_part()
            2.00000000000000
        """
        return new_Expression_from_GEx(g_real_part(self._gobj))

    def imag_part(self):
        """
        Return the imaginary part of self.

        EXAMPLES:
            sage: x = var('x', ns=1); SR = x.parent()
            sage: x.imag_part()
            imag_part(x)
            sage: SR(CC(2,3)).imag_part()
            3.00000000000000
            sage: SR(CDF(2,3)).imag_part()
            3.0
        """
        return new_Expression_from_GEx(g_imag_part(self._gobj))

    def sqrt(self):
        """
        EXAMPLES:
            sage: var('x, y', ns=1); S = parent(x)
            (x, y)
            sage: S(2).sqrt()
            sqrt(2)
            sage: (x^2+y^2).sqrt()
            sqrt(x^2 + y^2)
            sage: (x^2).sqrt()
            sqrt(x^2)
        """
        return new_Expression_from_GEx(g_sqrt(self._gobj))

    def sin(self):
        """
        EXAMPLES:
            sage: var('x, y', ns=1); S = parent(x)
            (x, y)
            sage: sin(x^2 + y^2)
            sin(x^2 + y^2)
            sage: sin(sage.symbolic.ring.pi)
            0
            sage: sin(S(1))
            sin(1)
            sage: sin(S(RealField(150)(1)))
            0.84147098480789650665250232163029899962256306
        """
        return new_Expression_from_GEx(g_sin(self._gobj))

    def cos(self):
        """
        Return the cosine of self.

        EXAMPLES:
            sage: var('x, y', ns=1); S = parent(x)
            (x, y)
            sage: cos(x^2 + y^2)
            cos(x^2 + y^2)
            sage: cos(sage.symbolic.ring.pi)
            -1
            sage: cos(S(1))
            cos(1)
            sage: cos(S(RealField(150)(1)))
            0.54030230586813971740093660744297660373231042
            sage: S(RR(1)).cos()
            0.540302305868140
            sage: S(float(1)).cos()
            0.54030230586813977
        """
        return new_Expression_from_GEx(g_cos(self._gobj))

    def tan(self):
        """
        EXAMPLES:
            sage: var('x, y', ns=1); S = parent(x)
            (x, y)
            sage: tan(x^2 + y^2)
            tan(x^2 + y^2)
            sage: tan(sage.symbolic.ring.pi/2)
            tan(1/2*Pi)
            sage: tan(S(1))
            tan(1)
            sage: tan(S(RealField(150)(1)))
            1.5574077246549022305069748074583601730872508
        """
        return new_Expression_from_GEx(g_tan(self._gobj))

    def arcsin(self):
        """
        Return the arcsin of x, i.e., the number y between -pi and pi
        such that sin(y) == x.

        EXAMPLES:
            sage: x = var('x', ns=1); SR = x.parent()
            sage: x.arcsin()
            arcsin(x)
            sage: SR(0.5).arcsin()
            0.523598775598299
            sage: SR(0.999).arcsin()
            1.52607123962616
            sage: SR(-0.999).arcsin()
            -1.52607123962616
        """
        return new_Expression_from_GEx(g_asin(self._gobj))

    def arccos(self):
        """
        Return the arc cosine of self.

        EXAMPLES:
            sage: x = var('x', ns=1); S = x.parent()
            sage: x.arccos()
            arccos(x)
            sage: S(1).arccos()
            0
            sage: S(1/2).arccos()
            1/3*Pi
            sage: S(0.4).arccos()
            1.15927948072741
            sage: plot(lambda x: S(x).arccos(), -1,1)
        """
        return new_Expression_from_GEx(g_acos(self._gobj))

    def arctan(self):
        """
        Return the arc tangent of self.

        EXAMPLES:
            sage: x = var('x', ns=1); S = x.parent()
            sage: x.arctan()
            arctan(x)
            sage: S(1).arctan()
            1/4*Pi
            sage: S(1/2).arctan()
            arctan(1/2)
            sage: S(0.5).arctan()
            0.463647609000806
            sage: plot(lambda x: S(x).arctan(), -20,20)
        """
        return new_Expression_from_GEx(g_atan(self._gobj))

    def arctan2(self, x):
        """
        Return the inverse of the 2-variable tan function on self and x.

        EXAMPLES:
            sage: var('x,y', ns=1); S = parent(x)
            (x, y)
            sage: x.arctan2(y)
            arctan2(x, y)
            sage: S(1/2).arctan2(1/2)
            1/4*Pi
            sage: maxima.eval('atan2(1/2,1/2)')
            '%pi/4'
        """
        cdef Expression nexp = self.coerce_in(x)
        return new_Expression_from_GEx(g_atan2(self._gobj, nexp._gobj))

    def sinh(self):
        r"""
        Return sinh of self.

        We have $\sinh(x) = (e^{x} - e^{-x})/2$.

        EXAMPLES:
            sage: x = var('x', ns=1); S = x.parent()
            sage: x.sinh()
            sinh(x)
            sage: S(1).sinh()
            sinh(1)
            sage: S(0).sinh()
            0
            sage: S(1.0).sinh()
            1.17520119364380
            sage: maxima('sinh(1.0)')
            1.175201193643801
            sage: S(1.0000000000000000000000000).sinh()
            1.1752011936438014568823819
            sage: S(RIF(1)).sinh()
            1.175201193643802?
            sage: plot(lambda x: S(x).sinh(), -1, 1)
        """
        return new_Expression_from_GEx(g_sinh(self._gobj))

    def cosh(self):
        """
        Return cosh of self.

        We have $\sinh(x) = (e^{x} + e^{-x})/2$.

        EXAMPLES:
            sage: x = var('x', ns=1); S = x.parent()
            sage: x.cosh()
            cosh(x)
            sage: S(1).cosh()
            cosh(1)
            sage: S(0).cosh()
            1
            sage: S(1.0).cosh()
            1.54308063481524
            sage: maxima('cosh(1.0)')
            1.543080634815244
            sage: S(1.0000000000000000000000000).cosh()
            1.5430806348152437784779056
            sage: S(RIF(1)).cosh()
            1.543080634815244?
            sage: plot(lambda x: S(x).cosh(), -1, 1)
        """
        return new_Expression_from_GEx(g_cosh(self._gobj))

    def tanh(self):
        """
        Return tanh of self.

        We have $\tanh(x) = \sinh(x) / \cosh(x)$.

        EXAMPLES:
            sage: x = var('x', ns=1); S = x.parent()
            sage: x.tanh()
            tanh(x)
            sage: S(1).tanh()
            tanh(1)
            sage: S(0).tanh()
            0
            sage: S(1.0).tanh()
            0.761594155955765
            sage: maxima('tanh(1.0)')
            .7615941559557649
            sage: plot(lambda x: S(x).tanh(), -1, 1)
        """
        return new_Expression_from_GEx(g_tanh(self._gobj))

    def arcsinh(self):
        """
        Return the inverse hyperbolic sine of self.

        EXAMPLES:
            sage: x = var('x', ns=1); S = x.parent()
            sage: x.arcsinh()
            arcsinh(x)
            sage: S(0).arcsinh()
            0
            sage: S(1).arcsinh()
            arcsinh(1)
            sage: S(1.0).arcsinh()
            1.17520119364380
            sage: maxima('sinh(1.0)')
            1.175201193643801

        Sage automatically applies certain identies:
            sage: S(3/2).arcsinh().cosh()
            sqrt(13/4)
        """
        return new_Expression_from_GEx(g_asinh(self._gobj))

    def arccosh(self):
        """
        Return the inverse hyperbolic cosine of self.

        EXAMPLES:
            sage: x = var('x', ns=1); S = x.parent()
            sage: x.arccosh()
            arccosh(x)
            sage: S(0).arccosh()
            (0.5*I)*Pi
            sage: S(1/2).arccosh()
            arccosh(1/2)
            sage: S(0.5).arccosh()    # TODO -- BUG; somehow in pynac this is now BROKEN!! it works in ginsh
            1.047197551196597746*I
            sage: maxima('acosh(0.5)')
            1.047197551196598*%i
            sage: plot(lambda x: S(x).arccosh(), 0.01, 2)
        """
        return new_Expression_from_GEx(g_acosh(self._gobj))

    def arctanh(self):
        """
        Return the inverse hyperbolic tangent of self.

        EXAMPLES:
            sage: x = var('x', ns=1); S = x.parent()
            sage: x.arctanh()
            arctanh(x)
            sage: S(0).arctanh()
            0
            sage: S(1/2).arctanh()
            arctanh(1/2)
            sage: S(0.5).arctanh()
            0.549306144334055
            sage: S(0.5).arctanh().tanh()
            0.500000000000000
            sage: maxima('atanh(0.5)')
            .5493061443340549
            sage: plot(lambda x: S(x).arctanh(), -0.95,0.95)
        """
        return new_Expression_from_GEx(g_atanh(self._gobj))

    def exp(self):
        """
        Return exponential function of self, i.e., e to the
        power of self.

        EXAMPLES:
            sage: x = var('x', ns=1); S = x.parent()
            sage: x.exp()
            exp(x)
            sage: S(0).exp()
            1
            sage: S(1/2).exp()
            exp(1/2)
            sage: S(0.5).exp()
            1.64872127070013
            sage: S(0.5).exp().log()
            0.500000000000000
            sage: math.exp(0.5)
            1.6487212707001282
            sage: plot(lambda x: S(x).exp(), -2,1)
        """
        return new_Expression_from_GEx(g_exp(self._gobj))

    def log(self):
        """
        Return the logarithm of self.

        EXAMPLES:
            sage: x, y = var('x, y', ns=1); S = x.parent()
            sage: x.log()
            log(x)
            sage: (x^y + y^x).log()
            log(x^y + y^x)
            sage: S(0).log()
            Traceback (most recent call last):
            ...
            RuntimeError: log_eval(): log(0)
            sage: S(1).log()
            0
            sage: S(1/2).log()
            log(1/2)
            sage: S(0.5).log()
            -0.693147180559945
            sage: S(0.5).log().exp()
            0.500000000000000
            sage: math.log(0.5)
            -0.69314718055994529
            sage: plot(lambda x: S(x).log(), 0.1,10)
        """
        return new_Expression_from_GEx(g_log(self._gobj))

    def zeta(self):
        """
        EXAMPLES:
            sage: x, y = var('x, y', ns=1); S = x.parent()
            sage: (x/y).zeta()
            zeta(x*y^(-1))
            sage: S(2).zeta()
            1/6*Pi^2
            sage: S(3).zeta()
            zeta(3)
            sage: S(CDF(0,1)).zeta()
            0.00330022368532-0.418155449141*I
            sage: CDF(0,1).zeta()
            0.00330022368532 - 0.418155449141*I
            sage: plot(lambda x: S(x).zeta(), -10,10).show(ymin=-3,ymax=3)
        """
        _sig_on
        cdef GEx x = g_zeta(self._gobj)
        _sig_off
        return new_Expression_from_GEx(x)

    def factorial(self):
        """
        Return the factorial of self.

        OUTPUT:
            symbolic expression

        EXAMPLES:
            sage: var('x, y', ns=1); S = parent(x)
            (x, y)
            sage: S(5).factorial()
            120
            sage: x.factorial()
            x!
            sage: (x^2+y^3).factorial()
            (x^2 + y^3)!
        """
        _sig_on
        cdef GEx x = g_factorial(self._gobj)
        _sig_off
        return new_Expression_from_GEx(x)

    def binomial(self, k):
        """
        Return binomial coefficient "self choose k".

        OUTPUT:
            symbolic expression

        EXAMPLES:
            sage: var('x, y', ns=1); S = parent(x)
            (x, y)
            sage: S(5).binomial(S(3))
            10
            sage: x.binomial(S(3))
            1/6*x^3 - 1/2*x^2 + 1/3*x
            sage: x.binomial(y)
            binomial(x,y)
        """
        cdef Expression nexp = self.coerce_in(k)
        _sig_on
        cdef GEx x = g_binomial(self._gobj, nexp._gobj)
        _sig_off
        return new_Expression_from_GEx(x)

    def Order(self):
        """
        Order, as in big oh notation.

        OUTPUT:
            symbolic expression

        EXAMPLES:
            sage: n = var('n', ns=1)
            sage: (17*n^3).Order()
            Order(n^3)
        """
        return new_Expression_from_GEx(g_Order(self._gobj))

    def gamma(self):
        """
        Return the Gamma function evaluated at self.

        EXAMPLES:
            sage: x = var('x', ns=1); S = x.parent()
            sage: x.gamma()
            gamma(x)
            sage: S(2).gamma()
            1
            sage: S(10).gamma()
            362880
            sage: S(10.0).gamma()
            362880.000000000
            sage: S(CDF(1,1)).gamma()
            0.498015668118-0.154949828302*I
            sage: gp('gamma(1+I)')
            0.4980156681183560427136911175 - 0.1549498283018106851249551305*I
            sage: set_verbose(-1); plot(lambda x: S(x).gamma(), -6,4).show(ymin=-3,ymax=3)
        """
        _sig_on
        cdef GEx x = g_tgamma(self._gobj)
        _sig_off
        return new_Expression_from_GEx(x)

    def lgamma(self):
        """
        Return the log-gamma function evaluated at self.
        This is the logarithm of gamma of self, where
        gamma is a complex function such that gamma(n)
        equals factorial(n-1).

        EXAMPLES:
            sage: x = var('x', ns=1); S = x.parent()
            sage: x.lgamma()
            lgamma(x)
            sage: S(2).lgamma()
            0
            sage: S(5).lgamma()
            log(24)
            sage: S(5-1).factorial().log()
            log(24)
            sage: set_verbose(-1); plot(lambda x: S(x).lgamma(), -7,8, plot_points=1000).show()
            sage: math.exp(0.5)
            1.6487212707001282
            sage: plot(lambda x: (S(x).exp() - S(-x).exp())/2 - S(x).sinh(), -1, 1)
        """
        _sig_on
        cdef GEx x = g_lgamma(self._gobj)
        _sig_off
        return new_Expression_from_GEx(x)

    # Functions to add later, maybe.  These were in Ginac mainly
    # implemented using a lot from cln, and I had to mostly delete
    # their implementations.   They are pretty specialized for
    # physics apps, maybe.
    # This doesn't work / isn't implemented yet / just segfaults.
    #def Li(self, x):
    #    """
    #    """
    #    cdef Expression nexp = self.coerce_in(x)
    #    return new_Expression_from_GEx(g_Li(self._gobj, nexp._gobj))
    #def Li2(self):
    #    return new_Expression_from_GEx(g_Li2(self._gobj))
    #def G(self, Expression y):
    #    return new_Expression_from_GEx(g_G(self._gobj, y._gobj))
    #def G2(self, Expression s, Expression y):
    #    return new_Expression_from_GEx(g_G2(self._gobj, s._gobj, y._gobj))
    #def S(self, Expression p, Expression x):
    #return new_Expression_from_GEx(g_S(self._gobj, p._gobj, x._gobj))
    #def H(self, Expression x):
    #return new_Expression_from_GEx(g_H(self._gobj, x._gobj))
    #def zeta2(self, Expression s):
    #    return new_Expression_from_GEx(g_zeta2(self._gobj, s._gobj))
    #def zetaderiv(self, Expression x):
    #    return new_Expression_from_GEx(g_zetaderiv(self._gobj, x._gobj))
    #def beta(self, Expression y):
    #    return new_Expression_from_GEx(g_beta(self._gobj, y._gobj))
    #def psi(self):
    #    return new_Expression_from_GEx(g_psi(self._gobj))
    #def psi2(self, Expression x):
    #    return new_Expression_from_GEx(g_psi2(self._gobj, x._gobj))



cdef Expression new_Expression_from_GEx(GEx juice):
    cdef Expression nex
    nex = <Expression>PY_NEW(Expression)
    GEx_construct_ex(&nex._gobj, juice)
    nex._parent = ring.NSR
    return nex


