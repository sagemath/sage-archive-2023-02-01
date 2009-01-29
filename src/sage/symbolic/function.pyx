###############################################################################
#   SAGE: Open Source Mathematical Software
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

r"""

Support for symbolic functions.

"""
include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"
include "../libs/ginac/decl.pxi"


from expression cimport new_Expression_from_GEx, Expression
from ring import NSR as SR

from pynac import get_ginac_serial

cdef class SFunction:
    """
    Return a formal symbolic function.

    EXAMPLES:
    We create a formal function of one variable, write down
    an expression that involves first and second derivatives,
    and extract off coefficients.
        sage: from sage.symbolic.function import function
        sage: var('r,kappa', ns=1)
        (r, kappa)
        sage: psi = function('psi', 1)(r); psi
        psi(r)
        sage: g = 1/r^2*(2*r*psi.diff(r,1) + r^2*psi.diff(r,2)); g
        (2*D[0](psi)(r)*r + D[0,0](psi)(r)*r^2)*r^(-2)
        sage: g.expand()
        D[0,0](psi)(r) + 2*D[0](psi)(r)*r^(-1)
        sage: g.coeff(psi.diff(r,2))
        1
        sage: g.coeff(psi.diff(r,1))
        2*r^(-1)
    """
    cdef unsigned int serial
    cdef object name
    cdef int nargs

    cdef object eval_f
    cdef object evalf_f
    cdef object conjugate_f
    cdef object real_part_f
    cdef object imag_part_f
    cdef object derivative_f
    cdef object power_f
    cdef object series_f

    def __init__(self, name, nargs=0, eval_func=None, evalf_func=None,
            conjugate_func=None, real_part_func=None, imag_part_func=None,
            derivative_func=None, power_func=None, series_func=None):
        """
        EXAMPLES:
            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: x,y,z = var("x y z", ns=1)
            sage: foo(x,y) + foo(y,z)^2
            foo(x,y) + foo(y,z)^2

            sage: def ev(x): return 2*x
            sage: foo = nfunction("foo", 1, eval_func=ev)
            sage: foo(x)
            2*x
            sage: foo = nfunction("foo", 1, eval_func=lambda x: 5)
            sage: foo(x)
            5

            sage: foo = nfunction("foo", 1, evalf_func=lambda x: 5)
            sage: foo(x)
            foo(x)
            sage: foo(x).n()
            5

            sage: foo = nfunction("foo", 1, conjugate_func=ev)
            sage: foo(x).conjugate()
            2*x

            sage: def deriv(*args,**kwds): print args, kwds; return args[kwds['diff_param']]^2
            sage: foo = nfunction("foo", 2, derivative_func=deriv)
            sage: foo(x,y).diff(y)
            (x, y) {'diff_param': 1}
            y^2

            sage: def pow(x, power_param=None): print x, power_param; return x*power_param
            sage: foo = nfunction("foo", 1, power_func=pow)
            sage: foo(y)^(x+y)
            y x + y
            (x + y)*y

            sage: def expand(*args, **kwds): print args, kwds; return sum(args[0]^i for i in range(kwds['order']))
            sage: foo = nfunction("foo", 1, series_func=expand)
            sage: foo(y).series(y, 5)
            (y,) {'var': y, 'options': 0, 'at': 0, 'order': 5}
            y^4 + y^3 + y^2 + y + 1


        TESTS:
            # eval_func raises exception
            sage: def ef(x): raise RuntimeError, "foo"
            sage: bar = nfunction("bar", 1, eval_func=ef)
            sage: bar(x)
            Traceback (most recent call last):
            ...
            RuntimeError: foo

            # eval_func returns None
            sage: def ef(x): pass
            sage: bar = nfunction("bar", 1, eval_func=ef)
            sage: bar(x)
            Traceback (most recent call last):
            ...
            TypeError: eval function returned None, expected return value of type sage.symbolic.expression.Expression

            # eval_func returns non coercable
            sage: def ef(x): return ZZ
            sage: bar = nfunction("bar", 1, eval_func=ef)
            sage: bar(x)
            Traceback (most recent call last):
            ...
            TypeError: eval function did not return a symbolic expression or an element that can be coerced into a symbolic expression

            # eval_func is not callable
            sage: bar = nfunction("bar", 1, eval_func=5)
            Traceback (most recent call last):
            ...
            ValueError: eval_func parameter must be callable

        """
        cdef GFunctionOpt opt
        serial = -1
        try:
            self.serial = find_function(name, nargs)
        except ValueError, err:
            pass

        opt = g_function_options_args(name, nargs)
        opt.set_python_func()
        if eval_func:
            if not callable(eval_func):
                raise ValueError, "eval_func parameter must be callable"
            self.eval_f = eval_func
            opt.eval_func(eval_func)
        if evalf_func:
            if not callable(evalf_func):
                raise ValueError, "evalf_func parameter must be callable"
            self.evalf_f = evalf_func
            opt.evalf_func(evalf_func)
        if conjugate_func:
            if not callable(conjugate_func):
                raise ValueError, "conjugate_func parameter must be callable"
            self.conjugate_f = conjugate_func
            opt.conjugate_func(conjugate_func)
        if real_part_func:
            if not callable(real_part_func):
                raise ValueError, "real_part_func parameter must be callable"
            self.real_part_f = real_part_func
            opt.real_part_func(real_part_func)
        if imag_part_func:
            if not callable(imag_part_func):
                raise ValueError, "imag_part_func parameter must be callable"
            self.imag_part_f = imag_part_func
            opt.imag_part_func(imag_part_func)
        if derivative_func:
            if not callable(derivative_func):
                raise ValueError, "derivative_func parameter must be callable"
            self.derivative_f = derivative_func
            opt.derivative_func(derivative_func)
        if power_func:
            if not callable(power_func):
                raise ValueError, "power_func parameter must be callable"
            self.power_f = power_func
            opt.power_func(power_func)
        if series_func:
            if not callable(series_func):
                raise ValueError, "series_func parameter must be callable"
            self.series_f = series_func
            opt.series_func(series_func)


        if serial is -1 or serial < get_ginac_serial():
            self.serial = g_register_new(opt)
        else:
            self.serial = serial

        g_foptions_assign(g_registered_functions().index(self.serial), opt)

        self.name = name
        if nargs:
            self.nargs = nargs
        else:
            self.nargs = 0

    def __repr__(self):
        """
        EXAMPLES:
            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2); foo
            foo
        """
        return self.name

    def __cmp__(self, other):
        """

        TESTS:
            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: foo == foo
            True
            sage: foo == 2
            False
            sage: foo(1,2).operator() == foo
            True

        """
        if PY_TYPE_CHECK(other, SFunction):
            return cmp(self.serial, (<SFunction>other).serial)
        return False

    def __call__(self, *args, coerce=True):
        """
        EXAMPLES:
            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: x,y,z = var("x y z", ns=1)
            sage: foo(x,y)
            foo(x,y)

            sage: foo(y)
            Traceback (most recent call last):
            ...
            ValueError: Symbolic function foo expects 2 arguments, got 1.

            sage: bar = nfunction("bar")
            sage: bar(x)
            bar(x)
            sage: bar(x,y)
            bar(x,y)

        TESTS:
            # test coercion
            sage: bar(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce arguments:

        """
        cdef GExVector vec
        if self.nargs > 0 and len(args) != self.nargs:
            raise ValueError, "Symbolic function %s expects %s arguments, got %s."%(self.name, self.nargs, len(args))

        if coerce:
            try:
                args = map(SR._coerce_, args)
            except TypeError, err:
                raise TypeError, "cannot coerce arguments: %s"%(err)

        if self.nargs == 0:
            for i from 0 <= i < len(args):
                vec.push_back((<Expression>args[i])._gobj)
            return new_Expression_from_GEx(g_function_evalv(self.serial, vec))

        elif self.nargs == 1:
            return new_Expression_from_GEx(g_function_eval1(self.serial,
                    (<Expression>args[0])._gobj))
        elif self.nargs == 2:
            return new_Expression_from_GEx(g_function_eval2(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj))
        elif self.nargs == 3:
            return new_Expression_from_GEx(g_function_eval3(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj))
        elif self.nargs == 4:
            return new_Expression_from_GEx(g_function_eval4(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj))
        elif self.nargs == 5:
            return new_Expression_from_GEx(g_function_eval5(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj))
        elif self.nargs == 6:
            return new_Expression_from_GEx(g_function_eval6(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj))
        elif self.nargs == 7:
            return new_Expression_from_GEx(g_function_eval7(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj))
        elif self.nargs == 8:
            return new_Expression_from_GEx(g_function_eval8(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj))
        elif self.nargs == 9:
            return new_Expression_from_GEx(g_function_eval9(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj))
        elif self.nargs == 10:
            return new_Expression_from_GEx(g_function_eval10(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj))
        elif self.nargs == 11:
            return new_Expression_from_GEx(g_function_eval11(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj))
        elif self.nargs == 12:
            return new_Expression_from_GEx(g_function_eval12(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj, (<Expression>args[11])._gobj))
        elif self.nargs == 13:
            return new_Expression_from_GEx(g_function_eval13(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj, (<Expression>args[11])._gobj,
                    (<Expression>args[12])._gobj))
        elif self.nargs == 14:
            return new_Expression_from_GEx(g_function_eval14(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj, (<Expression>args[11])._gobj,
                    (<Expression>args[12])._gobj, (<Expression>args[13])._gobj))

function = SFunction

cdef new_SFunction_from_serial(int serial, char* name, int nargs):
    cdef SFunction f = PY_NEW(SFunction)
    f.serial = serial
    f.name = name
    f.nargs = nargs
    return f

def init_sfunction_map():
    """
    Initializes a list mapping GiNaC function serials to the equivalent Sage
    functions.

    EXAMPLES:
        sage: from sage.symbolic.function import init_sfunction_map
        sage: f_list = init_sfunction_map()
        sage: f_list[8] # indices here depend on the GiNaC library
        gamma
        sage: f_list[12]
        exp
    """
    # colliding serials (variable number of arguments): psi2, G2

    f_list = [None]*get_ginac_serial()

    f_list[abs_serial] = abs
#    f_list[step_serial] = # step function
#    f_list[csgn_serial] = # complex sign
#    f_list[conjugate_serial] = # complex conjugation
#    f_list[real_part_serial] = # real part
#    f_list[imag_part_serial] = # imaginary part
    from sage.calculus.calculus import sin, cos, tan, asin, acos, \
            atan, sinh, cosh, tanh, asinh, acosh, atanh, exp, log, gamma, \
            factorial
    f_list[sin_serial] = sin # sine
    f_list[cos_serial] = cos # cosine
    f_list[tan_serial] = tan # tangent
    f_list[asin_serial] = asin # inverse sine
    f_list[acos_serial] = acos # inverse cosine
    f_list[atan_serial] = atan # inverse tangent
    f_list[atan2_serial] = atan # inverse tangent with two arguments
    f_list[sinh_serial] = sinh # hyperbolic sine
    f_list[cosh_serial] = cosh # hyperbolic cosine
    f_list[tanh_serial] = tanh # hyperbolic tangent
    f_list[asinh_serial] = asinh # inverse hyperbolic sine
    f_list[acosh_serial] = acosh # inverse hyperbolic cosine
    f_list[atanh_serial] = atanh # inverse hyperbolic tangent
    f_list[exp_serial] = exp # exponential function
    f_list[log_serial] = log # natural logarithm
#    f_list[Li2_serial] = # dilogarithm
#    f_list[Li_serial] = # classical polylogarithm as well as multiple polylogarithm
#    f_list[G_serial] = # multiple polylogarithm
        # G2_serial = # multiple polylogarithm with explicit signs for the imaginary parts
#    f_list[S_serial] = # Nielsen's generalized polylogarithm
#    f_list[H_serial] = # harmonic polylogarithm
    from sage.functions.transcendental import zeta
    f_list[zeta1_serial] = zeta # Riemann's zeta function as well as multiple zeta value
#    f_list[zeta2_serial] = # alternating Euler sum
#    f_list[zetaderiv_serial] = # derivatives of Riemann's zeta function
    f_list[tgamma_serial] = gamma # gamma function
#    f_list[lgamma_serial] = # logarithm of gamma function
#    f_list[beta_serial] = # beta function (tgamma*tgamma(y)/tgamma(x+y))
#    f_list[psi_serial] = # psi (digamma) function
        # psi2_serial = # derivatives of psi function (polygamma functions)
    from sage.rings.arith import binomial
    f_list[factorial_serial] = factorial # factorial function n!
    f_list[binomial_serial] = binomial # binomial coefficients
    from sage.rings.big_oh import O
    f_list[Order_serial] = O # order term function in truncated power series

    return f_list

sfunction_map = None
def get_sfunction_map():
    """
    Returns the mapping between GiNaC function serials and Sage functions,
    initializing it if necessary.

    EXAMPLES:
        sage: from sage.symbolic.function import get_sfunction_map
        sage: f_list = get_sfunction_map()
        sage: f_list[14]
        sin
        sage: f_list[15]
        cos
    """
    global sfunction_map
    if sfunction_map is None:
        sfunction_map = init_sfunction_map()
    return sfunction_map
