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

cdef dict sfunction_serial_dict = {}

from pynac import get_ginac_serial

from sage.misc.fpickle import pickle_function, unpickle_function

cdef class SFunction:
    """
    Return a formal symbolic function.

    EXAMPLES:

    We create a formal function of one variable, write down
    an expression that involves first and second derivatives,
    and extract off coefficients.::

        sage: from sage.symbolic.function import function
        sage: var('r,kappa', ns=1)
        (r, kappa)
        sage: psi = function('psi', 1)(r); psi
        psi(r)
        sage: g = 1/r^2*(2*r*psi.derivative(r,1) + r^2*psi.derivative(r,2)); g
        (r^2*D[0, 0](psi)(r) + 2*r*D[0](psi)(r))/r^2
        sage: g.expand()
        2*D[0](psi)(r)/r + D[0, 0](psi)(r)
        sage: g.coeff(psi.derivative(r,2))
        1
        sage: g.coeff(psi.derivative(r,1))
        2/r
    """
    def __init__(self, name, nargs=0, eval_func=None, evalf_func=None,
            conjugate_func=None, real_part_func=None, imag_part_func=None,
            derivative_func=None, power_func=None, series_func=None,
            latex_name=None, print_func=None, print_latex_func=None):
        """
        EXAMPLES::

            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: x,y,z = var("x y z", ns=1)
            sage: foo(x, y) + foo(y, z)^2
            foo(y, z)^2 + foo(x, y)

            sage: def ev(x): return 2*x
            sage: foo = nfunction("foo", 1, eval_func=ev)
            sage: foo(x)
            2*x
            sage: foo = nfunction("foo", 1, eval_func=lambda x: 5)
            sage: foo(x)
            5

            sage: def evalf_f(x, prec=0): return 6
            sage: foo = nfunction("foo", 1, evalf_func=evalf_f)
            sage: foo(x)
            foo(x)
            sage: foo(x).n()
            6

            sage: foo = nfunction("foo", 1, conjugate_func=ev)
            sage: foo(x).conjugate()
            2*x

            sage: def deriv(*args,**kwds): print args, kwds; return args[kwds['diff_param']]^2
            sage: foo = nfunction("foo", 2, derivative_func=deriv)
            sage: foo(x,y).derivative(y)
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

            sage: def my_print(*args): return "my args are: " + ', '.join(map(repr, args))
            sage: foo = nfunction('t', 2, print_func=my_print)
            sage: foo(x,y^z)
            my args are: x, y^z

            sage: latex(foo(x,y^z))
            \mbox{t}\left(x, y^{z}\right)
            sage: foo = nfunction('t', 2, print_latex_func=my_print)
            sage: foo(x,y^z)
            t(x, y^z)
            sage: latex(foo(x,y^z))
            my args are: x, y^z
            sage: foo = nfunction('t', 2, latex_name='foo')
            sage: latex(foo(x,y^z))
            \mbox{foo}\left(x, y^{z}\right)


        TESTS::

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
        self.name = name
        if nargs:
            self.nargs = nargs
        else:
            self.nargs = 0

        self.eval_f = eval_func
        if eval_func:
            if not callable(eval_func):
                raise ValueError, "eval_func parameter must be callable"

        self.evalf_f = evalf_func
        if evalf_func:
            if not callable(evalf_func):
                raise ValueError, "evalf_func parameter must be callable"

        self.conjugate_f = conjugate_func
        if conjugate_func:
            if not callable(conjugate_func):
                raise ValueError, "conjugate_func parameter must be callable"

        self.real_part_f = real_part_func
        if real_part_func:
            if not callable(real_part_func):
                raise ValueError, "real_part_func parameter must be callable"

        self.imag_part_f = imag_part_func
        if imag_part_func:
            if not callable(imag_part_func):
                raise ValueError, "imag_part_func parameter must be callable"

        self.derivative_f = derivative_func
        if derivative_func:
            if not callable(derivative_func):
                raise ValueError, "derivative_func parameter must be callable"

        self.power_f = power_func
        if power_func:
            if not callable(power_func):
                raise ValueError, "power_func parameter must be callable"

        self.series_f = series_func
        if series_func:
            if not callable(series_func):
                raise ValueError, "series_func parameter must be callable"

        # handle custom printing
        # if print_func is defined, it is used instead of name
        # latex printing can be customised either by setting a string latex_name
        # or giving a custom function argument print_latex_func
        if latex_name and print_latex_func:
            raise ValueError, "only one of latex_name and print_latex_name should be specified."

        self.print_latex_f = print_latex_func
        self.latex_name = latex_name
        self.print_f = print_func

        self._init_()

    # this is separated from the constructor since it is also called
    # during unpickling
    cdef _init_(self):
        # see if there is already an SFunction with the same state
        cdef SFunction sfunc
        cdef long myhash = self._hash_()
        for sfunc in sfunction_serial_dict.itervalues():
            if myhash == sfunc._hash_():
                # found one, set self.serial to be a copy
                self.serial = sfunc.serial
                return

        cdef GFunctionOpt opt

        opt = g_function_options_args(self.name, self.nargs)
        opt.set_python_func()

        if self.eval_f:
            opt.eval_func(self.eval_f)

        if self.evalf_f:
            opt.evalf_func(self.evalf_f)

        if self.conjugate_f:
            opt.conjugate_func(self.conjugate_f)

        if self.real_part_f:
            opt.real_part_func(self.real_part_f)

        if self.imag_part_f:
            opt.imag_part_func(self.imag_part_f)

        if self.derivative_f:
            opt.derivative_func(self.derivative_f)

        if self.power_f:
            opt.power_func(self.power_f)

        if self.series_f:
            opt.series_func(self.series_f)

        if self.print_latex_f:
            opt.set_print_latex_func(self.print_latex_f)

        if self.latex_name:
            opt.latex_name(self.latex_name)

        if self.print_f:
            opt.set_print_dflt_func(self.print_f)

        self.serial = g_register_new(opt)

        g_foptions_assign(g_registered_functions().index(self.serial), opt)
        global sfunction_serial_dict
        sfunction_serial_dict[self.serial] = self

        self.__hinit = False

    # cache the hash value of this function
    # this is used very often while unpickling to see if there is already
    # a function with the same properties
    cdef long _hash_(self):
        if not self.__hinit:
            # create a string representation of this SFunction
            slist = [self.nargs, self.name, str(self.latex_name)]
            for f in [self.eval_f, self.evalf_f, self.conjugate_f,
                    self.real_part_f, self.imag_part_f, self.derivative_f,
                    self.power_f, self.series_f, self.print_f,
                    self.print_latex_f]:
                if f:
                    slist.append(hash(f.func_code))
                else:
                    slist.append(' ')
            self.__hcache = hash(tuple(slist))
            self.__hinit = True
        return self.__hcache

    def __hash__(self):
        """
        TESTS::
            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: hash(foo)
            7313648655953480146

            sage: def ev(x): return 2*x
            sage: foo = nfunction("foo", 2, eval_func = ev)
            sage: hash(foo)
            4884169210301491732

        """
        return self.serial*self._hash_()

    def __getstate__(self):
        """
        Returns a tuple describing the state of this object for pickling.

        Pickling SFunction objects is limited by the ability to pickle
        functions in python. We use sage.misc.fpickle.pickle_function for
        this purpose, which only works if there are no nested functions.


        This should return all information that will be required to unpickle
        the object. The functionality for unpickling is implemented in
        __setstate__().

        In order to pickle SFunction objects, we return a tuple containing

         * 0  - as pickle version number
                in case we decide to change the pickle format in the feature
         * name of this function
         * number of arguments
         * latex_name
         * a tuple containing attempts to pickle the following optional
           functions, in the order below
           * eval_f
           * evalf_f
           * conjugate_f
           * real_part_f
           * imag_part_f
           * derivative_f
           * power_f
           * series_f
           * print_f
           * print_latex_f

        EXAMPLES::
            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: foo.__getstate__()
            (0, 'foo', 2, None, [None, None, None, None, None, None, None, None, None, None])
            sage: t = loads(dumps(foo))
            sage: t == foo
            True
            sage: var('x,y',ns=1)
            (x, y)
            sage: t(x,y)
            foo(x, y)

            sage: def ev(x,y): return 2*x
            sage: foo = nfunction("foo", 2, eval_func = ev)
            sage: foo.__getstate__()
            (0, 'foo', 2, None, ["...", None, None, None, None, None, None, None, None, None])

            sage: u = loads(dumps(foo))
            sage: u == foo
            True
            sage: t == u
            False
            sage: u(y,x)
            2*y

            sage: def evalf_f(x, prec=0): return int(6)
            sage: foo = nfunction("foo", 1, evalf_func=evalf_f)
            sage: foo.__getstate__()
            (0, 'foo', 1, None, [None, "...", None, None, None, None, None, None, None, None])

            sage: v = loads(dumps(foo))
            sage: v == foo
            True
            sage: v == u
            False
            sage: foo(y).n()
            6
            sage: v(y).n()
            6

        Test pickling expressions with symbolic functions::

            sage: u = loads(dumps(foo(x)^2 + foo(y) + x^y)); u
            x^y + foo(x)^2 + foo(y)
            sage: u.subs(y=0)
            foo(x)^2 + foo(0) + 1
            sage: u.subs(y=0).n()
            43.0000000000000
        """
        return (0, self.name, self.nargs, self.latex_name,
                map(pickle_wrapper, [self.eval_f, self.evalf_f,
                    self.conjugate_f, self.real_part_f, self.imag_part_f,
                    self.derivative_f, self.power_f, self.series_f,
                    self.print_f, self.print_latex_f]))

    def __setstate__(self, state):
        """
        Initializes the state of the object from data saved in a pickle.

        During unpickling __init__ methods of classes are not called, the saved
        data is passed to the class via this function instead.

        TESTS::

            sage: from sage.symbolic.function import function as nfunction
            sage: var('x,y', ns=1)
            (x, y)
            sage: foo = nfunction("foo", 2)
            sage: bar = nfunction("bar", 1)
            sage: bar.__setstate__(foo.__getstate__())

        Note that the other direction doesn't work here, since foo._hash_()
        hash already been initialized.::

            sage: bar
            foo
            sage: bar(x,y)
            foo(x, y)
        """
        # check input
        if state[0] != 0 or len(state) != 5:
            raise ValueError, "unknown state information"

        self.name = state[1]
        self.nargs = state[2]
        self.latex_name = state[3]
        self.eval_f = unpickle_wrapper(state[4][0])
        self.evalf_f = unpickle_wrapper(state[4][1])
        self.conjugate_f = unpickle_wrapper(state[4][2])
        self.real_part_f = unpickle_wrapper(state[4][3])
        self.imag_part_f = unpickle_wrapper(state[4][4])
        self.derivative_f = unpickle_wrapper(state[4][5])
        self.power_f = unpickle_wrapper(state[4][6])
        self.series_f = unpickle_wrapper(state[4][7])
        self.print_f = unpickle_wrapper(state[4][8])
        self.print_latex_f = unpickle_wrapper(state[4][9])

        self._init_()

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
        EXAMPLES::

            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: x,y,z = var("x y z", ns=1)
            sage: foo(x,y)
            foo(x, y)

            sage: foo(y)
            Traceback (most recent call last):
            ...
            ValueError: Symbolic function foo expects 2 arguments, got 1.

            sage: bar = nfunction("bar")
            sage: bar(x)
            bar(x)
            sage: bar(x,y)
            bar(x, y)

        TESTS::

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

def get_sfunction_from_serial(serial):
    global sfunction_serial_dict
    return sfunction_serial_dict.get(serial)

import base64
def pickle_wrapper(f):
    if f is None:
        return None
    return pickle_function(f)

def unpickle_wrapper(p):
    if p is None:
        return None
    return unpickle_function(p)

def init_sfunction_map():
    """
    Initializes a list mapping GiNaC function serials to the equivalent Sage
    functions.

    EXAMPLES:
        sage: from sage.symbolic.function import init_sfunction_map
        sage: f_list = init_sfunction_map()
        sage: gamma in f_list
        True
        sage: exp in f_list
        True
        sage: x = var('x',ns=1)
        sage: sin(x).operator()
        sin
        sage: gamma(x).operator()
        gamma
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
        sage: binomial in f_list
        True
        sage: gamma in f_list
        True
        sage: cos in f_list
        True
    """
    global sfunction_map
    if sfunction_map is None:
        sfunction_map = init_sfunction_map()
    return sfunction_map
