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


from sage.structure.sage_object cimport SageObject
from expression cimport new_Expression_from_GEx, Expression
from ring import SR

cdef dict sfunction_serial_dict = {}

from pynac import get_ginac_serial

from sage.misc.fpickle import pickle_function, unpickle_function

sfunctions_funcs = dict([(name, True) for name in
                         ['eval', 'evalf', 'conjugate', 'real_part',
                          'imag_part', 'derivative', 'power', 'series',
                          'print', 'print_latex']])

cdef class SFunction(SageObject):
    """
    Return a formal symbolic function.

    EXAMPLES:

    We create a formal function of one variable, write down
    an expression that involves first and second derivatives,
    and extract off coefficients.

    ::

        sage: from sage.symbolic.function import function
        sage: r, kappa = var('r,kappa')
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
            latex_name=None, print_func=None, print_latex_func=None,
            ginac_name=None, built_in_function=False):
        """
        EXAMPLES::

            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: x,y,z = var("x y z")
            sage: foo(x, y) + foo(y, z)^2
            foo(y, z)^2 + foo(x, y)

            sage: def ev(x): return 2*x
            sage: foo = nfunction("foo", 1, eval_func=ev)
            sage: foo(x)
            2*x
            sage: foo = nfunction("foo", 1, eval_func=lambda x: 5)
            sage: foo(x)
            5
            sage: def ef(x): pass
            sage: bar = nfunction("bar", 1, eval_func=ef)
            sage: bar(x)
            bar(x)

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
        self._name = name
        self._nargs = nargs
        self.__dict__ = {}

        # handle custom printing
        # if print_func is defined, it is used instead of name
        # latex printing can be customised either by setting a string latex_name
        # or giving a custom function argument print_latex_func
        if latex_name and print_latex_func:
            raise ValueError, "only one of latex_name and print_latex_name should be specified."

        self._ginac_name = ginac_name if ginac_name is not None else self._name
        self._latex_name = latex_name if latex_name else None
        self._built_in_function = built_in_function

        l = locals()
        for func_name in sfunctions_funcs:
            func = l.get(func_name+"_func", None)
            if func and not callable(func):
                raise ValueError, func_name + "_func" + " parameter must be callable"

        if eval_func:
            self._eval_ = eval_func
        if evalf_func:
            self._evalf_ = evalf_func
        if conjugate_func:
            self._conjugate_ = conjugate_func
        if real_part_func:
            self._real_part_ = real_part_func
        if imag_part_func:
            self._imag_part_ = imag_part_func
        if derivative_func:
            self._derivative_ = derivative_func
        if power_func:
            self._power_ = power_func
        if series_func:
            self._series_ = series_func
        if print_func:
            self._print_ = print_func
        if print_latex_func:
            self._print_latex_ = print_latex_func

        self._init_()

    # this is separated from the constructor since it is also called
    # during unpickling
    cdef _init_(self):
        # see if there is already an SFunction with the same state
        cdef SFunction sfunc
        cdef long myhash = self._hash_()
        for sfunc in sfunction_serial_dict.itervalues():
            if myhash == sfunc._hash_():
                # found one, set self._serial to be a copy
                self._serial = sfunc._serial
                return

        cdef GFunctionOpt opt
        opt = g_function_options_args(self._name, self._nargs)
        opt.set_python_func()

        if self._eval_:
            opt.eval_func(self._eval_)

        if self._evalf_:
            opt.evalf_func(self._evalf_)

        if self._conjugate_:
            opt.conjugate_func(self._conjugate_)

        if self._real_part_:
            opt.real_part_func(self._real_part_)

        if self._imag_part_:
            opt.imag_part_func(self._imag_part_)

        if self._derivative_:
            opt.derivative_func(self._derivative_)

        if self._power_:
            opt.power_func(self._power_)

        if self._series_:
            opt.series_func(self._series_)

        if self._print_latex_:
            opt.set_print_latex_func(self._print_latex_)

        if self._latex_name:
            opt.latex_name(self._latex_name)

        if self._print_:
            opt.set_print_dflt_func(self._print_)

        serial = -1
        try:
            if self._built_in_function:
                serial = find_function(self._ginac_name, self._nargs)
        except ValueError, err:
            pass

        if serial is -1 or serial > get_ginac_serial():
            self._serial = g_register_new(opt)
            g_foptions_assign(g_registered_functions().index(self._serial), opt)
        else:
            self._serial = serial

        global sfunction_serial_dict
        sfunction_serial_dict[self._serial] = self

        from sage.symbolic.pynac import symbol_table
        symbol_table['functions'][self._name] = self

        self.__hinit = False


    # cache the hash value of this function
    # this is used very often while unpickling to see if there is already
    # a function with the same properties
    cdef long _hash_(self):
        if not self.__hinit:
            # create a string representation of this SFunction
            slist = [self._nargs, self._name, str(self._latex_name)]
            for f in [self._eval_, self._evalf_, self._conjugate_,
                      self._real_part_, self._imag_part_, self._derivative_,
                      self._power_, self._series_, self._print_,
                      self._print_latex_]:
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
            sage: hash(foo)      # random output
            -6859868030555295348

            sage: def ev(x): return 2*x
            sage: foo = nfunction("foo", 2, eval_func = ev)
            sage: hash(foo)      # random output
            -6859868030555295348
        """
        return self._serial*self._hash_()

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
            sage: var('x,y')
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
        return (0, self._name, self._nargs, self._latex_name,
                map(pickle_wrapper, [self._eval_, self._evalf_,
                    self._conjugate_, self._real_part_, self._imag_part_,
                    self._derivative_, self._power_, self._series_,
                    self._print_, self._print_latex_]))

    def __setstate__(self, state):
        """
        Initializes the state of the object from data saved in a pickle.

        During unpickling __init__ methods of classes are not called, the saved
        data is passed to the class via this function instead.

        TESTS::

            sage: from sage.symbolic.function import function as function
            sage: var('x,y')
            (x, y)
            sage: foo = function("foo", 2)
            sage: bar = function("bar", 1)
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

        self.__dict__ = {}
        self._name = state[1]
        self._nargs = state[2]
        self._latex_name = state[3]
        self._eval_ = unpickle_wrapper(state[4][0])
        self._evalf_ = unpickle_wrapper(state[4][1])
        self._conjugate_ = unpickle_wrapper(state[4][2])
        self._real_part_ = unpickle_wrapper(state[4][3])
        self._imag_part_ = unpickle_wrapper(state[4][4])
        self._derivative_ = unpickle_wrapper(state[4][5])
        self._power_ = unpickle_wrapper(state[4][6])
        self._series_ = unpickle_wrapper(state[4][7])
        self._print_ = unpickle_wrapper(state[4][8])
        self._print_latex_ = unpickle_wrapper(state[4][9])
        self._built_in_function = False

        self._init_()


    def name(self):
        """
        Returns the name of this function.

        EXAMPLES::

            sage: from sage.symbolic.function import function as function
            sage: foo = function("foo", 2)
            sage: foo.name()
            'foo'
        """
        return self._name

    def number_of_arguments(self):
        """
        Returns the number of arguments that this function takes.

        EXAMPLES::

            sage: from sage.symbolic.function import function as function
            sage: foo = function("foo", 2)
            sage: foo.number_of_arguments()
            2
            sage: foo(x,x)
            foo(x, x)

            sage: foo(x)
            Traceback (most recent call last):
            ...
            TypeError: Symbolic function foo takes exactly 2 arguments (1 given)
        """
        return self._nargs

    def variables(self):
        """
        Returns the variables (of which there are none) present in
        this SFunction.

        EXAMPLES::

            sage: sqrt.variables()
            ()
        """
        return ()

    def default_variable(self):
        """
        Returns a default variable.

        EXAMPLES::

            sage: sin.default_variable()
            x
        """
        return SR.var('x')

    def serial(self):
        """
        Returns the "serial" of this symbolic function.  This is the
        integer used to identify this function in the Pynac C++
        library.  It is primarily useful lower level interactions with
        the library.

        EXAMPLES::

            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: foo.serial()               #random
            40
        """
        return int(self._serial)

    def __getattr__(self, attr):
        """
        This method allows us to access attributes set by
        :meth:`__setattr__`.

        EXAMPLES::

            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: foo.bar = 4
            sage: foo.bar
            4
        """
        try:
            return self.__dict__[attr]
        except KeyError:
            raise AttributeError, attr

    def __setattr__(self, attr, value):
        """
        This method allows us to store arbitrary Python attributes
        on symbolic functions which is normally not possible with
        Cython extension types.

        EXAMPLES::

            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: foo.bar = 4
            sage: foo.bar
            4
        """
        self.__dict__[attr] = value

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2); foo
            foo
        """
        return self.name()

    def _maxima_init_(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction('foo', 2)
            sage: print foo._maxima_init_()
            'foo
        """
        return "'%s"%self.name()

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
            return cmp(self._serial, (<SFunction>other)._serial)
        return False

    def __call__(self, *args, coerce=True):
        """
        EXAMPLES::

            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: x,y,z = var("x y z")
            sage: foo(x,y)
            foo(x, y)

            sage: foo(y)
            Traceback (most recent call last):
            ...
            TypeError: Symbolic function foo takes exactly 2 arguments (1 given)

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
            TypeError: cannot coerce arguments: ...

        """
        cdef GExVector vec
        if self._nargs > 0 and len(args) != self._nargs:
            raise TypeError, "Symbolic function %s takes exactly %s arguments (%s given)"%(self._name, self._nargs, len(args))

        if coerce:
            try:
                args = map(SR._coerce_, args)
            except TypeError, err:
                raise TypeError, "cannot coerce arguments: %s"%(err)

        if self._nargs == 0:
            for i from 0 <= i < len(args):
                vec.push_back((<Expression>args[i])._gobj)
            return new_Expression_from_GEx(SR, g_function_evalv(self._serial, vec))

        elif self._nargs == 1:
            return new_Expression_from_GEx(SR, g_function_eval1(self._serial,
                    (<Expression>args[0])._gobj))
        elif self._nargs == 2:
            return new_Expression_from_GEx(SR, g_function_eval2(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj))
        elif self._nargs == 3:
            return new_Expression_from_GEx(SR, g_function_eval3(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj))
        elif self._nargs == 4:
            return new_Expression_from_GEx(SR, g_function_eval4(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj))
        elif self._nargs == 5:
            return new_Expression_from_GEx(SR, g_function_eval5(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj))
        elif self._nargs == 6:
            return new_Expression_from_GEx(SR, g_function_eval6(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj))
        elif self._nargs == 7:
            return new_Expression_from_GEx(SR, g_function_eval7(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj))
        elif self._nargs == 8:
            return new_Expression_from_GEx(SR, g_function_eval8(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj))
        elif self._nargs == 9:
            return new_Expression_from_GEx(SR, g_function_eval9(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj))
        elif self._nargs == 10:
            return new_Expression_from_GEx(SR, g_function_eval10(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj))
        elif self._nargs == 11:
            return new_Expression_from_GEx(SR, g_function_eval11(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj))
        elif self._nargs == 12:
            return new_Expression_from_GEx(SR, g_function_eval12(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj, (<Expression>args[11])._gobj))
        elif self._nargs == 13:
            return new_Expression_from_GEx(SR, g_function_eval13(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj, (<Expression>args[11])._gobj,
                    (<Expression>args[12])._gobj))
        elif self._nargs == 14:
            return new_Expression_from_GEx(SR, g_function_eval14(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj, (<Expression>args[11])._gobj,
                    (<Expression>args[12])._gobj, (<Expression>args[13])._gobj))

    def _fast_float_(self, *vars):
        """
        Returns an object which provides fast floating point evaluation of
        self.

        See sage.ext.fast_eval? for more information.

        EXAMPLES::

            sage: sin._fast_float_()
            <sage.ext.fast_eval.FastDoubleFunc object at 0x...>
            sage: sin._fast_float_()(0)
            0.0

        ::

            sage: ff = cos._fast_float_(); ff
            <sage.ext.fast_eval.FastDoubleFunc object at 0x...>
            sage: ff.is_pure_c()
            True
            sage: ff(0)
            1.0

        ::

            sage: ff = erf._fast_float_()
            sage: ff.is_pure_c()
            False
            sage: ff(1.5)
            0.96610514647531076
            sage: erf(1.5)
            0.966105146475311
        """
        import sage.ext.fast_eval as fast_float

        args = [fast_float.fast_float_arg(n) for n in range(self.number_of_arguments())]
        try:
            return self(*args)
        except TypeError:
            return fast_float.fast_float_func(self, *args)

    def _fast_callable_(self, etb):
        r"""
        Given an ExpressionTreeBuilder, return an Expression representing
        this value.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x','y'])
            sage: sin._fast_callable_(etb)
            sin(v_0)
            sage: erf._fast_callable_(etb)
            {erf}(v_0)
        """
        args = [etb._var_number(n) for n in range(self.number_of_arguments())]
        return etb.call(self, *args)



function = SFunction

cdef class PrimitiveFunction(SFunction):
    def __init__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: from sage.symbolic.function import PrimitiveFunction
            sage: f = PrimitiveFunction('f', latex='f')
            sage: f
            f
            sage: loads(dumps(f))
            f
        """
        approx = kwds.pop('approx', None)
        conversions = kwds.pop('conversions', {})
        if 'ginac' in conversions:
            kwds['ginac_name'] = conversions['ginac']
        kwds['built_in_function'] = True
        kwds['latex_name'] = kwds.pop('latex', self.name())
        if 'nargs' not in kwds:
            kwds['nargs'] = 1

        for name in sfunctions_funcs:
            if hasattr(self, "_%s_"%name):
                kwds['%s_func'%name] = getattr(self, "_%s_"%name)

        SFunction.__init__(self, *args, **kwds)

        self._conversions = conversions
        self._approx_ = approx if approx is not None else self._generic_approx_

        from sage.symbolic.pynac import register_symbol
        register_symbol(self, self._conversions)

    def __call__(self, x, hold=False):
        """
        Evaluates this function at *x*.  First, it checks to see if *x*
        is a float in which case it calls :meth:`_approx_`.  Next,
        it checks to see if *x* has an attribute with the same name
        as.  If both of those fail, then it returns the symbolic version
        provided by :class:`SFunction`.

        EXAMPLES::

            sage: from sage.symbolic.function import PrimitiveFunction
            sage: s = PrimitiveFunction('sin'); s
            sin
            sage: s(float(1.0))
            0.8414709848078965
            sage: s(1.0)
            0.841470984807897
            sage: s(1)
            sin(1)
        """
        if isinstance(x, float):
            return self._approx_(x)

        try:
            obj = x.pyobject()
            return getattr(obj, self.name())()
        except (AttributeError, TypeError):
            pass

        if not hold:
            try:
                return getattr(x, self.name())()
            except AttributeError:
                pass
        return SFunction.__call__(self, x)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.symbolic.function import PrimitiveFunction
            sage: s = PrimitiveFunction('sin'); sin
            sin
            sage: latex(s)
            sin
            sage: s = PrimitiveFunction('sin', latex=r'\sin')
            sage: latex(s)
            \sin
        """
        if self._latex_name is not None:
            return self._latex_name
        else:
            return self.name()

    def _generic_approx_(self, x):
        """
        Returns the results of numerically evaluating this function at
        *x*.  This is a default implementation which tries to do the
        evaluation in Maxima.

        EXAMPLES::

            sage: from sage.symbolic.function import PrimitiveFunction
            sage: s = PrimitiveFunction('sin'); s
            sin
            sage: s._generic_approx_(0)
            0.0
        """
        from sage.rings.all import RR
        from sage.calculus.calculus import maxima
        try:
            return float(self(RR(x)))
        except TypeError:
            pass

        s = '%s(%s), numer'%(self._maxima_init_(), float(x))
        return float(maxima.eval(s))

    def _complex_approx_(self, x): # must be called with Python complex float as iput
        """
        Given a Python complex `x`, evaluate self and return a complex value.

        EXAMPLES::

            sage: complex(cos(3*I))
            (10.067661995777767+0j)

        The following fails because we and Maxima haven't implemented
        erf yet for complex values::

            sage: complex(erf(3*I))
            Traceback (most recent call last):
            ...
            TypeError: unable to simplify to complex approximation
        """
        from sage.calculus.calculus import maxima
        if x.imag == 0:
            return complex(self._approx_(x.real))
        s = '%s(%s+%s*%%i), numer'%(self._maxima_init_(), x.real, x.imag)
        a = maxima.eval(s).replace('%i', '1j')
        if '(' in a:
            # unable to simplify to a complex -- still function calls there.
            raise TypeError, "unable to simplify to complex approximation"
        return complex(eval(a))

    def _interface_init_(self, I=None):
        """
        EXAMPLES::

             sage: sin._maxima_init_()
             'sin'
        """
        if I is None:
            return repr(self)
        s = self._conversions.get(I.name(), None)
        return s if s is not None else repr(self)

    def _mathematica_init_(self):
        """
        EXAMPLES::

             sage: sin._mathematica_init_()
             'Sin'
             sage: exp._mathematica_init_()
             'Exp'
             sage: (exp(x) + sin(x) + tan(x))._mathematica_init_()
             '(Exp[x])+(Sin[x])+(Tan[x])'
        """
        s = self._conversions.get('mathematica', None)
        return s if s is not None else repr(self).capitalize()

    def _maxima_init_(self, I=None):
        """
        EXAMPLES::

            sage: exp._maxima_init_()
            'exp'
            sage: from sage.symbolic.function import PrimitiveFunction
            sage: f = PrimitiveFunction('f', latex='f', conversions=dict(maxima='ff'))
            sage: f._maxima_init_()
            'ff'
        """
        s = self._conversions.get('maxima', None)
        if s is None:
            return repr(self)
        else:
            return s


def get_sfunction_from_serial(serial):
    """
    Returns an already created SFunction given the serial.  These are
    stored in the dictionary
    :obj:`sage.symbolic.function.sfunction_serial_dict`.

    EXAMPLES::

        sage: from sage.symbolic.function import get_sfunction_from_serial
        sage: f = function('f'); f
        f
        sage: s = f.serial(); s #random
        65
        sage: get_sfunction_from_serial(s)
        f
    """
    global sfunction_serial_dict
    return sfunction_serial_dict.get(serial)

def pickle_wrapper(f):
    """
    Returns a pickled version of the function f if f is not None;
    otherwise, it returns None.  This is a wrapper around
    :func:`pickle_function`.

    EXAMPLES::

        sage: from sage.symbolic.function import pickle_wrapper
        sage: def f(x): return x*x
        sage: pickle_wrapper(f)
        "csage...."
        sage: pickle_wrapper(None) is None
        True
    """
    if f is None:
        return None
    return pickle_function(f)

def unpickle_wrapper(p):
    """
    Returns a unpickled version of the function defined by *p* if *p*
    is not None; otherwise, it returns None.  This is a wrapper around
    :func:`unpickle_function`.

    EXAMPLES::

        sage: from sage.symbolic.function import pickle_wrapper, unpickle_wrapper
        sage: def f(x): return x*x
        sage: s = pickle_wrapper(f)
        sage: g = unpickle_wrapper(s)
        sage: g(2)
        4
        sage: unpickle_wrapper(None) is None
        True
    """
    if p is None:
        return None
    return unpickle_function(p)

def init_sfunction_map():
    """
    Initializes a list mapping GiNaC function serials to the equivalent Sage
    functions.

    EXAMPLES::

        sage: from sage.symbolic.function import init_sfunction_map
        sage: f_list = init_sfunction_map()
        sage: gamma in f_list
        True
        sage: exp in f_list
        True
        sage: x = var('x')
        sage: sin(x).operator()
        sin
        sage: gamma(x).operator()
        gamma
    """
    # colliding serials (variable number of arguments): psi2, G2

    f_list = [None]*get_ginac_serial()

#    f_list[step_serial] = # step function
#    f_list[csgn_serial] = # complex sign
#    f_list[conjugate_serial] = # complex conjugation
#    f_list[real_part_serial] = # real part
#    f_list[imag_part_serial] = # imaginary part
    from sage.all import sin, cos, tan, asin, acos, \
            atan, sinh, cosh, tanh, asinh, acosh, atanh, exp, log, gamma, \
            factorial, abs_symbolic, polylog, dilog
    from sage.functions.log import function_log

    f_list[abs_serial] = abs_symbolic
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
    f_list[log_serial] = function_log # natural logarithm
    f_list[Li2_serial] = dilog # dilogarithm
    f_list[Li_serial] = polylog # classical polylogarithm as well as multiple polylogarithm
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

    EXAMPLES::

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
