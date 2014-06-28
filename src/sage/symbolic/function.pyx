###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2008 - 2010 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

r"""

Support for symbolic functions.

"""
include "sage/ext/interrupt.pxi"
include "sage/ext/cdefs.pxi"

from functools import reduce
from ginac cimport *

from sage.structure.sage_object cimport SageObject
from expression cimport new_Expression_from_GEx, Expression
from ring import SR

from sage.structure.parent cimport Parent
from sage.structure.coerce import parent
from sage.structure.element import get_coercion_model

# we keep a database of symbolic functions initialized in a session
# this also makes the .operator() method of symbolic expressions work
cdef dict sfunction_serial_dict = {}

from sage.misc.fpickle import pickle_function, unpickle_function
from sage.ext.fast_eval import FastDoubleFunc

# List of functions which ginac allows us to define custom behavior for.
# Changing the order of this list could cause problems unpickling old pickles.
sfunctions_funcs = ['eval', 'evalf', 'conjugate', 'real_part', 'imag_part',
        'derivative', 'power', 'series', 'print', 'print_latex', 'tderivative']

cdef class Function(SageObject):
    """
    Base class for symbolic functions defined through Pynac in Sage.

    This is an abstract base class, with generic code for the interfaces
    and a :meth:`__call__` method. Subclasses should implement the
    :meth:`_is_registered` and :meth:`_register_function` methods.

    This class is not intended for direct use, instead use one of the
    subclasses :class:`BuiltinFunction` or :class:`SymbolicFunction`.
    """
    def __init__(self, name, nargs, latex_name=None, conversions=None,
            evalf_params_first=True, alt_name=None):
        """
        This is an abstract base class. It's not possible to test it directly.

        EXAMPLES::

            sage: f = function('f', nargs=1, conjugate_func=lambda self,x: 2r*x) # indirect doctest
            sage: f(2)
            f(2)
            sage: f(2).conjugate()
            4

        TESTS::

            # eval_func raises exception
            sage: def ef(self, x): raise RuntimeError, "foo"
            sage: bar = function("bar", nargs=1, eval_func=ef)
            sage: bar(x)
            Traceback (most recent call last):
            ...
            RuntimeError: foo

            # eval_func returns non coercible
            sage: def ef(self, x): return ZZ
            sage: bar = function("bar", nargs=1, eval_func=ef)
            sage: bar(x)
            Traceback (most recent call last):
            ...
            TypeError: function did not return a symbolic expression or an element that can be coerced into a symbolic expression

            # eval_func is not callable
            sage: bar = function("bar", nargs=1, eval_func=5)
            Traceback (most recent call last):
            ...
            ValueError: eval_func parameter must be callable
        """
        self._name = name
        self._alt_name = alt_name
        self._nargs = nargs
        self._latex_name = latex_name
        self._evalf_params_first = evalf_params_first
        self._conversions = {} if conversions is None else conversions

        # handle custom printing
        # if print_func is defined, it is used instead of name
        # latex printing can be customised either by setting a string latex_name
        # or giving a custom function argument print_latex_func
        if latex_name and hasattr(self, '_print_latex_'):
            raise ValueError, "only one of latex_name or _print_latex_ should be specified."

        # only one of derivative and tderivative should be defined
        if hasattr(self, '_derivative_') and hasattr(self, '_tderivative_'):
            raise ValueError, "only one of _derivative_ or _tderivative_ should be defined."

        for fname in sfunctions_funcs:
            real_fname = '_%s_'%fname
            if hasattr(self, real_fname) and not \
                    callable(getattr(self, real_fname)):
                raise ValueError,  real_fname + " parameter must be callable"

        if not self._is_registered():
            self._register_function()

            global sfunction_serial_dict
            sfunction_serial_dict[self._serial] = self

            from sage.symbolic.pynac import symbol_table, register_symbol
            symbol_table['functions'][self._name] = self

            register_symbol(self, self._conversions)

    cdef _is_registered(self):
        """
        Check if this function is already registered. If it is, set
        `self._serial` to the right value.
        """
        raise NotImplementedError, "this is an abstract base class, it shouldn't be initialized directly"

    cdef _register_function(self):
        """

        TESTS:

        After :trac:`9240`, pickling and unpickling of symbolic
        functions was broken. We check here that this is fixed
        (:trac:`11919`)::

            sage: f = function('f', x)
            sage: s = dumps(f)
            sage: loads(s)
            f(x)
            sage: deepcopy(f)
            f(x)

        """
        cdef GFunctionOpt opt
        opt = g_function_options_args(self._name, self._nargs)

        if hasattr(self, '_eval_'):
            opt.eval_func(self)

        if not self._evalf_params_first:
            opt.do_not_evalf_params()

        if hasattr(self, '_subs_'):
            opt.subs_func(self)

        if hasattr(self, '_evalf_'):
            opt.evalf_func(self)

        if hasattr(self, '_conjugate_'):
            opt.conjugate_func(self)

        if hasattr(self, '_real_part_'):
            opt.real_part_func(self)

        if hasattr(self, '_imag_part_'):
            opt.imag_part_func(self)

        if hasattr(self, '_derivative_'):
            opt.derivative_func(self)

        if hasattr(self, '_tderivative_'):
            opt.do_not_apply_chain_rule()
            opt.derivative_func(self)

        if hasattr(self, '_power_'):
            opt.power_func(self)

        if hasattr(self, '_series_'):
            opt.series_func(self)

        # custom print functions are called from python
        # so we don't register them with the ginac function_options object

        if self._latex_name:
            opt.latex_name(self._latex_name)

        self._serial = g_register_new(opt)
        g_foptions_assign(g_registered_functions().index(self._serial), opt)

    def _eval_default(self, *args):
        """
        Default automatic evaluation function.

        Calls numeric evaluation if the argument is not exact.

        TESTS::

            sage: coth(5)  # indirect doctest
            coth(5)
            sage: coth(0.5)
            2.16395341373865
            sage: from sage.symbolic.function import BuiltinFunction
            sage: class Test(BuiltinFunction):
            ....:     def __init__(self):
            ....:         BuiltinFunction.__init__(self, 'test', nargs=2)
            ....:     def _evalf_(self, x, y, parent):
            ....:         return x + 1
            ....:     def _eval_(self, x, y):
            ....:         res = self._eval_default(x, y)
            ....:         if res:
            ....:             return res
            ....:         elif x == 2:
            ....:             return 3
            ....:         else:
            ....:             return
            sage: test = Test()
            sage: test(1.3, 4)
            2.30000000000000
            sage: test(pi, 4)
            test(pi, 4)
            sage: test(2, x)
            3
            sage: test(2., 4)
            3.00000000000000
            sage: test(1 + 1.0*I, 2)
            2.00000000000000 + 1.00000000000000*I
            sage: class Test2(BuiltinFunction):
            ....:     def __init__(self):
            ....:         BuiltinFunction.__init__(self, 'test', nargs=1)
            ....:     def _evalf_(self, x, parent):
            ....:         return 0.5
            ....:     def _eval_(self, x):
            ....:         res = self._eval_default(x)
            ....:         if res:
            ....:             return res
            ....:         else:
            ....:             return 3
            sage: test2 = Test2()
            sage: test2(1.3)
            0.500000000000000
            sage: test2(pi)
            3
        """
        if len(args) == 1:
            x = args[0]
            try:
                method = getattr(x, self.name())
            except AttributeError:
                pass
            else:
                return method()
            if is_inexact(x) and not parent_c(x) is SR:
                return self._evalf_(x, parent=parent(x))
            return
        else:
            cc = get_coercion_model().canonical_coercion
            coerced = reduce(lambda x, y: cc(x, y)[0], args)
            if is_inexact(coerced) and not parent_c(coerced) is SR:
                return self._evalf_(*args, parent=parent(coerced))
            else:
                return

    def __hash__(self):
        """
        EXAMPLES::

            sage: f = function('f', nargs=1, conjugate_func=lambda self,x: 2r*x)
            sage: f.__hash__() #random
            -2224334885124003860
            sage: hash(f(2)) #random
            4168614485
        """
        return hash(self._name)*(self._nargs+1)*self._serial

    def __repr__(self):
        """
        EXAMPLES::

            sage: foo = function("foo", nargs=2); foo
            foo
        """
        return self._name

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.symbolic.function import SymbolicFunction
            sage: s = SymbolicFunction('foo'); s
            foo
            sage: latex(s)
            foo
            sage: s = SymbolicFunction('foo', latex_name=r'{\rm foo}')
            sage: latex(s)
            {\rm foo}
            sage: s._latex_()
            '{\\rm foo}'
        """
        if self._latex_name is not None:
            return self._latex_name
        else:
            return self._name

    def __cmp__(self, other):
        """
        TESTS:
            sage: foo = function("foo", nargs=2)
            sage: foo == foo
            True
            sage: foo == 2
            False
            sage: foo(1,2).operator() == foo
            True

        """
        if PY_TYPE_CHECK(other, Function):
            return cmp(self._serial, (<Function>other)._serial)
        return False

    def __call__(self, *args, bint coerce=True, bint hold=False):
        """
        Evaluates this function at the given arguments.

        We coerce the arguments into symbolic expressions if coerce=True, then
        call the Pynac evaluation method, which in turn passes the arguments to
        a custom automatic evaluation method if ``_eval_()`` is defined.

        EXAMPLES::

            sage: foo = function("foo", nargs=2)
            sage: x,y,z = var("x y z")
            sage: foo(x,y)
            foo(x, y)

            sage: foo(y)
            Traceback (most recent call last):
            ...
            TypeError: Symbolic function foo takes exactly 2 arguments (1 given)

            sage: bar = function("bar")
            sage: bar(x)
            bar(x)
            sage: bar(x,y)
            bar(x, y)

        The `hold` argument prevents automatic evaluation of the function::

            sage: exp(log(x))
            x
            sage: exp(log(x), hold=True)
            e^log(x)

        We can also handle numpy types::

            sage: import numpy
            sage: sin(numpy.arange(5))
            array([ 0.        ,  0.84147098,  0.90929743,  0.14112001, -0.7568025 ])

        Symbolic functions evaluate non-exact input numerically, and return
        symbolic expressions on exact input::

            sage: arctan(1)
            1/4*pi
            sage: arctan(float(1))
            0.7853981633974483

        Precision of the result depends on the precision of the input::

            sage: arctan(RR(1))
            0.785398163397448
            sage: arctan(RealField(100)(1))
            0.78539816339744830961566084582

        Return types for non-exact input depends on the input type::

            sage: type(exp(float(0)))
            <type 'float'>
            sage: exp(RR(0)).parent()
            Real Field with 53 bits of precision


        TESTS:

        Test coercion::

            sage: bar(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce arguments: ...
            sage: exp(QQbar(I))
            0.540302305868140 + 0.841470984807897*I

        For functions with single argument, if coercion fails we try to call
        a method with the name of the function on the object::

            sage: M = matrix(SR, 2, 2, [x, 0, 0, I*pi])
            sage: exp(M)
            [e^x   0]
            [  0  -1]

        Make sure we can pass mpmath arguments (:trac:`13608`)::

            sage: import mpmath
            sage: with mpmath.workprec(128): sin(mpmath.mpc('0.5', '1.2'))
            mpc(real='0.86807452059118713192871150787046523179886', imag='1.3246769633571289324095313649562791720086')

        """
        if self._nargs > 0 and len(args) != self._nargs:
            raise TypeError, "Symbolic function %s takes exactly %s arguments (%s given)"%(self._name, self._nargs, len(args))

        # support fast_float
        if self._nargs == 1:
            if isinstance(args[0], FastDoubleFunc):
                try:
                    method = getattr(args[0], self._name)
                except AttributeError, err:
                    raise TypeError, "cannot handle fast float arguments"
                else:
                    return method()

        # support numpy arrays as arguments
        if any([type(arg).__module__ == 'numpy' for arg in args]): # avoid importing
            import numpy
            # check that at least one of the arguments is a numpy array
            if any([isinstance(arg, numpy.ndarray) for arg in args]):
                try:
                    modulefn = getattr(numpy, self.name())
                except AttributeError:
                    return self._eval_numpy_(*args)
                else:
                    return modulefn(*args)

        # support mpmath mpf and mpc numbers as arguments
        if any(['mpmath' in type(arg).__module__ for arg in args]): # avoid importing
            import mpmath
            # check that at least one of the arguments is an mpmath type
            if any([isinstance(arg, (mpmath.mpf, mpmath.mpc)) for arg in args]):
                try:
                    modulefn = getattr(mpmath, self.name())
                except AttributeError:
                    return self._eval_mpmath_(*args)
                else:
                    return modulefn(*args)

        # if the given input is a symbolic expression, we don't convert it back
        # to a numeric type at the end
        if len(args) == 1 and parent_c(args[0]) is SR:
            symbolic_input = True
        else:
            symbolic_input = False


        cdef Py_ssize_t i
        if coerce:
            try:
                args = map(SR.coerce, args)
            except TypeError as err:
                # If the function takes only one argument, we try to call
                # a method with the name of this function on the object.
                # This makes the following work:
                #     sage: M = matrix(SR, 2, 2, [x, 0, 0, I*pi])
                #     sage: exp(M)
                #     [e^x   0]
                #     [  0  -1]
                if len(args) == 1:
                    method = getattr(args[0], self._name, None)
                    if callable(method):
                        return method()

                # There is no natural coercion from QQbar to the symbolic ring
                # in order to support
                #     sage: QQbar(sqrt(2)) + sqrt(3)
                #     3.146264369941973?
                # to work around this limitation, we manually convert
                # elements of QQbar to symbolic expressions here
                from sage.rings.qqbar import QQbar, AA
                nargs = [None]*len(args)
                for i in range(len(args)):
                    carg = args[i]
                    if PY_TYPE_CHECK(carg, Element) and \
                            (<Element>carg)._parent is QQbar or \
                            (<Element>carg)._parent is AA:
                        nargs[i] = SR(carg)
                    else:
                        try:
                            nargs[i] = SR.coerce(carg)
                        except Exception:
                            raise TypeError, "cannot coerce arguments: %s"%(err)
                args = nargs
        else: # coerce == False
            for a in args:
                if not PY_TYPE_CHECK(a, Expression):
                    raise TypeError, "arguments must be symbolic expressions"

        cdef GEx res
        cdef GExVector vec
        if self._nargs == 0 or self._nargs > 3:
            for i from 0 <= i < len(args):
                vec.push_back((<Expression>args[i])._gobj)
            res = g_function_evalv(self._serial, vec, hold)
        elif self._nargs == 1:
            res = g_function_eval1(self._serial,
                    (<Expression>args[0])._gobj, hold)
        elif self._nargs == 2:
            res = g_function_eval2(self._serial, (<Expression>args[0])._gobj,
                    (<Expression>args[1])._gobj, hold)
        elif self._nargs == 3:
            res = g_function_eval3(self._serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, hold)


        if not symbolic_input and is_a_numeric(res):
            return py_object_from_numeric(res)

        return new_Expression_from_GEx(SR, res)

    def name(self):
        """
        Returns the name of this function.

        EXAMPLES::

            sage: foo = function("foo", nargs=2)
            sage: foo.name()
            'foo'
        """
        return self._name

    def number_of_arguments(self):
        """
        Returns the number of arguments that this function takes.

        EXAMPLES::

            sage: foo = function("foo", nargs=2)
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

            sage: sin.variables()
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

    def _interface_init_(self, I=None):
        """
        EXAMPLES::

             sage: sin._interface_init_(maxima)
             'sin'
        """
        if I is None:
            return self._name
        return self._conversions.get(I.name(), self._name)

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

    def _sympy_init_(self, I=None):
        """
        EXAMPLES::

            sage: arcsin._sympy_init_()
            'asin'
            sage: from sage.symbolic.function import SymbolicFunction
            sage: g = SymbolicFunction('g', conversions=dict(sympy='gg'))
            sage: g._sympy_init_()
            'gg'
            sage: g(x)._sympy_()
            Traceback (most recent call last):
            ...
            NotImplementedError: SymPy function 'gg' doesn't exist
        """
        return self._conversions.get('sympy', self._name)

    def _maxima_init_(self, I=None):
        """
        EXAMPLES::

            sage: exp._maxima_init_()
            'exp'
            sage: from sage.symbolic.function import SymbolicFunction
            sage: f = SymbolicFunction('f', latex_name='f', conversions=dict(maxima='ff'))
            sage: f._maxima_init_()
            'ff'
        """
        return self._conversions.get('maxima', self._name)

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
            0.9661051464753108
            sage: erf(1.5)
            0.966105146475311
        """
        import sage.ext.fast_eval as fast_float

        args = [fast_float.fast_float_arg(n) for n in range(self.number_of_arguments())]
        try:
            return self(*args)
        except TypeError as err:
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

    def _eval_numpy_(self, *args):
        r"""
        Evaluates this function at the given arguments.

        At least one of elements of args is supposed to be a numpy array.

        EXAMPLES::

            sage: import numpy
            sage: a = numpy.arange(5)
            sage: csc(a)
            doctest:270: RuntimeWarning: divide by zero encountered in divide
            array([        inf,  1.18839511,  1.09975017,  7.0861674 , -1.32134871])

            sage: factorial(a)
            Traceback (most recent call last):
            ...
            NotImplementedError: The Function factorial does not support numpy arrays as arguments
        """
        raise NotImplementedError("The Function %s does not support numpy arrays as arguments" % self.name())

    def _eval_mpmath_(self, *args):
        r"""
        Evaluates this function for arguments of mpmath types.

        The default implementation casts its arguments to sage reals
        of the appropriate precision.

        EXAMPLES::

        At the time of this writing, mpmath had no arcsin, only asin.
        So the following call would actually fall back to the default
        implementation, using sage reals instead of mpmath ones. This
        might change when aliases for these functions are established.

            sage: import mpmath
            sage: with mpmath.workprec(128): arcsin(mpmath.mpf('0.5'))
            mpf('0.52359877559829887307710723054658381403157')

        TESTS:

        To ensure that we actually can fall back to an implementation
        not using mpmath, we have to create a custom function which
        will certainly never get created in mpmath.

            sage: import mpmath
            sage: from sage.symbolic.function import BuiltinFunction
            sage: class NoMpmathFn(BuiltinFunction):
            ....:         def _eval_(self, arg):
            ....:                 parent = arg.parent()
            ....:                 prec = parent.prec()
            ....:                 assert parent == RealField(prec)
            ....:                 return prec
            sage: noMpmathFn = NoMpmathFn("noMpmathFn")
            sage: with mpmath.workprec(64): noMpmathFn(sqrt(mpmath.mpf('2')))
            mpf('64.0')
            sage: mpmath.noMpmathFn = lambda x: 123
            sage: with mpmath.workprec(64): noMpmathFn(sqrt(mpmath.mpf('2')))
            123
            sage: del mpmath.noMpmathFn

        """
        import mpmath
        from sage.libs.mpmath.utils import mpmath_to_sage, sage_to_mpmath
        prec = mpmath.mp.prec
        args = [mpmath_to_sage(x, prec)
                if isinstance(x, (mpmath.mpf, mpmath.mpc)) else x
                for x in args]
        res = self(*args)
        res = sage_to_mpmath(res, prec)
        return res

cdef class GinacFunction(BuiltinFunction):
    """
    This class provides a wrapper around symbolic functions already defined in
    Pynac/GiNaC.

    GiNaC provides custom methods for these functions defined at the C++ level.
    It is still possible to define new custom functionality or override those
    already defined.

    There is also no need to register these functions.
    """
    def __init__(self, name, nargs=1, latex_name=None, conversions=None,
            ginac_name=None, evalf_params_first=True):
        """
        TESTS::

            sage: from sage.functions.trig import Function_sin
            sage: s = Function_sin() # indirect doctest
            sage: s(0)
            0
            sage: s(pi)
            0
            sage: s(pi/2)
            1
        """
        self._ginac_name = ginac_name
        BuiltinFunction.__init__(self, name, nargs, latex_name, conversions,
                evalf_params_first=evalf_params_first)

    cdef _is_registered(self):
        # Since this is function is defined in C++, it is already in
        # ginac's function registry
        fname = self._ginac_name if self._ginac_name is not None else self._name
        # get serial
        try:
            self._serial = find_function(fname, self._nargs)
        except ValueError as err:
            raise ValueError, "cannot find GiNaC function with name %s and %s arguments"%(fname, self._nargs)

        global sfunction_serial_dict
        return self._serial in sfunction_serial_dict

    cdef _register_function(self):
        # We don't need to add anything to GiNaC's function registry
        # However, if any custom methods were provided in the python class,
        # we should set the properties of the function_options object
        # corresponding to this function
        cdef GFunctionOpt opt = g_registered_functions().index(self._serial)

        if hasattr(self, '_eval_'):
            opt.eval_func(self)

        if not self._evalf_params_first:
            opt.do_not_evalf_params()

        if hasattr(self, '_evalf_'):
            opt.evalf_func(self)

        if hasattr(self, '_conjugate_'):
            opt.conjugate_func(self)

        if hasattr(self, '_real_part_'):
            opt.real_part_func(self)

        if hasattr(self, '_imag_part_'):
            opt.imag_part_func(self)

        if hasattr(self, '_derivative_'):
            opt.derivative_func(self)

        if hasattr(self, '_tderivative_'):
            opt.do_not_apply_chain_rule()
            opt.derivative_func(self)

        if hasattr(self, '_power_'):
            opt.power_func(self)

        if hasattr(self, '_series_'):
            opt.series_func(self)

        # overriding print functions is not supported

        if self._latex_name:
            opt.latex_name(self._latex_name)

        g_foptions_assign(g_registered_functions().index(self._serial), opt)


cdef class BuiltinFunction(Function):
    """
    This is the base class for symbolic functions defined in Sage.

    If a function is provided by the Sage library, we don't need to pickle
    the custom methods, since we can just initialize the same library function
    again. This allows us to use Cython for custom methods.

    We assume that each subclass of this class will define one symbolic
    function. Make sure you use subclasses and not just call the initializer
    of this class.
    """
    def __init__(self, name, nargs=1, latex_name=None, conversions=None,
            evalf_params_first=True, alt_name=None):
        """
        TESTS::

            sage: from sage.functions.trig import Function_cot
            sage: c = Function_cot() # indirect doctest
            sage: c(pi/2)
            0
        """
        Function.__init__(self, name, nargs, latex_name, conversions,
                evalf_params_first, alt_name = alt_name)

    def __call__(self, *args, bint coerce=True, bint hold=False,
            bint dont_call_method_on_arg=False):
        r"""
        Evaluate this function on the given arguments and return the result.

        EXAMPLES::

            sage: exp(5)
            e^5
            sage: gamma(15)
            87178291200

        TESTS::

            sage: from sage.symbolic.function import BuiltinFunction
            sage: class A:
            ....:     def foo(self):
            ....:         return 'foo'
            sage: foo = BuiltinFunction(name='foo')
            sage: foo(A())
            'foo'
            sage: bar = BuiltinFunction(name='bar', alt_name='foo')
            sage: bar(A())
            'foo'
        """
        # if there is only one argument, and the argument has an attribute
        # with the same name as this function, try to call it to get the result
        # The argument dont_call_method_on_arg is used to prevent infinite loops
        # when .exp(), .log(), etc. methods call this symbolic function on
        # themselves
        if len(args) == 1 and not hold and not dont_call_method_on_arg:
            arg = args[0]
            method = getattr(arg, self._name, None)
            if callable(method):
                return method()
            elif self._alt_name is not None:
                method = getattr(arg, self._alt_name, None)
                if method is not None:
                    return method()

        res = super(BuiltinFunction, self).__call__(
                        *args, coerce=coerce, hold=hold)

        # we want to convert the result to the original parent if the input
        # is not exact, so we store the parent here
        org_parent = parent_c(args[0])

        # convert the result back to the original parent previously stored
        # otherwise we end up with
        #     sage: arctan(RR(1))
        #     1/4*pi
        # which is surprising, to say the least...
        if org_parent is not SR and \
                (org_parent is float or org_parent is complex or \
                (PY_TYPE_CHECK(org_parent, Parent) and \
                    not org_parent.is_exact())):
            try:
                return org_parent(res)
            except (TypeError, ValueError):
                pass

            # conversion to the original parent failed
            # we try if it works with the corresponding complex domain
            if org_parent is float:
                try:
                    return complex(res)
                except (TypeError, ValueError):
                    pass
            elif hasattr(org_parent, 'complex_field'):
                try:
                    return org_parent.complex_field()(res)
                except (TypeError, ValueError):
                    pass

        return res

    cdef _is_registered(self):
        """
        TESTS:

        Check if :trac:`13586` is fixed::

            sage: from sage.symbolic.function import BuiltinFunction
            sage: class AFunction(BuiltinFunction):
            ....:       def __init__(self, name, exp=1):
            ....:           self.exponent=exp
            ....:           BuiltinFunction.__init__(self, name, nargs=1)
            ....:       def _eval_(self, arg):
            ....:               return arg**self.exponent
            sage: p2 = AFunction('p2', 2)
            sage: p2(x)
            x^2
            sage: p3 = AFunction('p3', 3)
            sage: p3(x)
            x^3
            sage: loads(dumps(cot)) == cot    # :trac:`15138`
            True
        """
        # check if already defined
        cdef int serial = -1

        # search ginac registry for name and nargs
        try:
            serial = find_function(self._name, self._nargs)
        except ValueError as err:
            pass

        # if match, get operator from function table
        global sfunction_serial_dict
        if serial != -1 and serial in sfunction_serial_dict and \
                sfunction_serial_dict[serial].__class__ == self.__class__:
                    # if the returned function is of the same type
                    self._serial = serial
                    return True

        return False

    def __reduce__(self):
        """
        EXAMPLES::

            sage: cot.__reduce__()
            (<class 'sage.functions.trig.Function_cot'>, ())

            sage: f = loads(dumps(cot)) #indirect doctest
            sage: f(pi/2)
            0
        """
        return self.__class__, tuple()

    # this is required to read old pickles of erf, elliptic_ec, etc.
    def __setstate__(self, state):
        """
        EXAMPLES::

            sage: cot.__setstate__([1,0])
            Traceback (most recent call last):
            ...
            ValueError: cannot read pickle
            sage: cot.__setstate__([0]) #don't try this at home
        """
        if state[0] == 0:
            # old pickle data
            # we call __init__ since Python only allocates the class and does
            # not call __init__ before passing the pickled state to __setstate__
            self.__init__()
        else:
            # we should never end up here
            raise ValueError, "cannot read pickle"


cdef class SymbolicFunction(Function):
    """
    This is the basis for user defined symbolic functions. We try to pickle or
    hash the custom methods, so subclasses must be defined in Python not Cython.
    """
    def __init__(self, name, nargs=0, latex_name=None, conversions=None,
            evalf_params_first=True):
        """
        EXAMPLES::

            sage: from sage.symbolic.function import SymbolicFunction
            sage: class my_function(SymbolicFunction):
            ....:     def __init__(self):
            ....:         SymbolicFunction.__init__(self, 'foo', nargs=2)
            ....:     def _evalf_(self, x, y, parent=None, algorithm=None):
            ....:         return x*y*2r
            ....:     def _conjugate_(self, x, y):
            ....:         return x
            sage: foo = my_function()
            sage: foo
            foo
            sage: foo(2,3)
            foo(2, 3)
            sage: foo(2,3).n()
            12.0000000000000
            sage: foo(2,3).conjugate()
            2
        """
        self.__hinit = False
        Function.__init__(self, name, nargs, latex_name, conversions,
                evalf_params_first)


    cdef _is_registered(SymbolicFunction self):
        # see if there is already an SFunction with the same state
        cdef Function sfunc
        cdef long myhash = self._hash_()
        for sfunc in sfunction_serial_dict.itervalues():
            if PY_TYPE_CHECK(sfunc, SymbolicFunction) and \
                    myhash == (<SymbolicFunction>sfunc)._hash_():
                # found one, set self._serial to be a copy
                self._serial = sfunc._serial
                return True

        return False

    # cache the hash value of this function
    # this is used very often while unpickling to see if there is already
    # a function with the same properties
    cdef long _hash_(self) except -1:
        if not self.__hinit:
            # create a string representation of this SFunction
            slist = [self._nargs, self._name, str(self._latex_name),
                    self._evalf_params_first]
            for fname in sfunctions_funcs:
                real_fname = '_%s_'%fname
                if hasattr(self, '%s'%real_fname):
                    slist.append(hash(getattr(self, real_fname).__code__))
                else:
                    slist.append(' ')
            self.__hcache = hash(tuple(slist))
            self.__hinit = True
        return self.__hcache

    def __hash__(self):
        """
        TESTS::

            sage: foo = function("foo", nargs=2)
            sage: hash(foo)      # random output
            -6859868030555295348

            sage: def ev(self, x): return 2*x
            sage: foo = function("foo", nargs=2, eval_func = ev)
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

            sage: foo = function("foo", nargs=2)
            sage: foo.__getstate__()
            (2, 'foo', 2, None, {}, True, [None, None, None, None, None, None, None, None, None, None, None])
            sage: t = loads(dumps(foo))
            sage: t == foo
            True
            sage: var('x,y')
            (x, y)
            sage: t(x,y)
            foo(x, y)

            sage: def ev(self, x,y): return 2*x
            sage: foo = function("foo", nargs=2, eval_func = ev)
            sage: foo.__getstate__()
            (2, 'foo', 2, None, {}, True, ["...", None, None, None, None, None, None, None, None, None, None])

            sage: u = loads(dumps(foo))
            sage: u == foo
            True
            sage: t == u
            False
            sage: u(y,x)
            2*y

            sage: def evalf_f(self, x, **kwds): return int(6)
            sage: foo = function("foo", nargs=1, evalf_func=evalf_f)
            sage: foo.__getstate__()
            (2, 'foo', 1, None, {}, True, [None, "...", None, None, None, None, None, None, None, None, None])

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
            foo(x)^2 + x^y + foo(y)
            sage: u.subs(y=0)
            foo(x)^2 + foo(0) + 1
            sage: u.subs(y=0).n()
            43.0000000000000
        """
        return (2, self._name, self._nargs, self._latex_name, self._conversions,
                self._evalf_params_first,
                map(pickle_wrapper, [getattr(self, '_%s_'%fname) \
                        if hasattr(self, '_%s_'%fname) else None \
                        for fname in sfunctions_funcs]))

    def __setstate__(self, state):
        """
        Initializes the state of the object from data saved in a pickle.

        During unpickling __init__ methods of classes are not called, the saved
        data is passed to the class via this function instead.

        TESTS::

            sage: var('x,y')
            (x, y)
            sage: foo = function("foo", nargs=2)
            sage: bar = function("bar", nargs=1)
            sage: bar.__setstate__(foo.__getstate__())

        ::

            sage: g = function('g', nargs=1, conjugate_func=lambda y,x: 2*x)
            sage: st = g.__getstate__()
            sage: f = function('f')
            sage: f(x)
            f(x)
            sage: f(x).conjugate() # no special conjugate method
            conjugate(f(x))
            sage: f.__setstate__(st)
            sage: f(x+1).conjugate() # now there is a special method
            2*x + 2

        Note that the other direction doesn't work here, since foo._hash_()
        hash already been initialized.::

            sage: bar
            foo
            sage: bar(x,y)
            foo(x, y)
        """
        # check input
        if not ((state[0] == 1 and len(state) == 6) or \
                (state[0] == 2 and len(state) == 7)):
            raise ValueError, "unknown state information"

        name = state[1]
        nargs = state[2]
        latex_name = state[3]
        conversions = state[4]

        if state[0] == 1:
            evalf_params_first = True
            function_pickles = state[5]
        elif state[0] == 2:
            evalf_params_first = state[5]
            function_pickles = state[6]

        for pickle, fname in zip(function_pickles, sfunctions_funcs):
            if pickle:
                real_fname = '_%s_'%fname
                setattr(self, real_fname, unpickle_function(pickle))

        SymbolicFunction.__init__(self, name, nargs, latex_name,
                conversions, evalf_params_first)


cdef class DeprecatedSFunction(SymbolicFunction):
    cdef dict __dict__
    def __init__(self, name, nargs=0, latex_name=None):
        """
        EXAMPLES::

            sage: from sage.symbolic.function import DeprecatedSFunction
            sage: foo = DeprecatedSFunction("foo", 2)
            sage: foo
            foo
            sage: foo(x,2)
            foo(x, 2)
            sage: foo(2)
            Traceback (most recent call last):
            ...
            TypeError: Symbolic function foo takes exactly 2 arguments (1 given)
        """
        self.__dict__ = {}
        SymbolicFunction.__init__(self, name, nargs, latex_name)

    def __getattr__(self, attr):
        """
        This method allows us to access attributes set by
        :meth:`__setattr__`.

        EXAMPLES::

            sage: from sage.symbolic.function import DeprecatedSFunction
            sage: foo = DeprecatedSFunction("foo", 2)
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

            sage: from sage.symbolic.function import DeprecatedSFunction
            sage: foo = DeprecatedSFunction("foo", 2)
            sage: foo.bar = 4
            sage: foo.bar
            4
        """
        self.__dict__[attr] = value

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.function import DeprecatedSFunction
            sage: foo = DeprecatedSFunction("foo", 2)
            sage: foo.__reduce__()
            (<function unpickle_function at ...>, ('foo', 2, None, {}, True, [None, None, None, None, None, None, None, None, None, None, None]))
        """
        from sage.symbolic.function_factory import unpickle_function
        state = self.__getstate__()
        name = state[1]
        nargs = state[2]
        latex_name = state[3]
        conversions = state[4]
        evalf_params_first = state[5]
        pickled_functions = state[6]
        return (unpickle_function, (name, nargs, latex_name, conversions,
            evalf_params_first, pickled_functions))

    def __setstate__(self, state):
        """
        EXAMPLES::

            sage: from sage.symbolic.function import DeprecatedSFunction
            sage: foo = DeprecatedSFunction("foo", 2)
            sage: foo.__setstate__([0, 'bar', 1, '\\bar', [None]*10])
            sage: foo
            bar
            sage: foo(x)
            bar(x)
            sage: latex(foo(x))
            \bar\left(x\right)
        """
        name = state[1]
        nargs = state[2]
        latex_name = state[3]
        self.__dict__ = {}
        for pickle, fname in zip(state[4], sfunctions_funcs):
            if pickle:
                if fname == 'evalf':
                    from sage.symbolic.function_factory import \
                            deprecated_custom_evalf_wrapper
                    setattr(self, '_evalf_',
                            deprecated_custom_evalf_wrapper(
                                unpickle_function(pickle)))
                    continue
                real_fname = '_%s_'%fname
                setattr(self, real_fname, unpickle_function(pickle))

        SymbolicFunction.__init__(self, name, nargs, latex_name, None)

SFunction = DeprecatedSFunction
PrimitiveFunction = DeprecatedSFunction


def get_sfunction_from_serial(serial):
    """
    Returns an already created SFunction given the serial.  These are
    stored in the dictionary
    `sage.symbolic.function.sfunction_serial_dict`.

    EXAMPLES::

        sage: from sage.symbolic.function import get_sfunction_from_serial
        sage: get_sfunction_from_serial(65) #random
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

def is_inexact(x):
    """
    Returns True if the argument is an inexact object.

    TESTS::

        sage: from sage.symbolic.function import is_inexact
        sage: is_inexact(5)
        False
        sage: is_inexact(5.)
        True
        sage: is_inexact(pi)
        True
        sage: is_inexact(5r)
        False
        sage: is_inexact(5.4r)
        True
    """
    if isinstance(x, (float, complex)):
        return True
    if isinstance(x, Element):
        return not (<Element>x)._parent.is_exact()
    return False
