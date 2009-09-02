###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################
"""
Symbolic Expressions

RELATIONAL EXPRESSIONS:

We create a relational expression::

    sage: x = var('x')
    sage: eqn = (x-1)^2 <= x^2 - 2*x + 3
    sage: eqn.subs(x == 5)
    16 <= 18

Notice that squaring the relation squares both sides.

::

    sage: eqn^2
    (x - 1)^4 <= (x^2 - 2*x + 3)^2
    sage: eqn.expand()
    x^2 - 2*x + 1 <= x^2 - 2*x + 3

The can transform a true relational into a false one::

    sage: eqn = SR(-5) < SR(-3); eqn
    -5 < -3
    sage: bool(eqn)
    True
    sage: eqn^2
    25 < 9
    sage: bool(eqn^2)
    False

We can do arithmetic with relationals::

    sage: e = x+1 <= x-2
    sage: e + 2
    x + 3 <= x
    sage: e - 1
    x <= x - 3
    sage: e*(-1)
    -x - 1 <= -x + 2
    sage: (-2)*e
    -2*x - 2 <= -2*x + 4
    sage: e*5
    5*x + 5 <= 5*x - 10
    sage: e/5
    1/5*x + 1/5 <= 1/5*x - 2/5
    sage: 5/e
    5/(x + 1) <= 5/(x - 2)
    sage: e/(-2)
    -1/2*x - 1/2 <= -1/2*x + 1
    sage: -2/e
    -2/(x + 1) <= -2/(x - 2)

We can even add together two relations, so long as the operators are
the same::

    sage: (x^3 + x <= x - 17)  + (-x <= x - 10)
    x^3 <= 2*x - 27

Here they aren't::

    sage: (x^3 + x <= x - 17)  + (-x >= x - 10)
    Traceback (most recent call last):
    ...
    TypeError: incompatible relations


ARBITRARY SAGE ELEMENTS:

You can work symbolically with any Sage data type.  This can lead to
nonsense if the data type is strange, e.g., an element of a finite
field (at present).

We mix Singular variables with symbolic variables::

    sage: R.<u,v> = QQ[]
    sage: var('a,b,c')
    (a, b, c)
    sage: expand((u + v + a + b + c)^2)
    a^2 + 2*a*b + 2*a*c + 2*a*u + 2*a*v + b^2 + 2*b*c + 2*b*u + 2*b*v + c^2 + 2*c*u + 2*c*v + u^2 + 2*u*v + v^2

TESTS:

Test Jacobian on Pynac expressions. #5546 ::

    sage: var('x,y')
    (x, y)
    sage: f = x + y
    sage: jacobian(f, [x,y])
    [1 1]


Test if matrices work #5546 ::

    sage: var('x,y,z')
    (x, y, z)
    sage: M = matrix(2,2,[x,y,z,x])
    sage: v = vector([x,y])
    sage: M * v
    (x^2 + y^2, x*y + x*z)
    sage: v*M
    (x^2 + y*z, 2*x*y)

Test if comparison bugs from #6256 are fixed::

    sage: t = exp(sqrt(x)); u = 1/t
    sage: t*u
    1
    sage: t + u
    e^sqrt(x) + e^(-sqrt(x))
    sage: t
    e^sqrt(x)
"""

include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"

import operator
import ring
import sage.rings.integer
import sage.rings.rational
from sage.structure.element cimport ModuleElement, RingElement, Element
from sage.symbolic.function import get_sfunction_from_serial
from sage.rings.rational import Rational  # Used for sqrt.
from sage.misc.derivative import multi_derivative
from sage.rings.infinity import AnInfinity

def is_Expression(x):
    """
    Returns True if *x* is a symbolic Expression.

    EXAMPLES::

        sage: from sage.symbolic.expression import is_Expression
        sage: is_Expression(x)
        True
        sage: is_Expression(2)
        False
        sage: is_Expression(SR(2))
        True
    """
    return isinstance(x, Expression)

def is_SymbolicEquation(x):
    """
    Returns True if *x* is a symbolic equation.

    EXAMPLES:

    The following two examples are symbolic equations::

        sage: from sage.symbolic.expression import is_SymbolicEquation
        sage: is_SymbolicEquation(sin(x) == x)
        True
        sage: is_SymbolicEquation(sin(x) < x)
        True
        sage: is_SymbolicEquation(x)
        False

    This is not, since ``2==3`` evaluates to the boolean
    ``False``::

        sage: is_SymbolicEquation(2 == 3)
        False

    However here since both 2 and 3 are coerced to be symbolic, we
    obtain a symbolic equation::

        sage: is_SymbolicEquation(SR(2) == SR(3))
        True

    """
    return isinstance(x, Expression) and is_a_relational((<Expression>x)._gobj)

cdef class Expression(CommutativeRingElement):
    cpdef object pyobject(self):
        """
        Get the underlying Python object corresponding to this
        expression, assuming this expression is a single numerical
        value.   Otherwise, a TypeError is raised.

        EXAMPLES::

            sage: var('x')
            x
            sage: b = -17/3
            sage: a = SR(b)
            sage: a.pyobject()
            -17/3
            sage: a.pyobject() is b
            True
        """
        cdef GConstant* c
        if is_a_constant(self._gobj):
            from sage.symbolic.constants import constants_name_table
            return constants_name_table[GEx_to_str(&self._gobj)]
        if not is_a_numeric(self._gobj):
            raise TypeError, "self must be a numeric expression"
        return py_object_from_numeric(self._gobj)

    def __init__(self, SR, x=0):
        """
        Nearly all expressions are created by calling new_Expression_from_*,
        but we need to make sure this at least doesn't leave self._gobj
        uninitialized and segfault.

        TESTS::

            sage: sage.symbolic.expression.Expression(SR)
            0
            sage: sage.symbolic.expression.Expression(SR, 5)
            5

        We test subclassing ``Expression``::

            sage: from sage.symbolic.expression import Expression
            sage: class exp_sub(Expression): pass
            sage: f = function('f')
            sage: t = f(x)
            sage: u = exp_sub(SR, t)
            sage: u.operator()
            f
        """
        self._parent = SR
        cdef Expression exp = self.coerce_in(x)
        GEx_construct_ex(&self._gobj, exp._gobj)

    def __dealloc__(self):
        """
        Delete memory occupied by this expression.
        """
        GEx_destruct(&self._gobj)

    def __getstate__(self):
        """
        Returns a tuple describing the state of this expression for pickling.

        This should return all information that will be required to unpickle
        the object. The functionality for unpickling is implemented in
        __setstate__().

        In order to pickle Expression objects, we return a tuple containing

         * 0  - as pickle version number
                in case we decide to change the pickle format in the feature
         * names of symbols of this expression
         * a string representation of self stored in a Pynac archive.

        TESTS::
            sage: var('x,y,z')
            (x, y, z)
            sage: t = 2*x*y^z+3
            sage: s = dumps(t)

            sage: t.__getstate__()
            (0,
             ['x', 'y', 'z'],
             ...)

        """
        cdef GArchive ar
        ar.archive_ex(self._gobj, "sage_ex")
        ar_str = GArchive_to_str(&ar)
        return (0, map(repr, self.variables()), ar_str)

    def __setstate__(self, state):
        """
        Initializes the state of the object from data saved in a pickle.

        During unpickling __init__ methods of classes are not called, the saved
        data is passed to the class via this function instead.

        TESTS::
            sage: var('x,y,z')
            (x, y, z)
            sage: t = 2*x*y^z+3
            sage: u = loads(dumps(t)) # indirect doctest
            sage: u
            2*y^z*x + 3
            sage: bool(t == u)
            True
            sage: u.subs(x=z)
            2*y^z*z + 3

            sage: loads(dumps(x.parent()(2)))
            2
        """
        # check input
        if state[0] != 0 or len(state) != 3:
            raise ValueError, "unknown state information"
        # set parent
        self._set_parent(ring.SR)
        # get variables
        cdef GExList sym_lst
        for name in state[1]:
            sym_lst.append_sym(get_symbol(name))

        # initialize archive
        cdef GArchive ar
        GArchive_from_str(&ar, state[2], len(state[2]))

        # extract the expression from the archive
        GEx_construct_ex(&self._gobj, ar.unarchive_ex(sym_lst, <unsigned>0))

    # TODO: The keyword argument simplify is for compatibility with
    # old symbolics, once the switch is complete, it should be removed
    def _repr_(self, simplify=None):
        """
        Return string representation of this symbolic expression.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: repr(x+y)
            'x + y'

        TESTS::

            # printing of modular number equal to -1 as coefficient
            sage: k.<a> = GF(9); k(2)*x
            2*x

            sage: (x+1)*Mod(6,7)
            6*x + 6

            #printing rational functions
            sage: x/y
            x/y
            sage: x/2/y
            1/2*x/y
            sage: .5*x/y
            0.500000000000000*x/y
            sage: x^(-1)
            1/x
            sage: x^(-5)
            x^(-5)
            sage: x^(-y)
            x^(-y)
            sage: 2*x^(-1)
            2/x
            sage: i*x
            I*x
            sage: -x.parent(i)
            -I
            sage: y + 3*(x^(-1))
            y + 3/x

        Printing the exp function::

            sage: x.parent(1).exp()
            e
            sage: x.exp()
            e^x

        Powers::

            sage: _ = var('A,B,n'); (A*B)^n
            (A*B)^n
            sage: (A/B)^n
            (A/B)^n
            sage: n*x^(n-1)
            n*x^(n - 1)
            sage: (A*B)^(n+1)
            (A*B)^(n + 1)
            sage: (A/B)^(n-1)
            (A/B)^(n - 1)
            sage: n*x^(n+1)
            n*x^(n + 1)
            sage: n*x^(n-1)
            n*x^(n - 1)
            sage: n*(A/B)^(n+1)
            n*(A/B)^(n + 1)
            sage: (n+A/B)^(n+1)
            (n + A/B)^(n + 1)

        Powers where the base or exponent is a Python object::

            sage: (2/3)^x
            (2/3)^x
            sage: x^CDF(1,2)
            x^(1.0 + 2.0*I)
            sage: (2/3)^(2/3)
            (2/3)^(2/3)
            sage: (-x)^(1/4)
            (-x)^(1/4)
            sage: k.<a> = GF(9)
            sage: SR(a+1)^x
            (a + 1)^x
        """
        return self._parent._repr_element_(self)

    def _interface_(self, I):
        """
        EXAMPLES::

            sage: f = sin(e + 2)
            sage: f._interface_(sage.calculus.calculus.maxima)
            sin(%e+2)
        """
        if is_a_constant(self._gobj):
            return self.pyobject()._interface_(I)
        return super(Expression, self)._interface_(I)

    def _maxima_(self, session=None):
        """
        EXAMPLES::

            sage: f = sin(e + 2)
            sage: f._maxima_()
            sin(%e+2)
            sage: _.parent() is sage.calculus.calculus.maxima
            True
        """
        if session is None:
            from sage.calculus.calculus import maxima
            return super(Expression, self)._interface_(maxima)
        else:
            return super(Expression, self)._interface_(session)

    def _interface_init_(self, I):
        """
        EXAMPLES::

            sage: a = (pi + 2).sin()
            sage: a._maxima_init_()
            'sin((%pi)+(2))'

            sage: a = (pi + 2).sin()
            sage: a._maple_init_()
            'sin((pi)+(2))'

            sage: a = (pi + 2).sin()
            sage: a._mathematica_init_()
            'Sin[(Pi)+(2)]'

            sage: f = pi + I*e
            sage: f._pari_init_()
            '(Pi)+((exp(1))*(I))'
        """
        from sage.symbolic.expression_conversions import InterfaceInit
        return InterfaceInit(I)(self)

    def _gap_init_(self):
        """
        Conversion of symbolic object to GAP always results in a GAP
        string.

        EXAMPLES::

            sage: gap(e + pi^2 + x^3)
            pi^2 + x^3 + e
        """
        return '"%s"'%repr(self)

    def _singular_init_(self):
        """
        Conversion of a symbolic object to Singular always results in a
        Singular string.

        EXAMPLES::

            sage: singular(e + pi^2 + x^3)
            pi^2 + x^3 + e
        """
        return '"%s"'%repr(self)

    def _magma_init_(self, magma):
        """
        Return string representation in Magma of this symbolic expression.

        Since Magma has no notation of symbolic calculus, this simply
        returns something that evaluates in Magma to a a Magma string.

        EXAMPLES::

            sage: x = var('x')
            sage: f = sin(cos(x^2) + log(x))
            sage: f._magma_init_(magma)
            '"sin(log(x) + cos(x^2))"'
            sage: magma(f)                         # optional - magma
            sin(cos(x^2) + log(x))
            sage: magma(f).Type()                  # optional - magma
            MonStgElt
        """
        return '"%s"'%repr(self)

    def _latex_(self):
        r"""
        Return string representation of this symbolic expression.

        EXAMPLES:


        TESTS::

            sage: var('x,y,z')
            (x, y, z)
            sage: latex(y + 3*(x^(-1)))
            y + 3 \, \frac{1}{x}
            sage: latex(x^(y+z^(1/y)))
            x^{z^{\frac{1}{y}} + y}
            sage: latex(1/sqrt(x+y))
            \frac{1}{\sqrt{x + y}}
            sage: latex(sin(x*(z+y)^x))
            \sin\left({(y + z)}^{x} x\right)
            sage: latex(3/2*(x+y)/z/y)
            \frac{3}{2} \, \frac{{(x + y)}}{y z}
            sage: latex((2^(x^y)))
            2^{x^{y}}
            sage: latex(abs(x))
            {\left| x \right|}
            sage: latex((x*y).conjugate())
            \overline{x} \overline{y}

        Check spacing of coefficients of mul expressions (#3202)::

            sage: latex(2*3^x)
            2 \, 3^{x}

        Powers::

            sage: _ = var('A,B,n')
            sage: latex((n+A/B)^(n+1))
            {(n + \frac{A}{B})}^{n + 1}
            sage: latex((A*B)^n)
            {(A B)}^{n}
            sage: latex((A*B)^(n-1))
            {(A B)}^{n - 1}

        Powers where the base or exponent is a Python object::

            sage: latex((2/3)^x)
            \left(\frac{2}{3}\right)^{x}
            sage: latex(x^CDF(1,2))
            x^{1.0 + 2.0i}
            sage: latex((2/3)^(2/3))
            \left(\frac{2}{3}\right)^{\frac{2}{3}}
            sage: latex((-x)^(1/4))
            {(-x)}^{\frac{1}{4}}
            sage: k.<a> = GF(9)
            sage: latex(SR(a+1)^x)
            \left(a + 1\right)^{x}
        """
        return self._parent._latex_element_(self)

    def _mathml_(self):
        """
        Returns a MathML representation of this object.

        EXAMPLES::

            sage: mathml(pi)
            <mi>&pi;</mi>
            sage: mathml(pi+2)
            MATHML version of the string pi + 2

        """
        from sage.misc.all import mathml
        try:
            obj = self.pyobject()
        except TypeError:
            return mathml(repr(self))
        return mathml(obj)

    def _integer_(self, ZZ=None):
        """
        EXAMPLES::

            sage: f = x^3 + 17*x -3
            sage: ZZ(f.coeff(x^3))
            1
            sage: ZZ(f.coeff(x))
            17
            sage: ZZ(f.coeff(x,0))
            -3
            sage: type(ZZ(f.coeff(x,0)))
            <type 'sage.rings.integer.Integer'>

        Coercion is done if necessary::

            sage: f = x^3 + 17/1*x
            sage: ZZ(f.coeff(x))
            17
            sage: type(ZZ(f.coeff(x)))
            <type 'sage.rings.integer.Integer'>

        If the symbolic expression is just a wrapper around an integer,
        that very same integer is returned::

            sage: n = 17; SR(n)._integer_() is n
            True
        """
        try:
            n = self.pyobject()
        except TypeError:
            raise TypeError, "unable to convert x (=%s) to an integer"%(self)
        if isinstance(n, sage.rings.integer.Integer):
            return n
        return sage.rings.integer.Integer(n)

    def __int__(self):
        """
        EXAMPLES::

            sage: int(sin(2)*100)
            90
            sage: int(log(8)/log(2))
            3
        """
        try:
            return int(self.pyobject())
        except (ValueError, TypeError):
            #FIXME:  This should be fixed so that it does something
            #smarter to handle the log(8)/log(2) case
            from sage.functions.all import floor, ceil
            f = float(self)
            if f > 0:
                return int(floor(self))
            else:
                return int(ceil(self))

    def __long__(self):
        """
        EXAMPLES::

            sage: long(sin(2)*100)
            90L
        """
        return long(int(self))

    def _rational_(self):
        """
        EXAMPLES::

            sage: f = x^3 + 17/1*x - 3/8
            sage: QQ(f.coeff(x^2))
            0
            sage: QQ(f.coeff(x^3))
            1
            sage: a = QQ(f.coeff(x)); a
            17
            sage: type(a)
            <type 'sage.rings.rational.Rational'>
            sage: QQ(f.coeff(x,0))
            -3/8

        If the symbolic expression is just a wrapper around a rational,
        that very same rational is returned::

            sage: n = 17/1; SR(n)._rational_() is n
            True
        """
        try:
            n = self.pyobject()
        except TypeError:
            raise TypeError, "unable to convert %s to a rational"%self
        if isinstance(n, sage.rings.rational.Rational):
            return n
        return sage.rings.rational.Rational(n)

    cpdef _eval_self(self, R):
        """
        Evaluate this expression numerically, and try to coerce it to R.

        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: sin(x).subs(x=5)._eval_self(RR)
            -0.958924274663138
            sage: gamma(x).subs(x=I)._eval_self(CC)
            -0.154949828301811 - 0.498015668118356*I
            sage: x._eval_self(CC)
            Traceback (most recent call last):
            ...
            TypeError: Cannot evaluate symbolic expression to a numeric value.
        """
        cdef GEx res = self._gobj.evalf(0, R.prec())
        if is_a_numeric(res):
            return R(py_object_from_numeric(res))
        else:
            raise TypeError, "Cannot evaluate symbolic expression to a numeric value."

    def _convert(self, typ):
        """
        Convert self to the given type by converting each of the operands
        to that type and doing the arithmetic.

        FIXME: Make sure these docs are correct with the new symbolics.

        EXAMPLES::

            sage: f = sqrt(2) * cos(3); f
            sqrt(2)*cos(3)
            sage: f._convert(RDF)
            -1.40006081534
            sage: f._convert(float)
            -1.4000608153399503

        Converting to an int can have surprising consequences, since Python
        int is "floor" and one individual factor can floor to 0 but the
        product doesn't::

            sage: int(f)
            -1
            sage: f._convert(int)
            0
            sage: int(sqrt(2))
            1
            sage: int(cos(3))
            0

        TESTS: This illustrates how the conversion works even when a type
        exception is raised, since here one operand is still x (in the
        unsimplified form)::

            sage: f = sin(SR(0))*x
            sage: f._convert(CDF)
            0
        """
        operator = self.operator()
        if operator is None:
            return typ(self.pyobject())
        else:
            args = [typ(op) for op in self.operands()]
            if len(args) == 1:
                return operator(*args)
            else:
                return reduce(operator, args)

    def _mpfr_(self, R):
        """
        Return a numerical approximation to this expression in the RealField R.

        The precision of the approximation is determined by the precision of
        the input R.

        EXAMPLES::

            0.090909090909090909090909090909090909090909090909090909090909

            sage: a = sin(3); a
            sin(3)
            sage: RealField(200)(a)
            0.14112000805986722210074480280811027984693326425226558415188
            sage: a._mpfr_(RealField(100))
            0.14112000805986722210074480281
        """
        return self._eval_self(R)

    def _real_mpfi_(self, R):
        """
        Returns this expression as a real interval.

        EXAMPLES::

            sage: RIF(sqrt(2))
            1.414213562373095?
        """
        from sage.symbolic.expression_conversions import RingConverter
        return RingConverter(R)(self)

    def _complex_mpfi_(self, R):
        """
        Returns this expression as a complex interval.

        EXAMPLES::

            sage: CIF(pi)
            3.141592653589794?
        """
        from sage.symbolic.expression_conversions import RingConverter
        return RingConverter(R)(self)

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: RDF(sin(3))
            0.14112000806
        """
        return self._eval_self(R)

    def _complex_mpfr_field_(self, R):
        """
        Return a numerical approximation to this expression in the given
        ComplexField R.

        The precision of the approximation is determined by the precision of
        the input R.

        EXAMPLES::

            sage: ComplexField(200)(SR(1/11))
            0.090909090909090909090909090909090909090909090909090909090909
            sage: zeta(x).subs(x=I)._complex_mpfr_field_(ComplexField(70))
            0.0033002236853241028742 - 0.41815544914132167669*I
            sage: gamma(x).subs(x=I)._complex_mpfr_field_(ComplexField(60))
            -0.15494982830181069 - 0.49801566811835604*I
            sage: log(x).subs(x=I)._complex_mpfr_field_(ComplexField(50))
            1.5707963267949*I

            sage: CC(sqrt(2))
            1.41421356237310
            sage: a = sqrt(-2); a
            sqrt(-2)
            sage: CC(a).imag()
            1.41421356237309
            sage: ComplexField(200)(a).imag()
            1.4142135623730950488016887242096980785696718753769480731767
            sage: ComplexField(100)((-1)^(1/10))
            0.95105651629515357211643933338 + 0.30901699437494742410229341718*I
            sage: CC(x*sin(0))
            0
        """
        return self._eval_self(R)

    def _complex_double_(self, R):
        """
        Return a numerical approximation to this expression in the given
        Complex Double Field R.

        EXAMPLES::

            sage: CDF(SR(1/11))
            0.0909090909091
            sage: zeta(x).subs(x=I)._complex_double_(CDF)
            0.00330022368532 - 0.418155449141*I
            sage: gamma(x).subs(x=I)._complex_double_(CDF)
            -0.154949828302 - 0.498015668118*I
            sage: log(x).subs(x=I)._complex_double_(CDF)
            1.57079632679*I
            sage: CDF((-1)^(1/3))
            0.5 + 0.866025403784*I
        """
        return self._eval_self(R)

    def __float__(self):
        """
        Return float conversion of self, assuming self is constant.
        Otherwise, raise a TypeError.

        OUTPUT:

        - float - double precision evaluation of self

        EXAMPLES::

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
            TypeError: a float is required
        """
        cdef bint success
        cdef double ans = GEx_to_double(self._gobj, &success)
        if not success:
            raise TypeError, "float() argument must be a string or a number"
        return ans

    def __complex__(self):
        """
        EXAMPLES::

            sage: complex(I)
            1j
            sage: complex(erf(3*I))
            Traceback (most recent call last):
            ...
            TypeError: unable to simplify to complex approximation
        """
        try:
            return complex(self.n(53))
        except TypeError:
            raise TypeError, "unable to simplify to complex approximation"

    def _sympy_(self):
        """
        Returns a Sympy version of this object.

        EXAMPLES::

            sage: pi._sympy_()
            pi
            sage: type(_)
            <class 'sympy.core.numbers.Pi'>

        """
        from sage.symbolic.expression_conversions import sympy
        return sympy(self)

    def _algebraic_(self, field):
        """
        Convert a symbolic expression to an algebraic number.

        EXAMPLES::

            sage: QQbar(sqrt(2) + sqrt(8))
            4.242640687119285?
            sage: AA(sqrt(2) ^ 4) == 4
            True
            sage: AA(-golden_ratio)
            -1.618033988749895?
            sage: QQbar((2*I)^(1/2))
            1 + 1*I
            sage: QQbar(e^(pi*I/3))
            0.500000000000000? + 0.866025403784439?*I

            sage: QQbar(sqrt(2))
            1.414213562373095?
            sage: AA(abs(1+I))
            1.414213562373095?
            sage: golden_ratio._algebraic_(QQbar)
            1.618033988749895?
            sage: QQbar(golden_ratio)
            1.618033988749895?

            sage: AA(x*sin(0))
            0
            sage: QQbar(x*sin(0))
            0
        """
        from sage.symbolic.expression_conversions import algebraic
        return algebraic(self, field)

    def __hash__(self):
        """
        Return hash of this expression.

        EXAMPLES::

        The hash of an object in Python or its coerced version into
        the symbolic ring is the same::

            sage: hash(SR(3/1))
            3
            sage: hash(SR(19/23))
            4
            sage: hash(19/23)
            4

        The hash for symbolic expressions are unfortunately random. Here we
        only test that the hash() function returns without error, and that
        the return type is correct::

            sage: x, y = var("x y")
            sage: t = hash(x); type(t)
            <type 'int'>
            sage: t = hash(x^y); type(t)
            <type 'int'>
            sage: type(hash(x+y))
            <type 'int'>
            sage: d = {x+y: 5}
            sage: d
            {x + y: 5}

        In this example hashing is important otherwise the answer is
        wrong::

            sage: uniq([x-x, -x+x])
            [0]

        Test if exceptions during hashing are handled properly::

            sage: t = SR(matrix(2,2,range(4)))
            sage: hash(t)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        TESTS:

        Test if hashes for fderivatives with different parameters collide.
        #6243::

            sage: f = function('f'); t = f(x,y)
            sage: u = t.derivative(x); v = t.derivative(y)
            sage: hash(u) == hash(v)
            False
            sage: d = {u: 3, v: 5}; sorted(d.values())
            [3, 5]
        """
        return self._gobj.gethash()

    # Boilerplate code from sage/structure/element.pyx
    def __richcmp__(left, right, int op):
        """
        Create a formal symbolic inequality or equality.

        EXAMPLES::

            sage: var('x, y')
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
        return (<Element>left)._richcmp(right, op)

    cdef _richcmp_c_impl(left, Element right, int op):
        cdef Expression l, r

        l = left
        r = right

        if is_a_relational(l._gobj) or is_a_relational(r._gobj):
            c = cmp(hash(l), hash(r))
            if op == Py_NE:
                return c != 0
            elif op == Py_EQ:
                return c == 0
            else:
                return False

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
        else:
            raise TypeError
        return new_Expression_from_GEx(l._parent, e)

    def assume(self):
        r"""
        Assume that this equation holds. This is relevant for symbolic
        integration, among other things.

        EXAMPLES: We call the assume method to assume that `x>2`::

            sage: (x > 2).assume()

        Bool returns True below if the inequality is *definitely* known to
        be True.

        ::

            sage: bool(x > 0)
            True
            sage: bool(x < 0)
            False

        This may or may not be True, so bool returns False::

            sage: bool(x > 3)
            False

        TESTS::

            sage: v,c = var('v,c')
            sage: assume(c != 0)
            sage: integral((1+v^2/c^2)^3/(1-v^2/c^2)^(3/2),v)
            -75/8*sqrt(c^2)*arcsin(sqrt(c^2)*v/c^2) - 17/8*v^3/(sqrt(-v^2/c^2 + 1)*c^2) - 1/4*v^5/(sqrt(-v^2/c^2 + 1)*c^4) + 83/8*v/sqrt(-v^2/c^2 + 1)
        """
        from sage.symbolic.assumptions import _assumptions
        from sage.calculus.calculus import maxima
        if not self.is_relational():
            raise TypeError, "self (=%s) must be a relational expression"%self
        if not self in _assumptions:
            m = self._maxima_init_assume_()
            maxima.assume(m)
            _assumptions.append(self)

    def forget(self):
        """
        Forget the given constraint.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: forget()
            sage: assume(x>0, y < 2)
            sage: assumptions()
            [x > 0, y < 2]
            sage: forget(y < 2)
            sage: assumptions()
            [x > 0]
        """
        from sage.symbolic.assumptions import _assumptions
        from sage.calculus.calculus import maxima
        if not self.is_relational():
            raise TypeError, "self (=%s) must be a relational expression"%self
        m = self._maxima_init_assume_()
        maxima.forget(m)
        try:
            _assumptions.remove(self)
        except ValueError:
            pass

    def _maxima_init_assume_(self):
        """
        Return string that when evaluated in Maxima defines the assumption
        determined by this expression.

        EXAMPLES::

            sage: f = x+2 > sqrt(3)
            sage: f._maxima_init_assume_()
            '((x)+(2))>((3/1)^(1/2))'
        """
        from sage.calculus.calculus import maxima
        l = self.lhs()._maxima_init_()
        r = self.rhs()._maxima_init_()
        op = self.operator()
        if  op is operator.eq:
            m = 'equal(%s, %s)'%(l, r)
        elif op is operator.ne:
            m = 'notequal(%s, %s)'%(l, r)
        else:
            m = '(%s)%s(%s)' % (l, maxima._relation_symbols()[op], r)
        return m

    cpdef bint is_polynomial(self, var):
        """
        Return True if self is a polynomial in the given variable.

        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: t = x^2 + y; t
            x^2 + y
            sage: t.is_polynomial(x)
            True
            sage: t.is_polynomial(y)
            True
            sage: t.is_polynomial(z)
            True

            sage: t = sin(x) + y; t
            y + sin(x)
            sage: t.is_polynomial(x)
            False
            sage: t.is_polynomial(y)
            True
            sage: t.is_polynomial(sin(x))
            True
        """
        cdef Expression symbol0 = self.coerce_in(var)
        return self._gobj.is_polynomial(symbol0._gobj)

    cpdef bint is_relational(self):
        """
        Return True if self is a relational expression.

        EXAMPLES::

            sage: x = var('x')
            sage: eqn = (x-1)^2 == x^2 - 2*x + 3
            sage: eqn.is_relational()
            True
            sage: sin(x).is_relational()
            False
        """
        return is_a_relational(self._gobj)

    def left_hand_side(self):
        """
        If self is a relational expression, return the left hand side
        of the relation.  Otherwise, raise a ValueError.

        EXAMPLES::

            sage: x = var('x')
            sage: eqn = (x-1)^2 == x^2 - 2*x + 3
            sage: eqn.left_hand_side()
            (x - 1)^2
            sage: eqn.lhs()
            (x - 1)^2
            sage: eqn.left()
            (x - 1)^2
        """
        if not self.is_relational():
            raise ValueError, "self must be a relational expression"
        return new_Expression_from_GEx(self._parent, self._gobj.lhs())

    lhs = left = left_hand_side

    def right_hand_side(self):
        """
        If self is a relational expression, return the right hand side
        of the relation.  Otherwise, raise a ValueError.

        EXAMPLES::

            sage: x = var('x')
            sage: eqn = (x-1)^2 <= x^2 - 2*x + 3
            sage: eqn.right_hand_side()
            x^2 - 2*x + 3
            sage: eqn.rhs()
            x^2 - 2*x + 3
            sage: eqn.right()
            x^2 - 2*x + 3
        """
        if not self.is_relational():
            raise ValueError, "self must be a relation"
        return new_Expression_from_GEx(self._parent, self._gobj.rhs())

    rhs = right = right_hand_side

    def __nonzero__(self):
        """
        Return True if this element is definitely not zero.

        EXAMPLES::

            sage: x = var('x')
            sage: forget()
            sage: SR(0).__nonzero__()
            False
            sage: SR(1).__nonzero__()
            True
            sage: bool(abs(x))
            True
            sage: bool(x/x - 1)
            False

            sage: k = var('k')
            sage: pol = 1/(k-1) - 1/k -1/k/(k-1);
            sage: pol.is_zero()
            True

            sage: f = sin(x)^2 + cos(x)^2 - 1
            sage: f.is_zero()
            True

        A bunch of tests of nonzero (which is called by bool) for
        symbolic relations::

            sage: x = var('x')
            sage: bool((x-1)^2 == x^2 - 2*x + 1)
            True
            sage: bool(((x-1)^2 == x^2 - 2*x + 1).expand())
            True
            sage: bool(((x-1)^2 == x^2 - 2*x + 3).expand())
            False
            sage: bool(2 + x < 3 + x)
            True
            sage: bool(2 + x < 1 + x)
            False
            sage: bool(2 + x > 1 + x)
            True
            sage: bool(1 + x > 1 + x)
            False
            sage: bool(1 + x >= 1 + x)
            True
            sage: bool(1 + x < 1 + x)
            False
            sage: bool(1 + x <= 1 + x)
            True
            sage: bool(1 + x^2 != 1 + x*x)
            False
            sage: bool(1 + x^2 != 2 + x*x)
            True
            sage: bool(SR(oo) == SR(oo))
            True
            sage: bool(-SR(oo) == SR(oo))
            False
            sage: bool(-SR(oo) != SR(oo))
            True
        """
        if self.is_relational():
            #If both the left hand side and right hand side are wrappers
            #around Sage objects, then we should do the comparison directly
            #with those objects
            if is_a_constant(self._gobj.lhs()) and is_a_constant(self._gobj.rhs()):
                return self.operator()(self.lhs().pyobject(), self.rhs().pyobject())

            res = relational_to_bool(self._gobj)
            if res:
                return True

            # If assumptions are involved, falsification is more complicated...
            need_assumptions = False
            vars = self.variables()
            if len(vars) > 0:
                from sage.symbolic.assumptions import assumptions
                assumption_list = assumptions()
                if len(assumption_list) > 0:
                    assumption_vars = set(sum([eqn.variables() for eqn in assumptions()], ()))
                    if set(vars).intersection(assumption_vars):
                        need_assumptions = True

            # Use interval fields to try and falsify the relation
            if not need_assumptions:
                res = self.test_relation()
                if res is True:
                    return True
                elif res is False:
                    return False

            # we really have to do some work here...
            # I really don't like calling Maxima to test equality.  It
            # is SUPER SUPER SLOW, and it has all the problem
            # associated with different semantics, different
            # precision, etc., that can lead to subtle bugs.  Also, a
            # lot of basic Sage objects can't be put into maxima.
            from sage.symbolic.relation import test_relation_maxima
            return test_relation_maxima(self)

        self_is_zero = self._gobj.is_zero()
        if self_is_zero:
            return False
        else:
            return not bool(self == self._parent.zero_element())

    def test_relation(self, int ntests=20, domain=None, proof=True):
        """
        Test this relation at several random values, attempting to find
        a contradiction. If this relation has no variables, it will also
        test this relation after casting into the domain.

        Because the interval fields never return false positives, we can be
        assured that if True or False is returned (and proof is False) then
        the answer is correct.

        INPUT::

           ntests -- (default 20) the number of iterations to run
           domain -- (optional) the domain from which to draw the random values
                     defaults to CIF for equality testing and RIF for
                     order testing
           proof --  (default True) if False and the domain is an interval field,
                     regard overlapping (potentially equal) intervals as equal,
                     and return True if all tests succeeded.

        OUTPUT::

            True - this relation holds in the domain and has no variables
            False - a contradiction was found
            NotImplemented - no contradiction found

        EXAMPLES::

            sage: (3 < pi).test_relation()
            True
            sage: (0 >= pi).test_relation()
            False
            sage: (exp(pi) - pi).n()
            19.9990999791895
            sage: (exp(pi) - pi == 20).test_relation()
            False
            sage: (sin(x)^2 + cos(x)^2 == 1).test_relation()
            NotImplemented
            sage: (sin(x)^2 + cos(x)^2 == 1).test_relation(proof=False)
            True
            sage: (x == 1).test_relation()
            False
            sage: var('x,y')
            (x, y)
            sage: (x < y).test_relation()
            False

        TESTS::

            sage: all_relations = [op for name, op in sorted(operator.__dict__.items()) if len(name) == 2]
            sage: all_relations
            [<built-in function eq>, <built-in function ge>, <built-in function gt>, <built-in function le>, <built-in function lt>, <built-in function ne>]
            sage: [op(3, pi).test_relation() for op in all_relations]
            [False, False, False, True, True, True]
            sage: [op(pi, pi).test_relation() for op in all_relations]
            [True, True, False, True, False, False]

            sage: s = 'some_very_long_variable_name_which_will_definitely_collide_if_we_use_a_reasonable_length_bound_for_a_hash_that_respects_lexicographic_order'
            sage: t1, t2 = var(','.join([s+'1',s+'2']))
            sage: (t1 == t2).test_relation()
            False
        """
        cdef int k, eq_count = 0
        cdef bint is_interval
        if not self.is_relational():
            raise ValueError, "self must be a relation"
        cdef operators op = relational_operator(self._gobj)
        from sage.rings.real_mpfi import is_RealIntervalField
        from sage.rings.complex_interval_field import is_ComplexIntervalField
        from sage.rings.all import RIF, CIF
        if domain is None:
            is_interval = True
            if op == equal or op == not_equal:
                domain = CIF
            else:
                domain = RIF
        else:
            is_interval = is_RealIntervalField(domain) or is_ComplexIntervalField(domain)
        zero = domain(0)
        diff = self.lhs() - self.rhs()
        vars = diff.variables()
        if op == equal:
            falsify = operator.ne
        elif op == not_equal:
            falsify = operator.eq
        elif op == less:
            falsify = operator.ge
        elif op == less_or_equal:
            falsify = operator.gt
        elif op == greater:
            falsify = operator.le
        elif op == greater_or_equal:
            falsify = operator.lt
        cdef bint equality_ok = op in [equal, less_or_equal, greater_or_equal]
        cdef int errors = 0
        if len(vars) == 0:
            val = None
            try:
                val = domain(diff)
            except (TypeError, ValueError, ArithmeticError), ex:
                pass
            else:
                if self.operator()(val, zero):
                    return True
                elif falsify(val, zero):
                    return False
                if is_interval and not proof:
                    if val.contains_zero():
                        return equality_ok
                    else:
                        return not equality_ok
        else:
            for k in range(ntests):
                try:
                    if is_interval:
                        # Let's up the prec
                        if val and k > 4 and val.contains_zero() and domain.prec() < 1000:
                            domain = domain.to_prec(int(domain.prec() * 1.5))
                        # Uniform [-1,1] isn't the best distribution to use...
                        var_dict = dict([(v, domain.random_element() * domain.random_element(-2,6).exp()) for v in vars])
                    else:
                        var_dict = dict([(v, domain.random_element()) for v in vars])
                    val = domain(diff.subs(var_dict))
                    if falsify(val, zero):
                        return False
                    if is_interval:
                        eq_count += <bint>val.contains_zero()
                except (TypeError, ValueError, ArithmeticError), ex:
                    errors += 1
                    if k == errors > 3 and is_ComplexIntervalField(domain):
                        domain = RIF.to_prec(domain.prec())
                    # we are plugging in random values above, don't be surprised
                    # if something goes wrong...
                    eq_count += equality_ok

        if not proof:
            if not equality_ok:
                return eq_count == 0
            elif op == equal and is_interval:
                return eq_count == ntests
            else:
                return True
        # Nothing failed, so it *may* be True, but this method doesn't wasn't
        # able to find anything.
        return NotImplemented

    def is_unit(self):
        """
        Return True if this expression is a unit of the symbolic ring.

        EXAMPLES::

            sage: SR(1).is_unit()
            True
            sage: SR(-1).is_unit()
            True
            sage: SR(0).is_unit()
            False
        """
        if not not self:
            return True
        if self == 0:
            return False
        return NotImplementedError

    cdef Expression coerce_in(self, z):
        """
        Quickly coerce z to be an Expression.
        """
        cdef Expression w
        try:
            w = z
            return w
        except TypeError:
            return self._parent._coerce_(z)

    cpdef ModuleElement _add_(left, ModuleElement right):
        """
        Add left and right.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: x + y + y + x
            2*x + 2*y

            # adding relational expressions
            sage: ( (x+y) > x ) + ( x > y )
            2*x + y > x + y

            sage: ( (x+y) > x ) + x
            2*x + y > 2*x

        TESTS::

            sage: x + ( (x+y) > x )
            2*x + y > 2*x

            sage: ( x > y) + (y < x)
            Traceback (most recent call last):
            ...
            TypeError: incompatible relations

            sage: (x < 1) + (y <= 2)
            x + y < 3

            sage: x + oo
            +Infinity
            sage: x - oo
            -Infinity
            sage: x + unsigned_infinity
            Infinity
            sage: x - unsigned_infinity
            Infinity

            sage: nsr = x.parent()
            sage: nsr(oo) + nsr(oo)
            +Infinity
            sage: nsr(-oo) + nsr(-oo)
            -Infinity
            sage: nsr(oo) - nsr(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: Infinity - Infinity encountered.
            sage: nsr(-oo) - nsr(-oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: Infinity - Infinity encountered.

            sage: nsr(unsigned_infinity) + nsr(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: unsigned_infinity + x where x is Infinity, -Infinity or unsigned infinity encountered.
            sage: nsr(unsigned_infinity) - nsr(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: unsigned_infinity + x where x is Infinity, -Infinity or unsigned infinity encountered.
            sage: nsr(oo) + nsr(unsigned_infinity)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: unsigned_infinity + x where x is Infinity, -Infinity or unsigned infinity encountered.
            sage: nsr(oo) - nsr(unsigned_infinity)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: unsigned_infinity + x where x is Infinity, -Infinity or unsigned infinity encountered.
            sage: nsr(unsigned_infinity) + nsr(unsigned_infinity)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: unsigned_infinity + x where x is Infinity, -Infinity or unsigned infinity encountered.


        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        cdef operators op
        if is_a_relational(left._gobj):
            if is_a_relational(_right._gobj):
                op = compatible_relation(relational_operator(left._gobj),
                                         relational_operator(_right._gobj))
                x = relational(gadd(left._gobj.lhs(), _right._gobj.lhs()),
                               gadd(left._gobj.rhs(), _right._gobj.rhs()),
                               op)
            else:
                x = relational(gadd(left._gobj.lhs(), _right._gobj),
                               gadd(left._gobj.rhs(), _right._gobj),
                               relational_operator(left._gobj))
        elif is_a_relational(_right._gobj):
            x = relational(gadd(left._gobj, _right._gobj.lhs()),
                           gadd(left._gobj, _right._gobj.rhs()),
                           relational_operator(_right._gobj))
        else:
            x = gadd(left._gobj, _right._gobj)
        return new_Expression_from_GEx(left._parent, x)

    cpdef ModuleElement _sub_(left, ModuleElement right):
        """
        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: x - y
            x - y

            # subtracting relational expressions
            sage: ( (x+y) > x ) - ( x > y )
            y > x - y

            sage: ( (x+y) > x ) - x
            y > 0

        TESTS::

            sage: x - ( (x+y) > x )
            -y > 0

            sage: ( x > y) - (y < x)
            Traceback (most recent call last):
            ...
            TypeError: incompatible relations

            sage: x - oo
            -Infinity
            sage: oo - x
            +Infinity
        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        if is_a_relational(left._gobj):
            if is_a_relational(_right._gobj):
                op = compatible_relation(relational_operator(left._gobj),
                                         relational_operator(_right._gobj))
                x = relational(gsub(left._gobj.lhs(), _right._gobj.lhs()),
                               gsub(left._gobj.rhs(), _right._gobj.rhs()),
                               op)
            else:
                x = relational(gsub(left._gobj.lhs(), _right._gobj),
                               gsub(left._gobj.rhs(), _right._gobj),
                               relational_operator(left._gobj))
        elif is_a_relational(_right._gobj):
            x = relational(gsub(left._gobj, _right._gobj.lhs()),
                           gsub(left._gobj, _right._gobj.rhs()),
                           relational_operator(_right._gobj))
        else:
            x = gsub(left._gobj, _right._gobj)
        return new_Expression_from_GEx(left._parent, x)

    cpdef RingElement _mul_(left, RingElement right):
        """
        Multiply left and right.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: x*y*y
            x*y^2

            # multiplying relational expressions
            sage: ( (x+y) > x ) * ( x > y )
            (x + y)*x > x*y

            sage: ( (x+y) > x ) * x
            (x + y)*x > x^2

            sage: ( (x+y) > x ) * -1
            -x - y > -x

        TESTS::

            sage: x * ( (x+y) > x )
            (x + y)*x > x^2

            sage: ( x > y) * (y < x)
            Traceback (most recent call last):
            ...
            TypeError: incompatible relations

            sage: a = 1000 + 300*x + x^3 + 30*x^2
            sage: a*Mod(1,7)
            x^3 + 2*x^2 + 6*x + 6

            sage: var('z')
            z
            sage: 3*(x+y)/z
            3*(x + y)/z
            sage: (-x+z)*(3*x-3*z)
            -3*(x - z)^2

            # check if comparison of constant terms in Pynac add objects work
            sage: (y-1)*(y-2)
            (y - 2)*(y - 1)

        Check for simplifications when multiplying instances of exp::

            sage: exp(x)*exp(y)
            e^(x + y)
            sage: exp(x)^2*exp(y)
            e^(2*x + y)
            sage: x^y*exp(x+y)*exp(-y)
            x^y*e^x
            sage: x^y*exp(x+y)*(x+y)*(2*x+2*y)*exp(-y)
            2*(x + y)^2*x^y*e^x
            sage: x^y*exp(x+y)*(x+y)*(2*x+2*y)*exp(-y)*exp(z)^2
            2*(x + y)^2*x^y*e^(x + 2*z)
            sage: 1/exp(x)
            e^(-x)
            sage: exp(x)/exp(y)
            e^(x - y)
            sage: A = exp(I*pi/5)
            sage: t = A*A*A*A; t
            e^(4/5*I*pi)
            sage: t*A
            -1
            sage: b = -x*A; c = b*b; c
            x^2*e^(2/5*I*pi)
            sage: u = -t*A; u
            1

            sage: x*oo
            +Infinity
            sage: -x*oo
            -Infinity
            sage: x*unsigned_infinity
            Traceback (most recent call last):
            ...
            ValueError: oo times number < oo not defined

            sage: SR(oo)*SR(oo)
            +Infinity
            sage: SR(-oo)*SR(oo)
            -Infinity
            sage: SR(oo)*SR(-oo)
            -Infinity
            sage: SR(unsigned_infinity)*SR(oo)
            Infinity

        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        cdef operators o
        if is_a_relational(left._gobj):
            if is_a_relational(_right._gobj):
                op = compatible_relation(relational_operator(left._gobj),
                                         relational_operator(_right._gobj))
                x = relational(gmul(left._gobj.lhs(), _right._gobj.lhs()),
                               gmul(left._gobj.rhs(), _right._gobj.rhs()),
                               op)
            else:
                o = relational_operator(left._gobj)
                x = relational(gmul(left._gobj.lhs(), _right._gobj),
                               gmul(left._gobj.rhs(), _right._gobj),
                               o)
        elif is_a_relational(_right._gobj):
            o = relational_operator(_right._gobj)
            x = relational(gmul(left._gobj, _right._gobj.lhs()),
                           gmul(left._gobj, _right._gobj.rhs()),
                           o)
        else:
            x = gmul(left._gobj, _right._gobj)
        return new_Expression_from_GEx(left._parent, x)

    cpdef RingElement _div_(left, RingElement right):
        """
        Divide left and right.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: x/y/y
            x/y^2

            # dividing relational expressions
            sage: ( (x+y) > x ) / ( x > y )
            (x + y)/x > x/y

            sage: ( (x+y) > x ) / x
            (x + y)/x > 1

            sage: ( (x+y) > x ) / -1
            -x - y > -x

        TESTS::

            sage: x / ( (x+y) > x )
            x/(x + y) > 1

            sage: ( x > y) / (y < x)
            Traceback (most recent call last):
            ...
            TypeError: incompatible relations
            sage: x/oo
            0
            sage: oo/x
            +Infinity

            sage: SR(oo)/SR(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0*infinity encountered.

            sage: SR(-oo)/SR(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0*infinity encountered.

            sage: SR(oo)/SR(-oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0*infinity encountered.

            sage: SR(oo)/SR(unsigned_infinity)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0*infinity encountered.

            sage: SR(unsigned_infinity)/SR(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0*infinity encountered.

            sage: SR(0)/SR(oo)
            0

            sage: SR(0)/SR(unsigned_infinity)
            0

            sage: x/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Symbolic division by zero
        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        cdef operators o
        try:
            if is_a_relational(left._gobj):
                if is_a_relational(_right._gobj):
                    op = compatible_relation(relational_operator(left._gobj),
                                             relational_operator(_right._gobj))
                    x = relational(gdiv(left._gobj.lhs(), _right._gobj.lhs()),
                                   gdiv(left._gobj.rhs(), _right._gobj.rhs()),
                                   op)
                else:
                    o = relational_operator(left._gobj)
                    x = relational(gdiv(left._gobj.lhs(), _right._gobj),
                                   gdiv(left._gobj.rhs(), _right._gobj),
                                   o)
            elif is_a_relational(_right._gobj):
                o = relational_operator(_right._gobj)
                x = relational(gdiv(left._gobj, _right._gobj.lhs()),
                               gdiv(left._gobj, _right._gobj.rhs()),
                               o)
            else:
                x = gdiv(left._gobj, _right._gobj)
            return new_Expression_from_GEx(left._parent, x)
        except Exception, msg:
            # TODO: change this to maybe cleverly do something involving Cython C++ exception handling.
            # See http://docs.cython.org/docs/wrapping_CPlusPlus.html
            if 'division by zero' in str(msg):
                raise ZeroDivisionError, "Symbolic division by zero"
            else:
                raise

    def __invert__(self):
        """
        Return the inverse of this symbolic expression.

        EXAMPLES::

            sage: ~x
            1/x
            sage: ~SR(3)
            1/3
            sage: v1=var('v1'); a = (2*erf(2*v1*arcsech(0))/v1); ~a
            1/2*v1/erf(2*v1*arcsech(0))
        """
        return 1/self

    # Boilerplate code from sage/structure/element.pyx
    def __cmp__(left, right):
        """
        Compare self and right, returning -1, 0, or 1, depending on if
        self < right, self == right, or self > right, respectively.

        Use this instead of the operators <=, <, etc. to compare symbolic
        expressions when you do not want to get a formal inequality back.

        IMPORTANT: Both self and right *must* have the same type, or
        this function won't be called.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: x.__cmp__(y)
            -1
            sage: x < y
            x < y
            sage: cmp(x,y)
            -1
            sage: cmp(SR(0.5), SR(0.7))
            -1
            sage: SR(0.5) < SR(0.7)
            0.500000000000000 < 0.700000000000000
            sage: cmp(SR(0.5), 0.7)
            -1
            sage: cmp(sin(SR(2)), sin(SR(1)))
            1
            sage: float(sin(SR(2)))
            0.90929742682568171
            sage: float(sin(SR(1)))
            0.8414709848078965
        """
        return (<Element>left)._cmp(right)

    cdef int _cmp_c_impl(left, Element right) except -2:
        return left._gobj.compare((<Expression>right)._gobj)

    def __pow__(self, exp, ignored):
        """
        Return self raised to the power of exp.

        INPUT:
            self -- symbolic expression
            exp -- something that coerces to a symbolic expressions

        OUTPUT:
            symbolic expression

        EXAMPLES:
            sage: var('x,y')
            (x, y)
            sage: x.__pow__(y)
            x^y
            sage: x^(3/5)
            x^(3/5)
            sage: x^sin(x)^cos(y)
            x^(sin(x)^cos(y))

        TESTS::

            sage: (Mod(2,7)*x^2 + Mod(2,7))^7
            (2*x^2 + 2)^7
            sage: k = GF(7)
            sage: f = expand((k(1)*x^5 + k(1)*x^2 + k(2))^7); f
            x^35 + x^14 + 2

            sage: x^oo
            Traceback (most recent call last):
            ...
            RuntimeError: power::eval(): pow(x, Infinity) for non numeric x is not defined.
            sage: SR(oo)^2
            +Infinity
            sage: SR(-oo)^2
            +Infinity
            sage: SR(-oo)^3
            -Infinity
            sage: SR(unsigned_infinity)^2
            Infinity

        Test powers of exp::

            sage: exp(2)^5
            e^10
            sage: exp(x)^5
            e^(5*x)

        Test base a Python numeric type::

            sage: int(2)^x
            2^x
            sage: float(2.3)^(x^3 - x^2 + 1/3)
            2.30000000000000^(x^3 - x^2 + 1/3)
            sage: complex(1,3)^(sqrt(2))
            (1.00000000000000 + 3.00000000000000*I)^sqrt(2)
        """
        cdef Expression base, nexp

        try:
            # self is an Expression and exp might not be
            base = self
            nexp = base.coerce_in(exp)
        except TypeError:
            # exp is an Expression and self might not be
            nexp = exp
            base = nexp.coerce_in(self)
        cdef GEx x
        if is_a_relational(base._gobj):
            x = relational(g_pow(base._gobj.lhs(), nexp._gobj),
                           g_pow(base._gobj.rhs(), nexp._gobj),
                           relational_operator(base._gobj))
        else:
            x = g_pow(base._gobj, nexp._gobj)
        return new_Expression_from_GEx(base._parent, x)

    def derivative(self, *args):
        """
        Returns the derivative of this expressions with respect to the
        variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more
        details.

        .. seealso::

            :meth:`_derivative`

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: t = (x^2+y)^2
            sage: t.derivative(x)
            4*(x^2 + y)*x
            sage: t.derivative(x, 2)
            12*x^2 + 4*y
            sage: t.derivative(x, 2, y)
            4
            sage: t.derivative(y)
            2*x^2 + 2*y

        ::

            sage: t = sin(x+y^2)*tan(x*y)
            sage: t.derivative(x)
            (tan(x*y)^2 + 1)*y*sin(y^2 + x) + cos(y^2 + x)*tan(x*y)
            sage: t.derivative(y)
            (tan(x*y)^2 + 1)*x*sin(y^2 + x) + 2*y*cos(y^2 + x)*tan(x*y)

        ::

            sage: h = sin(x)/cos(x)
            sage: derivative(h,x,x,x)
            6*sin(x)^4/cos(x)^4 + 8*sin(x)^2/cos(x)^2 + 2
            sage: derivative(h,x,3)
            6*sin(x)^4/cos(x)^4 + 8*sin(x)^2/cos(x)^2 + 2

        ::

            sage: var('x, y')
            (x, y)
            sage: u = (sin(x) + cos(y))*(cos(x) - sin(y))
            sage: derivative(u,x,y)
            sin(x)*sin(y) - cos(x)*cos(y)
            sage: f = ((x^2+1)/(x^2-1))^(1/4)
            sage: g = derivative(f, x); g # this is a complex expression
            1/2*(x/(x^2 - 1) - (x^2 + 1)*x/(x^2 - 1)^2)/((x^2 + 1)/(x^2 - 1))^(3/4)
            sage: g.factor()
            -x/((x^2 - 1)^(5/4)*(x^2 + 1)^(3/4))

        ::

            sage: y = var('y')
            sage: f = y^(sin(x))
            sage: derivative(f, x)
            y^sin(x)*log(y)*cos(x)

        ::

            sage: g(x) = sqrt(5-2*x)
            sage: g_3 = derivative(g, x, 3); g_3(2)
            -3

        ::

            sage: f = x*e^(-x)
            sage: derivative(f, 100)
            x*e^(-x) - 100*e^(-x)

        ::

            sage: g = 1/(sqrt((x^2-1)*(x+5)^6))
            sage: derivative(g, x)
            -((x + 5)^6*x + 3*(x + 5)^5*(x^2 - 1))/((x + 5)^6*(x^2 - 1))^(3/2)

        TESTS::

            sage: t.derivative()
            Traceback (most recent call last):
            ...
            ValueError: No differentiation variable specified.
        """
        return multi_derivative(self, args)

    diff = differentiate = derivative

    def _derivative(self, symb=None, deg=1):
        """
        Return the deg-th (partial) derivative of self with respect to symb.

        EXAMPLES::

            sage: var("x y")
            (x, y)
            sage: b = (x+y)^5
            sage: b._derivative(x, 2)
            20*(x + y)^3

            sage: from sage.symbolic.function import function as myfunc
            sage: foo = myfunc('foo',2)
            sage: foo(x^2,x^2)._derivative(x)
            2*x*D[0](foo)(x^2, x^2) + 2*x*D[1](foo)(x^2, x^2)

            sage: SR(1)._derivative()
            0

        TESTS:

        Raise error if no variable is specified and there are multiple
        variables::

            sage: b._derivative()
            Traceback (most recent call last):
            ...
            ValueError: No differentiation variable specified.
        """
        if symb is None:
            # we specify a default value of None for symb and check for it here
            # to return more helpful error messages when no variable is
            # given by the multi_derivative framework
            vars = self.variables()
            if len(vars) == 1:
                symb = vars[0]
            elif len(vars) == 0:
                return self._parent(0)
            else:
                raise ValueError, "No differentiation variable specified."
        if not isinstance(deg, (int, long, sage.rings.integer.Integer)) \
                or deg < 1:
            raise TypeError, "argument deg should be an integer >1."
        cdef Expression symbol = self.coerce_in(symb)
        if not is_a_symbol(symbol._gobj):
            raise TypeError, "argument symb must be a symbol"
        _sig_on
        cdef GEx x = self._gobj.diff(ex_to_symbol(symbol._gobj), deg)
        _sig_off
        return new_Expression_from_GEx(self._parent, x)

    def gradient(self, variables=None):
        r"""
        Compute the gradient of a symbolic function.

        This function returns a vector whose components are the derivatives
        of the original function with respect to the arguments of the
        original function. Alternatively, you can specify the variables as
        a list.

        EXAMPLES::

            sage: x,y = var('x y')
            sage: f = x^2+y^2
            sage: f.gradient()
            (2*x, 2*y)
            sage: g(x,y) = x^2+y^2
            sage: g.gradient()
            ((x, y) |--> 2*x, (x, y) |--> 2*y)
            sage: n = var('n')
            sage: f(x,y) = x^n+y^n
            sage: f.gradient()
            ((x, y) |--> n*x^(n - 1), (x, y) |--> n*y^(n - 1))
            sage: f.gradient([y,x])
            ((x, y) |--> n*y^(n - 1), (x, y) |--> n*x^(n - 1))
        """
        from sage.modules.free_module_element import vector
        if variables is None:
            variables = self.arguments()
        return vector([self.derivative(x) for x in variables])

    def hessian(self):
        r"""
        Compute the hessian of a function. This returns a matrix components
        are the 2nd partial derivatives of the original function.

        EXAMPLES::

            sage: x,y = var('x y')
            sage: f = x^2+y^2
            sage: f.hessian()
            [2 0]
            [0 2]
            sage: g(x,y) = x^2+y^2
            sage: g.hessian()
            [(x, y) |--> 2 (x, y) |--> 0]
            [(x, y) |--> 0 (x, y) |--> 2]
        """
        from sage.matrix.constructor import matrix
        return matrix([[g.derivative(x) for x in self.arguments()]
                       for g in self.gradient()])


    def series(self, symbol, int order):
        r"""
        Return the power series expansion of self in terms of the variable
        symbol to the given order.

        INPUT:

        - ``symbol`` - a variable

        - ``order`` - an integer

        OUTPUT:

        - a power series

        To truncate the power series and obtain a normal expression, use the
        truncate command.

        EXAMPLES:

        We expand a polynomial in `x` about 0, about `1`, and also truncate
        it back to a polynomial::

            sage: var('x,y')
            (x, y)
            sage: f = (x^3 - sin(y)*x^2 - 5*x + 3); f
            x^3 - x^2*sin(y) - 5*x + 3
            sage: g = f.series(x, 4); g
            3 + (-5)*x + (-sin(y))*x^2 + 1*x^3
            sage: g.truncate()
            x^3 - x^2*sin(y) - 5*x + 3
            sage: g = f.series(x==1, 4); g
            (-sin(y) - 1) + (-2*sin(y) - 2)*(x - 1) + (-sin(y) + 3)*(x - 1)^2 + 1*(x - 1)^3
            sage: h = g.truncate(); h
            -(sin(y) - 3)*(x - 1)^2 + (x - 1)^3 - 2*(sin(y) + 1)*(x - 1) - sin(y) - 1
            sage: h.expand()
            x^3 - x^2*sin(y) - 5*x + 3

        We computer another series expansion of an analytic function::

            sage: f = sin(x)/x^2
            sage: f.series(x,7)
            1*x^(-1) + (-1/6)*x + 1/120*x^3 + (-1/5040)*x^5 + Order(x^7)
            sage: f.series(x==1,3)
            (sin(1)) + (-2*sin(1) + cos(1))*(x - 1) + (5/2*sin(1) - 2*cos(1))*(x - 1)^2 + Order((x - 1)^3)
            sage: f.series(x==1,3).truncate().expand()
            5/2*x^2*sin(1) - 2*x^2*cos(1) - 7*x*sin(1) + 5*x*cos(1) + 11/2*sin(1) - 3*cos(1)

        Following the GiNaC tutorial, we use John Machin's amazing
        formula `\pi = 16 \tan^{-1}(1/5) - 4 \tan^{-1}(1/239)` to compute
        digits of `\pi`. We expand the arc tangent around 0 and insert
        the fractions 1/5 and 1/239.

        ::

            sage: x = var('x')
            sage: f = atan(x).series(x, 10); f
            1*x + (-1/3)*x^3 + 1/5*x^5 + (-1/7)*x^7 + 1/9*x^9 + Order(x^10)
            sage: float(16*f.subs(x==1/5) - 4*f.subs(x==1/239))
            3.1415926824043994
        """
        cdef Expression symbol0 = self.coerce_in(symbol)
        _sig_on
        cdef GEx x = self._gobj.series(symbol0._gobj, order, 0)
        _sig_off
        return new_Expression_from_GEx(self._parent, x)

    def taylor(self, v, a, n):
        r"""
        Expands this symbolic expression in a truncated Taylor or
        Laurent series in the variable `v` around the point `a`,
        containing terms through `(x - a)^n`.

        INPUT:

        -  ``v`` - variable
        -  ``a`` - number
        -  ``n`` - integer

        EXAMPLES::

            sage: var('a, x, z')
            (a, x, z)
            sage: taylor(a*log(z), z, 2, 3)
            1/24*(z - 2)^3*a - 1/8*(z - 2)^2*a + 1/2*(z - 2)*a + a*log(2)
            sage: taylor(sqrt (sin(x) + a*x + 1), x, 0, 3)
            1/48*(3*a^3 + 9*a^2 + 9*a - 1)*x^3 - 1/8*(a^2 + 2*a + 1)*x^2 + 1/2*(a + 1)*x + 1
            sage: taylor (sqrt (x + 1), x, 0, 5)
            7/256*x^5 - 5/128*x^4 + 1/16*x^3 - 1/8*x^2 + 1/2*x + 1
            sage: taylor (1/log (x + 1), x, 0, 3)
            -19/720*x^3 + 1/24*x^2 - 1/12*x + 1/x + 1/2
            sage: taylor (cos(x) - sec(x), x, 0, 5)
            -1/6*x^4 - x^2
            sage: taylor ((cos(x) - sec(x))^3, x, 0, 9)
            -1/2*x^8 - x^6
            sage: taylor (1/(cos(x) - sec(x))^3, x, 0, 5)
            -15377/7983360*x^4 - 6767/604800*x^2 + 11/120/x^2 + 1/2/x^4 - 1/x^6 - 347/15120
        """
        from sage.all import SR, Integer
        l = self._maxima_().taylor(v, SR(a), Integer(n))
        return self.parent()(l)


    def truncate(self):
        """
        Given a power series or expression, return the corresponding
        expression without the big oh.

        INPUT:

        - a series as output by the series command

        OUTPUT:

        - expression

        EXAMPLES::

            sage: f = sin(x)/x^2
            sage: f.truncate()
            sin(x)/x^2
            sage: f.series(x,7)
            1*x^(-1) + (-1/6)*x + 1/120*x^3 + (-1/5040)*x^5 + Order(x^7)
            sage: f.series(x,7).truncate()
            -1/5040*x^5 + 1/120*x^3 - 1/6*x + 1/x
            sage: f.series(x==1,3).truncate().expand()
            5/2*x^2*sin(1) - 2*x^2*cos(1) - 7*x*sin(1) + 5*x*cos(1) + 11/2*sin(1) - 3*cos(1)
        """
        if not is_a_series(self._gobj):
            return self
        return new_Expression_from_GEx(self._parent, series_to_poly(self._gobj))

    def expand(Expression self, side=None):
        """
        Expand this symbolic expression. Products of sums and exponentiated
        sums are multiplied out, numerators of rational expressions which
        are sums are split into their respective terms, and multiplications
        are distributed over addition at all levels.

        EXAMPLES:

        We expand the expression `(x-y)^5` using both
        method and functional notation.

        ::

            sage: x,y = var('x,y')
            sage: a = (x-y)^5
            sage: a.expand()
            x^5 - 5*x^4*y + 10*x^3*y^2 - 10*x^2*y^3 + 5*x*y^4 - y^5
            sage: expand(a)
            x^5 - 5*x^4*y + 10*x^3*y^2 - 10*x^2*y^3 + 5*x*y^4 - y^5

        We expand some other expressions::

            sage: expand((x-1)^3/(y-1))
            x^3/(y - 1) - 3*x^2/(y - 1) + 3*x/(y - 1) - 1/(y - 1)
            sage: expand((x+sin((x+y)^2))^2)
            x^2 + 2*x*sin((x + y)^2) + sin((x + y)^2)^2

        We can expand individual sides of a relation::

            sage: a = (16*x-13)^2 == (3*x+5)^2/2
            sage: a.expand()
            256*x^2 - 416*x + 169 == 9/2*x^2 + 15*x + 25/2
            sage: a.expand('left')
            256*x^2 - 416*x + 169 == 1/2*(3*x + 5)^2
            sage: a.expand('right')
            (16*x - 13)^2 == 9/2*x^2 + 15*x + 25/2

        TESTS:

            sage: var('x,y')
            (x, y)
            sage: ((x + (2/3)*y)^3).expand()
            x^3 + 2*x^2*y + 4/3*x*y^2 + 8/27*y^3
            sage: expand( (x*sin(x) - cos(y)/x)^2 )
            x^2*sin(x)^2 - 2*sin(x)*cos(y) + cos(y)^2/x^2
            sage: f = (x-y)*(x+y); f
            (x - y)*(x + y)
            sage: f.expand()
            x^2 - y^2
        """
        if side is not None:
            if not is_a_relational(self._gobj):
                raise ValueError, "expansion on sides only makes sense for relations"
            if side == 'left':
                return self.operator()(self.lhs().expand(), self.rhs())
            elif side == 'right':
                return self.operator()(self.lhs(), self.rhs().expand())
            else:
                raise ValueError, "side must be 'left', 'right', or None"

        _sig_on
        cdef GEx x = self._gobj.expand(0)
        _sig_off
        return new_Expression_from_GEx(self._parent, x)

    expand_rational = rational_expand = expand

    def expand_trig(self, full=False, half_angles=False, plus=True, times=True):
        """
        Expands trigonometric and hyperbolic functions of sums of angles
        and of multiple angles occurring in self. For best results, self
        should already be expanded.

        INPUT:

        -  ``full`` - (default: False) To enhance user control
           of simplification, this function expands only one level at a time
           by default, expanding sums of angles or multiple angles. To obtain
           full expansion into sines and cosines immediately, set the optional
           parameter full to True.

        -  ``half_angles`` - (default: False) If True, causes
           half-angles to be simplified away.

        -  ``plus`` - (default: True) Controls the sum rule;
           expansion of sums (e.g. 'sin(x + y)') will take place only if plus
           is True.

        -  ``times`` - (default: True) Controls the product
           rule, expansion of products (e.g. sin(2\*x)) will take place only
           if times is True.


        OUTPUT: a symbolic expression

        EXAMPLES::

            sage: sin(5*x).expand_trig()
            sin(x)^5 - 10*sin(x)^3*cos(x)^2 + 5*sin(x)*cos(x)^4
            sage: cos(2*x + var('y')).expand_trig()
            -sin(2*x)*sin(y) + cos(2*x)*cos(y)

        We illustrate various options to this function::

            sage: f = sin(sin(3*cos(2*x))*x)
            sage: f.expand_trig()
            sin(-(sin(cos(2*x))^3 - 3*sin(cos(2*x))*cos(cos(2*x))^2)*x)
            sage: f.expand_trig(full=True)
            sin(((sin(sin(x)^2)*cos(cos(x)^2) - sin(cos(x)^2)*cos(sin(x)^2))^3 - 3*(sin(sin(x)^2)*cos(cos(x)^2) - sin(cos(x)^2)*cos(sin(x)^2))*(sin(sin(x)^2)*sin(cos(x)^2) + cos(sin(x)^2)*cos(cos(x)^2))^2)*x)
            sage: sin(2*x).expand_trig(times=False)
            sin(2*x)
            sage: sin(2*x).expand_trig(times=True)
            2*sin(x)*cos(x)
            sage: sin(2 + x).expand_trig(plus=False)
            sin(x + 2)
            sage: sin(2 + x).expand_trig(plus=True)
            sin(2)*cos(x) + sin(x)*cos(2)
            sage: sin(x/2).expand_trig(half_angles=False)
            sin(1/2*x)
            sage: sin(x/2).expand_trig(half_angles=True)
            1/2*sqrt(-cos(x) + 1)*sqrt(2)*(-1)^floor(1/2*x/pi)

        ALIASES:

        :meth:`trig_expand` and :meth:`expand_trig` are the same
        """
        from sage.calculus.calculus import maxima_options
        M = self._maxima_()
        P = M.parent()
        opt = maxima_options(trigexpand=full, halfangles=half_angles,
                             trigexpandplus=plus, trigexpandtimes=times)
        cmd = 'trigexpand(%s), %s'%(M.name(), opt)
        ans = P(cmd)
        return self.parent()(ans)

    trig_expand = expand_trig

    ############################################################################
    # Pattern Matching
    ############################################################################
    def match(self, pattern):
        """
        Check if self matches the given pattern.

        INPUT:

        -  ``pattern`` - a symbolic expression, possibly containing wildcards
           to match for

        OUTPUT:

        - None - if there is no match
        - a dictionary mapping the wildcards to the matching values if a match
          was found. Note that the dictionary is empty if there were no
          wildcards in the given pattern.

        See also http://www.ginac.de/tutorial/Pattern-matching-and-advanced-substitutions.html

        EXAMPLES::

            sage: var('x,y,z,a,b,c,d,e,f')
            (x, y, z, a, b, c, d, e, f)
            sage: w0 = SR.wild(0); w1 = SR.wild(1); w2 = SR.wild(2)
            sage: ((x+y)^a).match((x+y)^a)  # no wildcards, so empty dict
            {}
            sage: print ((x+y)^a).match((x+y)^b)
            None
            sage: t = ((x+y)^a).match(w0^w1)
            sage: t[w0], t[w1]
            (x + y, a)
            sage: print ((x+y)^a).match(w0^w0)
            None
            sage: ((x+y)^(x+y)).match(w0^w0)
            {$0: x + y}
            sage: t = ((a+b)*(a+c)).match((a+w0)*(a+w1))
            sage: t[w0], t[w1]
            (b, c)
            sage: ((a+b)*(a+c)).match((w0+b)*(w0+c))
            {$0: a}
            sage: print ((a+b)*(a+c)).match((w0+w1)*(w0+w2))    # surprising?
            None
            sage: t = (a*(x+y)+a*z+b).match(a*w0+w1)
            sage: t[w0], t[w1]
            (x + y, a*z + b)
            sage: print (a+b+c+d+e+f).match(c)
            None
            sage: (a+b+c+d+e+f).has(c)
            True
            sage: (a+b+c+d+e+f).match(c+w0)
            {$0: a + b + d + e + f}
            sage: (a+b+c+d+e+f).match(c+e+w0)
            {$0: a + b + d + f}
            sage: (a+b).match(a+b+w0)
            {$0: 0}
            sage: print (a*b^2).match(a^w0*b^w1)
            None
            sage: (a*b^2).match(a*b^w1)
            {$1: 2}
            sage: (x*x.arctan2(x^2)).match(w0*w0.arctan2(w0^2))
            {$0: x}

        Beware that behind-the-scenes simplification can lead to
        surprising results in matching::

            sage: print (x+x).match(w0+w1)
            None
            sage: t = x+x; t
            2*x
            sage: t.operator()
            <built-in function mul>

        Since asking to match w0+w1 looks for an addition operator,
        there is no match.
        """
        cdef Expression p = self.coerce_in(pattern)
        cdef GExList mlst
        cdef bint res = self._gobj.match(p._gobj, mlst)
        if not res:
            return None

        cdef dict rdict = {}
        cdef GExListIter itr = mlst.begin()
        cdef GExListIter lstend = mlst.end()
        while itr.is_not_equal(lstend):
            key = new_Expression_from_GEx(self._parent, itr.obj().lhs())
            val = new_Expression_from_GEx(self._parent, itr.obj().rhs())
            rdict[key] = val
            itr.inc()
        return rdict


    def find(self, pattern):
        """
        Find all occurrences of the given pattern in this expression.

        Note that once a subexpression matches the pattern, the search doesn't
        extend to subexpressions of it.

        EXAMPLES::

            sage: var('x,y,z,a,b')
            (x, y, z, a, b)
            sage: w0 = SR.wild(0); w1 = SR.wild(1)

            sage: (sin(x)*sin(y)).find(sin(w0))
            [sin(x), sin(y)]

            sage: ((sin(x)+sin(y))*(a+b)).expand().find(sin(w0))
            [sin(x), sin(y)]

            sage: (1+x+x^2+x^3).find(x)
            [x]
            sage: (1+x+x^2+x^3).find(x^w0)
            [x^3, x^2]

            sage: (1+x+x^2+x^3).find(y)
            []

            # subexpressions of a match are not listed
            sage: ((x^y)^z).find(w0^w1)
            [(x^y)^z]
        """
        cdef Expression p = self.coerce_in(pattern)
        cdef GExList found
        self._gobj.find(p._gobj, found)
        res = []
        cdef GExListIter itr = found.begin()
        while itr.is_not_equal(found.end()):
            res.append(new_Expression_from_GEx(self._parent, itr.obj()))
            itr.inc()
        return res

    def has(self, pattern):
        """
        EXAMPLES::

            sage: var('x,y,a'); w0 = SR.wild(); w1 = SR.wild()
            (x, y, a)
            sage: (x*sin(x + y + 2*a)).has(y)
            True

        Here "x+y" is not a subexpression of "x+y+2*a" (which has the
        subexpressions "x", "y" and "2*a")::

            sage: (x*sin(x + y + 2*a)).has(x+y)
            False
            sage: (x*sin(x + y + 2*a)).has(x + y + w0)
            True

        The following fails because "2*(x+y)" automatically gets converted to
        "2*x+2*y" of which "x+y" is not a subexpression::

            sage: (x*sin(2*(x+y) + 2*a)).has(x+y)
            False

        Although x^1==x and x^0==1, neither "x" nor "1" are actually of the
        form "x^something"::

            sage: (x+1).has(x^w0)
            False

        Here is another possible pitfall, where the first expression
        matches because the term "-x" has the form "(-1)*x" in GiNaC. To check
        whether a polynomial contains a linear term you should use the
        coeff() function instead.

        ::

            sage: (4*x^2 - x + 3).has(w0*x)
            True
            sage: (4*x^2 + x + 3).has(w0*x)
            False
            sage: (4*x^2 + x + 3).has(x)
            True
            sage: (4*x^2 - x + 3).coeff(x,1)
            -1
            sage: (4*x^2 + x + 3).coeff(x,1)
            1
        """
        cdef Expression p = self.coerce_in(pattern)
        return self._gobj.has(p._gobj)

    def substitute(self, in_dict=None, **kwds):
        """
        EXAMPLES::

            sage: var('x,y,z,a,b,c,d,e,f')
            (x, y, z, a, b, c, d, e, f)
            sage: w0 = SR.wild(0); w1 = SR.wild(1)
            sage: t = a^2 + b^2 + (x+y)^3

            # substitute with keyword arguments (works only with symbols)
            sage: t.subs(a=c)
            (x + y)^3 + b^2 + c^2

            # substitute with a dictionary argument
            sage: t.subs({a^2: c})
            (x + y)^3 + b^2 + c

            sage: t.subs({w0^2: w0^3})
            (x + y)^3 + a^3 + b^3

            # substitute with a relational expression
            sage: t.subs(w0^2 == w0^3)
            (x + y)^3 + a^3 + b^3

            sage: t.subs(w0==w0^2)
            (x^2 + y^2)^18 + a^16 + b^16

            # more than one keyword argument is accepted
            sage: t.subs(a=b, b=c)
            (x + y)^3 + b^2 + c^2

            # using keyword arguments with a dictionary is allowed
            sage: t.subs({a:b}, b=c)
            (x + y)^3 + b^2 + c^2

            # in this case keyword arguments override the dictionary
            sage: t.subs({a:b}, a=c)
            (x + y)^3 + b^2 + c^2

            sage: t.subs({a:b, b:c})
            (x + y)^3 + b^2 + c^2

        TESTS:
            # no arguments return the same expression
            sage: t.subs()
            (x + y)^3 + a^2 + b^2

            # similarly for an empty dictionary argument
            sage: t.subs({})
            (x + y)^3 + a^2 + b^2

            # non keyword or dictionary argument returns error
            sage: t.subs(5)
            Traceback (most recent call last):
            ...
            TypeError: subs takes either a set of keyword arguments, a dictionary, or a symbolic relational expression

            # substitutions with infinity
            sage: (x/y).subs(y=oo)
            0
            sage: (x/y).subs(x=oo)
            +Infinity
            sage: (x*y).subs(x=oo)
            +Infinity
            sage: (x^y).subs(x=oo)
            Traceback (most recent call last):
            ...
            RuntimeError: power::eval(): pow(Infinity, x) for non numeric x is not defined.
            sage: (x^y).subs(y=oo)
            Traceback (most recent call last):
            ...
            RuntimeError: power::eval(): pow(x, Infinity) for non numeric x is not defined.
            sage: (x+y).subs(x=oo)
            +Infinity
            sage: (x-y).subs(y=oo)
            -Infinity
            sage: gamma(x).subs(x=-1)
            Infinity
            sage: 1/gamma(x).subs(x=-1)
            0

            # verify that this operation does not modify the passed dictionary (#6622)
            sage: var('v t')
            (v, t)
            sage: f = v*t
            sage: D = {v: 2}
            sage: f(D, t=3)
            6
            sage: D
            {v: 2}
        """
        cdef dict sdict = {}
        if in_dict is not None:
            if isinstance(in_dict, Expression):
                return self._subs_expr(in_dict)
            if not isinstance(in_dict, dict):
                raise TypeError, "subs takes either a set of keyword arguments, a dictionary, or a symbolic relational expression"
            sdict.update(in_dict)

        if kwds:
            for k, v in kwds.iteritems():
                k = self._parent.var(k)
                sdict[k] = v

        cdef GExMap smap
        for k, v in sdict.iteritems():
            smap.insert(make_pair((<Expression>self.coerce_in(k))._gobj,
                                  (<Expression>self.coerce_in(v))._gobj))

        return new_Expression_from_GEx(self._parent, self._gobj.subs_map(smap))

    subs = substitute

    cpdef Expression _subs_expr(self, expr):
        """
        EXAMPLES::

            sage: var('x,y,z,a,b,c,d,e,f')
            (x, y, z, a, b, c, d, e, f)
            sage: w0 = SR.wild(0); w1 = SR.wild(1)
            sage: (a^2 + b^2 + (x+y)^2)._subs_expr(w0^2 == w0^3)
            (x + y)^3 + a^3 + b^3
            sage: (a^4 + b^4 + (x+y)^4)._subs_expr(w0^2 == w0^3)
            (x + y)^4 + a^4 + b^4
            sage: (a^2 + b^4 + (x+y)^4)._subs_expr(w0^2 == w0^3)
            (x + y)^4 + a^3 + b^4
            sage: ((a+b+c)^2)._subs_expr(a+b == x)
            (a + b + c)^2
            sage: ((a+b+c)^2)._subs_expr(a+b+w0 == x+w0)
            (c + x)^2
            sage: (a+2*b)._subs_expr(a+b == x)
            a + 2*b
            sage: (a+2*b)._subs_expr(a+b+w0 == x+w0)
            a + 2*b
            sage: (a+2*b)._subs_expr(a+w0*b == x)
            x
            sage: (a+2*b)._subs_expr(a+b+w0*b == x+w0*b)
            a + 2*b
            sage: (4*x^3-2*x^2+5*x-1)._subs_expr(x==a)
            4*a^3 - 2*a^2 + 5*a - 1
            sage: (4*x^3-2*x^2+5*x-1)._subs_expr(x^w0==a^w0)
            4*a^3 - 2*a^2 + 5*x - 1
            sage: (4*x^3-2*x^2+5*x-1)._subs_expr(x^w0==a^(2*w0))._subs_expr(x==a)
            4*a^6 - 2*a^4 + 5*a - 1
            sage: sin(1+sin(x))._subs_expr(sin(w0)==cos(w0))
            cos(cos(x) + 1)
            sage: (sin(x)^2 + cos(x)^2)._subs_expr(sin(w0)^2+cos(w0)^2==1)
            1
            sage: (1 + sin(x)^2 + cos(x)^2)._subs_expr(sin(w0)^2+cos(w0)^2==1)
            sin(x)^2 + cos(x)^2 + 1
            sage: (17*x + sin(x)^2 + cos(x)^2)._subs_expr(w1 + sin(w0)^2+cos(w0)^2 == w1 + 1)
            17*x + 1
            sage: ((x-1)*(sin(x)^2 + cos(x)^2)^2)._subs_expr(sin(w0)^2+cos(w0)^2 == 1)
            x - 1
            """
        cdef Expression p = self.coerce_in(expr)
        return new_Expression_from_GEx(self._parent, self._gobj.subs(p._gobj))

    def substitute_expression(self, *equations):
        """
        Given a dictionary of key:value pairs, substitute all occurrences
        of key for value in self.  The substitutions can also be given
        as a number of symbolic equalities key == value; see the
        examples.

        .. warning::

           This is a formal pattern substitution, which may or may not
           have any mathematical meaning. The exact rules used at
           present in Sage are determined by Maxima's subst
           command. Sometimes patterns are not replaced even though
           one would think they should be - see examples below.

        EXAMPLES::

            sage: f = x^2 + 1
            sage: f.subs_expr(x^2 == x)
            x + 1

        ::

            sage: var('x,y,z'); f = x^3 + y^2 + z
            (x, y, z)
            sage: f.subs_expr(x^3 == y^2, z == 1)
            2*y^2 + 1

        Or the same thing giving the substitutions as a dictionary::

            sage: f.subs_expr({x^3:y^2, z:1})
            2*y^2 + 1

            sage: f = x^2 + x^4
            sage: f.subs_expr(x^2 == x)
            x^4 + x
            sage: f = cos(x^2) + sin(x^2)
            sage: f.subs_expr(x^2 == x)
            sin(x) + cos(x)

        ::

            sage: f(x,y,t) = cos(x) + sin(y) + x^2 + y^2 + t
            sage: f.subs_expr(y^2 == t)
            (x, y, t) |--> x^2 + 2*t + sin(y) + cos(x)

        The following seems really weird, but it *is* what Maple does::

            sage: f.subs_expr(x^2 + y^2 == t)
            (x, y, t) |--> x^2 + y^2 + t + sin(y) + cos(x)
            sage: maple.eval('subs(x^2 + y^2 = t, cos(x) + sin(y) + x^2 + y^2 + t)')          # optional requires maple
            'cos(x)+sin(y)+x^2+y^2+t'
            sage: maxima.quit()
            sage: maxima.eval('cos(x) + sin(y) + x^2 + y^2 + t, x^2 + y^2 = t')
            'sin(y)+y^2+cos(x)+x^2+t'

        Actually Mathematica does something that makes more sense::

            sage: mathematica.eval('Cos[x] + Sin[y] + x^2 + y^2 + t /. x^2 + y^2 -> t')       # optional -- requires mathematica
            2 t + Cos[x] + Sin[y]
        """
        if isinstance(equations[0], dict):
            eq_dict = equations[0]
            equations = [ x == eq_dict[x] for x in eq_dict.keys() ]

        if not all([is_SymbolicEquation(eq) for eq in equations]):
            raise TypeError, "each expression must be an equation"

        d = dict([(eq.lhs(), eq.rhs()) for eq in equations])
        return self.subs(d)

    subs_expr = substitute_expression

    def substitute_function(self, original, new):
        """
        Returns this symbolic expressions all occurrences of the
        function *original* replaced with the function *new*.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: clear_functions()
            sage: foo = function('foo'); bar = function('bar')
            sage: f = foo(x) + 1/foo(pi*y)
            sage: f.substitute_function(foo, bar)
            1/bar(pi*y) + bar(x)
        """
        from sage.symbolic.expression_conversions import SubstituteFunction
        return SubstituteFunction(self, original, new)()

    def __call__(self, *args, **kwds):
        """
        Calls the :meth:`subs` on this expression.

        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: (x+y)(x=z^2, y=x^y)
            x^y + z^2
        """
        return self._parent._call_element_(self, *args, **kwds)

    def variables(self):
        """
        Return sorted list of variables that occur in this expression.

        EXAMPLES::

            sage: (x,y,z) = var('x,y,z')
            sage: (x+y).variables()
            (x, y)
            sage: (2*x).variables()
            (x,)
            sage: (x^y).variables()
            (x, y)
            sage: sin(x+y^z).variables()
            (x, y, z)

        """
        from sage.symbolic.ring import SR
        cdef GExSet sym_set
        g_list_symbols(self._gobj, sym_set)
        res = []
        cdef GExSetIter itr = sym_set.begin()
        while itr.is_not_equal(sym_set.end()):
            res.append(new_Expression_from_GEx(SR, itr.obj()))
            itr.inc()
        return tuple(res)

    def arguments(self):
        """
        EXAMPLES::

            sage: x,y = var('x,y')
            sage: f = x + y
            sage: f.arguments()
            (x, y)

            sage: g = f.function(x)
            sage: g.arguments()
            (x,)

        """
        try:
            return self._parent.arguments()
        except AttributeError:
            return self.variables()

    args = arguments

    def number_of_arguments(self):
        """
        EXAMPLES::

            sage: x,y = var('x,y')
            sage: f = x + y
            sage: f.number_of_arguments()
            2

            sage: g = f.function(x)
            sage: g.number_of_arguments()
            1

        ::

            sage: x,y,z = var('x,y,z')
            sage: (x+y).number_of_arguments()
            2
            sage: (x+1).number_of_arguments()
            1
            sage: (sin(x)+1).number_of_arguments()
            1
            sage: (sin(z)+x+y).number_of_arguments()
            3
            sage: (sin(x+y)).number_of_arguments()
            2

        ::

            sage: ( 2^(8/9) - 2^(1/9) )(x-1)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 0
        """
        return len(self.arguments())

    def number_of_operands(self):
        """
        Returns the number of arguments of this expression.

        EXAMPLES::

            sage: var('a,b,c,x,y')
            (a, b, c, x, y)
            sage: a.number_of_operands()
            0
            sage: (a^2 + b^2 + (x+y)^2).number_of_operands()
            3
            sage: (a^2).number_of_operands()
            2
            sage: (a*b^2*c).number_of_operands()
            3
        """
        return self._gobj.nops()

    nops = number_of_operands

    def __len__(self):
        """
        Returns the number of arguments of this expression.

        EXAMPLES::

            sage: var('a,b,c,x,y')
            (a, b, c, x, y)
            sage: len(a)
            0
            sage: len((a^2 + b^2 + (x+y)^2))
            3
            sage: len((a^2))
            2
            sage: len(a*b^2*c)
            3
        """
        return self.number_of_operands()

    def operands(self):
        """
        Returns a list containing the operands of this expression.

        EXAMPLES::

            sage: var('a,b,c,x,y')
            (a, b, c, x, y)
            sage: (a^2 + b^2 + (x+y)^2).operands()
            [(x + y)^2, a^2, b^2]
            sage: (a^2).operands()
            [a, 2]
            sage: (a*b^2*c).operands()
            [a, b^2, c]
        """
        from sage.symbolic.ring import SR
        return [new_Expression_from_GEx(SR, self._gobj.op(i)) \
                            for i from 0 <= i < self._gobj.nops()]

    def operator(self):
        """
        Returns the topmost operator in this expression.

        EXAMPLES::

            sage: x,y,z = var('x,y,z')
            sage: (x+y).operator()
            <built-in function add>
            sage: (x^y).operator()
            <built-in function pow>
            sage: (x^y * z).operator()
            <built-in function mul>
            sage: (x < y).operator()
            <built-in function lt>

            sage: abs(x).operator()
            abs
            sage: r = gamma(x).operator(); type(r)
            <class 'sage.functions.other.Function_gamma'>

            sage: from sage.symbolic.function import function
            sage: psi = function('psi', 1)
            sage: psi(x).operator()
            psi

            sage: r = psi(x).operator()
            sage: r == psi
            True

            sage: f = function('f', 1, conjugate_func=lambda x: 2*x)
            sage: nf = f(x).operator()
            sage: nf(x).conjugate()
            2*x

            sage: f = function('f')
            sage: a = f(x).diff(x); a
            D[0](f)(x)
            sage: a.operator()
            D[0](f)

        TESTS:
            sage: (x <= y).operator()
            <built-in function le>
            sage: (x == y).operator()
            <built-in function eq>
            sage: (x != y).operator()
            <built-in function ne>
            sage: (x > y).operator()
            <built-in function gt>
            sage: (x >= y).operator()
            <built-in function ge>
        """
        cdef operators o
        cdef unsigned serial
        import operator
        if is_a_add(self._gobj):
            return operator.add
        elif is_a_mul(self._gobj) or is_a_ncmul(self._gobj):
            return operator.mul
        elif is_a_power(self._gobj):
            return operator.pow
        elif is_a_relational(self._gobj):
            # find the operator and return it
            o = relational_operator(self._gobj)
            if o == equal:
                return operator.eq
            elif o == not_equal:
                return operator.ne
            elif o == less:
                return operator.lt
            elif o == less_or_equal:
                return operator.le
            elif o == greater:
                return operator.gt
            elif o == greater_or_equal:
                return operator.ge
            else:
                raise RuntimeError, "operator type not known, please report this as a bug"
        elif is_a_function(self._gobj):
            # get function id
            serial = ex_to_function(self._gobj).get_serial()

            from sage.symbolic.pynac import get_ginac_serial

            # if operator is a special function defined by us
            # find the python equivalent and return it
            res = get_sfunction_from_serial(serial)
            if res is None:
                raise RuntimeError, "cannot find SFunction in table"

            if is_a_fderivative(self._gobj):
                from sage.symbolic.pynac import paramset_from_Expression
                from sage.symbolic.operators import FDerivativeOperator
                parameter_set = paramset_from_Expression(self)
                res = FDerivativeOperator(res, parameter_set)

            return res

        # self._gobj is either a symbol, constant or numeric
        return None

    def __index__(self):
        """
        EXAMPLES::

            sage: a = range(10)
            sage: a[:SR(5)]
            [0, 1, 2, 3, 4]
        """
        return int(self._integer_())

    def iterator(self):
        """
        Return an iterator over the arguments of this expression.

        EXAMPLES::

            sage: x,y,z = var('x,y,z')
            sage: list((x+y+z).iterator())
            [x, y, z]
            sage: list((x*y*z).iterator())
            [x, y, z]
            sage: list((x^y*z*(x+y)).iterator())
            [x + y, x^y, z]
        """
        return new_ExpIter_from_Expression(self)

##     def __getitem__(self, ind):
##         """
##         EXAMPLES::

##             sage: x,y,z = var('x,y,z')
##             sage: e = x + x*y + z^y + 3*y*z; e
##             x*y + 3*y*z + z^y + x
##             sage: e[1]
##             3*y*z
##             sage: e[-1]
##             x
##             sage: e[1:]
##             [3*y*z, z^y, x]
##             sage: e[:2]
##             [x*y, 3*y*z]
##             sage: e[-2:]
##             [z^y, x]
##             sage: e[:-2]
##             [x*y, 3*y*z]

##         """
##         cdef int bind, eind, step, i
##         cdef int n_ops = self._gobj.nops()
##         if PY_TYPE_CHECK(ind, slice):
##             if ind.start:
##                 bind = ind.start
##                 if bind != ind.start:
##                     raise ValueError, "integer index expected"
##                 if bind < 0:
##                     bind = n_ops + bind
##             else:
##                 bind = 0
##             if ind.stop:
##                 eind = ind.stop
##                 if eind != ind.stop:
##                     raise ValueError, "integer index expected"
##                 if eind > n_ops:
##                     eind = n_ops
##                 if eind < 0:
##                     eind = n_ops + eind
##             else:
##                 eind = n_ops
##             if ind.step:
##                 step = ind.step
##                 if step != ind.step:
##                     raise ValueError, "step value must be an integer"
##             else:
##                 step = 1
##             return [new_Expression_from_GEx(self._parent, self._gobj.op(i))
##                     for i in xrange(bind, eind, step)]

##         try:
##             bind = ind
##             if bind != ind:
##                 raise ValueError, "integer index expected"
##         except TypeError:
##             raise TypeError, "index should either be a slice object, or an integer"
##         if bind < 0:
##             bind = n_ops + bind
##         return new_Expression_from_GEx(self._parent, self._gobj.op(bind))

    def n(self, prec=None, digits=None):
        """
        Return a numerical approximation this symbolic expression as
        either a real or complex number with at least the requested
        number of bits or digits of precision.

        EXAMPLES::

            sage: sin(x).subs(x=5).n()
            -0.958924274663138
            sage: sin(x).subs(x=5).n(100)
            -0.95892427466313846889315440616
            sage: sin(x).subs(x=5).n(digits=50)
            -0.95892427466313846889315440615599397335246154396460
            sage: zeta(x).subs(x=2).numerical_approx(digits=50)
            1.6449340668482264364724151666460251892189499012068

            sage: cos(3).numerical_approx(200)
            -0.98999249660044545727157279473126130239367909661558832881409
            sage: numerical_approx(cos(3), digits=10)
            -0.9899924966
            sage: (i + 1).numerical_approx(32)
            1.00000000 + 1.00000000*I
            sage: (pi + e + sqrt(2)).numerical_approx(100)
            7.2740880444219335226246195788

        TESTS:

        We test the evaluation of different infinities available in Pynac::

            sage: t = x - oo; t
            -Infinity
            sage: t.n()
            -infinity
            sage: t = x + oo; t
            +Infinity
            sage: t.n()
            +infinity
            sage: t = x - unsigned_infinity; t
            Infinity
            sage: t.n()
            +infinity
        """
        if prec is None:
            if digits is None:
                prec = 53
            else:
                prec = int((digits+1) * 3.32192) + 1
        x = new_Expression_from_GEx(self._parent, self._gobj.evalf(0, prec)).pyobject()
        # Important -- the  we get might not be a valid output for numerical_approx in
        # the case when one gets infinity.
        if isinstance(x, AnInfinity):
            return x.n(prec=prec,digits=digits)
        return x


    numerical_approx = n

    def function(self, *args):
        """
        Return a callable symbolic expression with the given variables.

        EXAMPLES:

        We will use several symbolic variables in the examples below::

            sage: var('x, y, z, t, a, w, n')
            (x, y, z, t, a, w, n)

        ::

            sage: u = sin(x) + x*cos(y)
            sage: g = u.function(x,y)
            sage: g(x,y)
            x*cos(y) + sin(x)
            sage: g(t,z)
            t*cos(z) + sin(t)
            sage: g(x^2, x^y)
            x^2*cos(x^y) + sin(x^2)

        ::

            sage: f = (x^2 + sin(a*w)).function(a,x,w); f
            (a, x, w) |--> x^2 + sin(a*w)
            sage: f(1,2,3)
            sin(3) + 4

        Using the :meth:`function` method we can obtain the above function
        `f`, but viewed as a function of different variables::

            sage: h = f.function(w,a); h
            (w, a) |--> x^2 + sin(a*w)

        This notation also works::

            sage: h(w,a) = f
            sage: h
            (w, a) |--> x^2 + sin(a*w)

        You can even make a symbolic expression `f` into a function
        by writing ``f(x,y) = f``::

            sage: f = x^n + y^n; f
            x^n + y^n
            sage: f(x,y) = f
            sage: f
            (x, y) |--> x^n + y^n
            sage: f(2,3)
            2^n + 3^n
        """
        # we override type checking in CallableSymbolicExpressionRing,
        # since it checks for old SymbolicVariable's
        # and do the check here instead
        from sage.symbolic.callable import CallableSymbolicExpressionRing
        from sage.symbolic.ring import is_SymbolicVariable
        for i in args:
            if not is_SymbolicVariable(i):
                break
        else:
            R = CallableSymbolicExpressionRing(args, check=False)
            return R(self)
        raise TypeError, "Must construct a function with a tuple (or list) of symbolic variables."

    ############################################################################
    # Polynomial functions
    ############################################################################
    def coefficient(self, s, int n=1):
        """
        Returns the coefficient of `s^n` in this symbolic expression.

        INPUT:

        - ``s`` - expression

        - ``n`` - integer, default 1

        OUTPUT:

        - coefficient of s^n

        Sometimes it may be necessary to expand or factor first, since this
        is not done automatically.

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + x/y + 2*sin(x*y)/x + 100
            sage: f.collect(x)
            x^3*sin(x*y) + (a + y + 1/y)*x + 2*sin(x*y)/x + 100
            sage: f.coefficient(x,0)
            100
            sage: f.coefficient(x,-1)
            2*sin(x*y)
            sage: f.coefficient(x,1)
            a + y + 1/y
            sage: f.coefficient(x,2)
            0
            sage: f.coefficient(x,3)
            sin(x*y)
            sage: f.coefficient(x^3)
            sin(x*y)
            sage: f.coefficient(sin(x*y))
            x^3 + 2/x
            sage: f.collect(sin(x*y))
            (x^3 + 2/x)*sin(x*y) + a*x + x*y + x/y + 100

            sage: var('a, x, y, z')
            (a, x, y, z)
            sage: f = (a*sqrt(2))*x^2 + sin(y)*x^(1/2) + z^z
            sage: f.coefficient(sin(y))
            sqrt(x)
            sage: f.coefficient(x^2)
            sqrt(2)*a
            sage: f.coefficient(x^(1/2))
            sin(y)
            sage: f.coefficient(1)
            0
            sage: f.coefficient(x, 0)
            sqrt(x)*sin(y) + z^z

        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._parent, self._gobj.coeff(ss._gobj, n))

    coeff = coefficient

    def coefficients(self, x=None):
        r"""
        Coefficients of this symbolic expression as a polynomial in x.

        INPUT:

        -  ``x`` - optional variable

        OUTPUT:

        - A list of pairs (expr, n), where expr is a symbolic
          expression and n is a power.

        EXAMPLES::

            sage: var('x, y, a')
            (x, y, a)
            sage: p = x^3 - (x-3)*(x^2+x) + 1
            sage: p.coefficients()
            [[1, 0], [3, 1], [2, 2]]
            sage: p = expand((x-a*sqrt(2))^2 + x + 1); p
            -2*sqrt(2)*a*x + 2*a^2 + x^2 + x + 1
            sage: p.coefficients(a)
            [[x^2 + x + 1, 0], [-2*sqrt(2)*x, 1], [2, 2]]
            sage: p.coefficients(x)
            [[2*a^2 + 1, 0], [-2*sqrt(2)*a + 1, 1], [1, 2]]

        A polynomial with wacky exponents::

            sage: p = (17/3*a)*x^(3/2) + x*y + 1/x + x^x
            sage: p.coefficients(x)
            [[1, -1], [x^x, 0], [y, 1], [17/3*a, 3/2]]
        """
        f = self._maxima_()
        maxima = f.parent()
        maxima._eval_line('load(coeflist)')
        if x is None:
            x = self.default_variable()
        x = self.parent().var(repr(x))
        G = f.coeffs(x)
        from sage.calculus.calculus import symbolic_expression_from_maxima_string
        S = symbolic_expression_from_maxima_string(repr(G))
        return S[1:]

    coeffs = coefficients

    def leading_coefficient(self, s):
        """
        Return the leading coefficient of s in self.

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + x/y + 2*sin(x*y)/x + 100
            sage: f.leading_coefficient(x)
            sin(x*y)
            sage: f.leading_coefficient(y)
            x
            sage: f.leading_coefficient(sin(x*y))
            x^3 + 2/x
        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._parent, self._gobj.lcoeff(ss._gobj))

    leading_coeff = leading_coefficient

    def trailing_coefficient(self, s):
        """
        Return the trailing coefficient of s in self, i.e., the coefficient
        of the smallest power of s in self.

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + x/y + 2*sin(x*y)/x + 100
            sage: f.trailing_coefficient(x)
            2*sin(x*y)
            sage: f.trailing_coefficient(y)
            x
            sage: f.trailing_coefficient(sin(x*y))
            a*x + x*y + x/y + 100
        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._parent, self._gobj.tcoeff(ss._gobj))

    trailing_coeff = trailing_coefficient

    def low_degree(self, s):
        """
        Return the exponent of the lowest nonpositive power of s in self.

        OUTPUT:

        - an integer <= 0.

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y^10 + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + 2*sin(x*y)/x + x/y^10 + 100
            sage: f.low_degree(x)
            -1
            sage: f.low_degree(y)
            -10
            sage: f.low_degree(sin(x*y))
            0
            sage: (x^3+y).low_degree(x)
            0
        """
        cdef Expression ss = self.coerce_in(s)
        return self._gobj.ldegree(ss._gobj)

    def degree(self, s):
        """
        Return the exponent of the highest nonnegative power of s in self.

        OUTPUT:

        - an integer >= 0.

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y^10 + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + 2*sin(x*y)/x + x/y^10 + 100
            sage: f.degree(x)
            3
            sage: f.degree(y)
            1
            sage: f.degree(sin(x*y))
            1
            sage: (x^-3+y).degree(x)
            0
        """
        cdef Expression ss = self.coerce_in(s)
        return self._gobj.degree(ss._gobj)

    def poly(self, x=None):
        r"""
        Express this symbolic expression as a polynomial in *x*. If
        this is not a polynomial in *x*, then some coefficients may be
        functions of *x*.

        .. warning::

           This is different from :meth:`polynomial` which returns
           a Sage polynomial over a given base ring.

        EXAMPLES::

            sage: var('a, x')
            (a, x)
            sage: p = expand((x-a*sqrt(2))^2 + x + 1); p
            -2*sqrt(2)*a*x + 2*a^2 + x^2 + x + 1
            sage: p.poly(a)
            -2*sqrt(2)*a*x + 2*a^2 + x^2 + x + 1
            sage: bool(expand(p.poly(a)) == p)
            True
            sage: p.poly(x)
            -(2*sqrt(2)*a - 1)*x + 2*a^2 + x^2 + 1
        """
        from sage.symbolic.ring import SR
        f = self._maxima_()
        P = f.parent()
        P._eval_line('load(coeflist)')
        if x is None:
            x = self.default_variable()
        x = self._parent.var(repr(x))
        G = f.coeffs(x)
        ans = None
        for i in range(1, len(G)):
            Z = G[i]
            coeff = SR(Z[0])
            n = SR(Z[1])
            if repr(coeff) != '0':
                if repr(n) == '0':
                    xpow = SR(1)
                elif repr(n) == '1':
                    xpow = x
                else:
                    xpow = x**n
                if ans is None:
                    ans = coeff*xpow
                else:
                    ans += coeff*xpow
        return ans

    def polynomial(self, base_ring, ring=None):
        r"""
        Return this symbolic expression as an algebraic polynomial
        over the given base ring, if possible.

        The point of this function is that it converts purely symbolic
        polynomials into optimised algebraic polynomials over a given
        base ring.

        .. warning::

           This is different from meth:`poly` which is used to rewrite
           self as a polynomial in terms of one of the variables.

        INPUT:

        -  ``base_ring`` - a ring

        EXAMPLES::

            sage: f = x^2 -2/3*x + 1
            sage: f.polynomial(QQ)
            x^2 - 2/3*x + 1
            sage: f.polynomial(GF(19))
            x^2 + 12*x + 1

        Polynomials can be useful for getting the coefficients of an
        expression::

            sage: g = 6*x^2 - 5
            sage: g.coefficients()
            [[-5, 0], [6, 2]]
            sage: g.polynomial(QQ).list()
            [-5, 0, 6]
            sage: g.polynomial(QQ).dict()
            {0: -5, 2: 6}

        ::

            sage: f = x^2*e + x + pi/e
            sage: f.polynomial(RDF)
            2.71828182846*x^2 + 1.0*x + 1.15572734979
            sage: g = f.polynomial(RR); g
            2.71828182845905*x^2 + 1.00000000000000*x + 1.15572734979092
            sage: g.parent()
            Univariate Polynomial Ring in x over Real Field with 53 bits of precision
            sage: f.polynomial(RealField(100))
            2.7182818284590452353602874714*x^2 + 1.0000000000000000000000000000*x + 1.1557273497909217179100931833
            sage: f.polynomial(CDF)
            2.71828182846*x^2 + 1.0*x + 1.15572734979
            sage: f.polynomial(CC)
            2.71828182845905*x^2 + 1.00000000000000*x + 1.15572734979092

        We coerce a multivariate polynomial with complex symbolic
        coefficients::

            sage: x, y, n = var('x, y, n')
            sage: f = pi^3*x - y^2*e - I; f
            pi^3*x - y^2*e - I
            sage: f.polynomial(CDF)
            (-2.71828182846)*y^2 + 31.0062766803*x - 1.0*I
            sage: f.polynomial(CC)
            (-2.71828182845905)*y^2 + 31.0062766802998*x - 1.00000000000000*I
            sage: f.polynomial(ComplexField(70))
            (-2.7182818284590452354)*y^2 + 31.006276680299820175*x - 1.0000000000000000000*I

        Another polynomial::

            sage: f = sum((e*I)^n*x^n for n in range(5)); f
            x^4*e^4 - I*x^3*e^3 - x^2*e^2 + I*x*e + 1
            sage: f.polynomial(CDF)
            54.5981500331*x^4 - 20.0855369232*I*x^3 - 7.38905609893*x^2 + 2.71828182846*I*x + 1.0
            sage: f.polynomial(CC)
            54.5981500331442*x^4 - 20.0855369231877*I*x^3 - 7.38905609893065*x^2 + 2.71828182845905*I*x + 1.00000000000000

        A multivariate polynomial over a finite field::

            sage: f = (3*x^5 - 5*y^5)^7; f
            (3*x^5 - 5*y^5)^7
            sage: g = f.polynomial(GF(7)); g
            3*x^35 + 2*y^35
            sage: parent(g)
            Multivariate Polynomial Ring in x, y over Finite Field of size 7
        """
        from sage.symbolic.expression_conversions import polynomial
        return polynomial(self, base_ring=base_ring, ring=ring)

    def _polynomial_(self, R):
        """
        Coerce this symbolic expression to a polynomial in `R`.

        EXAMPLES::

            sage: var('x,y,z,w')
            (x, y, z, w)

        ::

            sage: R = QQ[x,y,z]
            sage: R(x^2 + y)
            x^2 + y
            sage: R = QQ[w]
            sage: R(w^3 + w + 1)
            w^3 + w + 1
            sage: R = GF(7)[z]
            sage: R(z^3 + 10*z)
            z^3 + 3*z

        .. note::

           If the base ring of the polynomial ring is the symbolic ring,
           then a constant polynomial is always returned.

        ::

            sage: R = SR[x]
            sage: a = R(sqrt(2) + x^3 + y)
            sage: a
            y + sqrt(2) + x^3
            sage: type(a)
            <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field'>
            sage: a.degree()
            0

        We coerce to a double precision complex polynomial ring::

            sage: f = e*x^3 + pi*y^3 + sqrt(2) + I; f
            pi*y^3 + x^3*e + sqrt(2) + I
            sage: R = CDF[x,y]
            sage: R(f)
            2.71828182846*x^3 + 3.14159265359*y^3 + 1.41421356237 + 1.0*I

        We coerce to a higher-precision polynomial ring

        ::

            sage: R = ComplexField(100)[x,y]
            sage: R(f)
            2.7182818284590452353602874714*x^3 + 3.1415926535897932384626433833*y^3 + 1.4142135623730950488016887242 + 1.0000000000000000000000000000*I

        This shows that the issue at trac #4246 is fixed (attempting to
        coerce an expression containing at least one variable that's not in
        `R` raises an error)::

            sage: x, y = var('x y')
            sage: S = PolynomialRing(Integers(4), 1, 'x')
            sage: S(x)
            x
            sage: S(y)
            Traceback (most recent call last):
            ...
            TypeError: y is not a variable of Multivariate Polynomial Ring in x over Ring of integers modulo 4
            sage: S(x+y)
            Traceback (most recent call last):
            ...
            TypeError: y is not a variable of Multivariate Polynomial Ring in x over Ring of integers modulo 4
            sage: (x+y)._polynomial_(S)
            Traceback (most recent call last):
            ...
            TypeError: y is not a variable of Multivariate Polynomial Ring in x over Ring of integers modulo 4
        """
        from sage.symbolic.all import SR
        from sage.rings.all import is_MPolynomialRing

        base_ring = R.base_ring()
        if base_ring == SR:
            if is_MPolynomialRing(R):
                return R({tuple([0]*R.ngens()):self})
            else:
                return R([self])
        return self.polynomial(None, ring=R)

    def power_series(self, base_ring):
        """
        Return algebraic power series associated to this symbolic
        expression, which must be a polynomial in one variable, with
        coefficients coercible to the base ring.

        The power series is truncated one more than the degree.

        EXAMPLES::

            sage: theta = var('theta')
            sage: f = theta^3 + (1/3)*theta - 17/3
            sage: g = f.power_series(QQ); g
            -17/3 + 1/3*theta + theta^3 + O(theta^4)
            sage: g^3
            -4913/27 + 289/9*theta - 17/9*theta^2 + 2602/27*theta^3 + O(theta^4)
            sage: g.parent()
            Power Series Ring in theta over Rational Field
        """
        v = self.variables()
        if len(v) != 1:
            raise ValueError, "self must be a polynomial in one variable but it is in the variables %s"%tuple([v])
        f = self.polynomial(base_ring)
        from sage.rings.all import PowerSeriesRing
        R = PowerSeriesRing(base_ring, names=f.parent().variable_names())
        return R(f, f.degree()+1)

    def gcd(self, b):
        """
        Return the gcd of self and b, which must be integers or polynomials over
        the rational numbers.

        TODO: I tried the massive gcd from
        http://trac.sagemath.org/sage_trac/ticket/694 on Ginac dies
        after about 10 seconds.  Singular easily does that GCD now.
        Since Ginac only handles poly gcd over QQ, we should change
        ginac itself to use Singular.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: SR(10).gcd(SR(15))
            5
            sage: (x^3 - 1).gcd(x-1)
            x - 1
            sage: (x^3 - 1).gcd(x^2+x+1)
            x^2 + x + 1
            sage: (x^3 - sage.symbolic.constants.pi).gcd(x-sage.symbolic.constants.pi)
            Traceback (most recent call last):
            ...
            RuntimeError: gcd: arguments must be polynomials over the rationals
            sage: gcd(x^3 - y^3, x-y)
            -x + y
            sage: gcd(x^100-y^100, x^10-y^10)
            -x^10 + y^10
            sage: gcd(expand( (x^2+17*x+3/7*y)*(x^5 - 17*y + 2/3) ), expand((x^13+17*x+3/7*y)*(x^5 - 17*y + 2/3)) )
            1/7*x^5 - 17/7*y + 2/21
        """
        cdef Expression r = self.coerce_in(b)
        _sig_on
        cdef GEx x = g_gcd(self._gobj, r._gobj)
        _sig_off
        return new_Expression_from_GEx(self._parent, x)

    def collect(Expression self, s):
        """
        INPUT:

        - ``s`` - a symbol

        OUTPUT:

        - expression

        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: f = 4*x*y + x*z + 20*y^2 + 21*y*z + 4*z^2 + x^2*y^2*z^2
            sage: f.collect(x)
            x^2*y^2*z^2 + (4*y + z)*x + 20*y^2 + 21*y*z + 4*z^2
            sage: f.collect(y)
            (x^2*z^2 + 20)*y^2 + (4*x + 21*z)*y + x*z + 4*z^2
            sage: f.collect(z)
            (x^2*y^2 + 4)*z^2 + (x + 21*y)*z + 4*x*y + 20*y^2
        """
        cdef Expression s0 = self.coerce_in(s)
        _sig_on
        cdef GEx x = self._gobj.collect(s0._gobj, False)
        _sig_off
        return new_Expression_from_GEx(self._parent, x)

    def collect_common_factors(self):
        """
        EXAMPLES::

            sage: var('x')
            x
            sage: (x/(x^2 + x)).collect_common_factors()
            1/(x + 1)
        """
        _sig_on
        cdef GEx x = g_collect_common_factors(self._gobj)
        _sig_off
        return new_Expression_from_GEx(self._parent, x)

    def __abs__(self):
        """
        Return the absolute value of this expression.

        EXAMPLES::

            sage: var('x, y')
            (x, y)

        The absolute value of a symbolic expression::

            sage: abs(x^2+y^2)
            abs(x^2 + y^2)

        The absolute value of a number in the symbolic ring::

            sage: abs(SR(-5))
            5
            sage: type(abs(SR(-5)))
            <type 'sage.symbolic.expression.Expression'>
        """
        return new_Expression_from_GEx(self._parent, g_abs(self._gobj))

    def step(self):
        """
        Return the value of the Heaviside step function, which is 0 for
        negative x, 1/2 for 0, and 1 for positive x.

        EXAMPLES::

            sage: x = var('x')
            sage: SR(1.5).step()
            1
            sage: SR(0).step()
            1/2
            sage: SR(-1/2).step()
            0
            sage: SR(float(-1)).step()
            0
        """
        return new_Expression_from_GEx(self._parent, g_step(self._gobj))

    def csgn(self):
        """
        Return the sign of self, which is -1 if self < 0, 0 if self ==
        0, and 1 if self > 0, or unevaluated when self is a nonconstant
        symbolic expression.

        It can be somewhat arbitrary when self is not real.

        EXAMPLES:
            sage: x = var('x')
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
            sage: SR(I).csgn()
            1
        """
        return new_Expression_from_GEx(self._parent, g_csgn(self._gobj))

    def conjugate(self):
        """
        Return the complex conjugate of this symbolic expression.

        EXAMPLES::

            sage: a = 1 + 2*I
            sage: a.conjugate()
            -2*I + 1
            sage: a = sqrt(2) + 3^(1/3)*I; a
            sqrt(2) + I*3^(1/3)
            sage: a.conjugate()
            sqrt(2) - I*3^(1/3)

            sage: SR(CDF.0).conjugate()
            -1.0*I
            sage: x.conjugate()
            conjugate(x)
            sage: SR(RDF(1.5)).conjugate()
            1.5
            sage: SR(float(1.5)).conjugate()
            1.50000000000000
            sage: SR(I).conjugate()
            -I
            sage: ( 1+I  + (2-3*I)*x).conjugate()
            (3*I + 2)*conjugate(x) - I + 1
        """
        return new_Expression_from_GEx(self._parent, self._gobj.conjugate())

    def norm(self):
        r"""
        The complex norm of this symbolic expression, i.e.,
        the expression times its complex conjugate.

        EXAMPLES::

            sage: a = 1 + 2*I
            sage: a.norm()
            5
            sage: a = sqrt(2) + 3^(1/3)*I; a
            sqrt(2) + I*3^(1/3)
            sage: a.norm()
            3^(2/3) + 2
            sage: CDF(a).norm()
            4.08008382305
            sage: CDF(a.norm())
            4.08008382305
        """
        return (self*self.conjugate()).expand()

    def real_part(self):
        """
        Return the real part of this symbolic expression.

        EXAMPLES::

            sage: x = var('x')
            sage: x.real_part()
            real_part(x)
            sage: SR(2+3*I).real_part()
            2
            sage: SR(CDF(2,3)).real_part()
            2.0
            sage: SR(CC(2,3)).real_part()
            2.00000000000000

            sage: f = log(x)
            sage: f.real_part()
            log(abs(x))
        """
        return new_Expression_from_GEx(self._parent, self._gobj.real_part())

    real = real_part

    def imag_part(self):
        r"""
        Return the imaginary part of this symbolic expression.

        EXAMPLES::

            sage: sqrt(-2).imag_part()
            sqrt(2)

        We simplify `\ln(\exp(z))` to `z` for
        `-\pi<{\rm Im}(z)<=\pi`::

            sage: z = var('z')
            sage: f = log(exp(z))
            sage: assume(-pi < imag(z))
            sage: assume(imag(z) <= pi)
            sage: f
            log(e^z)
            sage: f.simplify()
            z
            sage: forget()

        A more symbolic example::

            sage: var('a, b')
            (a, b)
            sage: f = log(a + b*I)
            sage: f.imag_part()
            arctan2(real_part(b) + imag_part(a), real_part(a) - imag_part(b))

        TESTS::

            sage: x = var('x')
            sage: x.imag_part()
            imag_part(x)
            sage: SR(2+3*I).imag_part()
            3
            sage: SR(CC(2,3)).imag_part()
            3.00000000000000
            sage: SR(CDF(2,3)).imag_part()
            3.0
        """
        return new_Expression_from_GEx(self._parent, self._gobj.imag_part())

    imag = imag_part

    def sqrt(self):
        """
        EXAMPLES:
            sage: var('x, y')
            (x, y)
            sage: SR(2).sqrt()
            sqrt(2)
            sage: (x^2+y^2).sqrt()
            sqrt(x^2 + y^2)
            sage: (x^2).sqrt()
            sqrt(x^2)
        """
        return new_Expression_from_GEx(self._parent, g_sqrt(self._gobj))

    def sin(self):
        """
        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: sin(x^2 + y^2)
            sin(x^2 + y^2)
            sage: sin(sage.symbolic.constants.pi)
            0
            sage: sin(SR(1))
            sin(1)
            sage: sin(SR(RealField(150)(1)))
            0.84147098480789650665250232163029899962256306

        TESTS::

            sage: SR(oo).sin()
            Traceback (most recent call last):
            ...
            RuntimeError: sin_eval(): sin(infinity) encountered
            sage: SR(-oo).sin()
            Traceback (most recent call last):
            ...
            RuntimeError: sin_eval(): sin(infinity) encountered
            sage: SR(unsigned_infinity).sin()
            Traceback (most recent call last):
            ...
            RuntimeError: sin_eval(): sin(infinity) encountered
        """
        return new_Expression_from_GEx(self._parent, g_sin(self._gobj))

    def cos(self):
        """
        Return the cosine of self.

        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: cos(x^2 + y^2)
            cos(x^2 + y^2)
            sage: cos(sage.symbolic.constants.pi)
            -1
            sage: cos(SR(1))
            cos(1)
            sage: cos(SR(RealField(150)(1)))
            0.54030230586813971740093660744297660373231042


        In order to get a numeric approximation use .n()::

            sage: SR(RR(1)).cos().n()
            0.540302305868140
            sage: SR(float(1)).cos().n()
            0.540302305868140

        TESTS::

            sage: SR(oo).cos()
            Traceback (most recent call last):
            ...
            RuntimeError: cos_eval(): cos(infinity) encountered
            sage: SR(-oo).cos()
            Traceback (most recent call last):
            ...
            RuntimeError: cos_eval(): cos(infinity) encountered
            sage: SR(unsigned_infinity).cos()
            Traceback (most recent call last):
            ...
            RuntimeError: cos_eval(): cos(infinity) encountered
        """
        return new_Expression_from_GEx(self._parent, g_cos(self._gobj))

    def tan(self):
        """
        EXAMPLES::

            sage: var('x, y')
            (x, y)
            sage: tan(x^2 + y^2)
            tan(x^2 + y^2)
            sage: tan(sage.symbolic.constants.pi/2)
            Infinity
            sage: tan(SR(1))
            tan(1)
            sage: tan(SR(RealField(150)(1)))
            1.5574077246549022305069748074583601730872508

        TESTS::

            sage: SR(oo).tan()
            Traceback (most recent call last):
            ...
            RuntimeError: tan_eval(): tan(infinity) encountered
            sage: SR(-oo).tan()
            Traceback (most recent call last):
            ...
            RuntimeError: tan_eval(): tan(infinity) encountered
            sage: SR(unsigned_infinity).tan()
            Traceback (most recent call last):
            ...
            RuntimeError: tan_eval(): tan(infinity) encountered
        """
        return new_Expression_from_GEx(self._parent, g_tan(self._gobj))

    def arcsin(self):
        """
        Return the arcsin of x, i.e., the number y between -pi and pi
        such that sin(y) == x.

        EXAMPLES::

            sage: x.arcsin()
            arcsin(x)
            sage: SR(0.5).arcsin()
            arcsin(0.500000000000000)

        Use .n() to get a numerical approximation::

            sage: SR(0.5).arcsin().n()
            0.523598775598299

            sage: SR(0.999).arcsin()
            arcsin(0.999000000000000)
            sage: SR(-0.999).arcsin()
            -arcsin(0.999000000000000)
            sage: SR(0.999).arcsin().n()
            1.52607123962616

        TESTS::

            sage: SR(oo).arcsin()
            Traceback (most recent call last):
            ...
            RuntimeError: arcsin_eval(): arcsin(infinity) encountered
            sage: SR(-oo).arcsin()
            Traceback (most recent call last):
            ...
            RuntimeError: arcsin_eval(): arcsin(infinity) encountered
            sage: SR(unsigned_infinity).arcsin()
            Infinity
        """
        return new_Expression_from_GEx(self._parent, g_asin(self._gobj))

    def arccos(self):
        """
        Return the arc cosine of self.

        EXAMPLES::

            sage: x.arccos()
            arccos(x)
            sage: SR(1).arccos()
            0
            sage: SR(1/2).arccos()
            1/3*pi
            sage: SR(0.4).arccos()
            arccos(0.400000000000000)

        Use .n() to get a numerical approximation::

            sage: SR(0.4).arccos().n()
            1.15927948072741
            sage: plot(lambda x: SR(x).arccos(), -1,1)

        TESTS::

            sage: SR(oo).arccos()
            Traceback (most recent call last):
            ...
            RuntimeError: arccos_eval(): arccos(infinity) encountered
            sage: SR(-oo).arccos()
            Traceback (most recent call last):
            ...
            RuntimeError: arccos_eval(): arccos(infinity) encountered
            sage: SR(unsigned_infinity).arccos()
            Infinity
        """
        return new_Expression_from_GEx(self._parent, g_acos(self._gobj))

    def arctan(self):
        """
        Return the arc tangent of self.

        EXAMPLES::

            sage: x = var('x')
            sage: x.arctan()
            arctan(x)
            sage: SR(1).arctan()
            1/4*pi
            sage: SR(1/2).arctan()
            arctan(1/2)
            sage: SR(0.5).arctan()
            arctan(0.500000000000000)

        Use .n() to get a numerical approximation::

            sage: SR(0.5).arctan().n()
            0.463647609000806
            sage: plot(lambda x: SR(x).arctan(), -20,20)

        TESTS::

            sage: SR(oo).arctan()
            1/2*pi
            sage: SR(-oo).arctan()
            -1/2*pi
            sage: SR(unsigned_infinity).arctan()
            Traceback (most recent call last):
            ...
            RuntimeError: arctan_eval(): arctan(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent, g_atan(self._gobj))

    def arctan2(self, x):
        """
        Return the inverse of the 2-variable tan function on self and x.

        EXAMPLES::

            sage: var('x,y')
            (x, y)
            sage: x.arctan2(y)
            arctan2(x, y)
            sage: SR(1/2).arctan2(1/2)
            1/4*pi
            sage: maxima.eval('atan2(1/2,1/2)')
            '%pi/4'

            sage: SR(-0.7).arctan2(SR(-0.6))
            -pi + arctan(1.16666666666667)

        Use .n() to get a numerical approximation::

            sage: SR(-0.7).arctan2(SR(-0.6)).n()
            -2.27942259892257

        TESTS:

        We compare a bunch of different evaluation points between
        Sage and Maxima::

            sage: float(SR(0.7).arctan2(0.6))
            0.8621700546672264
            sage: maxima('atan2(0.7,0.6)')
            .862170054667226...
            sage: float(SR(0.7).arctan2(-0.6))
            2.2794225989225669
            sage: maxima('atan2(0.7,-0.6)')
            2.279422598922567
            sage: float(SR(-0.7).arctan2(0.6))
            -0.8621700546672264
            sage: maxima('atan2(-0.7,0.6)')
            -.862170054667226...
            sage: float(SR(-0.7).arctan2(-0.6))
            -2.2794225989225669
            sage: maxima('atan2(-0.7,-0.6)')
            -2.279422598922567
            sage: float(SR(0).arctan2(-0.6))
            3.1415926535897931
            sage: maxima('atan2(0,-0.6)')
            3.141592653589793
            sage: float(SR(0).arctan2(0.6))
            0.0
            sage: maxima('atan2(0,0.6)')
            0.0
            sage: SR(0).arctan2(0)
            0

            sage: SR(I).arctan2(1)
            arctan2(I, 1)
            sage: SR(CDF(0,1)).arctan2(1)
            arctan2(1.0*I, 1)
            sage: SR(1).arctan2(CDF(0,1))
            arctan2(1, 1.0*I)

            sage: SR(oo).arctan2(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: arctan2_eval(): arctan2(infinity, infinity) encountered
            sage: SR(oo).arctan2(0)
            0
            sage: SR(-oo).arctan2(0)
            pi
            sage: SR(-oo).arctan2(-2)
            -pi
            sage: SR(unsigned_infinity).arctan2(2)
            Traceback (most recent call last):
            ...
            RuntimeError: arctan2_eval(): arctan2(unsigned_infinity, x) encountered
            sage: SR(2).arctan2(oo)
            1/2*pi
            sage: SR(2).arctan2(-oo)
            -1/2*pi
            sage: SR(2).arctan2(SR(unsigned_infinity))
            Traceback (most recent call last):
            ...
            RuntimeError: arctan2_eval(): arctan2(x, unsigned_infinity) encountered
        """
        cdef Expression nexp = self.coerce_in(x)
        return new_Expression_from_GEx(self._parent, g_atan2(self._gobj, nexp._gobj))

    def sinh(self):
        r"""
        Return sinh of self.

        We have $\sinh(x) = (e^{x} - e^{-x})/2$.

        EXAMPLES::

            sage: x.sinh()
            sinh(x)
            sage: SR(1).sinh()
            sinh(1)
            sage: SR(0).sinh()
            0
            sage: SR(1.0).sinh()
            sinh(1.00000000000000)

        Use .n() to get a numerical approximation::

            sage: SR(1.0).sinh().n()
            1.17520119364380
            sage: maxima('sinh(1.0)')
            1.175201193643801

            sinh(1.0000000000000000000000000)
            sage: SR(1).sinh().n(90)
            1.1752011936438014568823819
            sage: SR(RIF(1)).sinh()
            sinh(1)
            sage: SR(RIF(1)).sinh().n()
            1.175201193643802?

        TESTS::

            sage: SR(oo).sinh()
            +Infinity
            sage: SR(-oo).sinh()
            -Infinity
            sage: SR(unsigned_infinity).sinh()
            Traceback (most recent call last):
            ...
            RuntimeError: sinh_eval(): sinh(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent, g_sinh(self._gobj))

    def cosh(self):
        """
        Return cosh of self.

        We have $\sinh(x) = (e^{x} + e^{-x})/2$.

        EXAMPLES::

            sage: x.cosh()
            cosh(x)
            sage: SR(1).cosh()
            cosh(1)
            sage: SR(0).cosh()
            1
            sage: SR(1.0).cosh()
            cosh(1.00000000000000)

        Use .n() to get a numerical approximation::

            sage: SR(1.0).cosh().n()
            1.54308063481524
            sage: maxima('cosh(1.0)')
            1.543080634815244
            sage: SR(1.0000000000000000000000000).cosh()
            cosh(1.000000000000000000000000)
            sage: SR(1).cosh().n(90)
            1.5430806348152437784779056
            sage: SR(RIF(1)).cosh()
            cosh(1)
            sage: SR(RIF(1)).cosh().n()
            1.543080634815244?

        TESTS::

            sage: SR(oo).cosh()
            +Infinity
            sage: SR(-oo).cosh()
            +Infinity
            sage: SR(unsigned_infinity).cosh()
            Traceback (most recent call last):
            ...
            RuntimeError: cosh_eval(): cosh(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent, g_cosh(self._gobj))

    def tanh(self):
        """
        Return tanh of self.

        We have $\tanh(x) = \sinh(x) / \cosh(x)$.

        EXAMPLES::

            sage: x.tanh()
            tanh(x)
            sage: SR(1).tanh()
            tanh(1)
            sage: SR(0).tanh()
            0
            sage: SR(1.0).tanh()
            tanh(1.00000000000000)

        Use .n() to get a numerical approximation::

            sage: SR(1.0).tanh().n()
            0.761594155955765
            sage: maxima('tanh(1.0)')
            .7615941559557649
            sage: plot(lambda x: SR(x).tanh(), -1, 1)

        TESTS::

            sage: SR(oo).tanh()
            1
            sage: SR(-oo).tanh()
            -1
            sage: SR(unsigned_infinity).tanh()
            Traceback (most recent call last):
            ...
            RuntimeError: tanh_eval(): tanh(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent, g_tanh(self._gobj))

    def arcsinh(self):
        """
        Return the inverse hyperbolic sine of self.

        EXAMPLES::

            sage: x.arcsinh()
            arcsinh(x)
            sage: SR(0).arcsinh()
            0
            sage: SR(1).arcsinh()
            arcsinh(1)
            sage: SR(1.0).arcsinh()
            arcsinh(1.00000000000000)

        Use .n() to get a numerical approximation::

            sage: SR(1.0).arcsinh().n()
            0.881373587019543
            sage: maxima('asinh(1.0)')
            0.881373587019543

        Sage automatically applies certain identities::
            sage: SR(3/2).arcsinh().cosh()
            1/2*sqrt(13)

        TESTS::

            sage: SR(oo).arcsinh()
            +Infinity
            sage: SR(-oo).arcsinh()
            -Infinity
            sage: SR(unsigned_infinity).arcsinh()
            Infinity
        """
        return new_Expression_from_GEx(self._parent, g_asinh(self._gobj))

    def arccosh(self):
        """
        Return the inverse hyperbolic cosine of self.

        EXAMPLES::

            sage: x.arccosh()
            arccosh(x)
            sage: SR(0).arccosh()
            1/2*I*pi
            sage: SR(1/2).arccosh()
            arccosh(1/2)
            sage: SR(CDF(1/2)).arccosh()
            arccosh(0.5)

        Use .n() to get a numerical approximation::

            sage: SR(CDF(1/2)).arccosh().n()
            1.0471975512*I
            sage: maxima('acosh(0.5)')
            1.047197551196598*%i

        TESTS::

            sage: SR(oo).arccosh()
            +Infinity
            sage: SR(-oo).arccosh()
            +Infinity
            sage: SR(unsigned_infinity).arccosh()
            +Infinity
        """
        return new_Expression_from_GEx(self._parent, g_acosh(self._gobj))

    def arctanh(self):
        """
        Return the inverse hyperbolic tangent of self.

        EXAMPLES::

            sage: x.arctanh()
            arctanh(x)
            sage: SR(0).arctanh()
            0
            sage: SR(1/2).arctanh()
            arctanh(1/2)
            sage: SR(0.5).arctanh()
            arctanh(0.500000000000000)

        Use .n() to get a numerical approximation::

            sage: SR(0.5).arctanh().n()
            0.549306144334055
            sage: SR(0.5).arctanh().tanh()
            0.500000000000000
            sage: maxima('atanh(0.5)')
            .5493061443340...

        TESTS::

            sage: SR(1).arctanh()
            +Infinity
            sage: SR(-1).arctanh()
            -Infinity

            sage: SR(oo).arctanh()
            -1/2*I*pi
            sage: SR(-oo).arctanh()
            1/2*I*pi
            sage: SR(unsigned_infinity).arctanh()
            Traceback (most recent call last):
            ...
            RuntimeError: arctanh_eval(): arctanh(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent, g_atanh(self._gobj))

    def exp(self):
        """
        Return exponential function of self, i.e., e to the
        power of self.

        EXAMPLES::

            sage: x.exp()
            e^x
            sage: SR(0).exp()
            1
            sage: SR(1/2).exp()
            e^(1/2)
            sage: SR(0.5).exp()
            e^0.500000000000000
            sage: (pi*I).exp()
            -1

        Use .n() to get a numerical approximation::

            sage: SR(0.5).exp().n()
            1.64872127070013
            sage: math.exp(0.5)
            1.6487212707001282

            sage: SR(0.5).exp().log()
            0.500000000000000

        TESTS:

        Test if #6377 is fixed::

            sage: SR(oo).exp()
            +Infinity
            sage: SR(-oo).exp()
            0
            sage: SR(unsigned_infinity).exp()
            Traceback (most recent call last):
            ...
            RuntimeError: exp_eval(): exp^(unsigned_infinity) encountered
        """
        return new_Expression_from_GEx(self._parent, g_exp(self._gobj))

    def log(self, b=None):
        """
        Return the logarithm of self.

        EXAMPLES::

            sage: x, y = var('x, y')
            sage: x.log()
            log(x)
            sage: (x^y + y^x).log()
            log(x^y + y^x)
            sage: SR(0).log()
            -Infinity
            sage: SR(-1).log()
            I*pi
            sage: SR(1).log()
            0
            sage: SR(1/2).log()
            log(1/2)
            sage: SR(0.5).log()
            log(0.500000000000000)

        Use .n() to get a numerical approximation::

            sage: SR(0.5).log().n()
            -0.693147180559945
            sage: SR(0.5).log().exp()
            0.500000000000000
            sage: math.log(0.5)
            -0.69314718055994529
            sage: plot(lambda x: SR(x).log(), 0.1,10)

        TESTS::

            sage: SR(oo).log()
            +Infinity
            sage: SR(-oo).log()
            +Infinity
            sage: SR(unsigned_infinity).log()
            +Infinity
        """
        res = new_Expression_from_GEx(self._parent, g_log(self._gobj))
        if b is None:
            return res
        else:
            return res/self._parent(b).log()

    def zeta(self):
        """
        EXAMPLES::

            sage: x, y = var('x, y')
            sage: (x/y).zeta()
            zeta(x/y)
            sage: SR(2).zeta()
            1/6*pi^2
            sage: SR(3).zeta()
            zeta(3)
            sage: SR(CDF(0,1)).zeta()
            zeta(1.0*I)
            sage: SR(CDF(0,1)).zeta().n()
            0.00330022368532 - 0.418155449141*I
            sage: CDF(0,1).zeta()
            0.00330022368532 - 0.418155449141*I
            sage: plot(lambda x: SR(x).zeta(), -10,10).show(ymin=-3,ymax=3)

        TESTS::

            sage: t = SR(1).zeta(); t
            zeta(1)
            sage: t.n()
            +infinity
        """
        _sig_on
        cdef GEx x = g_zeta(self._gobj)
        _sig_off
        return new_Expression_from_GEx(self._parent, x)

    def factorial(self):
        """
        Return the factorial of self.

        OUTPUT:
            symbolic expression

        EXAMPLES:
            sage: var('x, y')
            (x, y)
            sage: SR(5).factorial()
            120
            sage: x.factorial()
            factorial(x)
            sage: (x^2+y^3).factorial()
            factorial(x^2 + y^3)
        """
        _sig_on
        cdef GEx x = g_factorial(self._gobj)
        _sig_off
        return new_Expression_from_GEx(self._parent, x)

    def binomial(self, k):
        """
        Return binomial coefficient "self choose k".

        OUTPUT:
            symbolic expression

        EXAMPLES:
            sage: var('x, y')
            (x, y)
            sage: SR(5).binomial(SR(3))
            10
            sage: x.binomial(SR(3))
            1/6*x^3 - 1/2*x^2 + 1/3*x
            sage: x.binomial(y)
            binomial(x,y)
        """
        cdef Expression nexp = self.coerce_in(k)
        _sig_on
        cdef GEx x = g_binomial(self._gobj, nexp._gobj)
        _sig_off
        return new_Expression_from_GEx(self._parent, x)

    def Order(self):
        """
        Order, as in big oh notation.

        OUTPUT:
            symbolic expression

        EXAMPLES:
            sage: n = var('n')
            sage: (17*n^3).Order()
            Order(n^3)
        """
        return new_Expression_from_GEx(self._parent, g_Order(self._gobj))

    def gamma(self):
        """
        Return the Gamma function evaluated at self.

        EXAMPLES:
            sage: x = var('x')
            sage: x.gamma()
            gamma(x)
            sage: SR(2).gamma()
            1
            sage: SR(10).gamma()
            362880
            sage: SR(10.0r).gamma()
            gamma(10.0000000000000)

        Use .n() to get a numerical approximation::

            sage: SR(10.0r).gamma().n()
            362880.000000000
            sage: SR(CDF(1,1)).gamma()
            gamma(1.0 + 1.0*I)
            sage: SR(CDF(1,1)).gamma().n()
            0.498015668118 - 0.154949828302*I

            sage: gp('gamma(1+I)') # 32-bit
            0.4980156681183560427136911175 - 0.1549498283018106851249551305*I

            sage: gp('gamma(1+I)') # 64-bit
            0.49801566811835604271369111746219809195 - 0.15494982830181068512495513048388660520*I

            sage: set_verbose(-1); plot(lambda x: SR(x).gamma(), -6,4).show(ymin=-3,ymax=3)
        """
        _sig_on
        cdef GEx x = g_tgamma(self._gobj)
        _sig_off
        return new_Expression_from_GEx(self._parent, x)

    def lgamma(self):
        """
        Return the log-gamma function evaluated at self.
        This is the logarithm of gamma of self, where
        gamma is a complex function such that gamma(n)
        equals factorial(n-1).

        EXAMPLES:
            sage: x = var('x')
            sage: x.lgamma()
            lgamma(x)
            sage: SR(2).lgamma()
            0
            sage: SR(5).lgamma()
            log(24)
            sage: SR(5-1).factorial().log()
            log(24)
            sage: set_verbose(-1); plot(lambda x: SR(x).lgamma(), -7,8, plot_points=1000).show()
            sage: math.exp(0.5)
            1.6487212707001282
            sage: plot(lambda x: (SR(x).exp() - SR(-x).exp())/2 - SR(x).sinh(), -1, 1)
        """
        _sig_on
        cdef GEx x = g_lgamma(self._gobj)
        _sig_off
        return new_Expression_from_GEx(self._parent, x)

    def default_variable(self):
        """
        Return the default variable, which is by definition the first
        variable in self, or `x` is there are no variables in self.
        The result is cached.

        EXAMPLES::

            sage: sqrt(2).default_variable()
            x
            sage: x, theta, a = var('x, theta, a')
            sage: f = x^2 + theta^3 - a^x
            sage: f.default_variable()
            a

        Note that this is the first *variable*, not the first *argument*::

            sage: f(theta, a, x) = a + theta^3
            sage: f.default_variable()
            a
            sage: f.variables()
            (a, theta)
            sage: f.arguments()
            (theta, a, x)
        """
        v = self.variables()
        if len(v) == 0:
            return self.parent().var('x')
        else:
            return v[0]

    def combine(self):
        r"""
        Returns a simplified version of this symbolic expression
        by combining all terms with the same denominator into a single
        term.

        EXAMPLES::

            sage: var('x, y, a, b, c')
            (x, y, a, b, c)
            sage: f = x*(x-1)/(x^2 - 7) + y^2/(x^2-7) + 1/(x+1) + b/a + c/a; f
            (x - 1)*x/(x^2 - 7) + y^2/(x^2 - 7) + b/a + c/a + 1/(x + 1)
            sage: f.combine()
            ((x - 1)*x + y^2)/(x^2 - 7) + (b + c)/a + 1/(x + 1)
        """
        return self.parent()(self._maxima_().combine())

    def numerator(self):
        """
        Returns the numerator of this symbolic expression.  If the
        expression is not a quotient, then this will return the
        expression itself.

        EXAMPLES::

            sage: a, x, y = var('a,x,y')
            sage: f = x*(x-a)/((x^2 - y)*(x-a)); f
            x/(x^2 - y)
            sage: f.numerator()
            x
            sage: f.denominator()
            x^2 - y

            sage: y = var('y')
            sage: g = x + y/(x + 2); g
            x + y/(x + 2)
            sage: g.numerator()
            x + y/(x + 2)
            sage: g.denominator()
            1

        """
        return self.parent()(self._maxima_().num())

    def denominator(self):
        """
        Returns the denominator of this symbolic expression.  If the
        expression is not a quotient, then this will just return 1.

        EXAMPLES::

            sage: x, y, z, theta = var('x, y, z, theta')
            sage: f = (sqrt(x) + sqrt(y) + sqrt(z))/(x^10 - y^10 - sqrt(theta))
            sage: f.denominator()
            sqrt(theta) - x^10 + y^10

            sage: y = var('y')
            sage: g = x + y/(x + 2); g
            x + y/(x + 2)
            sage: g.numerator()
            x + y/(x + 2)
            sage: g.denominator()
            1
        """
        return self.parent()(self._maxima_().denom())

    def partial_fraction(self, var=None):
        r"""
        Return the partial fraction expansion of ``self`` with
        respect to the given variable.

        INPUT:


        -  ``var`` - variable name or string (default: first
           variable)


        OUTPUT: Symbolic expression

        EXAMPLES::

            sage: f = x^2/(x+1)^3
            sage: f.partial_fraction()
            1/(x + 1) - 2/(x + 1)^2 + 1/(x + 1)^3
            sage: f.partial_fraction()
            1/(x + 1) - 2/(x + 1)^2 + 1/(x + 1)^3

        Notice that the first variable in the expression is used by
        default::

            sage: y = var('y')
            sage: f = y^2/(y+1)^3
            sage: f.partial_fraction()
            1/(y + 1) - 2/(y + 1)^2 + 1/(y + 1)^3

            sage: f = y^2/(y+1)^3 + x/(x-1)^3
            sage: f.partial_fraction()
            y^2/(y^3 + 3*y^2 + 3*y + 1) + 1/(x - 1)^2 + 1/(x - 1)^3

        You can explicitly specify which variable is used::

            sage: f.partial_fraction(y)
            x/(x^3 - 3*x^2 + 3*x - 1) + 1/(y + 1) - 2/(y + 1)^2 + 1/(y + 1)^3
        """
        if var is None:
            var = self.default_variable()
        return self.parent()(self._maxima_().partfrac(var))

    def simplify(self):
        """
        Returns a simplified version of this symbolic expression.

        .. note::

           Currently, this does just sends the expression to Maxima
           and converts it back to Sage.

        .. seealso::

           :meth:`simplify_full`, :meth:`simplify_trig`,
           :meth:`simplify_rational`, :meth:`simplify_radical`

        EXAMPLES::

            sage: a = var('a'); f = x*sin(2)/(x^a); f
            x*sin(2)/x^a
            sage: f.simplify()
            x^(-a + 1)*sin(2)
        """
        return self._parent(self._maxima_())

    def simplify_full(self):
        """
        Applies simplify_trig, simplify_rational, and simplify_radical
        to self (in that order).

        ALIAS: simplify_full and full_simplify are the same.

        EXAMPLES::

            sage: a = log(8)/log(2)
            sage: a.simplify_full()
            3

        ::

            sage: f = sin(x)^2 + cos(x)^2
            sage: f.simplify_full()
            1

        ::

            sage: f = sin(x/(x^2 + x))
            sage: f.simplify_full()
            sin(1/(x + 1))
        """
        x = self
        x = x.simplify_trig()
        x = x.simplify_rational()
        x = x.simplify_radical()
        return x

    full_simplify = simplify_full

    def simplify_trig(self):
        r"""
        First expands using trig_expand, then employs the identities
        `\sin(x)^2 + \cos(x)^2 = 1` and
        `\cosh(x)^2 - \sin(x)^2 = 1` to simplify expressions
        containing tan, sec, etc., to sin, cos, sinh, cosh.

        ALIAS: trig_simplify and simplify_trig are the same

        EXAMPLES::

            sage: f = sin(x)^2 + cos(x)^2; f
            sin(x)^2 + cos(x)^2
            sage: f.simplify()
            sin(x)^2 + cos(x)^2
            sage: f.simplify_trig()
            1
        """
        # much better to expand first, since it often doesn't work
        # right otherwise!
        return self.parent()(self._maxima_().trigexpand().trigsimp())

    trig_simplify = simplify_trig

    def simplify_rational(self):
        """
        Simplify by expanding repeatedly rational expressions.

        ALIAS: rational_simplify and simplify_rational are the same

        EXAMPLES::

            sage: f = sin(x/(x^2 + x))
            sage: f
            sin(x/(x^2 + x))
            sage: f.simplify_rational()
            sin(1/(x + 1))

        ::

            sage: f = ((x - 1)^(3/2) - (x + 1)*sqrt(x - 1))/sqrt((x - 1)*(x + 1)); f
            ((x - 1)^(3/2) - sqrt(x - 1)*(x + 1))/sqrt((x - 1)*(x + 1))
            sage: f.simplify_rational()
            -2*sqrt(x - 1)/sqrt(x^2 - 1)
        """
        return self.parent()(self._maxima_().fullratsimp())

    rational_simplify = simplify_rational

    # TODO: come up with a way to intelligently wrap Maxima's way of
    # fine-tuning all simplificationsrational

    def simplify_radical(self):
        r"""
        Simplifies this symbolic expression, which can contain logs,
        exponentials, and radicals, by converting it into a form which is
        canonical over a large class of expressions and a given ordering of
        variables

        DETAILS: This uses the Maxima radcan() command. From the Maxima
        documentation: "All functionally equivalent forms are mapped into a
        unique form. For a somewhat larger class of expressions, produces a
        regular form. Two equivalent expressions in this class do not
        necessarily have the same appearance, but their difference can be
        simplified by radcan to zero. For some expressions radcan is quite
        time consuming. This is the cost of exploring certain relationships
        among the components of the expression for simplifications based on
        factoring and partial fraction expansions of exponents."

        ALIAS: radical_simplify, simplify_radical, simplify_log,
        log_simplify, exp_simplify, simplify_exp are all the same

        EXAMPLES::

            sage: var('x,y,a')
            (x, y, a)

        ::

            sage: f = log(x*y)
            sage: f.simplify_radical()
            log(x) + log(y)

        ::

            sage: f = (log(x+x^2)-log(x))^a/log(1+x)^(a/2)
            sage: f.simplify_radical()
            log(x + 1)^(1/2*a)

        ::

            sage: f = (e^x-1)/(1+e^(x/2))
            sage: f.simplify_exp()
            e^(1/2*x) - 1
        """
        from sage.calculus.calculus import maxima
        maxima.eval('domain: real$')
        res = self.parent()(self._maxima_().radcan())
        maxima.eval('domain: complex$')
        return res

    radical_simplify = simplify_log = log_simplify = simplify_radical
    simplify_exp = exp_simplify = simplify_radical


    def factor(self, dontfactor=[]):
        """
        Factors self, containing any number of variables or functions, into
        factors irreducible over the integers.

        INPUT:


        -  ``self`` - a symbolic expression

        -  ``dontfactor`` - list (default: []), a list of
           variables with respect to which factoring is not to occur.
           Factoring also will not take place with respect to any variables
           which are less important (using the variable ordering assumed for
           CRE form) than those on the 'dontfactor' list.


        EXAMPLES::

            sage: x,y,z = var('x, y, z')
            sage: (x^3-y^3).factor()
            (x - y)*(x^2 + x*y + y^2)
            sage: factor(-8*y - 4*x + z^2*(2*y + x))
            (z - 2)*(z + 2)*(x + 2*y)
            sage: f = -1 - 2*x - x^2 + y^2 + 2*x*y^2 + x^2*y^2
            sage: F = factor(f/(36*(1 + 2*y + y^2)), dontfactor=[x]); F
            1/36*(y - 1)*(x^2 + 2*x + 1)/(y + 1)

        If you are factoring a polynomial with rational coefficients (and
        dontfactor is empty) the factorization is done using Singular
        instead of Maxima, so the following is very fast instead of
        dreadfully slow::

            sage: var('x,y')
            (x, y)
            sage: (x^99 + y^99).factor()
            (x + y)*(x^2 - x*y + y^2)*(x^6 - x^3*y^3 + y^6)*...
        """
        from sage.calculus.calculus import symbolic_expression_from_maxima_string, symbolic_expression_from_string
        if len(dontfactor) > 0:
            m = self._maxima_()
            name = m.name()
            cmd = 'block([dontfactor:%s],factor(%s))'%(dontfactor, name)
            return symbolic_expression_from_maxima_string(cmd)
        else:
            try:
                from sage.rings.all import QQ
                f = self.polynomial(QQ)
                w = repr(f.factor())
                return symbolic_expression_from_string(w)
            except TypeError:
                pass
            return self.parent()(self._maxima_().factor())

    def factor_list(self, dontfactor=[]):
        """
        Returns a list of the factors of self, as computed by the
        factor command.

        INPUT:

        -  ``self`` - a symbolic expression

        -  ``dontfactor`` - see docs for :meth:`factor`

        .. note::

           If you already have a factored expression and just want to
           get at the individual factors, use :meth:`_factor_list`
           instead.

        EXAMPLES::

            sage: var('x, y, z')
            (x, y, z)
            sage: f = x^3-y^3
            sage: f.factor()
            (x - y)*(x^2 + x*y + y^2)

        Notice that the -1 factor is separated out::

            sage: f.factor_list()
            [(x - y, 1), (x^2 + x*y + y^2, 1)]

        We factor a fairly straightforward expression::

            sage: factor(-8*y - 4*x + z^2*(2*y + x)).factor_list()
            [(z - 2, 1), (z + 2, 1), (x + 2*y, 1)]

        A more complicated example::

            sage: var('x, u, v')
            (x, u, v)
            sage: f = expand((2*u*v^2-v^2-4*u^3)^2 * (-u)^3 * (x-sin(x))^3)
            sage: f.factor()
            -(x - sin(x))^3*(4*u^3 - 2*u*v^2 + v^2)^2*u^3
            sage: g = f.factor_list(); g
            [(x - sin(x), 3), (4*u^3 - 2*u*v^2 + v^2, 2), (u, 3), (-1, 1)]

        This function also works for quotients::

            sage: f = -1 - 2*x - x^2 + y^2 + 2*x*y^2 + x^2*y^2
            sage: g = f/(36*(1 + 2*y + y^2)); g
            1/36*(x^2*y^2 + 2*x*y^2 - x^2 + y^2 - 2*x - 1)/(y^2 + 2*y + 1)
            sage: g.factor(dontfactor=[x])
            1/36*(y - 1)*(x^2 + 2*x + 1)/(y + 1)
            sage: g.factor_list(dontfactor=[x])
            [(y - 1, 1), (y + 1, -1), (x^2 + 2*x + 1, 1), (1/36, 1)]

        This example also illustrates that the exponents do not have to be
        integers::

            sage: f = x^(2*sin(x)) * (x-1)^(sqrt(2)*x); f
            (x - 1)^(sqrt(2)*x)*x^(2*sin(x))
            sage: f.factor_list()
            [(x - 1, sqrt(2)*x), (x, 2*sin(x))]
        """
        return self.factor(dontfactor=dontfactor)._factor_list()

    def _factor_list(self):
        r"""
        Turn an expression already in factored form into a list of (prime,
        power) pairs.

        This is used, e.g., internally by the :meth:`factor_list`
        command.

        EXAMPLES::

            sage: g = factor(x^3 - 1); g
            (x - 1)*(x^2 + x + 1)
            sage: v = g._factor_list(); v
            [(x - 1, 1), (x^2 + x + 1, 1)]
            sage: type(v)
            <type 'list'>
        """
        op = self.operator()
        if op is operator.mul:
            return sum([f._factor_list() for f in self.operands()], [])
        elif op is operator.pow:
            return [tuple(self.operands())]
        else:
            return [(self, 1)]

    ###################################################################
    # solve
    ###################################################################
    def roots(self, x=None, explicit_solutions=True, multiplicities=True, ring=None):
        r"""
        Returns roots of ``self`` that can be found exactly,
        possibly with multiplicities.  Not all roots are guaranteed to
        be found.

        .. warning::

           This is *not* a numerical solver - use ``find_root`` to
           solve for self == 0 numerically on an interval.

        INPUT:

        - ``x`` - variable to view the function in terms of
          (use default variable if not given)

        - ``explicit_solutions`` - bool (default True); require that
          roots be explicit rather than implicit

        - ``multiplicities`` - bool (default True); when True, return
          multiplicities

        - ``ring`` - a ring (default None): if not None, convert
          self to a polynomial over ring and find roots over ring

        OUTPUT:

        list of pairs (root, multiplicity) or list of roots

        If there are infinitely many roots, e.g., a function like
        `\sin(x)`, only one is returned.

        EXAMPLES::

            sage: var('x, a')
            (x, a)

        A simple example::

            sage: ((x^2-1)^2).roots()
            [(-1, 2), (1, 2)]
            sage: ((x^2-1)^2).roots(multiplicities=False)
            [-1, 1]

        A complicated example.

        ::

            sage: f = expand((x^2 - 1)^3*(x^2 + 1)*(x-a)); f
            -a*x^8 + x^9 + 2*a*x^6 - 2*x^7 - 2*a*x^2 + 2*x^3 + a - x

        The default variable is `a`, since it is the first in
        alphabetical order::

            sage: f.roots()
            [(x, 1)]

        As a polynomial in `a`, `x` is indeed a root::

            sage: f.poly(a)
            x^9 - 2*x^7 + 2*x^3 - (x^8 - 2*x^6 + 2*x^2 - 1)*a - x
            sage: f(a=x)
            0

        The roots in terms of `x` are what we expect::

            sage: f.roots(x)
            [(a, 1), (-I, 1), (I, 1), (1, 3), (-1, 3)]

        Only one root of `\sin(x) = 0` is given::

            sage: f = sin(x)
            sage: f.roots(x)
            [(0, 1)]

        We derive the roots of a general quadratic polynomial::

            sage: var('a,b,c,x')
            (a, b, c, x)
            sage: (a*x^2 + b*x + c).roots(x)
            [(-1/2*(b + sqrt(-4*a*c + b^2))/a, 1), (-1/2*(b - sqrt(-4*a*c + b^2))/a, 1)]

        By default, all the roots are required to be explicit rather than
        implicit. To get implicit roots, pass
        ``explicit_solutions=False`` to
        ``.roots()``

        ::

            sage: var('x')
            x
            sage: f = x^(1/9) + (2^(8/9) - 2^(1/9))*(x - 1) - x^(8/9)
            sage: f.roots()
            Traceback (most recent call last):
            ...
            RuntimeError: no explicit roots found
            sage: f.roots(explicit_solutions=False)
            [((2^(8/9) - 2^(1/9) + x^(8/9) - x^(1/9))/(2^(8/9) - 2^(1/9)), 1)]

        Another example, but involving a degree 5 poly whose roots don't
        get computed explicitly::

            sage: f = x^5 + x^3 + 17*x + 1
            sage: f.roots()
            Traceback (most recent call last):
            ...
            RuntimeError: no explicit roots found
            sage: f.roots(explicit_solutions=False)
            [(x^5 + x^3 + 17*x + 1, 1)]
            sage: f.roots(explicit_solutions=False, multiplicities=False)
            [x^5 + x^3 + 17*x + 1]

        Now let's find some roots over different rings::

            sage: f.roots(ring=CC)
            [(-0.0588115223184495, 1), (-1.331099917875... - 1.52241655183732*I, 1), (-1.331099917875... + 1.52241655183732*I, 1), (1.36050567903502 - 1.51880872209965*I, 1), (1.36050567903502 + 1.51880872209965*I, 1)]
            sage: (2.5*f).roots(ring=RR)
            [(-0.058811522318449..., 1)]
            sage: f.roots(ring=CC, multiplicities=False)
            [-0.0588115223184495, -1.331099917875... - 1.52241655183732*I, -1.331099917875... + 1.52241655183732*I, 1.36050567903502 - 1.51880872209965*I, 1.36050567903502 + 1.51880872209965*I]
            sage: f.roots(ring=QQ)
            []
            sage: f.roots(ring=QQbar, multiplicities=False)
            [-0.05881152231844944?, -1.331099917875796? - 1.522416551837318?*I, -1.331099917875796? + 1.522416551837318?*I, 1.360505679035020? - 1.518808722099650?*I, 1.360505679035020? + 1.518808722099650?*I]

        Root finding over finite fields::

            sage: f.roots(ring=GF(7^2, 'a'))
            [(3, 1), (4*a + 6, 2), (3*a + 3, 2)]

        TESTS::

            sage: (sqrt(3) * f).roots(ring=QQ)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert sqrt(3) to a rational
        """
        if x is None:
            x = self.default_variable()
        if ring is not None:
            p = self.polynomial(ring)
            return p.roots(ring=ring, multiplicities=multiplicities)

        S, mul = self.solve(x, multiplicities=True, explicit_solutions=explicit_solutions)
        if len(mul) == 0 and explicit_solutions:
            raise RuntimeError, "no explicit roots found"
        else:
            rt_muls = [(S[i].rhs(), mul[i]) for i in range(len(mul))]
        if multiplicities:
            return rt_muls
        else:
            return [ rt for rt, mul in rt_muls ]

    def solve(self, x, multiplicities=False, solution_dict=False, explicit_solutions=None):
        r"""
        Analytically solve the equation ``self == 0`` for the
        variable `x`.

        .. warning::

           This is not a numerical solver - use ``find_root`` to solve
           for self == 0 numerically on an interval.

        INPUT:


        -  ``x`` - variable to solve for

        -  ``multiplicities`` - bool (default: False); if True,
           return corresponding multiplicities.

        - ``explicit_solutions`` - bool; if True, require that all
           solutions returned be explicit (rather than implicit)


        EXAMPLES::

            sage: z = var('z')
            sage: (z^5 - 1).solve(z)
            [z == e^(2/5*I*pi), z == e^(4/5*I*pi), z == e^(-4/5*I*pi), z == e^(-2/5*I*pi), z == 1]

            sage: var('Q')
            Q
            sage: solve(Q*sqrt(Q^2 + 2) - 1,Q)
            [Q == 1/sqrt(-sqrt(2) + 1), Q == 1/sqrt(sqrt(2) + 1)]
        """
        import operator
        cdef Expression ex
        if is_a_relational(self._gobj):
            if self.operator() is not operator.eq:
                raise NotImplementedError, "solving only implemented for equalities"
            ex = self
        else:
            ex = (self == 0)

        if x is None:
            v = ex.variables()
            if len(v) == 0:
                if multiplicities:
                    return [], []
                else:
                    return []
            x = v[0]

        if not isinstance(x, Expression):
            raise TypeError, "%s is not a valid variable."%x

        m = ex._maxima_()
        P = m.parent()
        if explicit_solutions:
            P.eval('solveexplicit: true')
        s = m.solve(x).str()
        if explicit_solutions:
            P.eval('solveexplicit: false')

        from sage.symbolic.relation import string_to_list_of_solutions

        X = string_to_list_of_solutions(s)


        #################
        # to_poly_solve #
        #################
        if explicit_solutions is not False:
            for eq in X:
                if repr(x) in map(repr, eq.rhs().variables()):
                    from sage.calculus.calculus import symbolic_expression_from_maxima_element
                    X = symbolic_expression_from_maxima_element(m.to_poly_solve(x))
                    X = [eq[0] for eq in X]
                    break

        if solution_dict is True:
            X=[dict([[sol.left(),sol.right()]]) for sol in X]

        if multiplicities:
            if len(X) == 0:
                return X, []
            else:
                return X, [int(e) for e in str(P.get('multiplicities'))[1:-1].split(',')]
        else:
            return X

    def find_root(self, a, b, var=None, xtol=10e-13, rtol=4.5e-16, maxiter=100, full_output=False):
        """
        Numerically find a root of self on the closed interval [a,b] (or
        [b,a]) if possible, where self is a function in the one variable.

        INPUT:

        -  ``a, b`` - endpoints of the interval

        -  ``var`` - optional variable

        -  ``xtol, rtol`` - the routine converges when a root
           is known to lie within xtol of the value return. Should be = 0. The
           routine modifies this to take into account the relative precision
           of doubles.

        -  ``maxiter`` - integer; if convergence is not
           achieved in maxiter iterations, an error is raised. Must be = 0.

        -  ``full_output`` - bool (default: False), if True,
           also return object that contains information about convergence.


        EXAMPLES:

        Note that in this example both f(-2) and f(3) are positive,
        yet we still find a root in that interval::

            sage: f = x^2 - 1
            sage: f.find_root(-2, 3)
            1.0
            sage: f.find_root(-2, 3, x)
            1.0
            sage: z, result = f.find_root(-2, 3, full_output=True)
            sage: result.converged
            True
            sage: result.flag
            'converged'
            sage: result.function_calls
            11
            sage: result.iterations
            10
            sage: result.root
            1.0

        More examples::

            sage: (sin(x) + exp(x)).find_root(-10, 10)
            -0.588532743981862...
            sage: sin(x).find_root(-1,1)
            0.0
            sage: (1/tan(x)).find_root(3,3.5)
            3.1415926535...

        An example with a square root::

            sage: f = 1 + x + sqrt(x+2); f.find_root(-2,10)
            -1.6180339887498949

        Some examples that Ted Kosan came up with::

            sage: t = var('t')
            sage: v = 0.004*(9600*e^(-(1200*t)) - 2400*e^(-(300*t)))
            sage: v.find_root(0, 0.002)
            0.001540327067911417...

        ::

            sage: a = .004*(8*e^(-(300*t)) - 8*e^(-(1200*t)))*(720000*e^(-(300*t)) - 11520000*e^(-(1200*t))) +.004*(9600*e^(-(1200*t)) - 2400*e^(-(300*t)))^2

        There is a 0 very close to the origin::

            sage: show(plot(a, 0, .002), xmin=0, xmax=.002)

        Using solve does not work to find it::

            sage: a.solve(t)
            []

        However ``find_root`` works beautifully::

            sage: a.find_root(0,0.002)
            0.0004110514049349...

        We illustrate that root finding is only implemented in one
        dimension::

            sage: x, y = var('x,y')
            sage: (x-y).find_root(-2,2)
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding currently only implemented in 1 dimension.

        TESTS:

        Test the special case that failed for the first attempt to fix
        #3980::

            sage: t = var('t')
            sage: find_root(1/t - x,0,2)
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding currently only implemented in 1 dimension.
        """
        if is_a_relational(self._gobj) and self.operator() is not operator.eq:
            raise ValueError, "Symbolic equation must be an equality."
        from sage.numerical.optimize import find_root
        if self.number_of_arguments() == 0:
            if bool(self == 0):
                return a
            else:
                raise RuntimeError, "no zero in the interval, since constant expression is not 0."
        elif self.number_of_arguments() == 1:
            f = self._fast_float_(self.default_variable())
            return find_root(f, a=a, b=b, xtol=xtol,
                             rtol=rtol,maxiter=maxiter,
                             full_output=full_output)
        else:
            raise NotImplementedError, "root finding currently only implemented in 1 dimension."

    def find_maximum_on_interval(self, a, b, var=None, tol=1.48e-08, maxfun=500):
        r"""
        Numerically find the maximum of the expression ``self``
        on the interval [a,b] (or [b,a]) along with the point at which the
        maximum is attained.

        See the documentation for
        ``self.find_minimum_on_interval`` for more details.

        EXAMPLES::

            sage: f = x*cos(x)
            sage: f.find_maximum_on_interval(0,5)
            (0.5610963381910451, 0.8603335890...)
            sage: f.find_maximum_on_interval(0,5, tol=0.1, maxfun=10)
            (0.561090323458081..., 0.857926501456...)
        """
        minval, x = (-self).find_minimum_on_interval(a, b, var=var, tol=tol,
                                                     maxfun=maxfun)
        return -minval, x

    def find_minimum_on_interval(self, a, b, var=None, tol=1.48e-08, maxfun=500):
        r"""
        Numerically find the minimum of the expression ``self``
        on the interval [a,b] (or [b,a]) and the point at which it attains
        that minimum. Note that ``self`` must be a function of
        (at most) one variable.

        INPUT:

        -  ``var`` - variable (default: first variable in
           self)

        -  ``a,b`` - endpoints of interval on which to minimize
           self.

        -  ``tol`` - the convergence tolerance

        -  ``maxfun`` - maximum function evaluations


        OUTPUT:

        - ``minval`` - (float) the minimum value that self takes on in the
          interval [a,b]

        - ``x`` - (float) the point at which self takes on the minimum value

        EXAMPLES::

            sage: f = x*cos(x)
            sage: f.find_minimum_on_interval(1, 5)
            (-3.288371395590..., 3.4256184695...)
            sage: f.find_minimum_on_interval(1, 5, tol=1e-3)
            (-3.288371361890..., 3.4257507903...)
            sage: f.find_minimum_on_interval(1, 5, tol=1e-2, maxfun=10)
            (-3.288370845983..., 3.4250840220...)
            sage: show(f.plot(0, 20))
            sage: f.find_minimum_on_interval(1, 15)
            (-9.477294259479..., 9.5293344109...)

        ALGORITHM:

        Uses ``scipy.optimize.fminbound`` which uses Brent's method.

        AUTHORS:

        - William Stein (2007-12-07)
        """
        from sage.numerical.optimize import find_minimum_on_interval

        if var is None:
            var = self.default_variable()
        return find_minimum_on_interval(self._fast_float_(var),
                                        a=a, b=b, tol=tol, maxfun=maxfun )

    ###################
    # Fast Evaluation #
    ###################
    def _fast_float_(self, *vars):
        """
        Returns an object which provides fast floating point
        evaluation of this symbolic expression.

        See :mod:`sage.ext.fast_eval` for more information.

        EXAMPLES::

            sage: f = sqrt(x+1)
            sage: ff = f._fast_float_('x')
            sage: ff(1.0)
            1.4142135623730951
            sage: type(_)
            <type 'float'>
        """
        from sage.symbolic.expression_conversions import fast_float
        return fast_float(self, *vars)

    def _fast_callable_(self, etb):
        """
        Given an ExpressionTreeBuilder *etb*, return an Expression representing
        this symbolic expression.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x','y'])
            sage: x,y = var('x,y')
            sage: f = y+2*x^2
            sage: f._fast_callable_(etb)
            add(mul(ipow(v_0, 2), 2), v_1)
        """
        from sage.symbolic.expression_conversions import fast_callable
        return fast_callable(self, etb)

    def show(self):
        """
        Show this symbolic expression, i.e., typeset it nicely.

        EXAMPLES::

            sage: (x^2 + 1).show()
            x^{2}  + 1
        """
        from sage.misc.functional import _do_show
        return _do_show(self)

    def plot(self, *args, **kwds):
        """
        Plot a symbolic expression. All arguments are passed onto the standard plot command.

        EXAMPLES:

        This displays a straight line::

            sage: sin(2).plot((x,0,3))

        This draws a red oscillatory curve::

            sage: sin(x^2).plot((x,0,2*pi), rgbcolor=(1,0,0))

        Another plot using the variable theta::

            sage: var('theta')
            theta
            sage: (cos(theta) - erf(theta)).plot((theta,-2*pi,2*pi))

        A very thick green plot with a frame::

            sage: sin(x).plot((x,-4*pi, 4*pi), thickness=20, rgbcolor=(0,0.7,0)).show(frame=True)

        You can embed 2d plots in 3d space as follows::

            sage: plot(sin(x^2), (x,-pi, pi), thickness=2).plot3d(z = 1)

        A more complicated family::

            sage: G = sum([plot(sin(n*x), (x,-2*pi, 2*pi)).plot3d(z=n) for n in [0,0.1,..1]])
            sage: G.show(frame_aspect_ratio=[1,1,1/2])

        A plot involving the floor function::

            sage: plot(1.0 - x * floor(1/x), (x,0.00001,1.0))

        Sage used to allow symbolic functions with "no arguments";
        this no longer works::

            sage: plot(2*sin, -4, 4)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Integer Ring' and '<class 'sage.functions.trig.Function_sin'>'

        You should evaluate the function first::

            sage: plot(2*sin(x), -4, 4)

        TESTS::

            sage: f(x) = x*(1 - x)
            sage: plot(f,0,1)
        """
        from sage.symbolic.callable import is_CallableSymbolicExpression
        from sage.symbolic.ring import is_SymbolicVariable
        from sage.plot.plot import plot

        # see if the user passed a variable in.
        if kwds.has_key('param'):
            param = kwds['param']
        else:
            param = None
            for i, arg in enumerate(args):
                if is_SymbolicVariable(arg):
                    param = arg
                    args = args[:i] + args[i+1:]
                    break

        if param is None:
            if is_CallableSymbolicExpression(self):
                A = self.arguments()
                if len(A) == 0:
                    raise ValueError, "function has no input arguments"
                else:
                    param = A[0]

                f = self._fast_float_(param)
            else:
                A = self.variables()
                if len(A) == 0:
                    #Here we handle the case where f is something
                    #like 2*sin, which has takes arguments which
                    #aren't explicitly given
                    n = self.number_of_arguments()
                    f = self._fast_float_()
                else:
                    param = A[0]
                    try:
                        f = self._fast_float_(param)
                    except NotImplementedError:
                        return self.function(param)
        else:
            try:
                f = self._fast_float_(param)
            except NotImplementedError:
                return self.function(param)
        return plot(f, *args, **kwds)

    ############
    # Calculus #
    ############
    def integral(self, *args, **kwds):
        """
        Compute the integral of self.  Please see
        :obj:`sage.calculus.calculus.integral` for more details.

        EXAMPLES::

            sage: sin(x).integral(x,0,3)
            -cos(3) + 1
            sage: sin(x).integral(x)
            -cos(x)
        """
        from sage.calculus.calculus import integral
        return integral(self, *args, **kwds)

    integrate = integral

    def nintegral(self, *args, **kwds):
        """
        Compute the numerical integral of self.  Please see
        :obj:`sage.calculus.calculus.nintegral` for more details.

        EXAMPLES::

            sage: sin(x).nintegral(x,0,3)
            (1.989992496600..., 2.209335488557...e-14, 21, 0)
        """
        from sage.calculus.calculus import nintegral
        return nintegral(self, *args, **kwds)

    nintegrate = nintegral

    def minpoly(self, *args, **kwds):
        """
        Return the minimal polynomial of this symbolic expression.

        EXAMPLES::

            sage: golden_ratio.minpoly()
            x^2 - x - 1
        """
        try:
            obj = self.pyobject()
            return obj.minpoly()
        except AttributeError:
            pass
        except TypeError:
            pass
        from sage.calculus.calculus import minpoly
        return minpoly(self, *args, **kwds)

    def limit(self, *args, **kwds):
        """
        Return a symbolic limit.  See
        :obj:`sage.calculus.calculus.limit`

        EXAMPLES::

            sage: (sin(x)/x).limit(x=0)
            1
        """
        from sage.calculus.calculus import limit
        return limit(self, *args, **kwds)

    def laplace(self, t, s):
        """
        Return Laplace transform of self.  See
        :obj:`sage.calculus.calculus.laplace`

        EXAMPLES::

            sage: var('x,s,z')
            (x, s, z)
            sage: (z + exp(x)).laplace(x, s)
            z/s + 1/(s - 1)
        """
        from sage.calculus.calculus import laplace
        return laplace(self, t, s)

    def inverse_laplace(self, t, s):
        """
        Return inverse Laplace transform of self.  See
        :obj:`sage.calculus.calculus.inverse_laplace`

        EXAMPLES::

            sage: var('w, m')
            (w, m)
            sage: f = (1/(w^2+10)).inverse_laplace(w, m); f
            1/10*sqrt(10)*sin(sqrt(10)*m)
        """
        from sage.calculus.calculus import inverse_laplace
        return inverse_laplace(self, t, s)

    def add_to_both_sides(self, x):
        """
        Returns a relation obtained by adding *x* to both sides of
        this relation.

        EXAMPLES::

            sage: var('x y z')
            (x, y, z)
            sage: eqn = x^2 + y^2 + z^2 <= 1
            sage: eqn.add_to_both_sides(-z^2)
            x^2 + y^2 <= -z^2 + 1
            sage: eqn.add_to_both_sides(I)
            x^2 + y^2 + z^2 + I <= (I + 1)
        """
        if not is_a_relational(self._gobj):
            raise TypeError, "this expression must be a relation"
        return self + x

    def subtract_from_both_sides(self, x):
        """
        Returns a relation obtained by subtracting *x* from both sides
        of this relation.

        EXAMPLES::

            sage: eqn = x*sin(x)*sqrt(3) + sqrt(2) > cos(sin(x))
            sage: eqn.subtract_from_both_sides(sqrt(2))
            sqrt(3)*x*sin(x) > -sqrt(2) + cos(sin(x))
            sage: eqn.subtract_from_both_sides(cos(sin(x)))
            sqrt(3)*x*sin(x) + sqrt(2) - cos(sin(x)) > 0
        """
        if not is_a_relational(self._gobj):
            raise TypeError, "this expression must be a relation"
        return self - x

    def multiply_both_sides(self, x, checksign=None):
        """
        Returns a relation obtained by multiplying both sides of this
        relation by *x*.

        .. note::

           The *checksign* keyword argument is currently ignored and
           is included for backward compatibility reasons only.

        EXAMPLES::

            sage: var('x,y'); f = x + 3 < y - 2
            (x, y)
            sage: f.multiply_both_sides(7)
            7*x + 21 < 7*y - 14
            sage: f.multiply_both_sides(-1/2)
            -1/2*x - 3/2 < -1/2*y + 1
            sage: f*(-2/3)
            -2/3*x - 2 < -2/3*y + 4/3
            sage: f*(-pi)
            -(x + 3)*pi < -(y - 2)*pi

        Since the direction of the inequality never changes when doing
        arithmetic with equations, you can multiply or divide the
        equation by a quantity with unknown sign::

            sage: f*(1+I)
            (I + 1)*x + 3*I + 3 < (I + 1)*y - 2*I - 2
            sage: f = sqrt(2) + x == y^3
            sage: f.multiply_both_sides(I)
            I*x + I*sqrt(2) == I*y^3
            sage: f.multiply_both_sides(-1)
            -x - sqrt(2) == -y^3

        Note that the direction of the following inequalities is
        not reversed::

            sage: (x^3 + 1 > 2*sqrt(3)) * (-1)
            -x^3 - 1 > -2*sqrt(3)
            sage: (x^3 + 1 >= 2*sqrt(3)) * (-1)
            -x^3 - 1 >= -2*sqrt(3)
            sage: (x^3 + 1 <= 2*sqrt(3)) * (-1)
            -x^3 - 1 <= -2*sqrt(3)
        """
        if not is_a_relational(self._gobj):
            raise TypeError, "this expression must be a relation"
        return self * x

    def divide_both_sides(self, x, checksign=None):
        """
        Returns a relation obtained by dividing both sides of this
        relation by *x*.

        .. note::

           The *checksign* keyword argument is currently ignored and
           is included for backward compatibility reasons only.

        EXAMPLES::

            sage: theta = var('theta')
            sage: eqn =   (x^3 + theta < sin(x*theta))
            sage: eqn.divide_both_sides(theta, checksign=False)
            (x^3 + theta)/theta < sin(theta*x)/theta
            sage: eqn.divide_both_sides(theta)
            (x^3 + theta)/theta < sin(theta*x)/theta
            sage: eqn/theta
            (x^3 + theta)/theta < sin(theta*x)/theta
        """
        if not is_a_relational(self._gobj):
            raise TypeError, "this expression must be a relation"
        return self / x


    # Functions to add later, maybe.  These were in Ginac mainly
    # implemented using a lot from cln, and I had to mostly delete
    # their implementations.   They are pretty specialized for
    # physics apps, maybe.
    # This doesn't work / isn't implemented yet / just segfaults.
    #def Li(self, x):
    #    """
    #    """
    #    cdef Expression nexp = self.coerce_in(x)
    #    return new_Expression_from_GEx(self._parent, g_Li(self._gobj, nexp._gobj))
    #def Li2(self):
    #    return new_Expression_from_GEx(self._parent, g_Li2(self._gobj))
    #def G(self, Expression y):
    #    return new_Expression_from_GEx(self._parent, g_G(self._gobj, y._gobj))
    #def G2(self, Expression s, Expression y):
    #    return new_Expression_from_GEx(self._parent, g_G2(self._gobj, s._gobj, y._gobj))
    #def SR(self, Expression p, Expression x):
    #return new_Expression_from_GEx(self._parent, g_SR(self._gobj, p._gobj, x._gobj))
    #def H(self, Expression x):
    #return new_Expression_from_GEx(self._parent, g_H(self._gobj, x._gobj))
    #def zeta2(self, Expression s):
    #    return new_Expression_from_GEx(self._parent, g_zeta2(self._gobj, s._gobj))
    #def zetaderiv(self, Expression x):
    #    return new_Expression_from_GEx(self._parent, g_zetaderiv(self._gobj, x._gobj))
    #def beta(self, Expression y):
    #    return new_Expression_from_GEx(self._parent, g_beta(self._gobj, y._gobj))
    #def psi(self):
    #    return new_Expression_from_GEx(self._parent, g_psi(self._gobj))
    #def psi2(self, Expression x):
    #    return new_Expression_from_GEx(self._parent, g_psi2(self._gobj, x._gobj))



cdef Expression new_Expression_from_GEx(parent, GEx juice):
    cdef Expression nex
    nex = <Expression>PY_NEW(Expression)
    GEx_construct_ex(&nex._gobj, juice)
    nex._parent = parent
    return nex

cdef Expression new_Expression_from_pyobject(parent, x):
    cdef GEx exp
    GEx_construct_pyobject(exp, x)
    return new_Expression_from_GEx(parent, exp)

cdef class ExpressionIterator:
    cdef Expression _ex
    cdef int _ind
    cdef int _len
    def __iter__(self):
        """
        Return this iterator object itself.

        EXAMPLES::

            sage: x,y,z = var('x,y,z')
            sage: i = (x+y).iterator()
            sage: iter(i) is i
            True
        """
        return self

    def __next__(self):
        """
        Return the next component of the expression.

        EXAMPLES::

            sage: x,y,z = var('x,y,z')
            sage: i = (x+y).iterator()
            sage: i.next()
            x
        """
        cdef GEx ex
        if self._ind == self._len:
            raise StopIteration
        ex = self._ex._gobj.op(self._ind)
        self._ind+=1
        return new_Expression_from_GEx(self._ex._parent, ex)

cdef inline ExpressionIterator new_ExpIter_from_Expression(Expression ex):
    """
    Construct a new iterator over a symbolic expression.

    EXAMPLES::

        sage: x,y,z = var('x,y,z')
        sage: i = (x+y).iterator() #indirect doctest
    """
    # The const_iterator in GiNaC just keeps an integer index to the current
    # subexpression. We do the same here, to avoid the trouble of having to
    # mess with C++ class constructors/destructors.
    cdef ExpressionIterator m = <ExpressionIterator>PY_NEW(ExpressionIterator)
    m._ex = ex
    m._ind = 0
    m._len  = ex._gobj.nops()
    return m


cdef operators compatible_relation(operators lop, operators rop) except <operators>-1:
    """
    TESTS::

        sage: var('a,b,x,y')
        (a, b, x, y)
        sage: (x < a) + (y <= b)     # indirect doctest
        x + y < a + b
        sage: (x >= 4) * (y > 7)
        x*y > 28
    """
    if lop == rop:
        return lop
    elif lop == not_equal or rop == not_equal:
        raise TypeError, "incompatible relations"
    elif lop == equal:
       return rop
    elif rop == equal:
       return lop
    elif (lop in [less, less_or_equal] and rop in [less, less_or_equal]):
       return less
    elif (lop in [greater, greater_or_equal] and rop in [greater, greater_or_equal]):
       return greater
    else:
        raise TypeError, "incompatible relations"
