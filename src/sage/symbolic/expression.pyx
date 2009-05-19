###############################################################################
#   SAGE: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

"""
EXAMPLES:

RELATIONAL EXPRESSIONS:

We create a relational expression:
    sage: x = var('x',ns=1); S = x.parent()
    sage: eqn = (x-1)^2 <= x^2 - 2*x + 3
    sage: eqn.subs(x == 5)
    16 <= 18

Notice that squaring the relation squares both sides.
    sage: eqn^2
    (x - 1)^4 <= (x^2 - 2*x + 3)^2
    sage: eqn.expand()
    x^2 - 2*x + 1 <= x^2 - 2*x + 3

The can transform a true relational into a false one:
    sage: eqn = S(-5) < S(-3); eqn
    -5 < -3
    sage: bool(eqn)
    True
    sage: eqn^2
    25 < 9
    sage: bool(eqn^2)
    False

We can do arithmetic with relationals:
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

We can even add together two relations, so long as the operators are the same:
    sage: (x^3 + x <= x - 17)  + (-x <= x - 10)
    x^3 <= 2*x - 27

Here they aren't:
    sage: (x^3 + x <= x - 17)  + (-x >= x - 10)
    Traceback (most recent call last):
    ...
    TypeError: incompatible relations


ANY SAGE ELEMENTS:

You can work symbolically with any Sage datatype.  This can lead to nonsense
if the data type is strange, e.g., an element of a finite field (at present).

We mix Singular variables with symbolic variables:
    sage: R.<u,v> = QQ[]
    sage: var('a,b,c', ns=1)
    (a, b, c)
    sage: expand((u + v + a + b + c)^2)
    a^2 + 2*a*b + 2*a*c + b^2 + 2*b*c + c^2 + (2*u + 2*v)*a + (2*u + 2*v)*b + (2*u + 2*v)*c + u^2 + 2*u*v + v^2

TESTS:

Test jacobian on pynac expressions. #5546 ::

    sage: var('x,y', ns=1)
    (x, y)
    sage: f = x + y
    sage: jacobian(f, [x,y])
    [1 1]


Test if matrices work #5546 ::

    sage: var('x,y,z', ns=1)
    (x, y, z)
    sage: M = matrix(2,2,[x,y,z,x])
    sage: v = vector([x,y])
    sage: M * v
    (x^2 + y^2, x*y + x*z)
    sage: v*M
    (x^2 + y*z, 2*x*y)

"""

include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"

import ring
import sage.rings.integer
import sage.rings.rational

from sage.structure.element cimport ModuleElement, RingElement, Element
from sage.symbolic.function import get_sfunction_from_serial


from sage.rings.rational import Rational  # Used for sqrt.

from sage.calculus.calculus import CallableSymbolicExpressionRing

from sage.misc.derivative import multi_derivative

cdef class Expression(CommutativeRingElement):
    cpdef object pyobject(self):
        """
        Get the underlying Python object corresponding to this
        expression, assuming this expression is a single numerical
        value.   Otherwise, a TypeError is raised.

        EXAMPLES:
            sage: var('x',ns=1); S = parent(x)
            x
            sage: b = -17/3
            sage: a = S(b)
            sage: a.pyobject()
            -17/3
            sage: a.pyobject() is b
            True
        """
        if not is_a_numeric(self._gobj):
            raise TypeError, "self must be a numeric expression"
        return py_object_from_numeric(self._gobj)

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
            sage: var('x,y,z',ns=1)
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
            sage: var('x,y,z',ns=1)
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
        self._set_parent(ring.NSR)
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

        EXAMPLES:
            sage: var("x y", ns=1)
            (x, y)
            sage: (x+y)._repr_()
            'x + y'

        TESTS:
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
            1/x^5
            sage: x^(-y)
            1/x^y
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
        """
        return GEx_to_str(&self._gobj)

    def _latex_(self):
        r"""
        Return string representation of this symbolic expression.

        EXAMPLES:


        TESTS:
            sage: var('x,y,z',ns=1)
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
            {\left|x\right|}
            sage: latex((x*y).conjugate())
            \bar{x} \bar{y}

        Check spacing of coefficients of mul expressions (#3202)::

            sage: latex(2*3^x)
            2 \, 3^{x}

        """
        return GEx_to_str_latex(&self._gobj)

    def _integer_(self, ZZ=None):
        """
        EXAMPLES:
            sage: var('x',ns=1); S = parent(x)
            x
            sage: f = x^3 + 17*x -3
            sage: ZZ(f.coeff(x^3))
            1
            sage: ZZ(f.coeff(x))
            17
            sage: ZZ(f.coeff(x,0))
            -3
            sage: type(ZZ(f.coeff(x,0)))
            <type 'sage.rings.integer.Integer'>

        Coercion is done if necessary:
            sage: f = x^3 + 17/1*x
            sage: ZZ(f.coeff(x))
            17
            sage: type(ZZ(f.coeff(x)))
            <type 'sage.rings.integer.Integer'>

        If the symbolic expression is just a wrapper around an integer,
        that very same integer is returned:
            sage: n = 17; S(n)._integer_() is n
            True
        """
        n = self.pyobject()
        if isinstance(n, sage.rings.integer.Integer):
            return n
        return sage.rings.integer.Integer(n)

    def _rational_(self):
        """
        EXAMPLES:
            sage: var('x',ns=1); S = parent(x)
            x
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
        that very same rational is returned:
            sage: n = 17/1; S(n)._rational_() is n
            True
        """
        n = self.pyobject()
        if isinstance(n, sage.rings.rational.Rational):
            return n
        return sage.rings.rational.Rational(n)

    cpdef _eval_self(self, R):
        """
        Evaluate this expression numerically, and try to coerce it to R.

        EXAMPLES::

            sage: var('x,y,z',ns=1)
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

    def _mpfr_(self, R):
        """
        Return a numerical approximation to this expression in the RealField R.

        The precision of the approximation is determined by the precision of
        the input R.

        EXAMPLES::

            sage: var('x',ns=1); S = parent(x)
            x
            sage: RealField(200)(S(1/11))
            0.090909090909090909090909090909090909090909090909090909090909

            sage: a = S(3).sin(); a
            sin(3)
            sage: RealField(200)(a)
            0.14112000805986722210074480280811027984693326425226558415188
            sage: a._mpfr_(RealField(100))
            0.14112000805986722210074480281
        """
        return self._eval_self(R)

    def _complex_mpfr_field_(self, R):
        """
        Return a numerical approximation to this expression in the given
        ComplexField R.

        The precision of the approximation is determined by the precision of
        the input R.

        EXAMPLES::

            sage: var('x',ns=1); S = parent(x)
            x
            sage: ComplexField(200)(S(1/11))
            0.090909090909090909090909090909090909090909090909090909090909
            sage: zeta(x).subs(x=I)._complex_mpfr_field_(ComplexField(70))
            0.0033002236853241028742 - 0.41815544914132167669*I
            sage: gamma(x).subs(x=I)._complex_mpfr_field_(ComplexField(60))
            -0.15494982830181069 - 0.49801566811835604*I
            sage: log(x).subs(x=I)._complex_mpfr_field_(ComplexField(50))
            1.5707963267949*I

        """
        return self._eval_self(R)

    def _complex_double_(self, R):
        """
        Return a numerical approximation to this expression in the given
        Complex Double Field R.

        EXAMPLES::

            sage: var('x',ns=1); S = parent(x)
            x
            sage: CDF(S(1/11))
            0.0909090909091
            sage: zeta(x).subs(x=I)._complex_double_(CDF)
            0.00330022368532 - 0.418155449141*I
            sage: gamma(x).subs(x=I)._complex_double_(CDF)
            -0.154949828302 - 0.498015668118*I
            sage: log(x).subs(x=I)._complex_double_(CDF)
            1.57079632679*I
        """
        return self._eval_self(R)

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
            2.0
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

        TESTS:

        Test if exceptions during hashing are handled properly::

            sage: t = S(matrix(2,2,range(4)))
            sage: hash(t)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
        """
        return self._gobj.gethash()

    # Boilerplate code from sage/structure/element.pyx
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
        return (<Element>left)._richcmp(right, op)

    cdef _richcmp_c_impl(left, Element right, int op):
        cdef Expression l, r

        l = left
        r = right

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

    cpdef bint is_polynomial(self, var):
        """
        Return True if self is a polynomial in the given variable.

        EXAMPLES::

            sage: var('x,y,z',ns=1)
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

        EXAMPLES:
            sage: x = var('x',ns=1); SR = x.parent()
            sage: eqn = (x-1)^2 == x^2 - 2*x + 3
            sage: eqn.is_relational()
            True
            sage: sin(x).is_relational()
            False
        """
        return is_a_relational(self._gobj)

    def lhs(self):
        """
        If self is a relational expression, return the left hand side
        of the relation.  Otherwise, raise a ValueError.

        EXAMPLES:
            sage: x = var('x',ns=1); SR = x.parent()
            sage: eqn = (x-1)^2 == x^2 - 2*x + 3
            sage: eqn.lhs()
            (x - 1)^2
        """
        if not self.is_relational():
            raise ValueError, "self must be a relational expression"
        return new_Expression_from_GEx(self._gobj.lhs())

    def rhs(self):
        """
        If self is a relational expression, return the right hand side of the relation.  Otherwise,
        raise a ValueError.

        EXAMPLES:
            sage: x = var('x',ns=1); SR = x.parent()
            sage: eqn = (x-1)^2 <= x^2 - 2*x + 3
            sage: eqn.rhs()
            x^2 - 2*x + 3
        """
        if not self.is_relational():
            raise ValueError, "self must be a relation"
        return new_Expression_from_GEx(self._gobj.rhs())

    def __nonzero__(self):
        """
        Return True if self is nonzero.

        EXAMPLES:
            sage: x = var('x',ns=1); SR = x.parent()
            sage: SR(0).__nonzero__()
            False
            sage: SR(1).__nonzero__()
            True
            sage: bool(abs(x))
            True
            sage: bool(x/x - 1)
            False

        A bunch of tests of nonzero (which is called by bool) for
        symbolic relations:
            sage: x = var('x',ns=1); SR = x.parent()
            sage: bool((x-1)^2 == x^2 - 2*x + 1)
            False
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
        """
        if self.is_relational():
            return relational_to_bool(self._gobj)
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

    cpdef ModuleElement _add_(left, ModuleElement right):
        """
        Add left and right.

        EXAMPLES;
            sage: var("x y", ns=1)
            (x, y)
            sage: x + y + y + x
            2*x + 2*y

            # adding relational expressions
            sage: ( (x+y) > x ) + ( x > y )
            2*x + y > x + y

            sage: ( (x+y) > x ) + x
            2*x + y > 2*x

        TESTS:
            sage: x + ( (x+y) > x )
            2*x + y > 2*x

            sage: ( x > y) + (y < x)
            Traceback (most recent call last):
            ...
            TypeError: incompatible relations

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
        if is_a_relational(left._gobj):
            if is_a_relational(_right._gobj):
                if relational_operator(left._gobj) != relational_operator(_right._gobj):
                    raise TypeError, "incompatible relations"
                x = relational(gadd(left._gobj.lhs(), _right._gobj.lhs()),
                               gadd(left._gobj.rhs(), _right._gobj.rhs()),
                               relational_operator(left._gobj))
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
        return new_Expression_from_GEx(x)

    cpdef ModuleElement _sub_(left, ModuleElement right):
        """
        EXAMPLES:
            sage: var("x y", ns=1)
            (x, y)
            sage: x - y
            x - y

            # subtracting relational expressions
            sage: ( (x+y) > x ) - ( x > y )
            y > x - y

            sage: ( (x+y) > x ) - x
            y > 0

        TESTS:
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
                if relational_operator(left._gobj) != relational_operator(_right._gobj):
                    raise TypeError, "incompatible relations"
                x = relational(gsub(left._gobj.lhs(), _right._gobj.lhs()),
                               gsub(left._gobj.rhs(), _right._gobj.rhs()),
                               relational_operator(left._gobj))
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
        return new_Expression_from_GEx(x)

    cpdef RingElement _mul_(left, RingElement right):
        """
        Multiply left and right.

        EXAMPLES:
            sage: var("x y", ns=1)
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

        TESTS:
            sage: x * ( (x+y) > x )
            (x + y)*x > x^2

            sage: ( x > y) * (y < x)
            Traceback (most recent call last):
            ...
            TypeError: incompatible relations

            sage: a = 1000 + 300*x + x^3 + 30*x^2
            sage: a*Mod(1,7)
            x^3 + 2*x^2 + 6*x + 6

            sage: var('z',ns=1)
            z
            sage: 3*(x+y)/z
            3*(x + y)/z
            sage: (-x+z)*(3*x-3*z)
            -3*(x - z)^2

            # check if comparison of constant terms in pynac add objects work
            sage: (y-1)*(y-2)
            (y - 2)*(y - 1)

            sage: x*oo
            +Infinity
            sage: x*unsigned_infinity
            Infinity
            sage: -x*oo
            -Infinity

            sage: nsr = x.parent()
            sage: nsr(oo)*nsr(oo)
            +Infinity
            sage: nsr(-oo)*nsr(oo)
            -Infinity
            sage: nsr(oo)*nsr(-oo)
            -Infinity
            sage: nsr(unsigned_infinity)*nsr(oo)
            Infinity

        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        cdef operators o
        if is_a_relational(left._gobj):
            if is_a_relational(_right._gobj):
                if relational_operator(left._gobj) != relational_operator(_right._gobj):
                    raise TypeError, "incompatible relations"
                x = relational(gmul(left._gobj.lhs(), _right._gobj.lhs()),
                               gmul(left._gobj.rhs(), _right._gobj.rhs()),
                               relational_operator(left._gobj))
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
        return new_Expression_from_GEx(x)

    cpdef RingElement _div_(left, RingElement right):
        """
        Divide left and right.

        EXAMPLES:
            sage: var("x y", ns=1)
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

        TESTS:
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

            sage: nsr = x.parent()
            sage: nsr(oo)/nsr(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0*infinity encountered.

            sage: nsr(-oo)/nsr(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0*infinity encountered.

            sage: nsr(oo)/nsr(-oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0*infinity encountered.

            sage: nsr(oo)/nsr(unsigned_infinity)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0*infinity encountered.

            sage: nsr(unsigned_infinity)/nsr(oo)
            Traceback (most recent call last):
            ...
            RuntimeError: indeterminate expression: 0*infinity encountered.

            sage: nsr(0)/nsr(oo)
            0

            sage: nsr(0)/nsr(unsigned_infinity)
            0


        """
        cdef GEx x
        cdef Expression _right = <Expression>right
        cdef operators o
        if is_a_relational(left._gobj):
            if is_a_relational(_right._gobj):
                if relational_operator(left._gobj) != relational_operator(_right._gobj):
                    raise TypeError, "incompatible relations"
                x = relational(gdiv(left._gobj.lhs(), _right._gobj.lhs()),
                               gdiv(left._gobj.rhs(), _right._gobj.rhs()),
                               relational_operator(left._gobj))
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
        return new_Expression_from_GEx(x)

    # Boilerplate code from sage/structure/element.pyx
    def __cmp__(left, right):
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
            sage: cmp(S(0.5), 0.7)
            -1
            sage: cmp(sin(S(2)), sin(S(1)))
            1
            sage: float(sin(S(2)))
            0.90929742682568171
            sage: float(sin(S(1)))
            0.8414709848078965
        """
        return (<Element>left)._cmp(right)

    cdef int _cmp_c_impl(left, Element right) except -2:
        return left._gobj.compare((<Expression>right)._gobj)

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

        TESTS:
            sage: (Mod(2,7)*x^2 + Mod(2,7))^7
            (2*x^2 + 2)^7
            sage: k = GF(7)
            sage: f = expand((k(1)*x^5 + k(1)*x^2 + k(2))^7); f
            x^35 + x^14 + 2

            sage: x^oo
            Traceback (most recent call last):
            ...
            RuntimeError: power::eval(): pow(x, Infinity) for non numeric x is not defined.
            sage: S(oo)^2
            +Infinity
            sage: S(-oo)^2
            +Infinity
            sage: S(-oo)^3
            -Infinity
            sage: S(unsigned_infinity)^2
            Infinity

        Test powers of exp::

            sage: S(2).exp()^5
            e^10
            sage: x.exp()^5
            e^(5*x)
        """
        cdef Expression nexp = self.coerce_in(exp)
        cdef GEx x
        if is_a_relational(self._gobj):
            x = relational(g_pow(self._gobj.lhs(), nexp._gobj),
                           g_pow(self._gobj.rhs(), nexp._gobj),
                           relational_operator(self._gobj))
        else:
            x = g_pow(self._gobj, nexp._gobj)
        return new_Expression_from_GEx(x)

    def derivative(self, *args):
        """
        Returns the derivative of this expressions with respect to the variables
        supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more details.

        .. seealso::

            :meth:`_derivative`

        EXAMPLES::
            sage: var("x y", ns=1)
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

        TESTS:
            sage: t.derivative()
            Traceback (most recent call last):
            ...
            ValueError: No differentiation variable specified.
        """
        return multi_derivative(self, args)

    def _derivative(self, symb=None, deg=1):
        """
        Return the deg-th (partial) derivative of self with respect to symb.

        EXAMPLES:
            sage: var("x y", ns=1)
            (x, y)
            sage: b = (x+y)^5
            sage: b._derivative(x, 2)
            20*(x + y)^3

            sage: from sage.symbolic.function import function as myfunc
            sage: foo = myfunc('foo',2)
            sage: foo(x^2,x^2)._derivative(x)
            2*x*D[0](foo)(x^2, x^2) + 2*x*D[1](foo)(x^2, x^2)

        TESTS:
            Raise error if no variable is specified::
            sage: b._derivative()
            Traceback (most recent call last):
            ...
            ValueError: No differentiation variable specified.
        """
        if symb is None:
            # we specify a default value of None for symb and check for it here
            # to return more helpful error messages when no variable is
            # given by the multi_derivative framework
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
        return new_Expression_from_GEx(x)

    def series(self, symbol, int order):
        r"""
        Return the power series expansion of self in terms of the variable
        symbol to the given order.

        INPUT:
            symbol -- a variable
            order -- an integer

        OUTPUT:
            a power series --

        To truncate the power series and obtain a normal expression, use the
        truncate command.

        EXAMPLES:
        We expand a polynomial in $x$ about 0, about $1$, and also truncate
        it back to a polynomial:
            sage: var('x,y',ns=1)
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

        We computer another series expansion of an analytic function:
            sage: f = sin(x)/x^2
            sage: f.series(x,7)
            1*x^(-1) + (-1/6)*x + 1/120*x^3 + (-1/5040)*x^5 + Order(x^7)
            sage: f.series(x==1,3)
            (sin(1)) + (-2*sin(1) + cos(1))*(x - 1) + (5/2*sin(1) - 2*cos(1))*(x - 1)^2 + Order((x - 1)^3)
            sage: f.series(x==1,3).truncate().expand()
            5/2*x^2*sin(1) - 2*x^2*cos(1) - 7*x*sin(1) + 5*x*cos(1) + 11/2*sin(1) - 3*cos(1)

        Following the GiNaC tutorial, we use John Machin's amazing
        formula $\pi = 16 \atan(1/5) - 4 \atan(1/239)$ to compute
        digits of $\pi$. We expand the arc tangent around 0 and insert
        the fractions 1/5 and 1/239.
            sage: x = var('x',ns=1)
            sage: f = atan(x).series(x, 10); f
            1*x + (-1/3)*x^3 + 1/5*x^5 + (-1/7)*x^7 + 1/9*x^9 + Order(x^10)
            sage: float(16*f.subs(x==1/5) - 4*f.subs(x==1/239))
            3.1415926824043994
        """
        cdef Expression symbol0 = self.coerce_in(symbol)
        _sig_on
        cdef GEx x = self._gobj.series(symbol0._gobj, order, 0)
        _sig_off
        return new_Expression_from_GEx(x)

    def truncate(self):
        """
        Given a power series or expression, return the corresponding
        expression without the big oh.

        INPUT:
            a series as output by the series command

        OUTPUT:
            expression

        EXAMPLES:
            sage: var('x,y',ns=1)
            (x, y)
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
        return new_Expression_from_GEx(series_to_poly(self._gobj))

    def expand(Expression self):
        """
        Return expanded form of this expression, obtained by multiplying out
        all products.

        OUTPUT:
            symbolic expression

        EXAMPLES:
            sage: var('x,y',ns=1)
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
        _sig_on
        cdef GEx x = self._gobj.expand(0)
        _sig_off
        return new_Expression_from_GEx(x)

    ############################################################################
    # Pattern Matching
    ############################################################################
    def match(self, pattern):
        """
        See http://www.ginac.de/tutorial/Pattern-matching-and-advanced-substitutions.html

        EXAMPLES:
            sage: var('x,y,z,a,b,c,d,e,f',ns=1); S = parent(x)
            (x, y, z, a, b, c, d, e, f)
            sage: w0 = S.wild(0); w1 = S.wild(1); w2 = S.wild(2)
            sage: ((x+y)^a).match((x+y)^a)
            True
            sage: ((x+y)^a).match((x+y)^b)
            False
            sage: ((x+y)^a).match(w0^w1)
            True
            sage: ((x+y)^a).match(w0^w0)
            False
            sage: ((x+y)^(x+y)).match(w0^w0)
            True
            sage: ((a+b)*(a+c)).match((a+w0)*(a+w1))
            True
            sage: ((a+b)*(a+c)).match((w0+b)*(w0+c))
            True
            sage: ((a+b)*(a+c)).match((w0+w1)*(w0+w2))    # surprising?
            False
            sage: (a*(x+y)+a*z+b).match(a*w0+w1)
            True
            sage: (a+b+c+d+e+f).match(c)
            False
            sage: (a+b+c+d+e+f).has(c)
            True
            sage: (a+b+c+d+e+f).match(c+w0)
            True
            sage: (a+b+c+d+e+f).match(c+e+w0)
            True
            sage: (a+b).match(a+b+w0)
            True
            sage: (a*b^2).match(a^w0*b^w1)
            False
            sage: (a*b^2).match(a*b^w1)
            True
            sage: (x*x.arctan2(x^2)).match(w0*w0.arctan2(w0^2))
            True
        """
        cdef Expression p = self.coerce_in(pattern)
        return self._gobj.match(p._gobj)

    def find(self, pattern):
        """
        Find all occurences of the given pattern in this expression.

        Note that once a subexpression matches the pattern, the search doesn't
        extend to subexpressions of it.

        EXAMPLES::

            sage: var('x,y,z,a,b',ns=1); S = parent(x)
            (x, y, z, a, b)
            sage: w0 = S.wild(0); w1 = S.wild(1)

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
            res.append(new_Expression_from_GEx(itr.obj()))
            itr.inc()
        return res

    def has(self, pattern):
        """
        EXAMPLES:
            sage: var('x,y,a', ns=1); S = x.parent(); w0 = S.wild(); w1 = S.wild()
            (x, y, a)
            sage: (x*sin(x + y + 2*a)).has(y)
            True

        Here "x+y" is not a subexpression of "x+y+2*a" (which has the
        subexpressions "x", "y" and "2*a"):
            sage: (x*sin(x + y + 2*a)).has(x+y)
            False
            sage: (x*sin(x + y + 2*a)).has(x + y + w0)
            True

        The following fails because "2*(x+y)" automatically gets converted to
        "2*x+2*y" of which "x+y" is not a subexpression:
            sage: (x*sin(2*(x+y) + 2*a)).has(x+y)
            False

        Although x^1==x and x^0==1, neither "x" nor "1" are actually of the
        form "x^something":
            sage: (x+1).has(x^w0)
            False

        Here is another possible pitfall, where the first expression
        matches because the term "-x" has the form "(-1)*x" in GiNaC. To check
        whether a polynomial contains a linear term you should use the
        coeff() function instead.
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

    def subs(self, in_dict=None, **kwds):
        """
        EXAMPLES:
            sage: var('x,y,z,a,b,c,d,e,f',ns=1); S = parent(x)
            (x, y, z, a, b, c, d, e, f)
            sage: w0 = S.wild(0); w1 = S.wild(1)
            sage: t = a^2 + b^2 + (x+y)^3

            # substitute with keyword arguments
            sage: t.subs(a=c)
            (x + y)^3 + b^2 + c^2

            sage: t.subs(w0 = w0^2)
            (x^2 + y^2)^18 + a^16 + b^16

            # substitute with a dictionary argument
            sage: t.subs({a^2: c})
            (x + y)^3 + b^2 + c

            sage: t.subs({w0^2: w0^3})
            (x + y)^3 + a^3 + b^3

            # substitute with a relational expression
            sage: t.subs(w0^2 == w0^3)
            (x + y)^3 + a^3 + b^3

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
            TypeError: subs takes either a single keyword argument, or a dictionary, or a symbolic relational expression

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
        """
        cdef dict sdict = {}
        if in_dict is not None:
            if isinstance(in_dict, Expression):
                return self._subs_expr(in_dict)
            if not isinstance(in_dict, dict):
                raise TypeError, "subs takes either a single keyword argument, or a dictionary, or a symbolic relational expression"
            sdict = in_dict

        if kwds:
            from sage.misc.sage_eval import sage_eval
            for k, v in kwds.iteritems():
                k = sage_eval(k, locals=globals())
                sdict[k] = v

        cdef GExMap smap
        for k, v in sdict.iteritems():
            smap.insert(make_pair((<Expression>self.coerce_in(k))._gobj,
                (<Expression>self.coerce_in(v))._gobj))

        return new_Expression_from_GEx(self._gobj.subs_map(smap))

    cpdef Expression _subs_expr(self, expr):
        """
        EXAMPLES:
            sage: var('x,y,z,a,b,c,d,e,f',ns=1); S = parent(x)
            (x, y, z, a, b, c, d, e, f)
            sage: w0 = S.wild(0); w1 = S.wild(1)
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
        return new_Expression_from_GEx(self._gobj.subs(p._gobj))

    def __call__(self, *args, **kwds):
        """
        Calls the .subs() method on this expression.

        EXAMPLES:
            sage: var('x,y,z',ns=1)
            (x, y, z)
            sage: (x+y)(x=z^2, y=x^y)
            x^y + z^2
        """
        return self.subs(*args, **kwds)

    def variables(self):
        """
        Return sorted list of variables that occur in this expression.

        EXAMPLES:
            sage: (x,y,z) = var('x,y,z', ns=1)
            sage: (x+y).variables()
            [x, y]
            sage: (2*x).variables()
            [x]
            sage: (x^y).variables()
            [x, y]
            sage: sin(x+y^z).variables()
            [x, y, z]

        """
        cdef GExSet sym_set
        g_list_symbols(self._gobj, sym_set)
        res = []
        cdef GExSetIter itr = sym_set.begin()
        while itr.is_not_equal(sym_set.end()):
            res.append(new_Expression_from_GEx(itr.obj()))
            itr.inc()
        return res

    def nargs(self):
        """
        Returns the number of arguments of this expression.

        EXAMPLES:
            sage: var('a,b,c,x,y',ns=1); S = parent(x)
            (a, b, c, x, y)
            sage: a.nargs()
            0
            sage: (a^2 + b^2 + (x+y)^2).nargs()
            3
            sage: (a^2).nargs()
            2
            sage: (a*b^2*c).nargs()
            3
        """
        return self._gobj.nops()

    def args(self):
        """
        Returns a list containing the arguments of this expression.

        EXAMPLES:
            sage: var('a,b,c,x,y',ns=1); S = parent(x)
            (a, b, c, x, y)
            sage: (a^2 + b^2 + (x+y)^2).args()
            [(x + y)^2, a^2, b^2]
            sage: (a^2).args()
            [a, 2]
            sage: (a*b^2*c).args()
            [a, b^2, c]
        """
        return [new_Expression_from_GEx(self._gobj.op(i)) \
                            for i from 0 <= i < self._gobj.nops()]

    def operator(self):
        """
        Returns the topmost operator in this expression.

        EXAMPLES:
            sage: x,y,z = var('x,y,z',ns=1)
            sage: (x+y).operator()
            <built-in function add>
            sage: (x^y).operator()
            <built-in function pow>
            sage: (x^y * z).operator()
            <built-in function mul>
            sage: (x < y).operator()
            <built-in function lt>

            sage: abs(x).operator()
            <built-in function abs>
            sage: r = gamma(x).operator(); type(r)
            <class 'sage.calculus.calculus.Function_gamma'>

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
            if serial < get_ginac_serial():
                from sage.symbolic.function import get_sfunction_map
                fn = get_sfunction_map()[serial]
                if fn is None:
                    raise NotImplementedError, "Sage equivalent of this special function is not implemented."
                return fn
            else: # else, look up the serial in the sfunction database
                res = get_sfunction_from_serial(serial)
                if res is None:
                    raise RuntimeError, "cannot find SFunction in table"
                return res

        # self._gobj is either a symbol, constant or numeric
        return None

    def __iter__(self):
        """
        Return an iterator over the arguments of this expression.

        EXAMPLES:
            sage: x,y,z = var('x,y,z',ns=1)
            sage: list(iter(x+y+z))
            [x, y, z]
            sage: list(iter(x*y*z))
            [x, y, z]
            sage: list(iter(x^y*z*(x+y)))
            [x + y, x^y, z]
        """
        return new_ExpIter_from_Expression(self)

    def __getitem__(self, ind):
        """
        EXAMPLES:
            sage: x,y,z = var('x,y,z',ns=1)
            sage: e = x + x*y + z^y + 3*y*z; e
            x*y + 3*y*z + z^y + x
            sage: e[1]
            3*y*z
            sage: e[-1]
            x
            sage: e[1:]
            [3*y*z, z^y, x]
            sage: e[:2]
            [x*y, 3*y*z]
            sage: e[-2:]
            [z^y, x]
            sage: e[:-2]
            [x*y, 3*y*z]

        """
        cdef int bind, eind, step, i
        cdef int n_ops = self._gobj.nops()
        if PY_TYPE_CHECK(ind, slice):
            if ind.start:
                bind = ind.start
                if bind != ind.start:
                    raise ValueError, "integer index expected"
                if bind < 0:
                    bind = n_ops + bind
            else:
                bind = 0
            if ind.stop:
                eind = ind.stop
                if eind != ind.stop:
                    raise ValueError, "integer index expected"
                if eind > n_ops:
                    eind = n_ops
                if eind < 0:
                    eind = n_ops + eind
            else:
                eind = n_ops
            if ind.step:
                step = ind.step
                if step != ind.step:
                    raise ValueError, "step value must be an integer"
            else:
                step = 1
            return [new_Expression_from_GEx(self._gobj.op(i))
                    for i in xrange(bind, eind, step)]

        try:
            bind = ind
            if bind != ind:
                raise ValueError, "integer index expected"
        except TypeError:
            raise TypeError, "index should either be a slice object, or an integer"
        if bind < 0:
            bind = n_ops + bind
        return new_Expression_from_GEx(self._gobj.op(bind))

    def n(self, prec=None, digits=None):
        """
        Return a numerical approximation of self.

        ALIAS: .numerical_approx()

        EXAMPLES::

            sage: var('x', ns=1)
            x
            sage: sin(x).subs(x=5).n()
            -0.958924274663138
            sage: sin(x).subs(x=5).n(100)
            -0.95892427466313846889315440616
            sage: sin(x).subs(x=5).n(digits=50)
            -0.95892427466313846889315440615599397335246154396460
            sage: zeta(x).subs(x=2).numerical_approx(digits=50)
            1.6449340668482264364724151666460251892189499012068

        TESTS::

        We test the evaluation of different infinities available in pynac::
            sage: t = x - oo; t
            -Infinity
            sage: t.n()
            -Infinity
            sage: t = x + oo; t
            +Infinity
            sage: t.n()
            +Infinity
            sage: t = x - unsigned_infinity; t
            Infinity
            sage: t.n()
            Infinity
        """
        if prec is None:
            if digits is None:
                prec = 53
            else:
                prec = int((digits+1) * 3.32192) + 1
        return new_Expression_from_GEx(self._gobj.evalf(0, prec))

    numerical_approx = n

    def function(self, *args):
        """
        Return a callable symbolic expression with the given variables.

        EXAMPLES:
            sage: var('x,y,z',ns=1)
            (x, y, z)
            sage: f = (x+2*y).function(x,y); f
            (x, y) |--> x + 2*y
            sage: f(1,2)
            5

            sage: f(1)
            2*y + 1

        """
        # we override type checking in CallableSymbolicExpressionRing,
        # since it checks for old SymbolicVariable's
        # and do the check here instead
        for i in args:
            if not PY_TYPE_CHECK(i, Expression):
                break
            elif not is_a_symbol((<Expression>i)._gobj):
                break
        else:
            R = CallableSymbolicExpressionRing(args, check=False)
            return R(self)
        raise TypeError, "Must construct a function with a tuple (or list) of symbolic variables."

    ############################################################################
    # Polynomial functions
    ############################################################################
    def coeff(self, s, int n=1):
        """
        INPUT:
            s -- expression
            n -- integer
        OUTPUT:
            coefficient of s^n

        EXAMPLES:
            sage: var('x,y,a', ns=1)
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + x/y + 2*sin(x*y)/x + 100
            sage: f.collect(x)
            x^3*sin(x*y) + (a + y + 1/y)*x + 2*sin(x*y)/x + 100
            sage: f.coeff(x,0)
            100
            sage: f.coeff(x,-1)
            2*sin(x*y)
            sage: f.coeff(x,1)
            a + y + 1/y
            sage: f.coeff(x,2)
            0
            sage: f.coeff(x,3)
            sin(x*y)
            sage: f.coeff(x^3)
            sin(x*y)
            sage: f.coeff(sin(x*y))
            x^3 + 2/x
            sage: f.collect(sin(x*y))
            (x^3 + 2/x)*sin(x*y) + a*x + x*y + x/y + 100
        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._gobj.coeff(ss._gobj, n))

    def leading_coeff(self, s):
        """
        Return the leading coefficient of s in self.

        EXAMPLES:
            sage: var('x,y,a', ns=1)
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + x/y + 2*sin(x*y)/x + 100
            sage: f.leading_coeff(x)
            sin(x*y)
            sage: f.leading_coeff(y)
            x
            sage: f.leading_coeff(sin(x*y))
            x^3 + 2/x
        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._gobj.lcoeff(ss._gobj))

    def trailing_coeff(self, s):
        """
        Return the trailing coefficient of s in self, i.e., the coefficient
        of the smallest power of s in self.

        EXAMPLES:
            sage: var('x,y,a', ns=1)
            (x, y, a)
            sage: f = 100 + a*x + x^3*sin(x*y) + x*y + x/y + 2*sin(x*y)/x; f
            x^3*sin(x*y) + a*x + x*y + x/y + 2*sin(x*y)/x + 100
            sage: f.trailing_coeff(x)
            2*sin(x*y)
            sage: f.trailing_coeff(y)
            x
            sage: f.trailing_coeff(sin(x*y))
            a*x + x*y + x/y + 100
        """
        cdef Expression ss = self.coerce_in(s)
        return new_Expression_from_GEx(self._gobj.tcoeff(ss._gobj))

    def low_degree(self, s):
        """
        Return the exponent of the lowest nonpositive power of s in self.

        OUTPUT:
            an integer <= 0.

        EXAMPLES:
            sage: var('x,y,a', ns=1)
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
           an integer >= 0.

        EXAMPLES:
            sage: var('x,y,a', ns=1)
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

    def gcd(self, b):
        """
        Return the gcd of self and b, which must be integers or polynomials over
        the rational numbers.

        TODO: I tried the massive gcd from
        http://trac.sagemath.org/sage_trac/ticket/694 on Ginac dies
        after about 10 seconds.  Singular easily does that GCD now.
        Since Ginac only handles poly gcd over QQ, we should change
        ginac itself to use Singular.

        EXAMPLES:
            sage: var('x,y',ns=1); S = parent(x)
            (x, y)
            sage: S(10).gcd(S(15))
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
        return new_Expression_from_GEx(x)

    def collect_common_factors(self):
        """

        EXAMPLES:
            sage: var('x', ns=1)
            x
            sage: (x/(x^2 + x)).collect_common_factors()
            1/(x + 1)
        """
        _sig_on
        cdef GEx x = g_collect_common_factors(self._gobj)
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
            sage: SR(I).csgn()
            1
        """
        return new_Expression_from_GEx(g_csgn(self._gobj))

    def conjugate(self):
        """
        Return the complex conjugate of self.

        EXAMPLES::

            sage: x = var('x', ns=1); SR = x.parent()
            sage: SR(CDF.0).conjugate()
            -1.0*I
            sage: x.conjugate()
            conjugate(x)
            sage: SR(RDF(1.5)).conjugate()
            1.5
            sage: SR(float(1.5)).conjugate()
            1.5
            sage: SR(I).conjugate()
            -I
            sage: ( 1+I  + (2-3*I)*x).conjugate()
            (3*I + 2)*conjugate(x) + 1 - I
        """
        return new_Expression_from_GEx(self._gobj.conjugate())

    def real_part(self):
        """
        Return the real part of self.

        EXAMPLES:
            sage: x = var('x', ns=1); SR = x.parent()
            sage: x.real_part()
            real_part(x)
            sage: SR(2+3*I).real_part()
            2
            sage: SR(CDF(2,3)).real_part()
            2.0
            sage: SR(CC(2,3)).real_part()
            2.00000000000000
        """
        return new_Expression_from_GEx(self._gobj.real_part())

    def imag_part(self):
        """
        Return the imaginary part of self.

        EXAMPLES:
            sage: x = var('x', ns=1); SR = x.parent()
            sage: x.imag_part()
            imag_part(x)
            sage: SR(2+3*I).imag_part()
            3
            sage: SR(CC(2,3)).imag_part()
            3.00000000000000
            sage: SR(CDF(2,3)).imag_part()
            3.0
        """
        return new_Expression_from_GEx(self._gobj.imag_part())

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
        EXAMPLES::

            sage: var('x, y', ns=1); S = parent(x)
            (x, y)
            sage: sin(x^2 + y^2)
            sin(x^2 + y^2)
            sage: sin(sage.symbolic.constants.pi)
            0
            sage: sin(S(1))
            sin(1)
            sage: sin(S(RealField(150)(1)))
            sin(1.0000000000000000000000000000000000000000000)

        In order to get a numeric approximation use .n()::

            sage: sin(S(1)).n(150)
            0.84147098480789650665250232163029899962256306

        """
        return new_Expression_from_GEx(g_sin(self._gobj))

    def cos(self):
        """
        Return the cosine of self.

        EXAMPLES::

            sage: var('x, y', ns=1); S = parent(x)
            (x, y)
            sage: cos(x^2 + y^2)
            cos(x^2 + y^2)
            sage: cos(sage.symbolic.constants.pi)
            -1
            sage: cos(S(1))
            cos(1)
            sage: cos(S(RealField(150)(1)))
            cos(1.0000000000000000000000000000000000000000000)

        In order to get a numeric approximation use .n()::

            sage: cos(S(1)).n(150)
            0.54030230586813971740093660744297660373231042
            sage: S(RR(1)).cos().n()
            0.540302305868140
            sage: S(float(1)).cos().n()
            0.540302305868140
        """
        return new_Expression_from_GEx(g_cos(self._gobj))

    def tan(self):
        """
        EXAMPLES::

            sage: var('x, y', ns=1); S = parent(x)
            (x, y)
            sage: tan(x^2 + y^2)
            tan(x^2 + y^2)
            sage: tan(sage.symbolic.constants.pi/2)
            Traceback (most recent call last):
            ...
            ValueError: simple pole at 1/2*pi
            sage: tan(S(1))
            tan(1)
            sage: tan(S(RealField(150)(1)))
            tan(1.0000000000000000000000000000000000000000000)

        Use .n() to get a numerical approximation::

            sage: tan(S(1)).n(150)
            1.5574077246549022305069748074583601730872508
        """
        try:
            return new_Expression_from_GEx(g_tan(self._gobj))
        except RuntimeError:
            raise ValueError, "simple pole at %s"%(self)

    def arcsin(self):
        """
        Return the arcsin of x, i.e., the number y between -pi and pi
        such that sin(y) == x.

        EXAMPLES::

            sage: x = var('x', ns=1); SR = x.parent()
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
        """
        return new_Expression_from_GEx(g_asin(self._gobj))

    def arccos(self):
        """
        Return the arc cosine of self.

        EXAMPLES::

            sage: x = var('x', ns=1); S = x.parent()
            sage: x.arccos()
            arccos(x)
            sage: S(1).arccos()
            0
            sage: S(1/2).arccos()
            1/3*pi
            sage: S(0.4).arccos()
            arccos(0.400000000000000)

        Use .n() to get a numerical approximation::

            sage: S(0.4).arccos().n()
            1.15927948072741
            sage: plot(lambda x: S(x).arccos(), -1,1)
        """
        return new_Expression_from_GEx(g_acos(self._gobj))

    def arctan(self):
        """
        Return the arc tangent of self.

        EXAMPLES::

            sage: x = var('x', ns=1); S = x.parent()
            sage: x.arctan()
            arctan(x)
            sage: S(1).arctan()
            1/4*pi
            sage: S(1/2).arctan()
            arctan(1/2)
            sage: S(0.5).arctan()
            arctan(0.500000000000000)

        Use .n() to get a numerical approximation::

            sage: S(0.5).arctan().n()
            0.463647609000806
            sage: plot(lambda x: S(x).arctan(), -20,20)
        """
        return new_Expression_from_GEx(g_atan(self._gobj))

    def arctan2(self, x):
        """
        Return the inverse of the 2-variable tan function on self and x.

        EXAMPLES::

            sage: var('x,y', ns=1); S = parent(x)
            (x, y)
            sage: x.arctan2(y)
            arctan2(x, y)
            sage: S(1/2).arctan2(1/2)
            1/4*pi
            sage: maxima.eval('atan2(1/2,1/2)')
            '%pi/4'

            sage: S(-0.7).arctan2(S(-0.6))
            -pi + arctan(1.16666666666667)

        Use .n() to get a numerical approximation::

            sage: S(-0.7).arctan2(S(-0.6)).n()
            -2.27942259892257

        TESTS:

        We compare a bunch of different evaluation points between
        Sage and Maxima::

            sage: float(S(0.7).arctan2(0.6))
            0.8621700546672264
            sage: maxima('atan2(0.7,0.6)')
            .8621700546672264
            sage: float(S(0.7).arctan2(-0.6))
            2.2794225989225669
            sage: maxima('atan2(0.7,-0.6)')
            2.279422598922567
            sage: float(S(-0.7).arctan2(0.6))
            -0.8621700546672264
            sage: maxima('atan2(-0.7,0.6)')
            -.8621700546672264
            sage: float(S(-0.7).arctan2(-0.6))
            -2.2794225989225669
            sage: maxima('atan2(-0.7,-0.6)')
            -2.279422598922567
            sage: float(S(0).arctan2(-0.6))
            3.1415926535897931
            sage: maxima('atan2(0,-0.6)')
            3.141592653589793
            sage: float(S(0).arctan2(0.6))
            0.0
            sage: maxima('atan2(0,0.6)')
            0.0
            sage: S(0).arctan2(0)
            0

            sage: S(I).arctan2(1)
            arctan2(I, 1)
            sage: S(CDF(0,1)).arctan2(1)
            arctan2(1.0*I, 1)
            sage: S(1).arctan2(CDF(0,1))
            arctan2(1, 1.0*I)
        """
        cdef Expression nexp = self.coerce_in(x)
        return new_Expression_from_GEx(g_atan2(self._gobj, nexp._gobj))

    def sinh(self):
        r"""
        Return sinh of self.

        We have $\sinh(x) = (e^{x} - e^{-x})/2$.

        EXAMPLES::

            sage: x = var('x', ns=1); S = x.parent()
            sage: x.sinh()
            sinh(x)
            sage: S(1).sinh()
            sinh(1)
            sage: S(0).sinh()
            0
            sage: S(1.0).sinh()
            sinh(1.00000000000000)

        Use .n() to get a numerical approximation::

            sage: S(1.0).sinh().n()
            1.17520119364380
            sage: maxima('sinh(1.0)')
            1.175201193643801

            sage: S(1.0000000000000000000000000).sinh()
            sinh(1.0000000000000000000000000)
            sage: S(1).sinh().n(90)
            1.1752011936438014568823819
            sage: S(RIF(1)).sinh()
            sinh(1)
            sage: S(RIF(1)).sinh().n()
            1.17520119364380
            sage: plot(lambda x: S(x).sinh(), -1, 1)
        """
        return new_Expression_from_GEx(g_sinh(self._gobj))

    def cosh(self):
        """
        Return cosh of self.

        We have $\sinh(x) = (e^{x} + e^{-x})/2$.

        EXAMPLES::

            sage: x = var('x', ns=1); S = x.parent()
            sage: x.cosh()
            cosh(x)
            sage: S(1).cosh()
            cosh(1)
            sage: S(0).cosh()
            1
            sage: S(1.0).cosh()
            cosh(1.00000000000000)

        Use .n() to get a numerical approximation::

            sage: S(1.0).cosh().n()
            1.54308063481524
            sage: maxima('cosh(1.0)')
            1.543080634815244
            sage: S(1.0000000000000000000000000).cosh()
            cosh(1.0000000000000000000000000)
            sage: S(1).cosh().n(90)
            1.5430806348152437784779056
            sage: S(RIF(1)).cosh()
            cosh(1)
            sage: S(RIF(1)).cosh().n()
            1.54308063481524
            sage: plot(lambda x: S(x).cosh(), -1, 1)
        """
        return new_Expression_from_GEx(g_cosh(self._gobj))

    def tanh(self):
        """
        Return tanh of self.

        We have $\tanh(x) = \sinh(x) / \cosh(x)$.

        EXAMPLES::

            sage: x = var('x', ns=1); S = x.parent()
            sage: x.tanh()
            tanh(x)
            sage: S(1).tanh()
            tanh(1)
            sage: S(0).tanh()
            0
            sage: S(1.0).tanh()
            tanh(1.00000000000000)

        Use .n() to get a numerical approximation::

            sage: S(1.0).tanh().n()
            0.761594155955765
            sage: maxima('tanh(1.0)')
            .7615941559557649
            sage: plot(lambda x: S(x).tanh(), -1, 1)
        """
        return new_Expression_from_GEx(g_tanh(self._gobj))

    def arcsinh(self):
        """
        Return the inverse hyperbolic sine of self.

        EXAMPLES::

            sage: x = var('x', ns=1); S = x.parent()
            sage: x.arcsinh()
            arcsinh(x)
            sage: S(0).arcsinh()
            0
            sage: S(1).arcsinh()
            arcsinh(1)
            sage: S(1.0).arcsinh()
            arcsinh(1.00000000000000)

        Use .n() to get a numerical approximation::

            sage: S(1.0).arcsinh().n()
            0.881373587019543
            sage: maxima('asinh(1.0)')
            0.881373587019543

        Sage automatically applies certain identies::
            sage: S(3/2).arcsinh().cosh()
            1/2*sqrt(13)
        """
        return new_Expression_from_GEx(g_asinh(self._gobj))

    def arccosh(self):
        """
        Return the inverse hyperbolic cosine of self.

        EXAMPLES::

            sage: x = var('x', ns=1); S = x.parent()
            sage: x.arccosh()
            arccosh(x)
            sage: S(0).arccosh()
            1/2*I*pi
            sage: S(1/2).arccosh()
            arccosh(1/2)
            sage: S(CDF(1/2)).arccosh()
            arccosh(0.5)

        Use .n() to get a numerical approximation::

            sage: S(CDF(1/2)).arccosh().n()
            1.0471975512*I
            sage: maxima('acosh(0.5)')
            1.047197551196598*%i
        """
        return new_Expression_from_GEx(g_acosh(self._gobj))

    def arctanh(self):
        """
        Return the inverse hyperbolic tangent of self.

        EXAMPLES::

            sage: x = var('x', ns=1); S = x.parent()
            sage: x.arctanh()
            arctanh(x)
            sage: S(0).arctanh()
            0
            sage: S(1/2).arctanh()
            arctanh(1/2)
            sage: S(0.5).arctanh()
            arctanh(0.500000000000000)

        Use .n() to get a numerical approximation::

            sage: S(0.5).arctanh().n()
            0.549306144334055
            sage: S(0.5).arctanh().tanh()
            0.500000000000000
            sage: maxima('atanh(0.5)')
            .5493061443340548
        """
        return new_Expression_from_GEx(g_atanh(self._gobj))

    def exp(self):
        """
        Return exponential function of self, i.e., e to the
        power of self.

        EXAMPLES::

            sage: x = var('x', ns=1); S = x.parent()
            sage: x.exp()
            e^x
            sage: S(0).exp()
            1
            sage: S(1/2).exp()
            e^(1/2)
            sage: S(0.5).exp()
            e^0.500000000000000


        Use .n() to get a numerical approximation::

            sage: S(0.5).exp().n()
            1.64872127070013
            sage: math.exp(0.5)
            1.6487212707001282

            sage: S(0.5).exp().log()
            0.500000000000000
            sage: plot(lambda x: S(x).exp(), -2,1)
        """
        return new_Expression_from_GEx(g_exp(self._gobj))

    def log(self):
        """
        Return the logarithm of self.

        EXAMPLES::

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
            log(0.500000000000000)

        Use .n() to get a numerical approximation::

            sage: S(0.5).log().n()
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
        EXAMPLES::
            sage: x, y = var('x, y', ns=1); S = x.parent()
            sage: (x/y).zeta()
            zeta(x/y)
            sage: S(2).zeta()
            1/6*pi^2
            sage: S(3).zeta()
            zeta(3)
            sage: S(CDF(0,1)).zeta()
            zeta(1.0*I)
            sage: S(CDF(0,1)).zeta().n()
            0.00330022368532 - 0.418155449141*I
            sage: CDF(0,1).zeta()
            0.00330022368532 - 0.418155449141*I
            sage: plot(lambda x: S(x).zeta(), -10,10).show(ymin=-3,ymax=3)

        TESTS::

            sage: t = S(1).zeta(); t
            zeta(1)
            sage: t.n()
            Infinity
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
            factorial(x)
            sage: (x^2+y^3).factorial()
            factorial(x^2 + y^3)
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
            sage: S(10.0r).gamma()
            gamma(10.0)

        Use .n() to get a numerical approximation::

            sage: S(10.0r).gamma().n()
            362880.000000000
            sage: S(CDF(1,1)).gamma()
            gamma(1.0 + 1.0*I)
            sage: S(CDF(1,1)).gamma().n()
            0.498015668118 - 0.154949828302*I

            sage: gp('gamma(1+I)') # 32-bit
            0.4980156681183560427136911175 - 0.1549498283018106851249551305*I

            sage: gp('gamma(1+I)') # 64-bit
            0.49801566811835604271369111746219809195 - 0.15494982830181068512495513048388660520*I

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

cdef class ExpressionIterator:
    cdef Expression _ex
    cdef int _ind
    cdef int _len
    def __iter__(self):
        """
        Return this iterator object itself.

        EXAMPLE:
            sage: x,y,z = var('x,y,z',ns=1)
            sage: i = iter(x+y)
            sage: iter(i) is i
            True
        """
        return self

    def __next__(self):
        """
        Return the next component of the expression.

        EXAMPLE:
            sage: x,y,z = var('x,y,z',ns=1)
            sage: i = iter(x+y)
            sage: i.next()
            x
        """
        cdef GEx ex
        if self._ind == self._len:
            raise StopIteration
        ex = self._ex._gobj.op(self._ind)
        self._ind+=1
        return new_Expression_from_GEx(ex)

cdef inline ExpressionIterator new_ExpIter_from_Expression(Expression ex):
    """
    Construct a new iterator over a symbolic expression.

    EXAMPLES:
        sage: x,y,z = var('x,y,z',ns=1)
        sage: i = iter(x+y) #indirect doctest
    """
    # The const_iterator in GiNaC just keeps an integer index to the current
    # subexpression. We do the same here, to avoid the trouble of having to
    # mess with C++ class constructors/destructors.
    cdef ExpressionIterator m = <ExpressionIterator>PY_NEW(ExpressionIterator)
    m._ex = ex
    m._ind = 0
    m._len  = ex._gobj.nops()
    return m
