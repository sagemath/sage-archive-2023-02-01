###############################################################################
#   SAGE: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"

#################################################################
# Initialize the library
#################################################################

#initialize_ginac()

from sage.rings.integer cimport Integer
from sage.rings.real_mpfr import RealNumber

from sage.symbolic.expression cimport Expression, new_Expression_from_GEx, new_Expression_from_pyobject
from sage.symbolic.expression import is_Expression

from sage.structure.element cimport RingElement, Element
from sage.structure.parent_base import ParentWithBase
from sage.rings.ring cimport CommutativeRing
from sage.categories.morphism cimport Morphism

cdef class SymbolicRing(CommutativeRing):
    """
    Symbolic Ring, parent object for all symbolic expressions.
    """
    def __init__(self):
        """
        Initialize the Symbolic Ring.

        EXAMPLES::

            sage: sage.symbolic.ring.SymbolicRing()
            Symbolic Ring
        """
        ParentWithBase.__init__(self, self)
        self._populate_coercion_lists_(convert_method_name='_symbolic_')

    def __reduce__(self):
        """
        EXAMPLES::

           sage: loads(dumps(SR)) == SR           # indirect doctest
           True
        """
        return the_SymbolicRing, tuple([])

    def __hash__(self):
        """
        EXAMPLES::

            sage: hash(SR)   #random
            139682705593888
        """
        return hash(SymbolicRing)

    def _repr_(self):
        """
        Return a string representation of self.

        EXAMPLES::

            sage: repr(SR)
            'Symbolic Ring'
        """
        return "Symbolic Ring"

    def _latex_(self):
        """
        Return latex representation of the symbolic ring.

        EXAMPLES::

            sage: latex(SR)
            \text{SR}
            sage: M = MatrixSpace(SR, 2); latex(M)
            \mathrm{Mat}_{2\times 2}(\text{SR})
        """
        return r'\text{SR}'

    cpdef _coerce_map_from_(self, R):
        """
        EXAMPLES::

            sage: SR.coerce(int(2))
            2
            sage: SR.coerce(pari(2/5))
            2/5
            sage: SR.coerce(-infinity)
            -Infinity
            sage: SR.has_coerce_map_from(ZZ['t'])
            True
            sage: SR.has_coerce_map_from(ZZ['t,u,v'])
            True
            sage: SR.has_coerce_map_from(Frac(ZZ['t,u,v']))
            True
            sage: SR.has_coerce_map_from(GF(5)['t'])
            True
            sage: SR.has_coerce_map_from(SR['t'])
            False
            sage: SR.has_coerce_map_from(Integers(8))
            True
            sage: SR.has_coerce_map_from(GF(9, 'a'))
            True
        """
        if isinstance(R, type):
            if R in [int, float, long, complex]:
                return True

            import numpy
            if R in [numpy.float, numpy.float32, numpy.float64, numpy.float128,
                     numpy.complex, numpy.complex64, numpy.complex128,
                     numpy.complex256]:
                return NumpyToSRMorphism(R, self)

            from sympy.core.basic import Basic
            if issubclass(R, Basic):
                return UnderscoreSageMorphism(R, self)
        else:
            from sage.rings.real_mpfr import mpfr_prec_min

            from sage.rings.all import ( ComplexField,
               is_PolynomialRing, is_MPolynomialRing,
               is_FractionField, RLF, CLF, AA, QQbar, InfinityRing,
               is_RealIntervalField, is_ComplexIntervalField,
               is_IntegerModRing, is_FiniteField)

            from sage.interfaces.maxima import Maxima
            from sage.libs.pari.gen import PariInstance

            if ComplexField(mpfr_prec_min()).has_coerce_map_from(R):
                # Anything with a coercion into any precicion of CC
                return R not in (RLF, CLF, AA, QQbar)
            elif is_PolynomialRing(R) or is_MPolynomialRing(R) or is_FractionField(R):
                base = R.base_ring()
                return base is not self and self.has_coerce_map_from(base)
            elif (R is InfinityRing
                  or is_RealIntervalField(R) or is_ComplexIntervalField(R)
                  or is_IntegerModRing(R) or is_FiniteField(R)):
                return True
            elif isinstance(R, (Maxima, PariInstance)):
                return True

    def _element_constructor_(self, x):
        """
        Coerce `x` into the symbolic expression ring SR.

        EXAMPLES::

            sage: a = SR(-3/4); a
            -3/4
            sage: type(a)
            <type 'sage.symbolic.expression.Expression'>
            sage: a.parent()
            Symbolic Ring
            sage: K.<a> = QuadraticField(-3)
            sage: a + sin(x)
            I*sqrt(3) + sin(x)
            sage: x=var('x'); y0,y1=PolynomialRing(ZZ,2,'y').gens()
            sage: x+y0/y1
            x + y0/y1
            sage: x.subs(x=y0/y1)
            y0/y1
            sage: x + long(1)
            x + 1

        If `a` is already in the symbolic expression ring, coercing returns
        `a` itself (not a copy)::

            sage: a = SR(-3/4); a
            -3/4
            sage: SR(a) is a
            True

        A Python complex number::

            sage: SR(complex(2,-3))
            2.00000000000000 - 3.00000000000000*I

        TESTS::

            sage: SR._coerce_(int(5))
            5
            sage: SR._coerce_(5)
            5
            sage: SR._coerce_(float(5))
            5.00000000000000
            sage: SR._coerce_(5.0)
            5.00000000000000

        An interval arithmetic number::

            sage: SR._coerce_(RIF(pi))
            3.141592653589794?

        A number modulo 7::

            sage: a = SR(Mod(3,7)); a
            3
            sage: a^2
            2

            sage: si = SR.coerce(I)
            sage: si^2
            -1
            sage: bool(si == CC.0)
            True

        """
        cdef GEx exp

        if is_Expression(x):
            if (<Expression>x)._parent is self:
                return x
            else:
                return new_Expression_from_GEx(self, (<Expression>x)._gobj)
        elif hasattr(x, '_symbolic_'):
            return x._symbolic_(self)
        elif isinstance(x, str):
            try:
                from sage.calculus.calculus import symbolic_expression_from_string
                return self(symbolic_expression_from_string(x))
            except SyntaxError, err:
                msg, s, pos = err.args
                raise TypeError, "%s: %s !!! %s" % (msg, s[:pos], s[pos:])

        from sage.interfaces.maxima import MaximaElement
        if isinstance(x, MaximaElement):
            from sage.calculus.calculus import symbolic_expression_from_maxima_element
            return self(symbolic_expression_from_maxima_element(x))

        from sage.rings.infinity import (infinity, minus_infinity,
                                         unsigned_infinity)

        if isinstance(x, (Integer, RealNumber)):
            GEx_construct_pyobject(exp, x)
        elif isinstance(x, int):
            GEx_construct_long(&exp, x)
        elif isinstance(x, float):
            from sage.rings.all import RR
            GEx_construct_pyobject(exp, RR(x))
        elif isinstance(x, long):
            GEx_construct_pyobject(exp, Integer(x))
        elif isinstance(x, complex):
            from sage.rings.all import CC
            GEx_construct_pyobject(exp, CC(x))
        elif x is infinity:
            return new_Expression_from_GEx(self, g_Infinity)
        elif x is minus_infinity:
            return new_Expression_from_GEx(self, g_mInfinity)
        elif x is unsigned_infinity:
            return new_Expression_from_GEx(self, g_UnsignedInfinity)
        elif isinstance(x, RingElement):
            GEx_construct_pyobject(exp, x)
        else:
            raise TypeError

        return new_Expression_from_GEx(self, exp)

    def wild(self, unsigned int n=0):
        """
        Return the n-th wild-card for pattern matching and substitution.

        INPUT:

        - ``n`` - a nonnegative integer

        OUTPUT:

        - `i^{th}` wildcard expression

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: w0 = SR.wild(0); w1 = SR.wild(1)
            sage: pattern = sin(x)*w0*w1^2; pattern
            $0*$1^2*sin(x)
            sage: f = atan(sin(x)*3*x^2); f
            arctan(3*x^2*sin(x))
            sage: f.has(pattern)
            True
            sage: f.subs(pattern == x^2)
            arctan(x^2)
        """
        return new_Expression_from_GEx(self, g_wild(n))

    def __cmp__(self, other):
        """
        Compare two symbolic expression rings. They are equal if and only
        if they have the same type. Otherwise their types are compared.

        EXAMPLES::

            sage: from sage.symbolic.ring import SymbolicRing
            sage: cmp(SR, RR) #random
            1
            sage: cmp(RR, SymbolicRing()) #random
            -1
            sage: cmp(SR, SymbolicRing()) #random
            0
        """
        return cmp(type(self), type(other))

    def __contains__(self, x):
        r"""
        True if there is an element of the symbolic ring that is equal to x
        under ``==``.

        EXAMPLES:

        The symbolic variable x is in the symbolic ring.::

            sage: x.parent()
            Symbolic Ring
            sage: x in SR
            True

        2 is also in the symbolic ring since it is equal to something in
        SR, even though 2's parent is not SR.

        ::

            sage: 2 in SR
            True
            sage: parent(2)
            Integer Ring
            sage: 1/3 in SR
            True
        """
        try:
            x2 = self(x)
            return bool(x2 == x)
        except TypeError:
            return False

    def characteristic(self):
        """
        Return the characteristic of the symbolic ring, which is 0.

        OUTPUT:

        - a Sage integer

        EXAMPLES::

            sage: c = SR.characteristic(); c
            0
            sage: type(c)
            <type 'sage.rings.integer.Integer'>
        """
        return Integer(0)

    def _an_element_(self):
        """
        Return an element of the symbolic ring, which is used by the
        coercion model.

        EXAMPLES::

            sage: SR._an_element_()
            some_variable
        """
        return self('some_variable')

    def is_field(self):
        """
        Returns True, since the symbolic expression ring is (for the most
        part) a field.

        EXAMPLES::

            sage: SR.is_field()
            True
        """
        return True

    def is_exact(self):
        """
        Return False, because there are approximate elements in the
        symbolic ring.

        EXAMPLES::

            sage: SR.is_exact()
            False

        Here is an inexact element.

        ::

            sage: SR(1.9393)
            1.93930000000000
        """
        return False

    def pi(self):
        """
        EXAMPLES::

            sage: SR.pi() is pi
            True
        """
        from sage.symbolic.constants import pi
        return self(pi)

    cpdef symbol(self, name):
        """
        EXAMPLES::

            sage: SR.symbol("asdfasdfasdf")
            asdfasdfasdf
        """
        cdef GSymbol symb = get_symbol(name)
        cdef Expression e
        e = <Expression>PY_NEW(Expression)
        GEx_construct_symbol(&e._gobj, symb)
        e._parent = SR
        return e

    def var(self, name):
        """
        Return the symbolic variable defined by x as an element of the
        symbolic ring.

        EXAMPLES::

            sage: zz = SR.var('zz'); zz
            zz
            sage: type(zz)
            <type 'sage.symbolic.expression.Expression'>
            sage: t = SR.var('theta2'); t
            theta2
        """
        if is_Expression(name):
            return name
        if not isinstance(name, str):
            name = repr(name)
        if ',' in name:
            return tuple([self.symbol(s.strip()) for s in name.split(',')])
        elif ' ' in name:
            return tuple([self.symbol(s.strip()) for s in name.split(' ')])
        else:
            return self.symbol(name)

    def _repr_element_(self, Expression x):
        """
        Returns the string representation of the element x.  This is
        used so that subclasses of the SymbolicRing (such the a
        CallableSymbolicExpressionRing) can provide their own
        implementations of how to print Expressions.

        EXAMPLES::

            sage: SR._repr_element_(x+2)
            'x + 2'
        """
        return GEx_to_str(&x._gobj)

    def _latex_element_(self, Expression x):
        """
        Returns the standard LaTeX version of the expression *x*.

        EXAMPLES::

            sage: latex(sin(x+2))
            \sin\left(x + 2\right)
            sage: latex(var('theta') + 2)
            \theta + 2
        """
        return GEx_to_str_latex(&x._gobj)

    def _call_element_(self, _the_element, *args, **kwds):
        """
        EXAMPLES::

            sage: x,y=var('x,y')
            sage: f = x+y
            sage: f.variables()
            (x, y)
            sage: f()
            x + y
            sage: f(3)
            doctest:...: DeprecationWarning: Substitution using function-call syntax and unnamed arguments is deprecated and will be removed from a future release of Sage; you can use named arguments instead, like EXPR(x=..., y=...)
            y + 3
            sage: f(x=3)
            y + 3
            sage: f(3,4)
            7
            sage: f(x=3,y=4)
            7
            sage: f(2,3,4)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 2
            sage: f(x=2,y=3,z=4)
            5

        ::

            sage: f({x:3})
            y + 3
            sage: f({x:3,y:4})
            7
            sage: f(x=3)
            y + 3
            sage: f(x=3,y=4)
            7

        ::

            sage: a = (2^(8/9))
            sage: a(4)
            Traceback (most recent call last):
            ...
            ValueError: the number of arguments must be less than or equal to 0


        Note that you make get unexpected results when calling
        symbolic expressions and not explicitly giving the variables::

            sage: f = function('Gamma', var('z'), var('w')); f
            Gamma(z, w)
            sage: f(2)
            Gamma(z, 2)
            sage: f(2,5)
            Gamma(5, 2)

        Thus, it is better to be explicit::

            sage: f(z=2)
            Gamma(2, w)
        """
        if len(args) == 0:
            d = None
        elif len(args) == 1 and isinstance(args[0], dict):
            d = args[0]
        else:
            from sage.misc.misc import deprecation
            vars = _the_element.operands()
            deprecation("Substitution using function-call syntax and unnamed arguments is deprecated and will be removed from a future release of Sage; you can use named arguments instead, like EXPR(x=..., y=...)")
            d = {}

            vars = _the_element.variables()
            for i, arg in enumerate(args):
                try:
                    d[ vars[i] ] = arg
                except IndexError:
                    raise ValueError, "the number of arguments must be less than or equal to %s"%len(vars)

        return _the_element.subs(d, **kwds)

SR = SymbolicRing()

cdef class NumpyToSRMorphism(Morphism):
    def __init__(self, numpy_type, R):
        """
        A Morphism which constructs Expressions from NumPy floats and
        complexes by converting them to elements of either RDF or CDF.

        EXAMPLES::

            sage: import numpy
            sage: from sage.symbolic.ring import NumpyToSRMorphism
            sage: f = NumpyToSRMorphism(numpy.float64, SR)
            sage: f(numpy.float64('2.0'))
            2.0
            sage: _.parent()
            Symbolic Ring
        """
        import sage.categories.homset
        from sage.structure.parent import Set_PythonType
        Morphism.__init__(self, sage.categories.homset.Hom(Set_PythonType(numpy_type), R))

    cpdef Element _call_(self, a):
        """
        EXAMPLES:

        This should be called when coercing or converting a NumPy
        float or complex to the Symbolic Ring::

            sage: import numpy
            sage: SR(numpy.float64('2.0')).pyobject().parent()
            Real Double Field
        """
        from sage.rings.all import RDF, CDF
        numpy_type = self._domain.object()
        if 'complex' in numpy_type.__name__:
            res = CDF(a)
        else:
            res = RDF(a)

        return new_Expression_from_pyobject(self._codomain, res)

cdef class UnderscoreSageMorphism(Morphism):
    def __init__(self, t, R):
        """
        A Morphism which constructs Expressions from an arbitrary Python
        object by calling the :meth:`_sage_` method on the object.

        EXAMPLES::

            sage: import sympy
            sage: from sage.symbolic.ring import UnderscoreSageMorphism
            sage: b = sympy.var('b')
            sage: f = UnderscoreSageMorphism(type(b), SR)
            sage: f(b)
            b
            sage: _.parent()
            Symbolic Ring
        """
        import sage.categories.homset
        from sage.structure.parent import Set_PythonType
        Morphism.__init__(self, sage.categories.homset.Hom(Set_PythonType(t), R))

    cpdef Element _call_(self, a):
        """
        EXAMPLES:

        This should be called when coercing or converting a SymPy
        object to the Symbolic Ring::

            sage: import sympy
            sage: b = sympy.var('b')
            sage: bool(SR(b) == SR(b._sage_()))
            True
        """
        return self._codomain(a._sage_())


def the_SymbolicRing():
    """
    Return the unique symbolic ring object.

    (This is mainly used for unpickling.)

    EXAMPES::

        sage: sage.symbolic.ring.the_SymbolicRing()
        Symbolic Ring
        sage: sage.symbolic.ring.the_SymbolicRing() is sage.symbolic.ring.the_SymbolicRing()
        True
        sage: sage.symbolic.ring.the_SymbolicRing() is SR
        True
    """
    return SR

def is_SymbolicExpressionRing(R):
    """
    Returns True if *R* is the symbolic expression ring.

    EXAMPLES::

        sage: from sage.symbolic.ring import is_SymbolicExpressionRing
        sage: is_SymbolicExpressionRing(ZZ)
        False
        sage: is_SymbolicExpressionRing(SR)
        True
    """
    return R is SR

def var(name):
    """
    EXAMPLES::

        sage: from sage.symbolic.ring import var
        sage: var("x y z")
        (x, y, z)
        sage: var("x,y,z")
        (x, y, z)
        sage: var("x , y , z")
        (x, y, z)
        sage: var("z")
        z
    """
    return SR.var(name)

def is_SymbolicVariable(x):
    """
    Returns True if x is a variable.

    EXAMPLES::

        sage: from sage.symbolic.ring import is_SymbolicVariable
        sage: is_SymbolicVariable(x)
        True
        sage: is_SymbolicVariable(x+2)
        False
    """
    return is_Expression(x) and is_a_symbol((<Expression>x)._gobj)
