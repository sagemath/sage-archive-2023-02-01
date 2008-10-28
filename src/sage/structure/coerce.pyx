r"""
The coercion model manages how elements of one parent get relate to elements
of another. For example, the integer 2 can canonically be viewed as an element
of the rational numbers. (The Parent of a non-element is its Python type.)

    sage: ZZ(2).parent()
    Integer Ring
    sage: QQ(2).parent()
    Rational Field

The most prominent role of the coercion model is to make sense of binary
operations between elements that have distinct parents. It does this by
finding a parent where both elements make sense, and doing the operation
there. For example

    sage: a = 1/2; a.parent()
    Rational Field
    sage: b = ZZ[x].gen(); b.parent()
    Univariate Polynomial Ring in x over Integer Ring
    sage: a+b
    x + 1/2
    sage: (a+b).parent()
    Univariate Polynomial Ring in x over Rational Field

If there is a coercion (see below) from one of the parents to the other,
the operation is always performed in the codomain of that coercion. Otherwise
a reasonable attempt to create a new parent with coercion maps from both
original parents is made. The results of these discoveries are cached.
On failure, a TypeError is always raised.

Some arithmetic operations (such as multiplication) can indicate an action
rather than arithmitic in a common parent. For example

    sage: E = EllipticCurve('37a')
    sage: P = E(0,0)
    sage: 5*P
    (1/4 : -5/8 : 1)

where there is action of $\Z$ on the points of $E$ given by the additive
group law. Parents can specify how they act on or are acted upon by other
parents.

There are two kinds of ways to get from one parent to another, coercions and
conversions.

Coercions are cannonical (possibly modulo a finite number of determanistic
choices) morphisms, and the set of all coercions between all parents forms a
comuting diagram (modulo possibly rounding issues). $\Z \rightarrow \Q$ is an
example of a coercion. These are invoked implictly by the coercion model.

Conversions try to construct an element out of their input if at all possible.
Examples include sections of coercions, creating an element from a string or
list, etc. and may fail on some inputs of a given type while succeeding on
others (i.e. they may not be defined on the whole domain). Conversions are
always explicitly invoked, and never used by the coercion model to resolve
binary operations.

For more information on how to specify coercions, conversions, and actions,
see the documentation for Parent.
"""


#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include "../ext/stdsage.pxi"
include "../ext/python_object.pxi"
include "coerce.pxi"

import operator

from sage_object cimport SageObject
from sage.categories.map cimport Map
import sage.categories.morphism
from sage.categories.morphism import IdentityMorphism
from sage.categories.action import InverseAction, PrecomposedAction
from parent import Set_PythonType

import sys

from coerce_actions import LeftModuleAction, RightModuleAction, IntegerMulAction

cpdef py_scalar_parent(py_type):
    """
    Returns the Sage equivalent of the given python type, if one exists.
    If there is no equivalent, return None.

    EXAMPLES:
        sage: from sage.structure.coerce import py_scalar_parent
        sage: py_scalar_parent(int)
        Integer Ring
        sage: py_scalar_parent(long)
        Integer Ring
        sage: py_scalar_parent(float)
        Real Double Field
        sage: py_scalar_parent(complex)
        Complex Double Field
        sage: py_scalar_parent(dict),
        (None,)
    """
    if py_type is int or py_type is long:
        import sage.rings.integer_ring
        return sage.rings.integer_ring.ZZ
    elif py_type is float:
        import sage.rings.real_double
        return sage.rings.real_double.RDF
    elif py_type is complex:
        import sage.rings.complex_double
        return sage.rings.complex_double.CDF
    else:
        return None

cdef object _Integer
cdef bint is_Integer(x):
    global _Integer
    if _Integer is None:
        from sage.rings.integer import Integer as _Integer
    return PY_TYPE_CHECK_EXACT(x, _Integer) or PY_TYPE_CHECK_EXACT(x, int)

cdef class CoercionModel_cache_maps(CoercionModel):
    """
    See also sage.categories.pushout

    EXAMPLES:
        sage: f = ZZ['t','x'].0 + QQ['x'].0 + CyclotomicField(13).gen(); f
        t + x + (zeta13)
        sage: f.parent()
        Multivariate Polynomial Ring in t, x over Cyclotomic Field of order 13 and degree 12
        sage: ZZ['x','y'].0 + ~Frac(QQ['y']).0
        (x*y + 1)/y
        sage: MatrixSpace(ZZ['x'], 2, 2)(2) + ~Frac(QQ['x']).0
        [(2*x + 1)/x           0]
        [          0 (2*x + 1)/x]
        sage: f = ZZ['x,y,z'].0 + QQ['w,x,z,a'].0; f
        w + x
        sage: f.parent()
        Multivariate Polynomial Ring in w, x, y, z, a over Rational Field
        sage: ZZ['x,y,z'].0 + ZZ['w,x,z,a'].1
        2*x

    AUTHOR:
        -- Robert Bradshaw
    """

    def __init__(self, lookup_dict_size=127, lookup_dict_threshold=.75):
        """
        INPUT:
            lookup_dict_size -- initial size of the coercion hashtables
            lookup_dict_threshold -- maximal density of the coercion hashtables before forcing a re-hash

        EXAMPLES:
            sage: from sage.structure.coerce import CoercionModel_cache_maps
            sage: cm = CoercionModel_cache_maps(4, .95)
            sage: A = cm.get_action(ZZ, NumberField(x^2-2, 'a'), operator.mul)
            sage: f, g = cm.coercion_maps(QQ, int)
            sage: f, g = cm.coercion_maps(ZZ, int)
            sage: cm.get_stats()
            ((0, 1.0, 4), (0, 0.25, 1))

        NOTE: In practice 4 would be a really bad number to choose, but it makes the
              hashing deterministic.
        """
        self.reset_cache(lookup_dict_size, lookup_dict_threshold)

    def reset_cache(self, lookup_dict_size=127, lookup_dict_threshold=.75):
        """
        Clear the coercion cache.

        This should have no impact on the result of arithmetic operations, as
        the exact same coercions and actions will be re-discovered when needed.

        It may be useful for debugging, and may also free some memory.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.get_stats()                # random
            ((0, 0.3307086614173229, 3), (0, 0.1496062992125984, 2))
            sage: cm.reset_cache()
            sage: cm.get_stats()
            ((0, 0.0, 0), (0, 0.0, 0))
        """
        # This MUST be a mapping of tuples, where each
        # tuple contains at least two elements that are either
        # None or of type Map.
        self._coercion_maps = TripleDict(lookup_dict_size, threshold=lookup_dict_threshold)
        # This MUST be a mapping to actions.
        self._action_maps = TripleDict(lookup_dict_size, threshold=lookup_dict_threshold)
        # This is a mapping from Parents to Parents, storing the result of division in the given parent.
        self._division_parents = TripleDict(lookup_dict_size, threshold=lookup_dict_threshold)

    def get_cache(self):
        """
        This returns the current cache of coercion maps and actions, primarily
        useful for debugging and introspection.

        EXAMPLES:
            sage: 1 + 1/2
            3/2
            sage: cm = sage.structure.element.get_coercion_model()
            sage: maps, actions = cm.get_cache()

        Now lets see what happens when we do a binary operations with
        an integer and a rational:
            sage: left_morphism, right_morphism = maps[ZZ, QQ]
            sage: print left_morphism
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
            sage: print right_morphism
            None

        We can see that it coerces the left operand from an integer to a
        rational, and doesn't do anything to the right.

        Now for some actions:

            sage: R.<x> = ZZ['x']
            sage: 1/2 * x
            1/2*x
            sage: maps, actions = cm.get_cache()
            sage: act = actions[QQ, R, operator.mul]; act
            Left scalar multiplication by Rational Field on Univariate Polynomial Ring in x over Integer Ring
            sage: act.actor()
            Rational Field
            sage: act.domain()
            Univariate Polynomial Ring in x over Integer Ring
            sage: act.codomain()
            Univariate Polynomial Ring in x over Rational Field
            sage: act(1/5, x+10)
            1/5*x + 2
        """
        return dict([((S, R), mors) for (S, R, op), mors in self._coercion_maps.iteritems()]), \
               dict(self._action_maps.iteritems())

    def get_stats(self):
        """
        This returns the state of the cache of coercion maps and actions,
        primarily useful for debugging and introspection.
    If a class is not part of the coercion system, we should call the
    __rmul__ method when it makes sense.

        The coercion maps are stored in a specialized TripleDict hashtable,
        and the stats returned are (min, avg, max) of the number of items
        per bucket. The lower the better.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.get_stats()                # random
            ((0, 0.16058394160583941, 2), (0, 0.13138686131386862, 3))
        """
        return self._coercion_maps.stats(), self._action_maps.stats()

    def _record_exception(self):
        r"""
        Pushes the last exception that occured onto the stack for later reference,
        for internal use.

        If the stack has not yet been flagged as cleared, we clear it now (rather
        than wasting time to do so for successful operations).

        TEST:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: 1+1/2+2 # make sure there aren't any errors hanging around
            7/2
            sage: cm.exception_stack()
            []
            sage: cm._test_exception_stack()
            sage: cm.exception_stack()
            [(<type 'exceptions.TypeError'>,  TypeError('just a test',),  <traceback object at ...>)]

            The function _test_exception_stack is executing the following code:
            try:
                raise TypeError, "just a test"
            except:
                cm._record_exception()
        """
        if not self._exceptions_cleared:
            self._exception_stack = []
            self._exceptions_cleared = True
        self._exception_stack.append(sys.exc_info())

    def _test_exception_stack(self):
        """
        A function to test the exception stack.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: 1 + 1/11 # make sure there aren't any errors hanging around
            12/11
            sage: cm.exception_stack()
            []
            sage: cm._test_exception_stack()
            sage: cm.exception_stack()
            [(<type 'exceptions.TypeError'>,  TypeError('just a test',),  <traceback object at ...>)]
        """
        try:
            raise TypeError, "just a test"
        except:
            self._record_exception()

    def exception_stack(self):
        r"""
        Returns the list of exceptions that were caught in the course of
        executing the last binary operation. Useful for diagnosis when
        user-defined maps or actions raise exceptions that are caught in
        the course of coercion detection.

        If all went well, this should be the empty list. If things aren't
        happening as you expect, this is a good place to check. See also
        \code{coercion_traceback}.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: 1/2 + 2
            5/2
            sage: cm.exception_stack()
            []
            sage: 1/2 + GF(3)(2)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': 'Rational Field' and 'Finite Field of size 3'

        Now see what the actual problem was:
            sage: import traceback
            sage: cm.exception_stack()
            [(<type 'exceptions.TypeError'>, TypeError("BUG: the base_extend method must be defined for 'Monoid of ideals of Integer Ring' (class '<class 'sage.rings.ideal_monoid.IdealMonoid_c'>')",), <traceback object at ...>), (<type 'exceptions.TypeError'>,  TypeError("no common canonical parent for objects with parents: 'Rational Field' and 'Finite Field of size 3'",),  <traceback object at ...>)]
            sage: print ''.join(sum([traceback.format_exception(*info) for info in cm.exception_stack()], []))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Rational Field' and 'Finite Field of size 3'

        This is typically accessed via the \code{coercion_traceback} function.
            sage: coercion_traceback()
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Rational Field' and 'Finite Field of size 3'
        """
        if not self._exceptions_cleared:
            self._exception_stack = []
            self._exceptions_cleared = True
        return self._exception_stack


    def explain(self, xp, yp, op=operator.mul, int verbosity=2):
        """
        This function can be used to understand what coercions will happen
        for an arithmetic operation between xp and yp (which may be either
        elements or parents). If the parent of the result can be determined
        then it will be returned.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()

            sage: cm.explain(ZZ, ZZ)
            Identical parents, arithmetic performed immediately.
            Result lives in Integer Ring
            Integer Ring

            sage: cm.explain(QQ, int)
            Coercion on right operand via
                Native morphism:
                  From: Set of Python objects of type 'int'
                  To:   Rational Field
            Arithmetic performed after coercions.
            Result lives in Rational Field
            Rational Field

            sage: cm.explain(ZZ['x'], QQ)
            Action discovered.
                Right scalar multiplication by Rational Field on Univariate Polynomial Ring in x over Integer Ring
            Result lives in Univariate Polynomial Ring in x over Rational Field
            Univariate Polynomial Ring in x over Rational Field

            sage: cm.explain(ZZ['x'], QQ, operator.add)
            Coercion on left operand via
                Call morphism:
                  From: Univariate Polynomial Ring in x over Integer Ring
                  To:   Univariate Polynomial Ring in x over Rational Field
            Coercion on right operand via
                Call morphism:
                  From: Rational Field
                  To:   Univariate Polynomial Ring in x over Rational Field
            Arithmetic performed after coercions.
            Result lives in Univariate Polynomial Ring in x over Rational Field
            Univariate Polynomial Ring in x over Rational Field

        Sometimes with non-sage types there is not enough information to deduce
        what will actually happen:

            sage: cm.explain(RealField(100), float, operator.add)
            Right operand is numeric, will attempt conversion in both directions.
            Unknown result parent.
            sage: parent(RealField(100)(1) + float(1))
            Real Field with 100 bits of precision
            sage: cm.explain(QQ, float, operator.add)
            Right operand is numeric, will attempt conversion in both directions.
            Unknown result parent.
            sage: parent(QQ(1) + float(1))
            <type 'float'>


        Special care is taken to deal with division:

            sage: cm.explain(ZZ, ZZ, operator.div)
            Identical parents, arithmetic performed immediately.
            Result lives in Rational Field
            Rational Field

            sage: cm.explain(ZZ['x'], QQ['x'], operator.div)
            Coercion on left operand via
                Call morphism:
                  From: Univariate Polynomial Ring in x over Integer Ring
                  To:   Univariate Polynomial Ring in x over Rational Field
            Arithmetic performed after coercions.
            Result lives in Fraction Field of Univariate Polynomial Ring in x over Rational Field
            Fraction Field of Univariate Polynomial Ring in x over Rational Field

            sage: cm.explain(int, ZZ, operator.div)
            Coercion on left operand via
                Native morphism:
                  From: Set of Python objects of type 'int'
                  To:   Integer Ring
            Arithmetic performed after coercions.
            Result lives in Rational Field
            Rational Field

            sage: cm.explain(ZZ['x'], ZZ, operator.div)
            Action discovered.
                Right inverse action by Rational Field on Univariate Polynomial Ring in x over Integer Ring
                with precomposition on right by Natural morphism:
                  From: Integer Ring
                  To:   Rational Field
            Result lives in Univariate Polynomial Ring in x over Rational Field
            Univariate Polynomial Ring in x over Rational Field

        NOTE: This function is accurate only in so far as analyse is kept in
              sync with the \code{bin_op} and \code{canonical_coercion} which
              are kept seperate for maximal efficiency.
        """
        all, res = self.analyse(xp, yp, op)
        indent = " "*4
        if verbosity >= 2:
            print "\n".join([s if isinstance(s, str) else indent+(repr(s).replace("\n", "\n"+indent)) for s in all])
        elif verbosity >= 1:
            print "\n".join([s for s in all if isinstance(s, str)])
        if verbosity >= 1:
            if res is None:
                print "Unknown result parent."
            else:
                print "Result lives in", res
        return res

    cpdef analyse(self, xp, yp, op=mul):
        """
        Emulate the process of doing arithmetic between xp and yp, returning
        a list of steps and the parent that the result will live in. The
        \code{explain} function is easier to use, but if one wants access to
        the acutal morphism and action objects (rather than their string
        representations) then this is the function to use.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: steps, res = cm.analyse(GF(7), ZZ)
            sage: print steps
            ['Coercion on right operand via', Natural morphism:
              From: Integer Ring
              To:   Finite Field of size 7, 'Arithmetic performed after coercions.']
            sage: print res
            Finite Field of size 7
            sage: f = steps[1]; type(f)
            <type 'sage.rings.integer_mod.Integer_to_IntegerMod'>
            sage: f(100)
            2
            """
        self._exceptions_cleared = False
        if not PY_TYPE_CHECK(xp, type) and not PY_TYPE_CHECK(xp, Parent):
            xp = parent_c(xp)
        if not PY_TYPE_CHECK(yp, type) and not PY_TYPE_CHECK(yp, Parent):
            yp = parent_c(yp)

        all = []
        if xp is yp:
            all.append("Identical parents, arithmetic performed immediately." % xp)
            if op is div and PY_TYPE_CHECK(xp, Parent):
                xp = self.division_parent(xp)
            return all, xp
        if xp == yp:
            all.append("Equal but distinct parents.")

        if (op is not sub) and (op is not isub):
            action = self.get_action(xp, yp, op)
            if action is not None:
                all.append("Action discovered.")
                all.append(action)
                return all, action.codomain()

        homs = self.discover_coercion(xp, yp)
        if homs is not None:
            x_mor, y_mor = homs
            if x_mor is not None:
                all.append("Coercion on left operand via")
                all.append(x_mor)
                res = x_mor.codomain()
            if y_mor is not None:
                all.append("Coercion on right operand via")
                all.append(y_mor)
                if res is not None and res is not y_mor.codomain():
                    raise RuntimeError, ("BUG in coercion model: codomains not equal!", x_mor, y_mor)
                res = y_mor.codomain()
            all.append("Arithmetic performed after coercions.")
            if op is div and PY_TYPE_CHECK(res, Parent):
                res = self.division_parent(res)
            return all, res

        if PY_TYPE_CHECK(yp, Parent) and xp in [int, long, float, complex, bool]:
            mor = yp.coerce_map_from(xp)
            if mor is not None:
                all.append("Coercion on numeric left operand via")
                all.append(mor)
                if op is div and PY_TYPE_CHECK(yp, Parent):
                    yp = self.division_parent(yp)
                return all, yp
            all.append("Left operand is numeric, will attempt conversion in both directions.")
        elif type(xp) is type:
            all.append("Left operand is not Sage element, will try _sage_.")

        if PY_TYPE_CHECK(xp, Parent) and yp in [int, long, float, complex, bool]:
            mor = xp.coerce_map_from(yp)
            if mor is not None:
                all.append("Coercion on numeric right operand via")
                all.append(mor)
                if op is div and PY_TYPE_CHECK(xp, Parent):
                    xp = self.division_parent(xp)
                return all, xp
            all.append("Right operand is numeric, will attempt conversion in both directions.")
        elif type(yp) is type:
            all.append("Right operand is not Sage element, will try _sage_.")

        if op is mul or op is imul:
            all.append("Will try _r_action and _l_action")

        return all, None

    cpdef Parent division_parent(self, Parent parent):
        r"""
        Deduces where the result of division in parent lies by calculating
        the inverse of \code{parent.one_element()} or \code{parent.an_element()}.

        The result is cached.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.division_parent(ZZ)
            Rational Field
            sage: cm.division_parent(QQ)
            Rational Field
            sage: cm.division_parent(ZZ['x'])
            Fraction Field of Univariate Polynomial Ring in x over Integer Ring
            sage: cm.division_parent(GF(41))
            Finite Field of size 41
            sage: cm.division_parent(Integers(100))
            Ring of integers modulo 100
            sage: cm.division_parent(SymmetricGroup(5))
            Symmetric group of order 5! as a permutation group
        """
        try:
            return self._division_parents.get(parent, None, None)
        except KeyError:
            pass
        try:
            ret = parent_c(~parent.one_element())
        except:
            self._record_exception()
            ret = parent_c(~parent.an_element())
        self._division_parents.set(parent, None, None, ret)
        return ret


    cpdef bin_op(self, x, y, op):
        """
        Execute the operation op on x and y. It first looks for an action
        corresponding to op, and failing that, it tries to coerces x and y
        into the a common parent and calls op on them.

        If it cannot make sense of the operation, a TypeError is raised.

        INPUT:
            x  -- the left operand
            y  -- the right operand
            op -- a python function taking 2 arguments
                  Note: op is often an arithmitic operation, but need
                  not be so.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.bin_op(1/2, 5, operator.mul)
            5/2

        The operator can be any callable:
            set Rational Field Integer Ring <function <lambda> at 0xc0b2270> None None
            (Rational Field, Rational Field)
            sage: R.<x> = ZZ['x']
            sage: cm.bin_op(x^2-1, x+1, gcd)
            x + 1

        Actions are detected and performed:
            sage: M = matrix(ZZ, 2, 2, range(4))
            sage: V = vector(ZZ, [5,7])
            sage: cm.bin_op(M, V, operator.mul)
            (7, 31)

        TESTS:
            sage: class Foo:
            ...      def __rmul__(self, left):
            ...          return 'hello'
            ...
            sage: H = Foo()
            sage: print int(3)*H
            hello
            sage: print Integer(3)*H
            hello
            sage: print H*3
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': '<type 'instance'>' and 'Integer Ring'

            sage: class Nonsense:
            ...       def __init__(self, s):
            ...           self.s = s
            ...       def __repr__(self):
            ...           return self.s
            ...       def __mul__(self, x):
            ...           return Nonsense(self.s + chr(x%256))
            ...       __add__ = __mul__
            ...       def __rmul__(self, x):
            ...           return Nonsense(chr(x%256) + self.s)
            ...       __radd__ = __rmul__
            ...
            sage: a = Nonsense('blahblah')
            sage: a*80
            blahblahP
            sage: 80*a
            Pblahblah
            sage: a+80
            blahblahP
            sage: 80+a
            Pblahblah

        """
        self._exceptions_cleared = False
        if (op is not sub) and (op is not isub):
            # Actions take preference over common-parent coercions.
            xp = parent_c(x)
            yp = parent_c(y)
            if xp is yp:
                return op(x,y)
            action = self.get_action(xp, yp, op)
            if action is not None:
                return (<Action>action)._call_(x, y)

        try:
            xy = self.canonical_coercion(x,y)
            return PyObject_CallObject(op, xy)
        except TypeError, err:
            self._record_exception()

        if op is mul or op is imul:

            # elements may also act on non-elements
            # (e.g. sequences or parents)
            try:
                return x._l_action(y)
            except (AttributeError, TypeError), err:
                self._record_exception()
            try:
                return y._r_action(x)
            except (AttributeError, TypeError), err:
                self._record_exception()

        if not isinstance(y, Element):
            try:
                op_name = op.__name__
                if op_name[0] == 'i':
                    op_name = op_name[1:]
                return getattr(y, '__r%s__'%op_name)(x)
            except (AttributeError, TypeError):
                pass

        # We should really include the underlying error.
        # This causes so much headache.
        raise TypeError, arith_error_message(x,y,op)

    cpdef canonical_coercion(self, x, y):
        """
        Given to elements x and y, with parents S and R respectively,
        find a common parent Z such that there are coercions
        $f: S \mapsto Z$ and $g: R \mapsto Z$ and return $f(x), g(y)$
        which will have the same parent.

        Raises a type error if no such Z can be found.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.canonical_coercion(mod(2, 10), 17)
            (2, 7)
            sage: x, y = cm.canonical_coercion(1/2, matrix(ZZ, 2, 2, range(4)))
            sage: x
            [1/2   0]
            [  0 1/2]
            sage: y
            [0 1]
            [2 3]
            sage: parent(x) is parent(y)
            True

        There is some support for non-Sage datatypes as well:
            sage: x, y = cm.canonical_coercion(int(5), 10)
            sage: type(x), type(y)
            (<type 'sage.rings.integer.Integer'>, <type 'sage.rings.integer.Integer'>)

            sage: class MyClass:
            ...       def _sage_(self):
            ...           return 13
            sage: a, b = cm.canonical_coercion(MyClass(), 1/3)
            sage: a, b
            (13, 1/3)
            sage: type(a)
            <type 'sage.rings.rational.Rational'>

        We also make an exception for 0, even if $\ZZ$ does not map in:
            sage: canonical_coercion(vector([1, 2, 3]), 0)
            ((1, 2, 3), (0, 0, 0))
        """
        xp = parent_c(x)
        yp = parent_c(y)
        if xp is yp:
            return x,y

        cdef Element x_elt, y_elt
        coercions = self.coercion_maps(xp, yp)
        if coercions is not None:
            x_map, y_map = coercions
            if x_map is not None:
                x_elt = (<Map>x_map)._call_(x)
            else:
                x_elt = x
            if y_map is not None:
                y_elt = (<Map>y_map)._call_(y)
            else:
                y_elt = y
            if x_elt is None:
                raise RuntimeError, "BUG in map, returned None %s %s %s" % (x, type(x_map), x_map)
            elif y_elt is None:
                raise RuntimeError, "BUG in map, returned None %s %s %s" % (y, type(y_map), y_map)
            if x_elt._parent is y_elt._parent:
                # We must verify this as otherwise we are prone to
                # getting into an infinite loop in c, and the above
                # maps may be written by (imperfect) users.
                return x_elt,y_elt
            elif x_elt._parent == y_elt._parent:
                # TODO: Non-uniqueness of parents strikes again!
                # print parent_c(x_elt), " is not ", parent_c(y_elt)
                y_elt = parent_c(x_elt)(y_elt)
                if x_elt._parent is y_elt._parent:
                    return x_elt,y_elt
            self._coercion_error(x, x_map, x_elt, y, y_map, y_elt)

        # Now handle the native python + sage object cases
        # that were not taken care of above.
        elif PY_IS_NUMERIC(x):
            try:
                x = yp(x)
                if PY_TYPE_CHECK(yp, type): return x,y
            except (TypeError, ValueError):
                self._record_exception()
                y = x.__class__(y)
                return x, y
            return _verify_canonical_coercion_c(x,y)

        elif PY_IS_NUMERIC(y):
            try:
                y = xp(y)
                if PY_TYPE_CHECK(xp, type): return x,y
            except (TypeError, ValueError):
                self._record_exception()
                x = y.__class__(x)
                return x, y
            return _verify_canonical_coercion_c(x,y)

        try:
            if not PY_TYPE_CHECK(x, SageObject) or not PY_TYPE_CHECK(y, SageObject):
                x = x._sage_()
                y = y._sage_()
                return self.canonical_coercion(x, y)
        except AttributeError:
            self._record_exception()

        # Allow coercion of 0 even if no coercion from Z
        if is_Integer(x) and not x and not PY_TYPE_CHECK_EXACT(yp, type):
            try:
                return yp(0), y
            except:
                self._record_exception()

        if is_Integer(y) and not y and not PY_TYPE_CHECK_EXACT(xp, type):
            try:
                return x, xp(0)
            except:
                self._record_exception()

        raise TypeError, "no common canonical parent for objects with parents: '%s' and '%s'"%(xp, yp)


    cpdef coercion_maps(self, R, S):
        """
        Give two parents R and S, return a pair of coercion maps
        $f: R \rightarrow Z$ and $g: S \rightarrow Z$, if such a $Z$
        can be found.

        In the (common) case that $R=Z$ or $S=Z$ then \code{None} is returned
        for $f$ or $g$ respectively rather than constructing (and subsequently
        calling) the identity morphism.

        If no suitable $f, g$ can be found, a single None is returned.
        This result is cached.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: f, g = cm.coercion_maps(ZZ, QQ)
            sage: print f
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
            sage: print g
            None

            sage: f, g = cm.coercion_maps(ZZ['x'], QQ)
            sage: print f
            Call morphism:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Univariate Polynomial Ring in x over Rational Field
            sage: print g
            Call morphism:
              From: Rational Field
              To:   Univariate Polynomial Ring in x over Rational Field

            sage: cm.coercion_maps(QQ, GF(7)) == None
            True

        Note that to break symmetry, if there is a coercion map in both
        directions, the parent on the left is used:
            sage: V = QQ^3
            sage: W = loads(dumps(V))
            sage: V == W
            True
            sage: V is W
            False
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.coercion_maps(V, W)
            (None,
             Call morphism:
              From: Vector space of dimension 3 over Rational Field
              To:   Vector space of dimension 3 over Rational Field)
            sage: cm.coercion_maps(W, V)
            (None,
             Call morphism:
              From: Vector space of dimension 3 over Rational Field
              To:   Vector space of dimension 3 over Rational Field)
            sage: v = V([1,2,3])
            sage: w = W([1,2,3])
            sage: parent(v+w) is V
            True
            sage: parent(w+v) is W
            True
        """
        try:
            return self._coercion_maps.get(R, S, None)
        except KeyError:
            homs = self.discover_coercion(R, S)
            if 0:
                # This breaks too many things that are going to change
                # in the new coercion model anyways.
                # COERCE TODO: Enable it then.
                homs = self.verify_coercion_maps(R, S, homs)
            else:
                if homs is not None:
                    x_map, y_map = homs
                    if x_map is not None and not isinstance(x_map, Map):
                        raise RuntimeError, "BUG in coercion model: coerce_map_from must return a Map"
                    if y_map is not None and not isinstance(y_map, Map):
                        raise RuntimeError, "BUG in coercion model: coerce_map_from must return a Map"
            if homs is None:
                swap = None
            else:
                R_map, S_map = homs
                if R_map is None and PY_TYPE_CHECK(S, Parent) and (<Parent>S).has_coerce_map_from(R):
                    swap = None, (<Parent>S).coerce_map_from(R)
                else:
                    swap = S_map, R_map
            self._coercion_maps.set(R, S, None, homs)
            self._coercion_maps.set(S, R, None, swap)
        return homs

    cpdef verify_coercion_maps(self, R, S, homs, bint fix=False):
        """
        Make sure this is a valid pair of homomorphisms from R and S to a common parent.
        This function is used to protect the user against buggy parents.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: homs = QQ.coerce_map_from(ZZ), None
            sage: cm.verify_coercion_maps(ZZ, QQ, homs) == homs
            True
            sage: homs = QQ.coerce_map_from(ZZ), RR.coerce_map_from(QQ)
            sage: cm.verify_coercion_maps(ZZ, QQ, homs) == homs
            Traceback (most recent call last):
            ...
            RuntimeError: ('BUG in coercion model, codomains must be identical', Natural morphism:
              From: Integer Ring
              To:   Rational Field, Generic map:
              From: Rational Field
              To:   Real Field with 53 bits of precision)
        """
        if homs is None:
            return None
        cdef Map x_map, y_map
        R_map, S_map = homs
        if PY_TYPE_CHECK(R, type):
            R = Set_PythonType(R)
        elif PY_TYPE_CHECK(S, type):
            S = Set_PythonType(S)
        if R_map is None:
            R_map = IdentityMorphism(R)
        elif S_map is None:
            S_map = IdentityMorphism(S)
        # Make sure the domains are correct
        if R_map.domain() is not R:
            if fix:
                connecting = R_map.domain().coerce_map_from(R)
                if connecting is not None:
                    R_map = R_map * connecting
            if R_map.domain() is not R:
                raise RuntimeError, ("BUG in coercion model, left domain must be original parent", R, R_map)
        if S_map is not None and S_map.domain() is not S:
            if fix:
                connecting = S_map.domain().coerce_map_from(S)
                if connecting is not None:
                    S_map = S_map * connecting
            if S_map.domain() is not S:
                raise RuntimeError, ("BUG in coercion model, right domain must be original parent", S, S_map)
        # Make sure the codomains are correct
        if R_map.codomain() is not S_map.codomain():
            if fix:
                connecting = R_map.codomain().coerce_map_from(S_map.codomain())
                if connecting is not None:
                    S_map = connecting * S_map
                else:
                    connecting = S_map.codomain().coerce_map_from(R_map.codomain())
                    if connecting is not None:
                        R_map = connecting * R_map
            if R_map.codomain() is not S_map.codomain():
                raise RuntimeError, ("BUG in coercion model, codomains must be identical", R_map, S_map)
        if PY_TYPE_CHECK(R_map, IdentityMorphism):
            R_map = None
        elif PY_TYPE_CHECK(S_map, IdentityMorphism):
            S_map = None
        return R_map, S_map


    cpdef discover_coercion(self, R, S):
        """
        This actually implements the finding of coercion maps as described in
        the \code{coercion_maps} method.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()

        If R is S, then two identity morphisms suffice:
            sage: cm.discover_coercion(SR, SR)
            (None, None)

        If there is a coercion map either direction, use that:
            sage: cm.discover_coercion(ZZ, QQ)
            (Natural morphism:
              From: Integer Ring
              To:   Rational Field, None)
            sage: cm.discover_coercion(RR, QQ)
            (None,
             Generic map:
              From: Rational Field
              To:   Real Field with 53 bits of precision)

        Otherwise, try and compute an appropriate cover:
            sage: cm.discover_coercion(ZZ['x,y'], RDF)
            (Call morphism:
              From: Multivariate Polynomial Ring in x, y over Integer Ring
              To:   Multivariate Polynomial Ring in x, y over Real Double Field,
             Call morphism:
              From: Real Double Field
              To:   Multivariate Polynomial Ring in x, y over Real Double Field)

        Sometimes there is a reasonable "cover," but no canonical coercion:
            sage: sage.categories.pushout.pushout(QQ, QQ^3)
            Vector space of dimension 3 over Rational Field
            sage: print cm.discover_coercion(QQ, QQ^3)
            None
        """
        from sage.categories.homset import Hom
        if R is S:
            return None, None

        # See if there is a natural coercion from R to S
        if PY_TYPE_CHECK(R, Parent):
            mor = (<Parent>R).coerce_map_from(S)
            if mor is not None:
                return None, mor

        # See if there is a natural coercion from S to R
        if PY_TYPE_CHECK(S, Parent):
            mor = (<Parent>S).coerce_map_from(R)
            if mor is not None:
                return mor, None

        # Try base extending
        if PY_TYPE_CHECK(R, Parent) and PY_TYPE_CHECK(S, Parent):
            from sage.categories.pushout import pushout
            try:
                Z = pushout(R, S)
                coerce_R = Z.coerce_map_from(R)
                coerce_S = Z.coerce_map_from(S)
                if coerce_R is None:
                    raise TypeError, "No coercion from %s to pushout %s" % (R, Z)
                if coerce_S is None:
                    raise TypeError, "No coercion from %s to pushout %s" % (S, Z)
                return coerce_R, coerce_S
            except:
                self._record_exception()

        return None


    cpdef get_action(self, R, S, op):
        """
        Get the action of R on S or S on R associated to the operation op.

        EXAMPLES:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.get_action(ZZ['x'], ZZ, operator.mul)
            Right scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Integer Ring
            sage: cm.get_action(ZZ['x'], ZZ, operator.imul)
            Right scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Integer Ring
            sage: cm.get_action(ZZ['x'], QQ, operator.mul)
            Right scalar multiplication by Rational Field on Univariate Polynomial Ring in x over Integer Ring
            sage: cm.get_action(QQ['x'], int, operator.mul)
            Right scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Rational Field
            with precomposition on right by Native morphism:
              From: Set of Python objects of type 'int'
              To:   Integer Ring

            sage: R.<x> = QQ['x']
            sage: A = cm.get_action(R, ZZ, operator.div); A
            Right inverse action by Rational Field on Univariate Polynomial Ring in x over Rational Field
            with precomposition on right by Natural morphism:
              From: Integer Ring
              To:   Rational Field
            sage: A(x+10, 5)
            1/5*x + 2

        """
        try:
            return self._action_maps.get(R, S, op)
        except KeyError:
            action = self.discover_action(R, S, op)
            action = self.verify_action(action, R, S, op)
            self._action_maps.set(R, S, op, action)
            return action

    cpdef verify_action(self, action, R, S, op, bint fix=True):
        r"""
        Verify that \code{action} takes an element of R on the left and S
        on the right, raising an error if not.

        This is used for consistency checking in the coercion model.

        EXAMPLES:
            sage: R.<x> = ZZ['x']
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.verify_action(R.get_action(QQ), R, QQ, operator.mul)
            Right scalar multiplication by Rational Field on Univariate Polynomial Ring in x over Integer Ring
            sage: cm.verify_action(R.get_action(QQ), RDF, R, operator.mul)
            Traceback (most recent call last):
            ...
            RuntimeError: There is a BUG in the coercion model:
                Action found for R <built-in function mul> S does not have the correct domains
                R = Real Double Field
                S = Univariate Polynomial Ring in x over Integer Ring
                (should be Univariate Polynomial Ring in x over Integer Ring, Rational Field)
                action = Right scalar multiplication by Rational Field on Univariate Polynomial Ring in x over Integer Ring (<type 'sage.structure.coerce_actions.RightModuleAction'>)
        """
        if action is None:
            return action
        elif PY_TYPE_CHECK(action, IntegerMulAction):
            return action
        cdef bint ok = True
        try:
            if action.left_domain() is not R:
                ok &= PY_TYPE_CHECK(R, type) and action.left_domain()._type is R
            if action.right_domain() is not S:
                ok &= PY_TYPE_CHECK(S, type) and action.right_domain()._type is S
        except AttributeError:
            ok = False
        if not ok:
            if PY_TYPE_CHECK(R, type):
                R = Set_PythonType(R)
            if PY_TYPE_CHECK(S, type):
                S = Set_PythonType(S)

            # Non-unique parents
            if fix and action.left_domain() is not R and action.left_domain() == R:
                action = PrecomposedAction(action, action.left_domain().coerce_map_from(R), None)
            if fix and action.right_domain() is not S and action.right_domain() == S:
                action = PrecomposedAction(action, None, action.right_domain().coerce_map_from(S))

            if action.left_domain() is not R or action.right_domain() is not S:
                raise RuntimeError, """There is a BUG in the coercion model:
                Action found for R %s S does not have the correct domains
                R = %s
                S = %s
                (should be %s, %s)
                action = %s (%s)
                """ % (op, R, S, action.left_domain(), action.right_domain(), action, type(action))

        return action

    cpdef discover_action(self, R, S, op):
        """
        INPUT
            R -- the left Parent (or type)
            S -- the right Parent (or type)
            op -- the operand, typically an element of the operator module.

        OUTPUT:
            An action A such that s op r is given by A(s,r).

        The steps taken are illistrated below.

        EXAMPLES:
            sage: P.<x> = ZZ['x']
            sage: P.get_action(ZZ)
            Right scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Integer Ring
            sage: ZZ.get_action(P) is None
            True
            sage: cm = sage.structure.element.get_coercion_model()

        If R or S is a Parent, ask it for an action by/on R:
            sage: cm.discover_action(ZZ, P, operator.mul)
            Left scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Integer Ring

        If R or S a type, recursively call get_action with the Sage versions of R and/or S:
            sage: cm.discover_action(P, int, operator.mul)
            Right scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Integer Ring
            with precomposition on right by Native morphism:
              From: Set of Python objects of type 'int'
              To:   Integer Ring

        If op in an inplace operation, look for the non-inplace action:
            sage: cm.discover_action(P, ZZ, operator.imul)
            Right scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Integer Ring

        If op is division, look for action on right by inverse:
            sage: cm.discover_action(P, ZZ, operator.div)
            Right inverse action by Rational Field on Univariate Polynomial Ring in x over Integer Ring
            with precomposition on right by Natural morphism:
              From: Integer Ring
              To:   Rational Field
        """
        #print "looking", R, <int><void *>R, op, S, <int><void *>S

        if PY_TYPE_CHECK(S, Parent):
            action = (<Parent>S).get_action(R, op, False)
            if action is not None:
                #print "found1", action
                return action

        if PY_TYPE_CHECK(R, Parent):
            action = (<Parent>R).get_action(S, op, True)
            if action is not None:
                #print "found2", action
                return action

        if PY_TYPE(R) == <void *>type:
            sageR = py_scalar_parent(R)
            if sageR is not None:
                action = self.discover_action(sageR, S, op)
                if action is not None:
                    if not PY_TYPE_CHECK(action, IntegerMulAction):
                        action = PrecomposedAction(action, sageR.coerce_map_from(R), None)
                    return action

        if PY_TYPE(S) == <void *>type:
            sageS = py_scalar_parent(S)
            if sageS is not None:
                action = self.discover_action(R, sageS, op)
                if action is not None:
                    if not PY_TYPE_CHECK(action, IntegerMulAction):
                        action = PrecomposedAction(action, None, sageS.coerce_map_from(S))
                    return action

        if op.__name__[0] == 'i':
            try:
                a = self.discover_action(R, S, no_inplace_op(op))
                if a is not None:
                    is_inverse = isinstance(a, InverseAction)
                    if is_inverse: a = ~a
                    if a is not None and PY_TYPE_CHECK(a, RightModuleAction):
                        # We want a new instance so that we don't alter the (potentially cached) original
                        a = RightModuleAction(S, R)
                        a.is_inplace = 1
                    if is_inverse: a = ~a
                return a
            except KeyError:
                self._record_exception()

        if op is div:
            # Division on right is the same acting on right by inverse, if it is so defined.
            # To return such an action, we need to verify that it would be an action for the mul
            # operator, but the action must be over a parent containing inverse elements.
            from sage.rings.ring import is_Ring
            if is_Ring(S):
                try:
                    K = S._pseudo_fraction_field()
                except (TypeError, AttributeError):
                    K = None
            elif PY_TYPE_CHECK(S, Parent):
                K = S
            else:
                # python scalar case handled recursively above
                K = None

            if K is not None:
                action = self.get_action(R, K, mul)
                if action is not None and action.actor() is K:
                    try:
                        action = ~action
                        if K is not S:
                            action = PrecomposedAction(action, None, K.coerce_map_from(S))
                        return action
                    except TypeError: # action may not be invertible
                        self._record_exception()

        return None

    def _coercion_error(self, x, x_map, x_elt, y, y_map, y_elt):
        """
        This function is only called when someone has incorrectly implemented
        a user-defined part of the coercion system (usually, a morphism).

        EXAMPLE:
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm._coercion_error('a', 'f', 'f(a)', 'b', 'g', 'g(b)')
            Traceback (most recent call last):
            ...
            RuntimeError: There is a bug in the coercion code in SAGE.
            Both x (='f(a)') and y (='g(b)') are supposed to have identical parents but they don't.
            In fact, x has parent '<type 'str'>'
            whereas y has parent '<type 'str'>'
            Original elements 'a' (parent <type 'str'>) and 'b' (parent <type 'str'>) and maps
            <type 'str'> 'f'
            <type 'str'> 'g'
        """
        raise RuntimeError, """There is a bug in the coercion code in SAGE.
Both x (=%r) and y (=%r) are supposed to have identical parents but they don't.
In fact, x has parent '%s'
whereas y has parent '%s'
Original elements %r (parent %s) and %r (parent %s) and maps
%s %r
%s %r"""%( x_elt, y_elt, parent_c(x_elt), parent_c(y_elt),
            x, parent_c(x), y, parent_c(y),
            type(x_map), x_map, type(y_map), y_map)


