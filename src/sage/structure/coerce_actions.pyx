"""
Coerce actions
"""
#*****************************************************************************
#     Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                    http://www.gnu.org/licenses/
#*****************************************************************************

import operator

include "../ext/stdsage.pxi"
include "../ext/interrupt.pxi"
include "../ext/python_int.pxi"
include "../ext/python_number.pxi"
include "coerce.pxi"

cdef extern from *:
    ctypedef struct RefPyObject "PyObject":
        int ob_refcnt



cdef class LAction(Action):
    """Action calls _l_action of the actor."""
    def __init__(self, G, S):
        Action.__init__(self, G, S, True, operator.mul)
    cpdef Element _call_(self, g, a):
        return g._l_action(a)  # a * g


cdef class RAction(Action):
    """Action calls _r_action of the actor."""
    def __init__(self, G, S):
        Action.__init__(self, G, S, False, operator.mul)
    cpdef Element _call_(self, a, g):
        return g._r_action(a)  # g * a



cdef class ModuleAction(Action):

    def __init__(self, G, S):
        """
        This creates an action of an element of a module by an element of its
        base ring. The simplest example to keep in mind is R acting on the
        polynomial ring R[x].

        The actual action is implemented by the _rmul_ or _lmul_ function on
        its elements. We must, however, be very particular about what we feed
        into these functions because they operate under the assumption that
        the inputs lie exactly in the base ring and may segfault otherwise.

        Thus we handle all possible base extensions manually here. This is
        an abstract class, one must actually instantiate a LeftModuleAction
        or a RightModuleAction

        EXAMPLES::

            sage: from sage.structure.coerce_actions import LeftModuleAction
            sage: LeftModuleAction(ZZ, ZZ['x'])
            Left scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Integer Ring
            sage: LeftModuleAction(ZZ, QQ['x'])
            Left scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Rational Field
            sage: LeftModuleAction(QQ, ZZ['x'])
            Left scalar multiplication by Rational Field on Univariate Polynomial Ring in x over Integer Ring
            sage: LeftModuleAction(QQ, ZZ['x']['y'])
            Left scalar multiplication by Rational Field on Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integer Ring
        """
        Action.__init__(self, G, S, not PY_TYPE_CHECK(self, RightModuleAction), operator.mul)
        if not isinstance(G, Parent):
            # only let Parents act
            raise TypeError, "Actor must be a parent."
        if S.base() is S:
            # The right thing to do is a normal multiplication
            raise TypeError
        # Objects are implemented with the assumption that
        # _rmul_ is given an element of the basering
        if G is not S.base():
            # first we try the easy case of coercing G to the basering of S
            self.connecting = S.base().coerce_map_from(G)
            if self.connecting is None:
                # otherwise, we try and find a base extension
                from sage.categories.pushout import pushout
                # this may raise a type error, which we propagate
                self.extended_base = pushout(G, S)
                # make sure the pushout actually gave correct a base extension of S
                if self.extended_base.base() != pushout(G, S.base()):
                    raise TypeError, "Actor must be coercible into base."
                else:
                    self.connecting = self.extended_base.base().coerce_map_from(G)
                    if self.connecting is None:
                        # this may happen if G is, say, int rather than a parent
                        # TODO: let python types be valid actions
                        raise TypeError

        # Don't waste time if our connecting morphisms happen to be the identity.
        if self.connecting is not None and self.connecting.codomain() is G:
            self.connecting = None

        if self.extended_base is not None and self.extended_base is S:
            self.extended_base = None

        # At this point, we can assert it is safe to call _Xmul_c
        the_ring = G if self.connecting is None else self.connecting.codomain()
        the_set = S if self.extended_base is None else self.extended_base
        assert the_ring is the_set.base(), "BUG in coercion model"

        cpdef RingElement g = G.an_element()
        cpdef ModuleElement a = S.an_element()
        res = self.act(g, a)

        if parent_c(res) is not S and parent_c(res) is not self.extended_base:
            raise TypeError


    def _repr_name_(self):
        """
        The default name of this action type, which is has a sane default.

        EXAMPLES::

            sage: from sage.structure.coerce_actions import LeftModuleAction, RightModuleAction
            sage: A = LeftModuleAction(ZZ, ZZ['x']); A
            Left scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Integer Ring
            sage: A._repr_name_()
            'scalar multiplication'
            sage: RightModuleAction(GF(5), GF(5)[['t']])
            Right scalar multiplication by Finite Field of size 5 on Power Series Ring in t over Finite Field of size 5
        """
        return "scalar multiplication"

    def codomain(self):
        """
        The codomain of self, which may or may not be equal to the domain.

        EXAMPLES::

            sage: from sage.structure.coerce_actions import LeftModuleAction
            sage: A = LeftModuleAction(QQ, ZZ['x,y,z'])
            sage: A.codomain()
            Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        if self.extended_base is not None:
            return self.extended_base
        return self.S

    def domain(self):
        """
        The domain of self, which is the module that is being acted on.

        EXAMPLES::

            sage: from sage.structure.coerce_actions import LeftModuleAction
            sage: A = LeftModuleAction(QQ, ZZ['x,y,z'])
            sage: A.domain()
            Multivariate Polynomial Ring in x, y, z over Integer Ring
        """
        return self.S



cdef class LeftModuleAction(ModuleAction):

    cpdef Element _call_(self, g, a):
        """
        A left module action is an action that takes the ring element as the
        first argument (the left side) and the module element as the second
        argument (the right side).

        EXAMPLES::

            sage: from sage.structure.coerce_actions import LeftModuleAction
            sage: R.<x> = QQ['x']
            sage: A = LeftModuleAction(ZZ, R)
            sage: A(5, x+1)
            5*x + 5
            sage: R.<x> = ZZ['x']
            sage: A = LeftModuleAction(QQ, R)
            sage: A(1/2, x+1)
            1/2*x + 1/2
            sage: A._call_(1/2, x+1) # safe only when arguments have exactly the correct parent
            1/2*x + 1/2
        """
        if self.connecting is not None:
            g = self.connecting._call_(g)
        if self.extended_base is not None:
            a = self.extended_base(a)
        return (<ModuleElement>a)._rmul_(<RingElement>g)  # a * g


cdef class RightModuleAction(ModuleAction):

    def __cinit__(self):
        self.is_inplace = 0

    cpdef Element _call_(self, a, g):
        """
        A right module action is an action that takes the module element as the
        first argument (the left side) and the ring element as the second
        argument (the right side).

        EXAMPLES::

            sage: from sage.structure.coerce_actions import RightModuleAction
            sage: R.<x> = QQ['x']
            sage: A = RightModuleAction(ZZ, R)
            sage: A(x+5, 2)
            2*x + 10
            sage: A._call_(x+5, 2) # safe only when arguments have exactly the correct parent
            2*x + 10
        """
        cdef PyObject* tmp
        if self.connecting is not None:
            g = self.connecting._call_(g)
        if self.extended_base is not None:
            a = self.extended_base(a)
            # TODO: figure out where/why the polynomial constructor is caching 'a'
            if (<RefPyObject *>a).ob_refcnt == 2:
                b = self.extended_base(0)
            if (<RefPyObject *>a).ob_refcnt == 1:
                # This is a truely new object, mutate it
                return (<ModuleElement>a)._ilmul_(<RingElement>g)  # a * g
            else:
                return (<ModuleElement>a)._lmul_(<RingElement>g)  # a * g
        else:
            # If we have few enough references to this object, then we know
            # it is safe to do a (mutating) inplace operation.
            if (<RefPyObject *>a).ob_refcnt < inplace_threshold + self.is_inplace:
                return (<ModuleElement>a)._ilmul_(<RingElement>g)  # a * g
            else:
                return (<ModuleElement>a)._lmul_(<RingElement>g)  # a * g



cdef class IntegerMulAction(Action):

    def __init__(self, ZZ, M, is_left=True):
        r"""
        This class implements the action `n \cdot a = a + a + \cdots + a` via
        repeated doubling.

        Both addition and negation must be defined on the set `A`.

        EXAMPLES::

            sage: from sage.structure.coerce_actions import IntegerMulAction
            sage: R.<x> = QQ['x']
            sage: act = IntegerMulAction(ZZ, R)
            sage: act(5, x)
            5*x
            sage: act(0, x)
            0
            sage: act(-3, x-1)
            -3*x + 3
        """
        if PY_TYPE_CHECK(ZZ, type):
            from sage.structure.parent import Set_PythonType
            ZZ = Set_PythonType(ZZ)
        test = M.an_element() + (-M.an_element()) # make sure addition and negation is allowed
        Action.__init__(self, ZZ, M, is_left, operator.mul)

    cpdef Element _call_(self, nn, a):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_actions import IntegerMulAction
            sage: act = IntegerMulAction(ZZ, GF(101))
            sage: act(3, 9)
            27
            sage: act(3^689, 9)
            42
            sage: 3^689 * mod(9, 101)
            42

        Use round off error to verify this is doing actual repeated addition
        instead of just multiplying::

            sage: act = IntegerMulAction(ZZ, RR)
            sage: act(49, 1/49) == 49*RR(1/49)
            False
        """
        if not self._is_left:
            a, nn = nn, a
        if not PyInt_CheckExact(nn):
            nn = PyNumber_Int(nn)
            if not PyInt_CheckExact(nn):
                return fast_mul(a, nn)

        return fast_mul_long(a, PyInt_AS_LONG(nn))

    def __invert__(self):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_actions import IntegerMulAction
            sage: act = IntegerMulAction(ZZ, CDF)
            sage: ~act
            Traceback (most recent call last):
            ...
            TypeError: No generic module division by Z.
        """
        raise TypeError, "No generic module division by Z."

    def _repr_name_(self):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_actions import IntegerMulAction
            sage: IntegerMulAction(ZZ, GF(5))
            Left Integer Multiplication by Integer Ring on Finite Field of size 5
        """
        return "Integer Multiplication"



cdef inline fast_mul(a, n):
    _sig_on
    if n < 0:
        n = -n
        a = -a
    pow2a = a
    while n & 1 == 0:
        pow2a += pow2a
        n = n >> 1
    sum = pow2a
    n = n >> 1
    while n != 0:
        pow2a += pow2a
        if n & 1:
            sum += pow2a
        n = n >> 1
    _sig_off
    return sum

cdef inline fast_mul_long(a, long n):
    if n < 0:
        n = -n
        a = -a
    if n < 4:
        if n == 0: return parent_c(a)(0)
        if n == 1: return a
        if n == 2: return a+a
        if n == 3: return a+a+a
    _sig_on
    pow2a = a
    while n & 1 == 0:
        pow2a += pow2a
        n = n >> 1
    sum = pow2a
    n = n >> 1
    while n != 0:
        pow2a += pow2a
        if n & 1:
            sum += pow2a
        n = n >> 1
    _sig_off
    return sum
