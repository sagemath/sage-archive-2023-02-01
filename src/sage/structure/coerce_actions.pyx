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

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"
from cpython.int cimport *
from cpython.number cimport *
include "coerce.pxi"

from coerce_exceptions import CoercionException

cdef extern from *:
    ctypedef struct RefPyObject "PyObject":
        int ob_refcnt

cdef _record_exception():
    from element import get_coercion_model
    get_coercion_model()._record_exception()

cdef inline an_element(R):
    if isinstance(R, Parent):
        return R.an_element()
    else:
        for x in ([(1, 2)], "abc", 10.5, 10):
            try:
                return R(x)
            except Exception:
                pass

cdef class LAction(Action):
    """Action calls _l_action of the actor."""
    def __init__(self, G, S):
        Action.__init__(self, G, S, True, operator.mul)
    cpdef _call_(self, g, a):
        return g._l_action(a)  # a * g


cdef class RAction(Action):
    """Action calls _r_action of the actor."""
    def __init__(self, G, S):
        Action.__init__(self, G, S, False, operator.mul)
    cpdef _call_(self, a, g):
        return g._r_action(a)  # g * a


# In the code below, I take the convention that g is acting on a.

cdef class GenericAction(Action):

    cdef _codomain

    def __init__(self, Parent G, S, is_left, bint check=True):
        """
        TESTS::

            sage: sage.structure.coerce_actions.ActedUponAction(MatrixSpace(ZZ, 2), Cusps, True)
            Left action by Full MatrixSpace of 2 by 2 dense matrices over Integer Ring on Set P^1(QQ) of all cusps

            sage: sage.structure.coerce_actions.GenericAction(QQ, Zmod(6), True)
            Traceback (most recent call last):
            ...
            NotImplementedError: Action not implemented.

        This will break if we tried to use it::

            sage: sage.structure.coerce_actions.GenericAction(QQ, Zmod(6), True, check=False)
            Left action by Rational Field on Ring of integers modulo 6

        """
        Action.__init__(self, G, S, is_left, operator.mul)
        if check:
            res = self.act(G.an_element(), S.an_element())
            if res is None:
                raise CoercionException
            _codomain = parent_c(res)

    def codomain(self):
        """
        Returns the "codomain" of this action, i.e. the Parent in which the
        result elements live. Typically, this should be the same as the
        acted upon set.

        EXAMPLES::

            sage: A = sage.structure.coerce_actions.ActedUponAction(MatrixSpace(ZZ, 2), Cusps, True)
            sage: A.codomain()
            Set P^1(QQ) of all cusps
            sage: A = sage.structure.coerce_actions.ActOnAction(SymmetricGroup(3), QQ['x,y,z'], False)
            sage: A.codomain()
            Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        if self._codomain is None:
            self._codomain = parent_c(self.act(an_element(self.G),
                                               an_element(self.underlying_set())))
        return self._codomain


cdef class ActOnAction(GenericAction):
    """
    Class for actions defined via the _act_on_ method.
    """
    cpdef _call_(self, a, b):
        """
        TESTS::

            sage: G = SymmetricGroup(3)
            sage: R.<x,y,z> = QQ[]
            sage: A = sage.structure.coerce_actions.ActOnAction(G, R, False)
            sage: A._call_(x^2 + y - z, G((1,2)))
            y^2 + x - z
            sage: A._call_(x+2*y+3*z, G((1,3,2)))
            2*x + 3*y + z

            sage: type(A)
            <type 'sage.structure.coerce_actions.ActOnAction'>
        """
        if self._is_left:
            return (<Element>a)._act_on_(b, True)
        else:
            return (<Element>b)._act_on_(a, False)

cdef class ActedUponAction(GenericAction):
    """
    Class for actions defined via the _acted_upon_ method.
    """
    cpdef _call_(self, a, b):
        """
        TESTS::

            sage: A = sage.structure.coerce_actions.ActedUponAction(MatrixSpace(ZZ, 2), Cusps, True)
            sage: A.act(matrix(ZZ, 2, [1,0,2,-1]), Cusp(1,2))
            Infinity
            sage: A._call_(matrix(ZZ, 2, [1,0,2,-1]), Cusp(1,2))
            Infinity

            sage: type(A)
            <type 'sage.structure.coerce_actions.ActedUponAction'>
        """
        if self._is_left:
            return (<Element>b)._acted_upon_(a, False)
        else:
            return (<Element>a)._acted_upon_(b, True)

def detect_element_action(Parent X, Y, bint X_on_left, X_el=None, Y_el=None):
    r"""
    Returns an action of X on Y or Y on X as defined by elements X, if any.

    EXAMPLES::

        sage: from sage.structure.coerce_actions import detect_element_action
        sage: detect_element_action(ZZ['x'], ZZ, False)
        Left scalar multiplication by Integer Ring on Univariate Polynomial Ring in x over Integer Ring
        sage: detect_element_action(ZZ['x'], QQ, True)
        Right scalar multiplication by Rational Field on Univariate Polynomial Ring in x over Integer Ring
        sage: detect_element_action(Cusps, MatrixSpace(ZZ, 2), False)
        Left action by Full MatrixSpace of 2 by 2 dense matrices over Integer Ring on Set P^1(QQ) of all cusps
        sage: detect_element_action(Cusps, MatrixSpace(ZZ, 2), True),
        (None,)
        sage: detect_element_action(ZZ, QQ, True),
        (None,)

    TESTS:

    This test checks that the issue in :trac:`7718` has been fixed::

        sage: class MyParent(Parent):
        ....:     def an_element(self):
        ....:         pass
        ....:
        sage: A = MyParent()
        sage: detect_element_action(A, ZZ, True)
        Traceback (most recent call last):
        ...
        RuntimeError: an_element() for <class '__main__.MyParent'> returned None
    """
    cdef Element x
    if X_el is None or (parent_c(X_el) is not X):
        x = an_element(X)
    else:
        x = X_el
    if x is None:
        raise RuntimeError("an_element() for %s returned None" % X)
    if Y_el is None or (parent_c(Y_el) is not Y):
        y = an_element(Y)
    else:
        y = Y_el
    if y is None:
        if isinstance(Y, Parent):
            raise RuntimeError("an_element() for %s returned None" % Y)
        else:
            return # don't know how to make elements of this type...
    if isinstance(x, ModuleElement) and isinstance(y, RingElement):
        # Elements defining _lmul_ and _rmul_
        try:
            return (RightModuleAction if X_on_left else LeftModuleAction)(Y, X, y, x)
        except CoercionException as msg:
            _record_exception()
    try:
        # Elements defining _act_on_
        if x._act_on_(y, X_on_left) is not None:
            return ActOnAction(X, Y, X_on_left, False)
    except CoercionException:
        _record_exception()
    if isinstance(Y, Parent):
        try:
            # Elements defining _acted_upon_
            if x._acted_upon_(y, X_on_left) is not None:
                return ActedUponAction(Y, X, not X_on_left, False)
        except CoercionException:
            _record_exception()


cdef class ModuleAction(Action):
    """
    Module action.

    .. SEEALSO::

        This is an abstract class, one must actually instantiate a
        :class:`LeftModuleAction` or a :class:`RightModuleAction`.

    INPUT:

    - ``G`` -- the actor, an instance of :class:`~sage.structure.parent.Parent`.
    - ``S`` -- the object that is acted upon.
    - ``g`` -- optional, an element of ``G``.
    - ``a`` -- optional, an element of ``S``.
    - ``check`` -- if True (default), then there will be no consistency tests
      performed on sample elements.

    NOTE:

    By default, the sample elements of ``S`` and ``G`` are obtained from
    :meth:`~sage.structure.parent.Parent.an_element`, which relies on the
    implementation of an ``_an_element_()`` method. This is not always
    awailable. But usually, the action is only needed when one already
    *has* two elements. Hence, by :trac:`14249`, the coercion model will
    pass these two elements the the :class:`ModuleAction` constructor.

    The actual action is implemented by the ``_rmul_`` or ``_lmul_``
    function on its elements. We must, however, be very particular about
    what we feed into these functions, because they operate under the
    assumption that the inputs lie exactly in the base ring and may
    segfault otherwise. Thus we handle all possible base extensions
    manually here.

    """
    def __init__(self, G, S, g=None, a=None, check=True):
        """
        This creates an action of an element of a module by an element of its
        base ring. The simplest example to keep in mind is R acting on the
        polynomial ring R[x].

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

        The following tests against a problem that was relevant during work on
        trac ticket #9944::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: S.<x> = PolynomialRing(ZZ, sparse=True)
            sage: 1/R.0
            1/x
            sage: 1/S.0
            1/x

        """
        Action.__init__(self, G, S, not PY_TYPE_CHECK(self, RightModuleAction), operator.mul)
        if not isinstance(G, Parent):
            # only let Parents act
            raise TypeError, "Actor must be a parent."
        base = S.base()
        if base is S or base is None:
            # The right thing to do is a normal multiplication
            raise CoercionException, "Best viewed as standard multiplication"
        # Objects are implemented with the assumption that
        # _rmul_ is given an element of the base ring
        if G is not base:
            # first we try the easy case of coercing G to the base ring of S
            self.connecting = base.coerce_map_from(G)
            if self.connecting is None:
                # otherwise, we try and find a base extension
                from sage.categories.pushout import pushout
                # this may raise a type error, which we propagate
                self.extended_base = pushout(G, S)
                # make sure the pushout actually gave correct a base extension of S
                if self.extended_base.base() != pushout(G, base):
                    raise CoercionException, "Actor must be coercible into base."
                else:
                    self.connecting = self.extended_base.base().coerce_map_from(G)
                    if self.connecting is None:
                        # this may happen if G is, say, int rather than a parent
                        # TODO: let python types be valid actions
                        raise CoercionException, "Missing connecting morphism"

        # Don't waste time if our connecting morphisms happen to be the identity.
        if self.connecting is not None and self.connecting.codomain() is G:
            self.connecting = None

        if self.extended_base is not None and self.extended_base is S:
            self.extended_base = None

        # At this point, we can assert it is safe to call _Xmul_c
        the_ring = G if self.connecting is None else self.connecting.codomain()
        the_set = S if self.extended_base is None else self.extended_base
        assert the_ring is the_set.base(), "BUG in coercion model\n    Apparently there are two versions of\n        %s\n    in the cache."%the_ring

        if not check:
            return
        if g is None:
            g = G.an_element()
        if parent_c(g) is not G:
            raise CoercionException("The parent of %s is not %s but %s"%(g,G,parent_c(g)))
        if a is None:
            a = S.an_element()
        if parent_c(a) is not S:
            raise CoercionException("The parent of %s is not %s but %s"%(a,S,parent_c(a)))
        if not isinstance(g, RingElement) or not isinstance(a, ModuleElement):
            raise CoercionException("Not a ring element acting on a module element.")
        res = self.act(g, a)
        if parent_c(res) is not the_set:
            # In particular we will raise an error if res is None
            raise CoercionException("Result is None or has wrong parent.")


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
        return self.underlying_set()

    def domain(self):
        """
        The domain of self, which is the module that is being acted on.

        EXAMPLES::

            sage: from sage.structure.coerce_actions import LeftModuleAction
            sage: A = LeftModuleAction(QQ, ZZ['x,y,z'])
            sage: A.domain()
            Multivariate Polynomial Ring in x, y, z over Integer Ring
        """
        return self.underlying_set()



cdef class LeftModuleAction(ModuleAction):

    cpdef _call_(self, g, a):
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

    cpdef _call_(self, a, g):
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
                # This is a truly new object, mutate it
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

    def __init__(self, ZZ, M, is_left=True, m=None):
        r"""
        This class implements the action `n \cdot a = a + a + \cdots + a` via
        repeated doubling.

        Both addition and negation must be defined on the set `M`.

        INPUT:

        - An integer ring, ``ZZ``
        - A ``ZZ`` module ``M``
        - Optional: An element ``m`` of ``M``

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
        if m is None:
            m = M.an_element()
        test = m + (-m) # make sure addition and negation is allowed
        Action.__init__(self, ZZ, M, is_left, operator.mul)

    cpdef _call_(self, nn, a):
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
    # We cannot use sig_on() here because this may call arbitrary Python
    # code raising exceptions. -- Jeroen Demeyer
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
    return sum

cdef inline fast_mul_long(a, long n):
    # We cannot use sig_on() here because this may call arbitrary Python
    # code raising exceptions. -- Jeroen Demeyer
    if n < 0:
        n = -n
        a = -a
    if n < 4:
        if n == 0: return parent_c(a)(0)
        elif n == 1: return a
        elif n == 2: return a+a
        elif n == 3: return a+a+a
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
    return sum
