"""
Coerce maps
"""
include "sage/ext/stdsage.pxi"

import re
import types

from parent import Set_PythonType
from sage.structure.parent cimport Parent
from sage.structure.element cimport Element

cdef object BuiltinMethodType = type(repr)

# COERCE TODO: remove or integrate better (as this bit is only checked on an error)
cdef bint print_warnings = 0


cdef class DefaultConvertMap(Map):
    """
    This morphism simply calls the codomain's element_constructor method,
    passing in the codomain as the first argument.
    """
    def __init__(self, domain, codomain, force_use=False):
        if not PY_TYPE_CHECK(domain, Parent):
            domain = Set_PythonType(domain)
        Map.__init__(self, domain, codomain)
        self._coerce_cost = 100
        self._force_use = force_use
        if self._codomain._element_constructor is None:
            raise RuntimeError, "BUG in coercion model, no element constructor for %s" % type(self._codomain)
        self._repr_type_str = "Coercion" if self._is_coercion else "Conversion"

    cpdef Element _call_(self, x):
        try:
            return self._codomain._element_constructor(self._codomain, x)
        except Exception:
            if print_warnings:
                print type(self._codomain), self._codomain
                print type(self._codomain._element_constructor), self._codomain._element_constructor
            raise

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        try:
            if len(args) == 0:
                if len(kwds) == 0:
                    # This line is apparently never used in any tests (hivert, 2009-04-28)
                    return self._codomain._element_constructor(self._codomain, x)
                else:
                    return self._codomain._element_constructor(self._codomain, x, **kwds)
            else:
                if len(kwds) == 0:
                    return self._codomain._element_constructor(self._codomain, x, *args)
                else:
                    return self._codomain._element_constructor(self._codomain, x, *args, **kwds)
        except Exception:
            if print_warnings:
                print type(self._codomain), self._codomain
                print type(self._codomain._element_constructor), self._codomain._element_constructor
            raise


cdef class DefaultConvertMap_unique(DefaultConvertMap):
    """
    This morphism simply defers action to the codomain's
    element_constructor method, WITHOUT passing in the codomain as the
    first argument.

    This is used for creating elements that don't take a parent as the
    first argument to their __init__ method, for example, Integers,
    Rationals, Algebraic Reals... all have a unique parent. It is also
    used when the element_constructor is a bound method (whose self
    argument is assumed to be bound to the codomain).
    """
    cpdef Element _call_(self, x):
        try:
            return self._codomain._element_constructor(x)
        except Exception:
            if print_warnings:
                print type(self._codomain), self._codomain
                print type(self._codomain._element_constructor), self._codomain._element_constructor
            raise

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        try:
            if len(args) == 0:
                if len(kwds) == 0:
                    return self._codomain._element_constructor(x)
                else:
                    return self._codomain._element_constructor(x, **kwds)
            else:
                if len(kwds) == 0:
                    return self._codomain._element_constructor(x, *args)
                else:
                    return self._codomain._element_constructor(x, *args, **kwds)
        except Exception:
            if print_warnings:
                print type(self._codomain), self._codomain
                print type(self._codomain._element_constructor), self._codomain._element_constructor
            raise


cdef class NamedConvertMap(Map):
    """
    This is used for creating a elements via the _xxx_ methods.

    For example, many elements implement an _integer_ method to
    convert to ZZ, or a _rational_ method to convert to QQ.
    """

    def __init__(self, domain, codomain, method_name, force_use=False):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_maps import NamedConvertMap
            sage: var('t')
            t
            sage: mor = NamedConvertMap(SR, QQ['t'], '_polynomial_')
            sage: mor(t^2/4+1)
            1/4*t^2 + 1
            sage: mor = NamedConvertMap(SR, GF(7)[['t']], '_polynomial_')
            sage: mor(t^2/4+1)
            1 + 2*t^2
        """
        if PY_TYPE_CHECK(domain, type):
            domain = Set_PythonType(domain)
        Map.__init__(self, domain, codomain)
        self._coerce_cost = 400
        self._force_use = force_use
        self.method_name = method_name
        self._repr_type_str = "Conversion via %s method" % self.method_name

    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_maps import NamedConvertMap
            sage: f = NamedConvertMap(GF(5), QQ, '_integer_'); f
            Conversion via _integer_ method map:
              From: Finite Field of size 5
              To:   Rational Field
            sage: f(19)
            4
            sage: f(19).parent()
            Rational Field
        """
        try:
            method = getattr(x, self.method_name)
        except AttributeError:
            if print_warnings:
                print type(x), x
                print type(self._codomain), self._codomain
                print self.method_name
            raise TypeError, "Cannot coerce %s to %s"%(x, self._codomain)
        cdef Map m
        cdef Element e = method(self._codomain)
        if e is None:
            raise RuntimeError, "BUG in coercion model: %s method of %s returned None" % (self.method_name, type(x))
        if e._parent is not self._codomain:
            m = self._codomain.convert_map_from(e._parent)
            if m is None or m is self:
                raise TypeError
            e = m._call_(e)
        return e

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_maps import NamedConvertMap
            sage: f = NamedConvertMap(SR, ZZ['x'], '_polynomial_')
            sage: f(x^2+1, check=True)
            x^2 + 1
        """
        return self._codomain._element_constructor(self._call_(x), *args, **kwds)


# Perhaps this could be a method, extracting (<PyMethodDescrObject *>(<object>Parent).coerce_map_from).d_method.ml_meth and/or PyCFunction_GET_FUNCTION(method)
# and constructing a CCallableConvertMap_class if it is bound to the codomain.

cdef class CallableConvertMap(Map):
    cdef bint _parent_as_first_arg
    cdef _func

    def __init__(self, domain, codomain, func, parent_as_first_arg=None):
        """
        This lets one easily create maps from any callable object.

        This is especially useful to create maps from bound methods.

        EXAMPLES::

            sage: from sage.structure.coerce_maps import CallableConvertMap
            sage: def foo(P, x): return x/2
            sage: f = CallableConvertMap(ZZ, QQ, foo)
            sage: f(3)
            3/2
            sage: f
            Conversion via foo map:
              From: Integer Ring
              To:   Rational Field

        Create a homomorphism from $\RR$ to $\RR^+$ viewed as additive groups.

        ::

            sage: f = CallableConvertMap(RR, RR, exp, parent_as_first_arg=False)
            sage: f(0)
            1.00000000000000
            sage: f(1)
            2.71828182845905
            sage: f(-3)
            0.0497870683678639
        """
        if PY_TYPE_CHECK(domain, type):
            domain = Set_PythonType(domain)
        Map.__init__(self, domain, codomain)
        self._coerce_cost = 100
        self._func = func
        if parent_as_first_arg is None:
            if PY_TYPE_CHECK(func, types.MethodType):
                # can't easily access self
                parent_as_first_arg = False
            elif PY_TYPE_CHECK(func, BuiltinMethodType):
                parent_as_first_arg = codomain is func.__self__
            else:
                parent_as_first_arg = True
        self._parent_as_first_arg = parent_as_first_arg
        try:
            self._repr_type_str = "Conversion via %s" % self._func.__name__
        except AttributeError:
            self._repr_type_str = "Conversion via %s" % self._func

    cpdef Element _call_(self, x):
        """
        Because self._func may be anything we do a little bit of sanity
        checking (the return value must be an element with the correct parent).

        TESTS::

            sage: from sage.structure.coerce_maps import CallableConvertMap
            sage: def foo(P, x): return x
            sage: f = CallableConvertMap(ZZ, ZZ, foo)
            sage: f(0)
            0
            sage: f = CallableConvertMap(ZZ, QQ, foo)
            sage: f(0)
            Traceback (most recent call last):
            ...
            RuntimeError: BUG in coercion model: <function foo at ...> returned element with wrong parent (expected Rational Field got Integer Ring)
            sage: def foo(P, x): return None
            sage: f = CallableConvertMap(ZZ, QQ, foo)
            sage: f(0)
            Traceback (most recent call last):
            ...
            RuntimeError: BUG in coercion model: <function foo at ...> returned None
        """
        cdef Element y
        try:
            if self._parent_as_first_arg:
                y = self._func(self._codomain, x)
            else:
                y = self._func(x)
        except Exception:
            if print_warnings:
                print self._func
                print self._codomain
            raise
        if y is None:
            raise RuntimeError, "BUG in coercion model: %s returned None" % (self._func)
        elif y._parent is not self._codomain:
            raise RuntimeError, "BUG in coercion model: %s returned element with wrong parent (expected %s got %s)" % (self._func, self._codomain, y._parent)
        return y

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        TESTS::

            sage: from sage.structure.coerce_maps import CallableConvertMap
            sage: def foo(P, x, y): return x or y
            sage: f = CallableConvertMap(ZZ, ZZ, foo)
            sage: f(0, 3)
            3
            sage: f = CallableConvertMap(ZZ, QQ, foo)
            sage: f(0, 3)
            Traceback (most recent call last):
            ...
            RuntimeError: BUG in coercion model: <function foo at ...> returned element with wrong parent (expected Rational Field got Integer Ring)
            sage: f(None, None)
            Traceback (most recent call last):
            ...
            RuntimeError: BUG in coercion model: <function foo at ...> returned None
        """
        cdef Element y
        try:
            if self._parent_as_first_arg:
                y = self._func(self._codomain, x, *args, **kwds)
            else:
                y = self._func(x, *args, **kwds)
        except Exception:
            if print_warnings:
                print self._func
                print self._codomain
            raise
        if y is None:
            raise RuntimeError, "BUG in coercion model: %s returned None" % (self._func)
        elif y._parent is not self._codomain:
            raise RuntimeError, "BUG in coercion model: %s returned element with wrong parent (expected %s got %s)" % (self._func, self._codomain, y._parent)
        return y


cdef class CCallableConvertMap_class(Map):
    cdef Element (*_func)(Parent, object)
    cdef public _name

    def __init__(self, domain, codomain, name):
        if PY_TYPE_CHECK(domain, type):
            domain = Set_PythonType(domain)
        Map.__init__(self, domain, codomain)
        self._coerce_cost = 10
        self._name = name

    cpdef Element _call_(self, x):
        """
        TESTS::

            sage: from sage.structure.coerce_maps import test_CCallableConvertMap
            sage: f = test_CCallableConvertMap(QQ, 'test')
            sage: f(1/3)
            -8/27
        """
        return self._func(self._codomain, x)

    def _repr_type(self):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_maps import test_CCallableConvertMap
            sage: test_CCallableConvertMap(ZZ, 'any name')
            Conversion via c call 'any name' map:
              From: Integer Ring
              To:   Integer Ring
            sage: test_CCallableConvertMap(ZZ, None)  # random address
            Conversion via c call at 0xc339000 map:
              From: Integer Ring
              To:   Integer Ring
        """
        if self._name is None:
            return "Conversion via c call at 0x%x" % <long>self._func
        else:
            return "Conversion via c call '%s'" % self._name


cdef Map CCallableConvertMap(domain, codomain, void* func, name):
    """
    Use this to create a map from domain to codomain by calling func
    (which must be a function pointer taking a Parent and object, and
    returning an Element in the given Parent).

    This is the c analogue of CallableConvertMap.
    """
    # Cython doesn't yet accept function pointers as arguments,
    # change this when it does.
    cdef CCallableConvertMap_class map = CCallableConvertMap_class(domain, codomain, name)
    map._func = <Element (*)(Parent, object)>func
    return map

cpdef Element _ccall_test_function(codomain, x):
    """
    For testing CCallableConvertMap_class. Returns x*x*x-x in the codomain.

    TESTS::

        sage: from sage.structure.coerce_maps import _ccall_test_function
        sage: _ccall_test_function(ZZ, 1)
        0
        sage: _ccall_test_function(ZZ, 2)
        6
        sage: _ccall_test_function(ZZ, -3)
        -24
    """
    return codomain(x*x*x-x)

def test_CCallableConvertMap(domain, name=None):
    """
    For testing CCallableConvertMap_class.

    TESTS::

        sage: from sage.structure.coerce_maps import test_CCallableConvertMap
        sage: f = test_CCallableConvertMap(ZZ, 'test'); f
        Conversion via c call 'test' map:
          From: Integer Ring
          To:   Integer Ring
        sage: f(3)
        24
        sage: f(9)
        720
    """
    return CCallableConvertMap(domain, domain, <void*>&_ccall_test_function, name)


cdef class ListMorphism(Map):

    cdef Map _real_morphism

    def __init__(self, domain, Map real_morphism):
        if not PY_TYPE_CHECK(domain, Parent):
            domain = Set_PythonType(domain)
        Map.__init__(self, domain, real_morphism.codomain())
        self._coerce_cost = real_morphism._coerce_cost + 3
        self._real_morphism = real_morphism
        self._repr_type_str = "List"

    cpdef Element _call_(self, x):
        try:
            x = x._data
        except AttributeError:
            x = list(x)
        return self._real_morphism._call_(x)

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        try:
            x = x._data
        except AttributeError:
            x = list(x)
        return self._real_morphism._call_with_args(x, args, kwds)


cdef class TryMap(Map):
    def __init__(self, Map morphism_preferred, Map morphism_backup, error_types=None):
        """
        TESTS::

            sage: sage.structure.coerce_maps.TryMap(RDF.coerce_map_from(QQ), RDF.coerce_map_from(ZZ))
            Traceback (most recent call last):
            ...
            TypeError: incorrectly matching parent
        """
        if (morphism_preferred.domain() is not morphism_backup.domain()
             or morphism_preferred.codomain() is not morphism_backup.codomain()):
            raise TypeError, "incorrectly matching parent"
        Map.__init__(self, morphism_preferred.parent())
        self._map_p = morphism_preferred
        self._map_b = morphism_backup
        if error_types is None:
            self._error_types = (ValueError, TypeError, AttributeError)
        else:
            self._error_types = error_types

    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: map1 = sage.structure.coerce_maps.CallableConvertMap(ZZ, QQ, lambda parent, x: 1/x)
            sage: map2 = QQ.coerce_map_from(ZZ)
            sage: map = sage.structure.coerce_maps.TryMap(map1, map2, error_types=(ZeroDivisionError,))
            sage: map(3)
            1/3
            sage: map(-7)
            -1/7
            sage: map(0)
            0
        """
        try:
            return self._map_p._call_(x)
        except self._error_types:
            return self._map_b._call_(x)

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        EXAMPLES::

            sage: map1 = sage.structure.coerce_maps.CallableConvertMap(ZZ, QQ, lambda parent, x, y:  y/x)
            sage: map2 = sage.structure.coerce_maps.CallableConvertMap(ZZ, QQ, lambda parent, x, y: 23/1)
            sage: map = sage.structure.coerce_maps.TryMap(map1, map2, error_types=(ZeroDivisionError,))
            sage: map._call_with_args(3, (2,))
            2/3
            sage: map._call_with_args(-7, (5,))
            -5/7
            sage: map._call_with_args(0, (1,))
            23
        """
        try:
            return self._map_p._call_with_args(x, args, kwds)
        except self._error_types:
            return self._map_b._call_with_args(x, args, kwds)
