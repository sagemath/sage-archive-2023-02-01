include "../ext/stdsage.pxi"

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

    cpdef Element _call_(self, x):
        try:
            return self._codomain._element_constructor(self._codomain, x)
        except:
            if print_warnings:
                print type(self._codomain), self._codomain
                print type(self._codomain._element_constructor), self._codomain._element_constructor
            raise

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        try:
            if len(args) == 0:
                if len(kwds) == 0:
                    return self._codomain.__element_constructor(self._codomain, x)
                else:
                    return self._codomain._element_constructor(self._codomain, x, **kwds)
            else:
                if len(kwds) == 0:
                    return self._codomain._element_constructor(self._codomain, x, *args)
                else:
                    return self._codomain._element_constructor(self._codomain, x, *args, **kwds)
        except:
            if print_warnings:
                print type(self._codomain), self._codomain
                print type(self._codomain._element_constructor), self._codomain._element_constructor
            raise

    def _repr_type(self):
        return "Conversion"

cdef class DefaultConvertMap_unique(DefaultConvertMap):
    """
    This morphism simply differs action to the codomain's element_constructor method,
    WITHOUT passing in the codomain as the first argument.

    This is used for creating elements that don't take a parent as the first argument
    to their __init__ method, for example, Integers, Rationals, Algebraic Reals... all
    have a unique parent. It is also used when the element_constructor is a bound method
    (whose self argument is assumed to be bound to the codomain).
    """
    cpdef Element _call_(self, x):
        try:
            return self._codomain._element_constructor(x)
        except:
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
        except:
            if print_warnings:
                print type(self._codomain), self._codomain
                print type(self._codomain._element_constructor), self._codomain._element_constructor
            raise

    def _repr_type(self):
        return "Coercion" if self._is_coercion else "Conversion"

cdef class NamedConvertMap(Map):
    """
    This is used for creating a elements via the _xxx_ methods.

    For example, many elements implement an _integer_ method to
    convert to ZZ, or a _rational_ method to convert to QQ.
    """

    def __init__(self, domain, codomain, method_name, force_use=False):
        if PY_TYPE_CHECK(domain, type):
            domain = Set_PythonType(domain)
        Map.__init__(self, domain, codomain)
        self._coerce_cost = 400
        self._force_use = force_use
        self.method_name = method_name

    cpdef Element _call_(self, x):
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

    def _repr_type(self):
        return "Conversion via %s" % self.method_name

cdef class CallableConvertMap(Map):
    cdef bint _parent_as_first_arg
    cdef _func

    def __init__(self, domain, codomain, func, parent_as_first_arg=None):
        """
        This lets one easily create maps from any callable object.

        This is especially useful to create maps from bound methods.

        EXAMPLES:
            sage: from sage.structure.coerce_maps import CallableConvertMap
            sage: def foo(P, x): return x/2
            sage: f = CallableConvertMap(ZZ, QQ, foo)
            sage: f(3)
            3/2
            sage: f
            Conversion via foo map:
              From: Integer Ring
              To:   Rational Field

        Create a homomorphism from $\R$ to $\R^+$ viewed as additave groups.
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

    cpdef Element _call_(self, x):
        """
        Because self._func may be anything we do a little bit of sanity
        checking (the return value must be an element with the correct parent).

        TESTS:
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
            sage: f(None)
            Traceback (most recent call last):
            ...
            RuntimeError: BUG in coercion model: <function foo at ...> returned None
        """
        cdef Element y
        if self._parent_as_first_arg:
            y = self._func(self._codomain, x)
        else:
            y = self._func(x)
        if y is None:
            raise RuntimeError, "BUG in coercion model: %s returned None" % (self._func)
        elif y._parent is not self._codomain:
            raise RuntimeError, "BUG in coercion model: %s returned element with wrong parent (expected %s got %s)" % (self._func, self._codomain, y._parent)
        return y

    def _repr_type(self):
        try:
            return "Conversion via %s" % self._func.__name__
        except AttributeError:
            return "Conversion via %s" % self._func

cdef class ListMorphism(Map):

    cdef Map _real_morphism

    def __init__(self, domain, Map real_morphism):
        if not PY_TYPE_CHECK(domain, Parent):
            domain = Set_PythonType(domain)
        Map.__init__(self, domain, real_morphism.codomain())
        self._coerce_cost = real_morphism._coerce_cost + 3
        self._real_morphism = real_morphism

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

    def _repr_type(self):
        return "List"

cdef class TryMap(Map):
    def __init__(self, morphism_preferred, morphism_backup, error_types=None):
        if morphism_preferred.parent() is not morphism_backup.parent():
            raise TypeError, "incorrectly matching parent"
        Map.__init__(self, morphism_preferred.parent())
        self._map_p = morphism_preferred
        self._map_b = morphism_backup
        if error_types is None:
            self._error_types = (ValueError, TypeError, AttributeError)
        else:
            self._error_types = error_types

    cpdef Element _call_(self, x):
        try:
            return self._map_p._call_(x)
        except self._error_types:
            return self._map_b._call_(x)

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        try:
            return self._map_p._call_with_args(x, args, kwds)
        except self._error_types:
            return self._map_b._call_with_args(x, args, kwds)
