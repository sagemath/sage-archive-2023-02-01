r"""
Examples of semigroups in cython
"""
from sage.structure.parent import Parent
from sage.structure.element cimport Element
from sage.categories.all import Category, Semigroups
from sage.misc.cachefunc import cached_method
from sage.categories.examples.semigroups import LeftZeroSemigroup as LeftZeroSemigroupPython

class DummyClass:
    def method(self):
        """
        TESTS::

            sage: from sage.categories.examples.semigroups_cython import DummyClass
            sage: DummyClass().method()
        """
        pass

cdef class DummyCClass:
    def method(self):
        """
        TESTS::

            sage: from sage.categories.examples.semigroups_cython import DummyCClass
            sage: DummyCClass().method()
        """
        pass

    cpdef cpmethod(self):
        """
        TESTS::

            sage: from sage.categories.examples.semigroups_cython import DummyCClass
            sage: DummyCClass().cpmethod()
        """
        pass

instancemethod     = type(DummyClass.method)
method_descriptor  = type(DummyCClass.method)

cdef class IdempotentSemigroupsElement(Element):
    def _pow_(self, i):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S(2)._pow_(3) # todo: not implemented (binding; see __getattr__)
            2
        """
        assert i > 0
        return self

    cpdef is_idempotent_cpdef(self):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S(2).is_idempotent_cpdef()  # todo: not implemented (binding; see __getattr__)
            True
        """
        return True

class IdempotentSemigroups(Category):
    #@cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import IdempotentSemigroups
            sage: IdempotentSemigroups().super_categories()
            [Category of semigroups]
        """
        return [Semigroups()]

    ElementMethods = IdempotentSemigroupsElement


cdef class LeftZeroSemigroupElement(Element):
    cdef object _value

    def __init__(self, value, parent):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: x = S(3)
            sage: TestSuite(x).run()
        """
        Element.__init__(self, parent = parent)
        self._value = value

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S(3)                    # indirect doctest
            3
        """
        return repr(self._value)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: x = S(3)
            sage: x.__reduce__()
            (<type 'sage.categories.examples.semigroups_cython.LeftZeroSemigroupElement'>,
             (3, An example of a semigroup: the left zero semigroup))
        """
        return LeftZeroSemigroupElement, (self._value, self._parent)

    def __cmp__(LeftZeroSemigroupElement self, LeftZeroSemigroupElement other):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S(3) == S(3)
            True
            sage: S(3) == S(2)
            False
        """
        return cmp(self._value, other._value)

    def __getattr__(self, name):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: x = S(3)
            sage: x.is_idempotent
            <bound method LeftZeroSemigroupElement.is_idempotent of 3>
            sage: x.is_idempotent()
            True
        """
        result = getattr(self.parent().category().element_class, name)
        if isinstance(result, instancemethod):
            return instancemethod(result.im_func, self, self.__class__)
        elif isinstance(result,  method_descriptor):
            return result # should bind the method descriptor to appropriate object
        else:
            return result

    # Apparently, python looks for __mul__, __pow__, ... in the
    # class of self rather than in self itself. No big deal, since
    # those will usually be defined in a cython super class of
    # this class
    def __mul__(self, other):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S(2) * S(3)
            2
        """
        return self._mul_(other)

    cpdef _mul_(self, other):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S(2)._mul_(S(3))
            2
        """
        return self.parent().product(self, other)

    def __pow__(self, i, dummy):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S(2)^3  # todo: not implemented (binding; see __getattr__)
        """
        return self._pow_(i)

class LeftZeroSemigroup(LeftZeroSemigroupPython):
    r"""
    An example of semigroup

    This class illustrates a minimal implementation of a semi-group
    where the element class is an extension type, and still gets code
    from the category. Also, the category itself includes some cython
    methods.

    This is purely a proof of concept. The code obviously needs refactorisation!

    Comments:
    - nested classes seem not to be currently supported by cython
    - one cannot play ugly use class surgery tricks (as with _mul_parent)
      available operations should really be declared to the coercion model!

    EXAMPLES::

        sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
        sage: S = LeftZeroSemigroup(); S
        An example of a semigroup: the left zero semigroup

    This is the semigroup which contains all sort of objects:

        sage: S.some_elements()
        [3, 42, 'a', 3.3999999999999999, 'raton laveur']

    with product rule is given by $a \times b = a$ for all $a,b$.

        sage: S('hello') * S('world')
        'hello'

        sage: S(3)*S(1)*S(2)
        3

        sage: S(3)^12312321312321         # todo: not implemented (see __getattr__)
        3

        sage: TestSuite(S).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_element_pickling() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass

        # That's really the only method which is obtained from the category ...
        sage: S(42).is_idempotent
        <bound method LeftZeroSemigroupElement.is_idempotent of 42>
        sage: S(42).is_idempotent()
        True

        sage: S(42)._pow_                 # how to bind it?
        <method '_pow_' of 'sage.categories.examples.semigroups_cython.IdempotentSemigroupsElement' objects>
        sage: S(42)^10                    # todo: not implemented (see __getattr__)
        42

        sage: S(42).is_idempotent_cpdef   # how to bind it?
        <method 'is_idempotent_cpdef' of 'sage.categories.examples.semigroups_cython.IdempotentSemigroupsElement' objects>
        sage: S(42).is_idempotent_cpdef() # todo: not implemented (see __getattr__)
        True
    """

    def __init__(self):
        """
        TESTS::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S.category()
            Category of idempotent semigroups
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category = IdempotentSemigroups())

    Element = LeftZeroSemigroupElement
