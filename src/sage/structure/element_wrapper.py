"""
A class for wrapping Sage or Python objects as Sage elements
"""

from sage.structure.element import Element

class ElementWrapper(Element):
    r"""
    A class for wrapping Sage or Python objects as Sage elements

    EXAMPLES::

        sage: o = ElementWrapper("bla", parent=ZZ); o
        'bla'
        sage: isinstance(o, sage.structure.element.Element)
        True
        sage: o.parent()
        Integer Ring
        sage: o.value
        'bla'

    This is of course a meaningless example: `bla` is not an element
    of ``ZZ``.

    Note that ``o`` is not *an instance of* ``str``, but rather
    *contains a* ``str``. Therefore, ``o`` does not inherit the string
    methods. On the other hand, it is provided with reasonable default
    implementations for equality testing, hashing, etc.

    The typical use case of `ElementWrapper` is for trivially
    constructing new element classes from preexisting Sage or Python
    classes, with a containment relation. Here we construct the
    tropical monoid of integers endowed with ``min`` as
    multiplication. There, it is desirable *not* to inherit the
    ``factor`` method from Integer::

        sage: class MinMonoid(Parent):
        ...       def _repr_(self):
        ...           return "The min monoid"
        ...
        sage: M = MinMonoid()
        sage: class MinMonoidElement(ElementWrapper):
        ...       wrapped_class = Integer
        ...
        ...       def __mul__(self, other):
        ...           return MinMonoidElement(min(self.value, other.value), parent = self.parent())
        sage: x = MinMonoidElement(5, parent = M); x
        5
        sage: x.parent()
        The min monoid
        sage: x.value
        5
        sage: y = MinMonoidElement(3, parent = M)
        sage: x * y
        3

    This example was voluntarily kept to a bare minimum. See the
    (upcoming) examples in the categories for several full featured
    applications.

    Caveat: the order between the value and the parent argument is
    likely to change shortly. At this point, all the code using it in
    the Sage library will be updated. There will be no transition period.
    """

    wrapped_class = object

    def __init__(self, value, parent):
        """
        EXAMPLES::

            sage: a = ElementWrapper(1, parent = ZZ)
            sage: a == loads(dumps(a))
            True
        """
        assert isinstance(value, self.wrapped_class)
        Element.__init__(self, parent = parent)
        self.value = value

    def __repr__(self):
        """
        EXAMPLES::

            sage: ElementWrapper(1, parent = ZZ)
            1
        """
        return repr(self.value)

    def __hash__(self):
        """
        Returns the same hash as for the wrapped element

        EXAMPLES::

            sage: hash(ElementWrapper(1, parent = ZZ))
            1
            sage: hash(ElementWrapper(1, parent = QQ))
            1

        TODO: should this take the parent and/or the class into account?
        """
        return hash(self.value)

    def __eq__(self, other):
        """
        Default implementation of equality testing: two elements are
        equal if they have the same value, same parent, and same class.

        EXAMPLES::

            sage: parent1 = ZZ
            sage: parent2 = QQ
            sage: l11 = ElementWrapper(1, parent = parent1)
            sage: l12 = ElementWrapper(2, parent = parent1)
            sage: l21 = ElementWrapper(1, parent = parent2)
            sage: l22 = ElementWrapper(2, parent = parent2)
            sage: l11 == l11
            True
            sage: l11 == l12
            False
            sage: l11 == l21
            False
        """
        return (self.__class__ is other.__class__ and
                self.parent() == other.parent() and
                self.value == other.value)

    def __cmp__(self, other):
        """
        Default implementation of comparison: compare first values,
        then parents, then class.

        EXAMPLES::

            sage: parent1 = ZZ
            sage: parent2 = QQ
            sage: l11 = ElementWrapper(1, parent = parent1)
            sage: l12 = ElementWrapper(2, parent = parent1)
            sage: l21 = ElementWrapper(1, parent = parent2)
            sage: l22 = ElementWrapper(2, parent = parent2)
            sage: cmp(l11, l11)
            0
            sage: cmp(l11, l12), cmp(l12, l11)   # values differ
            (-1, 1)
            sage: cmp(l11, l21) in [-1, 1]       # parents differ
            True
            sage: cmp(l21, l11) == -cmp(l11, l21)
            True
            sage: cmp(l11, 1) in [-1,1]          # class differ
            True
        """
        if self.__class__ != other.__class__:
            return cmp(self.__class__, other.__class__)
        if self.parent() != other.parent():
            return cmp(self.parent(), other.parent())
        return cmp(self.value, other.value)
