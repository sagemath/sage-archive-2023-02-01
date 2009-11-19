"""
This file contains a discussion, examples, and tests about nested
classes and parents. It is kept in a separate file to avoid import
loops.

EXAMPLES:

Currently pickling fails for parents using nested classes (typically
for categories), but deriving only from Parent::

    sage: from sage.misc.nested_class_test import TestParent1, TestParent2, TestParent3, TestParent4
    sage: P = TestParent1()
    sage: TestSuite(P).run()
    Traceback (most recent call last):
    ...
    PicklingError: Can't pickle <class 'sage.misc.nested_class_test.Element'>: attribute lookup sage.misc.nested_class_test.Element failed

They actually need to be in the NestedClassMetaclass. However, due to
a technical detail, this is currently not directly supported::

    sage: P = TestParent2()
    Traceback (most recent call last):
    ...
    TypeError: metaclass conflict: the metaclass of a derived class must be a (non-strict) subclass of the metaclasses of all its bases
    sage: TestSuite(P).run()  # todo: not implemented

Instead, the easiest is to inherit from UniqueRepresentation, which is
what you want to do anyway most of the time::

    sage: P = TestParent3()
    sage: TestSuite(P).run()

This is what all Sage's parents using categories currently do. An
alternative is to use ClasscallMetaclass as metaclass::

    sage: P = TestParent4()
    sage: TestSuite(P).run()

"""
#*****************************************************************************
#  Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.nested_class import NestedClassMetaclass

class TestParent1(Parent):
    def __init__(self):
        """
        EXAMPLES::

            sage: sage.misc.nested_class_test.TestParent1()
            <class 'sage.misc.nested_class_test.TestParent1_with_category'>
        """
        from sage.categories.all import Sets
        Parent.__init__(self, category = Sets())

    class Element(ElementWrapper):
        pass

class TestParent2(Parent):
    __metaclass__ = NestedClassMetaclass

    def __init__(self):
        """
        EXAMPLES::

            sage: sage.misc.nested_class_test.TestParent2()
            Traceback (most recent call last):
            TypeError: metaclass conflict: the metaclass of a derived class must be a (non-strict) subclass of the metaclasses of all its bases
        """
        from sage.categories.all import Sets
        Parent.__init__(self, category = Sets())

    class Element(ElementWrapper):
        pass

class TestParent3(UniqueRepresentation, Parent):

    def __init__(self):
        """
        EXAMPLES::

            sage: sage.misc.nested_class_test.TestParent3()
            <class 'sage.misc.nested_class_test.TestParent3_with_category'>
        """
        from sage.categories.all import Sets
        Parent.__init__(self, category = Sets())

    class Element(ElementWrapper):
        pass

class TestParent4(Parent):
    __metaclass__ = ClasscallMetaclass

    def __init__(self):
        """
        EXAMPLES::

            sage: sage.misc.nested_class_test.TestParent4()
            <class 'sage.misc.nested_class_test.TestParent4_with_category'>
        """
        from sage.categories.all import Sets
        Parent.__init__(self, category = Sets())

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.nested_class_test import TestParent4
            sage: TestParent4() == TestParent4()
            True
        """
        return self.__class__ == other.__class__

    class Element(ElementWrapper):
        pass
