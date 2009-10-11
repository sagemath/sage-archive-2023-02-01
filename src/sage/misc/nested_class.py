"""
As of Python 2.6, names for nested classes are set by Python in  a
way which is incompatible with the pickling of such classes (pickling by name)::

    sage: class A:
    ...       class B:
    ...            pass
    sage: A.B.__name__
    'B'

instead of the a priori more natural `A.B`.

Furthermore, upon pickling (here in save_global) *and* unpickling (in
load_global) a class with name "A.B" in a module mod, the standard
cPickle module searches for "A.B" in mod.__dict__ instead of looking
up "A" and then "B" in the result.

See: http://groups.google.com/group/sage-devel/browse_thread/thread/6c7055f4a580b7ae/

This module provides two utilities to workaround this issue:

 - :func:`nested_pickle` "fixes" recursively the name of the
   subclasses of a class and inserts their fullname 'A.B' in
   mod.__dict__

 - :class:`NestedClassMetaclass` is a metaclass ensuring that nested_pickle is
   called on a class upon creation.


EXAMPLES::

    sage: from sage.misc.nested_class import A1, nested_pickle

    sage: A1.A2.A3.__name__
    'A3'
    sage: A1.A2.A3
    <class sage.misc.nested_class.A3 at ...>

    sage: nested_pickle(A1)
    <class sage.misc.nested_class.A1 at ...>

    sage: A1.A2
    <class sage.misc.nested_class.A1.A2 at ...>

    sage: A1.A2.A3
    <class sage.misc.nested_class.A1.A2.A3 at ...>
    sage: A1.A2.A3.__name__
    'A1.A2.A3'

    sage: sage.misc.nested_class.__dict__['A1.A2'] is A1.A2
    True
    sage: sage.misc.nested_class.__dict__['A1.A2.A3'] is A1.A2.A3
    True

All of this is not perfect. In the following scenario::

    sage: class A1:
    ...       class A2:
    ...           pass
    sage: class B1:
    ...       A2 = A1.A2
    ...

The name for A1.A2 could potentially be set to B1.A2. But that will work anyway.
"""

import sys

def modify_for_nested_pickle(cls, name_prefix, module):
    r"""
    Modify the subclasses of the given class to be picklable, by
    giving them a mangled name and putting the mangled name in the
    module namespace.

    INPUTS:

    - ``cls`` - The class to modify.

    - ``name_prefix`` - The prefix to prepend to the class name.

    - ``module`` - The module object to modify with the mangled name.

    EXAMPLES::

        sage: from sage.misc.nested_class import *
        sage: class A(object):
        ...       class B(object):
        ...           pass
        ...
        sage: module = sys.modules['__main__']
        sage: A.B.__name__
        'B'
        sage: getattr(module, 'A.B', 'Not found')
        'Not found'
        sage: modify_for_nested_pickle(A, 'A', sys.modules['__main__'])
        sage: A.B.__name__
        'A.B'
        sage: getattr(module, 'A.B', 'Not found')
        <class '__main__.A.B'>

    """
    import types
    for (name, v) in cls.__dict__.iteritems():
        if isinstance(v, (type, types.ClassType)):
            if v.__name__ == name and v.__module__ == module.__name__ and getattr(module, name, None) is not v:
                # OK, probably this is a nested class.
                dotted_name = name_prefix + '.' + name
                v.__name__ = dotted_name
                setattr(module, dotted_name, v)
                modify_for_nested_pickle(v, dotted_name, module)

def nested_pickle(cls):
    r"""
    This decorator takes a class that potentially contains nested classes.
    For each such nested class, its name is modified to a new illegal
    identifier, and that name is set in the module.  For example, if
    you have::

        sage: from sage.misc.nested_class import nested_pickle
        sage: module = sys.modules['__main__']
        sage: class A(object):
        ...       class B:
        ...           pass
        sage: nested_pickle(A)
        <class '__main__.A'>

    then the name of class B will be modified to 'A.B', and the 'A.B'
    attribute of the module will be set to class B::

        sage: A.B.__name__
        'A.B'
        sage: getattr(module, 'A.B', 'Not found')
        <class __main__.A.B at ...>

    In Python 2.6, decorators work with classes; then @nested_pickle
    should work as a decorator::

        sage: @nested_pickle    # todo: not implemented
        ...   class A2(object):
        ...       class B:
        ...           pass
        sage: A2.B.__name__    # todo: not implemented
        'A2.B'
        sage: getattr(module, 'A2.B', 'Not found')    # todo: not implemented
        <class __main__.A2.B at ...>

    EXAMPLES::

        sage: from sage.misc.nested_class import *
        sage: loads(dumps(MainClass.NestedClass())) # indirect doctest
        <sage.misc.nested_class.MainClass.NestedClass object at 0x...>
    """
    modify_for_nested_pickle(cls, cls.__name__, sys.modules[cls.__module__])
    return cls


class NestedClassMetaclass(type):
    r"""
    A metaclass for nested pickling.

    Check that one can use a metaclass to ensure nested_pickle
    is called on any derived subclass::

        sage: from sage.misc.nested_class import NestedClassMetaclass
        sage: class ASuperClass(object):
        ...       __metaclass__ = NestedClassMetaclass
        ...
        sage: class A3(ASuperClass):
        ...       class B(object):
        ...           pass
        ...
        sage: A3.B.__name__
        'A3.B'
        sage: getattr(sys.modules['__main__'], 'A3.B', 'Not found')
        <class '__main__.A3.B'>
    """
    def __init__(self, *args):
        r"""
        This invokes the nested_pickle on construction.

        sage: from sage.misc.nested_class import NestedClassMetaclass
        sage: class A(object):
        ...       __metaclass__ = NestedClassMetaclass
        ...       class B(object):
        ...           pass
        ...
        sage: A.B
        <class '__main__.A.B'>
        sage: getattr(sys.modules['__main__'], 'A.B', 'Not found')
        <class '__main__.A.B'>
        """
        nested_pickle(self)


class MainClass(object):
    r"""
    A simple class to test nested_pickle.

    EXAMPLES::

        sage: from sage.misc.nested_class import *
        sage: loads(dumps(MainClass()))
        <sage.misc.nested_class.MainClass object at 0x...>
    """

    __metaclass__ = NestedClassMetaclass

    class NestedClass(object):
        r"""
        EXAMPLES::

            sage: from sage.misc.nested_class import *
            sage: loads(dumps(MainClass.NestedClass()))
            <sage.misc.nested_class.MainClass.NestedClass object at 0x...>
        """

        class NestedSubClass(object):
            r"""
            EXAMPLES::

                sage: from sage.misc.nested_class import *
                sage: loads(dumps(MainClass.NestedClass.NestedSubClass()))
                <sage.misc.nested_class.MainClass.NestedClass.NestedSubClass object at 0x...>
                sage: getattr(sage.misc.nested_class, 'MainClass.NestedClass.NestedSubClass')
                <class 'sage.misc.nested_class.MainClass.NestedClass.NestedSubClass'>
                sage: MainClass.NestedClass.NestedSubClass.__name__
                'MainClass.NestedClass.NestedSubClass'
            """
            pass

class SubClass(MainClass):
    r"""
    A simple class to test nested_pickle.

    EXAMPLES::

        sage: from sage.misc.nested_class import *
        sage: loads(dumps(SubClass.NestedClass()))
        <sage.misc.nested_class.MainClass.NestedClass object at 0x...>
    """
    pass

nested_pickle(SubClass)

class CopiedClass(object):
    r"""
    A simple class to test nested_pickle.

    EXAMPLES::

        sage: from sage.misc.nested_class import *
        sage: loads(dumps(CopiedClass.NestedClass()))
        <sage.misc.nested_class.MainClass.NestedClass object at 0x...>
        sage: loads(dumps(CopiedClass.NestedSubClass()))
        <sage.misc.nested_class.MainClass.NestedClass.NestedSubClass object at 0x...>
        sage: loads(dumps(SubClass()))
        <sage.misc.nested_class.SubClass object at 0x...>
    """
    NestedClass = MainClass.NestedClass
    NestedSubClass = MainClass.NestedClass.NestedSubClass
    SubClass = SubClass

nested_pickle(CopiedClass)

# Further classes for recursive tests

class A1:
    class A2:
        class A3:
            pass
