"""
Fixing Pickle for Nested Classes

As of Python 2.6, names for nested classes are set by Python in  a
way which is incompatible with the pickling of such classes (pickling by name)::

    sage: class A:
    ...       class B:
    ...            pass
    sage: A.B.__name__
    'B'

instead of the a priori more natural ``"A.B"``.

Furthermore, upon pickling (here in save_global) *and* unpickling (in
load_global) a class with name ``"A.B"`` in a module ``mod``, the standard
cPickle module searches for ``"A.B"`` in ``mod.__dict__`` instead of looking
up ``"A"`` and then ``"B"`` in the result.

See: http://groups.google.com/group/sage-devel/browse_thread/thread/6c7055f4a580b7ae/

This module provides two utilities to workaround this issue:

 - :func:`nested_pickle` "fixes" recursively the name of the
   subclasses of a class and inserts their fullname ``"A.B"`` in
   ``mod.__dict__``

 - :class:`NestedClassMetaclass` is a metaclass ensuring that
   :func:`nested_pickle` is called on a class upon creation.

See also :mod:`sage.misc.nested_class_test`.

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

The name for ``"A1.A2"`` could potentially be set to ``"B1.A2"``. But that will work anyway.
"""

import sys
cdef dict sys_modules = sys.modules

import types
from types import ClassType

__all__ = ['modify_for_nested_pickle', 'nested_pickle',
           'NestedClassMetaclass', 'MainClass'
           # Comment out to silence Sphinx warning about nested classes.
           #, 'SubClass', 'CopiedClass', 'A1'
           ]

cpdef modify_for_nested_pickle(cls, str name_prefix, module, first_run=True):
    r"""
    Modify the subclasses of the given class to be picklable, by
    giving them a mangled name and putting the mangled name in the
    module namespace.

    INPUTS:

    - ``cls`` - The class to modify.
    - ``name_prefix`` - The prefix to prepend to the class name.
    - ``module`` - The module object to modify with the mangled name.
    - ``first_run`` - optional bool (default True): Whether or not
      this function is run for the first time on ``cls``.

    NOTE:

    This function would usually not be directly called. It is internally used
    in :class:`NestedClassMetaclass`.

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

    Here we demonstrate the effect of the ``first_run`` argument::

        sage: modify_for_nested_pickle(A, 'X', sys.modules['__main__'])
        sage: A.B.__name__ # nothing changed
        'A.B'
        sage: modify_for_nested_pickle(A, 'X', sys.modules['__main__'], first_run=False)
        sage: A.B.__name__
        'X.A.B'

    Note that the class is now found in the module under both its old and
    its new name::

        sage: getattr(module, 'A.B', 'Not found')
        <class '__main__.X.A.B'>
        sage: getattr(module, 'X.A.B', 'Not found')
        <class '__main__.X.A.B'>


    TESTS:

    The following is a real life example, that was enabled by the internal
    use of the``first_run`` in :trac:`9107`::

        sage: cython_code = [
        ....:  "from sage.structure.unique_representation import UniqueRepresentation",
        ....:  "class A1(UniqueRepresentation):",
        ....:  "    class B1(UniqueRepresentation):",
        ....:  "        class C1: pass",
        ....:  "    class B2:",
        ....:  "        class C2: pass"]
        sage: import os
        sage: cython(os.linesep.join(cython_code))

    Before :trac:`9107`, the name of ``A1.B1.C1`` would have been wrong::

        sage: A1.B1.C1.__name__
        'A1.B1.C1'
        sage: A1.B2.C2.__name__
        'A1.B2.C2'
        sage: A_module = sys.modules[A1.__module__]
        sage: getattr(A_module, 'A1.B1.C1', 'Not found').__name__
        'A1.B1.C1'
        sage: getattr(A_module, 'A1.B2.C2', 'Not found').__name__
        'A1.B2.C2'

    """
    cdef str name, dotted_name
    cdef str mod_name = module.__name__
    cdef str cls_name = cls.__name__+'.'
    cdef str v_name
    if first_run:
        for (name, v) in cls.__dict__.iteritems():
            if isinstance(v, NestedClassMetaclass):
                v_name = v.__name__
                if v_name==name and v.__module__ == mod_name and getattr(module, v_name, None) is not v:
                    # OK, probably this is a nested class.
                    dotted_name = name_prefix + '.' + v_name
                    setattr(module, dotted_name, v)
                    modify_for_nested_pickle(v, name_prefix, module, False)
                    v.__name__ = dotted_name
            elif isinstance(v, (type, ClassType)):
                v_name = v.__name__
                if v_name==name and v.__module__ == mod_name and getattr(module, v_name, None) is not v:
                    # OK, probably this is a nested class.
                    dotted_name = name_prefix + '.' + v_name
                    setattr(module, dotted_name, v)
                    modify_for_nested_pickle(v, dotted_name, module)
                    v.__name__ = dotted_name
    else:
        for (name, v) in cls.__dict__.iteritems():
            if isinstance(v, (type, ClassType, NestedClassMetaclass)):
                v_name = v.__name__
                if v_name==cls_name+name and v.__module__ == mod_name:
                    # OK, probably this is a nested class.
                    dotted_name = name_prefix + '.' + v_name
                    setattr(module, dotted_name, v)
                    modify_for_nested_pickle(v, name_prefix, module, False)
                    v.__name__ = dotted_name
        

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

    then the name of class ``"B"`` will be modified to ``"A.B"``, and the ``"A.B"``
    attribute of the module will be set to class ``"B"``::

        sage: A.B.__name__
        'A.B'
        sage: getattr(module, 'A.B', 'Not found')
        <class __main__.A.B at ...>

    In Python 2.6, decorators work with classes; then ``@nested_pickle``
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
    modify_for_nested_pickle(cls, cls.__name__, sys_modules[cls.__module__])
    return cls

cdef class NestedClassMetaclass(type):
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
        modify_for_nested_pickle(self, self.__name__, sys_modules[self.__module__])

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

        sage: from sage.misc.nested_class import SubClass
        sage: loads(dumps(SubClass.NestedClass()))
        <sage.misc.nested_class.MainClass.NestedClass object at 0x...>
        sage: loads(dumps(SubClass()))
        <sage.misc.nested_class.SubClass object at 0x...>
    """
    pass

nested_pickle(SubClass)

def _provide_SubClass():
    return SubClass

class CopiedClass(object):
    r"""
    A simple class to test nested_pickle.

    EXAMPLES::

        sage: from sage.misc.nested_class import CopiedClass
        sage: loads(dumps(CopiedClass.NestedClass()))
        <sage.misc.nested_class.MainClass.NestedClass object at 0x...>
        sage: loads(dumps(CopiedClass.NestedSubClass()))
        <sage.misc.nested_class.MainClass.NestedClass.NestedSubClass object at 0x...>
    """
    NestedClass = MainClass.NestedClass
    NestedSubClass = MainClass.NestedClass.NestedSubClass
    SubClass = _provide_SubClass()

nested_pickle(CopiedClass)

# Further classes for recursive tests

class A1:
    class A2:
        class A3:
            pass

