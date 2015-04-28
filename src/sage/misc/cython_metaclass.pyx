"""
Metaclasses for Cython extension types

Cython does not support metaclasses, but this module can be used to
implement metaclasses for extension types.

.. WARNING::

    This module has many caveats and you can easily get segfaults if you
    make a mistake. It relies on undocumented Python and Cython
    behaviour, so things might break in future versions.

How to use
==========

To enable this metaclass mechanism, you need to put
``cimport sage.misc.cython_metaclass`` in your module (in the ``.pxd``
file if you are using one).

In the extension type (a.k.a. ``cdef class``) for which you want to
define a metaclass, define a method ``__getmetaclass__`` with a single
unused argument. This method should return a type to be used as
metaclass:

.. code-block:: cython

    cimport sage.misc.cython_metaclass
    cdef class MyCustomType(object):
        def __getmetaclass__(_):
            from foo import MyMetaclass
            return MyMetaclass

The above ``__getmetaclass__`` method is analogous to
``__metaclass__ = MyMetaclass`` in Python 2.

.. WARNING::

    ``__getmetaclass__`` must be defined as an ordinary method taking a
    single argument, but this argument should not be used in the
    method (it will be ``None``).

When a type ``cls`` is being constructed with metaclass ``meta``,
then ``meta.__init__(cls, None, None, None)`` is called from Cython.
In Python, this would be ``meta.__init__(cls, name, bases, dict)``.

.. WARNING::

    The ``__getmetaclass__`` method is called while the type is being
    created during the import of the module. Therefore,
    ``__getmetaclass__`` should not refer to any global objects,
    including the type being created or other types defined or imported
    in the module. Note that ``cdef`` functions are not Python objects,
    so those are safe to call.

    The same warning applies to the ``__init__`` method of the
    metaclass.

EXAMPLES::

    sage: cython('''
    ....: cimport sage.misc.cython_metaclass
    ....: cdef class MyCustomType(object):
    ....:     def __getmetaclass__(_):
    ....:         class MyMetaclass(type):
    ....:             def __init__(*args):
    ....:                 print "Calling MyMetaclass.__init__{}".format(args)
    ....:         return MyMetaclass
    ....:
    ....: cdef class MyDerivedType(MyCustomType):
    ....:     pass
    ....: ''')
    Calling MyMetaclass.__init__(<type '...MyCustomType'>, None, None, None)
    Calling MyMetaclass.__init__(<type '...MyDerivedType'>, None, None, None)
    sage: MyCustomType.__class__
    <class '...MyMetaclass'>
    sage: class MyPythonType(MyDerivedType):
    ....:     pass
    Calling MyMetaclass.__init__(<class '...MyPythonType'>, 'MyPythonType', (<type '...MyDerivedType'>,), {'__module__': '__main__'})

Application: patching extension types
=====================================

Sometimes it is needed to "patch" an extension type to implement
something using the Python C API which cannot be implemented in plain
Cython. For example, Cython does not allow defining a custom ``tp_new``
function for a type. Therefore, if we do want a custom ``tp_new``, we
need to manually change it in the ``PyTypeObject`` structure.

For this, we provide :class:`TypeInitMetaclass`, which is a metaclass
which calls ``cls.__typeinit__(cls)`` when ``cls`` is being created.

.. code-block:: cython

    cimport sage.misc.cython_metaclass
    from cpython.object cimport PyTypeObject

    cdef class Foo(object):
        def __getmetaclass__(_):
            from sage.misc.cython_metaclass import TypeInitMetaclass
            return TypeInitMetaclass

        def __typeinit__(cls):
            if not isinstance(cls, type):
                raise TypeError
            (<PyTypeObject*>cls).tp_new = ...

Note that ``__typeinit__`` is inherited as usual, so the class can also
be a subclass.

.. WARNING::

    ``__typeinit__`` must be defined as ordinary method taking a single
    argument, but it is called as if it was a class method (this is
    to work around the way how Cython initializes class methods).
    To be safe, check ``isinstance(cls, type)`` on the ``cls`` argument
    before casting it to ``PyTypeObject*``.

Implementation
==============

All this is implemented by defining

.. code-block:: c

    #define PyTypeReady(t)  Sage_PyType_Ready(t)

and then implementing the function ``Sage_PyType_Ready(t)`` which first
calls ``PyType_Ready(t)`` and then handles the metaclass stuff.

TESTS:

Check that a proper exception is raised if ``__getmetaclass__``
returns a non-type::

    sage: cython('''
    ....: cimport sage.misc.cython_metaclass
    ....: cdef class MyCustomType(object):
    ....:     def __getmetaclass__(_):
    ....:         return 2
    ....: ''')
    Traceback (most recent call last):
    ...
    TypeError: __getmetaclass__ did not return a type

Reference
=========
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef class TypeInitMetaclass(type):
    """
    Metaclass which calls ``cls.__typeinit__(cls)`` when ``cls``
    is initialized.

    For this to work from Cython, you need to add
    ``cimport sage.misc.cython_metaclass`` in the Cython module.

    EXAMPLES::

        sage: cython('''
        ....: cimport sage.misc.cython_metaclass
        ....: from cpython.object cimport PyTypeObject
        ....:
        ....: cdef class Foo(object):
        ....:     def __getmetaclass__(_):
        ....:         from sage.misc.cython_metaclass import TypeInitMetaclass
        ....:         return TypeInitMetaclass
        ....:
        ....:     def __typeinit__(cls):
        ....:         if not isinstance(cls, type):
        ....:             raise TypeError
        ....:         print "Calling __typeinit__ of", (<PyTypeObject*>cls).tp_name
        ....:
        ....: cdef class Bar(Foo):
        ....:     pass
        ....: ''')
        Calling __typeinit__ of ...Foo
        Calling __typeinit__ of ...Bar
        sage: class PythonFoo(Foo):
        ....:     pass
        Calling __typeinit__ of PythonFoo
    """
    def __init__(self, *args):
        """
        TESTS:

        We define a ``__typeinit__`` method with an extra argument and
        check that we get a proper exception::

            sage: cython('''
            ....: cimport sage.misc.cython_metaclass
            ....: cdef class BadTypeInit(object):
            ....:     def __getmetaclass__(_):
            ....:         from sage.misc.cython_metaclass import TypeInitMetaclass
            ....:         return TypeInitMetaclass
            ....:     def __typeinit__(cls, arg):
            ....:         pass
            ....: ''')
            Traceback (most recent call last):
            ...
            TypeError: PyMethodDescr_CallSelf requires a method without arguments
        """
        PyMethodDescr_CallSelf(self.__typeinit__, self)
