r"""
Dynamic documentation for instances of classes

The functionality in this module allows to define specific docstrings
of *instances* of a class, which are different from the class docstring.
A typical use case is given by cached methods: the documentation of a
cached method should not be the documentation of the class
:class:`CachedMethod`; it should be the documentation of the underlying
method.

In order to use this, define a class docstring as usual. Also define a
method ``def _instancedoc_(self)`` which should return the docstring of
the instance ``self``. Finally, add the decorator ``@instancedoc`` to
the class.

.. WARNING::

    Since the ``__doc__`` attribute is never inherited, the decorator
    ``@instancedoc`` must be added to all subclasses of the class
    defining ``_instancedoc_``. Doing it on the base class is not
    sufficient.

EXAMPLES::

    sage: from sage.docs.instancedoc import instancedoc
    sage: @instancedoc
    ....: class X(object):
    ....:     "Class docstring"
    ....:     def _instancedoc_(self):
    ....:         return "Instance docstring"
    sage: X.__doc__
    'Class docstring'
    sage: X().__doc__
    'Instance docstring'

For a Cython ``cdef class``, a decorator cannot be used. Instead, call
:func:`instancedoc` as a function after defining the class::

    sage: cython('''
    ....: from sage.docs.instancedoc import instancedoc
    ....: cdef class Y:
    ....:     "Class docstring"
    ....:     def _instancedoc_(self):
    ....:         return "Instance docstring"
    ....: instancedoc(Y)
    ....: ''')
    sage: Y.__doc__
    'File:...\nClass docstring'
    sage: Y().__doc__
    'Instance docstring'

TESTS:

Check that inheritance works (after passing the subclass to
:func:`instancedoc`)::

    sage: @instancedoc
    ....: class A(object):
    ....:     "Class A docstring"
    ....:     def _instancedoc_(self):
    ....:         return "Instance docstring"
    sage: class B(A):
    ....:     "Class B docstring"
    sage: B.__doc__
    'Class B docstring'
    sage: B().__doc__  # Ideally, this would return the instance docstring
    'Class B docstring'
    sage: B = instancedoc(B)
    sage: B.__doc__
    'Class B docstring'
    sage: B().__doc__
    'Instance docstring'
"""

#*****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.object cimport PyObject, PyTypeObject

cdef extern from *:
    cdef int PyDict_SetItemString(PyObject*, const char*, object) except -1
    cdef void PyType_Modified(PyTypeObject*)

cdef inline PyTypeObject* TypeObject(cls) except NULL:
    if not isinstance(cls, type):
        raise TypeError(f"expected type, got {cls!r}")
    return <PyTypeObject*>cls


cdef class InstanceDocDescriptor:
    """
    Descriptor for dynamic documentation, to be installed as the
    ``__doc__`` attribute.

    INPUT:

    - ``classdoc`` -- (string) class documentation

    - ``instancedoc`` -- (method) documentation for an instance

    EXAMPLES::

        sage: from sage.docs.instancedoc import InstanceDocDescriptor
        sage: def instancedoc(self):
        ....:     return "Instance doc"
        sage: docattr = InstanceDocDescriptor("Class doc", instancedoc)
        sage: class Z(object):
        ....:     __doc__ = InstanceDocDescriptor("Class doc", instancedoc)
        sage: Z.__doc__
        'Class doc'
        sage: Z().__doc__
        'Instance doc'
    """
    cdef classdoc
    cdef instancedoc

    def __init__(self, classdoc, instancedoc):
        """
        TESTS::

            sage: from sage.docs.instancedoc import InstanceDocDescriptor
            sage: InstanceDocDescriptor(None, None)
            <sage.docs.instancedoc.InstanceDocDescriptor object at ...>
        """
        self.classdoc = classdoc
        self.instancedoc = instancedoc

    def __get__(self, obj, typ):
        """
        TESTS::

            sage: from sage.docs.instancedoc import InstanceDocDescriptor
            sage: def instancedoc(self):
            ....:     return "Doc for {!r}".format(self)
            sage: descr = InstanceDocDescriptor("Class doc", instancedoc)
            sage: descr.__get__(None, object)
            'Class doc'
            sage: descr.__get__(42, type(42))
            'Doc for 42'
        """
        if obj is None:
            return self.classdoc
        else:
            return self.instancedoc(obj)


def instancedoc(cls):
    """
    Add support for ``_instancedoc_`` to the class ``cls``.

    Typically, this will be used as decorator.

    INPUT:

    - ``cls`` -- a new-style class

    OUTPUT: ``cls``

    .. WARNING::

        ``instancedoc`` mutates the given class. So you are *not* supposed
        to use it as ``newcls = instancedoc(cls)`` because that would
        mutate ``cls`` (and ``newcls`` would be the same object as ``cls``)

    TESTS:

    We get a useful error message if ``_instancedoc_`` is not defined::

        sage: from sage.docs.instancedoc import instancedoc
        sage: class X(object): pass
        sage: instancedoc(X)
        Traceback (most recent call last):
        ...
        TypeError: instancedoc requires <class '__main__.X'> to have an '_instancedoc_' attribute

    This does not work on old-style classes or things which are not a
    class at all::

        sage: instancedoc(7)
        Traceback (most recent call last):
        ...
        TypeError: expected type, got 7
        sage: class OldStyle: pass
        sage: instancedoc(OldStyle)
        Traceback (most recent call last):
        ...
        TypeError: expected type, got <class __main__.OldStyle at ...>
    """
    cdef PyTypeObject* tp = TypeObject(cls)
    try:
        instdoc = cls._instancedoc_
    except AttributeError:
        raise TypeError(f"instancedoc requires {cls!r} to have an '_instancedoc_' attribute")
    docattr = InstanceDocDescriptor(cls.__doc__, instdoc)
    PyDict_SetItemString(tp.tp_dict, "__doc__", docattr)
    tp.tp_doc = NULL
    PyType_Modified(tp)
    return cls
