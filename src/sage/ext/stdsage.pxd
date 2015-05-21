"""
Standard C helper code for Cython modules
"""
#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include 'python_rich_object.pxi'
from cpython.ref cimport Py_TYPE


cdef inline PY_NEW(type t):
    """
    Return ``t.__new__(t)``.  This works even for types like
    :class:`Integer` where we change ``tp_new`` at runtime (Cython
    optimizations assume that ``tp_new`` doesn't change).
    """
    return (<RichPyTypeObject*>t).tp_new(<PyTypeObject*>t, NULL, NULL)


cdef inline void PY_SET_TP_NEW(type dst, type src):
    """
    Manually set ``dst.__new__`` to ``src.__new__``.  This is used to
    speed up Cython's boilerplate object construction code by skipping
    irrelevant base class ``tp_new`` methods.
    """
    (<RichPyTypeObject*>dst).tp_new = (<RichPyTypeObject*>src).tp_new


cdef inline bint HAS_DICTIONARY(obj):
    """
    Test whether the given object has a Python dictionary.
    """
    return (<RichPyTypeObject*>Py_TYPE(obj)).tp_dictoffset != 0
