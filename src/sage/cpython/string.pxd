# -*- encoding: utf-8 -*-
#*****************************************************************************
#       Copyright (C) 2017 Erik M. Bray <erik.bray@lri.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from libc.string cimport strlen

from cpython.bytes cimport PyBytes_AS_STRING, PyBytes_FromString
from cpython.unicode cimport PyUnicode_Decode, PyUnicode_AsEncodedString

IF PY_MAJOR_VERSION >= 3:
    cdef extern from "Python.h":
        # Missing from cpython.unicode in Cython 0.27.3
        char* PyUnicode_AsUTF8(object unicode)
        unicode PyUnicode_DecodeLocale(const char* str, const char* errors)
        bytes PyUnicode_EncodeLocale(object unicode, const char* errors)


cdef inline str char_to_str(char* c, encoding=None, errors=None):
    IF PY_MAJOR_VERSION <= 2:
        return <str>PyBytes_FromString(c)
    ELSE:
        cdef char* err
        if errors is None:
            err = NULL  # implies "strict"
        else:
            err = PyUnicode_AsUTF8(errors)

        if encoding is None:
            return PyUnicode_DecodeLocale(c, err)

        return PyUnicode_Decode(c, strlen(c), PyUnicode_AsUTF8(encoding), err)


cpdef inline bytes_to_str(b, encoding=None, errors=None):
    """
    Convert ``bytes`` to ``str``.

    On Python 2 this is a no-op since ``bytes is str``.  On Python 3
    this decodes the given ``bytes`` to a Python 3 unicode ``str`` using
    the specified encoding.

    EXAMPLES::

        sage: import six
        sage: from sage.cpython.string import bytes_to_str
        sage: s = bytes_to_str(b'\xe2\x98\x83')
        sage: if six.PY2:
        ....:     s == b'\xe2\x98\x83'
        ....: else:
        ....:     s == u'☃'
        True
        sage: bytes_to_str([])
        Traceback (most recent call last):
        ...
        TypeError: expected bytes, list found
    """
    if not isinstance(b, bytes):
        raise TypeError(f"expected bytes, {type(b).__name__} found")

    IF PY_MAJOR_VERSION <= 2:
        return b
    ELSE:
        return char_to_str(PyBytes_AS_STRING(b), encoding=encoding,
                           errors=errors)


cpdef inline str_to_bytes(s, encoding=None, errors=None):
    """
    Convert ``str`` to ``bytes``.

    On Python 2 this is a no-op since ``str is bytes``.  On Python 3
    this encodes the given ``str`` to a Python 3 ``bytes`` using the
    specified encoding.

    EXAMPLES::

        sage: import six
        sage: from sage.cpython.string import str_to_bytes
        sage: if six.PY2:
        ....:     b = str_to_bytes('\xe2\x98\x83')
        ....: else:
        ....:     b = str_to_bytes(u'☃')
        sage: b == b'\xe2\x98\x83'
        True
        sage: str_to_bytes([])
        Traceback (most recent call last):
        ...
        TypeError: expected str, list found
    """
    # Make this check explicit to avoid obscure error message below
    if not isinstance(s, str):
        raise TypeError(f"expected str, {type(s).__name__} found")

    IF PY_MAJOR_VERSION <= 2:
        return s
    ELSE:
        cdef char* err
        if errors is None:
            err = NULL  # implies "strict"
        else:
            err = PyUnicode_AsUTF8(errors)

        if encoding is None:
            return PyUnicode_EncodeLocale(s, err)

        return PyUnicode_AsEncodedString(s, PyUnicode_AsUTF8(encoding), err)
