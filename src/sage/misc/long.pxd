r"""
Fast conversion macro to long
"""

from libc.limits cimport LONG_MIN

from cpython.int cimport PyInt_CheckExact, PyInt_AS_LONG
from cpython.long cimport PyLong_CheckExact, PyLong_AsLong
from cpython.number cimport PyNumber_Index

from sage.libs.gmp.mpz cimport mpz_fits_slong_p, mpz_get_si

from sage.rings.integer cimport Integer

cdef inline long pyobject_to_long(x) except LONG_MIN:
    r"""
    Given a Python object ``x`` cast it quickly to a C long.

    If it fails, it either raises a ``TypeError`` or a ``OverflowError`` (and
    the exceptional value is `LONG_MIN`)
    """
    if PyInt_CheckExact(x):
        return PyInt_AS_LONG(x)
    elif type(x) is Integer:
        if mpz_fits_slong_p((<Integer>x).value):
            return mpz_get_si((<Integer>x).value)
        else:
            raise OverflowError
    elif PyLong_CheckExact(x):
        return PyLong_AsLong(x)

    return PyNumber_Index(x)
