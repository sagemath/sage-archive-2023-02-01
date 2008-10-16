"""
Polynomial Template for C/C++ Library Interfaces
"""

from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class Polynomial_template(Polynomial):
    cdef celement x

