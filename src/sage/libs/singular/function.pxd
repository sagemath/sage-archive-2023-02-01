"""
Direct Access to Singular's Functions via libSingular.

AUTHORS:

- Michael Brickenstein (2009-07): initial implementation, overall design
- Martin Albrecht (2009-07): clean up, enhancements, etc.
"""
#*****************************************************************************
#       Copyright (C) 2009 Michael Brickenstein <brickenstein@mfo.de>
#       Copyright (C) 2009 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object cimport SageObject
from sage.libs.singular.decl cimport leftv, idhdl
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular

cdef class Converter(SageObject)

cdef class BaseCallHandler:
    cdef leftv* handle_call(self, Converter argument_list) except NULL
    cdef bint free_res(self)

cdef class LibraryCallHandler(BaseCallHandler):
    cdef idhdl * proc_idhdl

cdef class KernelCallHandler(BaseCallHandler):
    cdef long arity
    cdef long cmd_n

cdef class SingularFunction(SageObject):
    cdef object _name
    cdef MPolynomialRing_libsingular _ring
    cdef BaseCallHandler call_handler

    cdef BaseCallHandler get_call_handler(self)
    cdef bint function_exists(self)
    cdef MPolynomialRing_libsingular common_ring(self, tuple args, ring=?)

cdef class SingularLibraryFunction(SingularFunction):
    pass

cdef class SingularKernelFunction(SingularFunction):
    pass

# the most direct function call interface
cdef inline call_function(SingularFunction self, tuple args, MPolynomialRing_libsingular R, bint signal_handler=?)
