"""
Parent objects with generators
"""

###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

cimport parent

cdef class ParentWithGens(parent.Parent):
    cdef public object _gens
    cdef public object _gens_dict
    cdef public object _names
    cdef public object _latex_names
    cdef public object _list
    cdef public object _base

cdef class ParentWithMultiplicativeAbelianGens(ParentWithGens):
    cdef public object _generator_orders

cdef class ParentWithAdditiveAbelianGens(ParentWithGens):
    cdef public object _generator_orders
