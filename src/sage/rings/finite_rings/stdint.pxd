"""
C Integer Types Used in Finite Rings
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2013 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cimport libc.stdint

cdef extern from "sage/rings/finite_rings/stdint.h":
    ctypedef libc.stdint.int_fast32_t int_fast32_t
    ctypedef libc.stdint.int_fast64_t int_fast64_t
    int_fast32_t INTEGER_MOD_INT32_LIMIT
    int_fast64_t INTEGER_MOD_INT64_LIMIT

