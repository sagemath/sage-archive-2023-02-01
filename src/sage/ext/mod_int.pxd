"""
The `mod_int` Data Type

* In C/C++ headers, you can `#include "mod_int.h"`
* In Cython files, use `from sage.ext.mod_int cimport *`
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2013 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


cdef extern from "sage/ext/mod_int.h":
    ctypedef long mod_int
    mod_int MOD_INT_MAX
    mod_int MOD_INT_OVERFLOW
