###############################################################################
#       Copyright (C) 2009, William Stein <wstein@gmail.com>
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################


cdef extern from "<gap/system.h>":
    ctypedef char libGAP_Char
    ctypedef int libGAP_Int
    ctypedef unsigned char libGAP_UChar

cdef extern from "<gap/code.h>":
    ctypedef unsigned int libGAP_Stat
    ctypedef libGAP_Stat* libGAP_PtrBody

cdef extern from "<gap/gap.h>":
    ctypedef unsigned int libGAP_UInt
    ctypedef void* libGAP_ExecStatus

cdef extern from "<gap/objects.h>":
    ctypedef void* libGAP_Obj
