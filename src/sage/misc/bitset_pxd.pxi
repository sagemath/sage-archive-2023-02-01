#*****************************************************************************
#     Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************



# This is a .pxi file so that one can inline functions. Doctests in misc_c.

include "../ext/stdsage.pxi"

cdef extern from *:
    void *memset(void *, int, size_t)
    void *memcpy(void *, void *, size_t)
    int memcmp(void *, void *, size_t)
    size_t strlen(char *)

    # constant literals
    int index_shift "(sizeof(unsigned long)==8 ? 6 : 5)"
    unsigned long offset_mask "(sizeof(unsigned long)==8 ? 0x3F : 0x1F)"

    # Given an element index n in a set, (n >> index_shift) gives the
    # corresponding limb number, while (n & offset_mask) gives the bit
    # number inside of the limb.

cdef struct bitset_s:
    long size
    long limbs
    unsigned long *bits

ctypedef bitset_s bitset_t[1]
