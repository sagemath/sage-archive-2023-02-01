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
    # The size of a bitset B counts the maximum number of bits that B can
    # hold. This size is independent of how many elements of B are toggled to
    # 1. For example, say B is the bitset 1001. Then B has size 4, with the
    # first and fourth elements toggled to 1, reading from left to right.
    # We can also think of the size of a bitset as its capacity.
    long size

    # A limb is that part of a bitset that can fit into a machine word. On an
    # n-bit machine, a word is usually taken to be n bits in length. Thus on a
    # 32-bit machine, a word usually is 32 bits in length. Similarly, on a
    # 64-bit machine, a word is usually 64 bits in length. If a bitset has
    # size <= n, then the whole bitset fits into a limb and we only require
    # one limb to represent the bitset. However, if the bitset has size > n,
    # we require more than one limb to represent the bitset. For example, if
    # a limb is 64 bits in length and the bitset has size 96 bits, then we
    # require at most two limbs to represent the bitset.
    long limbs

    # The individual bits of a bitset.
    unsigned long *bits

ctypedef bitset_s bitset_t[1]
