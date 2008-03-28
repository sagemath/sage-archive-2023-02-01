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

# This should hopefully be doable in a future release
# of Cython. Basically, I need to be able to do
#
#   ctypedef bitset_s bitset_t[1]
#
# to have stack-allocated pass-by-reference.

cdef extern from "bitset.h":
    struct bitset_s:
        long size
        long limbs
        unsigned long *bits

    ctypedef bitset_s* bitset_t

    int index_shift "(sizeof(unsigned long)==8 ? 6 : 5)"
    unsigned long offset_mask "(sizeof(unsigned long)==8 ? 0x3F : 0x1F)"

#############################################################################
# Bitset Initalization
#############################################################################

cdef inline bint bitset_init(bitset_t bits, unsigned long size) except -1:
    bits.size = size
    bits.limbs = (size - 1)/(8*sizeof(unsigned long)) + 1
    bits.bits = <unsigned long*>sage_malloc(bits.limbs * sizeof(unsigned long))
    if bits.bits == NULL:
        raise MemoryError
    bits.bits[bits.limbs-1] = 0

cdef inline void bitset_clear(bitset_t bits):
    sage_free(bits.bits)

cdef inline void bitset_zero(bitset_t bits):
    memset(bits.bits, 0, (bits.size+7) >> 3)

cdef inline void bitset_copy(bitset_t dst, bitset_t src):
    memcpy(dst.bits, src.bits, (dst.size+7) >> 3)

#############################################################################
# Bitset Comparison
#############################################################################

cdef inline bint bitset_is_zero(bitset_t bits):
    cdef long i
    for i from 0 <= i < bits.limbs:
        if bits.bits[i] != 0:
            return False
    return True

cdef inline bint bitset_is_equal(bitset_t a, bitset_t b):
    return memcmp(a.bits, b.bits, (a.size+7) >> 3) == 0

cdef inline int bitset_cmp(bitset_t a, bitset_t b):
    cdef long i
    for i from a.limbs > i >= 0:
        if a.bits[i] != b.bits[i]:
            if a.bits[i] < b.bits[i]:
                return -1
            else:
                return 1
    return 0

#############################################################################
# Bitset Bit Manipulation
#############################################################################

cdef inline bint bitset_check(bitset_t bits, unsigned long n):
    return (bits.bits[n >> index_shift] >> (n & offset_mask)) & 1

cdef inline void bitset_unset(bitset_t bits, unsigned long n):
    bits.bits[n >> index_shift] &= ~((<unsigned long>1) << (n & offset_mask))

cdef inline void bitset_set(bitset_t bits, unsigned long n):
    bits.bits[n >> index_shift] |= (<unsigned long>1) << (n & offset_mask)

cdef inline void bitset_set_to(bitset_t bits, unsigned long n, bint b):
    bitset_unset(bits, n)
    bits.bits[n >> index_shift] |= (<unsigned long>b) << (n & offset_mask)

cdef inline void bitset_flip(bitset_t bits, unsigned long n):
    bits.bits[n >> index_shift] ^= (<unsigned long>1) << (n & offset_mask)

#############################################################################
# Bitset Searching
#############################################################################

cdef inline long _bitset_first_in_limb(unsigned long limb):
    cdef long j
    if limb & (((<unsigned long>1) << 4*sizeof(unsigned long)) - 1):
        for j from 0 <= j < 4*sizeof(unsigned long):
            if limb & ((<unsigned long>1) << j):
                return j
    else:
        for j from 4*sizeof(unsigned long) <= j < 8*sizeof(unsigned long):
            if limb & ((<unsigned long>1) << j):
                return j
    return -1

cdef inline long bitset_first(bitset_t a):
    cdef long i
    for i from 0 <= i < a.limbs:
        if a.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb(a.bits[i])
    return -1

cdef inline long bitset_first_diff(bitset_t a, bitset_t b):
    cdef long i
    for i from 0 <= i < a.limbs:
        if a.bits[i] != b.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb(a.bits[i] ^ b.bits[i])
    return -1

cdef inline long bitset_next(bitset_t a, long n):
    if n >= a.size:
        return -1
    cdef long i
    cdef long limb = a.bits[n >> index_shift] & ~(((<unsigned long>1) << (n & offset_mask)) - 1)
    cdef long ret = _bitset_first_in_limb(limb)
    if ret != -1:
        return (n & ~offset_mask) | ret
    for i from (n >> index_shift) < i < a.limbs:
        if a.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb(a.bits[i])
    return -1

cdef inline long bitset_next_diff(bitset_t a, bitset_t b, long n):
    if n >= a.size:
        return -1
    cdef long i
    cdef long limb = (a.bits[n >> index_shift] ^ b.bits[n >> index_shift]) \
        & ~(((<unsigned long>1) << (n & offset_mask)) - 1)
    cdef long ret = _bitset_first_in_limb(limb)
    if ret != -1:
        return (n & ~offset_mask) | ret
    for i from (n >> index_shift) < i < a.limbs:
        if a.bits[i] != b.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb(a.bits[i] ^ b.bits[i])
    return -1

#############################################################################
# Bitset Arithmatic
#############################################################################

cdef inline void bitset_not(bitset_t r, bitset_t a):
    cdef long i
    for i from 0 <= i < r.limbs:
        r.bits[i] = ~a.bits[i]

cdef inline void bitset_and(bitset_t r, bitset_t a, bitset_t b):
    cdef long i
    for i from 0 <= i < r.limbs:
        r.bits[i] = a.bits[i] & b.bits[i]

cdef inline void bitset_or(bitset_t r, bitset_t a, bitset_t b):
    cdef long i
    for i from 0 <= i < r.limbs:
        r.bits[i] = a.bits[i] | b.bits[i]

cdef inline void bitset_xor(bitset_t r, bitset_t a, bitset_t b):
    cdef long i
    for i from 0 <= i < r.limbs:
        r.bits[i] = a.bits[i] ^ b.bits[i]

cdef void bitset_rshift(bitset_t r, bitset_t a, long n):
    if n <= 0:
        if n != 0:
            bitset_lshift(r, a, -n)
        return
    elif n >= a.size:
        bitset_zero(r)
        return
    cdef long i
    cdef long off = n >> index_shift
    cdef int shift = offset_mask & n
    cdef int shift2 = 8*sizeof(unsigned long) - shift
    if shift == 0:
        for i from 0 <= i < r.limbs - off:
            r.bits[i] = a.bits[i+off]

    else:
        for i from 0 <= i < r.limbs - off - 1:
            r.bits[i] = (a.bits[i+off] >> shift) | (a.bits[i+off+1] << shift2)
        r.bits[r.limbs - off - 1] = a.bits[r.limbs - 1] >> shift

    if off > 0:
        memset(r.bits + r.limbs - off, 0, off * sizeof(unsigned long))

cdef void bitset_lshift(bitset_t r, bitset_t  a, long n):
    if n <= 0:
        if n != 0:
            bitset_rshift(r, a, -n)
        return
    elif n >= a.size:
        bitset_zero(r)
        return
    cdef long i
    cdef long off = n >> index_shift
    cdef int shift = offset_mask & n
    cdef int shift2 = 8*sizeof(unsigned long) - shift
    if shift == 0:
        for i from r.limbs - off > i >= 0:
            r.bits[i+off] = a.bits[i]

    else:
        for i from r.limbs - off > i >= 1:
            r.bits[i+off] = (a.bits[i] << shift) | (a.bits[i] >> shift2)
        r.bits[off] = a.bits[0] << shift

    if off > 0:
        memset(r.bits, 0, off * sizeof(unsigned long))

#############################################################################
# Hamming Weights
#############################################################################

cdef enum:
    _bitset_hamming_table_bits = 8

cdef int _bitset_hamming_table[1 << _bitset_hamming_table_bits]

cdef void _bitset_fill_hamming_table():
    cdef int i, j
    for i from 0 <= i < (1 << _bitset_hamming_table_bits):
        _bitset_hamming_table[i] = 0
        for j from 0 <= j < _bitset_hamming_table_bits:
            _bitset_hamming_table[i] += (i >> j) & 1

_bitset_fill_hamming_table()

cdef inline int bitset_hamming_weight(bitset_t a):
    cdef long i, j
    cdef long w = 0
    cdef unsigned long limb
    for i from 0 <= i < a.limbs:
        limb = a.bits[i]
        for j from 0 <= j < (8*sizeof(unsigned long) + _bitset_hamming_table_bits - 1) / _bitset_hamming_table_bits:
            w += _bitset_hamming_table[(limb >> (j*_bitset_hamming_table_bits)) & ((1 << _bitset_hamming_table_bits) - 1)]
    return w

cdef inline long bitset_hamming_weight_sparse(bitset_t a):
    cdef long i
    cdef long w = 0
    cdef unsigned long limb
    for i from 0 <= i < a.limbs:
        limb = a.bits[i]
        while limb:
            w += _bitset_hamming_table[limb & ((1 << _bitset_hamming_table_bits) - 1)]
            limb = limb >> _bitset_hamming_table_bits
    return w

#############################################################################
# Bitset Conversion
#############################################################################

cdef char* bitset_chars(char* s, bitset_t bits, char zero=c'0', char one=c'1'):
    cdef long i
    if s == NULL:
        s = <char *>sage_malloc(bits.size+1)
    for i from 0 <= i < bits.size:
        s[i] = one if bitset_check(bits, i) else zero
    s[bits.size] = 0
    return s

cdef void bitset_from_str(bitset_t bits, char* s, char zero=c'0', char one=c'1'):
    bitset_init(bits, strlen(s))
    cdef long i
    for i from 0 <= i < bits.size:
        bitset_set_to(bits, i, s[i] == one)

cdef bitset_string(bitset_t bits):
    cdef char* s = bitset_chars(NULL, bits)
    py_s = s
    sage_free(s)
    return py_s

