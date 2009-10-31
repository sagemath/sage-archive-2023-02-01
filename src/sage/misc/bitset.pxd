include 'bitset_pxd.pxi'

cdef class FrozenBitset:
    cdef bitset_t _bitset
    cdef FrozenBitset _new(self,long int capacity)
    cpdef FrozenBitset _larger_capacity_(self, long size)
    cpdef long capacity(self)
    cpdef bint isempty(self)
    cpdef bint issubset(self, FrozenBitset other)
    cpdef bint issuperset(self, FrozenBitset other)
    cpdef bint isdisjoint(self, FrozenBitset other)
    cpdef _union(self, FrozenBitset other)
    cpdef intersection(self, FrozenBitset other)
    cpdef difference(self, FrozenBitset other)
    cpdef symmetric_difference(self, FrozenBitset other)
    cpdef __copy__(self)

cdef class Bitset(FrozenBitset):
    cpdef __copy__(self)
    cpdef update(self, FrozenBitset other)
    cpdef intersection_update(self, FrozenBitset other)
    cpdef difference_update(self, FrozenBitset other)
    cpdef symmetric_difference_update(self, FrozenBitset other)
    cpdef add(self, unsigned long n)
    cpdef remove(self, unsigned long n)
    cpdef discard(self, unsigned long n)
    cpdef pop(self)
    cpdef clear(self)

