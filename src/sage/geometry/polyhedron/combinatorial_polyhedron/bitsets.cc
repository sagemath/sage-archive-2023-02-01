/*
#############################################################################
# More or less taken from ``src/sage/data_structures/bitset.pxi``
#############################################################################
*/

/*
#############################################################################
# Creating limb patterns
#############################################################################
*/

inline uint64_t limb_one_set_bit(size_t n){
    /*
    Return a limb with only bit n set.
    */
    return (uint64_t) 1 << (n % 64);
}

inline uint64_t limb_one_zero_bit(mp_bitcnt_t n){
    /*
    Return a limb with all bits set, except for bit n.
    */
    return !((uint64_t) 1 << (n % 64));
}

inline uint64_t limb_lower_bits_down(mp_bitcnt_t n){
    /*
    Return a limb with the lower n bits set, where n is interpreted
    in [0 .. 64-1].
    */
    return ((uint64_t) 1 << (n % 64)) - 1;
}

/*
#############################################################################
# Bitset Bit Manipulation
#############################################################################
*/

inline int bitset_in(uint64_t* bits, size_t n){
    /*
    Check if n is in bits.  Return True (i.e., 1) if n is in the
    set, False (i.e., 0) otherwise.
    */
    return bits[n/64] & limb_one_set_bit(n);
}

inline int bitset_discard(uint64_t* bits, size_t n){
    /*
    Remove n from bits.
    */
    bits[n/64] &= limb_one_zero_bit(n);
}

