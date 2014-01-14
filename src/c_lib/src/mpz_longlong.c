/* long long -> mpz conversion
*
* Author:  Robert Bradshaw <robertwb@math.washington.edu>
* Date:    November 2008
* License: GPL v2 or later
*
* this is free software: if it breaks, you get to keep all the pieces
*/

#include "mpz_longlong.h"
#include <assert.h>

void mpz_set_longlong(mpz_ptr z, long long val) {
    if (sizeof(mp_limb_t) == sizeof(long long)) {
        mpz_set_si(z, val);
    }
    else if (2*sizeof(mp_limb_t) == sizeof(long long)) {
        mp_limb_t* limbs = (mp_limb_t*)&val;
        int sign = (val > 0) - (val < 0);
        val *= sign;
        _mpz_realloc (z, 2);
        const int endianness = 1;
        if (((char*)&endianness)[0]) {
            // little endian
            z->_mp_d[0] = limbs[0];
            z->_mp_d[1] = limbs[1];
        }
        else {
            // big endian
            z->_mp_d[0] = limbs[1];
            z->_mp_d[1] = limbs[0];
        }
        z->_mp_size = sign * (2 - !z->_mp_d[1]);
    }
    else {
        assert(sizeof(mp_limb_t) == sizeof(long long) || 2*sizeof(mp_limb_t) == sizeof(long long));
    }
}


void mpz_set_ulonglong(mpz_ptr z, unsigned long long val) {
    if (sizeof(mp_limb_t) == sizeof(unsigned long long)) {
        mpz_set_ui(z, (mp_limb_t)val);
    }
    else if (2*sizeof(mp_limb_t) == sizeof(unsigned long long)) {
        mp_limb_t* limbs = (mp_limb_t*)&val;
        _mpz_realloc (z, 2);
        const int endianness = 1;
        if (((char*)&endianness)[0]) {
            // little endian
            z->_mp_d[0] = limbs[0];
            z->_mp_d[1] = limbs[1];
        }
        else {
            // big endian
            z->_mp_d[0] = limbs[1];
            z->_mp_d[1] = limbs[0];
        }
        z->_mp_size = (!!z->_mp_d[1] + (z->_mp_d[1] || z->_mp_d[0])); // 0, 1, or 2
    }
    else {
        assert(sizeof(mp_limb_t) == sizeof(unsigned long long) || 2*sizeof(mp_limb_t) == sizeof(unsigned long long));
    }
}
