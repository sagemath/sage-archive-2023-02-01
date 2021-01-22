# ****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************


# Operations on digits (intended to be small elements in the exact subring)
###########################################################################

cdef cdigit digit_zero
digit_init(digit_zero)

cdef inline void digit_init(cdigit a):
    pass

cdef inline void digit_clear(cdigit a):
    pass

# get and set

cdef inline Integer digit_get_sage(cdigit a):
    pass

cdef inline void digit_set(cdigit a, cdigit b):
    pass

cdef inline void digit_set_ui(cdigit a, slong b):
    pass

cdef inline void digit_set_sage(cdigit a, Integer elt):
    pass

# comparisons

cdef inline bint digit_equal(cdigit a, cdigit b):
    pass

cdef inline bint digit_equal_ui(cdigit a, slong b):
    pass

cdef inline bint digit_is_zero(cdigit a):
    pass

# operations

cdef inline void digit_random(cdigit res, PowComputer_class prime_pow):
    pass

cdef inline void digit_add(cdigit res, cdigit a, cdigit b):
    pass

cdef inline void digit_sub(cdigit res, cdigit a, cdigit b):
    pass

cdef inline void digit_mul(cdigit res, cdigit a, cdigit b):
    pass

cdef inline void digit_mod(cdigit res, cdigit a, PowComputer_class prime_pow):
    pass

cdef inline void digit_quorem(cdigit quo, cdigit rem, cdigit a, PowComputer_class prime_pow):
    pass

cdef inline void digit_inv(cdigit res, cdigit a, PowComputer_class prime_pow):
    pass

cdef bint digit_sqrt(cdigit_ptr ans, cdigit_ptr x, PowComputer_class prime_pow):
    pass


# Operations on elements (represented as series of digits)
##########################################################

cdef inline void element_init(celement x):
    pass

cdef inline void element_clear(celement x):
    pass

# get and set

cdef inline element_get_sage(celement x, PowComputer_class prime_pow):
    pass

cdef inline void element_set(celement x, celement y):
    pass

# get and set digits

cdef inline cdigit_ptr element_get_digit(celement x, slong i):
    pass

cdef inline Integer element_get_digit_sage(celement x, slong i):
    pass

cdef inline void element_get_slice(celement res, celement x, slong start, slong length):
    pass

cdef inline void element_set_digit(celement x, cdigit a, slong i):
    pass

cdef inline void element_set_digit_ui(celement x, slong a, slong i):
    pass

cdef inline void element_set_digit_sage(celement x, Integer a, slong i):
    pass

# operations

cdef inline void element_iadd_digit(celement x, cdigit a, slong i):
    pass

cdef inline void element_isub_digit(celement x, cdigit a, slong i):
    pass

cdef inline void element_iadd_slice(celement x, celement slice, slong start):
    pass

cdef inline void element_isub_slice(celement x, celement slice, slong start):
    pass

cdef inline void element_scalarmul(celement res, celement x, cdigit a):
    pass

cdef inline void element_mul(celement res, celement x, celement y):
    pass

cdef inline void element_reduce_digit(celement x, slong i, PowComputer_class prime_pow):
    pass

cdef inline void element_shift_right(celement x):
    pass
