"""
Some arithmetic functions that use GMP
"""

"""
To use these functions in your pyrex extension:
  1. Put
         cimport arith_gmp
     at the top of your extension file.

  2. Put the following at the top of your file:

         cimport arith_gmp
         import arith_gmp
         cdef arith_gmp.functions ag
         ag = arith_gmp.functions()

     Then the functions defined below will be calleable by typing
             ag.function_name.

  3. Put "arith_gmp" in the first list in file_setup.py.

  See polynomial_pyx.pyx for an example.
"""

include "gmp.pxi"
include "interrupt.pxi"

cdef class functions:
    cdef public int mpz_crt(self, mpz_t z, mpz_t x, mpz_t y, mpz_t m, mpz_t n) except -1:
        """
        Use the CRT to find the unique z such that z = x (mod m) and z = y (mod n),
        and z lies between 0 and n*m-1, inclusive.

        The modulo m and n are assumed to be coprime.
        The mpz_t z is assumed to have been initialized
        """
        cdef mpz_t g, s, t, mn
        _sig_on
        mpz_init(g); mpz_init(s); mpz_init(t); mpz_init(mn)
        mpz_gcdext(g, s, t, m, n)
        if mpz_cmp_si(g,1) != 0:
            mpz_clear(g); mpz_clear(s); mpz_clear(t); mpz_clear(mn)
            _sig_off
            raise ArithmeticError, "Moduli must be coprime (but gcd=%s)."%mpz_to_str(g)
        # Now s*m + t*n = 1, so mod n, s*m = 1.
        # The answer is x + (y-x)*s*m, since mod m this is x, and mod n this is y
        mpz_sub(z, y, x)     # z = y-x
        mpz_mul(z, z, s)     # z = (y-x)*s
        mpz_mul(z, z, m)     # z = (y-x)*s*m
        mpz_add(z, z, x)     # z = z + (y-x)*s*m
        mpz_mul(mn, m, n)
        mpz_mod(z, z, mn)
        # Finally
        mpz_clear(g); mpz_clear(s); mpz_clear(t); mpz_clear(mn)
        _sig_off
        return 0

    cdef public int mpz_vec(self, mpz_t **z, int *v, int n) except -1:
        """
        Given a C array v of ints, and an mpz_t ** that has not been initialized
        in any way, construct a C array of elements mpz_t elements.
        """
        cdef int i
        z[0] = <mpz_t*> PyMem_Malloc(sizeof(mpz_t)*n)
        if z[0] == <mpz_t*> 0:
            raise MemoryError
        for i from 0 <= i < n:
            mpz_init_set_si(z[0][i], v[i])
        return 0

    cdef public int mpzvec_to_intmod(self, unsigned long **z, mpz_t *v, int n, unsigned long p) except -1:
        """
        Take each element of v and reduce it modulo p.  The vector v must
        have length n.  z should not have been initialized in any way, and
        the calling function is responsible for de-allocating z.
        """
        cdef int i
        z[0] = <unsigned long*> PyMem_Malloc(sizeof(unsigned long)*n)
        if z[0] == <unsigned long *>0:
            raise MemoryError, "Error allocating memory"
        _sig_on
        for i from 0 <= i < n:
            z[0][i] = mpz_fdiv_ui(v[i], p)
        _sig_off
        return 0

    cdef public int intmodvec_to_mpz(self, mpz_t **z, unsigned long *v, int n) except -1:
        """
        Take each element of v and create the corresponding mpz_t (GMP integer).
        The vector v must have length n.
        z should not have been initialized in any way, and
        the calling function is responsible for de-allocating z.
        """
        cdef int i
        z[0] = <mpz_t *> PyMem_Malloc(sizeof(mpz_t) * n)
        if z[0] == <mpz_t *>0:
            raise MemoryError
        _sig_on
        for i from 0 <= i < n:
            mpz_init_set_ui(z[0][i], v[i])
        _sig_off
        return 0

    cdef public int allocate_mpz_zero_array(self, mpz_t **z, int n) except -1:
        """
        Set z, which should not have been initialized in any way,
        to point to an array of n 0's.
        """
        cdef int i
        z[0] = <mpz_t *> PyMem_Malloc(sizeof(mpz_t) * n)
        if z[0] == <mpz_t *>0:
            raise MemoryError
        for i from 0 <= i < n:
            mpz_init_set_si(z[0][i], 0)
        return 0

    cdef public int mpzvec_clear(self, mpz_t *z, int n) except -1:
        """
        z is an allocated array of mpz_t elements.
        This function goes through and clears all those
        mpz_t elements.
        The caller should then delete the array of mpz_t elements.
        """
        cdef int i
        for i from 0 <= i < n:
            mpz_clear(z[i])
        return 0

    cdef public int mpz_height_vec(self, mpz_t H, mpz_t *v, int n) except -1:
        """
        Computes the max of the absolute values of the entries of the C array v of length n.
        """
        cdef int i
        cdef mpz_t x
        mpz_init(x)

        _sig_on
        mpz_set_si(H, 0)
        for i from 0 <= i < n:
            mpz_abs(x, v[i])
            if mpz_cmp(x, H) > 0:
                mpz_set(H, x)
        mpz_clear(x)
        _sig_off
        return 0


