# This file was (manually) generated from FLINT's fmpz_vec.h.
#*****************************************************************************
#        Copyright (C) 2010 William Hart
#        Copyright (C) 2010 Sebastian Pancratz
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef extern from "flint/fmpz_vec.h":
    #  Memory management  ******************************************************
    fmpz * _fmpz_vec_init(long len)
    void _fmpz_vec_clear(fmpz * vec, long len)

    #  Randomisation  **********************************************************
    void _fmpz_vec_randtest(fmpz * f, flint_rand_t state, long len, mp_bitcnt_t bits)
    void _fmpz_vec_randtest_unsigned(fmpz * f, flint_rand_t state, long len, mp_bitcnt_t bits)

    #  Norms  ******************************************************************
    long _fmpz_vec_max_bits(fmpz * vec, long len)
    long _fmpz_vec_max_bits_ref(fmpz * vec, long len)
    mp_size_t _fmpz_vec_max_limbs(fmpz * vec, long len)
    void _fmpz_vec_height(fmpz_t height, fmpz * vec, long len)
    long _fmpz_vec_height_index(fmpz * vec, long len)

    #  Input and output  *******************************************************
    int _fmpz_vec_fprint(FILE * file, fmpz * vec, long len)
    int _fmpz_vec_print(fmpz * vec, long len)
    int _fmpz_vec_fread(FILE * file, fmpz ** vec, long * len)
    int _fmpz_vec_read(fmpz ** vec, long * len)

    #  Conversions  ************************************************************
    void _fmpz_vec_set_nmod_vec(fmpz * res, mp_srcptr poly, long len, nmod_t mod)
    void _fmpz_vec_get_nmod_vec(mp_ptr res, fmpz * poly, long len, nmod_t mod)
    long _fmpz_vec_get_fft(mp_limb_t ** coeffs_f, fmpz * coeffs_m, long l, long length)
    void _fmpz_vec_set_fft(fmpz * coeffs_m, long length, mp_limb_t ** coeffs_f, long limbs, long sign)

    #  Assignment and basic manipulation  **************************************
    void _fmpz_vec_set(fmpz * vec1, fmpz * vec2, long len2)
    void _fmpz_vec_swap(fmpz * vec1, fmpz * vec2, long len2)
    void _fmpz_vec_zero(fmpz * vec, long len)
    void _fmpz_vec_neg(fmpz * vec1, fmpz * vec2, long len2)

    #  Comparison  *************************************************************
    int _fmpz_vec_equal(fmpz * vec1, fmpz * vec2, long len)
    int _fmpz_vec_is_zero(fmpz * vec, long len)

    # Sorting  *****************************************************************
    void _fmpz_vec_sort(fmpz * vec, long len)

    #  Addition  ***************************************************************
    void _fmpz_vec_add(fmpz * res, fmpz * vec1, fmpz * vec2, long len2)
    void _fmpz_vec_sub(fmpz * res, fmpz * vec1, fmpz * vec2, long len2)

    #  Scalar multiplication and division  *************************************
    void _fmpz_vec_scalar_mul_si(fmpz * vec1, fmpz * vec2, long len2, long c)
    void _fmpz_vec_scalar_mul_ui(fmpz * vec1, fmpz * vec2, long len2, unsigned long c)
    void _fmpz_vec_scalar_mul_fmpz(fmpz * vec1, fmpz * vec2, long len2, fmpz_t x)
    void _fmpz_vec_scalar_mul_2exp(fmpz * vec1, fmpz * vec2, long len2, unsigned long exp)
    void _fmpz_vec_scalar_divexact_fmpz(fmpz * vec1, fmpz * vec2, long len2, fmpz_t x)
    void _fmpz_vec_scalar_divexact_si(fmpz * vec1, fmpz * vec2, long len2, long c)
    void _fmpz_vec_scalar_divexact_ui(fmpz * vec1, fmpz * vec2, long len2, unsigned long c)
    void _fmpz_vec_scalar_fdiv_q_fmpz(fmpz * vec1, fmpz * vec2, long len2, fmpz_t c)
    void _fmpz_vec_scalar_fdiv_q_si(fmpz * vec1, fmpz * vec2, long len2, long c)
    void _fmpz_vec_scalar_fdiv_q_ui(fmpz * vec1, fmpz * vec2, long len2, unsigned long c)
    void _fmpz_vec_scalar_fdiv_q_2exp(fmpz * vec1, fmpz * vec2, long len2, unsigned long exp)
    void _fmpz_vec_scalar_tdiv_q_fmpz(fmpz * vec1, fmpz * vec2, long len2, fmpz_t c)
    void _fmpz_vec_scalar_tdiv_q_si(fmpz * vec1, fmpz * vec2, long len2, long c)
    void _fmpz_vec_scalar_tdiv_q_ui(fmpz * vec1, fmpz * vec2, long len2, unsigned long c)
    void _fmpz_vec_scalar_tdiv_q_2exp(fmpz * vec1, fmpz * vec2, long len2, unsigned long exp)
    void _fmpz_vec_scalar_addmul_si(fmpz * vec1, fmpz * vec2, long len2, long c)
    void _fmpz_vec_scalar_addmul_fmpz(fmpz * poly1, fmpz * poly2, long len2, fmpz_t x)
    void _fmpz_vec_scalar_addmul_si_2exp(fmpz * vec1, fmpz * vec2, long len2, long c, unsigned long exp)
    void _fmpz_vec_scalar_submul_si(fmpz * vec1, fmpz * vec2, long len2, long c)
    void _fmpz_vec_scalar_submul_fmpz(fmpz * vec1, fmpz * vec2, long len2, fmpz_t x)
    void _fmpz_vec_scalar_submul_si_2exp(fmpz * vec1, fmpz * vec2, long len2, long c, unsigned long exp)

    #  Vector sum and product  *************************************************
    void _fmpz_vec_sum(fmpz_t res, fmpz * vec, long len)
    void _fmpz_vec_prod(fmpz_t res, fmpz * vec, long len)

    #  Reduction mod p *********************************************************
    void _fmpz_vec_scalar_mod_fmpz(fmpz *res, fmpz *vec, long len, fmpz_t p)
    void _fmpz_vec_scalar_smod_fmpz(fmpz *res, fmpz *vec, long len, fmpz_t p)

    #  Gaussian content  *******************************************************
    void _fmpz_vec_content(fmpz_t res, fmpz * vec, long len)
    void _fmpz_vec_lcm(fmpz_t res, fmpz * vec, long len)
