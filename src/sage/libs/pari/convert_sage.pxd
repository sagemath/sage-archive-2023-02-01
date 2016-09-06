from .types cimport GEN
from .gen cimport gen
from sage.libs.flint.types cimport fmpz_t, fmpz_mat_t
from sage.libs.gmp.types cimport mpz_t, mpq_t, mpz_ptr, mpq_ptr

cdef gen new_gen_from_mpz_t(mpz_t value)
cdef GEN _new_GEN_from_mpz_t(mpz_t value)
cdef gen new_gen_from_mpq_t(mpq_t value)
cdef GEN _new_GEN_from_mpq_t(mpq_t value)
cdef GEN _new_GEN_from_fmpz_t(fmpz_t value)
cdef gen new_gen_from_padic(long ordp, long relprec, mpz_t prime, mpz_t p_pow, mpz_t unit)
cdef GEN _new_GEN_from_fmpz_mat_t(fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc)
cdef GEN _new_GEN_from_fmpz_mat_t_rotate90(fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc)
cdef gen integer_matrix(fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf)
cdef GEN _new_GEN_from_mpq_t_matrix(mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)
cdef gen rational_matrix(mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)
cdef void INT_to_mpz(mpz_ptr value, GEN g)
cdef void INTFRAC_to_mpq(mpq_ptr value, GEN g)
