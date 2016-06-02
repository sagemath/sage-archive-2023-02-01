r"""
Fast binary operations for basic types
"""

from .types cimport mpz_t, mpq_t
from .mpz cimport mpz_set, mpz_sgn, mpz_sub, mpz_add, mpz_divexact, mpz_mul, mpz_neg
from .mpq cimport mpq_init, mpq_canonicalize, mpq_add, mpq_sub, mpq_mul, mpq_numref, mpq_denref

cdef inline void mpq_add_z(mpq_t res, mpq_t op1, mpz_t op2):
    mpz_mul(mpq_numref(res), mpq_denref(op1), op2)
    mpz_add(mpq_numref(res), mpq_numref(res), mpq_numref(op1))
    mpz_set(mpq_denref(res), mpq_denref(op1))

cdef inline void mpq_sub_z(mpq_t res, mpq_t op1, mpz_t op2):
    mpz_mul(mpq_numref(res), mpq_denref(op1), op2)
    mpz_sub(mpq_numref(res), mpq_numref(op1), mpq_numref(res))
    mpz_set(mpq_denref(res), mpq_denref(op1))

cdef inline void mpq_mul_z(mpq_t res, mpq_t op1, mpz_t op2):
    # (A/B) * C = (A/C) * B
    mpz_set(mpq_numref(res), op2)
    mpz_set(mpq_denref(res), mpq_denref(op1))
    mpq_canonicalize(res)
    mpz_mul(mpq_numref(res), mpq_numref(res), mpq_numref(op1))

cdef inline void mpq_div_z(mpq_t res, mpq_t op1, mpz_t op2):
    # (A/B) / C = (A/C) / B
    mpz_set(mpq_numref(res), mpq_numref(op1))
    mpz_set(mpq_denref(res), op2)
    mpq_canonicalize(res)
    mpz_mul(mpq_denref(res), mpq_denref(res), mpq_denref(op1))
    if mpz_sgn(mpq_denref(res)) == -1:
        mpz_neg(mpq_numref(res), mpq_numref(res))
        mpz_neg(mpq_denref(res), mpq_denref(res))



