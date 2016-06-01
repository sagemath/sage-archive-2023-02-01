r"""
Fast binary operations for basic types
"""

from .types cimport mpz_t, mpq_t
from .mpz cimport mpz_init, mpz_set, mpz_sgn, mpz_sub, mpz_add, mpz_clear, mpz_divexact, mpz_gcd, mpz_mul, mpz_neg
from .mpq cimport mpq_init, mpq_clear, mpq_set_z, mpq_add, mpq_sub, mpq_mul, mpq_numref, mpq_denref

cdef inline void mpq_add_z(mpq_t res, mpq_t op1, mpz_t op2):
    mpz_mul(mpq_numref(res), mpq_denref(op1), op2)
    mpz_add(mpq_numref(res), mpq_numref(res), mpq_numref(op1))
    mpz_set(mpq_denref(res), mpq_denref(op1))

cdef inline void mpq_sub_z(mpq_t res, mpq_t op1, mpz_t op2):
    mpz_mul(mpq_numref(res), mpq_denref(op1), op2)
    mpz_sub(mpq_numref(res), mpq_numref(op1), mpq_numref(res))
    mpz_set(mpq_denref(res), mpq_denref(op1))

cdef inline void mpq_mul_z(mpq_t res, mpq_t op1, mpz_t op2):
    # op1.num/op1.den * c = op1.num(c/gcd(op1.den,c))/ (op1.den/gcd(op1.den,c))
    cdef mpz_t z
    mpz_init(z)
    mpz_gcd(z, mpq_denref(op1), op2)
    mpz_divexact(mpq_numref(res), op2, z)
    mpz_mul(mpq_numref(res), mpq_numref(res), mpq_numref(op1))
    mpz_divexact(mpq_denref(res), mpq_denref(op1), z)
    mpz_clear(z)

cdef inline void mpq_div_z(mpq_t res, mpq_t op1, mpz_t op2):
    # op1.num/op1.den / c = (op1.num/gcd(op1.num, c)) / (op1.den * c/gcd(op1.num, c))
    cdef mpz_t z
    mpz_init(z)
    mpz_gcd(z, mpq_numref(op1), op2)
    mpz_divexact(mpq_denref(res), op2, z)
    mpz_mul(mpq_denref(res), mpq_denref(res), mpq_denref(op1))
    mpz_divexact(mpq_numref(res), mpq_numref(op1), z)
    if mpz_sgn(mpq_denref(res)) == -1:
        mpz_neg(mpq_numref(res), mpq_numref(res))
        mpz_neg(mpq_denref(res), mpq_denref(res))
    mpz_clear(z)

