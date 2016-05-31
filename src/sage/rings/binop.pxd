from sage.libs.gmp.types cimport mpz_t, mpq_t

cdef inline void mpq_add_z(mpq_t res, mpq_t op1, mpz_t op2)
cdef inline void mpq_sub_z(mpq_t res, mpq_t op1, mpz_t op2)
cdef void mpq_mul_z(mpq_t res, mpq_t op1, mpz_t op2)
cdef void mpq_div_z(mpq_t res, mpq_t op1, mpz_t op2)
