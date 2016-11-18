from sage.libs.arb.types cimport acb_t, acb_mat_t

cdef extern from "acb_mat.h":
    unsigned int acb_mat_nrows(acb_mat_t mat)
    unsigned int acb_mat_ncols(acb_mat_t mat)
    acb_t acb_mat_entry(acb_mat_t mat, unsigned long i, unsigned long j)

    void acb_mat_init(acb_mat_t mat, long r, long c)
    void acb_mat_clear(acb_mat_t mat)

    void acb_mat_set(acb_mat_t dest, const acb_mat_t src)

    void acb_mat_zero(acb_mat_t mat)
    void acb_mat_one(acb_mat_t mat)

    void acb_mat_add(acb_mat_t res, const acb_mat_t mat1,
                     const acb_mat_t mat2, long prec)
    void acb_mat_sub(acb_mat_t res, const acb_mat_t mat1,
                     const acb_mat_t mat2, long prec);
    void acb_mat_mul(acb_mat_t res, const acb_mat_t mat1,
                     const acb_mat_t mat2, long prec);

    void acb_mat_scalar_addmul_acb(acb_mat_t B, const acb_mat_t A,
                                   const acb_t c, long prec)
