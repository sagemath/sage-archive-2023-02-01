cdef extern from "flint/fmpz_vec.h":

    void _fmpz_vec_scalar_mod_fmpz(fmpz *res, fmpz *vec, long len, fmpz_t p)
