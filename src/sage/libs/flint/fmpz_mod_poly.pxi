include "fmpz.pxi"

cdef extern from "flint/fmpz_mod_poly.h":

    int _fmpz_mod_poly_invmod(fmpz *A, fmpz *B, long lenB, fmpz *P, long lenP, fmpz_t p)
