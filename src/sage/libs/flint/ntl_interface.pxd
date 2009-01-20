include "../../ext/cdefs.pxi"
include "fmpz.pxi"

from sage.libs.flint.fmpz_poly cimport fmpz_poly_t
from sage.libs.ntl.ntl_ZZ_decl cimport ZZ_c
from sage.libs.ntl.ntl_ZZX_decl cimport ZZX_c

cdef extern from "FLINT/NTL-interface.h":
    unsigned long ZZ_limbs(ZZ_c z)

    void fmpz_poly_to_ZZX(ZZX_c output, fmpz_poly_t poly)
    void ZZX_to_fmpz_poly(fmpz_poly_t output, ZZX_c poly)

    void fmpz_to_mpz(mpz_t res, fmpz_t f)
    void ZZ_to_fmpz(fmpz_t output, ZZ_c z)

