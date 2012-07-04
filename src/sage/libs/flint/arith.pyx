include "../../ext/interrupt.pxi"
include "fmpz.pxi"

cdef extern from "flint/arith.h":
    void arith_number_of_partitions(fmpz_t x, unsigned long n)

from sage.rings.integer cimport Integer


def number_of_partitions(unsigned long n):
    """
    Returns the number of partitions of the integer n.

    EXAMPLES::

        sage: from sage.libs.flint.arith import number_of_partitions
        sage: number_of_partitions(10)
        42
    """
    cdef fmpz_t ans_fmpz
    cdef Integer ans

    fmpz_init(ans_fmpz)

    if n > 1000:
        sig_on()

    arith_number_of_partitions(ans_fmpz, n)

    if n > 1000:
        sig_off()

    ans = Integer(0)
    fmpz_get_mpz(ans.value, ans_fmpz)
    fmpz_clear(ans_fmpz)
    return ans

