"""
FLINT Arithmetic Functions
"""
###########################################################################
#       Copyright (C) 2013 Fredrik Johansson <fredrik.johansson@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

include "sage/ext/interrupt.pxi"
include "fmpz.pxi"
include "fmpq_poly.pxi"

cdef extern from "flint/arith.h":
    void arith_number_of_partitions(fmpz_t x, unsigned long n)
    void arith_dedekind_sum(fmpq_t, fmpz_t, fmpz_t)

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

def number_of_partitions(unsigned long n):
    """
    Returns the number of partitions of the integer ``n``.

    EXAMPLES::

        sage: from sage.libs.flint.arith import number_of_partitions
        sage: number_of_partitions(3)
        3
        sage: number_of_partitions(10)
        42
        sage: number_of_partitions(40)
        37338
        sage: number_of_partitions(100)
        190569292
        sage: number_of_partitions(100000)
        27493510569775696512677516320986352688173429315980054758203125984302147328114964173055050741660736621590157844774296248940493063070200461792764493033510116079342457190155718943509725312466108452006369558934464248716828789832182345009262853831404597021307130674510624419227311238999702284408609370935531629697851569569892196108480158600569421098519

    TESTS::

        sage: n = 500 + randint(0,500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1500 + randint(0,1500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 100000000 + randint(0,100000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0  # long time
        True
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

def dedekind_sum(p, q):
    """
    Return the Dedekind sum `s(p, q)` where `p` and `q` are arbitrary integers.

    EXAMPLES::

        sage: from sage.libs.flint.arith import dedekind_sum
        sage: dedekind_sum(4, 5)
        -1/5
    """
    p = Integer(p)
    q = Integer(q)
    s = Rational(0)

    cdef fmpz_t p_fmpz, q_fmpz
    cdef fmpq_t s_fmpq

    fmpz_init(p_fmpz)
    fmpz_init(q_fmpz)
    fmpq_init(s_fmpq)

    fmpz_set_mpz(p_fmpz, (<Integer>p).value)
    fmpz_set_mpz(q_fmpz, (<Integer>q).value)

    arith_dedekind_sum(s_fmpq, p_fmpz, q_fmpz)

    fmpq_get_mpq((<Rational>s).value, s_fmpq)

    fmpz_clear(p_fmpz)
    fmpz_clear(q_fmpz)
    fmpq_clear(s_fmpq)

    return s

