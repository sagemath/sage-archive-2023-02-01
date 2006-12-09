"""
Basic RSA commands

Created 11-24-2005 by wdj. Last updated 12-02-2005.
"""

###########################################################################
#  Copyright (C) 2006 David Joyner <wdj@usna.edu> and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
###########################################################################

def rsa_encrypt(m, N, e):
    """
    INPUT:
        N -- typically a product of distinct primes.
        e -- element of Z/phiN Z.
        m -- integer relatively prime to e.

    OUTPUT:
        RSA encryption $m^e \pmod N$.

    EXAMPLES:
        sage: p = 53
        sage: q = 61
        sage: N = p*q
        sage: phiN = euler_phi(N)
        sage: R = IntegerModRing(int(phiN))
        sage: e = R(17)
        sage: rsa_encrypt(123,N,e)
        855

    AUTHOR: David Joyner (11-2005)
    """
    return (m**e)%N

def rsa_decrypt(c,N,e):
    """
    Inverse function to rsa_encrypt.

    INPUT:
        N -- typically a product of distinct primes.
        e -- in Z/phiN Z.
        c -- an integer relatively prime to e.

    EXAMPLES:
        sage: p = 53
        sage: q = 61
        sage: N = p*q
        sage: phiN = euler_phi(N)
        sage: R = IntegerModRing(int(phiN))
        sage: e = R(17)
        sage: rsa_decrypt(855,N,e)
        123

    AUTHOR:
        -- David Joyner (11-2005)
    """
    d = 1/e;
    return (c**d)%N

