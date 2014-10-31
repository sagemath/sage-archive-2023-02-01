/*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
#*****************************************************************************/

/*  Author:  Joel B. Mohler <joel@kiwistrawberry.us> (2007-06-17) */

#include "ZZ_pylong.h"
#include "ntl_wrap.h"
extern "C" {
    #include "mpz_pylong.h"
}

using namespace NTL;

/* ZZ -> pylong conversion */
PyObject * ZZ_get_pylong(ZZ &z)
{
    mpz_t temp;
    PyObject *val;
    mpz_init(temp);
    ZZ_to_mpz( temp, &z );
    val = mpz_get_pylong( temp );
    mpz_clear( temp );
    return val;
}

/* pylong -> ZZ conversion */
int ZZ_set_pylong(ZZ &z, PyObject * ll)
{
    mpz_t temp;
    mpz_init(temp);
    mpz_set_pylong( temp, ll );
    mpz_to_ZZ( &z, temp );
    mpz_clear( temp );
    return 0;
}
