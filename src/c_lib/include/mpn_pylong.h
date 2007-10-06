#ifndef MPN_PYLONG_H
#define MPN_PYLONG_H

#include <Python.h>
#include <gmp.h>

/************************************************************/

/* Python internals for pylong */

#include <longintrepr.h>
typedef int py_size_t; /* what python uses for ob_size */

/************************************************************/

/* mpn -> pylong conversion */

int mpn_pylong_size (mp_ptr up, mp_size_t un);

/* Assume digits points to a chunk of size size
 * where size >= mpn_pylong_size(up, un)
 */
void mpn_get_pylong (digit *digits, py_size_t size, mp_ptr up, mp_size_t un);

/************************************************************/

/* pylong -> mpn conversion */

mp_size_t mpn_size_from_pylong (digit *digits, py_size_t size);

/* Assume up points to a chunk of size un
 * where un >= mpn_size_from_pylong(digits, size)
 */
void mpn_set_pylong (mp_ptr up, mp_size_t un, digit *digits, py_size_t size);

/************************************************************/

/* Python hashing */

/*
 * for an mpz, this number has to be multiplied by the sign
 * also remember to catch -1 and map it to -2 !
 */

long mpn_pythonhash (mp_ptr up, mp_size_t un);

/************************************************************/

#endif
