#ifndef MPZ_PYLONG_H
#define MPZ_PYLONG_H

#include <Python.h>
#include <longintrepr.h>
#include <gmp.h>
#include "mpn_pylong.h"

/* mpz python hash */
long
mpz_pythonhash (mpz_srcptr z);

/* mpz -> pylong conversion */
PyObject *
mpz_get_pylong(mpz_srcptr z);

/* pylong -> mpz conversion */
int
mpz_set_pylong(mpz_ptr z, PyObject * ll);

#endif
