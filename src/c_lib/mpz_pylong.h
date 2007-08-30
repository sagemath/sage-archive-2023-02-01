#ifndef MPZ_PYLONG_H
#define MPZ_PYLONG_H

#include <Python.h>
#include <gmp.h>


#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

/* mpz -> pylong conversion */
EXTERN PyObject * mpz_get_pylong(mpz_srcptr z);

/* pylong -> mpz conversion */
EXTERN int mpz_set_pylong(mpz_ptr z, PyObject * ll);

/* mpz python hash */
EXTERN long mpz_pythonhash (mpz_srcptr z);

#endif
