/******************************************************************************
       Copyright (C) 2006 William Stein <wstein@gmail.com>
                     2006 Martin Albrecht <malb@informatik.uni-bremen.de>

   Distributed under the terms of the GNU General Public License (GPL)
   as published by the Free Software Foundation; either version 2 of
   the License, or (at your option) any later version.
                   http://www.gnu.org/licenses/

******************************************************************************/

/**
 * @file stdsage.h
 *
 * @author William Stein <wstein@gmail.com>
 * @auhtor Martin Albrecht <malb@informatik.uni-bremen.de>
 *
 * @brief General C (.h) code this is useful to include in any pyrex module.
 *
 * Put
 @verbatim

  include 'relative/path/to/stdsage.pxi'

 @endverbatim
 *
 * at the top of your Pyrex file.
 *
 * These are mostly things that can't be done in Pyrex.
 */

#ifndef STDSAGE_H
#define STDSAGE_H

#include "Python.h"
#include "memory.h"

/* Building with this not commented out causes
   serious problems on RHEL5 64-bit for Kiran Kedlaya... i.e., it doesn't work. */
/* #include "ccobject.h" */

#ifdef __cplusplus
extern "C" {
#endif


/*****************************************
   PARI array element assignment
 *****************************************/
#define set_gel(x, n, z)         (gel(x,n) = z)
#define set_gmael(x, i, j, z)    (gmael(x,i,j) = z)
#define set_gcoeff(x, i, j, z)   (gcoeff(x,i,j) = z)


/******************************************
 Some macros exported for Pyrex in cdefs.pxi
 ****************************************/

/** Tests whether zzz_obj is of type zzz_type. The zzz_type must be a
 * built-in or extension type. This is just a C++-compatible wrapper
 * for PyObject_TypeCheck.
 */
#define PY_TYPE_CHECK(zzz_obj, zzz_type) \
    (PyObject_TypeCheck((PyObject*)(zzz_obj), (PyTypeObject*)(zzz_type)))

/** Tests whether zzz_obj is exactly of type zzz_type. The zzz_type must be a
  * built-in or extension type.
  */
#define PY_TYPE_CHECK_EXACT(zzz_obj, zzz_type) \
  ((PyTypeObject*)PY_TYPE(zzz_obj) == (PyTypeObject*)(zzz_type))

/** Returns the type field of a python object, cast to void*. The
 *  returned value should only be used as an opaque object e.g. for
 *  type comparisons.
 */
#define PY_TYPE(zzz_obj) ((void*)((zzz_obj)->ob_type))

/** Constructs a new object of type zzz_type by calling tp_new
 *  directly, with no arguments.
 */

#define PY_NEW(zzz_type) \
    (((PyTypeObject*)(zzz_type))->tp_new((PyTypeObject*)(zzz_type), global_empty_tuple, NULL))


  /** Constructs a new object of type the same type as zzz_obj by calling tp_new
   *  directly, with no arguments.
   */

#define PY_NEW_SAME_TYPE(zzz_obj) \
  PY_NEW(PY_TYPE(zzz_obj))

/** Resets the tp_new slot of zzz_type1 to point to the tp_new slot of
 *  zzz_type2. This is used in SAGE to speed up Pyrex's boilerplate
 *  object construction code by skipping irrelevant base class tp_new
 *  methods.
 */
#define PY_SET_TP_NEW(zzz_type1, zzz_type2) \
    (((PyTypeObject*)zzz_type1)->tp_new = ((PyTypeObject*)zzz_type2)->tp_new)


/**
 * Tests whether the given object has a python dictionary.
 */
#define HAS_DICTIONARY(zzz_obj) \
    (((PyObject*)(zzz_obj))->ob_type->tp_dictoffset != NULL)

/**
 * Very very unsafe access to the list of pointers to PyObject*'s
 * underlying a list / sequence.  This does error checking of any kind
 * -- make damn sure you hand it a list or sequence!
 */
#define FAST_SEQ_UNSAFE(zzz_obj) \
    PySequence_Fast_ITEMS(PySequence_Fast(zzz_obj, "expected sequence type"))

/** Returns the type field of a python object, cast to void*. The
 *  returned value should only be used as an opaque object e.g. for
 *  type comparisons.
 */
#define PY_IS_NUMERIC(zzz_obj) \
     (PyInt_Check(zzz_obj) ||  PyBool_Check(zzz_obj) || PyLong_Check(zzz_obj) || \
       PyFloat_Check(zzz_obj) || PyComplex_Check(zzz_obj))


/** This is exactly the same as isinstance (and does return a Python
 *  bool), but the second argument must be a C-extension type -- so it
 *  can't be a Python class or a list.  If you just want an int return
 *  value, i.e., aren't going to pass this back to Python, just use
 *  PY_TYPE_CHECK.
 */
#define IS_INSTANCE(zzz_obj, zzz_type) \
    PyBool_FromLong(PY_TYPE_CHECK(zzz_obj, zzz_type))


/**
 * A global empty python tuple object. This is used to speed up some
 * python API calls where we want to avoid constructing a tuple every
 * time.
 */

extern PyObject* global_empty_tuple;


/**
 * Initialisation of signal handlers, global variables, etc. Called
 * exactly once at Sage start-up.
 */
void init_csage(void);



/**
 * a handy macro to be placed at the top of a function definition
 * below the variable declarations to ensure a function is called once
 * at maximum.
 */
#define _CALLED_ONLY_ONCE static int ncalls = 0; if (ncalls>0) return; else ncalls++

#ifdef __cplusplus
}
#endif

#endif
