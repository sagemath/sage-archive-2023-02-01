/*
General C (.h) code this is useful to include in any pyrex module.

These are mostly things that can't be done in Pyrex.
*/

/******************************************************************************
       Copyright (C) 2006 William Stein <wstein@gmail.com>

  Distributed under the terms of the GNU General Public License (GPL), Version 2.

  The full text of the GPL is available at:
                  http://www.gnu.org/licenses/

******************************************************************************/


/*****************************************
          For PARI
 *****************************************/

#define set_gel(x, n, z)  gel(x,n)=z;


/******************************************
 Some macros exported for Pyrex in cdefs.pxi
 ****************************************/

/* Tests whether zzz_obj is of type zzz_type. The zzz_type must be a built-in
   or extension type. This is just a C++-compatible wrapper for
   PyObject_TypeCheck. */
#define PY_TYPE_CHECK(zzz_obj, zzz_type) \
    (PyObject_TypeCheck((PyObject*)(zzz_obj), (PyTypeObject*)(zzz_type)))

/* Returns the type field of a python object, cast to void*. The returned
   value should only be used as an opaque object e.g. for type comparisons. */
#define PY_TYPE(zzz_obj) ((void*)((zzz_obj)->ob_type))

/* Constructs a new object of type zzz_type by calling tp_new directly,
   with no arguments. */
#define PY_NEW(zzz_type) \
    (((PyTypeObject*)(zzz_type))->tp_new((zzz_type), global_empty_tuple, NULL))

/* Resets the tp_new slot of zzz_type1 to point to the tp_new slot of
   zzz_type2. This is used in SAGE to speed up Pyrex's boilerplate object
   construction code by skipping irrelevant base class tp_new methods. */
#define PY_SET_TP_NEW(zzz_type1, zzz_type2) \
    (((PyTypeObject*)zzz_type1)->tp_new = ((PyTypeObject*)zzz_type2)->tp_new)

extern PyObject* global_empty_tuple;    /* from stdsage.c */
