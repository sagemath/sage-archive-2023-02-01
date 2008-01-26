/*
 * $Id: rpymodule.c 384 2007-11-30 00:40:17Z warnes $
 * Implementation of the module '_rpy' and the 'Robj' type.
 */

/* ***** BEGIN LICENSE BLOCK *****
 * Version: MPL 1.1/GPL 2.0/LGPL 2.1
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * The Original Code is the RPy python module.
 *
 * The Initial Developer of the Original Code is Walter Moreira.
 * Portions created by the Initial Developer are Copyright (C) 2002
 * the Initial Developer. All Rights Reserved.
 *
 * Contributor(s):
 *    Gregory R. Warnes <greg@warnes.net> (Maintainer)
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU General Public License Version 2 or later (the "GPL"), or
 * the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
 * in which case the provisions of the GPL or the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of either the GPL or the LGPL, and not to allow others to
 * use your version of this file under the terms of the MPL, indicate your
 * decision by deleting the provisions above and replace them with the notice
 * and other provisions required by the GPL or the LGPL. If you do not delete
 * the provisions above, a recipient may use your version of this file under
 * the terms of any one of the MPL, the GPL or the LGPL.
 *
 * ***** END LICENSE BLOCK ***** */

#define CSTACK_DEFNS
#include "RPy.h"

#define NONAMELESSUNION
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Flag indicating whether Numpy/Numeric is available in this session
 *
 * This is necessary since Numpy/Numeric may not available at run time, even if
 * it was available at compile time.
*/
static int use_numeric=0;


/* Local function definitions */
DL_EXPORT(void) INIT_RPY(void);           /* Module initializer */
static PyObject *r_init(PyObject *self,   /* Class initializer */
                        PyObject *args);
static PyObject *r_cleanup(void);         /* Clean up R & release resources */

#ifdef _WIN32
static void init_embedded_win32( void );
#endif

/* Global objects */
static SEXP get_item;
static SEXP set_item;
static SEXP length;
static SEXP aperm;
static PyObject *class_table;
static PyObject *proc_table;
static int default_mode;
static PyObject *r_lock;
PyObject *RPy_Exception;
PyObject *RPy_TypeConversionException;
PyObject *RPy_RException;

static char RHOME[BUFSIZ];
static char RVERSION[BUFSIZ];
static char RVER[BUFSIZ];
static char RUSER[BUFSIZ];

/* Global interpreter */
PyInterpreterState *my_interp;

/* Signal whether R is running interactively */
int R_interact;

/* RPy namespace */
PyObject *rpy;
PyObject *rpy_dict;


#ifdef WITH_NUMERIC
static PyObject *Py_transpose;
#endif

/* Global list to protect R objects from garbage collection */
/* This is inspired in $R_SRC/src/main/memory.c */
static SEXP R_References;

static SEXP
RecursiveRelease(SEXP obj, SEXP list)
{
  if (!isNull(list)) {
    if (obj == CAR(list))
      return CDR(list);
    else
      SETCDR(list, RecursiveRelease(obj, CDR(list)));
  }
  return list;
}

/* Robj methods. Following xxmodule.c from Python distro. */

static void
Robj_dealloc(RobjObject *self)
{
  /* Remove the object from the list of protected objects */
  R_References = RecursiveRelease(self->R_obj, R_References);
  SET_SYMVALUE(install("R.References"), R_References);

  PyObject_Del(self);
}

RobjObject *
Robj_new(SEXP robj, int conversion)
{
  RobjObject *self;
  self = PyObject_New(RobjObject, &Robj_Type);
  if (!self)
    return NULL;

  if (!robj)
    return NULL;

  /* Protect the R object */
  R_References = CONS(robj, R_References);
  SET_SYMVALUE(install("R.References"), R_References);

  self->R_obj = robj;
  self->conversion = conversion;
  return self;
}

#ifndef PRE_2_2
static PyObject *
Robj_tpnew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyObject *self;

  self = type->tp_alloc(type, 0);
  return self;
}
#endif

/* Type conversion routines. See documentation for details */

/* These are auxiliaries for a state machine for converting Python
   list to the coarsest R vector type */
#define ANY_T 0
#define BOOL_T 1
#define INT_T 2
#define FLOAT_T 3
#define COMPLEX_T 4
#define STRING_T 5
#define ROBJ_T 6

static int
type_to_int(PyObject *obj)
{
  if (PyBool_Check(obj))
    return BOOL_T;
  else if (PyInt_Check(obj))
    return INT_T;
  else if (PyFloat_Check(obj))
    return FLOAT_T;
  else if (PyComplex_Check(obj))
    return COMPLEX_T;
  else if (PyNumber_Check(obj))
    return ANY_T;
  else if (PyString_Check(obj))
    return STRING_T;
  else if (PyUnicode_Check(obj))
    return STRING_T;
  else if (Robj_Check(obj))
    return ROBJ_T;
  else
    return ANY_T;
}

/* Make a R list or vector from a Python sequence */
static SEXP
seq_to_R(PyObject *obj)
{
  PyObject *it;
  SEXP robj, rit;
  int i, len, state;

  /* This matrix defines what mode a vector should take given what
     it already contains and a new item

     E.g. Row 0 indicates that if we've seen an any, the vector will
     always remain an any.  Row 3 indicates that if we've seen a
     float, then seeing an boolean, integer, or float will preserve
     the vector as a float vector, while seeing a string or an Robj will
     convert it into an any vector.
  */
  int fsm[7][7] = {
    {0, 0, 0, 0, 0, 0, 0}, // any
    {0, 1, 2, 3, 4, 0, 0}, // bool
    {0, 2, 2, 3, 4, 0, 0}, // int
    {0, 3, 3, 3, 4, 0, 0}, // float
    {0, 4, 4, 4, 4, 0, 0}, // complex
    {0, 0, 0, 0, 0, 5, 0}, // string
    {0, 0, 0, 0, 0, 0, 6}  // RObj
  };

  len = PySequence_Length(obj);
  if (len == 0)
    return R_NilValue;

  PROTECT(robj = NEW_LIST(len));

  state = -1;
  for (i=0; i<len; i++)
    {
      it = PySequence_GetItem(obj, i);
      if (it==NULL)
	{
	  UNPROTECT(1);
	  return NULL;
	}

      if (state < 0)
	state = type_to_int(it);
      else
	state = fsm[state][type_to_int(it)];

      rit = to_Robj(it);
      if (rit==NULL || PyErr_Occurred()	)
	{
	  Py_XDECREF(it);
	  UNPROTECT(1);
	  return NULL;
	}

      SET_VECTOR_ELT(robj, i, rit);
      Py_XDECREF(it);
    }

  switch(state)
    {
    case INT_T:
      robj = AS_INTEGER(robj);
      break;
    case BOOL_T:
      robj = AS_LOGICAL(robj);
      break;
    case FLOAT_T:
      robj = AS_NUMERIC(robj);
      break;
    case COMPLEX_T:
      robj = AS_COMPLEX(robj);
      break;
    case STRING_T:
      robj = AS_CHARACTER(robj);
      break;
    default:;
      /* Otherwise, it's either an ANY_T or ROBJ_T - we want ANY */
    }

  UNPROTECT(1);
  return robj;
}

/* Make a R named list or vector from a Python dictionary */
static SEXP
dict_to_R(PyObject *obj)
{
  int len;
  PyObject *keys, *values;
  SEXP robj, names;

  len = PyMapping_Length(obj);
  if (len == 0)
    return R_NilValue;

  /* If 'keys' succeed and 'values' fails this leaks */
  if (!(keys = PyMapping_Keys(obj)))
    return NULL;
  if (!(values = PyMapping_Values(obj)))
    return NULL;

  robj = seq_to_R(values);
  names = seq_to_R(keys);
  if ( robj==NULL || robj==NULL )
    {
      Py_DECREF(keys);
      Py_DECREF(values);
      return NULL;
    }

  PROTECT(robj);
  SET_NAMES(robj, names);

  Py_DECREF(keys);
  Py_DECREF(values);
  UNPROTECT(1);

  return robj;
}

#ifdef WITH_NUMERIC
/* Convert a Numpy/Numeric array to a R array */
SEXP
to_Rarray(PyObject *o)
{
  PyObject *pytl, *nobj;
  PyArrayObject *obj;
  SEXP Rdims, tRdims, Rarray, e;
  n_intp *dims;
  int i;
  int type;
  long size;

  obj = (PyArrayObject *)o;
  dims = obj->dimensions;
  type = obj->descr->type_num;
  size = PyArray_Size( (PyObject*) obj);

  /* Handle a vector without dimensions, just length */
  if(obj->nd==0)
    {
      PROTECT(Rdims = allocVector(INTSXP, 1));
      PROTECT(tRdims = allocVector(INTSXP, 1));
      INTEGER(Rdims)[0] = size;
      INTEGER(tRdims)[0] = size;
    }
  else
    {
      PROTECT(Rdims = allocVector(INTSXP, obj->nd));
      PROTECT(tRdims = allocVector(INTSXP, obj->nd));

      for (i=0; i<obj->nd; i++)
	{
	  if (dims[i] == 0)
	    {
	      UNPROTECT(2);
	      return R_NilValue;
	    }
	  INTEGER(Rdims)[i] = dims[(obj->nd)-i-1];
	  INTEGER(tRdims)[i] = (obj->nd)-i;
	}
    }

    switch(type)
    {

    /*******************/
    /* String Variants */
    /*******************/
    /* TODO: Add proper handling of NumPy character arrays.
             The following code DOES NOT WORK:

	     #if WITH_NUMERIC==3
	     case PyArray_UNICODE:
	     case PyArray_STRING:
	     case PyArray_CHAR:
	       obj = (PyArrayObject *)PyArray_ContiguousFromObject((PyObject *)obj,
							  PyArray_STRING, 0, 0);
	     #endif

      The problem is that the PyArray call throws an exception,
      presumably because we haven't given a width specifier.

      NumPy strings are fixed-width, and may not be null terminated.  R only handles
      null terminated (varying width) strings.  We need a separate
      code path to handle this, as it requires quite different
      handling than the numeric arrays dealt with below.
    */


    /******************************************/
    /* All complex to (double,double) complex */
    /******************************************/

#if WITH_NUMERIC==1       /* Numeric */
    case PyArray_CFLOAT:
    case PyArray_CDOUBLE:
#else                        /* NumPy */
    case PyArray_COMPLEX64:
    case PyArray_COMPLEX128:
#endif
      obj = (PyArrayObject *)PyArray_ContiguousFromObject((PyObject *)obj,
							  PyArray_CDOUBLE, 0, 0);
      break;


    /**********************************************************************************/
    /* Convert all integers to platform integer (except 64 bit int on 32 bit platforms) */
    /************************************************************************************/

#if WITH_NUMERIC==1      /* Numeric */
    case PyArray_UBYTE:
    case PyArray_SBYTE:
    case PyArray_SHORT:
    case PyArray_INT:
    case PyArray_LONG:
      obj = (PyArrayObject *)PyArray_ContiguousFromObject((PyObject *)obj,
							  PyArray_INT, 0, 0);
      break;
#else                    /* NumPy */
    case PyArray_BOOL:
    case PyArray_INT8:
    case PyArray_UINT8:
    case PyArray_INT16:
    case PyArray_UINT16:
    case PyArray_INT32:
    case PyArray_UINT32:
#if PyArray_INT==PyArray_INT64 /* 64 bit platform */
    case PyArray_INT64:
    case PyArray_UINT64:
#else
      obj = (PyArrayObject *)PyArray_ContiguousFromObject((PyObject *)obj,
							  PyArray_INT, 0, 0);
      break;
#endif
#endif

    /**************************************************/
    /* All floats (and over-sized integers) to double */
    /**************************************************/
#if WITH_NUMERIC==1       /* Numeric */
    case PyArray_FLOAT:
    case PyArray_DOUBLE:
#else                     /* NumPy */
    case PyArray_FLOAT32:
    case PyArray_FLOAT64:
#if PyArray_INT!=PyArray_INT64 /* 32 bit platform */
    case PyArray_INT64:
    case PyArray_UINT64:
#endif
#endif
      obj = (PyArrayObject *)PyArray_ContiguousFromObject((PyObject *)obj,
							  PyArray_DOUBLE, 0, 0);
      break;

    default:
      UNPROTECT(2);
      PyErr_Format(RPy_TypeConversionException,
		   "Numeric/NumPy arrays containing %s are not supported.",
		   obj->ob_type->tp_name);
      return R_NilValue;
      break;
    }


  pytl = Py_BuildValue("[i]", size);
  nobj = PyArray_Reshape(obj, pytl);
  Py_XDECREF(pytl);
  Py_XDECREF(obj);

  if (nobj == NULL)
    {
      UNPROTECT(2);
      return R_NilValue;
    }


  PROTECT(Rarray = seq_to_R(nobj));
  if (Rarray == NULL)
    {
      UNPROTECT(3);
      return R_NilValue;
    }


  Py_XDECREF(nobj);
  SET_DIM(Rarray, Rdims);

  PROTECT(e = allocVector(LANGSXP, 3));
  SETCAR(e, aperm);
  SETCAR(CDR(e), Rarray);
  SETCAR(CDR(CDR(e)), tRdims);
  PROTECT(Rarray = do_eval_expr(e));

  UNPROTECT(5);
  return Rarray;
}
#endif

/* Convert a Python object to a R object. An Robj is passed w/o
 * modifications, an object which provides a '.as_r()' method, is
 * passed as the result of that method */
SEXP
to_Robj(PyObject *obj)
{
  SEXP robj;
  Py_complex c;
  PyObject *to_r_meth;
  PyObject *tempObj;
  int do_decref = 0;

  if (obj==NULL)
    return NULL;

  if (obj == Py_None) {
    return R_NilValue;
  }

  to_r_meth = PyObject_GetAttrString(obj, "as_r");
  if (to_r_meth) {
    obj = PyObject_CallObject(to_r_meth, NULL);
    Py_DECREF(to_r_meth);
    if (obj==NULL)
      return NULL;
    do_decref = 1;
  }
  PyErr_Clear();

  to_r_meth = PyObject_GetAttrString(obj, "_rpy_");
  if (to_r_meth) {
    obj = PyObject_CallObject(to_r_meth, NULL);
    Py_DECREF(to_r_meth);
    if (obj==NULL)
      return NULL;
    do_decref = 1;
  }
  PyErr_Clear();


  if (Robj_Check(obj))
    {
      PROTECT(robj = ((RobjObject *)obj)->R_obj);
    }
  else if (PyBool_Check(obj))
    {
      PROTECT(robj = NEW_LOGICAL(1));
      LOGICAL_DATA(robj)[0] = (Py_True==obj);
    }
  else if (PyInt_Check(obj))
    {
      PROTECT(robj = NEW_INTEGER(1));
      INTEGER_DATA(robj)[0] = (int) PyInt_AsLong(obj);
    }
  else if (PyFloat_Check(obj))
    {
      PROTECT(robj = NEW_NUMERIC(1));
      NUMERIC_DATA(robj)[0] = PyFloat_AsDouble(obj);
    }
  else if (PyComplex_Check(obj))
    {
      PROTECT(robj = NEW_COMPLEX(1));
      c = PyComplex_AsCComplex(obj);
      COMPLEX_DATA(robj)[0].r = c.real;
      COMPLEX_DATA(robj)[0].i = c.imag;
    }
  else if (PyUnicode_Check(obj))
    {
      /** Handle Unicode Strings.
       *
       * Ideally:  Python Unicode -> R Unicode,
       *
       * Unfortunately, the R documentation is not forthcoming on how
       * to accomplish this
       *
       * So, for the moment:
       *     python Unicode -> Python ASCII -> ordinary string -> R string
       *
       */
      PROTECT(robj = NEW_STRING(1));
      SET_STRING_ELT(robj, 0, COPY_TO_USER_STRING(PyString_AsString(PyUnicode_AsASCIIString(obj))));
    }
  else if (PyString_Check(obj))
    {
      PROTECT(robj = NEW_STRING(1));
      SET_STRING_ELT(robj, 0, COPY_TO_USER_STRING(PyString_AsString(obj)));
    }
#ifdef WITH_NUMERIC
  else if (use_numeric && PyArray_Check(obj))
    {
      PROTECT(robj = to_Rarray(obj));
    }
#endif
  else if ((PySequence_Check(obj)) &&
           (PySequence_Size(obj) >= 0))
    {
      PROTECT(robj = seq_to_R(obj));      /* No labels */
    }
  else if ((PyMapping_Check(obj)) &&
            (PyMapping_Size(obj) >= 0))
    {
      PROTECT(robj = dict_to_R(obj));
    }
  else if (PyNumber_Check(obj)) /* generic number interface */
    {
      tempObj = PyNumber_Float(obj);
      if(!tempObj) goto error;

      PROTECT(robj = NEW_NUMERIC(1));
      NUMERIC_DATA(robj)[0] = PyFloat_AsDouble(tempObj);
      Py_DECREF(tempObj);
    }
  else
    {
    error:
      PyErr_Format(RPy_TypeConversionException,
		   "cannot convert from type '%s'",
                   obj->ob_type->tp_name);
      PROTECT(robj = NULL);    /* Protected to avoid stack inbalance */
    }

  if (do_decref)
    {
      Py_DECREF(obj);
    }
  UNPROTECT(1);
  return robj;
}

/* Convert a R named vector or list to a Python dictionary */
static PyObject *
to_PyDict(PyObject *obj, SEXP names)
{
  int len, i;
  PyObject *it, *dict;
  const char *name;

  if ((len = PySequence_Length(obj)) < 0)
    return NULL;

  dict = PyDict_New();
  for (i=0; i<len; i++) {
    if (!(it = PyList_GetItem(obj, i)))

      return NULL;
    name = CHAR(STRING_ELT(names, i));
    if ((PyDict_SetItemString(dict, name, it)) < 0) {
      return NULL;
    }
  }

  return dict;
}

/* We need to transpose the list because R makes array by the
 * fastest index */
static PyObject *
ltranspose(PyObject *list, int *dims, int *strides,
             int pos, int shift, int len)
{
  PyObject *nl, *it;
  int i;

  if (!(nl = PyList_New(dims[pos])))
    return NULL;

  if (pos == len-1) {
    for (i=0; i<dims[pos]; i++) {
      if (!(it = PyList_GetItem(list, i*strides[pos]+shift)))
        return NULL;
      Py_INCREF(it);
      if (PyList_SetItem(nl, i, it) < 0)
        return NULL;
    }
    return nl;
  }

  for (i=0; i<dims[pos]; i++) {
    if (!(it = ltranspose(list, dims, strides, pos+1, shift, len)))
      return NULL;
    if (PyList_SetItem(nl, i, it) < 0)
      return NULL;
    shift += strides[pos];
  }

  return nl;
}

/* Convert a Python list to a Python array (in the form of
 * list of lists of ...) */
static PyObject *
to_PyArray(PyObject *obj, int *dims, int l)
{
  PyObject *list;
  int i, c, *strides;

  strides = (int *)PyMem_Malloc(l*sizeof(int));
  if (!strides)
    PyErr_NoMemory();

  c = 1;
  for (i=0; i<l; i++) {
    strides[i] = c;
    c *= dims[i];
  }

  list = ltranspose(obj, dims, strides, 0, 0, l);
  PyMem_Free(strides);

  return list;
}

/* Convert a Python sequence to a Numeric array */
#ifdef WITH_NUMERIC
static PyObject *
to_PyNumericArray(PyObject *seq, SEXP dim)
{
  PyObject *array, *ret, *dims, *it;
  int l, i, j;

  array = PyArray_ContiguousFromObject(seq, PyArray_DOUBLE, 0,0);
  if (!array)
    return NULL;

  l = GET_LENGTH(dim);
  dims = PyList_New(l);
  for (i=0; i<l; i++) {
    j = INTEGER(dim)[l-i-1];
    if (j == 0) {
      Py_DECREF(array);
      Py_DECREF(dims);
      Py_INCREF(Py_None);
      return Py_None;
    }
    if (!(it = PyInt_FromLong(j)))
      return NULL;
    if (PyList_SetItem(dims, i, it) < 0)
      return NULL;
  }

  ret = PyArray_Reshape((PyArrayObject *)array, dims);
  Py_DECREF(array);
  Py_DECREF(dims);
  if (!ret)
    return NULL;

  array = PyObject_CallFunction(Py_transpose, "O", ret);
  Py_XDECREF(ret);
  return array;
}
#endif

/* Convert an R object to a 'basic' Python object (mode 2) */
/* NOTE: R vectors of length 1 will yield a python scalar */
int
to_Pyobj_basic(SEXP robj, PyObject **obj)
{
  int status;
  PyObject *tmp;

  status = to_Pyobj_vector(robj, &tmp, BASIC_CONVERSION);

  if(status==1 && PyList_Check(tmp) && PyList_Size(tmp) == 1)
    {
      *obj = PyList_GetItem(tmp, 0);
      Py_XINCREF(*obj);
      Py_DECREF(tmp);
    }
  else
    *obj = tmp;

  return status;
}


/* Convert an R object to a 'vector' Python object (mode 1) */
/* NOTE: R vectors of length 1 will yield a python list of length 1*/
int
to_Pyobj_vector(SEXP robj, PyObject **obj, int mode)
{
  PyObject *it, *tmp;
  SEXP names, dim;
  int len, *integers, i, type;
  const char *strings, *thislevel;
  double *reals;
  Rcomplex *complexes;
#ifdef WITH_NUMERIC
  PyObject *array;
#endif

  if (!robj)
    return -1;                  /* error */

  if (robj == R_NilValue) {
    Py_INCREF(Py_None);
    *obj = Py_None;
    return 1;                   /* succeed */
  }

  len = GET_LENGTH(robj);
  tmp = PyList_New(len);
  type = TYPEOF(robj);

  for (i=0; i<len; i++) {
    switch (type)
      {
      case LGLSXP:
         integers = INTEGER(robj);
         if(integers[i]==NA_INTEGER) /* watch out for NA's */
           {
             if (!(it = PyInt_FromLong(integers[i])))
             return -1;
             //it = Py_None;
           }
         else if (!(it = PyBool_FromLong(integers[i])))
           return -1;
         break;
      case INTSXP:
        integers = INTEGER(robj);
        if(isFactor(robj)) {
          /* Watch for NA's! */
          if(integers[i]==NA_INTEGER)
            it = PyString_FromString(CHAR(NA_STRING));
          else
            {
              thislevel = CHAR(STRING_ELT(GET_LEVELS(robj), integers[i]-1));
              if (!(it = PyString_FromString(thislevel)))
                return -1;
            }
        }
        else {
          if (!(it = PyInt_FromLong(integers[i])))
            return -1;
        }
        break;
      case REALSXP:
        reals = REAL(robj);
        if (!(it = PyFloat_FromDouble(reals[i])))
          return -1;
        break;
      case CPLXSXP:
        complexes = COMPLEX(robj);
        if (!(it = PyComplex_FromDoubles(complexes[i].r,
                                         complexes[i].i)))
          return -1;
        break;
      case STRSXP:
        if(STRING_ELT(robj, i)==R_NaString)
          it = PyString_FromString(CHAR(NA_STRING));
        else
          {
            strings = CHAR(STRING_ELT(robj, i));
            if (!(it = PyString_FromString(strings)))
              return -1;
          }
        break;
      case LISTSXP:
        if (!(it = to_Pyobj_with_mode(elt(robj, i), mode)))
          return -1;
        break;
      case VECSXP:
        if (!(it = to_Pyobj_with_mode(VECTOR_ELT(robj, i), mode)))
          return -1;
        break;
      default:
        Py_DECREF(tmp);
        return 0;                 /* failed */
    }

    if (PyList_SetItem(tmp, i, it) < 0)
      return -1;
  }

  dim = GET_DIM(robj);
  if (dim != R_NilValue) {
#ifdef WITH_NUMERIC
    if(use_numeric)
      {
        array = to_PyNumericArray(tmp, dim);
        if (array) {                /* If the conversion to Numeric succeed.. */
          *obj = array;             /* we are done */
          Py_DECREF(tmp);
          return 1;
        }
        PyErr_Clear();
      }
#endif
    len = GET_LENGTH(dim);
    *obj = to_PyArray(tmp, INTEGER(dim), len);
    Py_DECREF(tmp);
    return 1;
  }

  names = GET_NAMES(robj);
  if (names == R_NilValue)
    *obj = tmp;
  else {
    *obj = to_PyDict(tmp, names);
    Py_DECREF(tmp);
  }
  return 1;
}

/* Search a conversion procedure from the class attribute */
PyObject *
from_class_table(SEXP robj)
{
  SEXP rclass;
  PyObject *lkey, *key, *fun;
  int i;

  PROTECT(rclass = GET_CLASS(robj));

  fun = NULL;
  if (rclass != R_NilValue) {

    lkey = to_Pyobj_with_mode(rclass, BASIC_CONVERSION);
    key = PyList_AsTuple(lkey);
    if (key) {
      Py_DECREF(lkey);
    } else {
      PyErr_Clear();
      key = lkey;
    }
    fun = PyDict_GetItem(class_table, key);
    Py_DECREF(key);

    if (!fun) {
      PyErr_Clear();
      for (i=0; i<GET_LENGTH(rclass); i++)
        if ((fun = PyDict_GetItemString(class_table,
                                        CHAR(STRING_ELT(rclass, i)))))
          break;
    }
    else
      Py_INCREF(fun);
  }
  UNPROTECT(1);
  return fun;
}

/* Search a conversion procedure from the proc table */
int
from_proc_table(SEXP robj, PyObject **fun)
{
  PyObject *procs, *proc, *funs, *res, *obj;
  int i, l, k, error;

  proc = NULL;
  procs = PyDict_Keys(proc_table);
  funs = PyDict_Values(proc_table);
  l = PyMapping_Size(proc_table);

  obj = (PyObject *)Robj_new(robj, TOP_MODE);

  error = 0;
  for (i=0; i<l; i++) {
    proc = PyList_GetItem(procs, i);
    Py_XINCREF(proc);
    res = PyObject_CallFunction(proc, "O", obj);
    if (!res) {
      error = -1;
      break;
    }
    k = PyObject_IsTrue(res);
    Py_DECREF(res);
    if (k) {
      *fun = PyList_GetItem(funs, i);
      Py_XINCREF(*fun);
      break;
    }
  }

  Py_DECREF(obj);
  Py_XDECREF(proc);
  Py_XDECREF(procs);
  Py_XDECREF(funs);
  return error;
}

int
to_Pyobj_proc(SEXP robj, PyObject **obj)
{
  PyObject *fun=NULL, *tmp;
  int i;

  i = from_proc_table(robj, &fun);
  if (i < 0)
    return -1;                  /* an error occurred */

  if (!fun)
    return 0;                   /* conversion failed */

  tmp = (PyObject *)Robj_new(robj, TOP_MODE);
  *obj = PyObject_CallFunction(fun, "O", tmp);
  Py_DECREF(fun);
  Py_DECREF(tmp);
  return 1;                     /* conversion succeed */
}

/* Convert a Robj to a Python object via the class table (mode 3) */
/* See the docs for conversion rules */
int
to_Pyobj_class(SEXP robj, PyObject **obj)
{
  PyObject *fun, *tmp;

  fun = from_class_table(robj);

  if (!fun)
    return 0;                   /* conversion failed */

  tmp = (PyObject *)Robj_new(robj, TOP_MODE);
  *obj = PyObject_CallFunction(fun, "O", tmp);
  Py_DECREF(fun);
  Py_DECREF(tmp);
  return 1;                     /* conversion succeed */
}

PyObject *
to_Pyobj_with_mode(SEXP robj, int mode)
{
  PyObject *obj;
  int i;

  switch (mode)
    {
    case PROC_CONVERSION:
      i = to_Pyobj_proc(robj, &obj);
      if (i<0) return NULL;
      if (i==1) break;
    case CLASS_CONVERSION:
      i = to_Pyobj_class(robj, &obj);
      if (i<0) return NULL;
      if (i==1) break;
    case BASIC_CONVERSION:
      i = to_Pyobj_basic(robj, &obj);
      if (i<0) return NULL;
      if (i==1) break;
    case VECTOR_CONVERSION:
      i = to_Pyobj_vector(robj, &obj, mode=VECTOR_CONVERSION);
      if (i<0) return NULL;
      if (i==1) break;
    default:
      obj = (PyObject *)Robj_new(robj, TOP_MODE);
  }

  return obj;
}

/* Convert a tuple to arguments for a R function */
int
make_args(int largs, PyObject *args, SEXP *e)
{
  SEXP r;
  int i;

  for (i=0; i<largs; i++) {
    r = to_Robj(PyTuple_GetItem(args, i));
    if (!r || PyErr_Occurred()	)
      return 0;
    SETCAR(*e, r);
    *e = CDR(*e);
  }
  return 1;
}

/* Implements the conversion rules for names. See the 'USING' file. We
   don't care about '<-' because it doesn't appear in keywords. */
const char *
dotter(char *s)
{
  char *r, *res;
  int l;

  if (!s)
    return NULL;                /* assume prev PyString_AsString has failed */
  l = strlen(s);
  r = (char *)PyMem_Malloc(l+1);
  if (!r) {
    PyErr_NoMemory();
    return NULL;
  }
  res = strcpy(r, s);

  if ((l > 1) && (res[l-1] == '_') && (res[l-2] != '_'))
    res[l-1]=0;

  while ((r=strchr(r, '_')))
    *r = '.';

  return res;
}

/* Convert a dict to keywords arguments for a R function */
int
make_kwds(int lkwds, PyObject *kwds, SEXP *e)
{
  SEXP r;
  const char *s;
  int i;
  PyObject *citems=NULL, *it;
  PyObject *kwname;

  if (kwds) {
    citems = PyMapping_Items(kwds);
  }

  for (i=0; i<lkwds; i++) {
    it = PySequence_GetItem(citems, i);
    if (!it)
      goto fail;
    r = to_Robj(PyTuple_GetItem(it, 1));
    Py_DECREF(it);
    if (!r || PyErr_Occurred())
      goto fail;

    SETCAR(*e, r);
    kwname = PyTuple_GetItem(it, 0);
    if (!kwname)
      goto fail;
    if (!PyString_Check(kwname)) {
      PyErr_SetString(PyExc_TypeError, "keywords must be strings");
      goto fail;
    }
    s = dotter(PyString_AsString(kwname));
    if (!s)
      goto fail;

    SET_TAG(*e, Rf_install(s));
    PyMem_Free( (void*) s);
    *e = CDR(*e);
  }
  Py_XDECREF(citems);
  return 1;

 fail:
  Py_XDECREF(citems);
  return 0;
}


/* This is the method to call when invoking an 'Robj' */
static PyObject *
Robj_call(PyObject *self, PyObject *args, PyObject *kwds)
{
  SEXP exp, e, res;
  int largs, lkwds, conv;
  PyObject *obj;

  largs = lkwds = 0;
  if (args)
    largs = PyObject_Length(args);
  if (kwds)
    lkwds = PyObject_Length(kwds);
  if ((largs<0) || (lkwds<0))
    return NULL;

  /* A SEXP with the function to call and the arguments and keywords. */
  PROTECT(exp = allocVector(LANGSXP, largs+lkwds+1));
  e = exp;
  SETCAR(e, ((RobjObject *)self)->R_obj);
  e = CDR(e);

  if (!make_args(largs, args, &e)) {
    UNPROTECT(1);
    return NULL;
  }
  if (!make_kwds(lkwds, kwds, &e)) {
    UNPROTECT(1);
    return NULL;
  }

  PROTECT(res = do_eval_expr(exp));
  if (!res) {
    UNPROTECT(2);
    return NULL;
  }

  if (default_mode < 0)
    conv = ((RobjObject *)self)->conversion;
  else
    conv = default_mode;

  obj = to_Pyobj_with_mode(res, conv);
  UNPROTECT(2);

  PrintWarnings(); /* show any warning messages */

  return obj;
}

/* Convert a sequence of (name, value) pairs to arguments to an R
   function call */
int
make_argl(int largl, PyObject *argl, SEXP *e)
{
  SEXP rvalue;
  const char *name;
  int i;
  PyObject *it, *nobj, *value;

  if( !PySequence_Check(argl) ) goto fail_arg;

  for (i=0; i<largl; i++) {
    it = PySequence_GetItem(argl, i);
    if(!it) goto fail_arg;
    if( PySequence_Size(it) != 2 )
      {
        Py_DECREF(it);
        goto fail_arg;
      }
    nobj = PySequence_GetItem(it, 0);

    /* Name can be a string, None, or NULL, error otherwise. */
    if (PyString_Check(nobj))
      {
        name = dotter(PyString_AsString(nobj));
        Py_DECREF(nobj);
      }
    else if (nobj == Py_None)
      {
        name = NULL;
        Py_DECREF(nobj);
      }
    else if(nobj == NULL)
      {
        name = NULL;
      }
    else
      {
        Py_DECREF(nobj);
        goto fail_arg;
      }

    /* Value can be anything. */
    value = PySequence_GetItem(it, 1);
    if (!value || PyErr_Occurred() )
      {
        PyMem_Free( (void*) name);
        goto fail;
      }

    rvalue =  to_Robj(value);
    Py_DECREF(value);
    Py_DECREF(it);

    if(PyErr_Occurred())
      goto fail;

    /* Add parameter value to call */
    SETCAR(*e, rvalue);

    /* Add name (if present) */
    if (name && strlen(name)>0)
      {
        SET_TAG(*e, Rf_install(name));
        PyMem_Free((void*) name);
      }

    /* Move index to new end of call */
    *e = CDR(*e);
  }
  return 1;

 fail_arg:
   PyErr_SetString(PyExc_ValueError,
                "Argument must be a sequence of (\"name\", value) pairs.\n");
 fail:
   return 0;
}

/* Methods for the 'Robj' type */

/* Explicitly call an R object with a list containing (name, value) *
 * argument pairs.  'name' can be None or '' to provide unnamed
 * arguments.  This function is necessary when the *order* of named
 * arguments needs to be preserved.
 */

static PyObject *
Robj_lcall(PyObject *self, PyObject *args)
{
  SEXP exp, e, res;
  int largs, largl, conv;
  PyObject *obj, *argl;

  /* Check arguments, there should be *exactly one* unnamed sequence. */
  largs = 0;
  if (args)
    largs = PyObject_Length(args);
  if (largs<0)
    return NULL;

  if(largs != 1 || !PySequence_Check(args) )
    {
      PyErr_SetString(PyExc_ValueError,
                "Argument must be a sequence of (\"name\", value) pairs.\n");
      return NULL;
    }

  // extract our one argument
  argl = PySequence_GetItem(args, 0);
  Py_DECREF(args);

  largl = 0;
  if (argl)
    largl = PyObject_Length(argl);
  if (largl<0)
    return NULL;

  // A SEXP with the function to call and the arguments
  PROTECT(exp = allocVector(LANGSXP, largl+1));
  e = exp;
  SETCAR(e, ((RobjObject *)self)->R_obj);
  e = CDR(e);

  // Add the arguments to the SEXP
  if (!make_argl(largl, argl, &e)) {
    UNPROTECT(1);
    return NULL;
  }

  // Evaluate
  PROTECT(res = do_eval_expr(exp));
  if (!res) {
    UNPROTECT(2);
    return NULL;
  }

  // Convert
  if (default_mode < 0)
    conv = ((RobjObject *)self)->conversion;
  else
    conv = default_mode;

  obj = to_Pyobj_with_mode(res, conv);
  UNPROTECT(2);

  // Return
  return obj;
}


/* Without args return the value of the conversion flag. With an
   argument set the conversion flag to the truth value of the argument. */
static PyObject *
Robj_autoconvert(PyObject *self, PyObject *args, PyObject *kwds)
{
  PyObject *obj;
  int conversion=-2;
  char *kwlist[] = {"val", 0};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|i:autoconvert", kwlist,
                                  &conversion))
    return NULL;

  if (conversion > TOP_MODE) {
    PyErr_SetString(PyExc_ValueError, "wrong mode");
    return NULL;
  }

  if (conversion == -2) {
    obj = PyInt_FromLong((long)((RobjObject *)self)->conversion);
  } else {
    ((RobjObject *)self)->conversion = conversion;
    obj = Py_None;
    Py_XINCREF(obj);
  }

  return obj;
}

static PyObject *
Robj_as_py(PyObject *self, PyObject *args, PyObject *kwds)
{
  PyObject *obj;
  char *kwlist[] = {"mode", 0};
  int conv=default_mode;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|i:as_py", kwlist,
                                   &conv))
    return NULL;

  if (conv <= -2 || conv > TOP_MODE) {
    PyErr_SetString(PyExc_ValueError, "wrong mode");
    return NULL;
  }

  if (conv < 0)
    conv = TOP_MODE;

  obj = to_Pyobj_with_mode(((RobjObject *)self)->R_obj, conv);
  return obj;
}

static PyMethodDef Robj_methods[] = {
  {"autoconvert", (PyCFunction)Robj_autoconvert, METH_VARARGS|METH_KEYWORDS},
  {"local_mode", (PyCFunction)Robj_autoconvert, METH_VARARGS|METH_KEYWORDS},
  {"as_py", (PyCFunction)Robj_as_py, METH_VARARGS|METH_KEYWORDS},
  {"lcall", (PyCFunction)Robj_lcall, METH_VARARGS},
  {NULL, NULL}          /* sentinel */
};

/* Sequence protocol implementation */

/* len(a) */
static int
Robj_len(PyObject *a)
{
  SEXP e, robj;

  PROTECT(e = allocVector(LANGSXP, 2));
  SETCAR(e, length);
  SETCAR(CDR(e), ((RobjObject *)a)->R_obj);

  if (!(robj = do_eval_expr(e))) {
    UNPROTECT(1);
    return -1;
  }

  UNPROTECT(1);
  return INTEGER_DATA(robj)[0];
}

/* a[i] = v */
static int
Robj_ass_item(PyObject *a, int i, PyObject *v)
{
  SEXP e, ri, robj;

  PROTECT(e = allocVector(LANGSXP, 4));
  ri = NEW_INTEGER(1);
  INTEGER_DATA(ri)[0] = i+1;
  SETCAR(e, set_item);
  SETCAR(CDR(e), ((RobjObject *)a)->R_obj);
  SETCAR(CDR(CDR(e)), ri);
  SETCAR(CDR(CDR(CDR(e))), to_Robj(v));

  if(PyErr_Occurred())
    return -1;

  if (!(robj = do_eval_expr(e))) {
    UNPROTECT(1);
    return -1;
  }

  ((RobjObject *)a)->R_obj = robj;
  UNPROTECT(1);
  return 0;
}

/* a[i] */
static PyObject *
Robj_item(PyObject *a, int i)
{
  SEXP ri, robj, e;
  PyObject *obj;
  int len, c;

  if ((len = Robj_len(a)) < 0)
    return NULL;
  if (i >= len || i < 0) {
    PyErr_SetString(PyExc_IndexError, "R object index out of range");
    return NULL;
  }

  PROTECT(ri = NEW_INTEGER(1));
  INTEGER_DATA(ri)[0] = i+1;
  PROTECT(e = allocVector(LANGSXP, 3));
  SETCAR(e, get_item);
  SETCAR(CDR(e), ((RobjObject *)a)->R_obj);
  SETCAR(CDR(CDR(e)), ri);

  if (!(robj = do_eval_expr(e))) {
    UNPROTECT(2);
    return NULL;
  }

  UNPROTECT(2);

  /* If there is a default mode, use it; otherwise, use the top mode. */
  if (default_mode < 0)
    c = TOP_MODE;
  else
    c = default_mode;
  obj = to_Pyobj_with_mode(robj, c);
  return obj;
}

/* We should implement sq_slice, sq_contains ... */
static PySequenceMethods Robj_as_sequence = {
  (inquiry)Robj_len,              /* sq_length */
  0,                            /* sq_concat */
  0,                            /* sq_repeat */
  (intargfunc)Robj_item,           /* sq_item */
  0,                            /* sq_slice */
  (intobjargproc)Robj_ass_item,     /* sq_ass_item */
  0,                            /* sq_ass_slice */
  0,                            /* sq_contains */
};


/* The 'Robj' table. When compiled under Python 2.2, the type 'Robj'
   is subclassable. */

#ifdef PRE_2_2
static PyObject *
Robj_getattr(RobjObject *self, char *name)
{
  return Py_FindMethod(Robj_methods, (PyObject *)self, name);
}
#endif

PyTypeObject Robj_Type = {
  /* The ob_type field must be initialized in the module init function
   * to be portable to Windows without using C++. */
#if defined(PRE_2_2) || defined(_WIN32)    // Matjaz
  PyObject_HEAD_INIT(NULL)
#else
  PyObject_HEAD_INIT(&PyType_Type)
#endif
  0,                    /*ob_size*/
  "Robj",               /*tp_name*/
  sizeof(RobjObject),   /*tp_basicsize*/
  0,                    /*tp_itemsize*/
  /* methods */
  (destructor)Robj_dealloc, /*tp_dealloc*/
  0,                    /*tp_print*/
#ifdef PRE_2_2
  (getattrfunc)Robj_getattr,
#else
  0,
#endif
  0,
  0,                    /*tp_compare*/
  0,                    /*tp_repr*/
  0,                    /*tp_as_number*/
  &Robj_as_sequence,    /*tp_as_sequence*/
  0,                    /*tp_as_mapping*/
  0,                    /*tp_hash*/
  (ternaryfunc)Robj_call,  /*tp_call*/
  0,                    /*tp_str*/
#if defined(PRE_2_2) || defined(_WIN32)
  0,
#else
  PyObject_GenericGetAttr, /*tp_getattro*/
#endif
  0,                      /*tp_setattro*/
  0,                      /*tp_as_buffer*/
#ifdef PRE_2_2
  0,
#else
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,     /*tp_flags*/
#endif
  0,                      /*tp_doc*/
  0,                      /*tp_traverse*/
#ifndef PRE_2_2
  0,                      /*tp_clear*/
  0,                      /*tp_richcompare*/
  0,                      /*tp_weaklistoffset*/
  0,                      /*tp_iter*/
  0,                      /*tp_iternext*/
  Robj_methods,           /*tp_methods*/
  0,                      /*tp_members*/
  0,                      /*tp_getset*/
  0,                      /*tp_base*/
  0,                      /*tp_dict*/
  0,                      /*tp_descr_get*/
  0,                      /*tp_descr_set*/
  0,                      /*tp_dictoffset*/
  0,                      /*tp_init*/
#ifdef _WIN32
  0,                      /*tp_alloc*/
#else
  PyType_GenericAlloc,    /*tp_alloc*/
#endif
  Robj_tpnew,             /*tp_new*/
  0,                      /*tp_free*/
  0,                      /*tp_is_gc*/
#endif
};


/* Module functions */

/* Obtain an R object via its name. 'autoconvert' is the keyword to
   set the autoconversion flag. */
static PyObject *
get_fun(PyObject *self, PyObject *args, PyObject *kwds)
{
  char *obj_str;
  int conversion=TOP_MODE;
  SEXP robj;

  static char *kwlist[] = {"name", "autoconvert", 0};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|i:get", kwlist,
                                   &obj_str, &conversion))
    return NULL;

  robj = get_fun_from_name(obj_str);
  if (!robj)
    return NULL;

  return (PyObject *)Robj_new(robj, conversion);
}

static PyObject *
set_mode(PyObject *self, PyObject *args)
{
  int i=-1;

  if (!PyArg_ParseTuple(args, "i:set_mode", &i))
    return NULL;

  if (i<-1 || i>TOP_MODE) {
    PyErr_SetString(PyExc_ValueError, "wrong mode");
    return NULL;
  }

  default_mode = i;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
get_mode(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, ":get_mode"))
    return NULL;

  return PyInt_FromLong(default_mode);
}

static PyObject *
r_events(PyObject *self, PyObject *args, PyObject *kwds)
#ifdef _WIN32
{
  return NULL;
}
#else
{
  fd_set *what;
  int usec=10000;

  static char *kwlist[] = {"usec", 0};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|i:r_events",
                                   kwlist, &usec))
    return NULL;

  if (R_interact) {
    Py_BEGIN_ALLOW_THREADS
    what = R_checkActivity(usec, 0);
    R_runHandlers(R_InputHandlers, what);
    Py_END_ALLOW_THREADS
  }

  Py_INCREF(Py_None);
  return Py_None;
}
#endif

void
stop_events(void)
{
  PyObject *o;

  if (!rpy_dict)
    return;

  if (!r_lock)
    r_lock = PyDict_GetItemString(rpy_dict, "_r_lock");

  o = PyObject_CallMethod(r_lock, "acquire", NULL);
  Py_XDECREF(o);
}

void
start_events(void)
{
  PyObject *o;

  if (!rpy_dict)
    return;

  if (!r_lock)
    r_lock = PyDict_GetItemString(rpy_dict, "_r_lock");

  o = PyObject_CallMethod(r_lock, "release", NULL);
  Py_XDECREF(o);
}


/*
 * Based on code from Rstd_CleanUp();
 * from src/unix/sys-std.c
 */

void r_finalize(void)
{
    unsigned char buf[1024];
    char * tmpdir;

    R_dot_Last();
    R_RunExitFinalizers();
    CleanEd();
    KillAllDevices();
    if((tmpdir = getenv("R_SESSION_TMPDIR"))) {

#ifdef _WIN32
      snprintf((char *)buf, 1024, "rmdir /S /Q %s", tmpdir);
#else
      snprintf((char *)buf, 1024, "rm -rf %s", tmpdir);
#endif

      R_system((char *)buf);
    }

    PrintWarnings();	/* from device close and .Last */
    R_gc();  /* Remove any remaining R objects from memory */
}


static PyObject *
r_cleanup(void)
{
  r_cleanup();
  Py_INCREF(Py_None);
  return Py_None;
}

#ifdef WITH_NUMERIC
static void
init_numeric(void)
{
  PyObject *multiarray, *dict;

  if(use_numeric)
    {
      import_array();
      multiarray = PyImport_ImportModule(PY_ARRAY_MODULE_NAME);
      if (multiarray) {
        dict = PyModule_GetDict(multiarray);
        if (dict)
          Py_transpose = PyDict_GetItemString(dict, "transpose");
      }
    }
}
#endif

static PyObject *
r_init(PyObject *self, PyObject *args)
{
  static int first=1;
  int i;

  if (!PyArg_ParseTuple(args, "i:r_init", &i))
    return NULL;
  use_numeric = i;

#ifdef WITH_NUMERIC
  if(use_numeric)
    init_numeric();
#endif

  if(first==1)
    {
      first=0;
      Py_INCREF(Py_None);
      return Py_None;
    }
  else
    {
      PyErr_SetString(PyExc_RuntimeError, "Only one R object may be instantiated per session");
      return NULL;
    }
}

/* List of functions defined in the module */
static PyMethodDef rpy_methods[] = {
  {"get_fun",       (PyCFunction)get_fun,           METH_VARARGS | METH_KEYWORDS},
  {"set_mode",      (PyCFunction)set_mode,      METH_VARARGS},
  {"get_mode",      (PyCFunction)get_mode,      METH_VARARGS},
  {"set_output",    (PyCFunction)set_output,    METH_VARARGS},
  {"set_input",     (PyCFunction)set_input,     METH_VARARGS},
  {"set_showfiles", (PyCFunction)set_showfiles, METH_VARARGS},
  {"get_output",    (PyCFunction)get_output,    METH_VARARGS},
  {"get_input",     (PyCFunction)get_input,     METH_VARARGS},
  {"get_showfiles", (PyCFunction)get_showfiles, METH_VARARGS},
  {"r_events",      (PyCFunction)r_events,      METH_VARARGS | METH_KEYWORDS},
  {"r_cleanup",     (PyCFunction)r_cleanup,     METH_NOARGS},
  {"r_init",        (PyCFunction)r_init,        METH_VARARGS},
  {NULL, NULL}          /* sentinel */
};

#ifdef _WIN32
static void char_message( char *s )
{
  if (!s) return;
  R_WriteConsole(s, strlen(s));
}

static int char_yesnocancel( char *s )
{
  return 1;
}

static void
RPyBusy( int which )
{
  /* set a busy cursor ... in which = 1, unset if which = 0 */
}

static void
RPyDoNothing( void )
{
}

/* initialise embedded R; based on rproxy_impl.c from the R distribution */
static void
init_embedded_win32( void ) {
  structRstart rp;
  Rstart Rp = &rp;
  char Rversion[25];
  int index;


  snprintf( Rversion, 25, "%s.%s", R_MAJOR, R_MINOR );
  if( strcmp( getDLLVersion(), Rversion ) != 0 ) {
    PyErr_SetString( PyExc_ImportError, "R.DLL version does not match" );
    return;
  }

  R_DefParams(Rp);

  /* set R_HOME */
  Rp->rhome = RHOME;

  index = strlen(RUSER) - 1;

  if (RUSER[index] == '/' || RUSER[index] == '\\')
    RUSER[index] = '\0';

  Rp->home = RUSER;
  Rp->CharacterMode = LinkDLL;

  Rp->ReadConsole = (blah1) RPy_ReadConsole;    // Matjaz
  Rp->WriteConsole = (blah2) RPy_WriteConsole;  // Matjaz

  Rp->CallBack = (blah3) RPyDoNothing;
#if R_VERSION < 0x20100
  Rp->message = char_message;
  Rp->yesnocancel = char_yesnocancel;
  Rp->busy = RPyBusy;
#else
  Rp->ShowMessage = char_message;
  Rp->YesNoCancel = char_yesnocancel;
  Rp->Busy = RPyBusy;
#endif

  Rp->R_Quiet = TRUE;

  /* run as "interactive", so server won't be killed after an error */
  Rp->R_Slave = Rp->R_Verbose = 0;
  Rp->R_Interactive = TRUE;
  Rp->RestoreAction = SA_NORESTORE; /* no restore */
  Rp->SaveAction    = SA_NOSAVE;  /* no save */

#if R_VERSION < 0x20000   // pre-R-2.0.0

  Rp->CommandLineArgs = NULL;
  Rp->NumCommandLineArgs = 0;
#else
  R_set_command_line_arguments(0, NULL);
#endif
  R_SetParams(Rp); /* so R_ShowMessage is set */
  R_SizeFromEnv(Rp);

  R_SetParams(Rp);

  setup_term_ui();
  setup_Rmainloop();
}
#endif

/* Initialization function for the module */
DL_EXPORT(void)
INIT_RPY(void)
{
  PyObject *m, *d;
  PyOS_sighandler_t old_int;
#ifndef _WIN32
  char *defaultargv[] = {"rpy", "-q", "--vanilla"};
  PyOS_sighandler_t old_usr1, old_usr2;
#endif
  SEXP interact;

  /* Get path and version information from environment */
  strncpy(RHOME,    getenv("RPY_RHOME"),    BUFSIZ);
  strncpy(RVERSION, getenv("RPY_RVERSION"), BUFSIZ);
  strncpy(RVER,     getenv("RPY_RVER"),     BUFSIZ);
  strncpy(RUSER,    getenv("RPY_RUSER"),    BUFSIZ);

  if( !strlen(RHOME) || !strlen(RVERSION) || !strlen(RVER) || !strlen(RUSER))
    {
      PyErr_Format(RPy_Exception,
                   "Unable to load R path or version information");
      return;
    }

  Robj_Type.ob_type = &PyType_Type;
#if defined( _WIN32 ) && ! defined( PRE_2_2 )
  Robj_Type.tp_getattro = PyObject_GenericGetAttr;
  Robj_Type.tp_alloc = PyType_GenericAlloc;
#endif

  m = Py_InitModule(xstr(RPY_SHNAME), rpy_methods);
  d = PyModule_GetDict(m);

  /* Save this interpreter */
  PyEval_InitThreads();
  my_interp = PyThreadState_Get()->interp;

  /* Save the Python signal handlers. If R inserts its handlers, we
     cannot return to the Python interpreter. */
  old_int = PyOS_getsig(SIGINT);
  python_sigint = old_int;
#ifndef _WIN32
  old_usr1 = PyOS_getsig(SIGUSR1);
  old_usr2 = PyOS_getsig(SIGUSR2);
#endif

#ifdef _WIN32
  init_embedded_win32();
#else
  Rf_initEmbeddedR( sizeof(defaultargv) / sizeof(defaultargv[0]),
                    defaultargv);
#endif


#ifndef CSTACK_DEFNS
  /* Disable C stack checking, which is incompatible with use as a
     shared library. */
  R_CStackLimit = (uintptr_t)-1;
#endif

  /* Restore Python handlers */
  PyOS_setsig(SIGINT, old_int);
#ifndef _WIN32
  PyOS_setsig(SIGUSR1, old_usr1);
  PyOS_setsig(SIGUSR2, old_usr2);
#endif

  /* Several new exceptions: */
  RPy_Exception               = PyErr_NewException("rpy.RPy_Exception",               NULL,          NULL);
  RPy_TypeConversionException = PyErr_NewException("rpy.RPy_TypeConversionException", RPy_Exception, NULL);
  RPy_RException              = PyErr_NewException("rpy.RPy_RException",              RPy_Exception, NULL);

  if (!RPy_Exception || !RPy_TypeConversionException || !RPy_RException )
    {
      PyErr_Format(RPy_Exception, "Unable create RPy exceptions");
      return;
    }

  PyDict_SetItemString(d, "RPy_Exception",               RPy_Exception);
  PyDict_SetItemString(d, "RPy_TypeConversionException", RPy_TypeConversionException);
  PyDict_SetItemString(d, "RPy_RException",              RPy_RException);

  // The conversion table
  class_table = PyDict_New();
  proc_table = PyDict_New();
  PyDict_SetItemString(d, "__class_table__", class_table);
  PyDict_SetItemString(d, "__proc_table__", proc_table);

  // The globals R objects for the sequence protocol
  get_item = get_fun_from_name("[");
  set_item = get_fun_from_name("[<-");
  length = get_fun_from_name("length");

  // Function to transpose arrays
  aperm = get_fun_from_name("aperm");

  // Initialize the list of protected objects
  R_References = R_NilValue;
  SET_SYMVALUE(install("R.References"), R_References);

  // Initialize the default mode
  default_mode = -1;

  // Check whether R is interactive or no
  interact = do_eval_fun("interactive");
  R_interact = INTEGER(interact)[0];

  // I/O routines
  init_io_routines();

  rpy = PyImport_ImportModule("rpy");
  rpy_dict = PyModule_GetDict(rpy);
  //  r_lock = PyDict_GetItemString(rpy_dict, "_r_lock");
  //  PyObject_Print(r_lock, stderr, Py_PRINT_RAW);
  r_lock = NULL;

  if( Py_AtExit( r_finalize ) )
    {
      fprintf(stderr, "Warning: Unable to set R finalizer.");
      fflush(stderr);
    }


}


