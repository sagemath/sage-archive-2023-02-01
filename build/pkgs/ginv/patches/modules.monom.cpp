/***************************************************************************
 *   Copyright (C) 2004 by Blinkov Yu.A.                                   *
 *   BlinkovUA@info.sgu.ru                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "./systemtype.h"
#include "./monom.h"

#include <string.h>
#include <sstream>

void monomInterface_dealloc(MonomInterface *self) {
//   delete self->monomInterface;
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* monomInterface_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  MonomInterface *self = (MonomInterface *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->monomInterface = NULL;
  }

  return (PyObject *)self;
}

int monomInterface_init(MonomInterface *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"order", "systemType", "independ", "depend", "varSep", NULL};

  const char *py_order = NULL;
  PyObject *py_systemType = NULL;
  PyObject *py_independ = NULL;
  PyObject *py_depend = NULL;
  int varSep = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "sO!O!|O!i:MonomInterface", kwlist,
                                    &py_order,
                                    &systemType_type, &py_systemType,
                                    &PyList_Type, &py_independ,
                                    &PyList_Type, &py_depend,
                                    &varSep))
      return -1;


  ISystemType *systemType = ((SystemType *)py_systemType)->systemType;

  IVariables *independ = new IVariables;
  int n = PyList_Size(py_independ);
  for(int i=0; i < n; i++) {
    independ->add(PyString_AsString(PyList_GetItem(py_independ, i)));
  }

  IVariables *depend = NULL;
  if (py_depend) {
    depend = new IVariables;
    int n = PyList_Size(py_depend);
    for(int i=0; i < n; i++) {
      depend->add(PyString_AsString(PyList_GetItem(py_depend, i)));
    }
  }

  IMonomInterface::Order orderType = IMonomInterface::Indefinitely;
  if (strcmp(py_order, "DegRevLex") == 0)
    orderType = IMonomInterface::DegRevLex;
  else if (strcmp(py_order, "DegRevLexByte") == 0)
    orderType = IMonomInterface::DegRevLexByte;
  else if (strcmp(py_order, "Lex") == 0)
    orderType = IMonomInterface::Lex;
  else if (strcmp(py_order, "Elim") == 0)
    orderType = IMonomInterface::Elim;
  else if (strcmp(py_order, "TopDegRevLexByte") == 0)
    orderType = IMonomInterface::TopDegRevLexByte;
  else if (strcmp(py_order, "TopDegRevLex") == 0)
    orderType = IMonomInterface::TopDegRevLex;
  else if (strcmp(py_order, "TopElim") == 0)
    orderType = IMonomInterface::TopElim;
  else if (strcmp(py_order, "TopLex") == 0)
    orderType = IMonomInterface::TopLex;
  else if (strcmp(py_order, "PotDegRevLex") == 0)
    orderType = IMonomInterface::PotDegRevLex;
  else if (strcmp(py_order, "PotLex") == 0)
    orderType = IMonomInterface::PotLex;
  else if (strcmp(py_order, "PosElim") == 0)
    orderType = IMonomInterface::PosElim;
  else {
    PyErr_SetString(PyExc_TypeError, "unknown value 'order' in 'MonomInterface'");
    return -1;
  }

  if (!(orderType & (IMonomInterface::Pot | IMonomInterface::Top)) &&
        (depend != NULL || systemType->module() > 0)) {
    PyErr_SetString(PyExc_TypeError, "error value 'depend' in 'MonomInterface' or 'module' in 'SystemType'");
    return -1;
  }

  if ((orderType & (IMonomInterface::Pot | IMonomInterface::Top)) &&
       (depend == NULL && systemType->module() == 0)) {
    PyErr_SetString(PyExc_TypeError, "unknown value 'depend' in 'MonomInterface' or 'module=0' in 'SystemType'");
    return -1;
  }

  if (orderType & IMonomInterface::Elim) {
    if (varSep <= 0) {
      PyErr_SetString(PyExc_TypeError, "value 'varSep > 0' in 'MonomInterface");
      return -1;
    }
    else {
      if (orderType == IMonomInterface::Elim) {
        if (0 >= varSep || varSep >= independ->dim()) {
          PyErr_SetString(PyExc_TypeError, "'0 < varSep < len(independ)' in MonomInterface");
          return -1;
        }
      }
      else {
        if (depend) {
          if (0 >= varSep || varSep >= depend->dim()) {
            PyErr_SetString(PyExc_TypeError, "'0 < varSep < len(depend)' in MonomInterface");
            return -1;
          }
        }
        else {
          if (0 >= varSep || varSep >= systemType->module()) {
            PyErr_SetString(PyExc_TypeError, "'0 < varSep < systemType->module()' in MonomInterface");
            return -1;
          }
        }
      }
    }
  }
  self->monomInterface = IMonomInterface::create(orderType, systemType, independ, depend, varSep);

  return 0;
}

PyObject *monomInterface_getattr(MonomInterface *self, char *name) {
  PyObject *r = NULL;
  static PyMethodDef methods[] = {
    {"order", (PyCFunction)monomInterface_order, METH_NOARGS , ""},
    {"dimIndepend", (PyCFunction)monomInterface_dimIndepend, METH_NOARGS , ""},
    {"independ", (PyCFunction)monomInterface_independ, METH_NOARGS , ""},
    {"dimDepend", (PyCFunction)monomInterface_dimDepend, METH_NOARGS , ""},
    {"depend", (PyCFunction)monomInterface_depend, METH_NOARGS , ""},
    {"varSep", (PyCFunction)monomInterface_varSep, METH_NOARGS , ""},
    {NULL}
  };
  r = Py_FindMethod(methods, (PyObject*)self, name);

  if (r == NULL)
    PyErr_SetString(PyExc_TypeError, "Cannot definition operation in 'MonomInterface'");

  return r;
}

PyObject* monomInterface_order(MonomInterface *self) {
  PyObject *r = NULL;
  switch(self->monomInterface->order()) {
    case IMonomInterface::DegRevLex:
      r = PyString_FromString("DegRevLex");
      break;
    case IMonomInterface::DegRevLexByte:
      r = PyString_FromString("DegRevLexByte");
      break;
    case IMonomInterface::Lex:
      r = PyString_FromString("Lex");
      break;
    case IMonomInterface::Elim:
      r = PyString_FromString("Elim");
      break;
    case IMonomInterface::TopDegRevLexByte:
      r = PyString_FromString("TopDegRevLexByte");
      break;
    case IMonomInterface::TopDegRevLex:
      r = PyString_FromString("TopDegRevLex");
      break;
    case IMonomInterface::TopElim:
      r = PyString_FromString("TopElim");
      break;
    case IMonomInterface::TopLex:
      r = PyString_FromString("TopLex");
      break;
    case IMonomInterface::PotDegRevLex:
      r = PyString_FromString("PotDegRevLex");
      break;
    case IMonomInterface::PotLex:
      r = PyString_FromString("PotLex");
      break;
    case IMonomInterface::PosElim:
      r = PyString_FromString("PosElim");
      break;
  }
  return r;
}

PyObject* monomInterface_dimIndepend(MonomInterface *self) {
  return PyInt_FromLong(self->monomInterface->dimIndepend());
}

PyObject* monomInterface_independ(MonomInterface *self) {
  PyObject *r = NULL;
  int n = self->monomInterface->dimIndepend();
  r = PyList_New(n);
  IVariables::ConstIterator i = self->monomInterface->independ()->begin();
  int k = 0;
  while (k < n) {
    PyObject *var = PyString_FromString(*i);
    PyList_SetItem(r, k, var);
    ++i;
    ++k;
  }
  return r;
}

PyObject* monomInterface_dimDepend(MonomInterface *self) {
  return PyInt_FromLong(self->monomInterface->dimDepend());
}

PyObject* monomInterface_depend(MonomInterface *self) {
  PyObject *r = NULL;
  int n = self->monomInterface->dimDepend();
  r = PyList_New(n);
  if (n) {
    IVariables::ConstIterator i = self->monomInterface->depend()->begin();
    int k = 0;
    while (k < n) {
      PyObject *var = PyString_FromString(*i);
      PyList_SetItem(r, k, var);
      ++i;
      ++k;
    }
  }
  return r;
}

PyObject* monomInterface_varSep(MonomInterface *self) {
  return PyInt_FromLong(self->monomInterface->elim());
}

PyTypeObject monomInterface_type = {
  PyObject_HEAD_INIT(NULL)
  0,                                            /*ob_size*/
  "MonomInterface",                             /*tp_name*/
  sizeof(MonomInterface),                       /*tp_size*/
  0,                                            /*tp_itemsize*/
  (destructor)monomInterface_dealloc,           /*tp_dealloc*/
  0,                                            /*tp_print*/
  (getattrfunc)monomInterface_getattr,          /*tp_getattr*/
  0,                                            /*tp_setattr*/
  0,                                            /*tp_compare*/
  0,                                            /*tp_repr*/
  0,                                            /*tp_as_number*/
  0,                                            /*tp_as_sequence*/
  0,                                            /*tp_as_mapping*/
  0,                                            /*tp_hash*/
  0,                                            /*tp_call*/
  0,                                            /*tp_str*/
  0,                                            /*tp_getattro*/
  0,                                            /*tp_setattro*/
  0,                                            /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,     /*tp_flags*/
  "MonomInterface object",                      /*tp_doc */
  0,                                            /*tp_traverse */
  0,                                            /*tp_clear */
  0,                                            /*tp_richcompare */
  0,                                            /*tp_weaklistoffset */
  0,                                            /*tp_iter */
  0,                                            /*tp_iternext */
  0,                                            /*tp_methods */
  0,                                            /*tp_members */
  0,                                            /*tp_getset */
  0,                                            /*tp_base */
  0,                                            /*tp_dict */
  0,                                            /*tp_descr_get */
  0,                                            /*tp_descr_set */
  0,                                            /*tp_dictoffset */
  (initproc)monomInterface_init,                /*tp_init */
  0,                                            /*tp_alloc */
  (newfunc)monomInterface_new,                  /*tp_new */
};

void monom_dealloc(Monom *self) {
//   delete self->monom;
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* monom_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  Monom *self = (Monom *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->monom = NULL;
  }

  return (PyObject *)self;
}

int monom_init(Monom *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"monomInterface", "monom", NULL};

  PyObject *py_monomInterface = NULL;
  const char *py_in = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|s!:Monom", kwlist,
                                    &monomInterface_type, &py_monomInterface,
                                    &py_in))
      return -1;

  self->monom = ((MonomInterface*)py_monomInterface)->monomInterface->create();
  self->monom->setZero();
  if (py_in) {
    std::istringstream in(py_in);
    in >> *self->monom >> std::ws;
    if (!in.eof()) {
      PyErr_SetString(PyExc_TypeError, "error parsing second argument");
      return -1;
    }
  }

  return 0;
}

int monom_print(Monom *m, FILE *fp, int flags) {
  std::ostringstream out;
  out << *m->monom;
  fputs(out.str().c_str(), fp);
  return 0;
}

PyObject* monom_str(Monom *m) {
  std::ostringstream out;
  out << *m->monom;
  return PyString_FromString(out.str().c_str());
}

PyObject* monom_getattr(Monom *self, char *name) {
  PyObject *r = NULL;

  static PyMethodDef methods[] = {
    {"degree", (PyCFunction)monom_degree, METH_NOARGS, ""},
    {"dependVar", (PyCFunction)monom_dependVar, METH_NOARGS, ""},
    {"setZero", (PyCFunction)monom_setZero, METH_NOARGS, ""},
    {"prolong", (PyCFunction)monom_prolong, METH_VARARGS, ""},
    {"gcd", (PyCFunction)monom_gcd, METH_VARARGS, ""},
    {"lcm", (PyCFunction)monom_lcm, METH_VARARGS, ""},
    {"divisibility", (PyCFunction)monom_divisibility, METH_VARARGS, ""},
    {"divisibilityTrue", (PyCFunction)monom_divisibilityTrue, METH_VARARGS, ""},
    {NULL}
  };
  r = Py_FindMethod(methods, (PyObject*)self, name);

  if (r == NULL)
    PyErr_SetString(PyExc_TypeError, "Cannot definition operation in 'Monom'");

  return r;
}

PyObject* monom_degree(Monom *self) {
  return PyInt_FromLong(self->monom->degree());
}

PyObject* monom_dependVar(Monom *self) {
  return PyInt_FromLong(self->monom->dependVar());
}

int monom_compare(Monom *o1, PyObject *o2) {
  int r = -1;
  if (!PyObject_TypeCheck(o2, &monom_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Monomial");
  else {
    Monom *monom = (Monom*)o2;
    if (o1->monom->monomInterface() != monom->monom->monomInterface())
      PyErr_SetString(PyExc_TypeError, "Monomial with difference MonomInterface");
    else
      r = o1->monom->compare(*monom->monom);
  }
  return r;
}

PyObject* monom_richCompare(Monom *o1, PyObject *o2, int opid) {
  PyObject *r = NULL;
  if (!PyObject_TypeCheck(o2, &monom_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Monomial");
  else {
    Monom *monom = (Monom*)o2;
    if (o1->monom->monomInterface() != monom->monom->monomInterface())
      PyErr_SetString(PyExc_TypeError, "Monomial with difference MonomInterface");
    else {
      int cmp = o1->monom->compare(*monom->monom);
      bool res = -1;
      switch (opid) {
        case Py_LT:
          res = (cmp == -1);
          break;
        case Py_LE:
          res = (cmp == -1) || (cmp == 0);
          break;
        case Py_EQ:
          res = (cmp == 0);
          break;
        case Py_NE:
          res = (cmp != 0);
          break;
        case Py_GT:
          res = (cmp == 1);
          break;
        case Py_GE:
          res = (cmp == 1) || (cmp == 0);
          break;
      }
      if (res)
        r = Py_True;
      else
        r = Py_False;
      Py_INCREF(r);
    }
  }
  return r;
}

PyObject* monom_setZero(Monom *self) {
  PyObject *r = NULL;
  self->monom->setZero();
  Py_INCREF(self);
  r = (PyObject*)self;

  return r;
}

PyObject* monom_prolong(Monom *self, PyObject *args) {
  int var = 0;
  unsigned deg = 1;
  if (!PyArg_ParseTuple(args, "i|I:Monom",
                        &var, &deg))
      return NULL;

  PyObject *r = NULL;
  if (0 > var || var >= self->monom->dimIndepend())
    PyErr_SetString(PyExc_TypeError, "'0 <= var < len(independ)' in MonomInterface");
  else {
    self->monom->prolong(var, deg);
    Py_INCREF(self);
    r = (PyObject*)self;
  }

  return r;
}

PyObject* monom_gcd(Monom *self, PyObject *args) {
  PyObject *py_monom = NULL;
  if (!PyArg_ParseTuple(args, "O!:Monom",
                        &monom_type, &py_monom))
      return NULL;

  PyObject *r = NULL;
  Monom *monom = (Monom*)py_monom;
  if (self->monom->monomInterface() != monom->monom->monomInterface())
    PyErr_SetString(PyExc_TypeError, "Monomial with difference MonomInterface");
  else
    r = PyInt_FromLong(self->monom->gcd(*monom->monom));

  return r;
}

PyObject* monom_lcm(Monom *self, PyObject *args) {
  PyObject *py_monom = NULL;
  if (!PyArg_ParseTuple(args, "O!:Monom",
                        &monom_type, &py_monom))
      return NULL;

  PyObject *r = NULL;
  Monom *monom = (Monom*)py_monom;
  if (self->monom->monomInterface() != monom->monom->monomInterface())
    PyErr_SetString(PyExc_TypeError, "Monomial with difference MonomInterface");
  else {
    IMonom *m = self->monom->monomInterface()->create();
    m->lcm(*self->monom, *monom->monom);
    Monom* res = PyObject_New(Monom, &monom_type);
    res->monom = m;
    r = (PyObject*)res;
  }

  return r;
}

PyObject* monom_divisibility(Monom *self, PyObject *args) {
  PyObject *py_monom = NULL;
  if (!PyArg_ParseTuple(args, "O!:Monom",
                        &monom_type, &py_monom))
      return NULL;

  PyObject *r = NULL;
  Monom *monom = (Monom*)py_monom;
  if (self->monom->monomInterface() != monom->monom->monomInterface())
    PyErr_SetString(PyExc_TypeError, "Monomial with difference MonomInterface");
  else {
    if (self->monom->divisibility(*monom->monom))
      r = Py_True;
    else
      r = Py_False;
    Py_INCREF(r);
  }

  return r;
}

PyObject* monom_divisibilityTrue(Monom *self, PyObject *args) {
  PyObject *py_monom = NULL;
  if (!PyArg_ParseTuple(args, "O!:Monom",
                        &monom_type, &py_monom))
      return NULL;

  PyObject *r = NULL;
  Monom *monom = (Monom*)py_monom;
  if (self->monom->monomInterface() != monom->monom->monomInterface())
    PyErr_SetString(PyExc_TypeError, "Monomial with difference MonomInterface");
  else {
    if (self->monom->divisibilityTrue(*monom->monom))
      r = Py_True;
    else
      r = Py_False;
    Py_INCREF(r);
  }

  return r;
}

PyObject* monom_getiter(Monom *self) {
  MonomIterator* iter = PyObject_New(MonomIterator, &monomIterator_type);
  iter->monom = self->monom;
  iter->iter = 0;
  return (PyObject*)iter;
}

PyObject* monom_multiply(Monom *o1, PyObject* o2) {
  PyObject *r = NULL;
  if (!PyObject_TypeCheck(o2, &monom_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Monomial");
  else {
    Monom *monom = (Monom*)o2;
    if (o1->monom->monomInterface() != monom->monom->monomInterface())
      PyErr_SetString(PyExc_TypeError, "Monomial with difference MonomInterface");
    else {
      IMonom *m = o1->monom->monomInterface()->create();
      m->mult(*o1->monom, *monom->monom);
      Monom* res = PyObject_New(Monom, &monom_type);
      res->monom = m;
      r = (PyObject*)res;
    }
  }
  return r;
}

PyObject* monom_divide(Monom *o1, PyObject* o2) {
  PyObject *r = NULL;
  if (!PyObject_TypeCheck(o2, &monom_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Monomial");
  else {
    Monom *monom = (Monom*)o2;
    if (o1->monom->monomInterface() != monom->monom->monomInterface())
      PyErr_SetString(PyExc_TypeError, "Monomial with difference MonomInterface");
    else if (!o1->monom->divisibility(*monom->monom))
      PyErr_SetString(PyExc_TypeError, "not divisibility");
    else {
      IMonom *m = o1->monom->monomInterface()->create();
      m->divide(*o1->monom, *monom->monom);
      Monom* res = PyObject_New(Monom, &monom_type);
      res->monom = m;
      r = (PyObject*)res;
    }
  }
  return r;
}

int monom_nonzero(Monom *o1) {
  return o1->monom->degree() > 0;
}

PyObject* monom_inplace_multiply(Monom *o1, PyObject* o2) {
  PyObject *r = NULL;
  if (!PyObject_TypeCheck(o2, &monom_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Monomial");
  else {
    Monom *monom = (Monom*)o2;
    if (o1->monom->monomInterface() != monom->monom->monomInterface())
      PyErr_SetString(PyExc_TypeError, "Monomial with difference MonomInterface");
    else {
      o1->monom->mult(*monom->monom);
      r = (PyObject*)o1;
      Py_INCREF(r);
    }
  }
  return r;
}

#if PY_VERSION_HEX < 0x02050000
int monom_length(Monom *self) {
#else
Py_ssize_t monom_length(Monom *self) {
#endif
  return self->monom->dimIndepend();
}

#if PY_VERSION_HEX < 0x02050000
PyObject* monom_item(Monom *self, int i) {
#else
PyObject* monom_item(Monom *self, Py_ssize_t i) {
#endif
  PyObject* r = NULL;
  if (0 > i || i >= self->monom->dimIndepend())
    PyErr_SetString(PyExc_ValueError, "Index out of range.");
  else
    r = PyInt_FromLong(self->monom->deg(i));
  return r;
}

static PyNumberMethods monom_as_number = {
  0,                                    /*nb_add*/
  0,                                    /*nb_subtract*/
  (binaryfunc)monom_multiply,           /*nb_multiply*/
  (binaryfunc)monom_divide,             /*nb_divide*/
  0,                                    /*nb_remainder*/
  0,                                    /*nb_divmod*/
  0,                                    /*nb_power*/
  0,                                    /*nb_negative*/
  0,                                    /*nb_positive*/
  0,                                    /*nb_absolute*/
  (inquiry)monom_nonzero,               /*nb_nonzero*/
  0,                                    /*nb_invert*/
  0,                                    /*nb_lshift*/
  0,                                    /*nb_rshift*/
  0,                                    /*nb_and*/
  0,                                    /*nb_xor*/
  0,                                    /*nb_or*/
  0,                                    /*nb_coerce*/
  0,                                    /*nb_int*/
  0,                                    /*nb_long*/
  0,                                    /*nb_float*/
  0,                                    /*nb_oct*/
  0,                                    /*nb_hex*/
  0,                                    /*nb_inplace_add*/
  0,                                    /*nb_inplace_subtract*/
  (binaryfunc)monom_inplace_multiply,   /*nb_inplace_multiply*/
  0,                                    /*nb_inplace_divide*/
  0,                                    /*nb_inplace_remainder*/
  0,                                    /*nb_inplace_power*/
  0,                                    /*nb_inplace_lshift*/
  0,                                    /*nb_inplace_rshift*/
  0,                                    /*nb_inplace_and*/
  0,                                    /*nb_inplace_xor*/
  0,                                    /*nb_inplace_or*/
  0,                                    /*nb_floor_divide*/
  0,                                    /*nb_true_divide*/
  0,                                    /*nb_inplace_floor_divide*/
  0,                                    /*nb_inplace_true_divide*/
};

static PySequenceMethods monom_sequence = {
#if PY_VERSION_HEX < 0x02050000
  (inquiry)monom_length,                /*sq_length*/
#else
  (lenfunc)monom_length,                /*sq_length*/
#endif
  0,                                    /*sq_concat*/
  0,                                    /*sq_repeat*/
#if PY_VERSION_HEX < 0x02050000
  (intargfunc)monom_item,               /*sq_item*/
#else
  (ssizeargfunc)monom_item,               /*sq_item*/
#endif
  0,                                    /*sq_slice*/
  0,                                    /*sq_ass_item*/
  0,                                    /*sq_ass_slice*/
  0,                                    /*sq_contains */
  0,                                    /*sq_inplace_concat*/
  0,                                    /*sq_inplace_repeat*/
};

PyTypeObject monom_type = {
  PyObject_HEAD_INIT(NULL)
  0,                                            /*ob_size*/
  "Monom",                                      /*tp_name*/
  sizeof(Monom),                                /*tp_size*/
  0,                                            /*tp_itemsize*/
  (destructor)monom_dealloc,                    /*tp_dealloc*/
  (printfunc)monom_print,                       /*tp_print*/
  (getattrfunc)monom_getattr,                   /*tp_getattr*/
  0,                                            /*tp_setattr*/
  (cmpfunc)monom_compare,                       /*tp_compare*/
  (reprfunc)monom_str,                          /*tp_repr*/
  &monom_as_number,                             /*tp_as_number*/
  &monom_sequence,                              /*tp_as_sequence*/
  0,                                            /*tp_as_mapping*/
  0,                                            /*tp_hash*/
  0,                                            /*tp_call*/
  (reprfunc)monom_str,                          /*tp_str*/
  0,                                            /*tp_getattro*/
  0,                                            /*tp_setattro*/
  0,                                            /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,     /*tp_flags*/
  "Monom object",                               /*tp_doc */
  0,                                            /*tp_traverse */
  0,                                            /*tp_clear */
  (richcmpfunc)monom_richCompare,               /*tp_richcompare */
  0,                                            /*tp_weaklistoffset */
  (getiterfunc)monom_getiter,                   /*tp_iter */
  0,                                            /*tp_iternext */
  0,                                            /*tp_methods */
  0,                                            /*tp_members */
  0,                                            /*tp_getset */
  0,                                            /*tp_base */
  0,                                            /*tp_dict */
  0,                                            /*tp_descr_get */
  0,                                            /*tp_descr_set */
  0,                                            /*tp_dictoffset */
  (initproc)monom_init,                         /*tp_init */
  0,                                            /*tp_alloc */
  (newfunc)monom_new,                           /*tp_new */
};

void monomIterator_dealloc(MonomIterator *self) {
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* monomIterator_next(MonomIterator *self) {
  PyObject* r = NULL;
  if (self->iter < self->monom->monomInterface()->dimIndepend()) {
/*    r = PyInt_FromLong((*self->monom)[self->iter]);
    r = PyInt_FromLong(self->iter);*/
    r = PyInt_FromLong(self->monom->deg(self->iter));
    ++self->iter;
  }
  return r;
}

PyTypeObject monomIterator_type = {
  PyObject_HEAD_INIT(NULL)
  0,                                            /*ob_size*/
  "MonomIterator",                              /*tp_name*/
  sizeof(MonomIterator),                        /*tp_size*/
  0,                                            /*tp_itemsize*/
  (destructor)monomIterator_dealloc,            /*tp_dealloc*/
  0,                                            /*tp_print*/
  0,                                            /*tp_getattr*/
  0,                                            /*tp_setattr*/
  0,                                            /*tp_compare*/
  0,                                            /*tp_repr*/
  0,                                            /*tp_as_number*/
  0,                                            /*tp_as_sequence*/
  0,                                            /*tp_as_mapping*/
  0,                                            /*tp_hash*/
  0,                                            /*tp_call*/
  0,                                            /*tp_str*/
  0,                                            /*tp_getattro*/
  0,                                            /*tp_setattro*/
  0,                                            /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,     /*tp_flags*/
  "MonomIterator object",                       /*tp_doc */
  0,                                            /*tp_traverse */
  0,                                            /*tp_clear */
  0,                                            /*tp_richcompare */
  0,                                            /*tp_weaklistoffset */
  0,                                            /*tp_iter */
  (iternextfunc)monomIterator_next,             /*tp_iternext */
  0,                                            /*tp_methods */
  0,                                            /*tp_members */
  0,                                            /*tp_getset */
  0,                                            /*tp_base */
  0,                                            /*tp_dict */
  0,                                            /*tp_descr_get */
  0,                                            /*tp_descr_set */
  0,                                            /*tp_dictoffset */
  0,                                            /*tp_init */
  0,                                            /*tp_alloc */
  0,                                            /*tp_new */
};
