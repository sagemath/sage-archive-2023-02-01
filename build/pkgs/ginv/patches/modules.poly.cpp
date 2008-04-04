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
#include "./coeff.h"
#include "./poly.h"

#include <string.h>
#include <sstream>

#include "../ginv/poly/iexpression.h"
#include "../ginv/poly/ipolylist.h"

void polyInterface_dealloc(PolyInterface *self) {
//   delete self->polyInterface;
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* polyInterface_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  PolyInterface *self = (PolyInterface*)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->polyInterface = NULL;
  }

  return (PyObject *)self;
}

int polyInterface_init(PolyInterface *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"type", "systemType", "monomInterface", "coeffInterface", NULL};

  const char *py_type = NULL;
  PyObject *py_systemType = NULL;
  PyObject *py_monomInterface = NULL;
  PyObject *py_coeffInterface = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "sO!O!O!:PolyInterface", kwlist,
                                    &py_type,
                                    &systemType_type, &py_systemType,
                                    &monomInterface_type, &py_monomInterface,
                                    &coeffInterface_type, &py_coeffInterface))
      return -1;

  IPolyInterface::Type type = IPolyInterface::Indefinitely;
  if (strcmp("PolyList", py_type) == 0)
    type = IPolyInterface::PolyList;
  else {
    PyErr_SetString(PyExc_TypeError, "unknown type in 'PolyInterface'");
    return -1;
  }

  ISystemType *systemType = ((SystemType *)py_systemType)->systemType;
  IMonomInterface *monomInterface = ((MonomInterface *)py_monomInterface)->monomInterface;
  ICoeffInterface *coeffInterface = ((CoeffInterface *)py_coeffInterface)->coeffInterface;

  self->polyInterface = IPolyInterface::create(type, systemType, monomInterface, coeffInterface);

  return 0;
}

PyObject *polyInterface_getattr(PolyInterface *self, char *name) {
 PyObject *r = NULL;
  static PyMethodDef methods[] = {
    {"type", (PyCFunction)polyInterface_getType, METH_NOARGS , ""},
    {NULL}
  };
  r = Py_FindMethod(methods, (PyObject*)self, name);

  if (r == NULL)
    PyErr_SetString(PyExc_TypeError, "Cannot definition operation in 'PolyInterface'");

  return r;
}

PyObject* polyInterface_getType(PolyInterface *self) {
  PyObject *r = NULL;
  switch(self->polyInterface->type()) {
    case IPolyInterface::PolyList:
      r = PyString_FromString("PolyList");
      break;
  }
  return r;
}

PyTypeObject polyInterface_type = {
  PyObject_HEAD_INIT(NULL)
  0,                                            /*ob_size*/
  "PolyInterface",                              /*tp_name*/
  sizeof(PolyInterface),                        /*tp_size*/
  0,                                            /*tp_itemsize*/
  (destructor)polyInterface_dealloc,            /*tp_dealloc*/
  0,                                            /*tp_print*/
  (getattrfunc)polyInterface_getattr,           /*tp_getattr*/
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
  "PolyInterface object",                       /*tp_doc */
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
  (initproc)polyInterface_init,                 /*tp_init */
  0,                                            /*tp_alloc */
  (newfunc)polyInterface_new,                   /*tp_new */
};

void poly_dealloc(Poly *self) {
//   delete self->poly;
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* poly_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  Poly *self = (Poly*)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->poly = NULL;
  }

  return (PyObject *)self;
}

int poly_init(Poly *self, PyObject *args, PyObject *kwds) {
  static char *kwlist[] = {"polyInterface", "poly", NULL};

  PyObject *py_polyInterface = NULL;
  const char *py_in = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|s:Poly", kwlist,
       &polyInterface_type, &py_polyInterface,
       &py_in))
    return -1;

  if (py_in == NULL)
    self->poly = ((PolyInterface*)py_polyInterface)->polyInterface->create();
  else {
    std::istringstream in(py_in);
    IExpression expr(((PolyInterface*)py_polyInterface)->polyInterface);
    in >> expr >> std::ws;
    if (!in.eof()) {
      PyErr_SetString(PyExc_TypeError, "error parsing second argument");
      return -1;
    }
    self->poly = expr.toPoly();
  }

  return 0;
}

int poly_print(Poly *p, FILE *fp, int flags) {
  std::ostringstream out;
  out << *p->poly;
  fputs(out.str().c_str(), fp);
  return 0;
}

PyObject* poly_str(Poly *p) {
  std::ostringstream out;
  out << *p->poly;
  return PyString_FromString(out.str().c_str());
}

PyObject* poly_getattr(Poly *self, char *name) {
  PyObject *r = NULL;
  static PyMethodDef methods[] = {
    {"length", (PyCFunction)poly_getLength, METH_NOARGS, ""},
    {"degree", (PyCFunction)poly_degree, METH_NOARGS, ""},
    {"norm", (PyCFunction)poly_norm, METH_NOARGS, ""},
    {"lm", (PyCFunction)poly_lm, METH_NOARGS, ""},
    {"lc", (PyCFunction)poly_lc, METH_NOARGS, ""},
    {"setZero", (PyCFunction)poly_setZero, METH_NOARGS, ""},
    {"pp", (PyCFunction)poly_pp, METH_NOARGS, ""},
    {"isPp", (PyCFunction)poly_isPp, METH_NOARGS, ""},
    {"prolong", (PyCFunction)poly_prolong, METH_VARARGS, ""},
    {"mult", (PyCFunction)poly_mult, METH_VARARGS, ""},
    {"diff", (PyCFunction)poly_diff, METH_VARARGS, ""},
    {"shift", (PyCFunction)poly_shift, METH_VARARGS, ""},
    {"reduction", (PyCFunction)poly_reduction, METH_VARARGS, ""},
    {"spoly", (PyCFunction)poly_spoly, METH_VARARGS, ""},
    {NULL}
  };
  r = Py_FindMethod(methods, (PyObject*)self, name);
  if (r == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot definition operation in 'Poly'");
  }

  return r;
}

int poly_compare(Poly *o1, PyObject *o2) {
  int r = -1;
  if (!PyObject_TypeCheck(o2, &poly_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Polynomial");
  else {
    Poly *poly = (Poly*)o2;
    if (o1->poly->polyInterface() != poly->poly->polyInterface())
      PyErr_SetString(PyExc_TypeError, "Polynomial with difference PolyInterface");
    else {
      if (o1->poly->isZero())
        if (poly->poly->isZero())
          r = 0;
        else
          r = -1;
      else if (poly->poly->isZero())
        r = 1;
      else
        r = o1->poly->lm().compare(poly->poly->lm());
    }
  }
  return r;
}

PyObject* poly_richCompare(Poly *o1, PyObject *o2, int opid) {
  PyObject *r = NULL;
  if (!PyObject_TypeCheck(o2, &poly_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Polynomial");
  else {
    Poly *poly = (Poly*)o2;
    if (o1->poly->polyInterface() != poly->poly->polyInterface())
      PyErr_SetString(PyExc_TypeError, "Polynomial with difference PolyInterface");
    else {
      int cmp = -1;
      if (o1->poly->isZero())
        if (poly->poly->isZero())
          cmp = 0;
        else
          cmp = -1;
      else if (poly->poly->isZero())
        cmp = 1;
      else
        cmp = o1->poly->lm().compare(poly->poly->lm());

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

PyObject* poly_isZero(Poly *self) {
  PyObject *r = NULL;
  if (self->poly->isZero())
    r = Py_True;
  else
    r = Py_False;
  Py_INCREF(r);
  return r;
}

PyObject* poly_getLength(Poly *self) {
  return PyInt_FromLong(self->poly->length());
}

PyObject* poly_degree(Poly *self) {
  return PyInt_FromLong(self->poly->degree());
}

PyObject* poly_norm(Poly *self) {
  return PyInt_FromLong(self->poly->norm());
}

PyObject* poly_lm(Poly *self) {
  PyObject *r = NULL;
  if (self->poly->isZero())
    PyErr_SetString(PyExc_TypeError, "Polynomial is zero");
  else {
    IMonom *m = self->poly->monomInterface()->copy(self->poly->lm());
    Monom* res = PyObject_New(Monom, &monom_type);
    res->monom = m;
    r = (PyObject*)res;
  }
  return r;
}

PyObject* poly_lc(Poly *self) {
  PyObject *r = NULL;
  if (self->poly->isZero())
    PyErr_SetString(PyExc_TypeError, "Polynomial is zero");
  else {
    ICoeff *c = self->poly->coeffInterface()->copy(self->poly->lc());
    Coeff* res = PyObject_New(Coeff, &coeff_type);
    res->coeff = c;
    r = (PyObject*)res;
  }
  return r;
}

PyObject* poly_setZero(Poly *self) {
 PyObject *r = NULL;
  self->poly->setZero();
  Py_INCREF(self);
  r = (PyObject*)self;

  return r;
}

PyObject* poly_pp(Poly *self) {
  PyObject *r = NULL;
  if (!self->poly->isZero())
    self->poly->pp();
  Py_INCREF(self);
  r = (PyObject*)self;

  return r;
}

PyObject* poly_isPp(Poly *self) {
  PyObject *r = NULL;
  if (self->poly->isZero() || self->poly->isPp())
    r = Py_True;
  else
    r = Py_False;
  Py_INCREF(r);

  return r;
}

PyObject* poly_prolong(Poly *self, PyObject *args) {
  int var = 0;
  unsigned deg = 1;
  if (!PyArg_ParseTuple(args, "i|I:Poly",
                        &var, &deg))
      return NULL;

  PyObject *r = NULL;
  if (0 > var || var >= self->poly->monomInterface()->dimIndepend())
    PyErr_SetString(PyExc_TypeError, "'0 <= var < len(independ)' in MonomInterface");
  else {
    self->poly->mult(var, deg);
    Py_INCREF(self);
    r = (PyObject*)self;
  }

  return r;
}

PyObject* poly_mult(Poly *self, PyObject *args) {
  PyObject *py_mult = NULL;
  if (!PyArg_ParseTuple(args, "O:Poly",
                        &py_mult))
      return NULL;

  PyObject *r = NULL;
  if (PyObject_TypeCheck(py_mult, &monom_type)) {
    Monom *monom = (Monom*)py_mult;
    if (self->poly->monomInterface() != monom->monom->monomInterface())
      PyErr_SetString(PyExc_TypeError, "Monomial with difference MonomInterface");
    else {
      self->poly->mult(*monom->monom);
      Py_INCREF(self);
      r = (PyObject*)self;
    }
  }
  else if (PyObject_TypeCheck(py_mult, &coeff_type)) {
    Coeff *coeff = (Coeff*)py_mult;
    if (self->poly->coeffInterface() != coeff->coeff->coeffInterface())
      PyErr_SetString(PyExc_TypeError, "Coeff with difference CoeffInterface");
    else {
      self->poly->mult(*coeff->coeff);
      Py_INCREF(self);
      r = (PyObject*)self;
    }
  }

  return r;
}

PyObject* poly_diff(Poly *self, PyObject *args) {
  PyObject *r = NULL;
  PyErr_SetString(PyExc_TypeError, "not implementation");
  return r;
}

PyObject* poly_shift(Poly *self, PyObject *args) {
  PyObject *r = NULL;
  PyErr_SetString(PyExc_TypeError, "not implementation");
  return r;
}

PyObject* poly_reduction(Poly *self, PyObject *args) {
  PyObject *py_poly = NULL;
  if (!PyArg_ParseTuple(args, "O!:Poly",
                        &poly_type, &py_poly))
      return NULL;

  PyObject *r = NULL;
  Poly *poly = (Poly*)py_poly;
  if (self->poly->polyInterface() != poly->poly->polyInterface())
    PyErr_SetString(PyExc_TypeError, "Polynomial with difference PolyInterface");
  else if (self->poly->isZero() || poly->poly->isZero()) {
    Py_INCREF(self);
    r = (PyObject*)self;
  }
  else if (self->poly->lm().divisibility(poly->poly->lm())) {
    self->poly->reduction(*poly->poly);
    Py_INCREF(self);
    r = (PyObject*)self;
  }
  else
    PyErr_SetString(PyExc_TypeError, "reduction impossible");

  return r;
}

PyObject* poly_spoly(Poly *self, PyObject *args) {
  PyObject *r = NULL;
  PyErr_SetString(PyExc_TypeError, "not implementation");
  return r;
}

PyObject *poly_getiter(Poly *self) {
  PolyIterator* iter = PyObject_New(PolyIterator, &polyIterator_type);
  iter->iter = ((const IPoly*)self->poly)->begin();
  return (PyObject*)iter;
}


PyObject* poly_add(Poly *o1, PyObject* o2) {
  PyObject *r = NULL;
  if (!PyObject_TypeCheck(o2, &poly_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Poly");
  else {
    Poly *poly = (Poly*)o2;
    if (o1->poly->polyInterface() != poly->poly->polyInterface())
      PyErr_SetString(PyExc_TypeError, "Poly with difference PolyInterface");
    else {
      IPoly *p1 = o1->poly->polyInterface()->copy(*o1->poly);
      IPoly *p2 = o1->poly->polyInterface()->copy(*poly->poly);
      p1->add(*p2);
      Poly* res = PyObject_New(Poly, &poly_type);
      res->poly = p1;
      r = (PyObject*)res;
      delete p2;
    }
  }
  return r;
}

PyObject* poly_sub(Poly *o1, PyObject* o2) {
  PyObject *r = NULL;
  if (!PyObject_TypeCheck(o2, &poly_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Poly");
  else {
    Poly *poly = (Poly*)o2;
    if (o1->poly->polyInterface() != poly->poly->polyInterface())
      PyErr_SetString(PyExc_TypeError, "Poly with difference PolyInterface");
    else {
      IPoly *p1 = o1->poly->polyInterface()->copy(*o1->poly);
      IPoly *p2 = o1->poly->polyInterface()->copy(*poly->poly);
      p1->sub(*p2);
      Poly* res = PyObject_New(Poly, &poly_type);
      res->poly = p1;
      r = (PyObject*)res;
      delete p2;
    }
  }
  return r;
}

PyObject* poly_multiply(Poly *o1, PyObject* o2) {
  PyObject *r = NULL;
  if (!PyObject_TypeCheck(o2, &poly_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Poly");
  else {
    Poly *poly = (Poly*)o2;
    if (o1->poly->polyInterface() != poly->poly->polyInterface())
      PyErr_SetString(PyExc_TypeError, "Poly with difference PolyInterface");
    else {
      IPoly *p = o1->poly->polyInterface()->copy(*o1->poly);
      p->mult(*poly->poly);
      Poly* res = PyObject_New(Poly, &poly_type);
      res->poly = p;
      r = (PyObject*)res;
    }
  }
  return r;
}

PyObject* poly_negative(Poly *o1) {
  PyObject *r = NULL;
  IPoly *p = o1->poly->polyInterface()->copy(*o1->poly);
  p->minus();
  Poly* res = PyObject_New(Poly, &poly_type);
  res->poly = p;
  r = (PyObject*)res;
  return r;
}

int poly_nonzero(Poly *o1) {
  return !o1->poly->isZero();
}

PyObject* poly_inplace_add(Poly *o1, PyObject* o2) {
  PyObject *r = NULL;
  if (!PyObject_TypeCheck(o2, &poly_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Poly");
  else {
    Poly *poly = (Poly*)o2;
    if (o1->poly->polyInterface() != poly->poly->polyInterface())
      PyErr_SetString(PyExc_TypeError, "Poly with difference PolyInterface");
    else {
      IPoly *p = o1->poly->polyInterface()->copy(*poly->poly);
      o1->poly->add(*p);
      r = (PyObject*)o1;
      Py_INCREF(r);
      delete p;
    }
  }
  return r;
}

PyObject* poly_inplace_sub(Poly *o1, PyObject* o2) {
  PyObject *r = NULL;
  if (!PyObject_TypeCheck(o2, &poly_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Poly");
  else {
    Poly *poly = (Poly*)o2;
    if (o1->poly->polyInterface() != poly->poly->polyInterface())
      PyErr_SetString(PyExc_TypeError, "Poly with difference PolyInterface");
    else {
      IPoly *p = o1->poly->polyInterface()->copy(*poly->poly);
      o1->poly->sub(*p);
      r = (PyObject*)o1;
      Py_INCREF(r);
      delete p;
    }
  }
  return r;
}

PyObject* poly_inplace_multiply(Poly *o1, PyObject* o2) {
  PyObject *r = NULL;
  if (!PyObject_TypeCheck(o2, &poly_type))
    PyErr_SetString(PyExc_TypeError, "second argument not Poly");
  else {
    Poly *poly = (Poly*)o2;
    if (o1->poly->polyInterface() != poly->poly->polyInterface())
      PyErr_SetString(PyExc_TypeError, "Poly with difference PolyInterface");
    else {
      o1->poly->mult(*poly->poly);
      r = (PyObject*)o1;
      Py_INCREF(r);
    }
  }
  return r;
}


#if PY_VERSION_HEX < 0x02050000
int poly_length(Poly *self) {
#else
Py_ssize_t poly_length(Poly *self) {
#endif
  return self->poly->length();
}

#if PY_VERSION_HEX < 0x02050000
PyObject* poly_item(Poly *self, int i) {
#else
PyObject* poly_item(Poly *self, Py_ssize_t i) {
#endif
  PyObject* r = NULL;
  if (0 > i || i >= self->poly->length())
    PyErr_SetString(PyExc_ValueError, "Index out of range.");
  else {
    IPoly::ConstIterator *k = ((const IPoly*)self->poly)->begin();
    for(; i > 0; i--)
     ++(*k);

    Monom* m = PyObject_New(Monom, &monom_type);
    m->monom = (*k).monom().monomInterface()->copy((*k).monom());
    Coeff* c = PyObject_New(Coeff, &coeff_type);
    c->coeff = (*k).coeff().coeffInterface()->copy((*k).coeff());
    r = PyTuple_New(2);
    PyTuple_SetItem(r, 0, (PyObject*)m);
    PyTuple_SetItem(r, 1, (PyObject*)c);
  }
  return r;
}

static PyNumberMethods poly_as_number = {
  (binaryfunc)poly_add,                 /*nb_add*/
  (binaryfunc)poly_sub,                 /*nb_subtract*/
  (binaryfunc)poly_multiply,            /*nb_multiply*/
  0,                                    /*nb_divide*/
  0,                                    /*nb_remainder*/
  0,                                    /*nb_divmod*/
  0,                                    /*nb_power*/
  (unaryfunc)poly_negative,             /*nb_negative*/
  0,                                    /*nb_positive*/
  0,                                    /*nb_absolute*/
  (inquiry)poly_nonzero,                /*nb_nonzero*/
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
  (binaryfunc)poly_inplace_add,         /*nb_inplace_add*/
  (binaryfunc)poly_inplace_sub,         /*nb_inplace_subtract*/
  (binaryfunc)poly_inplace_multiply,    /*nb_inplace_multiply*/
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

static PySequenceMethods poly_sequence = {
#if PY_VERSION_HEX < 0x02050000
  (inquiry)poly_length,                 /*sq_length*/
#else
  (lenfunc)poly_length,                 /*sq_length*/
#endif
  0,                                    /*sq_concat*/
  0,                                    /*sq_repeat*/
#if PY_VERSION_HEX < 0x02050000
  (intargfunc)poly_item,                /*sq_item*/
#else
  (ssizeargfunc)poly_item,                /*sq_item*/
#endif
  0,                                    /*sq_slice*/
  0,                                    /*sq_ass_item*/
  0,                                    /*sq_ass_slice*/
  0,                                    /*sq_contains */
  0,                                    /*sq_inplace_concat*/
  0,                                    /*sq_inplace_repeat*/
};

PyTypeObject poly_type = {
  PyObject_HEAD_INIT(NULL)
  0,                                            /*ob_size*/
  "Poly",                                       /*tp_name*/
  sizeof(Poly),                                 /*tp_size*/
  0,                                            /*tp_itemsize*/
  (destructor)poly_dealloc,                     /*tp_dealloc*/
  (printfunc)poly_print,                        /*tp_print*/
  (getattrfunc)poly_getattr,                    /*tp_getattr*/
  0,                                            /*tp_setattr*/
  (cmpfunc)poly_compare,                        /*tp_compare*/
  (reprfunc)poly_str,                           /*tp_repr*/
  &poly_as_number,                              /*tp_as_number*/
  &poly_sequence,                               /*tp_as_sequence*/
  0,                                            /*tp_as_mapping*/
  0,                                            /*tp_hash*/
  0,                                            /*tp_call*/
  (reprfunc)poly_str,                           /*tp_str*/
  0,                                            /*tp_getattro*/
  0,                                            /*tp_setattro*/
  0,                                            /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,     /*tp_flags*/
  "Poly object",                                /*tp_doc */
  0,                                            /*tp_traverse */
  0,                                            /*tp_clear */
  (richcmpfunc)poly_richCompare,                /*tp_richcompare */
  0,                                            /*tp_weaklistoffset */
  (getiterfunc)poly_getiter,                    /*tp_iter */
  0,                                            /*tp_iternext */
  0,                                            /*tp_methods */
  0,                                            /*tp_members */
  0,                                            /*tp_getset */
  0,                                            /*tp_base */
  0,                                            /*tp_dict */
  0,                                            /*tp_descr_get */
  0,                                            /*tp_descr_set */
  0,                                            /*tp_dictoffset */
  (initproc)poly_init,                          /*tp_init */
  0,                                            /*tp_alloc */
  (newfunc)poly_new,                            /*tp_new */
};

void polyIterator_dealloc(PolyIterator *self) {
  self->ob_type->tp_free((PyObject*)self);
}

PyObject *polyIterator_getiter(Poly *self) {
  PolyIterator* iter = PyObject_New(PolyIterator, &polyIterator_type);
  iter->iter = ((const IPoly*)(self->poly))->begin();
  return (PyObject*)iter;
}

PyObject* polyIterator_next(PolyIterator *self) {
  PyObject* r = NULL;
  if (*self->iter) {
    Monom* m = PyObject_New(Monom, &monom_type);
    m->monom = (*self->iter).monom().monomInterface()->copy((*self->iter).monom());
    Coeff* c = PyObject_New(Coeff, &coeff_type);
    c->coeff = (*self->iter).coeff().coeffInterface()->copy((*self->iter).coeff());
    ++(*self->iter);
    r = PyTuple_New(2);
    PyTuple_SetItem(r, 0, (PyObject*)m);
    PyTuple_SetItem(r, 1, (PyObject*)c);
  }
  return r;
}

PyTypeObject polyIterator_type = {
  PyObject_HEAD_INIT(NULL)
  0,                                            /*ob_size*/
  "PolyIterator",                               /*tp_name*/
  sizeof(PolyIterator),                         /*tp_size*/
  0,                                            /*tp_itemsize*/
  (destructor)polyIterator_dealloc,             /*tp_dealloc*/
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
  "PolyIterator object",                        /*tp_doc */
  0,                                            /*tp_traverse */
  0,                                            /*tp_clear */
  0,                                            /*tp_richcompare */
  0,                                            /*tp_weaklistoffset */
  (getiterfunc)polyIterator_getiter,            /*tp_iter */
  (iternextfunc)polyIterator_next,              /*tp_iternext */
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
