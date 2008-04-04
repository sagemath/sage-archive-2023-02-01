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

#include "./monom.h"
#include "./poly.h"
#include "./wrap.h"
#include "./division.h"
#include "./basis.h"

#include "../ginv/poly/iexpression.h"
#include "../ginv/algorithm/ialgorithm.h"

#include <string.h>
#include <sstream>

void basis_dealloc(Basis *self) {
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* basis_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  Basis *self = (Basis*)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->algorithm = NULL;
  }

  return (PyObject *)self;
}

PyObject* basis_getattr(Basis *self, char *name) {
  PyObject *r = NULL;
  static PyMethodDef methods[] = {
    {"lengthIB", (PyCFunction)basis_lengthIB, METH_NOARGS , ""},
    {"lengthGB", (PyCFunction)basis_lengthGB, METH_NOARGS , ""},
    {"numCriterion", (PyCFunction)basis_numCriterion, METH_NOARGS , ""},
    {"numSpoly", (PyCFunction)basis_numSpoly, METH_NOARGS , ""},
    {"numReduction", (PyCFunction)basis_numReduction, METH_NOARGS , ""},
    {"userTime", (PyCFunction)basis_userTime, METH_NOARGS , ""},
    {"sysTime", (PyCFunction)basis_sysTime, METH_NOARGS , ""},
    {"realTime", (PyCFunction)basis_realTime, METH_NOARGS , ""},
    {"hilbertPolynomial", (PyCFunction)basis_hilbertPolynomial, METH_NOARGS , ""},
    {"iterStatistics", (PyCFunction)basis_iterStatistics, METH_NOARGS , ""},
    {"iterIB", (PyCFunction)basis_iterIB, METH_NOARGS , ""},
    {"iterGB", (PyCFunction)basis_iterGB, METH_NOARGS , ""},
    {"iterLmIB", (PyCFunction)basis_iterLmIB, METH_NOARGS , ""},
    {"iterLmGB", (PyCFunction)basis_iterLmGB, METH_NOARGS , ""},
    {"find", (PyCFunction)basis_find, METH_VARARGS , ""},
    {"nfHead", (PyCFunction)basis_nf, METH_VARARGS , ""},
    {"isNfHead", (PyCFunction)basis_isNf, METH_VARARGS , ""},
    {"nfTail", (PyCFunction)basis_nfTail, METH_VARARGS , ""},
    {"isNfTail", (PyCFunction)basis_isNfTail, METH_VARARGS , ""},
    {NULL}
  };
  r = Py_FindMethod(methods, (PyObject*)self, name);

  if (r == NULL)
    PyErr_SetString(PyExc_TypeError, "Cannot definition operation in 'Basis'");

  return r;
}

PyObject* basis_getiter(Basis *self) {
  BasisIterator* iter = PyObject_New(BasisIterator, &basisIterator_type);
  iter->iter = self->algorithm->setT().begin();
  return (PyObject*)iter;
}

PyObject* basis_lengthIB(Basis *self) {
  return PyInt_FromLong(self->algorithm->statistics().mLengthIB);
}

PyObject* basis_lengthGB(Basis *self) {
  return PyInt_FromLong(self->algorithm->statistics().mLengthGB);
}

PyObject* basis_numCriterion(Basis *self) {
  return PyInt_FromLong(self->algorithm->statistics().mCriterion);
}

PyObject* basis_numSpoly(Basis *self) {
  return PyInt_FromLong(self->algorithm->statistics().mSpoly);
}

PyObject* basis_numReduction(Basis *self) {
  return PyInt_FromLong(self->algorithm->statistics().mReduction);
}

PyObject* basis_userTime(Basis *self) {
  return PyFloat_FromDouble(self->algorithm->statistics().mTimer.userTime());
}

PyObject* basis_sysTime(Basis *self) {
  return PyFloat_FromDouble(self->algorithm->statistics().mTimer.sysTime());
}

PyObject* basis_realTime(Basis *self) {
  return PyFloat_FromDouble(self->algorithm->statistics().mTimer.realTime());
}

PyObject* basis_hilbertPolynomial(Basis *self) {
  std::ostringstream out;
  out << self->algorithm->statistics().mHilbertPoly;
  return PyString_FromString(out.str().c_str());
}

PyObject* basis_iterStatistics(Basis *self) {
  StatisticsIterator* iter = PyObject_New(StatisticsIterator, &statisticsIterator_type);
  iter->iter = self->algorithm->statistics().begin();
  return (PyObject*)iter;
}

PyObject* basis_iterIB(Basis *self) {
  IBIterator* iter = PyObject_New(IBIterator, &iBIterator_type);
  iter->iter = self->algorithm->setT().begin();
  return (PyObject*)iter;
}

PyObject* basis_iterGB(Basis *self) {
  GBIterator* iter = PyObject_New(GBIterator, &gBIterator_type);
  iter->iter = self->algorithm->setT().begin();
  while(iter->iter && (*iter->iter)->isProlong())
    ++iter->iter;
  return (PyObject*)iter;
}

PyObject* basis_iterLmIB(Basis *self) {
  LmIBIterator* iter = PyObject_New(LmIBIterator, &lMiBIterator_type);
  iter->iter = self->algorithm->setT().begin();
  return (PyObject*)iter;
}

PyObject* basis_iterLmGB(Basis *self) {
  LmGBIterator* iter = PyObject_New(LmGBIterator, &lMgBIterator_type);
  iter->iter = self->algorithm->setT().begin();
  while(iter->iter && (*iter->iter)->isProlong())
    ++iter->iter;
  return (PyObject*)iter;
}

PyObject* basis_find(Basis *self, PyObject *args) {
  PyObject *py_monom = NULL;
  if (!PyArg_ParseTuple(args, "O:Basis",
                        &py_monom))
      return NULL;

  IMonom *monom = NULL;
  IMonomInterface *monomInterface = self->algorithm->setT().wrapInterface()->monomInterface();
  if (PyObject_TypeCheck(py_monom, &PyString_Type)) {
    IMonom *m = monomInterface->create();
    std::istringstream in(PyString_AsString(py_monom));
    in >> *m >> std::ws;
    if (!in.eof()) {
      PyErr_SetString(PyExc_TypeError, "error parsing string in Basis.basis_find()");
      delete m;
    }
    else
      monom = m;
  }
  else if (PyObject_TypeCheck(py_monom, &monom_type)) {
    IMonom *m = ((Monom*)py_monom)->monom;
    if (monomInterface != m->monomInterface())
      PyErr_SetString(PyExc_TypeError, "in Basis.basis_find() Monomiall with difference MonomInterface");
    else
      monom = m;
  }
  else
    PyErr_SetString(PyExc_TypeError, "nf element type string or 'Monom' in 'Basis.basis_find()'");

  PyObject *r = NULL;
  if (monom != NULL) {
    if (self->algorithm->division() == NULL) {
    }
    else {
      Wrap* res = PyObject_New(Wrap, &wrap_type);
      res->wrap = self->algorithm->division()->find(*monom);
      r = (PyObject*)res;
      if (PyObject_TypeCheck(py_monom, &PyString_Type))
        delete monom;
    }
  }

  return r;
}

PyObject* basis_nf(Basis *self, PyObject *args) {
  PyObject *py_poly = NULL;
  if (!PyArg_ParseTuple(args, "O:Basis",
                        &py_poly))
      return NULL;

  IPoly *poly = NULL;
  IPolyInterface *polyInterface = self->algorithm->setT().wrapInterface()->polyInterface();
  if (PyObject_TypeCheck(py_poly, &PyString_Type)) {
    std::istringstream in(PyString_AsString(py_poly));
    IExpression expr(polyInterface);
    in >> expr >> std::ws;
    if (!in.eof())
      PyErr_SetString(PyExc_TypeError, "error parsing string in Basis.nf()");
    else
      poly = expr.toPoly();
  }
  else if (PyObject_TypeCheck(py_poly, &poly_type)) {
    IPoly *p = ((Poly*)py_poly)->poly;
    if (polyInterface != p->polyInterface())
      PyErr_SetString(PyExc_TypeError, "in Basis.nf() Polynomial with difference PolyInterface");
    else
      poly = p->polyInterface()->copy(*p);
  }
  else
    PyErr_SetString(PyExc_TypeError, "nf element type string or 'Poly' in 'Basis.nf()'");

  PyObject *r = NULL;
  if (poly != NULL) {
    if (self->algorithm->division() == NULL) {

    }
    else {
      self->algorithm->division()->nf(*poly);
      if (!poly->isZero())
        self->algorithm->division()->nfTail(*poly);
      Poly* res = PyObject_New(Poly, &poly_type);
      res->poly = poly;
      r = (PyObject*)res;
    }
  }

  return r;
}

PyObject* basis_isNf(Basis *self, PyObject *args) {
  PyObject *py_poly = NULL;
  if (!PyArg_ParseTuple(args, "O:Basis",
                        &py_poly))
      return NULL;

  IPoly *poly = NULL;
  IPolyInterface *polyInterface = self->algorithm->setT().wrapInterface()->polyInterface();
  if (PyObject_TypeCheck(py_poly, &PyString_Type)) {
    std::istringstream in(PyString_AsString(py_poly));
    IExpression expr(polyInterface);
    in >> expr >> std::ws;
    if (!in.eof())
      PyErr_SetString(PyExc_TypeError, "error parsing string in Basis.nf()");
    else
      poly = expr.toPoly();
  }
  else if (PyObject_TypeCheck(py_poly, &poly_type)) {
    IPoly *p = ((Poly*)py_poly)->poly;
    if (polyInterface != p->polyInterface())
      PyErr_SetString(PyExc_TypeError, "in Basis.nf() Polynomial with difference PolyInterface");
    else
      poly = p;
  }
  else
    PyErr_SetString(PyExc_TypeError, "nf element type string or 'Poly' in 'Basis.nf()'");

  PyObject *r = NULL;
  if (poly != NULL) {
    if (self->algorithm->division() == NULL) {
    }
    else {
      if (!poly->isZero() &&
          self->algorithm->division()->isNf(*poly) &&
          self->algorithm->division()->isNfTail(*poly))
        r = Py_True;
      else
        r = Py_False;
      Py_INCREF(r);
      if (PyObject_TypeCheck(py_poly, &PyString_Type))
        delete poly;
    }
  }

  return r;
}

PyObject* basis_nfTail(Basis *self, PyObject *args) {
  PyObject *py_poly = NULL;
  if (!PyArg_ParseTuple(args, "O:Basis",
                        &py_poly))
      return NULL;

  IPoly *poly = NULL;
  IPolyInterface *polyInterface = self->algorithm->setT().wrapInterface()->polyInterface();
  if (PyObject_TypeCheck(py_poly, &PyString_Type)) {
    std::istringstream in(PyString_AsString(py_poly));
    IExpression expr(polyInterface);
    in >> expr >> std::ws;
    if (!in.eof())
      PyErr_SetString(PyExc_TypeError, "error parsing string in Basis.nfTail()");
    else
      poly = expr.toPoly();
  }
  else if (PyObject_TypeCheck(py_poly, &poly_type)) {
    IPoly *p = ((Poly*)py_poly)->poly;
    if (polyInterface != p->polyInterface())
      PyErr_SetString(PyExc_TypeError, "in Basis.nfTail() Polynomial with difference PolyInterface");
    else
      poly = p->polyInterface()->copy(*p);
  }
  else
    PyErr_SetString(PyExc_TypeError, "nfTail element type string or 'Poly' in 'Basis.nfTail()'");

  PyObject *r = NULL;
  if (poly != NULL) {
    if (self->algorithm->division() == NULL) {
    }
    else {
      if (!poly->isZero())
        self->algorithm->division()->nfTail(*poly);
      Poly* res = PyObject_New(Poly, &poly_type);
      res->poly = poly;
      r = (PyObject*)res;
    }
  }

  return r;
}

PyObject* basis_isNfTail(Basis *self, PyObject *args) {
  PyObject *py_poly = NULL;
  if (!PyArg_ParseTuple(args, "O:Basis",
                        &py_poly))
      return NULL;

  IPoly *poly = NULL;
  IPolyInterface *polyInterface = self->algorithm->setT().wrapInterface()->polyInterface();
  if (PyObject_TypeCheck(py_poly, &PyString_Type)) {
    std::istringstream in(PyString_AsString(py_poly));
    IExpression expr(polyInterface);
    in >> expr >> std::ws;
    if (!in.eof())
      PyErr_SetString(PyExc_TypeError, "error parsing string in Basis.nf()");
    else
      poly = expr.toPoly();
  }
  else if (PyObject_TypeCheck(py_poly, &poly_type)) {
    IPoly *p = ((Poly*)py_poly)->poly;
    if (polyInterface != p->polyInterface())
      PyErr_SetString(PyExc_TypeError, "in Basis.nf() Polynomial with difference PolyInterface");
    else
      poly = p;
  }
  else
    PyErr_SetString(PyExc_TypeError, "nf element type string or 'Poly' in 'Basis.nf()'");

  PyObject *r = NULL;
  if (poly != NULL) {
    if (self->algorithm->division() == NULL) {
    }
    else {
      if (self->algorithm->division()->isNfTail(*poly))
        r = Py_True;
      else
        r = Py_False;
      Py_INCREF(r);
      if (PyObject_TypeCheck(py_poly, &PyString_Type))
        delete poly;
    }
  }

  return r;
}

#if PY_VERSION_HEX < 0x02050000
int basis_length(Basis *self) {
#else
Py_ssize_t basis_length(Basis *self) {
#endif
  return self->algorithm->setT().length();
}

#if PY_VERSION_HEX < 0x02050000
PyObject* basis_item(Basis *self, int i) {
#else
PyObject* basis_item(Basis *self, Py_ssize_t i) {
#endif
  PyObject* r = NULL;
  if (0 > i || i >= self->algorithm->setT().length())
    PyErr_SetString(PyExc_ValueError, "Index out of range.");
  else {
    ISetT::ConstIterator k = self->algorithm->setT().begin();
    for(; i > 0; i--)
     ++k;

    Wrap* wrap = PyObject_New(Wrap, &wrap_type);
    wrap->wrap = *k;
    r = (PyObject*)wrap;
  }
  return r;
}

static PySequenceMethods basis_sequence = {
#if PY_VERSION_HEX < 0x02050000
  (inquiry)basis_length,		/*sq_length*/
#else
  (lenfunc)basis_length,		/*sq_length*/
#endif
  0,					/*sq_concat*/
  0,					/*sq_repeat*/
#if PY_VERSION_HEX < 0x02050000
  (intargfunc)basis_item,		/*sq_item*/
#else
  (ssizeargfunc)basis_item,		/*sq_item*/
#endif
  0,					/*sq_slice*/
  0,					/*sq_ass_item*/
  0,					/*sq_ass_slice*/
  0,					/*sq_contains */
  0,					/*sq_inplace_concat*/
  0,					/*sq_inplace_repeat*/
};

PyTypeObject basis_type = {
  PyObject_HEAD_INIT(NULL)
  0,						/*ob_size*/
  "Basis",					/*tp_name*/
  sizeof(Basis),				/*tp_size*/
  0,						/*tp_itemsize*/
  (destructor)basis_dealloc,			/*tp_dealloc*/
  0,						/*tp_print*/
  (getattrfunc)basis_getattr,			/*tp_getattr*/
  0,						/*tp_setattr*/
  0,						/*tp_compare*/
  0,						/*tp_repr*/
  0,						/*tp_as_number*/
  &basis_sequence,				/*tp_as_sequence*/
  0,						/*tp_as_mapping*/
  0,						/*tp_hash*/
  0,						/*tp_call*/
  0,						/*tp_str*/
  0,						/*tp_getattro*/
  0,						/*tp_setattro*/
  0,						/*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/*tp_flags*/
  "Basis object",				/*tp_doc */
  0,						/*tp_traverse */
  0,						/*tp_clear */
  0,						/*tp_richcompare */
  0,						/*tp_weaklistoffset */
  (getiterfunc)basis_getiter,			/*tp_iter */
  0,						/*tp_iternext */
  0,						/*tp_methods */
  0,						/*tp_members */
  0,						/*tp_getset */
  0,						/*tp_base */
  0,						/*tp_dict */
  0,						/*tp_descr_get */
  0,						/*tp_descr_set */
  0,						/*tp_dictoffset */
  0,						/*tp_init */
  0,						/*tp_alloc */
  (newfunc)basis_new,				/*tp_new */
};

void basisIterator_dealloc(BasisIterator *self) {
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* basisIterator_next(BasisIterator *self) {
  PyObject* r = NULL;
  if (self->iter) {
    Wrap* wrap = PyObject_New(Wrap, &wrap_type);
    wrap->wrap = *self->iter;
    ++self->iter;
    r = (PyObject*)wrap;
  }
  return r;
}

PyTypeObject basisIterator_type = {
  PyObject_HEAD_INIT(NULL)
  0,						/*ob_size*/
  "BasisIterator",				/*tp_name*/
  sizeof(BasisIterator),			/*tp_size*/
  0,						/*tp_itemsize*/
  (destructor)basisIterator_dealloc,		/*tp_dealloc*/
  0,						/*tp_print*/
  0,						/*tp_getattr*/
  0,						/*tp_setattr*/
  0,						/*tp_compare*/
  0,						/*tp_repr*/
  0,						/*tp_as_number*/
  0,						/*tp_as_sequence*/
  0,						/*tp_as_mapping*/
  0,						/*tp_hash*/
  0,						/*tp_call*/
  0,						/*tp_str*/
  0,						/*tp_getattro*/
  0,						/*tp_setattro*/
  0,						/*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/*tp_flags*/
  "BasisIterator object",			/*tp_doc */
  0,						/*tp_traverse */
  0,						/*tp_clear */
  0,						/*tp_richcompare */
  0,						/*tp_weaklistoffset */
  0,						/*tp_iter */
  (iternextfunc)basisIterator_next,		/*tp_iternext */
  0,						/*tp_methods */
  0,						/*tp_members */
  0,						/*tp_getset */
  0,						/*tp_base */
  0,						/*tp_dict */
  0,						/*tp_descr_get */
  0,						/*tp_descr_set */
  0,						/*tp_dictoffset */
  0,						/*tp_init */
  0,						/*tp_alloc */
  0,						/*tp_new */
};

void statisticsIterator_dealloc(StatisticsIterator *self) {
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* statisticsIterator_getiter(StatisticsIterator *self) {
  StatisticsIterator* iter = PyObject_New(StatisticsIterator, &statisticsIterator_type);
  iter->iter = self->iter;
  return (PyObject*)iter;
}

PyObject* statisticsIterator_next(StatisticsIterator *self) {
  PyObject* r = NULL;
  if (self->iter) {
    r = PyTuple_New(4);
    PyTuple_SetItem(r, 0, PyFloat_FromDouble((*self->iter).mUserTime));
    PyTuple_SetItem(r, 1, PyInt_FromLong((*self->iter).mTlength));
    PyTuple_SetItem(r, 2, PyInt_FromLong((*self->iter).mQlength));
    PyTuple_SetItem(r, 3, PyInt_FromLong((*self->iter).mEqLmLength));
    ++self->iter;
  }
  return r;
}

PyTypeObject statisticsIterator_type = {
  PyObject_HEAD_INIT(NULL)
  0,						/*ob_size*/
  "StatisticsIterator",				/*tp_name*/
  sizeof(StatisticsIterator),			/*tp_size*/
  0,						/*tp_itemsize*/
  (destructor)statisticsIterator_dealloc,	/*tp_dealloc*/
  0,						/*tp_print*/
  0,						/*tp_getattr*/
  0,						/*tp_setattr*/
  0,						/*tp_compare*/
  0,						/*tp_repr*/
  0,						/*tp_as_number*/
  0,						/*tp_as_sequence*/
  0,						/*tp_as_mapping*/
  0,						/*tp_hash*/
  0,						/*tp_call*/
  0,						/*tp_str*/
  0,						/*tp_getattro*/
  0,						/*tp_setattro*/
  0,						/*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/*tp_flags*/
  "StatisticsIterator object",			/*tp_doc */
  0,						/*tp_traverse */
  0,						/*tp_clear */
  0,						/*tp_richcompare */
  0,						/*tp_weaklistoffset */
  (getiterfunc)statisticsIterator_getiter,	/*tp_iter */
  (iternextfunc)statisticsIterator_next,	/*tp_iternext */
  0,						/*tp_methods */
  0,						/*tp_members */
  0,						/*tp_getset */
  0,						/*tp_base */
  0,						/*tp_dict */
  0,						/*tp_descr_get */
  0,						/*tp_descr_set */
  0,						/*tp_dictoffset */
  0,						/*tp_init */
  0,						/*tp_alloc */
  0,						/*tp_new */
};

void iBIterator_dealloc(IBIterator *self) {
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* iBIterator_getiter(IBIterator *self) {
  IBIterator* iter = PyObject_New(IBIterator, &iBIterator_type);
  iter->iter = self->iter;
  return (PyObject*)iter;
}

PyObject* iBIterator_next(IBIterator *self) {
  PyObject* r = NULL;
  if (self->iter) {
    Poly* poly = PyObject_New(Poly, &poly_type);
    poly->poly = (*self->iter)->poly().polyInterface()->copy((*self->iter)->poly());
    ++self->iter;
    r = (PyObject*)poly;
  }
  return r;
}

PyTypeObject iBIterator_type = {
  PyObject_HEAD_INIT(NULL)
  0,						/*ob_size*/
  "IBIterator",					/*tp_name*/
  sizeof(IBIterator),				/*tp_size*/
  0,						/*tp_itemsize*/
  (destructor)iBIterator_dealloc,		/*tp_dealloc*/
  0,						/*tp_print*/
  0,						/*tp_getattr*/
  0,						/*tp_setattr*/
  0,						/*tp_compare*/
  0,						/*tp_repr*/
  0,						/*tp_as_number*/
  0,						/*tp_as_sequence*/
  0,						/*tp_as_mapping*/
  0,						/*tp_hash*/
  0,						/*tp_call*/
  0,						/*tp_str*/
  0,						/*tp_getattro*/
  0,						/*tp_setattro*/
  0,						/*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/*tp_flags*/
  "IBIterator object",				/*tp_doc */
  0,						/*tp_traverse */
  0,						/*tp_clear */
  0,						/*tp_richcompare */
  0,						/*tp_weaklistoffset */
  (getiterfunc)iBIterator_getiter,		/*tp_iter */
  (iternextfunc)iBIterator_next,		/*tp_iternext */
  0,						/*tp_methods */
  0,						/*tp_members */
  0,						/*tp_getset */
  0,						/*tp_base */
  0,						/*tp_dict */
  0,						/*tp_descr_get */
  0,						/*tp_descr_set */
  0,						/*tp_dictoffset */
  0,						/*tp_init */
  0,						/*tp_alloc */
  0,						/*tp_new */
};

void gBIterator_dealloc(GBIterator *self) {
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* gBIterator_getiter(GBIterator *self) {
  GBIterator* iter = PyObject_New(GBIterator, &gBIterator_type);
  iter->iter = self->iter;
  return (PyObject*)iter;
}

PyObject* gBIterator_next(GBIterator *self) {
  PyObject* r = NULL;
  if (self->iter) {
    Poly* poly = PyObject_New(Poly, &poly_type);
    poly->poly = (*self->iter)->poly().polyInterface()->copy((*self->iter)->poly());
    do {
      ++self->iter;
    } while (self->iter && (*self->iter)->isProlong());
    r = (PyObject*)poly;
  }
  return r;
}

PyTypeObject gBIterator_type = {
  PyObject_HEAD_INIT(NULL)
  0,						/*ob_size*/
  "GBIterator",					/*tp_name*/
  sizeof(GBIterator),				/*tp_size*/
  0,						/*tp_itemsize*/
  (destructor)gBIterator_dealloc,		/*tp_dealloc*/
  0,						/*tp_print*/
  0,						/*tp_getattr*/
  0,						/*tp_setattr*/
  0,						/*tp_compare*/
  0,						/*tp_repr*/
  0,						/*tp_as_number*/
  0,						/*tp_as_sequence*/
  0,						/*tp_as_mapping*/
  0,						/*tp_hash*/
  0,						/*tp_call*/
  0,						/*tp_str*/
  0,						/*tp_getattro*/
  0,						/*tp_setattro*/
  0,						/*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/*tp_flags*/
  "GBIterator object",				/*tp_doc */
  0,						/*tp_traverse */
  0,						/*tp_clear */
  0,						/*tp_richcompare */
  0,						/*tp_weaklistoffset */
  (getiterfunc)gBIterator_getiter,		/*tp_iter */
  (iternextfunc)gBIterator_next,		/*tp_iternext */
  0,						/*tp_methods */
  0,						/*tp_members */
  0,						/*tp_getset */
  0,						/*tp_base */
  0,						/*tp_dict */
  0,						/*tp_descr_get */
  0,						/*tp_descr_set */
  0,						/*tp_dictoffset */
  0,						/*tp_init */
  0,						/*tp_alloc */
  0,						/*tp_new */
};

void lMiBIterator_dealloc(LmIBIterator *self) {
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* lMiBIterator_getiter(LmIBIterator *self) {
  LmIBIterator* iter = PyObject_New(LmIBIterator, &lMiBIterator_type);
  iter->iter = self->iter;
  return (PyObject*)iter;
}

PyObject* lMiBIterator_next(LmIBIterator *self) {
  PyObject* r = NULL;
  if (self->iter) {
    Monom* monom = PyObject_New(Monom, &monom_type);
    monom->monom = (*self->iter)->lm().monomInterface()->copy((*self->iter)->lm());
    ++self->iter;
    r = (PyObject*)monom;
  }
  return r;
}

PyTypeObject lMiBIterator_type = {
  PyObject_HEAD_INIT(NULL)
  0,						/*ob_size*/
  "LmIBIterator",				/*tp_name*/
  sizeof(LmIBIterator),				/*tp_size*/
  0,						/*tp_itemsize*/
  (destructor)lMiBIterator_dealloc,		/*tp_dealloc*/
  0,						/*tp_print*/
  0,						/*tp_getattr*/
  0,						/*tp_setattr*/
  0,						/*tp_compare*/
  0,						/*tp_repr*/
  0,						/*tp_as_number*/
  0,						/*tp_as_sequence*/
  0,						/*tp_as_mapping*/
  0,						/*tp_hash*/
  0,						/*tp_call*/
  0,						/*tp_str*/
  0,						/*tp_getattro*/
  0,						/*tp_setattro*/
  0,						/*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/*tp_flags*/
  "LmIBIterator object",			/*tp_doc */
  0,						/*tp_traverse */
  0,						/*tp_clear */
  0,						/*tp_richcompare */
  0,						/*tp_weaklistoffset */
  (getiterfunc)lMiBIterator_getiter,		/*tp_iter */
  (iternextfunc)lMiBIterator_next,		/*tp_iternext */
  0,						/*tp_methods */
  0,						/*tp_members */
  0,						/*tp_getset */
  0,						/*tp_base */
  0,						/*tp_dict */
  0,						/*tp_descr_get */
  0,						/*tp_descr_set */
  0,						/*tp_dictoffset */
  0,						/*tp_init */
  0,						/*tp_alloc */
  0,						/*tp_new */
};

void lMgBIterator_dealloc(LmGBIterator *self) {
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* lMgBIterator_getiter(LmGBIterator *self) {
  LmGBIterator* iter = PyObject_New(LmGBIterator, &lMgBIterator_type);
  iter->iter = self->iter;
  return (PyObject*)iter;
}

PyObject* LmgBIterator_next(LmGBIterator *self) {
  PyObject* r = NULL;
  if (self->iter) {
    Monom* monom = PyObject_New(Monom, &monom_type);
    monom->monom = (*self->iter)->lm().monomInterface()->copy((*self->iter)->lm());
    do {
      ++self->iter;
    } while(self->iter && (*self->iter)->isProlong());
    r = (PyObject*)monom;
  }
  return r;
}

PyTypeObject lMgBIterator_type = {
  PyObject_HEAD_INIT(NULL)
  0,						/*ob_size*/
  "LmGBIterator",				/*tp_name*/
  sizeof(LmGBIterator),				/*tp_size*/
  0,						/*tp_itemsize*/
  (destructor)lMgBIterator_dealloc,		/*tp_dealloc*/
  0,						/*tp_print*/
  0,						/*tp_getattr*/
  0,						/*tp_setattr*/
  0,						/*tp_compare*/
  0,						/*tp_repr*/
  0,						/*tp_as_number*/
  0,						/*tp_as_sequence*/
  0,						/*tp_as_mapping*/
  0,						/*tp_hash*/
  0,						/*tp_call*/
  0,						/*tp_str*/
  0,						/*tp_getattro*/
  0,						/*tp_setattro*/
  0,						/*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,	/*tp_flags*/
  "LmGBIterator object",			/*tp_doc */
  0,						/*tp_traverse */
  0,						/*tp_clear */
  0,						/*tp_richcompare */
  0,						/*tp_weaklistoffset */
  (getiterfunc)lMgBIterator_getiter,		/*tp_iter */
  (iternextfunc)LmgBIterator_next,		/*tp_iternext */
  0,						/*tp_methods */
  0,						/*tp_members */
  0,						/*tp_getset */
  0,						/*tp_base */
  0,						/*tp_dict */
  0,						/*tp_descr_get */
  0,						/*tp_descr_set */
  0,						/*tp_dictoffset */
  0,						/*tp_init */
  0,						/*tp_alloc */
  0,						/*tp_new */
};
