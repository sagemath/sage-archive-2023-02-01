//
//  farey.cpp
//
//  Implementation of FareySymbol
//
//****************************************************************************
//       Copyright (C) 2011 Hartmut Monien <monien@th.physik.uni-bonn.de>
//                          Stefan Krämer <skraemer@th.physik.uni-bonn.de>
//
//  Distributed under the terms of the GNU General Public License (GPL)
//
//    This code is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    General Public License for more details.
//
//  The full text of the GPL is available at:
//
//                  http://www.gnu.org/licenses/
//****************************************************************************

#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <gmpxx.h>
#include <Python.h>

#include "farey.hpp"
#include "farey_symbol.h"

using namespace std;

template <class T>
inline
ostream& operator<<(ostream& os, const vector<T>& v) {
  os << v.size() << " ";
  for(typename vector<T>::const_iterator i=v.begin(); i!=v.end(); i++) {
    os << *i << " ";
  }
  return os;
}

template <class T>
inline
istream& operator>>(istream& is, vector<T>& v) {
  size_t n;
  is >> n;
  for(size_t i=0; i<n; i++) {
    T tmp;
    is >> tmp;
    v.push_back(tmp);
  }
  return is;
}

template <>
inline
istream& operator>>(istream& is, vector<SL2Z>& v) {
  size_t n;
  is >> n;
  for(size_t i=0; i<n; i++) {
    SL2Z tmp(1, 0, 0, 1);
    is >> tmp;
    v.push_back(tmp);
  }
  return is;
}

inline
istream& operator>>(istream& is, FareySymbol& F) {
  is >> F.pairing_max
     >> F.pairing
     >> F.cusp_classes
     >> F.a
     >> F.b
     >> F.x
     >> F.coset
     >> F.generators
     >> F.cusps
     >> F.cusp_widths;
  return is;
}

inline
ostream& operator<<(ostream& os, const FareySymbol& F) {
  os << F.pairing_max << " "
     << F.pairing
     << F.cusp_classes
     << F.a
     << F.b
     << F.x
     << F.coset
     << F.generators
     << F.cusps
     << F.cusp_widths;
  return os;
}

inline
mpq_class operator/(const mpz_class& a, const mpz_class& b) {
  return mpq_class(a, b);
}

inline
mpz_class lcm(const mpz_class& a, const mpz_class& b) {
  mpz_class result;
  mpz_lcm(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return result;
}

inline
mpz_class lcm(const vector<mpz_class>& v) {
  mpz_class q(1);
  for(size_t i=0; i<v.size(); i++) {
    q = lcm(q, v[i]);
  }
  return q;
}

//--- Helper class for checking membership of SL2Z matrix in group GammaH ----

is_element_GammaH::is_element_GammaH(int p_, PyObject* gen_list) : p(p_) {
  typedef vector<long>::const_iterator const_iterator;
  // list of generators
  vector<long> gen;
  Py_ssize_t ngen = PyList_Size(gen_list);
  for(Py_ssize_t i=0; i<ngen; i++) {
    PyObject* item = PyList_GetItem(gen_list, i);
    gen.push_back(convert_to_long(item));
  }
  // generate H from generators
  H = gen;
  for(;;) {
    vector<long> m;
    for(const_iterator i=gen.begin(); i!=gen.end(); i++) {
      for(const_iterator j=H.begin(); j!=H.end(); j++) {
	long q = ((*i)*(*j))%p;
	if( find(H.begin(), H.end(), q) == H.end() and
	    find(m.begin(), m.end(), q) == m.end() ) {
	  m.push_back(q);
	}
      }
    }
    if( m.size() == 0 ) break;
    for(const_iterator i=m.begin(); i!=m.end(); i++) H.push_back(*i);
  }
  // sort for binary searches
  sort(H.begin(), H.end());
}

is_element_GammaH::~is_element_GammaH() {
}

bool is_element_GammaH::is_member(const SL2Z& m) const {
  mpz_class a = m.a()%p; if(a < 0) a+=p;
  mpz_class d = m.d()%p; if(d < 0) d+=p;
  if( m.c()%p != 0 ) return false;
  if( not binary_search(H.begin(), H.end(), a.get_si()) ) return false;
  if( not binary_search(H.begin(), H.end(), d.get_si()) ) return false;
  return true;
}

//--- Helper class for checking membership of SL2Z matrix in group -------------
// group defined by the python object.

is_element_general::is_element_general(PyObject* group_) : group(group_) {
  if( PyObject_HasAttrString(group, "__contains__") ) {
    method = PyObject_GetAttrString(group, "__contains__");
  } else {
    cerr << "group has to define __contains__" << endl;
    throw string(__FUNCTION__) + ": error.";
  }
}

is_element_general::~is_element_general() {
  Py_DECREF(method);
}

bool is_element_general::is_member(const SL2Z& m) const {
  PyObject* arg = convert_to_SL2Z(m);
  PyObject* tuple = PyTuple_New(1);
  PyTuple_SetItem(tuple, 0, arg);
  PyObject *result = PyEval_CallObject(method, tuple);
  Py_DECREF(tuple);
  if( not PyBool_Check(result) ) {
    cerr << "__contains__ does not return bool." << endl;
    throw string(__FUNCTION__) + ": error.";
  }
  bool value = (result == Py_True);
  Py_DECREF(result);
  return value;
}

//--- FareySymbol ------------------------------------------------------------

// SL2Z

FareySymbol::FareySymbol() {
  pairing = vector<int>(2);
  pairing[0] = EVEN;
  pairing[1] = ODD;
  pairing_max = NO;
  a.push_back(0);
  b.push_back(1);
  cusp_widths.push_back(1);
  coset.push_back(SL2Z::E);
  generators.push_back(SL2Z::S);
  generators.push_back(SL2Z::S*SL2Z::R);
  cusp_classes.push_back(1);
  for(size_t i=0; i<a.size(); i++) x.push_back(a[i]/b[i]);
}

// User defined group and restoring from pickle

FareySymbol::FareySymbol(PyObject* o) {
  if( PyString_Check(o) ) {
    // restoration from data
    istringstream is(PyString_AsString(o));
    is >> (*this);
  } else {
    // init with user defined group
    is_element_general *group = new is_element_general(o);
    // check for user defined SL2Z
    if( group->is_member(SL2Z::S) and
	group->is_member(SL2Z::T) ) {
      pairing = vector<int>(2);
      pairing[0] = EVEN;
      pairing[1] = ODD;
      pairing_max = NO;
      a.push_back(0);
      b.push_back(1);
      cusp_widths.push_back(1);
      coset.push_back(SL2Z::E);
      generators.push_back(SL2Z::S);
      generators.push_back(SL2Z::S*SL2Z::R);
      cusp_classes.push_back(1);
      for(size_t i=0; i<a.size(); i++) x.push_back(a[i]/b[i]);
    }
    // check for index two subgroup
    else if( group->is_member(SL2Z( 0,  1, -1, -1)) and
	     group->is_member(SL2Z(-1,  1, -1,  0)) ) {
      pairing = vector<int>(2);
      pairing[0] = ODD;
      pairing[1] = ODD;
      a.push_back(0);
      b.push_back(1);
      for(size_t i=0; i<a.size(); i++) x.push_back(a[i]/b[i]);
      generators.push_back(SL2Z(1,  2,  0,  1));
      generators.push_back(SL2Z(0, -1,  1, -1));
      generators.push_back(SL2Z(1, -1,  1,  0));
      coset.push_back(SL2Z::E);
      coset.push_back(SL2Z::T);
      cusp_classes.push_back(1);
    } else {
      // everything else
      init_pairing(group);
      cusp_widths  = init_cusp_widths();
      coset        = init_coset_reps();
      generators   = init_generators(group);
      cusp_classes = init_cusp_classes();
      for(size_t i=0; i<a.size(); i++) x.push_back(a[i]/b[i]);
      cusps        = init_cusps();
    }
    delete group;
  }
}

// Predefined subgroups of SL2Z

FareySymbol::FareySymbol(PyObject* o, const is_element_group* group) {
  init_pairing(group);
  cusp_widths  = init_cusp_widths();
  coset        = init_coset_reps();
  generators   = init_generators(group);
  cusp_classes = init_cusp_classes();
  for(size_t i=0; i<a.size(); i++) x.push_back(a[i]/b[i]);
  cusps        = init_cusps();
}

FareySymbol::~FareySymbol() {
}

void FareySymbol::init_pairing(const is_element_group* group) {
  pairing = vector<int>(3, NO);
  const mpq_class infinity(10000000);
  pairing_max = NO;
  if( group->is_member(SL2Z(-1, 1, -1, 0)) ) {
    a.push_back(-1);
    a.push_back(0);
    b.push_back(1);
    b.push_back(1);
  } else {
    a.push_back(0);
    a.push_back(1);
    b.push_back(1);
    b.push_back(1);
  }
  check_pair(group, 0);
  check_pair(group, 1);
  check_pair(group, 2);
  for(;;) {
    int missing_pair(-1);
    mpq_class largest_diameter(0);
    for(size_t i=0; i<pairing.size(); i++) {
      if( pairing[i] == NO ) {
	if( i+1 != pairing.size() ) {
	  if( i != 0 ) {
	    mpq_class d = a[i]/b[i] - a[i-1]/b[i-1];
	    if( d > largest_diameter ) {
	      largest_diameter = d;
	      missing_pair = (int)(i);
	    }
	  } else {
	    largest_diameter = infinity;
	    missing_pair = 0;
	  }
	} else {
	  largest_diameter = infinity;
	  missing_pair = (int)(pairing.size()-1);
	  break;
	}
      }
    }
    if( missing_pair == -1 ) {
      break;
    } else {
      mpz_class A, B;
      if( missing_pair+1 == pairing.size() ) {
	A = a[missing_pair-1] + 1;
	B = b[missing_pair-1] + 0;
      } else {
	if( missing_pair == 0 ) {
	  A = a[0] - 1;
	  B = b[0] + 0;
	} else {
	  A = a[missing_pair-1]+a[missing_pair];
	  B = b[missing_pair-1]+b[missing_pair];
	}
      }
      add_term(missing_pair, A/B);
    }
    check_pair(group, missing_pair);
    check_pair(group, missing_pair+1);
  }
}

void FareySymbol::add_term(const int i, const mpq_class& r) {
  a.insert(a.begin()+i, r.get_num());
  b.insert(b.begin()+i, r.get_den());
  pairing.insert(pairing.begin()+i, NO);
}

void FareySymbol::check_pair(const is_element_group* group, const int i) {
  if( pairing[i] == NO ) {
    std::vector<int> even(pairing), odd(pairing);
    even[i] = EVEN;
    odd [i] = ODD;
    SL2Z A = pairing_matrix(even, i);
    SL2Z B = pairing_matrix(odd , i);
    if( group->is_member(A) or group->is_member(-A) ) {
      pairing[i] = EVEN;
      return;
    } else if( group->is_member(B) or group->is_member(-B) ) {
      pairing[i] = ODD;
      return;
    }
  }
  if( pairing[i] == NO ) {
    for(size_t j=0; j<pairing.size(); j++) {
      if( pairing[j] == NO and i != j ) {
	vector<int> p(pairing);
	p[i] = pairing_max+1;
	p[j] = pairing_max+1;
	SL2Z C = pairing_matrix(p, i);
	if( group->is_member(C) or group->is_member(-C) ) {
	  pairing_max++;
	  pairing[i] = pairing_max;
	  pairing[j] = pairing_max;
	  return;
	}
      }
    }
  }
}

SL2Z FareySymbol::pairing_matrix(const vector<int>& p, const size_t i) const {
  mpz_class ai, ai1, bi, bi1, aj, aj1, bj, bj1;
  if( i == 0 ) {
    ai = -1; bi = 0; ai1 = a[0], bi1 = b[0];
  } else if( i+1 == p.size() ) {
    ai = a[i-1]; bi = b[i-1]; ai1 = 1; bi1 = 0;
  } else {
    ai = a[i-1]; bi = b[i-1]; ai1 = a[i]; bi1 = b[i];
  }
  if( p[i] == NO ) {
    throw(string(__FUNCTION__)+string(": error"));
  } else if( p[i] == EVEN ) {
    return SL2Z(ai1*bi1+ai*bi, -ai*ai-ai1*ai1,
		bi*bi+bi1*bi1, -ai1*bi1-ai*bi);
  } else if( p[i] == ODD ) {
    return SL2Z(ai1*bi1+ai*bi1+ai*bi, -ai*ai-ai*ai1-ai1*ai1,
		bi*bi+bi*bi1+bi1*bi1, -ai1*bi1-ai1*bi-ai*bi);
  } else if( p[i] > NO ) {
    const size_t j = paired_side(p, i);
    if( j == 0 ) {
      aj = -1; bj = 0; aj1 = a[0]; bj1 = b[0];
    } else if( j == a.size() ) {
      aj = a[j-1]; bj = b[j-1]; aj1 = 1; bj1 = 0;
    } else {
      aj = a[j-1]; bj = b[j-1]; aj1 = a[j]; bj1 = b[j];
    }
    return SL2Z(aj1*bi1+aj*bi, -aj*ai-aj1*ai1,
		bj*bi+bj1*bi1, -ai1*bj1-ai*bj);
  }
  return SL2Z::E;
}

SL2Z FareySymbol::pairing_matrix(const size_t i) const {
  return pairing_matrix(pairing, i);
}

size_t FareySymbol::paired_side(const vector<int>& p, const size_t n) const {
  if( p[n] == EVEN or p[n] == ODD ) {
    return n;
  } else if( p[n] > NO ) {
    vector<int>::const_iterator i = find(p.begin(), p.end(), p[n]);
    if( i-p.begin() != n ) {
      return i-p.begin();
    } else {
      vector<int>::const_iterator j = find(i+1, p.end(), p[n]);
      return j-p.begin();
    }
  }
  throw(string(__FUNCTION__)+string(": error"));
  return 0;
}

vector<SL2Z> FareySymbol::init_generators(const is_element_group *group) const {
  const SL2Z I(-1, 0, 0, -1);
  vector<SL2Z> gen;
  vector<int> p;
  for(size_t i=0; i<pairing.size(); i++) {
    if( find(p.begin(), p.end(), pairing[i]) == p.end() ) {
      SL2Z m = pairing_matrix(i);
      if( not group->is_member(m) ) m = I*m;
      gen.push_back(m);
      if( pairing[i] > NO ) p.push_back(pairing[i]);
    }
  }
  if( nu2() == 0 and group->is_member(I) ) gen.push_back(I);
  return gen;
}

vector<mpq_class> FareySymbol::init_cusp_widths() const {
  static const mpq_class one_half(1, 2);
  vector<mpz_class> A(a), B(b);
  A.push_back(1);
  B.push_back(0);
  vector<mpq_class> w(A.size(), 0);
  for(size_t i=0; i<w.size(); i++) {
    size_t im = (i==0 ? A.size()-1 : i-1);
    size_t ip = (i+1==A.size() ? 0 : i+1);
    w[i] = abs(A[im]*B[ip]-A[ip]*B[im]);
    if( pairing[i ] == ODD ) w[i] += one_half;
    if( pairing[ip] == ODD ) w[i] += one_half;
  }
  return w;
}

vector<SL2Z> FareySymbol::init_coset_reps() const {
  static const mpq_class one_half(1, 2);
  vector<mpz_class> A(a), B(b);
  A.insert(A.begin(), -1);
  B.insert(B.begin(),  0);
  vector<mpq_class> cw(cusp_widths);
  rotate(cw.begin(), cw.end()-1, cw.end());
  vector<int> p(pairing);
  rotate(p.begin(), p.end()-1, p.end());
  vector<SL2Z> reps;
  for(size_t i=0; i<p.size(); i++) {
    size_t j = (i+1) % p.size();
    mpz_class c(A[j]), d(B[j]);
    if( d == 0 ) c = 1;
    mpq_class upper_bound(cw[i]);
    if( p[i] == ODD ) upper_bound += one_half;
    if( p[j] == ODD ) upper_bound -= one_half;
    for(size_t k=0; k<upper_bound; k++) {
      reps.push_back(SL2Z(1, -(int)(k), 0, 1)/SL2Z(-A[i], c, -B[i], d));
    }
  }
  return reps;
}

// init cusp classes is a class of the pairing alone !!!

vector<int> FareySymbol::init_cusp_classes() const {
  vector<int> c(pairing.size(), 0);
  int cusp_class(1);
  for(size_t m=0; m<c.size(); m++) {
    if( c[m] != 0 ) {
      continue;
    }
    c[m] = cusp_class;
    size_t i(m), I, J;
    for(;;) {
      if( pairing[i] == NO ) {
	I = i;
	J = (i==0? pairing.size()-1 : (i-1)%c.size());
      } else {
	I = (i+1)%c.size();
	J = I;
      }
      if( pairing[I] == ODD or pairing[I] == EVEN ) {
	if( c[I] == cusp_class ) {
	  cusp_class++;
	  break;
	}
	c[J] = cusp_class;
	i = J;
	continue;
      } else if( pairing[I] > NO ) {
	size_t j;
	for(size_t k=0; k<c.size(); k++) {
	  if( pairing[k] == pairing[I] and k != I ) j = k;
	}
	if( I != i ) {
	  if( c[j] == cusp_class ) {
	    cusp_class++;
	    break;
	  }
	  c[j] = cusp_class;
	  i = j;
	  continue;
	} else {
	  if( c[j-1] == cusp_class ) {
	    cusp_class++;
	    break;
	  }
	  c[j-1] = cusp_class;
	  i = j-1;
	  continue;
	}
      }
    }
  }
  return c;
}

vector<mpq_class> FareySymbol::init_cusps() const {
  // initialize cusps by identifying fractions using the cusp_class number
  vector<mpq_class> c;
  for(int i=1; i<=number_of_cusps(); i++) {
    for(size_t j=0; j<cusp_classes.size(); j++) {
      if( cusp_classes[j] == i and cusp_classes[j] != cusp_classes.back() ) {
	c.push_back(x[j]);
	break;
      }
    }
  }
  // width of cusp at infinity
  mpq_class W(0);
  for(size_t i=0; i<cusp_classes.size(); i++) {
    if( cusp_classes[i] == cusp_classes.back() ) W += cusp_widths[i];
  }
  // shift negative cusps to positve
  for(size_t i=0; i<c.size(); i++) while(c[i]<0) c[i] += W;
  sort(c.begin(), c.end());
  return c;
}

size_t FareySymbol::nu2() const {
  return count(pairing.begin(), pairing.end(), EVEN);
}

size_t FareySymbol::nu3() const {
  return count(pairing.begin(), pairing.end(), ODD);
}

size_t FareySymbol::rank_pi() const {
  if( index() == 2 ) return 1;
  return count_if(pairing.begin(), pairing.end(),
		  bind2nd(greater<int>(), 0))/2;
}

size_t FareySymbol::number_of_cusps() const {
  return size_t(*max_element(cusp_classes.begin(), cusp_classes.end()));
}

size_t FareySymbol::genus() const {
  return (rank_pi()-number_of_cusps()+1)/2;
}

size_t FareySymbol::level() const {
  if( index() == 2 ) return 2;
  vector<mpz_class> A(a), B(b);
  A.push_back(1);
  B.push_back(0);
  vector<mpz_class> width;
  for(size_t i=1; i<=number_of_cusps(); i++) {
    mpq_class cusp_width(0);
    for(size_t j=0; j<cusp_widths.size(); j++) {
      if( cusp_classes[j] == i ) {
	cusp_width += cusp_widths[j];
      }
    }
    width.push_back(cusp_width.get_num());
  }
  return lcm(width).get_ui();
}

//--- communication with sage ------------------------------------------------

size_t FareySymbol::index() const {
  return coset.size();
}

PyObject* FareySymbol::get_coset() const {
  PyObject* coset_list = PyList_New(coset.size());
  for(size_t i=0; i<coset.size(); i++) {
    PyObject* m = convert_to_SL2Z(coset[i]);
    PyList_SetItem(coset_list, i, m);
  }
  return coset_list;
}

PyObject* FareySymbol::get_generators() const {
  PyObject* generators_list = PyList_New(generators.size());
  for(size_t i=0; i<generators.size(); i++) {
    PyObject* m = convert_to_SL2Z(generators[i]);
    PyList_SetItem(generators_list, i, m);
  }
  return generators_list;
}

PyObject* FareySymbol::get_cusp_widths() const {
  PyObject* cusp_widths_list = PyList_New(cusp_widths.size());
  for(size_t i=0; i<cusp_widths.size(); i++) {
    PyObject* m = convert_to_rational(cusp_widths[i]);
    PyList_SetItem(cusp_widths_list, i, m);
  }
  return cusp_widths_list;
}

PyObject* FareySymbol::get_cusps() const {
  PyObject* cusps_list = PyList_New(cusps.size());
  for(size_t i=0; i<cusps.size(); i++) {
    PyObject* m = convert_to_cusp(cusps[i]);
    PyList_SetItem(cusps_list, i, m);
  }
  return cusps_list;
}

PyObject* FareySymbol::get_fractions() const {
  PyObject* x_list = PyList_New(x.size());
  for(size_t i=0; i<x.size(); i++) {
    PyObject* m = convert_to_rational(x[i]);
    PyList_SetItem(x_list, i, m);
  }
  return x_list;
}

PyObject* FareySymbol::get_pairings() const {
  PyObject* pairing_list = PyList_New(pairing.size());
  for(size_t i=0; i<pairing.size(); i++) {
    PyObject* m = PyInt_FromLong(long(pairing[i]));
    PyList_SetItem(pairing_list, i, m);
  }
  return pairing_list;
}

PyObject* FareySymbol::get_paired_sides() const {
  vector<int> p;
  for(size_t i=0; i<pairing.size(); i++) {
    if( pairing[i] > NO and
	p.end() == find(p.begin(), p.end(), pairing[i]) ) {
      p.push_back(pairing[i]);
    }
  }
  sort(p.begin(), p.end());
  PyObject* pairing_list = PyList_New(p.size());
  for(vector<int>::const_iterator i=p.begin(); i!=p.end(); i++) {
    vector<int>::const_iterator j = find(pairing.begin(), pairing.end(), *i);
    vector<int>::const_iterator k = find(j+1, pairing.end(), *i);
    PyObject* J = PyInt_FromLong(long(j-pairing.begin()));
    PyObject* K = PyInt_FromLong(long(k-pairing.begin()));
    PyObject* tuple = PyTuple_New(2);
    PyTuple_SetItem(tuple, 0, J);
    PyTuple_SetItem(tuple, 1, K);
    PyList_SetItem(pairing_list, i-p.begin(), tuple);
  }
  return pairing_list;
}

PyObject* FareySymbol::get_pairing_matrices() const {
  PyObject* pairing_matrix_list = PyList_New(pairing.size());
  for(size_t i=0; i<pairing.size(); i++) {
    PyObject* pm = convert_to_SL2Z(pairing_matrix(i));
    PyList_SetItem(pairing_matrix_list, i, pm);
  }
  return pairing_matrix_list;
}

PyObject* FareySymbol::dumps() const {
  ostringstream os(ostringstream::out|ostringstream::binary);
  os << (*this);
  PyObject* d = PyString_FromString(os.str().c_str());
  return d;
}

