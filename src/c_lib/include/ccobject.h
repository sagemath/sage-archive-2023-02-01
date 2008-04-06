/******************************************************************************
       Copyright (C) 2007 Joel B. Mohler <joel@kiwistrawberry.us>

  Distributed under the terms of the GNU General Public License (GPL), Version 2.

  The full text of the GPL is available at:
                  http://www.gnu.org/licenses/

******************************************************************************/

/**
 * @file ccobject.h
 *
 * @author Joel B. Mohler <joel@kiwistrawberry.us>
 *
 * @brief These functions provide assistance for constructing and destructing
 *  C++ objects from pyrex.
 *
 */

#ifdef __cplusplus

#include <iostream>
#include <sstream>

/* Ok, we are in C++ mode.  Lets define some templated functions for construction
and destruction of allocated memory chunks. */

/* Allocate and Construct */
template <class T>
T* New(){
  return new T();
}

/* Construct */
template <class T>
T* Construct(void* mem){
  return new(mem) T();
}

/* Construct with one parameter */
template <class T, class P>
T* Construct_p(void* mem, const P &d){
  return new(mem) T(d);
}

/* Construct with two parameters */
template <class T, class P, class Q>
T* Construct_pp(void* mem, const P &d, const Q &e){
  return new(mem) T(d, e);
}

/* Construct with three parameters */
template <class T, class P, class Q, class R>
T* Construct_ppp(void* mem, const P &d, const Q &e, const R &f){
  return new(mem) T(d, e, f);
}

/* Construct with four parameters */
template <class T, class P, class Q, class R, class S>
T* Construct_pppp(void* mem, const P &d, const Q &e, const R &f, const S &g){
  return new(mem) T(d, e, f, g);
}

/* Destruct */
template <class T>
void Destruct(T* mem){
  mem->~T();
}

/* Deallocate and Destruct */
template <class T>
void Delete(T* mem){
  delete mem;
}

template <class T>
bool _equal(T lhs, T rhs)
{
    return lhs == rhs;
}

template <class T>
void _from_str(T* dest, char* src){
  std::istringstream out(src);
  out >> *dest;
}

template <class T>
PyObject* _to_PyString(const T *x)
{
  std::ostringstream instore;
  instore << (*x);
  return PyString_FromString(instore.str().data());
}

#endif
