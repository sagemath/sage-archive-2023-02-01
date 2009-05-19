/*******************************************************************************
    SAGE: Open Source Mathematical Software
        Copyright (C) 2008 William Stein <wstein@gmail.com>
        Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
   Distributed under the terms of the GNU General Public License (GPL),
   version 2 or any later version.  The full text of the GPL is available at:
                   http://www.gnu.org/licenses/
*******************************************************************************/

#include "ccobject.h"
#include <pynac/ginac.h>

#include <string>

using namespace GiNaC;

void list_symbols(const ex& e, std::set<ex, ex_is_less> &s)
{
    if (is_a<symbol>(e)) {
        s.insert(e);
    } else {
        for (size_t i=0; i<e.nops(); i++)
            list_symbols(e.op(i), s);
    }
}


ex g_function_evalv(unsigned serial, exvector& vec)
{
    return ex(function(serial, vec));
}

ex g_function_eval0(unsigned serial)
{
    return ex(function(serial));
}

ex g_function_eval1(unsigned serial, const ex& arg1)
{
    return ex(function(serial, arg1));
}

ex g_function_eval2(unsigned serial, const ex& arg1, const ex& arg2)
{
    return ex(function(serial, arg1, arg2));
}

ex g_function_eval3(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3)
{
    return ex(function(serial, arg1, arg2, arg3));
}

ex g_function_eval4(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, const ex& arg4)
{
    return ex(function(serial, arg1, arg2, arg3, arg4));
}

ex g_function_eval5(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, const ex& arg4, const ex& arg5)
{
    return ex(function(serial, arg1, arg2, arg3, arg4, arg5));
}

ex g_function_eval6(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, const ex& arg4, const ex& arg5, const ex& arg6)
{
    return ex(function(serial, arg1, arg2, arg3, arg4, arg5, arg6));
}

ex g_function_eval7(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, const ex& arg4, const ex& arg5, const ex& arg6,
        const ex& arg7)
{
    return ex(function(serial, arg1, arg2, arg3, arg4, arg5, arg6, arg7));
}

ex g_function_eval8(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, const ex& arg4, const ex& arg5, const ex& arg6,
        const ex& arg7, const ex& arg8)
{
    return ex(function(serial, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8));
}

ex g_function_eval9(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, const ex& arg4, const ex& arg5, const ex& arg6,
        const ex& arg7, const ex& arg8, const ex& arg9)
{
    return ex(function(serial, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8,
                arg9));
}

ex g_function_eval10(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, const ex& arg4, const ex& arg5, const ex& arg6,
        const ex& arg7, const ex& arg8, const ex& arg9, const ex& arg10)
{
    return ex(function(serial, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8,
                arg9, arg10));
}

ex g_function_eval11(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, const ex& arg4, const ex& arg5, const ex& arg6,
        const ex& arg7, const ex& arg8, const ex& arg9, const ex& arg10,
        const ex& arg11)
{
    return ex(function(serial, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8,
                arg9, arg10, arg11));
}

ex g_function_eval12(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, const ex& arg4, const ex& arg5, const ex& arg6,
        const ex& arg7, const ex& arg8, const ex& arg9, const ex& arg10,
        const ex& arg11, const ex& arg12)
{
    return ex(function(serial, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8,
                arg9, arg10, arg11, arg12));
}

ex g_function_eval13(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, const ex& arg4, const ex& arg5, const ex& arg6,
        const ex& arg7, const ex& arg8, const ex& arg9, const ex& arg10,
        const ex& arg11, const ex& arg12, const ex& arg13)
{
    return ex(function(serial, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8,
                arg9, arg10, arg11, arg12, arg13));
}

ex g_function_eval14(unsigned serial, const ex& arg1, const ex& arg2,
        const ex& arg3, const ex& arg4, const ex& arg5, const ex& arg6,
        const ex& arg7, const ex& arg8, const ex& arg9, const ex& arg10,
        const ex& arg11, const ex& arg12, const ex& arg13, const ex& arg14)
{
    return ex(function(serial, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8,
                arg9, arg10, arg11, arg12, arg13, arg14));
}

double GEx_to_double(ex& e, int* success) {
  /*
     Convert an expression to a double if possible.

     INPUT:
         e -- an expression
         success -- pointer to int
     OUTPUT:
         double -- the answer,
         *AND* sets success to true on success, and false on failure
         if e is not coercible to a double
  */
  ex f = e.evalf(0, 53);
  if (is_a<numeric>(f)) {
    *success = true;
    return (ex_to<numeric>(f)).to_double();
  } else {
    *success = false;
    return 0;
  }
}

bool relational_to_bool(const ex& e) {
    if (ex_to<relational>(e))
	return 1;
    else
	return 0;
}

relational::operators relational_operator(const ex& e) {
    // unsafe cast -- be damn sure the input is a relational.
    return (ex_to<relational>(e)).the_operator();
}

relational::operators switch_operator(relational::operators o) {
    switch(o) {
    case relational::less:
	return relational::greater;
    case relational::less_or_equal:
	return relational::greater_or_equal;
    case relational::greater:
	return relational::less;
    case relational::greater_or_equal:
	return relational::less_or_equal;
    default:
	return o;
    }
}

bool is_negative(ex x) {
    if (is_a<numeric>(x)) {
	return (ex_to<numeric>(x)).is_negative();
    }
    return false;
}

PyObject* py_object_from_numeric(ex x) {
    return (ex_to<numeric>(x)).to_pyobject();
}

template <class T>
PyObject* _to_PyString_latex(const T *x)
{
  std::ostringstream instore;
  instore << latex << (*x);
  return PyString_FromString(instore.str().data());
}

#define ASSIGN_WRAP(x,y) x = y

#define ADD_WRAP(x,y) (x)+(y)
#define SUB_WRAP(x,y) (x)-(y)
#define MUL_WRAP(x,y) (x)*(y)
#define DIV_WRAP(x,y) (x)/(y)

#define LT_WRAP(x,y)  (x)<(y)
#define EQ_WRAP(x,y)  (x)==(y)
#define GT_WRAP(x,y)  (x)>(y)
#define LE_WRAP(x,y)  (x)<=(y)
#define NE_WRAP(x,y)  (x)!=(y)
#define GE_WRAP(x,y)  (x)>=(y)

