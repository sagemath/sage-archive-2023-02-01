/*

Some global C stuff that gets imported into pyrex modules.

*/

/******************************************************************************
       Copyright (C) 2006 William Stein <wstein@gmail.com>
                     2006 David Harvey <dmharvey@math.harvard.edu>

  Distributed under the terms of the GNU General Public License (GPL), Version 2.

  The full text of the GPL is available at:
                  http://www.gnu.org/licenses/

******************************************************************************/

#include "Python.h"

/* A global empty python tuple object. This is used to speed up some python
   API calls where we want to avoid constructing a tuple every time. */
PyObject* global_empty_tuple;

/* This is called exactly once at startup from sage_object.pyx */
void init_global_empty_tuple(void)
{
  global_empty_tuple = PyTuple_New(0);
}
