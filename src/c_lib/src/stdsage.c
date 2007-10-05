/**
 * @file stdsage.c
 *
 * Some global C stuff that gets imported into pyrex modules.
 *
 */

/******************************************************************************
       Copyright (C) 2006 William Stein <wstein@gmail.com>
                     2006 David Harvey <dmharvey@math.harvard.edu>
                     2006 Martin Albrecht <malb@informatik.uni-bremen.de>

  Distributed under the terms of the GNU General Public License (GPL), Version 2.

  The full text of the GPL is available at:
                  http://www.gnu.org/licenses/

******************************************************************************/


#include "stdsage.h"
#include "interrupt.h"

PyObject* global_empty_tuple;

void init_global_empty_tuple(void) {
  _CALLED_ONLY_ONCE;

  global_empty_tuple = PyTuple_New(0);
}


void init_csage(void) {
  //_CALLED_ONLY_ONCE;

  init_global_empty_tuple();
  setup_signal_handler();
}
