/**
 * @file stdsage.c
 *
 * Some global C stuff that gets imported into pyrex modules.
 *
 */

/******************************************************************************
 *     Copyright (C) 2006 William Stein <wstein@gmail.com>
 *                   2006 David Harvey <dmharvey@math.harvard.edu>
 *                   2006 Martin Albrecht <malb@informatik.uni-bremen.de>
 *                   2011 Jeroen Demeyer <jdemeyer@cage.ugent.be>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/


#include "stdsage.h"
#include "interrupt.h"

PyObject* global_empty_tuple;

void init_global_empty_tuple(void) {
  _CALLED_ONLY_ONCE;

  global_empty_tuple = PyTuple_New(0);
}


/* This is called once during Sage startup. On some platforms like
 * Cygwin, this is also called from init_csage_module(). */
void init_csage() {
    init_global_empty_tuple();
    setup_sage_signal_handler();
}
