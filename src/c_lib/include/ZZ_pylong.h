/*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL) as
#  published by the Free Software Foundation; either version 2 of the
#  License, or (at your option) any later version.
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************/

/*  Author:  Joel B. Mohler <joel@kiwistrawberry.us>
			 2007-06-17                                                      */

#ifndef ZZ_PYLONG_H
#define ZZ_PYLONG_H

/* Yes, this is kind of weird.  I'm only wrapping this for C++ */
#ifdef __cplusplus

#include <Python.h>
#include <NTL/ZZ.h>
using namespace NTL;
#include "gmp.h"

/* ZZ -> pylong conversion */
PyObject * ZZ_get_pylong(ZZ &z);

/* pylong -> ZZ conversion */
int ZZ_set_pylong(ZZ &z, PyObject * ll);
#endif

#endif
