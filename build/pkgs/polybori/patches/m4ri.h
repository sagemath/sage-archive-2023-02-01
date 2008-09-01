/**
 * \file m4ri.h
 * \brief Main include file for the M4RI library.
 *
 * \author Gregory Bard <bard@fordham.edu>
 * \author Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
 */
/******************************************************************************
*
*            M4RI: Method of the Four Russians Inversion
*
*       Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org>
*       Copyright (C) 2007,2008 Martin Albrecht <malb@informatik.uni-bremen.de>
*
*  Distributed under the terms of the GNU General Public License (GPL)
*
*    This code is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    General Public License for more details.
*
*  The full text of the GPL is available at:
*
*                  http://www.gnu.org/licenses/
******************************************************************************/
#ifndef M4RI_H
#define M4RI_H

/**
 * \mainpage
 *
 * M4RI is a library to do fast arithmetic with dense matrices over
 * \f$F_2\f$. M4RI is available under the GPLv2+ and used by the Sage
 * mathematics software and the PolyBoRi library. See
 * http://m4ri.sagemath.org for details.
 *
 * \example testsuite/test_multiplication.c
 * \example testsuite/test_reduction.c
 */

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "watch.h"
#include "packedmatrix.h"
#include "brilliantrussian.h"
#include "strassen.h"
#include "grayflex.h"
#include "parity.h"

#ifdef __cplusplus
}
#endif //__cplusplus

#endif //M4RI_H
