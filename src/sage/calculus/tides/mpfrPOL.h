/****************************************************************************
 libTIDES. 
 This file is part of TIDES.
 
 Copyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
 Grupo de Mecanica Espacial
 University of Zaragoza
 SPAIN
 
 http://gme.unizar.es/software/tides
 Contact: <tides@unizar.es>
 
 TIDES is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 TIDES is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with TIDES.  If not, see <http://www.gnu.org/licenses/>.
 
 *****************************************************************************/

/*  V-20 */

#ifndef _mpfrPOL_H
#define _mpfrPOL_H

#include <stdio.h>
#include <stdlib.h>
#include "mpfr.h"
#include <math.h>

#include "mpfrNUM.h"

void mp_pol_derivative(mpfr_t *pol, int grado, mpfr_t *dpol);

void mp_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);
void mp_horner_der(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);

void mp_classic_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);

void mp_twosum(mpfr_t a, mpfr_t b, mpfr_t *x, mpfr_t *y);
void mp_split_real(mpfr_t a, mpfr_t *ah, mpfr_t *al);
void mp_twoproduct(mpfr_t a, mpfr_t b, mpfr_t *x, mpfr_t *y);
void mp_compensated_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);

void mp_ch_evaluate_poly_der(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *vpol, mpfr_t *vder);

int  mp_one_zero(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *raiz);
int  mp_one_extremum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr);
int  mp_one_maximum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr);
int  mp_one_minimum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr);
int  mp_brent(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *raiz);
int  mp_rtsafe(mpfr_t *pol, int grado, mpfr_t x1, mpfr_t x2, mpfr_t tol, mpfr_t *raiz);


#endif

