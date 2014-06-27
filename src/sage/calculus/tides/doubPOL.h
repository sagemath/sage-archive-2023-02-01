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

/*  V-20  */

#ifndef _doublePOL_H
#define _doublePOL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void dp_pol_derivative(double *pol, int grado, double *dpol);

void dp_horner(double *pol, int grado, double t, double *eval);
void dp_horner_der(double *pol, int grado, double t, double *eval);

void dp_classic_horner(double *pol, int grado, double t, double *eval);

void dp_twosum(double a, double b, double *x, double *y);
void dp_split_real(double a, double *ah, double *al);
void dp_twoproduct(double a, double b, double *x, double *y);
void dp_compensated_horner(double *pol, int grado, double t, double *eval);

void dp_ch_evaluate_poly_der(double *pol, int grado, double t, double *vpol, double *vder);



int  dp_one_zero(double *pol, int grado, double ta, double tb, double tol, double *raiz);
int  dp_one_extremum(double *pol, int grado, double ta, double tb, double tol, double *extr);
int  dp_one_maximum(double *pol, int grado, double ta, double tb, double tol, double *extr);
int  dp_one_minimum(double *pol, int grado, double ta, double tb, double tol, double *extr);
int  dp_brent(double *pol, int grado, double ta, double tb, double tol, double *raiz);
int  dp_rtsafe(double *pol, int grado, double x1, double x2, double tol, double *raiz);



#endif

