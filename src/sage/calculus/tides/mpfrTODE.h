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

#ifndef _taylorODE_H
#define _taylorODE_H
#ifndef real_MP
#define real_MP
#endif
#include <stdio.h>
#include <stdlib.h>
#include "mpfr.h"
#include <math.h>

#include "mpfrNUM.h"
#include "mpfrPOL.h"
#include "commonITER.h"


typedef long (*MPLinkedFunction)(iteration_data *itd, mpfr_t t, mpfr_t *v, mpfr_t *p, int orden, mpfr_t *coef);

typedef struct mp_DM {
	int rows;
	int columns;
	mpfr_t **data;
} mp_data_matrix;


void mp_set_relminstep(mpfr_t val);

void mp_set_info_taylor(void);
void mp_unset_info_taylor(void);
void mp_str_info_taylor(void);
void mp_add_info_step(mpfr_t tstep);

int getOrder (void);
int getNsteps (void);

int mp_taylor_order(mpfr_t eps);

void mp_norm_inf_vec(mpfr_t *rop, int nvar, mpfr_t *y) ;
void  mp_norm_inf(mpfr_t *rop, int n, int k, mpfr_t *coef, int MAX_ORDER);

void mp_compute_step(mpfr_t *rop, double *hant, mpfr_t tol, int n, int ord, mpfr_t *coef, int MAX_ORDER);

void mp_compute_tol (mpfr_t *tol, mpfr_t tolrel, mpfr_t tolabs, int nvar, mpfr_t *y); 

void mp_taylor_horner(int n, int ord, mpfr_t *coef, mpfr_t t, mpfr_t *x, int MAX_ORDER);
void mp_taylor_horner_der(int n, int ord, mpfr_t *coef, mpfr_t t, mpfr_t *x, int MAX_ORDER);

void mp_write_taylor_solution( int n, int j, mpfr_t tini, 
								  mpfr_t *x,  mpfr_t*** mat, FILE* fileout);
	

int mp_valid_step (MPLinkedFunction fcn, iteration_data *itd,  mpfr_t *step, mpfr_t tip, 
					   mpfr_t eps, int nvar, int ncol, int order,
					   mpfr_t *cvfd, mpfr_t *p, int MAX_ORDER);

int  mp_tides(MPLinkedFunction fcn, int *pdd,
			  int nvar, int npar, int nfun, 
			  mpfr_t *x, mpfr_t *p,
			  mpfr_t *lt, int ntes, 	
			  mpfr_t tolrel, mpfr_t tolabs, 
			  mp_data_matrix *mat, FILE* fileout);

int  mp_tides_point(MPLinkedFunction fcn, int *pdd,
					int nvar, int npar, int nfun, 
					mpfr_t *x, mpfr_t *p,
					mpfr_t t0, mpfr_t tf, mpfr_t dt, 	
					mpfr_t tolrel, mpfr_t tolabs,   
					mp_data_matrix *mat, FILE* fileout);

int  mp_tides_list(MPLinkedFunction fcn, int *pdd,
				   int nvar, int npar, int nfun, 
				   mpfr_t *x, mpfr_t *p,
				   mpfr_t *lt, int ntes, 	
				   mpfr_t tolrel, mpfr_t tolabs,   
				   mp_data_matrix *mat, FILE* fileout) ;

int  mp_tides_delta(MPLinkedFunction fcn, int *pdd,
					int nvar, int npar, int nfun, 
					mpfr_t *x, mpfr_t *p,
					mpfr_t t0, mpfr_t dt, int ntot,	
					mpfr_t tolrel, mpfr_t tolabs,   
					mp_data_matrix *mat, FILE* fileout);

int  mp_tides_delta_ft(MPLinkedFunction fcn, int *pdd, 
					int nvar, int npar, int nfun, 
					mpfr_t *x, mpfr_t *p,
					mpfr_t t0, mpfr_t tf, mpfr_t dt, 	
					mpfr_t tolrel, mpfr_t tolabs,   
					mp_data_matrix *mat, FILE* fileout);

int  mp_tides_kernel(MPLinkedFunction fcn,  int *pdd,
					 int nvar, int npar, int nfun, 
					 mpfr_t *x, mpfr_t *p,
					 mpfr_t t0, mpfr_t dt,  mpfr_t *dtl, int ntot,	
					 mpfr_t tolrel, mpfr_t tolabs,   
					 mp_data_matrix *mat, FILE* fileout, mp_data_matrix *der);

long mp_number_of_columns(MPLinkedFunction fcn);


void init_mp_data_matrix(mp_data_matrix *dm, int r, int c);
void delete_mp_data_matrix(mp_data_matrix *dm);

#endif
