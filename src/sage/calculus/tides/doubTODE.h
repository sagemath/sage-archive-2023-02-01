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
#ifndef real_Double
#define real_Double
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "doubNUM.h"
#include "doubPOL.h"
#include "commonITER.h"



typedef long (*DBLinkedFunction)(iteration_data *itd, double t, double *v, double *p, int orden, double *cvfd);

typedef struct dp_DM {
	int rows;
	int columns;
	double **data;
} dp_data_matrix;

void	check_relminstep(void);
void	uncheck_relminstep(void);
void	dp_set_relminstep(double val);


void use_default_step_estimator (void);
void dp_set_info_taylor(void);
void dp_unset_info_taylor(void);
void dp_str_info_taylor(void);
void dp_add_info_step(double tstep);


int  dp_taylor_order(double eps);

void dp_norm_inf_var(double *rop, int nvar, double *y) ;
void dp_norm_inf(double *rop, int n, int k, double *coef, int MAX_ORDER);

void dp_compute_step(double *rop, double *hant, double tol, int n, int ord, double *cvfd, int MAX_ORDER);

void dp_compute_tol(double *tol, double tolrel, double tolabs, int nvar, double *y) ;

void dp_taylor_horner(int n, int ord, double *cvfd, double t, double *x, int MAX_ORDER);
void dp_taylor_horner_der (int n, int ord, double *cvfd, double t, double *x, int MAX_ORDER);

void dp_write_taylor_solution(int n,  int j, double tini, double *x, double*** mat, FILE* fileout);

int  dp_valid_step (DBLinkedFunction fcn, iteration_data *itd,  double *step, double tip, 
				double tol, int nvar, int ncol, int order, double *cvfd,
				double *p, int MAX_ORDER);

int  dp_tides(DBLinkedFunction fcn, int *pdd,
			  int nvar, int npar, int nfun, 
			  double *x, double *p,
			  double *lt, int ntes, 	
			  double tolrel, double tolabs,   
			  dp_data_matrix *mat, FILE* fileout);

int  dp_tides_point(DBLinkedFunction fcn, int *pdd,
					int nvar, int npar, int nfun, 
					double *x, double *p,
					double t0, double tf, double dt, 	
					double tolrel, double tolabs,   
					dp_data_matrix *mat, FILE* fileout);

int  dp_tides_list(DBLinkedFunction fcn, int *pdd,
				   int nvar, int npar, int nfun, 
				   double *x, double *p,
				   double *lt, int ntes, 	
				   double tolrel, double tolabs,   
				   dp_data_matrix *mat, FILE* fileout) ;

int  dp_tides_delta(DBLinkedFunction fcn, int *pdd,
					int nvar, int npar, int nfun, 
					double *x, double *p,
					double t0, double dt, int ntot,	
					double tolrel, double tolabs,   
					dp_data_matrix *mat, FILE* fileout);

int  dp_tides_delta_ft(DBLinkedFunction fcn, int *pdd,
					int nvar, int npar, int nfun, 
					double *x, double *p,
					double t0, double tf, double dt, 	
					double tolrel, double tolabs,   
					dp_data_matrix *mat, FILE* fileout);


int  dp_tides_kernel(DBLinkedFunction fcn, int *pdd,
					 int nvar, int npar, int nfun, 
					 double *x, double *p,
					 double t0, double dt, double *dtl, int ntot, 	
					 double tolrel, double tolabs,   
					 dp_data_matrix *mat, FILE* fileout, dp_data_matrix *der) ;

long dp_number_of_columns(DBLinkedFunction fcn);


void init_dp_data_matrix(dp_data_matrix *dm, int r, int c);
void delete_dp_data_matrix(dp_data_matrix *dm);



#endif
