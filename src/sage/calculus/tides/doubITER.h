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

#ifndef _mpfrITER_H
#define _mpfrITER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "doubNUM.h"
#include "commonITER.h"


void	varDB_init(iteration_data *itd, double *var,  double *v, double t);
void	parDB_init(iteration_data *itd, double *par, double *p);
void	linkDB_init(iteration_data *itd, double *lk);
void	derDB_init(iteration_data *itd, double *var, double *par, double *v);
void	write_sol_DB(iteration_data *itd, double *cvfd, double *var, double *link);

void	double_htilde(iteration_data *itd, double *h, long j, long v, long i, double *ht, int ORDER_INDEX);

void	double_var_t(iteration_data *itd, double *f, double *u, int ORDER_INDEX);
void	double_var_t_c(iteration_data *itd, char* cs, double *w, int ORDER_INDEX);
void	double_var_t_cc(iteration_data *itd, double c, double *w, int ORDER_INDEX);

void	double_add_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX);
void	double_sub_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX);
void	double_add_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX);
void	double_sub_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX);
void	double_add_t_cc(iteration_data *itd, double c, double *u,  double *w, int ORDER_INDEX);
void	double_sub_t_cc(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX);


void	double_mul_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX);
void	double_mul_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX);
void	double_mul_t_cc(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX);
void	double_div_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_div_t_vc(iteration_data *itd, double *u, double c, double *w, int ORDER_INDEX);
void	double_div_t_cv(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX);

void	double_inv_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX);
void	double_exp_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX);
void	double_pow_t_c(iteration_data *itd, double *u, char* cs, double *w, int ORDER_INDEX);
void	double_pow_t_cc(iteration_data *itd, double *u, double c, double *w, int ORDER_INDEX);

void	double_sct_0(iteration_data *itd, double *f, double *g, double *h, long i);
void	double_sct_i(iteration_data *itd, double *f, double *g, double *h, long i, int ORDER_INDEX);
void	double_sin_t(iteration_data *itd, double *s, double *c, double *w, int ORDER_INDEX);
void	double_cos_t(iteration_data *itd, double *c, double *s, double *w, int ORDER_INDEX);
void    double_sin_cos_t (iteration_data *itd, double *f, double *s, double *c,  int ORDER_INDEX);
void	double_sinh_t(iteration_data *itd, double *s, double *c, double *w, int ORDER_INDEX);
void	double_cosh_t(iteration_data *itd, double *c, double *s, double *w, int ORDER_INDEX);
void    double_sinh_cosh_t (iteration_data *itd, double *f, double *s, double *c, int ORDER_INDEX);

void	double_fgt_0(iteration_data *itd, double *f, double *g, double *h, long i);
void	double_fgt_i(iteration_data *itd, double *f, double *g, double *h, long i, int ORDER_INDEX);
void	double_log_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX);
void	double_asin_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_acos_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_atan_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_asinh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_acosh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_atanh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);


#endif

