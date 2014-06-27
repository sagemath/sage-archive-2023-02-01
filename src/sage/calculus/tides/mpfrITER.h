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
#include "mpfr.h"
#include "mpfrNUM.h"



		
void	varMP_init(iteration_data *itd, mpfr_t *var, mpfr_t *v, mpfr_t t);
void	parMP_init(iteration_data *itd, mpfr_t *par, mpfr_t *p);
void	linkMP_init(iteration_data *itd, mpfr_t *lk);
void	derMP_init(iteration_data *itd, mpfr_t *var, mpfr_t *par, mpfr_t *v);
void 	clear_mpfr_vpl (iteration_data *itd, mpfr_t *var, mpfr_t *par, mpfr_t *link);
void	write_sol_MP(iteration_data *itd, mpfr_t *cvfd, mpfr_t *var, mpfr_t *link);

void	mpfrts_htilde(iteration_data *itd, mpfr_t *h, long j, long v, long i, mpfr_t *ht, int ORDER_INDEX);

void	mpfrts_var_t(iteration_data *itd, mpfr_t *f, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_var_t_c(iteration_data *itd, char* cs, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_var_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *w, int ORDER_INDEX);

void	mpfrts_add_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_sub_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_add_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_sub_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_add_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_sub_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);


void	mpfrts_mul_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_mul_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_mul_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void    mpfrts_div_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_div_t_vc(iteration_data *itd, mpfr_t *u, mpfr_t c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_div_t_cv(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);


void	mpfrts_inv_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_exp_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_pow_t_c(iteration_data *itd, mpfr_t *u, char* cs, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_pow_t_cc(iteration_data *itd, mpfr_t *u, mpfr_t c, mpfr_t *w, int ORDER_INDEX);

void	mpfrts_sct_0(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i);
void	mpfrts_sct_i(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i, int ORDER_INDEX);
void	mpfrts_sin_t(iteration_data *itd, mpfr_t *s, mpfr_t *c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_cos_t(iteration_data *itd, mpfr_t *c, mpfr_t *s, mpfr_t *w, int ORDER_INDEX);
void    mpfrts_sin_cos_t (iteration_data *itd, mpfr_t *f, mpfr_t *s, mpfr_t *c, int ORDER_INDEX);
void	mpfrts_sinh_t(iteration_data *itd, mpfr_t *s, mpfr_t *c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_cosh_t(iteration_data *itd, mpfr_t *c, mpfr_t *s, mpfr_t *w, int ORDER_INDEX);
void    mpfrts_sinh_cosh_t (iteration_data *itd, mpfr_t *f,  mpfr_t *s, mpfr_t *c, int ORDER_INDEX);


void	mpfrts_fgt_0(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i);
void	mpfrts_fgt_i(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i, int ORDER_INDEX);
void	mpfrts_log_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_asin_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_acos_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_atan_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_asinh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_acosh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_atanh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);


#endif

