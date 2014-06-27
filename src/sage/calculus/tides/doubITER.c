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

#include "doubITER.h"

extern int	printError;

void varDB_init(iteration_data *itd, double *var,  double *v, double t)
{
	int i, j, k;
	for (i=0; i<=itd->NVARS; i++) 
		for (j=0; j<itd->NDER; j++)
			for (k=0; k< itd->MAX_ORDER; k++)
				*(var + i*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k) = 0.;
	double_set (var, t);
	double_set (var+1, 1.);
	for(i=1; i<=itd->NVARS; i++) 
		double_set (var + i*itd->LINK_ELEMENTS , v[i-1]);
}
void parDB_init(iteration_data *itd, double *par, double *p)
{
	int i, j, k;
	for (i=0; i<itd->NPARS; i++) 
		for (j=0; j<itd->NDER; j++)
			for (k=0; k< itd->MAX_ORDER; k++)
				*(par + i*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k) = 0.;
	for(i=0; i<itd->NPARS; i++) 
		double_set (par + i*itd->LINK_ELEMENTS , p[i]);
	for(i=0; i <itd->NPARTIALS; i++)
		if(!is_variable(itd, itd->PARTIAL_LIST[i]))
			*(par + (itd->PARTIAL_LIST[i]-itd->NVARS-1)*itd->LINK_ELEMENTS +(i+1)*itd->MAX_ORDER ) = 1.;
	
}
void linkDB_init(iteration_data *itd, double *lk)
{
	int i, j, k;
	for (i=0; i<itd->NLINKS; i++) 
		for (j=0; j<itd->NDER; j++)
			for (k=0; k< itd->MAX_ORDER; k++)
				*(lk + i*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k) = 0.;
}


void derDB_init(iteration_data *itd, double *var, double *par, double *v)
{
	int i; 
	if (itd->clearPartials) {
		for(i=0; i <itd->NPARTIALS; i++) {
			if(is_variable(itd, itd->PARTIAL_LIST[i])) 
				*(var + itd->PARTIAL_LIST[i]*itd->LINK_ELEMENTS +(i+1)*itd->MAX_ORDER ) = 1.;
		}
	}else {
		int j, nvf;
		nvf = itd->NVARS+itd->NFUNS;
		for( i = 0; i < itd->NVARS; i++) 
			for(j = 1; j < itd->NDER; j++) 
				*(var + (i+1)*itd->LINK_ELEMENTS +j*itd->MAX_ORDER ) = v[i+(j*nvf)];
	}

	itd->clearPartials = 0;
}


void	write_sol_DB(iteration_data *itd, double *cvfd, double *var, double *link)
{
	int i,j,k, nvf;
	double vfun;
	nvf = itd->NVARS+itd->NFUNS;
	for( i = 0; i < itd->NVARS; i++) 
		for(j = 0; j < itd->NDER; j++) 
			for(k = 0; k < itd->MAX_ORDER; k++)
				double_set (&cvfd[(i+(j*nvf))*itd->MAX_ORDER + k], *(var + (i+1)*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k));
	for( i = 0; i < itd->NFUNS; i++) 
		for(j = 0; j < itd->NDER; j++) 
			for(k = 0; k < itd->MAX_ORDER; k++) {
				if (itd->FUNCTION_LIST[i] >= 0) vfun = *(link + itd->FUNCTION_LIST[i]*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k);
				else vfun = *(var -itd->FUNCTION_LIST[i]*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k);
				double_set (&cvfd[(itd->NVARS+i+(j*nvf))*itd->MAX_ORDER +  k], vfun);
			}
}

/************************************************************************/

void double_htilde(iteration_data *itd, double *h, long j, long v, long i, double *ht, int ORDER_INDEX)
{
	double cero;
	double_init (&cero);
	if(ORDER_INDEX >= 0) {
		if( j == ORDER_INDEX && i==v) double_set(ht, cero);
		else double_set(ht, h[v*itd->MAX_ORDER+j]);
	}
}

/************************************************************************/

void	double_var_t(iteration_data *itd, double *f, double *w, int ORDER_INDEX)
{
	double val; double_init (&val);
	int i;
	if(ORDER_INDEX > 0) {
		for(i = 0; i < itd->NDER; i++) {
			double_div_i(&val, f[i*itd->MAX_ORDER+ORDER_INDEX-1], ORDER_INDEX);
			double_set(&w[i*itd->MAX_ORDER+ORDER_INDEX], val);
		}
	}

}

void	double_var_t_c(iteration_data *itd, char* cs, double *w, int ORDER_INDEX)
{
	if(ORDER_INDEX == 1) {
		double_set_str(&w[1], cs);
	}
}
void	double_var_t_cc(iteration_data *itd, double c, double *w, int ORDER_INDEX)
{
	if(ORDER_INDEX == 1) {
		double_set(&w[1], c);
	}
}

void	double_add_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < itd->NDER; i++) 
			double_add(&w[i*itd->MAX_ORDER+ORDER_INDEX], u[i*itd->MAX_ORDER+ORDER_INDEX], v[i*itd->MAX_ORDER+ORDER_INDEX]);
	}
}

void	double_sub_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < itd->NDER; i++) 
			double_sub(&w[i*itd->MAX_ORDER+ORDER_INDEX], u[i*itd->MAX_ORDER+ORDER_INDEX], v[i*itd->MAX_ORDER+ORDER_INDEX]);
	}
}



void	double_add_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX)
{
	long i;
	double c; 
	double_init (&c);
	if(ORDER_INDEX >= 0) {
		double_set_str(&c, cs);
		if(ORDER_INDEX == 0) 
			double_add(&w[ORDER_INDEX], c, u[ORDER_INDEX]);
		else 
			double_set(&w[ORDER_INDEX], u[ORDER_INDEX]);
		for(i = 1; i < itd->NDER; i++) 
			double_set(&w[i*itd->MAX_ORDER+ORDER_INDEX], u[i*itd->MAX_ORDER+ORDER_INDEX]);
	}

}
void	double_add_t_cc(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) 
			double_add(&w[ORDER_INDEX], c, u[ORDER_INDEX]);
		else 
			double_set(&w[ORDER_INDEX], u[ORDER_INDEX]);
		for(i = 1; i < itd->NDER; i++) 
			double_set(&w[i*itd->MAX_ORDER+ORDER_INDEX], u[i*itd->MAX_ORDER+ORDER_INDEX]);
	}
}


void	double_sub_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX)
{
	long i;
	double c; double_init (&c);
	if(ORDER_INDEX >= 0) {
		double_set_str(&c, cs);
		if(ORDER_INDEX == 0) 
			double_sub(&w[ORDER_INDEX], c, u[ORDER_INDEX]);
		else
			double_mul_i(&w[ORDER_INDEX], u[ORDER_INDEX],-1);
		for(i = 1; i < itd->NDER; i++) 
			double_mul_i(&w[i*itd->MAX_ORDER+ORDER_INDEX], u[i*itd->MAX_ORDER+ORDER_INDEX],-1);
	}

}
void	double_sub_t_cc(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) 
			double_sub(&w[ORDER_INDEX], c, u[ORDER_INDEX]);
		else
			double_mul_i(&w[ORDER_INDEX], u[ORDER_INDEX],-1);
		for(i = 1; i < itd->NDER; i++) 
			double_mul_i(&w[i*itd->MAX_ORDER+ORDER_INDEX], u[i*itd->MAX_ORDER+ORDER_INDEX],-1);
	}	
}



void	double_mul_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX)
{
	long i, vi, j;
	double sumint, sumext, partial;
	if(ORDER_INDEX >= 0) {
		sumint = sumext = partial = 0.;
		for(i = 0; i < itd->NDER; i++) {
			sumext = 0.;
			for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
				sumint = 0.;
				for(j = 0; j <= ORDER_INDEX; j++)
					sumint = sumint + u[itd->PREV_VI[vi]*itd->MAX_ORDER+ORDER_INDEX-j] * v[itd->PREV_IV[vi]*itd->MAX_ORDER+j];
				sumext = sumext + sumint*itd->PREV_COEF[vi];
			}
			w[i*itd->MAX_ORDER+ORDER_INDEX] = sumext;	
		}
	}
}


void	double_mul_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX)
{
	long i;
	double c; double_init (&c);
	if(ORDER_INDEX >= 0) {
		double_set_str(&c, cs);
		for(i = 0; i < itd->NDER; i++) 
			double_mul(&w[i*itd->MAX_ORDER+ORDER_INDEX], c, u[i*itd->MAX_ORDER+ORDER_INDEX]);
	}

}
void	double_mul_t_cc(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < itd->NDER; i++) 
			double_mul(&w[i*itd->MAX_ORDER+ORDER_INDEX], c, u[i*itd->MAX_ORDER+ORDER_INDEX]);
	}
	
}



void double_div_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX)
{
	long i, j, vi;
	double sumint, sumext, partial, zero, ht;
	if(g[0] == 0 ) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********   double_div_t   ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		g[0]= 1;
	}
	if(ORDER_INDEX >= 0) {
		double_init(&sumint);
		double_init(&sumext);
		double_init(&partial);
		double_init(&zero);
		for(i = 0; i < itd->NDER; i++) {
			double_set(&sumext, zero);
			for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
				double_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					double_htilde(itd, h, ORDER_INDEX-j, itd->PREV_VI[vi], i, &ht, ORDER_INDEX);
					double_mul(&partial, ht, g[itd->PREV_IV[vi]*itd->MAX_ORDER+j]);
					double_add(&sumint, sumint, partial);
				}
				double_mul(&sumint, sumint, itd->PREV_COEF[vi]);
				double_add(&sumext, sumext, sumint);
			}
			double_sub(&sumext, f[i*itd->MAX_ORDER+ORDER_INDEX], sumext);
			double_div(&sumext, sumext, g[0]);
			double_set(&h[i*itd->MAX_ORDER+ORDER_INDEX], sumext);
		}
	}
}


void	double_div_t_vc(iteration_data *itd, double *u, double c, double *w, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < itd->NDER; i++) 
			double_mul(&w[i*itd->MAX_ORDER+ORDER_INDEX], 1./c, u[i*itd->MAX_ORDER+ORDER_INDEX]);
	}
}

void	double_div_t_cv(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX)
{
	long i, j, vi;
	double sumint, sumext, partial, zero, f, wt;
	double_init (&zero);
	if(double_equal (u[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********  double_div_t_cv ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		double_set_d (&u[0], 1.);
	}
	if(ORDER_INDEX >= 0) {
		double_init(&sumint);
		double_init(&sumext);
		double_init(&partial);
		double_init(&f);
		double_init(&wt);
		for(i = 0; i < itd->NDER; i++) {
			if(ORDER_INDEX == 0 && i == 0 ) double_set_d(&f, 1.);
			else double_set(&f, zero);
			double_set(&sumext, zero);
			for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
				double_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					double_htilde(itd, w, ORDER_INDEX-j, itd->PREV_VI[vi], i, &wt, ORDER_INDEX);
					double_div(&wt, wt, c);
					double_mul(&partial, wt, u[itd->PREV_IV[vi]*itd->MAX_ORDER+j]);
					double_add(&sumint, sumint, partial);
				}
				double_mul_i(&sumint, sumint, itd->PREV_COEF[vi]);
				double_add(&sumext, sumext, sumint);
			}
			double_sub(&sumext, f, sumext);
			double_div(&sumext, sumext, u[0]);
			double_mul(&sumext, sumext, c);
			double_set(&w[i*itd->MAX_ORDER+ORDER_INDEX], sumext);
		}
	}
}

void double_inv_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX)
{
	long i, j, vi;
	double sumint, sumext, partial, zero, f, wt;
	double_init (&zero);
	if(double_equal (u[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********   double_inv_t   ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		double_set_d (&u[0], 1.);
	}
	if(ORDER_INDEX >= 0) {
		double_init(&sumint);
		double_init(&sumext);
		double_init(&partial);
		double_init(&f);
		double_init(&wt);
		for(i = 0; i < itd->NDER; i++) {
			if(ORDER_INDEX == 0 && i == 0 ) double_set_d(&f, 1.);
			else double_set(&f, zero);
			double_set(&sumext, zero);
			for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
				double_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					double_htilde(itd, w, ORDER_INDEX-j, itd->PREV_VI[vi], i, &wt, ORDER_INDEX);
					double_mul(&partial, wt, u[itd->PREV_IV[vi]*itd->MAX_ORDER+j]);
					double_add(&sumint, sumint, partial);
				}
				double_mul_i(&sumint, sumint, itd->PREV_COEF[vi]);
				double_add(&sumext, sumext, sumint);
			}
			double_sub(&sumext, f, sumext);
			double_div(&sumext, sumext, u[0]);
			double_set(&w[i*itd->MAX_ORDER+ORDER_INDEX], sumext);
		}
	}
}

void double_exp_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX)
{
	long i, j, vi;
	double sumint,sumext, partial, zero;
	if(ORDER_INDEX >= 0) {
		zero = 0.;
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_exp(&sumext, u[0]);
					double_set(&w[0], sumext);
				} else {
					double_set(&sumint, zero);					
					for( vi = itd->PREVSTAR_ACCUM[i]; vi < itd->PREVSTAR_ACCUM[i+1]; vi++) {
						double_mul(&partial, w[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER], u[itd->PREVSTAR_IV[vi]*itd->MAX_ORDER]);
						double_mul_i(&partial, partial, itd->PREVSTAR_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_set(&w[i*itd->MAX_ORDER], sumint);
				}
			}
		} else {		
			for(i = 0; i < itd->NDER; i++) {
				double_set(&sumext, zero);
				for(j = 0; j < ORDER_INDEX; j++) {
					double_set(&sumint, zero);
					for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
						double_mul(&partial, w[itd->PREV_VI[vi]*itd->MAX_ORDER+j], u[itd->PREV_IV[vi]*itd->MAX_ORDER+ORDER_INDEX-j]);
						double_mul_i(&partial, partial, itd->PREV_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_mul_i(&sumint, sumint, ORDER_INDEX-j);
					double_add(&sumext, sumext, sumint);
				}
				double_div_i(&sumext, sumext, ORDER_INDEX);
				double_set(&w[i*itd->MAX_ORDER+ORDER_INDEX], sumext);
			}
			
		}
	}
}


void double_pow_t_c(iteration_data *itd, double *u, char* cs, double *w, int ORDER_INDEX)
{
	long i, j, vi;
	double sumint, sumext, partial, partialb, zero, wt, c;
	double_init (&zero);
	if(double_equal (u[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********  double_pow_t_c  ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		double_set_d (&u[0], 1.);
	}
	
	if(ORDER_INDEX >= 0) {
		double_init(&sumint);
		double_init(&sumext);
		double_init(&partial);
		double_init(&partialb);
		double_init(&wt);
		double_init(&c);
		double_set_str(&c, cs);
	
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_pow(&sumext, u[0], c);
					double_set(&w[0], sumext);
				} else {
					double_set(&sumint, zero);					
					for( vi = itd->PREVSTAR_ACCUM[i]; vi < itd->PREVSTAR_ACCUM[i+1]; vi++) {
						double_mul(&partial, w[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER], u[itd->PREVSTAR_IV[vi]*itd->MAX_ORDER]);
						double_mul(&partial, partial, c);
						double_htilde(itd, w, 0, itd->PREVSTAR_IV[vi], i, &wt, 0);
						double_mul(&partialb, wt, u[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER]);
						double_sub(&partial, partial, partialb);
						double_mul_i(&partial, partial, itd->PREVSTAR_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_div(&sumint, sumint, u[0]);
					double_set(&w[i*itd->MAX_ORDER], sumint);
				}
			}
		} else {		
			for(i = 0; i < itd->NDER; i++) {
				double_set(&sumext, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					double_set(&sumint, zero);
					for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
						double_htilde(itd, w, j, itd->PREV_VI[vi], i, &wt, ORDER_INDEX);
						double_mul(&partial, wt, u[itd->PREV_IV[vi]*itd->MAX_ORDER+ORDER_INDEX-j]);
						double_mul_i(&partial, partial, itd->PREV_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_add_i(&partial, c, 1);
					double_mul_i(&partial, partial, j);
					double_mul_i(&partialb, c, ORDER_INDEX);
					double_sub(&partial, partialb, partial);
					double_mul(&sumint, sumint, partial);
					double_add(&sumext, sumext, sumint);
				}
				double_div(&sumext, sumext, u[0]);
				double_div_i(&sumext, sumext, ORDER_INDEX);
				double_set(&w[i*itd->MAX_ORDER+ORDER_INDEX], sumext);
			}
			
		}
		
	}

	
	
}
void double_pow_t_cc(iteration_data *itd, double *u, double c, double *w, int ORDER_INDEX)
{
	long i, j, vi;
	double sumint, sumext, partial, partialb, zero, wt;
	double_init (&zero);
	if(double_equal (u[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********  double_pow_t_cc ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		double_set_d (&u[0], 1.);
	}
	
	if(ORDER_INDEX >= 0) {
		double_init(&sumint);
		double_init(&sumext);
		double_init(&partial);
		double_init(&partialb);
		double_init(&wt);
		
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_pow(&sumext, u[0], c);
					double_set(&w[0], sumext);
				} else {
					double_set(&sumint, zero);					
					for( vi = itd->PREVSTAR_ACCUM[i]; vi < itd->PREVSTAR_ACCUM[i+1]; vi++) {
						double_mul(&partial, w[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER], u[itd->PREVSTAR_IV[vi]*itd->MAX_ORDER]);
						double_mul(&partial, partial, c);
						double_htilde(itd, w, 0, itd->PREVSTAR_IV[vi], i, &wt, 0);
						double_mul(&partialb, wt, u[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER]);
						double_sub(&partial, partial, partialb);
						double_mul_i(&partial, partial, itd->PREVSTAR_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_div(&sumint, sumint, u[0]);
					double_set(&w[i*itd->MAX_ORDER], sumint);
				}
			}
		} else {		
			for(i = 0; i < itd->NDER; i++) {
				double_set(&sumext, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					double_set(&sumint, zero);
					for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
						double_htilde(itd, w, j, itd->PREV_VI[vi], i, &wt, ORDER_INDEX);
						double_mul(&partial, wt, u[itd->PREV_IV[vi]*itd->MAX_ORDER+ORDER_INDEX-j]);
						double_mul_i(&partial, partial, itd->PREV_COEF[vi]);
						double_add(&sumint, sumint, partial);
					}
					double_add_i(&partial, c, 1);
					double_mul_i(&partial, partial, j);
					double_mul_i(&partialb, c, ORDER_INDEX);
					double_sub(&partial, partialb, partial);
					double_mul(&sumint, sumint, partial);
					double_add(&sumext, sumext, sumint);
				}
				double_div(&sumext, sumext, u[0]);
				double_div_i(&sumext, sumext, ORDER_INDEX);
				double_set(&w[i*itd->MAX_ORDER+ORDER_INDEX], sumext);
			}
			
		}
		
	}
	
	
	
}

/************************************************************************/

void double_sct_0(iteration_data *itd, double *f, double *g, double *h, long i)
{
	long vi;
	double sumint, partial, zero;
	double_init(&zero);
	double_init(&partial);
	double_init(&sumint);
	double_set(&sumint, zero);					
	for( vi = itd->PREVSTAR_ACCUM[i]; vi < itd->PREVSTAR_ACCUM[i+1]; vi++) {
		double_mul(&partial, g[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER], f[itd->PREVSTAR_IV[vi]*itd->MAX_ORDER]);
		double_mul_i(&partial, partial, itd->PREVSTAR_COEF[vi]);		
		double_add(&sumint, sumint, partial);
	}
	double_set(&h[i*itd->MAX_ORDER], sumint);
}

void double_sct_i(iteration_data *itd, double *f, double *g, double *h, long i, int ORDER_INDEX)
{
	long vi,j;
	double sumint, sumext, partial, zero;
	double_init(&zero);
	double_init(&partial);
	double_init(&sumint);
	double_init(&sumext);
	double_set(&sumext, zero);
	for(j = 1; j <= ORDER_INDEX; j++) {
		double_set(&sumint, zero);
		for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
			double_mul(&partial, g[itd->PREV_VI[vi]*itd->MAX_ORDER+ORDER_INDEX-j], f[itd->PREV_IV[vi]*itd->MAX_ORDER+j]);
			double_mul_i(&partial, partial, itd->PREV_COEF[vi]);		
			double_add(&sumint, sumint, partial);
		}
		double_mul_i(&sumint, sumint, j);		
		double_add(&sumext, sumext, sumint);
	}
	double_div_i(&sumext, sumext, ORDER_INDEX);		
	double_set(&h[i*itd->MAX_ORDER+ORDER_INDEX], sumext);
}

void    double_sin_cos_t (iteration_data *itd, double *f, double *s, double *c, int ORDER_INDEX) {
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_sin (&s[0], f[0]);
					double_cos (&c[0], f[0]);
				} else {
					double_sct_0 (itd, f,c,s,i);
					double_sct_0 (itd, f,s,c,i);
					double_mul_i (&c[i*itd->MAX_ORDER], c[i*itd->MAX_ORDER], -1);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_sct_i(itd, f,s,c,i,ORDER_INDEX);
				double_sct_i(itd, f,c,s,i,ORDER_INDEX);
				double_mul_i (&c[i*itd->MAX_ORDER+ORDER_INDEX],
					 c[i*itd->MAX_ORDER+ORDER_INDEX], -1);
			}
		}
	}
}



void double_sin_t(iteration_data *itd, double *s, double *c, double *w, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_sin(&w[0], s[0]);
				} else {
					double_sct_0(itd, s,c,w,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_sct_i(itd, s,c,w,i,ORDER_INDEX);
			}
		}
	}
}

void double_cos_t(iteration_data *itd, double *c, double *s, double *w, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_cos(&w[0], c[0]);
				} else {
					double_sct_0(itd, c,s,w,i);
					double_mul_i(&w[i*itd->MAX_ORDER], w[i*itd->MAX_ORDER], -1);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_sct_i(itd, c,s,w,i,ORDER_INDEX);
				double_mul_i(&w[i*itd->MAX_ORDER+ORDER_INDEX], w[i*itd->MAX_ORDER+ORDER_INDEX], -1);
			}
		}
	}
}



void    double_sinh_cosh_t (iteration_data *itd, double *f, double *s, double *c, int ORDER_INDEX) {
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_sinh (&s[0], f[0]);
					double_cosh (&c[0], f[0]);
				} else {
					double_sct_0 (itd, f,c,s,i);
					double_sct_0 (itd, f,s,c,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_sct_i(itd, f,s,c,i,ORDER_INDEX);
				double_sct_i(itd, f,c,s,i,ORDER_INDEX);
			}
		}
	}
}




void double_sinh_t(iteration_data *itd, double *s, double *c, double *w, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_sinh(&w[0], s[0]);
				} else {
					double_sct_0(itd, s,c,w,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_sct_i(itd, s,c,w,i,ORDER_INDEX);
			}
		}
	}
}
void double_cosh_t(iteration_data *itd, double *c, double *s, double *w, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_cosh(&w[0], c[0]);
				} else {
					double_sct_0(itd, c,s,w,i);
					double_mul_i(&w[i*itd->MAX_ORDER], w[i*itd->MAX_ORDER], -1);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_sct_i(itd, c,s,w,i,ORDER_INDEX);
			}
		}
	}
}

/***************************************************************/

void double_fgt_0(iteration_data *itd, double *f, double *g, double *h, long i)
{
	long vi;
	double sumint, partial, zero;
	double_init(&zero);
	double_init(&partial);
	double_init(&sumint);
	double_set(&sumint, zero);					
	if(double_equal (g[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********   double_fgt_0   ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		double_set_d (&g[0], 1.);
	}
	
	for( vi = itd->PREVSTAR_ACCUM[i]; vi < itd->PREVSTAR_ACCUM[i+1]; vi++) {
		if(itd->PREVSTAR_VI[vi]>0) {
			double_mul(&partial, h[itd->PREVSTAR_IV[vi]*itd->MAX_ORDER], g[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER]);
			double_mul_i(&partial, partial, itd->PREVSTAR_COEF[vi]);
			double_add(&sumint, sumint, partial);
		}
	}
	double_sub(&sumint, f[i*itd->MAX_ORDER], sumint);
	double_div(&sumint, sumint, g[0]);
	double_set(&h[i*itd->MAX_ORDER], sumint);
}
void double_fgt_i(iteration_data *itd, double *f, double *g, double *h, long i, int ORDER_INDEX)
{
	long vi,j;
	double sumint, sumext, partial, zero, ht;
	double_init(&zero);
	double_init(&partial);
	double_init(&sumint);
	double_init(&sumext);
	double_init(&ht);
	double_set(&sumext, zero);
	if(double_equal (g[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********   double_fgt_i   ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		double_set_d (&g[0], 1.);
	}
	for(j = 0; j < ORDER_INDEX; j++) {
		double_set(&sumint, zero);
		for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
			double_htilde(itd, h, ORDER_INDEX-j, itd->PREV_IV[vi], i, &ht, ORDER_INDEX);
			double_mul(&partial, ht, g[itd->PREV_VI[vi]*itd->MAX_ORDER+j]);
			double_mul_i(&partial, partial, itd->PREV_COEF[vi]);
			double_add(&sumint, sumint, partial);
		}
		double_mul_i(&sumint, sumint, ORDER_INDEX-j);
		double_add(&sumext, sumext, sumint);
	}
	double_div_i(&sumext, sumext, ORDER_INDEX);
	double_sub(&sumext, f[i*itd->MAX_ORDER+ORDER_INDEX], sumext);
	double_div(&sumext, sumext, g[0]);
	double_set(&h[i*itd->MAX_ORDER+ORDER_INDEX], sumext);
}


void	double_log_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_log(&w[0], u[0]);
				} else {
					double_fgt_0(itd,u,u,w,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_fgt_i(itd, u,u,w,i,ORDER_INDEX);
			}
		}
	}
}
void	double_asin_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_asin(&h[0], f[0]);
				} else {
					double_fgt_0(itd,f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_fgt_i(itd, f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	double_acos_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_acos(&h[0], f[0]);
				} else {
					double_fgt_0(itd, f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_fgt_i(itd, f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	double_atan_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_atan(&h[0], f[0]);
				} else {
					double_fgt_0(itd, f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_fgt_i(itd, f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	double_asinh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_asinh(&h[0], f[0]);
				} else {
					double_fgt_0(itd, f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_fgt_i(itd, f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	double_acosh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_acosh(&h[0], f[0]);
				} else {
					double_fgt_0(itd, f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_fgt_i(itd, f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	double_atanh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					double_atanh(&h[0], f[0]);
				} else {
					double_fgt_0(itd, f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				double_fgt_i(itd, f,g,h,i,ORDER_INDEX);
			}
		}
	}
}

