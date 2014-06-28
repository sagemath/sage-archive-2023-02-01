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

#include "commonITER.h"
#include "mpfrITER.h"

extern int	printError;


void	varMP_init(iteration_data *itd, mpfr_t *var, mpfr_t *v, mpfr_t t)
{
	int i, j, k;
	for (i=0; i<=itd->NVARS; i++) 
		for (j=0; j<itd->NDER; j++)
			for (k=0; k<itd->MAX_ORDER; k++)
				mpfrts_init (var + i*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k);
	mpfrts_set (var, t);
	mpfrts_set_i (var+1, 1);
	for(i=1; i<=itd->NVARS; i++) 
		mpfrts_set (var + i*itd->LINK_ELEMENTS, v[i-1]);
}
void parMP_init(iteration_data *itd, mpfr_t *par, mpfr_t *p)
{
	int i, j, k;
	for (i=0; i<itd->NPARS; i++) 
		for (j=0; j<itd->NDER; j++)
			for (k=0; k<itd->MAX_ORDER; k++)
				mpfrts_init (par + i*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k);
	for(i=0; i<itd->NPARS; i++) mpfrts_set (par + i*itd->LINK_ELEMENTS, p[i]);
	for(i=0; i <itd->NPARTIALS; i++)
		if(!is_variable(itd, itd->PARTIAL_LIST[i]))
			mpfrts_set_i (par + (itd->PARTIAL_LIST[i]-itd->NVARS-1)*itd->LINK_ELEMENTS +(i+1)*itd->MAX_ORDER, 1);	
}
void linkMP_init(iteration_data *itd, mpfr_t *lk)
{
	int i, j, k;
	for (i=0; i<itd->NLINKS; i++) 
		for (j=0; j<itd->NDER; j++)
			for (k=0; k<itd->MAX_ORDER; k++)
				mpfrts_init (lk + i*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k);
}


void derMP_init(iteration_data *itd, mpfr_t *var, mpfr_t *par, mpfr_t *v)
{
	int i;
	if (itd->clearPartials)
	for(i=0; i <itd->NPARTIALS; i++) {
		if(is_variable(itd, itd->PARTIAL_LIST[i])) 
			mpfrts_set_i (var + itd->PARTIAL_LIST[i]*itd->LINK_ELEMENTS +(i+1)*itd->MAX_ORDER , 1);
	}
	else {
		int j, nvf;
		nvf = itd->NVARS+itd->NFUNS;
		for( i = 0; i < itd->NVARS; i++) 
			for(j = 1; j < itd->NDER; j++) 
				mpfrts_set (var + (i+1)*itd->LINK_ELEMENTS +j*itd->MAX_ORDER, v[i+(j*nvf)]);
	}

	itd->clearPartials = 0;
}


void clear_mpfr_vpl (iteration_data *itd,  mpfr_t *var, mpfr_t *par, mpfr_t *link) 
{

	int i, j, k;
	for (i=0; i<itd->NVARS+1; i++)
		for (j=0; j<itd->NDER; j++)
			for (k=0; k<itd->MAX_ORDER; k++)
				mpfr_clear (*(var + i*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k));
	for (i=0; i<itd->NPARS; i++)
		for (j=0; j<itd->NDER; j++)
			for (k=0; k<itd->MAX_ORDER; k++)
				mpfr_clear (*(par + i*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k));
	for (i=0; i<itd->NLINKS; i++)
		for (j=0; j<itd->NDER; j++)
			for (k=0; k<itd->MAX_ORDER; k++)
				mpfr_clear (*(link + i*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k));
	mpfr_free_cache ();

}


void	write_sol_MP(iteration_data *itd, mpfr_t *cvfd, mpfr_t *var, mpfr_t *link)
{
	int i,j,k, nvf;
	nvf = itd->NVARS+itd->NFUNS;
	mpfr_t vfun;
	mpfrts_init(&vfun);
	for( i = 0; i < itd->NVARS; i++) 
		for(j = 0; j < itd->NDER; j++) 
			for(k = 0; k < itd->MAX_ORDER; k++)
				mpfrts_set (&cvfd[(i+(j*nvf))*itd->MAX_ORDER + k], *(var + (i+1)*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k));
	for( i = 0; i < itd->NFUNS; i++) 
		for(j = 0; j < itd->NDER; j++) 
			for(k = 0; k < itd->MAX_ORDER; k++) {
				if (itd->FUNCTION_LIST[i] >= 0) mpfrts_set(&vfun, *(link + itd->FUNCTION_LIST[i]*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k));
				else mpfrts_set(&vfun, *(var -itd->FUNCTION_LIST[i]*itd->LINK_ELEMENTS +j*itd->MAX_ORDER +k));
				mpfrts_set (&cvfd[(itd->NVARS+i+(j*nvf))*itd->MAX_ORDER +  k], vfun);
			}
	mpfr_clear (vfun);
	mpfr_free_cache ();
}

/************************************************************************/

void mpfrts_htilde(iteration_data *itd, mpfr_t *h, long j, long v, long i, mpfr_t *ht, int ORDER_INDEX)
{
	mpfr_t cero;
	mpfrts_init (&cero);
	if(ORDER_INDEX >= 0) {
		if( j == ORDER_INDEX && i==v) mpfrts_set(ht, cero);
		else mpfrts_set(ht, h[v*itd->MAX_ORDER +j]);
	}
	mpfr_clear (cero);
	mpfr_free_cache ();
}

/************************************************************************/

void	mpfrts_var_t(iteration_data *itd, mpfr_t *f, mpfr_t *w, int ORDER_INDEX)
{
	mpfr_t val; mpfrts_init (&val);
	int i;
	if(ORDER_INDEX > 0) {
		for(i = 0; i < itd->NDER; i++) {
			mpfrts_div_i(&val, f[i*itd->MAX_ORDER +ORDER_INDEX-1], ORDER_INDEX);
			mpfrts_set(&w[i*itd->MAX_ORDER +ORDER_INDEX], val);
		}
	}
	mpfr_clear (val);
	mpfr_free_cache ();

}

void	mpfrts_var_t_c(iteration_data *itd, char* cs, mpfr_t *w, int ORDER_INDEX)
{
	if(ORDER_INDEX == 1) {
		mpfrts_set_str(&w[1], cs);
	}
}
void	mpfrts_var_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *w, int ORDER_INDEX)
{
	if(ORDER_INDEX == 1) {
		mpfrts_set(&w[1], c);
	}
}


void	mpfrts_add_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < itd->NDER; i++) 
			mpfrts_add(&w[i*itd->MAX_ORDER +ORDER_INDEX], u[i*itd->MAX_ORDER +ORDER_INDEX], v[i*itd->MAX_ORDER +ORDER_INDEX]);
	}
}

void	mpfrts_sub_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < itd->NDER; i++) 
			mpfrts_sub(&w[i*itd->MAX_ORDER +ORDER_INDEX], u[i*itd->MAX_ORDER +ORDER_INDEX], v[i*itd->MAX_ORDER +ORDER_INDEX]);
	}
}

void	mpfrts_add_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX)
{
	long i;
	mpfr_t c; mpfrts_init (&c);
	if(ORDER_INDEX >= 0) {
		mpfrts_set_str(&c, cs);
		if(ORDER_INDEX == 0) 
			mpfrts_add(&w[ORDER_INDEX], c, u[ORDER_INDEX]);
		else
			mpfrts_set(&w[ORDER_INDEX], u[ORDER_INDEX]);
		for(i = 1; i < itd->NDER; i++) 
			mpfrts_set(&w[i*itd->MAX_ORDER +ORDER_INDEX], u[i*itd->MAX_ORDER +ORDER_INDEX]);
	}
	mpfr_clear (c);
	mpfr_free_cache ();
}
void	mpfrts_add_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) 
			mpfrts_add(&w[ORDER_INDEX], c, u[ORDER_INDEX]);
		else
			mpfrts_set(&w[ORDER_INDEX], u[ORDER_INDEX]);
		for(i = 1; i < itd->NDER; i++) 
			mpfrts_set(&w[i*itd->MAX_ORDER +ORDER_INDEX], u[i*itd->MAX_ORDER +ORDER_INDEX]);
	}
	mpfr_free_cache ();	
}

void	mpfrts_sub_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX)
{
	long i;
	mpfr_t c; mpfrts_init (&c);
	if(ORDER_INDEX >= 0) {
		mpfrts_set_str(&c, cs);
		if(ORDER_INDEX == 0) 
			mpfrts_sub(&w[ORDER_INDEX], c, u[ORDER_INDEX]);
		else
			mpfrts_mul_i(&w[ORDER_INDEX], u[ORDER_INDEX],-1);
		for(i = 1; i < itd->NDER; i++) 
			mpfrts_mul_i(&w[i*itd->MAX_ORDER +ORDER_INDEX], u[i*itd->MAX_ORDER +ORDER_INDEX],-1);
	}
	mpfr_clear (c);
	mpfr_free_cache ();
}
void	mpfrts_sub_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) 
			mpfrts_sub(&w[ORDER_INDEX], c, u[ORDER_INDEX]);
		else
			mpfrts_mul_i(&w[ORDER_INDEX], u[ORDER_INDEX],-1);
		for(i = 1; i < itd->NDER; i++) 
			mpfrts_mul_i(&w[i*itd->MAX_ORDER +ORDER_INDEX], u[i*itd->MAX_ORDER +ORDER_INDEX],-1);
	}
	mpfr_free_cache ();
}


void	mpfrts_mul_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX)
{
	long i, vi, j;
	mpfr_t sumint, sumext, partial, zero;
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&zero);

		for(i = 0; i < itd->NDER; i++) {
			mpfrts_set(&sumext, zero);
			for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
				mpfrts_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {

					mpfrts_mul(&partial, u[itd->PREV_VI[vi]*itd->MAX_ORDER +ORDER_INDEX-j], v[itd->PREV_IV[vi]*itd->MAX_ORDER +j]);
					mpfrts_add(&sumint, sumint, partial);
				}
				mpfrts_mul_i(&sumint, sumint, itd->PREV_COEF[vi]);
				mpfrts_add(&sumext, sumext, sumint);
			}
			mpfrts_set(&w[i*itd->MAX_ORDER +ORDER_INDEX], sumext);
		}
	}
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (partial);
	mpfr_clear (zero);
	mpfr_free_cache ();

}


void	mpfrts_mul_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX)
{
	long i;
	mpfr_t c; mpfrts_init (&c);
	if(ORDER_INDEX >= 0) {
		mpfrts_set_str(&c, cs);
		for(i = 0; i < itd->NDER; i++) 
			mpfrts_mul(&w[i*itd->MAX_ORDER +ORDER_INDEX], c, u[i*itd->MAX_ORDER +ORDER_INDEX]);
	}
	mpfr_clear (c);
	mpfr_free_cache ();

}
void	mpfrts_mul_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX)
{
	long i;
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < itd->NDER; i++) 
			mpfrts_mul(&w[i*itd->MAX_ORDER +ORDER_INDEX], c, u[i*itd->MAX_ORDER +ORDER_INDEX]);
	}
	mpfr_free_cache ();	
}

void mpfrts_div_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint, sumext, partial, zero, ht;
	if(g[0] == 0 ) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********   mpfrts_div_t   ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		mpfrts_set_str (&g[0], "1");
	}
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&zero);
		mpfrts_init(&ht);
		for(i = 0; i < itd->NDER; i++) {
			mpfrts_set(&sumext, zero);
			for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
				mpfrts_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					mpfrts_htilde(itd, h, ORDER_INDEX-j, itd->PREV_VI[vi], i, &ht, ORDER_INDEX);
					mpfrts_mul(&partial, ht, g[itd->PREV_IV[vi]*itd->MAX_ORDER +j]);
					mpfrts_add(&sumint, sumint, partial);
				}
				mpfrts_mul_i(&sumint, sumint, itd->PREV_COEF[vi]);
				mpfrts_add(&sumext, sumext, sumint);
			}
			mpfrts_sub(&sumext, f[i*itd->MAX_ORDER +ORDER_INDEX], sumext);
			mpfrts_div(&sumext, sumext, g[0]);
			mpfrts_set(&h[i*itd->MAX_ORDER +ORDER_INDEX], sumext);
		}
		mpfr_clear (sumint);
		mpfr_clear (sumext);
		mpfr_clear (partial);
		mpfr_clear (zero);
		mpfr_clear (ht);
	}
	
}

void	mpfrts_div_t_vc(iteration_data *itd, mpfr_t *u, mpfr_t c, mpfr_t *w, int ORDER_INDEX)
{
	long i;
	mpfr_t ci;
	mpfrts_init(&ci);
	mpfrts_i_div(&ci, 1, c);
	if(ORDER_INDEX >= 0) {
		for(i = 0; i < itd->NDER; i++) 
			mpfrts_mul(&w[i*itd->MAX_ORDER +ORDER_INDEX], ci, u[i*itd->MAX_ORDER +ORDER_INDEX]);
	}
	mpfr_free_cache ();	
}

void	mpfrts_div_t_cv(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint, sumext, partial, zero, f, wt;
	mpfrts_init (&zero);
	if(mpfrts_equal (u[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("*********   mpfrts_div_t_cv  **********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		mpfrts_set_d (&u[0], 1.);
	}
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&f);
		mpfrts_init(&wt);
		for(i = 0; i < itd->NDER; i++) {
			if(ORDER_INDEX == 0 && i == 0 ) mpfrts_set_d(&f, 1.);
			else mpfrts_set(&f, zero);
			mpfrts_set(&sumext, zero);
			for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
				mpfrts_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					mpfrts_htilde(itd, w, ORDER_INDEX-j, itd->PREV_VI[vi], i, &wt, ORDER_INDEX);
					mpfrts_div(&wt, wt, c);
					mpfrts_mul(&partial, wt, u[itd->PREV_IV[vi]*itd->MAX_ORDER +j]);
					mpfrts_add(&sumint, sumint, partial);
				}
				mpfrts_mul_i(&sumint, sumint, itd->PREV_COEF[vi]);
				mpfrts_add(&sumext, sumext, sumint);
			}
			mpfrts_sub(&sumext, f, sumext);
			mpfrts_div(&sumext, sumext, u[0]);
			mpfrts_mul(&sumext, sumext, c);
			mpfrts_set(&w[i*itd->MAX_ORDER +ORDER_INDEX], sumext);
		}
	}
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (partial);
	mpfr_clear (zero);
	mpfr_clear (f);
	mpfr_clear (wt);
	mpfr_free_cache ();
	
}


void mpfrts_inv_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint, sumext, partial, zero, f, wt;
	mpfrts_init (&zero);
	if(mpfrts_equal (u[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********  mpfrts_inv_t    ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		mpfrts_set_d (&u[0], 1.);
	}
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&f);
		mpfrts_init(&wt);
		for(i = 0; i < itd->NDER; i++) {
			if(ORDER_INDEX == 0 && i == 0 ) mpfrts_set_d(&f, 1.);
			else mpfrts_set(&f, zero);
			mpfrts_set(&sumext, zero);
			for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
				mpfrts_set(&sumint, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					mpfrts_htilde(itd, w, ORDER_INDEX-j, itd->PREV_VI[vi], i, &wt, ORDER_INDEX);
					mpfrts_mul(&partial, wt, u[itd->PREV_IV[vi]*itd->MAX_ORDER +j]);
					mpfrts_add(&sumint, sumint, partial);
				}
				mpfrts_mul_i(&sumint, sumint, itd->PREV_COEF[vi]);
				mpfrts_add(&sumext, sumext, sumint);
			}
			mpfrts_sub(&sumext, f, sumext);
			mpfrts_div(&sumext, sumext, u[0]);
			mpfrts_set(&w[i*itd->MAX_ORDER +ORDER_INDEX], sumext);
		}
	}
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (partial);
	mpfr_clear (zero);
	mpfr_clear (f);
	mpfr_clear (wt);
	mpfr_free_cache ();

}

void mpfrts_exp_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint,sumext, partial, zero;
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&zero);
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_exp(&sumext, u[0]);
					mpfrts_set(&w[0], sumext);
				} else {
					mpfrts_set(&sumint, zero);					
					for( vi = itd->PREVSTAR_ACCUM[i]; vi < itd->PREVSTAR_ACCUM[i+1]; vi++) {
						mpfrts_mul(&partial, w[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER ], u[itd->PREVSTAR_IV[vi]*itd->MAX_ORDER ]);
						mpfrts_mul_i(&partial, partial, itd->PREVSTAR_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_set(&w[i*itd->MAX_ORDER ], sumint);
				}
			}
		} else {		
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_set(&sumext, zero);
				for(j = 0; j < ORDER_INDEX; j++) {
					mpfrts_set(&sumint, zero);
					for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
						mpfrts_mul(&partial, w[itd->PREV_VI[vi]*itd->MAX_ORDER +j], u[itd->PREV_IV[vi]*itd->MAX_ORDER +ORDER_INDEX-j]);
						mpfrts_mul_i(&partial, partial, itd->PREV_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_mul_i(&sumint, sumint, ORDER_INDEX-j);
					mpfrts_add(&sumext, sumext, sumint);
				}
				mpfrts_div_i(&sumext, sumext, ORDER_INDEX);
				mpfrts_set(&w[i*itd->MAX_ORDER +ORDER_INDEX], sumext);
			}
			
		}
	}
}


void mpfrts_pow_t_c(iteration_data *itd, mpfr_t *u, char* cs, mpfr_t *w, int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint, sumext, partial, partialb, zero, wt, c;
	mpfrts_init (&zero);
	if(mpfrts_equal (u[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********  mpfrts_pow_t_c  ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		mpfrts_set_d (&u[0], 1.);
	}
	
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&partialb);
		mpfrts_init(&wt);
		mpfrts_init(&c);
		mpfrts_set_str(&c, cs);
		mpfrts_set (&sumext, zero);
	
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_pow(&sumext, u[0], c);
					mpfrts_set(&w[0], sumext);
				} else {
					mpfrts_set(&sumint, zero);					
					for( vi = itd->PREVSTAR_ACCUM[i]; vi < itd->PREVSTAR_ACCUM[i+1]; vi++) {
						mpfrts_mul(&partial, w[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER ], u[itd->PREVSTAR_IV[vi]*itd->MAX_ORDER ]);
						mpfrts_mul(&partial, partial, c);
						mpfrts_htilde(itd, w, 0, itd->PREVSTAR_IV[vi], i, &wt, 0);
						mpfrts_mul(&partialb, wt, u[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER ]);
						mpfrts_sub(&partial, partial, partialb);
						mpfrts_mul_i(&partial, partial, itd->PREVSTAR_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_div(&sumint, sumint, u[0]);
					mpfrts_set(&w[i*itd->MAX_ORDER ], sumint);
				}
			}
		} else {		
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_set(&sumext, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					mpfrts_set(&sumint, zero);
					for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
						mpfrts_htilde(itd, w, j, itd->PREV_VI[vi], i, &wt, ORDER_INDEX);
						mpfrts_mul(&partial, wt, u[itd->PREV_IV[vi]*itd->MAX_ORDER +ORDER_INDEX-j]);
						mpfrts_mul_i(&partial, partial, itd->PREV_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_add_i(&partial, c, 1);
					mpfrts_mul_i(&partial, partial, j);
					mpfrts_mul_i(&partialb, c, ORDER_INDEX);
					mpfrts_sub(&partial, partialb, partial);
					mpfrts_mul(&sumint, sumint, partial);
					mpfrts_add(&sumext, sumext, sumint);
				}
				mpfrts_div(&sumext, sumext, u[0]);
				mpfrts_div_i(&sumext, sumext, ORDER_INDEX);
				mpfrts_set(&w[i*itd->MAX_ORDER +ORDER_INDEX], sumext);
			}
			
		}
		
	}
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (partial);
	mpfr_clear (partialb);
	mpfr_clear (zero);
	mpfr_clear (c);
	mpfr_clear (wt);
	mpfr_free_cache ();

	
	
}

void mpfrts_pow_t_cc(iteration_data *itd, mpfr_t *u, mpfr_t c, mpfr_t *w, int ORDER_INDEX)
{
	long i, j, vi;
	mpfr_t sumint, sumext, partial, partialb, zero, wt;
	mpfrts_init (&zero);
	if(mpfrts_equal (u[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("*********  mpfrts_pow_t_cc  ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		mpfrts_set_d (&u[0], 1.);
	}
	
	if(ORDER_INDEX >= 0) {
		mpfrts_init(&sumint);
		mpfrts_init(&sumext);
		mpfrts_init(&partial);
		mpfrts_init(&partialb);
		mpfrts_init(&wt);
		mpfrts_set (&sumext, zero);
		
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_pow(&sumext, u[0], c);
					mpfrts_set(&w[0], sumext);
				} else {
					mpfrts_set(&sumint, zero);					
					for( vi = itd->PREVSTAR_ACCUM[i]; vi < itd->PREVSTAR_ACCUM[i+1]; vi++) {
						mpfrts_mul(&partial, w[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER ], u[itd->PREVSTAR_IV[vi]*itd->MAX_ORDER ]);
						mpfrts_mul(&partial, partial, c);
						mpfrts_htilde(itd, w, 0, itd->PREVSTAR_IV[vi], i, &wt, 0);
						mpfrts_mul(&partialb, wt, u[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER ]);
						mpfrts_sub(&partial, partial, partialb);
						mpfrts_mul_i(&partial, partial, itd->PREVSTAR_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_div(&sumint, sumint, u[0]);
					mpfrts_set(&w[i*itd->MAX_ORDER ], sumint);
				}
			}
		} else {		
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_set(&sumext, zero);
				for(j = 0; j <= ORDER_INDEX; j++) {
					mpfrts_set(&sumint, zero);
					for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
						mpfrts_htilde(itd, w, j, itd->PREV_VI[vi], i, &wt, ORDER_INDEX);
						mpfrts_mul(&partial, wt, u[itd->PREV_IV[vi]*itd->MAX_ORDER +ORDER_INDEX-j]);
						mpfrts_mul_i(&partial, partial, itd->PREV_COEF[vi]);
						mpfrts_add(&sumint, sumint, partial);
					}
					mpfrts_add_i(&partial, c, 1);
					mpfrts_mul_i(&partial, partial, j);
					mpfrts_mul_i(&partialb, c, ORDER_INDEX);
					mpfrts_sub(&partial, partialb, partial);
					mpfrts_mul(&sumint, sumint, partial);
					mpfrts_add(&sumext, sumext, sumint);
				}
				mpfrts_div(&sumext, sumext, u[0]);
				mpfrts_div_i(&sumext, sumext, ORDER_INDEX);
				mpfrts_set(&w[i*itd->MAX_ORDER +ORDER_INDEX], sumext);
			}
			
		}
		
	}
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (partial);
	mpfr_clear (partialb);
	mpfr_clear (zero);
	mpfr_clear (wt);
	mpfr_free_cache ();
}

/************************************************************************/

void mpfrts_sct_0(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i)
{
	long vi;
	mpfr_t sumint, partial, zero;
	mpfrts_init(&zero);
	mpfrts_init(&partial);
	mpfrts_init(&sumint);
	mpfrts_set(&sumint, zero);					
	for( vi = itd->PREVSTAR_ACCUM[i]; vi < itd->PREVSTAR_ACCUM[i+1]; vi++) {
		mpfrts_mul(&partial, g[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER ], f[itd->PREVSTAR_IV[vi]*itd->MAX_ORDER ]);
		mpfrts_mul_i(&partial, partial, itd->PREVSTAR_COEF[vi]);		
		mpfrts_add(&sumint, sumint, partial);
	}
	mpfrts_set(&h[i*itd->MAX_ORDER ], sumint);
	mpfr_clear (sumint);
	mpfr_clear (partial);
	mpfr_clear (zero);
	mpfr_free_cache ();
}

void mpfrts_sct_i(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i, int ORDER_INDEX)
{
	long vi,j;
	mpfr_t sumint, sumext, partial, zero;
	mpfrts_init(&zero);
	mpfrts_init(&partial);
	mpfrts_init(&sumint);
	mpfrts_init(&sumext);
	mpfrts_set(&sumext, zero);
	for(j = 1; j <= ORDER_INDEX; j++) {
		mpfrts_set(&sumint, zero);
		for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
			mpfrts_mul(&partial, g[itd->PREV_VI[vi]*itd->MAX_ORDER +ORDER_INDEX-j], f[itd->PREV_IV[vi]*itd->MAX_ORDER +j]);
			mpfrts_mul_i(&partial, partial, itd->PREV_COEF[vi]);		
			mpfrts_add(&sumint, sumint, partial);
		}
		mpfrts_mul_i(&sumint, sumint, j);		
		mpfrts_add(&sumext, sumext, sumint);
	}
	mpfrts_div_i(&sumext, sumext, ORDER_INDEX);		
	mpfrts_set(&h[i*itd->MAX_ORDER +ORDER_INDEX], sumext);
	mpfr_clear (partial);
	mpfr_clear (sumint);
	mpfr_clear (zero);
	mpfr_clear (sumext);
}

void mpfrts_sin_cos_t (iteration_data *itd, mpfr_t *f, mpfr_t *s, mpfr_t *c, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_sin (&s[0], f[0]);
					mpfrts_cos (&c[0], f[0]);
				} else {
					mpfrts_sct_0 (itd,f,c,s,i);
					mpfrts_sct_0 (itd,f,s,c,i);
						mpfrts_mul_i (&c[i*itd->MAX_ORDER], c[i*itd->MAX_ORDER], -1);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_sct_i(itd,f,s,c,i,ORDER_INDEX);
				mpfrts_sct_i(itd,f,c,s,i,ORDER_INDEX);
				mpfrts_mul_i (&c[i*itd->MAX_ORDER +ORDER_INDEX],
					 c[i*itd->MAX_ORDER +ORDER_INDEX], -1);
			}
		}
	}
}


void mpfrts_sin_t(iteration_data *itd, mpfr_t *s, mpfr_t *c, mpfr_t *w, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_sin(&w[0], s[0]);
				} else {
					
					mpfrts_sct_0(itd,s,c,w,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_sct_i(itd,s,c,w,i,ORDER_INDEX);
			}
		}
	}
}

void mpfrts_cos_t(iteration_data *itd, mpfr_t *c, mpfr_t *s, mpfr_t *w, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_cos(&w[0], c[0]);
				} else {
					mpfrts_sct_0(itd,c,s,w,i);
					mpfrts_mul_i(&w[i*itd->MAX_ORDER ], w[i*itd->MAX_ORDER ], -1);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_sct_i(itd,c,s,w,i,ORDER_INDEX);
				mpfrts_mul_i(&w[i*itd->MAX_ORDER +ORDER_INDEX], w[i*itd->MAX_ORDER +ORDER_INDEX], -1);
			}
		}
	}
}


void    mpfrts_sinh_cosh_t (iteration_data *itd, mpfr_t *f,  mpfr_t *s, mpfr_t *c, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_sinh (&s[0], f[0]);
					mpfrts_cosh (&c[0], f[0]);
				} else {
					mpfrts_sct_0 (itd,f,c,s,i);
					mpfrts_sct_0 (itd,f,s,c,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_sct_i(itd,f,s,c,i,ORDER_INDEX);
				mpfrts_sct_i(itd,f,c,s,i,ORDER_INDEX);
			}
		}
	}
}



void mpfrts_sinh_t(iteration_data *itd, mpfr_t *s, mpfr_t *c, mpfr_t *w, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_sinh(&w[0], s[0]);
				} else {
					mpfrts_sct_0(itd,s,c,w,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_sct_i(itd,s,c,w,i,ORDER_INDEX);
			}
		}
	}
}
void mpfrts_cosh_t(iteration_data *itd, mpfr_t *c, mpfr_t *s, mpfr_t *w, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_cosh(&w[0], c[0]);
				} else {
					mpfrts_sct_0(itd,c,s,w,i);
					mpfrts_mul_i(&w[i*itd->MAX_ORDER ], w[i*itd->MAX_ORDER ], -1);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_sct_i(itd,c,s,w,i,ORDER_INDEX);
			}
		}
	}
}

/***************************************************************/

void mpfrts_fgt_0(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i)
{
	long vi;
	mpfr_t sumint, partial, zero;
	mpfrts_init(&zero);
	mpfrts_init(&partial);
	mpfrts_init(&sumint);
	mpfrts_set(&sumint, zero);					
	if(mpfrts_equal (g[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********   mpfrts_fgt_0   ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		mpfrts_set_d (&g[0], 1.);
	}
	
	for( vi = itd->PREVSTAR_ACCUM[i]; vi < itd->PREVSTAR_ACCUM[i+1]; vi++) {
		if(itd->PREVSTAR_VI[vi]>0) {
			mpfrts_mul(&partial, h[itd->PREVSTAR_IV[vi]*itd->MAX_ORDER ], g[itd->PREVSTAR_VI[vi]*itd->MAX_ORDER ]);
			mpfrts_mul_i(&partial, partial, itd->PREVSTAR_COEF[vi]);
			mpfrts_add(&sumint, sumint, partial);
		}
	}
	mpfrts_sub(&sumint, f[i*itd->MAX_ORDER ], sumint);
	mpfrts_div(&sumint, sumint, g[0]);
	mpfrts_set(&h[i*itd->MAX_ORDER], sumint);
	mpfr_clear (zero);
	mpfr_clear (partial);
	mpfr_clear (sumint);
}

void mpfrts_fgt_i(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i, int ORDER_INDEX)
{
	long vi,j;
	mpfr_t sumint, sumext, partial, zero, ht;
	mpfrts_init(&zero);
	mpfrts_init(&partial);
	mpfrts_init(&sumint);
	mpfrts_init(&sumext);
	mpfrts_init(&ht);
	mpfrts_set(&sumext, zero);
	if(mpfrts_equal (g[0], zero)) {
		if(printError) {
			printf("**********  Divide by cero  ***********\n");
			printf("**********   mpfrts_fgt_i   ***********\n");
			printf("**********   bad   result   ***********\n");
		}
		itd->errorITER = 1;
		mpfrts_set_d (&g[0], 1.);
	}
	for(j = 0; j < ORDER_INDEX; j++) {
		mpfrts_set(&sumint, zero);
		for( vi = itd->PREV_ACCUM[i]; vi < itd->PREV_ACCUM[i+1]; vi++) {
			mpfrts_htilde(itd, h, ORDER_INDEX-j, itd->PREV_IV[vi], i, &ht, ORDER_INDEX);
			mpfrts_mul(&partial, ht, g[itd->PREV_VI[vi]*itd->MAX_ORDER +j]);
			mpfrts_mul_i(&partial, partial, itd->PREV_COEF[vi]);
			mpfrts_add(&sumint, sumint, partial);
		}
		mpfrts_mul_i(&sumint, sumint, ORDER_INDEX-j);
		mpfrts_add(&sumext, sumext, sumint);
	}
	mpfrts_div_i(&sumext, sumext, ORDER_INDEX);
	mpfrts_sub(&sumext, f[i*itd->MAX_ORDER +ORDER_INDEX], sumext);
	mpfrts_div(&sumext, sumext, g[0]);
	mpfrts_set(&h[i*itd->MAX_ORDER +ORDER_INDEX], sumext);
	mpfr_clear (zero);
	mpfr_clear (partial);
	mpfr_clear (sumint);
	mpfr_clear (sumext);
	mpfr_clear (ht);
}


void	mpfrts_log_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_log(&w[0], u[0]);
				} else {
					mpfrts_fgt_0(itd,u,u,w,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_fgt_i(itd,u,u,w,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_asin_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_asin(&h[0], f[0]);
				} else {
					mpfrts_fgt_0(itd,f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_fgt_i(itd,f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_acos_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_acos(&h[0], f[0]);
				} else {
					mpfrts_fgt_0(itd,f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_fgt_i(itd,f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_atan_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_atan(&h[0], f[0]);
				} else {
					mpfrts_fgt_0(itd,f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_fgt_i(itd,f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_asinh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_asinh(&h[0], f[0]);
				} else {
					mpfrts_fgt_0(itd,f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_fgt_i(itd,f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_acosh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_acosh(&h[0], f[0]);
				} else {
					mpfrts_fgt_0(itd,f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_fgt_i(itd,f,g,h,i,ORDER_INDEX);
			}
		}
	}
}
void	mpfrts_atanh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX)
{
	long i; 
	if(ORDER_INDEX >= 0) {
		if(ORDER_INDEX == 0) {
			for(i = 0; i < itd->NDER; i++) {
				if(i == 0) {
					mpfrts_atanh(&h[0], f[0]);
				} else {
					mpfrts_fgt_0(itd,f,g,h,i);
				}
			}
		} else {
			for(i = 0; i < itd->NDER; i++) {
				mpfrts_fgt_i(itd,f,g,h,i,ORDER_INDEX);
			}
		}
	}
}

