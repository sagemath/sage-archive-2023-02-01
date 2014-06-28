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

#include "mpfrTODE.h"


extern int	_info_steps_taylor_;
mpfr_t	etapa_mp_minima;
mpfr_t	etapa_mp_maxima;
mpfr_t	etapa_mp_total;
extern int	num_etapas;
extern int	order_series;

extern int		order_estimator;
extern double 	fac1;
extern double	fac2;
extern double	fac3;
extern double	rmaxstep;
extern double	rminstep;
extern double	nitermax;
extern int		nordinc;
extern int		minord;
extern int		defect_error_control;
extern int		kahan_summation;

extern int		DECIMAL_PRECISION;

extern double	hant;

mpfr_t  *mp_derivative_initial_point = NULL;
mpfr_t  *mp_derivative_final_point   = NULL;

extern int	test_relminstep;
extern int	printError;
mpfr_t  mp_relminstep ;
extern unsigned long    max_num_steps;

void	mp_set_relminstep(mpfr_t val)
{
	test_relminstep = 1;
	mpfrts_set(&mp_relminstep ,val);
}



int getNsteps (void) {return num_etapas;}
int getOrder (void) {return order_series;}

	
void mp_set_info_taylor(void) {  
	_info_steps_taylor_ = 1;
	mpfrts_init (&etapa_mp_minima);
	mpfrts_init (&etapa_mp_maxima);
	mpfrts_init (&etapa_mp_total);

	mpfrts_set_str (&etapa_mp_minima, "1.e30");
	mpfrts_set_str (&etapa_mp_maxima, "0");
	mpfrts_set_str (&etapa_mp_total, "0.");
	num_etapas = 0; 
	order_series = 0;
}

void mp_unset_info_taylor(void) { 
	_info_steps_taylor_ = 0;
	mpfrts_init (&etapa_mp_minima);
	mpfrts_init (&etapa_mp_maxima);
	mpfrts_init (&etapa_mp_total);
	mpfrts_set_str (&etapa_mp_minima, "1.e30");
	mpfrts_set_str (&etapa_mp_maxima, "0");
	mpfrts_set_str (&etapa_mp_total, "0.");
	num_etapas = 0; 
	order_series = 0;
}

void mp_str_info_taylor(void) {
	double aux;
	mpfr_t tot, num, average;

	mpfrts_init (&tot); mpfrts_init (&num); mpfrts_init (&average);
	mpfrts_set (&tot, etapa_mp_total);
	aux = num_etapas;
	mpfrts_set_d (&num, aux);
	mpfrts_div (&average, tot, num);

	if(_info_steps_taylor_ == 1) {
		printf("============================================================\n");
		printf ("Number of integration steps in Taylor method: %d\n", num_etapas);
		mpfrts_write ("Minimum   integration step:", etapa_mp_minima);
		mpfrts_write ("Maximum   integration step:", etapa_mp_maxima);
		mpfrts_write ("Averaged  integration step:",average);
		printf("Order of the Taylor series  : %d\n",  order_series);
		printf("============================================================\n\n");
	}
	mpfr_clear (tot); 
	mpfr_clear (num); 
	mpfr_clear (average);
	mpfr_free_cache ();
}


void mp_add_info_step(mpfr_t tstep) {
	if (mpfrts_less (tstep, etapa_mp_minima)) 
		mpfrts_set (&etapa_mp_minima, tstep);
	if (mpfrts_greater (tstep, etapa_mp_maxima)) 
		mpfrts_set (&etapa_mp_maxima, tstep);
	num_etapas++;
	mpfrts_add (&etapa_mp_total, etapa_mp_total, tstep);
}


int mp_taylor_order(mpfr_t eps) {
	mpfr_t rop, coc;
	double aux, nrop;

	mpfrts_init (&rop);
	mpfrts_log (&rop, eps);
	mpfrts_init (&coc);
	mpfrts_set_str(&coc, "2.0");
	mpfrts_div (&rop, rop, coc);
	mpfrts_i_sub (&rop, nordinc, rop);
	aux = mpfrts_get_d (rop);
	mpfr_clear (rop);
	mpfr_clear (coc);
	mpfr_free_cache ();
	nrop = (int) ceil (aux);
	if(nrop < minord) nrop = minord;
	return nrop;
}

void mp_norm_inf_vec(mpfr_t *rop, int nvar, mpfr_t *y) 
{
	int i;
	mpfr_t ninf, ncoef;
	mpfrts_init (&ninf); 
	mpfrts_init (&ncoef);
	
	for (i=0; i<nvar; i++) {
		mpfrts_abs(&ncoef, y[i]);
		if (mpfrts_greater (ncoef, ninf))  mpfrts_set (&ninf, ncoef);
	}
	mpfrts_set(rop, ninf);
	
	mpfr_clear (ninf);
	mpfr_clear (ncoef);
	mpfr_free_cache ();
}


void  mp_norm_inf(mpfr_t *rop, int n, int k, mpfr_t *coef, int MAX_ORDER) 
{
	int i;
	mpfr_t ninf, ncoef;
	mpfrts_init (&ninf); 
	mpfrts_init (&ncoef);

	for (i=0; i<n; i++) {
		mpfrts_abs(&ncoef, coef[i*MAX_ORDER+k]);
		if (mpfrts_greater (ncoef, ninf))  mpfrts_set (&ninf, ncoef);
	}
	mpfrts_set(rop, ninf);

	mpfr_clear (ninf);
	mpfr_clear (ncoef);
	mpfr_free_cache ();
}

void mp_compute_step(mpfr_t *rop, double *hant, mpfr_t tol, int n, int ord, mpfr_t *coef, int MAX_ORDER) 
{
	mpfr_t tult, tpen, nor, zero, aux;
	int i;
	mpfrts_init (&zero); mpfrts_init (&aux);
	mpfrts_init (&tult); mpfrts_init (&tpen);
	mpfrts_init (&nor);
	mp_norm_inf (&nor, n, ord, coef, MAX_ORDER);
	if(mpfrts_equal(nor, zero)) mpfrts_set_d (&tult, 0.);
	else {
		mpfrts_set_d (&aux, 1./(ord+1));
		mpfrts_pow (&aux, tol, aux);
		mpfrts_set_d (&tult, -1./ord);
		mpfrts_pow (&tult, nor, tult);
		mpfrts_mul (&tult, aux, tult);	
	}

	mp_norm_inf (&nor, n, ord-1, coef, MAX_ORDER);
	if(mpfrts_equal(nor, zero)) mpfrts_set_d (&tpen, 0.);
	else {
		mpfrts_set_d (&aux, 1./(ord));
		mpfrts_pow (&aux, tol, aux);
		mpfrts_set_d (&tpen, -1./(ord-1));
		mpfrts_pow (&tpen, nor, tpen);
		mpfrts_mul (&tpen, aux, tpen);	
	} 

	if (mpfrts_equal (tult, zero) && (mpfrts_equal (tpen, zero))) {
		i = ord -1;
		mpfrts_set_d (rop, 0);

		while ((i > 0) && mpfrts_equal (nor, zero)) {
			mp_norm_inf (&nor, n, i-1, coef, MAX_ORDER);
			mpfrts_set_d (&aux, 1./i);
			mpfrts_pow (&aux, tol, aux);
			mpfrts_set_d (rop, -1./(i-1));
			mpfrts_pow (rop, tol, *rop);
			i--;
		}
	} else if (mpfrts_equal (tult, zero)) mpfrts_set (rop, tpen);
	else if (mpfrts_equal (tpen, zero)) mpfrts_set (rop, tult);
	else {
		if (mpfrts_less (tult, tpen)) mpfrts_set (rop, tult);
		else mpfrts_set (rop, tpen);
	}
	mpfr_t mp_fac; mpfrts_init (&mp_fac); mpfrts_set_d (&mp_fac, fac1);
	mpfrts_mul (rop, *rop, mp_fac);

	mpfr_t mp_hant; mpfrts_init (&mp_hant); mpfrts_set_d (&mp_hant, *hant);	
	mpfr_t mp_rmax; mpfrts_init (&mp_rmax); 
	mpfrts_set_d (&mp_rmax, rmaxstep);	
	mpfr_t mp_rmin; mpfrts_init (&mp_rmin); 
	mpfrts_set_d (&mp_rmin, rminstep);	

	if ((*hant) != 0.) {
		mpfrts_div (&aux, *rop, mp_hant);
		if (mpfrts_greater(aux, mp_rmax))
			mpfrts_mul (rop, mp_hant, mp_rmax);
		if (mpfrts_less (aux, mp_rmin))
			mpfrts_mul (rop, mp_hant, mp_rmin);

	}

	*hant = mpfrts_get_d (*rop);

	if(_info_steps_taylor_ == 1) mp_add_info_step(*rop);
	mpfr_clear (tult);
	mpfr_clear (tpen);
	mpfr_clear (nor);
	mpfr_clear (zero);
	mpfr_clear (aux);
	mpfr_clear (mp_hant);
	mpfr_clear (mp_rmax);
	mpfr_clear (mp_rmin);
	mpfr_clear (mp_fac);
	mpfr_free_cache ();	
}



void mp_compute_tol (mpfr_t *tol, mpfr_t tolrel, mpfr_t tolabs, int nvar, mpfr_t *y) 
{
	static mpfr_t norm_y1 ;
	static int first = 1;
	mpfr_t norm_y0, max_norm;
	
	if(first) {
		mpfrts_init (&norm_y1);
		first = 0;
	}
	mpfrts_init (&norm_y0);  
	mpfrts_init (&max_norm);


	mp_norm_inf_vec(&norm_y0, nvar, y);

	if (mpfrts_greater (norm_y0, norm_y1)) 
		mpfrts_set (&max_norm, norm_y0);
	else 
		mpfrts_set (&max_norm, norm_y1);


	mpfrts_mul (tol, tolrel, max_norm);
	mpfrts_add( tol, *tol, tolabs);
	
	mpfr_clear (norm_y0);
	mpfr_clear (max_norm);

	mpfr_free_cache ();

}


void mp_taylor_horner (int n, int ord, mpfr_t *coef, mpfr_t t, mpfr_t x[], int MAX_ORDER) 
{
	int i;
	for(i = 0; i < n; i++ )
		mp_horner(&coef[i*MAX_ORDER], ord, t, &x[i]);
}

void mp_taylor_horner_der(int n, int ord, mpfr_t *coef, mpfr_t t, mpfr_t x[], int MAX_ORDER) 
{
	int i;
	for(i = 0; i < n; i++ )
		mp_horner_der(&coef[i*MAX_ORDER], ord, t, &x[i]);
}


void mp_write_taylor_solution (int n, int j, mpfr_t tini, 
	mpfr_t x[], mpfr_t*** mat, FILE* fileout) 
{
	int i, precision;
	mpfr_t zero;
	mpfrts_init (&zero);
	precision = mpfrts_get_prec ();
	if (fileout != NULL) {
		mpfrts_fwrite(fileout,tini, precision);
		for (i=0; i<n; i++){ 
			mpfrts_fwrite(fileout, x[i], precision);
		}
		fprintf (fileout, "\n");
	}
	if (mat != NULL) {
		mpfrts_set (&(*mat)[j][0], tini);

		for (i=0; i<n; i++)
			mpfrts_set (&(*mat)[j][i+1], x[i]);
	}
	mpfr_clear (zero);
	mpfr_free_cache ();	
}




int mp_valid_step (MPLinkedFunction fcn,  iteration_data *itd, mpfr_t *step, mpfr_t tip, 
		mpfr_t tol, int nvar, int ncol, int order, 
		mpfr_t *cvfd, mpfr_t p[], int MAX_ORDER) {
	int i, j, iter = 1;
	int accepted = 0;
	mpfr_t b[nvar*order], y[nvar], yp[nvar], t, cn[ncol*MAX_ORDER],
		dif[nvar], nor, aux, up, fac;

	mpfrts_init (&t); 
	mpfrts_init (&nor); 
	mpfrts_init (&aux);
	mpfrts_init (&up); 
	mpfrts_init (&fac); 
	for (i=0; i<nvar; i++) {
		mpfrts_init (&y[i]);
		mpfrts_init (&dif[i]);
		mpfrts_init (&yp[i]);
		for (j=0; j<order; j++)
			mpfrts_init (&b[i*order+j]);
		for (j=0; j<=order; j++)
			mpfrts_init (&cn[i*MAX_ORDER+j]);
	}
	
	mpfrts_set_d (&fac, fac3);
	mpfrts_set_d (&up, fac2); 
	mpfrts_mul (&up, up, tol);
	
	while(iter <= nitermax) { 
//	

		for (i=0; i<nvar; i++) for (j=0; j<=order-1; j++)
			mpfrts_mul_i (&b[i*order+j], cvfd[i*MAX_ORDER+j+1], j+1);

		mpfrts_add (&t, tip, *step);
		mp_taylor_horner (nvar, order, cvfd, *step, y, MAX_ORDER);
		mp_taylor_horner (nvar, order-1, b, *step, yp, MAX_ORDER);
		fcn (itd, t, y, p, 1, cn);
		for (i=0; i<nvar; i++) {
			mpfrts_sub (&dif[i], yp[i], cn[i*MAX_ORDER+1]);
			mpfrts_abs (&aux, dif[i]);
			if (mpfrts_greater (aux, nor)) mpfrts_set (&nor, aux);
		}

		if (mpfrts_greater (nor, up)){
			mpfrts_mul (step, *step, fac);
		}
		
		if (mpfrts_greater (nor, up)){
			mpfrts_mul (step, *step, fac);
		}
		
//	
		iter++;
	}

	
	if (iter > nitermax) {
		printf("The maximum number of iterations in defect_error_control has been reached.\n");
		printf("The program continues with a stepsize that does not achieve the defect_error_control condition.\n");
	}
	
	mpfr_clear (t);
	mpfr_clear (nor);
	mpfr_clear (aux);
	mpfr_clear (up);
	mpfr_clear (fac);
	for (i=0; i<nvar; i++) {
		for (j=0; j<order; j++) 
			mpfr_clear (b[i*order+j]);
		mpfr_clear (y[i]);
		mpfr_clear (yp[i]);
		mpfr_clear (dif[i]);
	}
	for (i=0; i<ncol; i++) {
		for (j=0; j<=order; j++)
			mpfr_clear (cn[i*MAX_ORDER+j]);
	}

	mpfr_free_cache ();	

	return accepted;
}

/*****************************************************************************/
int  mp_tides_point(MPLinkedFunction fcn,  int *pdd,
					int nvar, int npar, int nfun, 
					mpfr_t x[], mpfr_t p[],
					mpfr_t t0, mpfr_t tf, mpfr_t dt, 	
					mpfr_t tolrel, mpfr_t tolabs,   
					mp_data_matrix *mat, FILE* fileout) 
{
	int err;
	mpfr_t aux;
	mpfrts_init(&aux);
	mpfrts_set(&aux, tf);
	mpfrts_sub(&aux, aux, t0);
	mpfrts_div(&aux, aux, dt);
	double auxd;
	auxd = mpfrts_get_d(aux);
	int ntot = (int) floor (auxd) ;
	err = mp_tides_delta(fcn,pdd,nvar,npar,nfun,x,p,t0,dt,ntot,tolrel,tolabs,mat,fileout);
	mpfr_clear (aux); 
	mpfr_free_cache ();	
	return err;
}	

int  mp_tides(MPLinkedFunction fcn,  int *pdd,
			  int nvar, int npar, int nfun, 
			  mpfr_t x[], mpfr_t p[],
			  mpfr_t lt[], int ntes, 	
			  mpfr_t tolrel, mpfr_t tolabs,   
			  mp_data_matrix *mat, FILE* fileout) 
{
	return mp_tides_list(fcn,pdd,nvar,npar,nfun,x,p,lt,ntes,tolrel,tolabs,mat,fileout);
}
/********************************************************************************/

int  mp_tides_list(MPLinkedFunction fcn,  int *pdd,
				   int nvar, int npar, int nfun, 
				   mpfr_t x[], mpfr_t p[],
				   mpfr_t lt[], int ntes, 	
				   mpfr_t tolrel, mpfr_t tolabs,   
				   mp_data_matrix *mat, FILE* fileout) 
{
	int i, ntot = ntes-1;
	mpfr_t delta[ntot];
	for(i = 0; i < ntot; i++) {
		mpfrts_init(&delta[i]);
		mpfrts_sub (&delta[i], lt[i+1], lt[i]);
	}
	return mp_tides_kernel(fcn,pdd,nvar,npar,nfun,x,p,lt[0],delta[0],delta,ntot,tolrel,tolabs,mat,fileout, NULL);
}

int  mp_tides_delta(MPLinkedFunction fcn,  int *pdd,
					int nvar, int npar, int nfun, 
					mpfr_t x[], mpfr_t p[],
					mpfr_t t0, mpfr_t dt, int ntot,	
					mpfr_t tolrel, mpfr_t tolabs,   
					mp_data_matrix *mat, FILE* fileout)
{
	return mp_tides_kernel(fcn,pdd,nvar,npar,nfun,x,p,t0,dt,NULL,ntot,tolrel,tolabs,mat,fileout, NULL);
}

int  mp_tides_delta_ft(MPLinkedFunction fcn,  int *pdd,
					int nvar, int npar, int nfun, 
					mpfr_t x[], mpfr_t p[],
					mpfr_t t0, mpfr_t tf, mpfr_t dt, 	
					mpfr_t tolrel, mpfr_t tolabs,   
					mp_data_matrix *mat, FILE* fileout) 
{
	int err;
	mpfr_t aux;
	mpfrts_init(&aux);
	mpfrts_set(&aux, tf);
	mpfrts_sub(&aux, aux, t0);
	mpfrts_div(&aux, aux, dt);
	double auxd;
	auxd = mpfrts_get_d(aux);
	int ntot = (int) floor (auxd) ;
	err = mp_tides_delta(fcn,pdd,nvar,npar,nfun,x,p,t0,dt,ntot,tolrel,tolabs,mat,fileout);
	mpfr_clear (aux); 
	mpfr_free_cache ();	
	return err;
}	


/********************************************************************************/
int  mp_tides_kernel(MPLinkedFunction fcn,  int *pdd,
					 int nvar, int npar, int nfun, 
					 mpfr_t x[], mpfr_t p[],
					 mpfr_t t0, mpfr_t dt, mpfr_t *dtl, int ntot,	
					 mpfr_t tolrel, mpfr_t tolabs,   
					 mp_data_matrix *mat, FILE* fileout, mp_data_matrix *der)
{
	int order, i, j,  istep=0, signo, first_time = 1, errorSTEP = 0,ncol, MAX_ORDER;
    int nstep = 0;
	
	mpfr_t tini, tip, tdense, ndt, dtc, sdtc, tstep, ststep, eps, zero, aux, auxb, intT;
	mpfr_t kahan_cg , kahan_cd1 , kahan_cd2 , kahan_ca ; 
	mpfr_t kahan_t, kahan_d;
	mpfrts_init (&tini); 
	mpfrts_init (&tip); 
	mpfrts_init (&tdense); 
	mpfrts_init (&ndt); 
	mpfrts_init (&dtc); 
	mpfrts_init (&sdtc); 
	mpfrts_init (&tstep); 
	mpfrts_init (&ststep); 
	mpfrts_init (&eps); 
	mpfrts_init (&zero); 
	mpfrts_init (&aux); 
	mpfrts_init (&auxb); 
	mpfrts_init (&intT); 
	mpfrts_init (&kahan_cg); 
	mpfrts_init (&kahan_cd1); 
	mpfrts_init (&kahan_cd2); 
	mpfrts_init (&kahan_ca); 
	mpfrts_init (&kahan_t); 
	mpfrts_init (&kahan_d); 

	
	if(_info_steps_taylor_ == 1) mp_set_info_taylor();
	else mp_unset_info_taylor();
	double hant = 0.; 
	
	if (dtl != NULL) {
		mpfrts_set(&ndt, dtl[istep]);
		for(i = 0; i < ntot; i++) mpfrts_add(&intT, intT, dtl[i]);
	} else {
		mpfrts_set (&ndt, dt);
		mpfrts_mul_i(&intT, dt, ntot);
	}
	if (mpfrts_greaterequal(ndt, zero)) signo = 1;
	else signo = -1;
	mpfrts_set (&tini, t0);
	mpfrts_set (&tip, tini);
	mpfrts_set (&dtc, ndt);
	mpfrts_mul_i(&sdtc, dtc, signo); 
	
	order = mp_taylor_order (tolabs); 
	order_series = order;
	MAX_ORDER = order + 1;

	iteration_data itd;
	set_iteration_parameters(&itd, nvar, npar, nfun, order, pdd);

	ncol = (int)fcn (&itd, zero, NULL, NULL, -1, NULL);
	mpfr_t cvfd[ncol*MAX_ORDER];
	mpfr_t y[ncol],yn[ncol];
	for (i=0; i<ncol; i++) {
		mpfrts_init (&y[i]);
		mpfrts_init (&yn[i]);
		for (j=0; j<=order; j++) mpfrts_init (&cvfd[i*MAX_ORDER+j]);			
	}
	for (i=0; i<nvar; i++) mpfrts_set(&y[i], x[i]);
	if(mat != NULL) init_mp_data_matrix(mat, ntot+1, ncol+1);
	if(der != NULL) init_mp_data_matrix(der, ntot+1, nvar);
	
	while (istep <  ntot && errorSTEP == 0 )
	{
		mp_compute_tol (&eps, tolrel, tolabs, nvar, y); 
		fcn (&itd, tip, y, p, order, cvfd);
		if (itd.errorITER)  errorSTEP = 1;
		if (first_time  ) {
			mp_taylor_horner(ncol,order,cvfd, zero, yn, MAX_ORDER);
			if (mat != NULL)  mp_write_taylor_solution(ncol, istep ,tip , yn , &(mat->data), fileout);
			else mp_write_taylor_solution(ncol, istep ,tip , yn , NULL, fileout);
			
			// Derivada en el punto inicial
			if (der != NULL) {
				for(i=0; i<nvar; i++) 
					mpfrts_set (&(der->data)[0][i], cvfd[i*MAX_ORDER+1]);
			}
			first_time =0;
		}
		mp_compute_step(&tstep, &hant, eps, nvar, order, cvfd, MAX_ORDER);
		mpfrts_mul_i (&ststep, tstep, signo);			
		if (defect_error_control) 
			mp_valid_step (fcn, &itd, &ststep, tip, eps, nvar, ncol, order, cvfd, p, MAX_ORDER);
		
		if (test_relminstep) {
			mpfrts_abs(&aux, ststep); 
			mpfrts_mul(&auxb, intT, mp_relminstep); 
			mpfrts_abs(&auxb, auxb); 
		}
		if (test_relminstep && mpfrts_less(aux, auxb)) {
			errorSTEP = 1;
			if(printError) printf("Step reaches a very small value: %.5le.\n Integration ended before the last integration point.\n", ststep);
		} else {
			//Busqueda del siguiente origen de la intgracion
						
			if (kahan_summation) {
				mpfrts_sub (&ststep, ststep, kahan_cg);
				mpfrts_add (&kahan_t, tip, ststep);
				mpfrts_sub (&kahan_cg, kahan_t, tip);
				mpfrts_sub (&kahan_cg, kahan_cg, ststep);
				mpfrts_set (&tip, kahan_t);
			} else mpfrts_add (&tip, tip, ststep);
		
			// Salida densa
			while ( mpfrts_lessequal(sdtc, tstep) && istep < ntot) {
				
				//Nuevo punto de escritura densa
				if (dtl != NULL ) {
					if (kahan_summation) {
						mpfrts_sub (&kahan_d, ndt, kahan_cd2);
						mpfrts_add (&kahan_t, tdense, kahan_d);
						mpfrts_sub (&kahan_cd2, kahan_t, tdense);
						mpfrts_sub (&kahan_cd2, kahan_cd2, kahan_d);
						mpfrts_set (&tdense, kahan_t);
					} else mpfrts_add (&tdense, tdense, ndt);
				} else {
					mpfrts_mul_i(&aux, ndt, istep +1);
					mpfrts_add(&tdense, t0, aux);
				}
							
				mp_taylor_horner(ncol,order,cvfd,dtc,yn, MAX_ORDER);
				if (mat != NULL)  mp_write_taylor_solution(ncol, istep+1 ,tdense ,yn , &(mat->data), fileout);
				else mp_write_taylor_solution(ncol, istep+1 ,tdense ,yn , NULL, fileout);
				// Derivada en el punto final
				if (istep+1 <= ntot && der != NULL) {
					mp_taylor_horner_der(nvar,order,cvfd,dtc,(der->data)[istep+1], MAX_ORDER);
				}
				
				//Incremento para el punto siguiente
				istep++;
				if (dtl != NULL && istep < ntot) mpfrts_set (&ndt, dtl[istep]); 
				if (kahan_summation) {
					mpfrts_sub (&kahan_d, ndt, kahan_cd1);
					mpfrts_add (&kahan_t, dtc, kahan_d);
					mpfrts_sub (&kahan_cd1, kahan_t, dtc);
					mpfrts_sub (&kahan_cd1, kahan_cd1, kahan_d);
					mpfrts_set (&dtc, kahan_t);
				} else mpfrts_add (&dtc, dtc, ndt);
				mpfrts_mul_i(&sdtc, dtc, signo); 

			}
			
			// Continuacion
			
			mp_taylor_horner(ncol, order, cvfd, ststep, y, MAX_ORDER);
			
			//Kahan retroceso
			if (kahan_summation) {
	//				mpfrts_mul_i (&aux, ststep, -1);	
				mpfrts_neg(&aux, ststep);
				mpfrts_sub (&kahan_d, aux, kahan_ca);
				mpfrts_add (&kahan_t, dtc, kahan_d);
				mpfrts_sub (&kahan_ca, kahan_t, dtc);
				mpfrts_sub (&kahan_ca, kahan_ca, kahan_d);
				mpfrts_set (&dtc, kahan_t);
			} else mpfrts_sub (&dtc, dtc, ststep);
			mpfrts_mul_i(&sdtc, dtc, signo); 
		}
        nstep++;
        if (nstep == max_num_steps) {
  			errorSTEP = 1;
			if(printError) printf("Step reaches after: %d steps.\n Integration ended before the last integration point.\n", nstep);          
        }

	}
	
	if(_info_steps_taylor_==1)  mp_str_info_taylor(); 
	delete_iteration_parameters(&itd);
	for (i=0; i<nvar; i++) mpfrts_set (&x[i], yn[i]);
	mpfr_clear (tini); 
	mpfr_clear (tip); 
	mpfr_clear (tdense); 
	mpfr_clear (ndt); 
	mpfr_clear (dtc); 
	mpfr_clear (sdtc); 
	mpfr_clear (tstep); 
	mpfr_clear (ststep); 
	mpfr_clear (eps); 
	mpfr_clear (zero); 
	mpfr_clear (aux); 
	mpfr_clear (auxb); 
	mpfr_clear (intT); 
	mpfr_clear (kahan_cg); 
	mpfr_clear (kahan_cd1); 
	mpfr_clear (kahan_cd2); 
	mpfr_clear (kahan_ca); 
	mpfr_clear (kahan_t); 
	mpfr_clear (kahan_d); 
	for (i=0; i<ncol; i++) {
		mpfr_clear (y[i]);
		mpfr_clear (yn[i]);
		for (j=0; j<=order; j++) mpfr_clear (cvfd[i*MAX_ORDER+j]);			
	}
	mpfr_free_cache ();	
	return errorSTEP;	
}

/*
long mp_number_of_columns(MPLinkedFunction fcn)
{
	int ncol;
	mpfr_t zero;
	mpfrts_init(&zero);
	ncol = (int)fcn (NULL, zero, NULL, NULL, -1, NULL);
	mpfr_clear(zero);
	return (ncol+1);
}
 */

/*****************************************************************/
void init_mp_data_matrix(mp_data_matrix *dm, int r, int c)
{
	dm->rows = r;
	dm->columns = c;
	Array2MP_init(&(dm->data), r, c);
}
void delete_mp_data_matrix(mp_data_matrix *dm)
{
	Array2MP_clear(&(dm->data), dm->rows, dm->columns);
	dm->rows = 0;
	dm->columns = 0;
	dm->data = NULL;
}

