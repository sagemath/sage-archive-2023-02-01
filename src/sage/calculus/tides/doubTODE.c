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

#include "doubTODE.h"


int		_info_steps_taylor_ = 0;
double	etapa_dp_minima;
double	etapa_dp_maxima;
double	etapa_dp_total;
int		num_etapas = 0;
int		order_series = 0;
int		order_estimator = 0;

double 	fac1 = 0.95;
double	fac2 = 10.;
double	fac3 = 0.8;
double	rmaxstep = 1.e2;
double	rminstep = 1.e-2;
double	nitermax = 5;
int		nordinc = 5;
int		minord  = 6;
int		defect_error_control = 0;
int		kahan_summation = 1;


int		test_relminstep = 0;
double  dp_relminstep = 1.e-12;
extern int	printError;
extern unsigned long    max_num_steps;



void	check_relminstep(void)
{
	test_relminstep = 1;
}
void	uncheck_relminstep(void)
{
	test_relminstep = 0;
}

void	dp_set_relminstep(double val)
{
	test_relminstep = 1;
	dp_relminstep = val;
}


void use_default_step_estimator (void) {
	defect_error_control = 1;
}
	
void dp_set_info_taylor(void) {  
	_info_steps_taylor_ = 1;

	etapa_dp_minima = 1.e30;
	etapa_dp_maxima = 0.;
	etapa_dp_total = 0.;
	num_etapas = 0; 
	order_series = 0;
}

void dp_unset_info_taylor(void) {  
	_info_steps_taylor_ = 0;
	etapa_dp_minima = 1.e30;
	etapa_dp_maxima = 0.;
	etapa_dp_total = 0.;
	num_etapas = 0; 
	order_series = 0;
}

void dp_str_info_taylor(void) {
	double tot, num, average;

	tot = etapa_dp_total;
	num = num_etapas;
	average = tot / num;

	if(_info_steps_taylor_ == 1) {
		printf("============================================================\n");
		printf ("Number of integration steps in Taylor method: %d\n", num_etapas);
		double_write ("Minimum   integration step:", etapa_dp_minima);
		double_write ("Maximum   integration step:", etapa_dp_maxima);
		double_write ("Averaged  integration step:",average);
		printf("Order of the Taylor series  : %d\n",  order_series);
		printf("============================================================\n\n");
	}
}


void dp_add_info_step(double tstep) {
	if (tstep < etapa_dp_minima) 
		etapa_dp_minima = tstep;
	if (tstep > etapa_dp_maxima)
		etapa_dp_maxima = tstep;
	num_etapas++;
	etapa_dp_total = etapa_dp_total + tstep;
}


int dp_taylor_order(double tol) {
	double rop;
	int		nrop;
	rop = nordinc - log(tol)/2.0;
	nrop = (int) ceil (rop);
	if(nrop < minord) nrop = minord;
	return nrop;
}

void dp_norm_inf_var(double *rop, int nvar, double *y) 
{
	int i;
	double ninf = 0, ncoef;
	
	for (i=0; i<nvar; i++) {
		ncoef = fabs(y[i]);
		if (ninf < ncoef) ninf = ncoef;
	}
	double_set (rop, ninf);
	
}

void dp_norm_inf(double *rop, int n, int k, double *coef, int MAX_ORDER) {
	int i;
	double ninf = 0, ncoef;
	
	for (i=0; i<n; i++) {
		ncoef = fabs(coef[i*MAX_ORDER+k]);
		if (ninf < ncoef) ninf = ncoef;
	}
	double_set (rop, ninf);
	
}

void dp_compute_step(double *rop, double *hant, double tol, int n, int ord, double *coef, int MAX_ORDER) {
	double tult, tpen, nor, aux;
	int i;
	dp_norm_inf (&nor, n, ord, coef, MAX_ORDER);

	if (nor == 0.) tult = 0.;
	else tult = pow(tol, 1./(ord+1)) * pow(nor, -1./ord);
	
	dp_norm_inf (&nor, n, ord-1, coef, MAX_ORDER);
	if (nor == 0.) tpen = 0.;
	else tpen = pow(tol, 1./ord) * pow(nor, -1./(ord-1));

	if ((tult == 0.) && (tpen == 0.)) {
		i = ord -1;
		*rop = 0.;

		while ((i > 0) && (nor == 0.)) {
			dp_norm_inf (&nor, n, i-1, coef, MAX_ORDER);
			double_set_d (&aux, 1./i);
			double_pow (&aux, tol, aux);
			double_set_d (rop, -1./(i-1));
			double_pow (rop, tol, *rop);
			i--;
		}
	} else if (tult == 0.) *rop = tpen;
	else if (tpen == 0.) *rop = tult;
	else if (tult < tpen) *rop = tult;
	else *rop = tpen;

	*rop = (*rop) * fac1;
	if ((*hant) != 0.) {
		if ((*rop)/(*hant) > rmaxstep) *rop = (*hant) * rmaxstep;
		else if ((*rop)/(*hant) < rminstep) *rop = (*hant) * rminstep;	
	}
	(*hant) = *rop;
	
	if(_info_steps_taylor_ == 1) dp_add_info_step(*rop);
}


void dp_compute_tol (double *tol, double tolrel, double tolabs, int nvar, double *y) 
{
	static double norm_y1 =0.e0;
	double norm_y0, max_norm;

	
	dp_norm_inf_var(&norm_y0, nvar, y);

	if (norm_y0 > norm_y1) max_norm = norm_y0;
	else max_norm = norm_y1;

	*tol = tolrel * max_norm +tolabs;
}


void dp_taylor_horner (int n, int ord, double *coef, double t, double *x, int MAX_ORDER) 
{
	int i;
	for(i = 0; i < n; i++ )
		dp_horner(&coef[i*MAX_ORDER], ord, t, &x[i]);
}
void dp_taylor_horner_der (int n, int ord, double *coef, double t, double *x, int MAX_ORDER) 
{
	int i;
	for(i = 0; i < n; i++ )
		dp_horner_der(&coef[i*MAX_ORDER], ord, t, &x[i]);
}


void dp_write_taylor_solution( int n, int j, double tini, 
	double *x, double*** mat, FILE* fileout) 
{
	int i;
	char formato[] = "%23.15le ";
	if(fileout != NULL) {
		fprintf (fileout, formato, tini);
		for(i = 0; i < n; i++) {
			fprintf (fileout, formato, x[i]);
		}
		fprintf (fileout, "\n");
	}
	if(mat != NULL) {
		(*mat)[j][0] = tini;
		for(i = 0; i < n; i++)
			(*mat)[j][i+1] = x[i];
	}
}




int dp_valid_step (DBLinkedFunction fcn, iteration_data *itd, double *step, double tip, 
		double tol, int nvar, int ncol, int order, double *cvfd, double *p, int MAX_ORDER) {
	int i, j, iter = 1;
	int accepted = 0;
	double b[nvar*order], y[nvar], yp[nvar], t, cn[ncol*MAX_ORDER];
	double dif[nvar], nor, aux, up;

	up  = fac2*tol; 

	while(iter <= nitermax) { 
//	
		for (i=0; i<nvar; i++) 
			for (j=0; j<=order-1; j++)
				b[i*order+j] = cvfd[i*MAX_ORDER+j+1] * (j+1);

		t = tip + *step;
		dp_taylor_horner (nvar, order, cvfd, *step, y, MAX_ORDER);
		dp_taylor_horner (nvar, order-1, b, *step, yp, MAX_ORDER);
		fcn (itd, t, y, p, 1, cn);

		for (i=0; i<nvar; i++) {
			double_sub (&dif[i], yp[i], cn[i*MAX_ORDER+1]);
			aux = fabs(dif[i]);
			if (double_greater (aux, nor)) double_set (&nor, aux);
		}

		if( nor > up) {
			*step = *step * fac3;
		} else break;
//	
		iter++;
	}
	
	if (iter > nitermax) {
		printf("The maximum number of iterations in defect_error_control has been reached.\n");
		printf("The program continues with a stepsize that does not achieve the defect_error_control condition.\n");
	}
	
	return accepted;
}




/********************************************************************************/

int dp_tides_point(DBLinkedFunction fcn, int *pdd,
					int nvar, int npar, int nfun, 
					double *x, double *p,
					double t0, double tf, double dt, 	
					double tolrel, double tolabs,   
					dp_data_matrix *mat, FILE* fileout) 
{
	int ntot = (int) floor ((tf - t0)/dt) ;
	return dp_tides_delta(fcn,pdd,nvar,npar,nfun,x,p,t0,dt,ntot,tolrel,tolabs,mat,fileout);
}

int dp_tides(DBLinkedFunction fcn, int *pdd,
			  int nvar, int npar, int nfun, 
			  double *x, double *p,
			  double *lt, int ntes, 	
			  double tolrel, double tolabs,   
			  dp_data_matrix *mat, FILE* fileout) 
{
	return dp_tides_list(fcn,pdd,nvar,npar,nfun,x,p,lt,ntes,tolrel,tolabs,mat,fileout);
}
	
 
/********************************************************************************/

int dp_tides_list(DBLinkedFunction fcn, int *pdd,
			  int nvar, int npar, int nfun, 
			  double *x, double *p,
			  double *lt, int ntes, 	
			  double tolrel, double tolabs,   
			  dp_data_matrix *mat, FILE* fileout) 
{
	int i, ntot = ntes-1;
	double delta[ntot];
	for(i = 0; i < ntot; i++) delta[i] = lt[i+1]-lt[i];
	return dp_tides_kernel(fcn,pdd,nvar,npar,nfun,x,p,lt[0],delta[0],delta,ntot,tolrel,tolabs,mat,fileout, NULL);
}

int  dp_tides_delta(DBLinkedFunction fcn, int *pdd,
					 int nvar, int npar, int nfun, 
					 double *x, double *p,
					 double t0, double dt, int ntot,	
					 double tolrel, double tolabs,   
					 dp_data_matrix *mat, FILE* fileout)
{
	return dp_tides_kernel(fcn,pdd,nvar,npar,nfun,x,p,t0,dt,NULL,ntot,tolrel,tolabs,mat,fileout, NULL);
}

int  dp_tides_delta_ft(DBLinkedFunction fcn, int *pdd,
					int nvar, int npar, int nfun, 
					double *x, double *p,
					double t0, double tf, double dt, 	
					double tolrel, double tolabs,   
					dp_data_matrix *mat, FILE* fileout) 
{
	int ntot = (int) floor ((tf - t0)/dt) ;
	return dp_tides_delta(fcn,pdd,nvar,npar,nfun,x,p,t0,dt,ntot,tolrel,tolabs,mat,fileout);
}


/********************************************************************************/


int dp_tides_kernel(DBLinkedFunction fcn, int *pdd,
					 int nvar, int npar, int nfun, 
					 double x[], double p[],
					 double t0, double dt, double *dtl, int ntot,	
					 double tolrel, double tolabs,   
					 dp_data_matrix *mat, FILE* fileout, dp_data_matrix *der)
{
	int order, i, istep=0, signo, first_time = 1, errorSTEP = 0,ncol, MAX_ORDER;
    int nstep = 0;
	
	double tini, tip, tipo, tdense, ndt, dtc, tstep, ststep, eps, hant, intT;
	volatile double kahan_cg = 0., kahan_cd1 = 0., kahan_cd2 = 0., kahan_ca = 0.0; 
	volatile double kahan_t,  kahan_d; 
	
	if(_info_steps_taylor_ == 1) dp_set_info_taylor();
	else dp_unset_info_taylor();
	hant = 0.; 
	
	if (dtl != NULL) {
		ndt = dtl[istep];
		intT = 0;
		for(i = 0; i < ntot; i++) intT +=dtl[i];
	} else {
		ndt = dt;
		intT = ntot * dt;
	}
	if (ndt >= 0) signo = 1;
	else signo = -1;
	tini = t0; 
	tip = tini;
	dtc = ndt;
	
	order = dp_taylor_order (tolabs); 
	order_series = order;
	MAX_ORDER = order + 1;
	
	iteration_data itd;
	set_iteration_parameters(&itd, nvar, npar, nfun, order, pdd);
	double ***dmat;
	if(mat == NULL) dmat = NULL;
	else dmat = &(mat->data);
	
	ncol = (int)fcn (&itd, 0., NULL, NULL, -1, NULL);
	double cvfd[ncol*MAX_ORDER];
	double y[ncol],yn[ncol];
	for(i=0; i<nvar; i++) y[i] = x[i];	
	if(mat != NULL) init_dp_data_matrix(mat, ntot+1, ncol+1);
	if(der != NULL) init_dp_data_matrix(der, ntot+1, nvar);
		
	while (istep <  ntot && errorSTEP == 0 )
	{
		dp_compute_tol (&eps, tolrel, tolabs, nvar, y);
		fcn (&itd, tip, y, p, order, cvfd);
		if (itd.errorITER)  errorSTEP = 1;
		if (first_time) {
			dp_taylor_horner(ncol,order,cvfd, 0., yn, MAX_ORDER);
 			dp_write_taylor_solution(ncol, istep ,tip , yn , dmat, fileout);
			// Derivada en el punto inicial
			if (der != NULL) {
				for(i=0; i<nvar; i++) der->data[0][i] = cvfd[i*MAX_ORDER+1];
			}
			first_time =0;
		}
		dp_compute_step(&tstep, &hant, eps, nvar, order, cvfd, MAX_ORDER);
		ststep = tstep * signo;
		if (defect_error_control)
			dp_valid_step (fcn, &itd, &ststep, tip, eps, nvar, ncol, order, cvfd, p, MAX_ORDER);
		
		if (test_relminstep && fabs(ststep) < fabs(intT*dp_relminstep)) {
			errorSTEP = 1;
			if(printError) printf("Step reaches a very small value: %.5le.\n Integration ended before the last integration point.\n", ststep);
		} else {
			//Busqueda del siguiente origen de la integracion
			
			tipo = tip;
			if(kahan_summation) {
				ststep = ststep - kahan_cg;
				kahan_t = tip + ststep;
				kahan_cg = (kahan_t - tip) - ststep;
				tip = kahan_t;
			} else tip += ststep;
			
			
			// Salida densa
			while ((signo*dtc) <= (signo*ststep) && istep < ntot) {

				//Nuevo punto de escritura densa
				if (dtl != NULL ) {
					if(kahan_summation) {
						kahan_d = ndt - kahan_cd2;
						kahan_t = tdense +  kahan_d;
						kahan_cd2 = (kahan_t - tdense) - kahan_d;
						tdense = kahan_t;
					} else tdense += ndt;
				} else tdense = t0 + (istep +1) * ndt;
							
				dp_taylor_horner(ncol,order,cvfd,dtc,yn, MAX_ORDER);
				dp_write_taylor_solution(ncol, istep+1 ,tdense ,yn , dmat, fileout);
				// Derivada en el punto final
				if (istep+1 <= ntot && der != NULL) {
					dp_taylor_horner_der(nvar,order,cvfd,dtc,der->data[istep+1], MAX_ORDER);
				}
				
				istep++;
				if (dtl != NULL && istep < ntot) ndt = dtl[istep];
				//Incremento para el punto siguiente
				if(kahan_summation) {
					kahan_d = ndt - kahan_cd1;
					kahan_t = dtc +  kahan_d;
					kahan_cd1 = (kahan_t - dtc) - kahan_d;
					dtc = kahan_t;
				} else dtc += ndt;
			}
			
			// Continuacion
			
			dp_taylor_horner(ncol, order, cvfd, ststep, y, MAX_ORDER);
			
			//Kahan retroceso
			
			if(kahan_summation) {
				kahan_d = -ststep - kahan_ca;
				kahan_t = dtc +  kahan_d;
				kahan_ca = (kahan_t - dtc) - kahan_d;
				dtc = kahan_t;
			} else dtc -= ststep;
		}
        
        nstep++;
        if (nstep == max_num_steps) {
  			errorSTEP = 1;
			if(printError) printf("Step reaches after: %d steps.\n Integration ended before the last integration point.\n", nstep);          
        }

	}
	if(_info_steps_taylor_==1)  dp_str_info_taylor();
	delete_iteration_parameters(&itd);
	for (i=0; i<nvar; i++) x[i] = yn[i];
	return errorSTEP;
}
/*
long dp_number_of_columns(DBLinkedFunction fcn)
{
	int ncol;
	ncol = (int)fcn (NULL, 0., NULL, NULL, -1, NULL);
	return (ncol+1);
}
*/
/*****************************************************************/
void init_dp_data_matrix(dp_data_matrix *dm, int r, int c)
{
	dm->rows = r;
	dm->columns = c;
	Array2DB_init(&(dm->data), r, c);
}
void delete_dp_data_matrix(dp_data_matrix *dm)
{
	Array2DB_clear(&(dm->data), dm->rows, dm->columns);
	dm->rows = 0;
	dm->columns = 0;
	dm->data = NULL;
}

