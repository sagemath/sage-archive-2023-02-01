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

#include "mpfrEVENTS.h"


extern int	_info_steps_taylor_;
extern int	order_series;
extern int	defect_error_control;
extern int	kahan_summation;
extern int	DECIMAL_PRECISION;
extern int  noPartialDerivative[20];
extern double	hant;

extern int	printError;
extern int	test_relminstep;
extern mpfr_t	mp_relminstep;
extern unsigned long    max_num_steps;


void mp_init_event_list(mp_event_list *lista, int ncol)
{
	lista->total = 0;
	lista->dim = ncol;
	lista->first = NULL;
}

void mp_add_to_event_list(mp_event_list *lista, mpfr_t *val)
{
	int i;
	mp_event_node *nodo, *lnodo;
	nodo = (mp_event_node *) malloc(sizeof(mp_event_node));

	nodo->data = (mpfr_t *) calloc(lista->dim,sizeof(mpfr_t));
	for(i = 0; i < lista->dim; i++) {
		mpfrts_init(&(nodo->data[i]));
		mpfrts_set(&(nodo->data[i]), val[i]);
	}
	nodo->next = NULL;
	if(lista->first == NULL) lista->first = nodo;
	else {
		lnodo = lista->first;
		while (lnodo->next != NULL) lnodo = lnodo->next;
		lnodo->next = nodo;
	}
	lista->total++;
}

void mp_event_list_to_array(mp_event_list lista, mp_data_matrix *array)
{
	int i,j;
	mp_event_node *lnodo;
	lnodo = lista.first;
	for (i=0; i< lista.total; i++) {
		for (j=0; j< lista.dim; j++)  
			mpfrts_set(&(array->data[i][j]), lnodo->data[j]);
		lnodo = lnodo->next;
	}
}


void mp_delete_event_list(mp_event_list *lis)
{
	mp_event_node *noda, *nodp;
	noda = lis->first ;
	while (noda != NULL) {
		nodp = noda->next;
		free(noda->data);
		free(noda);
		noda = nodp;
	}
	lis->first = NULL;
	lis->dim = 0;
	lis->total = 0;
}

/*****************************************************************************/

int  mp_tides_find_zeros(MPLinkedFunction fcn, 
						 int nvar, int npar,mpfr_t *x, mpfr_t *p,
						 mpfr_t tini, mpfr_t tend, mpfr_t tol, 
						 int *numevents, mp_data_matrix *events, FILE* fileout) 
{
	return mp_tides_events(fcn,nvar,npar,x,p,tini,tend,tol,numevents,events,fileout,1);
}

int  mp_tides_find_extrema(MPLinkedFunction fcn, 
						   int nvar, int npar,mpfr_t *x, mpfr_t *p,
						   mpfr_t tini, mpfr_t tend, mpfr_t tol, 
						   int *numevents, mp_data_matrix *events, FILE* fileout) 
{
	return mp_tides_events(fcn,nvar,npar,x,p,tini,tend,tol,numevents,events,fileout,2);
}


int  mp_tides_find_minimum(MPLinkedFunction fcn, 
						   int nvar, int npar,mpfr_t *x, mpfr_t *p,
						   mpfr_t tini, mpfr_t tend, mpfr_t tol, 
						   int *numevents, mp_data_matrix *events, FILE* fileout) 
{
	return mp_tides_events(fcn,nvar,npar,x,p,tini,tend,tol,numevents,events,fileout,3);
}

int  mp_tides_find_maximum(MPLinkedFunction fcn, 
						   int nvar, int npar,mpfr_t *x, mpfr_t *p,
						   mpfr_t tini, mpfr_t tend, mpfr_t tol, 
						   int *numevents, mp_data_matrix *events, FILE* fileout)
{
	return mp_tides_events(fcn,nvar,npar,x,p,tini,tend,tol,numevents,events,fileout,4);
}

int  mp_tides_events(MPLinkedFunction fcn, 
					 int nvar, int npar,mpfr_t *x, mpfr_t *p,
					 mpfr_t tini, mpfr_t tend, mpfr_t tol, 
					 int *numevents, mp_data_matrix *events, 
					 FILE* fileout, int evcase) 
{
	int order, i, j, num, nz = 0, end = 0;
	int cond1, cond2, signo, MAX_ORDER, errorSTEP=0;
	int ncol;
    int nstep = 0;
	double hant = 0.;
	
	mpfr_t tip, tstep, ststep, ta,tb, zero, aux;
	mpfr_t dlt, eps, traiz, tolrel, tolabs;	
	mpfr_t temp, kahan_c, kahan_t;
	mpfrts_init (&tip); 
	mpfrts_init (&tstep); 
	mpfrts_init (&ststep);
	mpfrts_init (&ta); 
	mpfrts_init (&tb); 
	mpfrts_init (&zero); 
	mpfrts_init (&aux); 
	mpfrts_init (&dlt); 
	mpfrts_init (&eps); 
	mpfrts_init (&traiz);
	mpfrts_init (&tolrel); 
	mpfrts_init (&tolabs); 
	mpfrts_init (&temp); 
	mpfrts_init (&kahan_c); 
	mpfrts_init (&kahan_t);
	
	mp_event_list tlist;	
	
	mpfrts_set_i(&temp, 10);
	mpfrts_pow_i(&tolrel, temp, - DECIMAL_PRECISION);
	mpfrts_set(&tolabs, tolrel);
	
	if(_info_steps_taylor_ == 1) mp_set_info_taylor();
	else mp_unset_info_taylor();
	
	order = mp_taylor_order (tolabs);
	order_series = order;
	MAX_ORDER = order + 1;

	iteration_data itd;
	set_iteration_parameters(&itd, nvar, npar, 1, order, noPartialDerivative);

	if (mpfrts_greater(tend, tini)) signo = 1;
	else signo = -1;
	ncol = (int)fcn (&itd, zero, NULL, NULL, -1, NULL);
	mpfr_t cvfd[ncol*MAX_ORDER];
	mpfr_t y[ncol], yn[ncol], ynt[ncol+1];;
	for (i=0; i<ncol; i++) {
		mpfrts_init (&y[i]);
		mpfrts_init (&yn[i]);
		mpfrts_init (&ynt[i]);
		for (j=0; j<=order; j++) mpfrts_init (&cvfd[i*MAX_ORDER+j]);			
	}
	mpfrts_init (&ynt[ncol]);
	for (i=0; i<nvar; i++) mpfrts_set(&y[i], x[i]);
	mp_init_event_list(&tlist,ncol+1);
	
	mpfrts_set(&tip, tini);
	
	cond1 = ((signo ==  1 && mpfrts_lessequal(tip, tend)) || 
			 (signo == -1 && mpfrts_greaterequal(tip, tend)));
	if(*numevents) cond2 = (nz < *numevents);
	else cond2 = 1;
	
	while (cond1 && cond2 &&  errorSTEP == 0) {
		
		mp_compute_tol (&eps, tolrel, tolabs, nvar, y);
		fcn (&itd, tip, y, p, order, cvfd);
		if(itd.errorITER) errorSTEP = 1;
		mp_compute_step(&tstep, &hant, eps, nvar, order, cvfd, MAX_ORDER);
		mpfrts_mul_i (&ststep, tstep, signo);			
		if (defect_error_control)
			mp_valid_step (fcn, NULL, &ststep, tip, eps, nvar, ncol, order, cvfd, p, MAX_ORDER);
		
		mpfrts_sub(&dlt, tend, tip);
		
		mpfrts_abs(&temp, dlt);
		if(mpfrts_less(temp, tstep)) end = 1;
		else mpfrts_set(&dlt, ststep);
		
		if(signo == 1) {
			mpfrts_set(&ta, zero);
			mpfrts_set(&tb, dlt);
		} else {
			mpfrts_set(&tb, zero);
			mpfrts_set(&ta, dlt);
		}
		

		switch (evcase) {
			case 1:
				num = mp_one_zero(&cvfd[(ncol-1)*MAX_ORDER], order, ta, tb, tol, &traiz);
				break;
			case 2:
				num = mp_one_extremum(&cvfd[(ncol-1)*MAX_ORDER], order, ta, tb, tol, &traiz);
				break;
			case 3:
				num = mp_one_minimum(&cvfd[(ncol-1)*MAX_ORDER], order, ta, tb, tol, &traiz);
				break;
			case 4:
				num = mp_one_maximum(&cvfd[(ncol-1)*MAX_ORDER], order, ta, tb, tol, &traiz);
				break;
			default:
				break;
		}
		
		if(num) {
			mp_taylor_horner(ncol,order,cvfd,traiz,yn, MAX_ORDER);
			for(i=0; i < ncol; i++) mpfrts_set(&ynt[i+1], yn[i]);
			mpfrts_add(&ynt[0], tip, traiz);
			mp_add_to_event_list(&tlist, ynt);
			if(fileout !=NULL) mp_write_taylor_solution(ncol,nz,ynt[0],yn,NULL,fileout);
			nz++;
		}
		
		if (end) {
			mp_taylor_horner(ncol,order,cvfd,dlt,y, MAX_ORDER);
			if(_info_steps_taylor_==1)  mp_str_info_taylor(); 
			for(i=0; i<nvar; i++) mpfrts_set(&x[i], y[i]);
			*numevents = tlist.total;
			if (events != NULL) {
				init_mp_data_matrix(events, *numevents, tlist.dim);
				mp_event_list_to_array(tlist, events);
			}
			mp_delete_event_list(&tlist);
			delete_iteration_parameters(&itd);
			mpfr_clear (tip); 
			mpfr_clear (tstep); 
			mpfr_clear (ststep);
			mpfr_clear (ta); 
			mpfr_clear (tb); 
			mpfr_clear (zero); 
			mpfr_clear (dlt); 
			mpfr_clear (eps); 
			mpfr_clear (traiz);
			mpfr_clear (tolrel); 
			mpfr_clear (tolabs); 
			mpfr_clear (temp); 
			mpfr_clear (kahan_c); 
			mpfr_clear (kahan_t);
			for (i=0; i<ncol; i++) {
				mpfr_clear (y[i]);
				mpfr_clear (yn[i]);
				for (j=0; j<=order; j++) mpfr_clear (cvfd[i*MAX_ORDER+j]);			
			}
			mpfr_free_cache ();	
			return errorSTEP; 
		} 
		
		if (kahan_summation) {
			mpfrts_sub (&ststep, ststep, kahan_c);
			mpfrts_add (&kahan_t, tip, ststep);
			mpfrts_sub (&kahan_c, kahan_t, tip);
			mpfrts_sub (&kahan_c, kahan_c, ststep);
			mpfrts_set (&tip, kahan_t);
		} else mpfrts_add (&tip, tip, ststep);
		mp_taylor_horner(ncol, order, cvfd, ststep, y, MAX_ORDER);
		
		cond1 = ((signo ==  1 && mpfrts_lessequal(tip, tend)) || 
				 (signo == -1 && mpfrts_greaterequal(tip, tend)));
		if(*numevents) cond2 = (nz < *numevents);
		else cond2 = 1;
		
        nstep++;
        if (nstep == max_num_steps) {
  			errorSTEP = 1;
            end = 1;
			if(printError) printf("Step reaches after: %d steps.\n Integration ended before the last integration point.\n", nstep);          
        }
	}
	
	mpfr_clear (tip); 
	mpfr_clear (tstep); 
	mpfr_clear (ststep);
	mpfr_clear (ta); 
	mpfr_clear (tb); 
	mpfr_clear (zero); 
	mpfr_clear (aux); 
	mpfr_clear (dlt); 
	mpfr_clear (eps); 
	mpfr_clear (traiz);
	mpfr_clear (tolrel); 
	mpfr_clear (tolabs); 
	mpfr_clear (temp); 
	mpfr_clear (kahan_c); 
	mpfr_clear (kahan_t);
	for (i=0; i<ncol; i++) {
		mpfr_clear (y[i]);
		mpfr_clear (yn[i]);
		for (j=0; j<=order; j++) mpfr_clear (cvfd[i*MAX_ORDER+j]);			
	}
	mpfr_free_cache ();	
	
	
	if(_info_steps_taylor_==1)  mp_str_info_taylor(); 
	for(i=0; i<nvar; i++) mpfrts_set(&x[i], y[i]);
    *numevents = tlist.total;
	if (events != NULL) {
		init_mp_data_matrix(events, *numevents, tlist.dim);
		mp_event_list_to_array(tlist, events);
	}
	mp_delete_event_list(&tlist);
	delete_iteration_parameters(&itd);
	return errorSTEP; 
	
}
