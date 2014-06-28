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

#include "doubEVENTS.h"

extern int	_info_steps_taylor_ ;
extern int 	order_series;
extern int	defect_error_control;
extern int	kahan_summation;
extern int  noPartialDerivative[20];

extern int	printError;
extern int	test_relminstep;
extern double	dp_relminstep;
extern unsigned long    max_num_steps;


void dp_init_event_list(dp_event_list *lista, int ncol)
{
	lista->total = 0;
	lista->dim = ncol;
	lista->first = NULL;
}

void dp_add_to_event_list(dp_event_list *lista, double *val)
{
	int i;
	dp_event_node *nodo, *lnodo;
	nodo = (dp_event_node *) malloc(sizeof(dp_event_node));
	nodo->data = (double *) calloc(lista->dim,sizeof(double));
	for(i = 0; i < lista->dim; i++) nodo->data[i] = val[i];
	nodo->next = NULL;
	if(lista->first == NULL) lista->first = nodo;
	else {
		lnodo = lista->first;
		while (lnodo->next != NULL) lnodo = lnodo->next;
		lnodo->next = nodo;
	}
	lista->total++; 
}


void dp_event_list_to_array(dp_event_list lista, dp_data_matrix *array)
{
	int i,j;
	dp_event_node *lnodo;
	lnodo = lista.first;
	for (i=0; i< lista.total; i++) {
		for (j=0; j< lista.dim; j++)  (array->data)[i][j]= lnodo->data[j];
		lnodo = lnodo->next;
	}
}

void dp_delete_event_list(dp_event_list *lis)
{
	dp_event_node *noda, *nodp;
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

int  dp_tides_find_zeros(DBLinkedFunction fcn, 
						 int nvar, int npar,double *x, double *p,
						 double tini, double tend, double tol, 
						 int *numevents, dp_data_matrix *events, FILE* fileout) 
{
	return dp_tides_events(fcn,nvar,npar,x,p,tini,tend,tol,numevents,events,fileout,1);
}

int  dp_tides_find_extrema(DBLinkedFunction fcn, 
						   int nvar, int npar,double *x, double *p,
						   double tini, double tend, double tol, 
						   int *numevents, dp_data_matrix *events, FILE* fileout) 
{
	return dp_tides_events(fcn,nvar,npar,x,p,tini,tend,tol,numevents,events,fileout,2);
}


int  dp_tides_find_minimum(DBLinkedFunction fcn, 
						   int nvar, int npar,double *x, double *p,
						   double tini, double tend, double tol, 
						   int *numevents, dp_data_matrix *events, FILE* fileout) 
{
	return dp_tides_events(fcn,nvar,npar,x,p,tini,tend,tol,numevents,events,fileout,3);
}

int  dp_tides_find_maximum(DBLinkedFunction fcn, 
						   int nvar, int npar,double *x, double *p,
						   double tini, double tend, double tol, 
						   int *numevents, dp_data_matrix *events, FILE* fileout)
{
	return dp_tides_events(fcn,nvar,npar,x,p,tini,tend,tol,numevents,events,fileout,4);
}

int  dp_tides_events(DBLinkedFunction fcn, 
					 int nvar, int npar,double *x, double *p,
					 double tini, double tend, double tol, 
					 int *numevents, dp_data_matrix *events,  
					 FILE* fileout, int evcase) 
{
	int order, i, num, nz = 0, end = 0;
	int cond1, cond2, signo, MAX_ORDER, errorSTEP=0;
	int ncol;
    int nstep = 0;
	
	double tip, tstep, ststep, ta,tb;
	double dlt, eps, traiz, tolrel, tolabs;	
	double kahan_c = 0., kahan_t, hant;
	dp_event_list tlist;
	hant = 0.;
	tolrel = tolabs = 1.e-16;
	
	if(_info_steps_taylor_ == 1) dp_set_info_taylor();
	else dp_unset_info_taylor();
	
	order = dp_taylor_order (tolabs); 
	order_series = order;
	MAX_ORDER = order + 1;

	iteration_data itd;
	set_iteration_parameters(&itd, nvar, npar, 1, order, noPartialDerivative);

	if (tend > tini) signo = 1;
	else signo = -1;
	ncol = (int)fcn (&itd, 0., NULL, NULL, -1, NULL);
	double cvfd[ncol*MAX_ORDER];
	double y[ncol], yn[ncol],ynt[ncol+1];
	dp_init_event_list(&tlist,ncol+1);
	tip = tini;
	for(i=0; i<nvar; i++) y[i] = x[i];
	
	cond1 = ((signo == 1 && tip <= tend) || (signo == -1 && tip >= tend));
	if(*numevents) cond2 = (nz < *numevents);
	else cond2 = 1;
	
	while (cond1 && cond2 &&  errorSTEP == 0) {
		
		dp_compute_tol (&eps, tolrel, tolabs, nvar, y);
		fcn (&itd, tip, y, p, order, cvfd);
		if(itd.errorITER) errorSTEP = 1;
		dp_compute_step(&tstep, &hant, eps, nvar, order, cvfd, MAX_ORDER);
		ststep = tstep * signo;
		if (defect_error_control)
			dp_valid_step (fcn, &itd, &ststep, tip, eps, nvar, ncol, order, cvfd, p, MAX_ORDER);
		
		dlt = tend - tip;
		
		if(fabs(dlt) < tstep) end = 1;
		else dlt = ststep;
		
		if(signo == 1) {
			ta = 0;
			tb = dlt;
		} else {
			tb = 0;
			ta = dlt;
		}
		if (test_relminstep && fabs(dlt) < tip*dp_relminstep) {
			end = 1;
			num = 0;
			errorSTEP = 1;
			if(printError) printf("Step reaches a very small value: %.5le.\n Integration ended before the last integration point.\n", dlt);
		} else {
			switch (evcase) {
				case 1:
					num = dp_one_zero(&cvfd[(ncol-1)*MAX_ORDER], order, ta, tb, tol, &traiz);
					break;
				case 2:
					num = dp_one_extremum(&cvfd[(ncol-1)*MAX_ORDER], order, ta, tb, tol, &traiz);
					break;
				case 3:
					num = dp_one_minimum(&cvfd[(ncol-1)*MAX_ORDER], order, ta, tb, tol, &traiz);
					break;
				case 4:
					num = dp_one_maximum(&cvfd[(ncol-1)*MAX_ORDER], order, ta, tb, tol, &traiz);
					break;
				default:
					break;
			}
		}
		
		if(num) {
			dp_taylor_horner(ncol,order,cvfd,traiz,yn, MAX_ORDER);
			for(i=0; i < ncol; i++) ynt[i+1] = yn[i];
			ynt[0] = tip+traiz;
			dp_add_to_event_list(&tlist, ynt);
			if(fileout !=NULL) dp_write_taylor_solution(ncol,nz,ynt[0],yn,NULL,fileout);
			nz++;
		}
		
		if (end) {
			dp_taylor_horner(ncol,order,cvfd,dlt,y, MAX_ORDER);
			if(_info_steps_taylor_==1)  dp_str_info_taylor(); 
			for(i=0; i<nvar; i++) x[i] = y[i];
			*numevents = tlist.total;
			if (events != NULL) {
				init_dp_data_matrix(events, *numevents, tlist.dim);
				dp_event_list_to_array(tlist, events);
			}
			dp_delete_event_list(&tlist);
			delete_iteration_parameters(&itd);
			return errorSTEP; 
		} 
		
		if(kahan_summation) {
			ststep = ststep - kahan_c;
			kahan_t = tip + ststep;
			kahan_c = (kahan_t - tip) - ststep;
			tip = kahan_t;
		} else tip +=ststep;
		dp_taylor_horner(ncol, order, cvfd, ststep, y, MAX_ORDER);
		
		cond1 = ((signo == 1 && tip <= tend) || (signo == -1 && tip >= tend));
		if(*numevents) cond2 = (nz < *numevents);
		else cond2 = 1;
        
        nstep++;
        if (nstep == max_num_steps) {
  			errorSTEP = 1;
            end = 1;
			if(printError) printf("Step reaches after: %d steps.\n Integration ended before the last integration point.\n", nstep);          
        }

		
	}
	if(_info_steps_taylor_==1)  dp_str_info_taylor(); 
	for (i=0; i<nvar; i++) x[i] = y[i];
    *numevents = tlist.total;
	if (events != NULL) {
		init_dp_data_matrix(events, *numevents, tlist.dim);
		dp_event_list_to_array(tlist, events);
	}
	dp_delete_event_list(&tlist);
	delete_iteration_parameters(&itd);
	return errorSTEP; 
	
}


