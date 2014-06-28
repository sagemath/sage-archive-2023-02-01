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

#ifndef _doubleEVENTS_H
#define _doubleEVENTS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "doubTODE.h"
#include "doubPOL.h"
#include "commonITER.h"


typedef struct dp_node_event {
	double *data;
	struct dp_node_event *next;
} dp_event_node;

typedef struct dp_linked_events_nodes {
	int total;
	int dim;
	dp_event_node *first;
} dp_event_list;

void dp_init_event_list(dp_event_list *lista, int ncol);
void dp_add_to_event_list(dp_event_list *lista, double *val);
void dp_event_list_to_array(dp_event_list lista, dp_data_matrix *array);

void dp_delete_event_list(dp_event_list *lis);

/*****************************************************************************/


int  dp_tides_find_zeros(DBLinkedFunction fcn, 
						 int nvar, int npar,double *x, double *p,
						 double tini, double tend, double tol, 
						 int *numevents, dp_data_matrix *events, FILE* fileout) ;

int  dp_tides_find_extrema(DBLinkedFunction fcn, 
						   int nvar, int npar,double *x, double *p,
						   double tini, double tend, double tol,
						   int *numevents, dp_data_matrix *events, FILE* fileout) ;

int  dp_tides_find_minimum(DBLinkedFunction fcn, 
						   int nvar, int npar,double *x, double *p,
						   double tini, double tend, double tol, 
						   int *numevents, dp_data_matrix *events, FILE* fileout) ;

int  dp_tides_find_maximum(DBLinkedFunction fcn, 
						   int nvar, int npar,double *x, double *p,
						   double tini, double tend, double tol, 
						   int *numevents, dp_data_matrix *events, FILE* fileout);

int  dp_tides_events(DBLinkedFunction fcn, 
					 int nvar, int npar,double *x, double *p,
					 double tini, double tend, double tol, 
					 int *numevents, dp_data_matrix *events,  
					 FILE* fileout, int evcase) ;
#endif

