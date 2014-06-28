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

/*  V-20 */

#ifndef _mpfrEVENTS_H
#define _mpfrEVENTS_H

#include <stdio.h>
#include <stdlib.h>
#include "mpfr.h"
#include <math.h>
#include "mpfrTODE.h"
#include "mpfrPOL.h"
#include "commonITER.h"


typedef struct mp_node_event {
	mpfr_t *data;
	struct mp_node_event *next;
} mp_event_node;

typedef struct mp_linked_events_nodes {
	int total;
	int dim;
	mp_event_node *first;
} mp_event_list;

void mp_init_event_list(mp_event_list *lista, int ncol);
void mp_add_to_event_list(mp_event_list *lista, mpfr_t *val);
void mp_event_list_to_array(mp_event_list lista, mp_data_matrix *array);
void mp_delete_event_list(mp_event_list *lis);


int  mp_tides_find_zeros(MPLinkedFunction fcn, 
						 int nvar, int npar,mpfr_t *x, mpfr_t *p,
						 mpfr_t tini, mpfr_t tend, mpfr_t tol, 
						 int *numevents, mp_data_matrix *events, FILE* fileout) ;

int  mp_tides_find_extrema(MPLinkedFunction fcn, 
						   int nvar, int npar,mpfr_t *x, mpfr_t *p,
						   mpfr_t tini, mpfr_t tend, mpfr_t tol, 
						   int *numevents, mp_data_matrix *events, FILE* fileout);

int  mp_tides_find_minimum(MPLinkedFunction fcn, 
						   int nvar, int npar,mpfr_t *x, mpfr_t *p,
						   mpfr_t tini, mpfr_t tend, mpfr_t tol, 
						   int *numevents, mp_data_matrix *events, FILE* fileout) ;

int  mp_tides_find_maximum(MPLinkedFunction fcn, 
						   int nvar, int npar,mpfr_t *x, mpfr_t *p,
						   mpfr_t tini, mpfr_t tend, mpfr_t tol, 
						   int *numevents,mp_data_matrix *events, FILE* fileout);

int  mp_tides_events(MPLinkedFunction fcn, 
					 int nvar, int npar,mpfr_t *x, mpfr_t *p,
					 mpfr_t tini, mpfr_t tend, mpfr_t tol, 
					 int *numevents, mp_data_matrix *events, 
					 FILE* fileout, int evcase) ;


#endif

