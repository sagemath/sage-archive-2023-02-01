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

#ifndef _commonITER_H
#define _commonITER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>


typedef struct it_dt {
	int		MAX_ORDER;
	long	NDER;
	int		NVARS, NPARS, NFUNS, NLINKS, NPARTIALS; 
	int		*PARTIAL_LIST, *FUNCTION_LIST;
	long	*PREV_ACCUM, *PREV_VI, *PREV_IV, *PREV_COEF; 
	long	*PREVSTAR_ACCUM, *PREVSTAR_VI, *PREVSTAR_IV, *PREVSTAR_COEF;
	long	LINK_ELEMENTS;
	int		clearPartials;
	int		errorITER;
} iteration_data;

void	set_info_error(void);
void	unset_info_error(void);
void	set_maxnumsteps(unsigned long val);


int  same_der(int dim, int *der1, int *der2);
void string_to_der(int dim, char *sder, int *der);
long position_derivative(char* sder, int *pdd);
long position_variable(int v, char* der, int NVARS, int NFUNS, int *pdd);
long position_function(int f, char* der, int NVARS, int NFUNS, int *pdd);

int	 is_variable(iteration_data *itd, int num);

void set_iteration_parameters(iteration_data *itd, int v, int p, int f, int o, int *pdd);
void delete_iteration_parameters(iteration_data *itd);
void set_links(iteration_data *itd, int l, int *flst);
void check_iteration_data_parameters(int t, int pa, int pb);

#endif


