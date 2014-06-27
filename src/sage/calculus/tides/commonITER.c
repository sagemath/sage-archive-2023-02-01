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

int   noPartialDerivative[20]  = {0, 1, 2, 0, 1, 1, 1, 1, 0, 1, 0, 2, 0, 1, 1, 1, 1, 0, 1, 0};

int	  printError = 1;

unsigned long    max_num_steps = ULONG_MAX;


void	set_maxnumsteps(unsigned long val)
{
	max_num_steps = val;
}

void set_info_error(void)
{ 
	printError = 1;
}

void unset_info_error(void)
{ 
	printError = 0;
}


int is_variable(iteration_data *itd, int num)
{
	return (num <= itd->NVARS);
}

int same_der(int dim, int *der1, int *der2)
{
	int i;
	for (i = 0; i < dim; i++) if(der1[i] != der2[i]) return  0;
	return 1; 
}

void string_to_der(int dim, char *sder, int *der)
{
	int i=0, k=0, j; 
	char *pder, nder[dim];
	for(j = 0; j < dim; j++) der[j] = 0;
	strcpy(nder, sder); 
	j=0;
	while (sder[j] != '\0') if(sder[j++] =='/') k++;
	if(k != dim -1) {
		der[0] = -1;
		return;
	}

	pder = strtok (nder,"/");
	while (pder != NULL && i < dim)
	{
		der[i++] = atoi(pder);
		pder= strtok (NULL, "/");
	}
}


long position_derivative(char* sder, int *pdd)
{
	long i=0, j, k, num, sal = -1;
	int NPARTIALS, NDER;
	if (pdd != NULL)  {
		NPARTIALS = pdd[i++];
		for(j=0; j<NPARTIALS; j++) i++;
		NDER = pdd[i++]; 
		for(k = 0; k <8; k++) {
			num = pdd[i++];
			for(j=0; j<num; j++) i++;
		}
		int PD_NAMES[NDER][NPARTIALS];
		for(j = 0; j <NDER; j++) {
			for(k = 0; k < NPARTIALS; k++) PD_NAMES[j][k] = pdd[i++];
		}
		int der[NPARTIALS];
		string_to_der(NPARTIALS, sder, der);
		for(i=0; i < NDER; i++)
			if(same_der(NPARTIALS, der, PD_NAMES[i])) sal =  i;
	}
	return sal;
}

long position_variable(int v, char* der, int NVARS, int NFUNS, int *pdd)
{
	long nvf, pd;
	nvf = NVARS+NFUNS;
	pd  = position_derivative(der, pdd);
	if(pd <= -1 || v >= NVARS || v < 0 ) return -1;
	return (v+(pd*nvf)+1);
}
long position_function(int f, char* der, int NVARS, int NFUNS, int *pdd)
{
	long nvf, pd;
	nvf = NVARS+NFUNS;
	pd  = position_derivative(der,pdd);
	if(pd <= -1 || f >= NFUNS || f < 0 ) return -1;
	return (NVARS+f+(pd*nvf)+1);
}

 
/*************************************************************/
void set_iteration_parameters(iteration_data *itd, int v, int p, int f, int o, int *pdd)
{
	itd->NVARS			= v;
	itd->NPARS			= p;
	itd->NFUNS			= f;
	itd->MAX_ORDER		= o + 1;
	itd->clearPartials  =  1;
	itd->errorITER      =  0;

	
	int i=0, j, num, *pddata;
	
	if (pdd == NULL) pddata	= noPartialDerivative;
	else pddata = pdd;
	
	itd->NPARTIALS = pddata[i++];
	
	itd->PARTIAL_LIST = (int*) calloc(itd->NPARTIALS, sizeof(int));
	for(j=0; j<itd->NPARTIALS; j++) itd->PARTIAL_LIST[j] = pddata[i++];
	
	itd->NDER = pddata[i++]; 
	
	num = pddata[i++];
	itd->PREV_ACCUM = (long*) calloc(num, sizeof(long));
	for(j=0; j<num; j++) itd->PREV_ACCUM[j] = pddata[i++];

	num = pddata[i++];
	itd->PREV_COEF = (long*) calloc(num, sizeof(long));
	for(j=0; j<num; j++) itd->PREV_COEF[j] = pddata[i++];

	num = pddata[i++];
	itd->PREV_VI = (long*) calloc(num, sizeof(long));
	for(j=0; j<num; j++) itd->PREV_VI[j] = pddata[i++];
	
	num = pddata[i++];
	itd->PREV_IV = (long*) calloc(num, sizeof(long));
	for(j=0; j<num; j++) itd->PREV_IV[j] = pddata[i++];

	num = pddata[i++];
	itd->PREVSTAR_ACCUM = (long*) calloc(num, sizeof(long));
	for(j=0; j<num; j++) itd->PREVSTAR_ACCUM[j] = pddata[i++];
	
	num = pddata[i++];
	itd->PREVSTAR_COEF = (long*) calloc(num, sizeof(long));
	for(j=0; j<num; j++) itd->PREVSTAR_COEF[j] = pddata[i++];
	
	num = pddata[i++];
	itd->PREVSTAR_VI = (long*) calloc(num, sizeof(long));
	for(j=0; j<num; j++) itd->PREVSTAR_VI[j] = pddata[i++];
	
	num = pddata[i++];
	itd->PREVSTAR_IV = (long*) calloc(num, sizeof(long));
	for(j=0; j<num; j++) itd->PREVSTAR_IV[j] = pddata[i++];
}

void delete_iteration_parameters(iteration_data *itd)
{
	free(itd->PREV_ACCUM);
	free(itd->PREV_COEF);
	free(itd->PREV_VI);
	free(itd->PREV_IV);
	free(itd->PREVSTAR_ACCUM);
	free(itd->PREVSTAR_COEF);
	free(itd->PREVSTAR_VI);
	free(itd->PREVSTAR_IV);
}

void set_links(iteration_data *itd, int l, int *flst)
{
	itd->NLINKS			= l;
	itd->FUNCTION_LIST	= flst;
}


void check_iteration_data_parameters(int t, int pa, int pb)
{
	if (pa != pb) {
		switch (t) {
			case 0:
				printf("Different number of variables in driver and ODE function\n");
				break;
			case 1:
				printf("Different number of parameters in driver and ODE function\n\n");
				break;
			case 2:
				printf("Different number of functions in driver and ODE function\n\n");
				break;
			default:
				break;
		}
		exit(0);
	}
}


/************************************************************************/
