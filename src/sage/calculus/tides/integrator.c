/****************************************************************************
	This file has been created by SageTIDES (1.00)

	Copyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Marco, M. Rodriguez
	Grupo de Mecanica Espacial
	University of Zaragoza
	SPAIN

	http://gme.unizar.es/software/tides
	Contact: <tides@unizar.es>

	This file is part of TIDES.

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



#include "minc_tides.h"

void    mincseries(double t,double *v, double *p, double **XVAR,int ORDER, int MO)
{
	int VAR,PAR,TT,i,j, inext;

	VAR = 3;
	PAR = 1;
	TT = 7;
	double XX[TT+1][MO+1];

	for(j=0; j<=TT; j++)
		for(i=0; i<=ORDER; i++)
			XX[j][i] = 0.e0;
	XX[0][0] = t;
	XX[0][1] = 1.e0;
	for(i=1;i<=VAR;i++) {
		XX[i][0] = v[i-1];
	}

	for(i=0;i<ORDER;i++) {

		XX[4][i] = -1.00000000000000*XX[1][i];
		XX[5][i] = (i==0)? 1.09861228866811+XX[4][0] : XX[4][i];
		XX[6][i] = sin_mc(XX[0],XX[7], i);
		XX[7][i] = cos_mc(XX[0],XX[6], i);
		XX[1][i+1] = XX[2][i] / (i+1.0);
		XX[2][i+1] = XX[5][i] / (i+1.0);
		XX[3][i+1] = XX[7][i] / (i+1.0);
	}

	for(j=0; j<=VAR; j++)
		for(i=0; i<=ORDER; i++)
			XVAR[i][j] = XX[j][i];
}

