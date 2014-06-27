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

#include "doubPOL.h"

int compensated_horner = 0;

/******************************************************
 Calcula la derivada del polinomio pol
 *****************************************************/
void dp_pol_derivative(double *pol, int grado, double *dpol)
{
	int i;
	for(i = 1; i <= grado ; i++)  dpol[i-1] = i*pol[i];
}

/******************************************************
 Horner TIDES
 *****************************************************/
void dp_horner(double *pol, int grado, double t, double *eval)
{
	if (compensated_horner) 
		dp_compensated_horner(pol, grado, t, eval);
	else 
		dp_classic_horner(pol, grado, t, eval);
}

void dp_horner_der(double *pol, int grado, double t, double *eval)
{
	double dpol[grado];
	int i; 
	for(i=1; i<=grado; i++) dpol[i-1] = i*pol[i];
	dp_horner(dpol, grado-1, t, eval);
}

/******************************************************
 Horner clasico
 *****************************************************/

void dp_classic_horner(double *pol, int grado, double t, double *eval)
{
	int i;
	double vpol= 0;
	for(i = grado; i >= 0 ; i--)  vpol = vpol*t + pol[i];
	*eval = vpol;
}

/******************************************************
 Horner compensado
 *****************************************************/


void dp_twosum(double a, double b, double *x, double *y)
{
	volatile double z;
	*x = a+b;
	z = (*x)-a;
	*y = (a -((*x)-z))+(b-z);
}

void dp_split_real(double a, double *ah, double *al)
{
	volatile double z; 
	double factor;
	factor = exp2(27) + 1;
	z = a*factor;
	*ah = z -(z-a);
	*al = a-(*ah);
}

void dp_twoproduct(double a, double b, double *x, double *y)
{
	double ah,al,bh,bl;
	dp_split_real(a, &ah, &al);
	dp_split_real(b, &bh, &bl);
	*x = a*b;
	*y = al*bl - ((((*x) - ah*bh) - al*bh) - ah*bl);
}

void dp_compensated_horner(double *pol, int grado, double t, double *eval)
{
	double exceso, p, pe,se;
	*eval = pol[grado];
	exceso = 0.;
	int i; 
	for(i= grado-1; i>=0; i--) {
		dp_twoproduct(*eval, t, &p, &pe);
		dp_twosum(p, pol[i], eval, &se);
		exceso = exceso*t+pe+se;
	}
	*eval +=exceso;	
}

/******************************************************
 Evalua el polinomio pol y la derivada  en el punto t
 *****************************************************/
void dp_ch_evaluate_poly_der(double *pol, int grado, double t, double *vpol, double *vder)
{
	double dpol[grado];
	dp_pol_derivative(pol, grado, dpol);
	dp_compensated_horner(pol, grado, t, vpol);
	dp_compensated_horner(dpol, grado-1, t, vder);
}

#define ITMAX 1000


/******************************************************
 Encuentra un cero en un intervalo.
 Debe haber como máximo un unico cero en ese intervalo y ese cero
 debe ser simple (cambio de signo) o doble 
 (minimo o maximo proximo al cero).
 Devuelve 1 si lo encuentra y cero si no lo encuentra.
 *****************************************************/
int dp_one_zero(double *pol, int grado, double ta, double tb, double tol, double *raiz)
{
	int dgrado = grado-1, iter;
	double dpol[grado];
	dp_pol_derivative(pol, grado, dpol);
	double fa,fb;
	dp_compensated_horner(pol, grado, ta, &fa);
	dp_compensated_horner(pol, grado, tb, &fb);
 	if (fa*fb <= 0.) {
		iter = dp_brent(pol,grado,ta,tb,tol, raiz);
//		iter = dp_rtsafe(pol,grado,ta,tb,tol, raiz);
		if(iter > ITMAX) 
			printf("Bad approximation of the zero   %.15le (Maximum number of iterations exceeded) ", *raiz);
		return 1;
	}
	dp_compensated_horner(dpol, dgrado, ta, &fa);
	dp_compensated_horner(dpol, dgrado, tb, &fb);
 	if (fa*fb <= 0.) {
		iter = dp_brent(dpol,dgrado,ta,tb,tol, raiz); 
//		iter = dp_rtsafe(dpol,dgrado,ta,tb,tol, raiz); 
		dp_compensated_horner(pol, grado, *raiz, &fa);
		if(fabs(fa) > tol) return 0;
		if(iter > ITMAX) 
			printf("Bad approximation of the zero   %.15le (Maximum number of iterations exceeded) ", *raiz);
		return 1;
	}
	return 0;
}
/******************************************************
 Encuentra un extremo en un intervalo.
 Debe haber como máximo un unico extremo en ese intervalo 
 y este debe corresponder a un cero simple de la derivada.
 Devuelve -1 si se encuentra un minimo
 Devuelve  1 si se encuentra un maximo
 Devuelve  0 si se no se encuentra nada
 *****************************************************/

int  dp_one_extremum(double *pol, int grado, double ta, double tb, double tol, double *extr)
{
	int dgrado = grado-1, ddgrado = grado-2, iter, sal=0;
	double dpol[grado], ddpol[dgrado];
	dp_pol_derivative(pol, grado, dpol);
	dp_pol_derivative(dpol, dgrado, ddpol);
	double fa,fb;
	dp_compensated_horner(dpol, dgrado, ta, &fa);
	dp_compensated_horner(dpol, dgrado, tb, &fb);
 	if (fa*fb <= 0.) {
		iter = dp_brent(dpol,dgrado,ta,tb,tol, extr); 
//		iter = dp_rtsafe(dpol,dgrado,ta,tb,tol, extr); 
		dp_compensated_horner(ddpol, ddgrado, *extr, &fa);
		if(fa < 0) sal =  1;
		else if(fa > 0) sal = -1;
		if(iter > ITMAX) 
			printf("Bad approximation of the extremum  %.15le (Maximum number of iterations exceeded) ", *extr);
	}
	return sal;
}

int  dp_one_maximum(double *pol, int grado, double ta, double tb, double tol, double *extr)
{
	int caso;
	caso = dp_one_extremum(pol, grado, ta, tb, tol, extr);
	if(caso == -1) caso = 0;
	return caso;
}
int  dp_one_minimum(double *pol, int grado, double ta, double tb, double tol, double *extr)
{
	int caso;
	caso = dp_one_extremum(pol, grado, ta, tb, tol, extr);
	if(caso == -1 ) caso = 1;
	else caso = 0;
	return caso;
}


/******************************************************
 Metodo de Brent para encontrar el cero (corte) de un
 polinomio pol en un intervalo con signo distinto de 
 la  funcion en los extremos del intervalo. 
 No se comprueba que solo hay una raiz y que es simple.
 *****************************************************/

int dp_brent(double *pol, int grado, double ta, double tb, double tol, double *raiz)
{
	int iter;
	double a=ta,b=tb,c=tb,d,e,min1,min2;
	double fa,fb,fc,p,q,r,s,tol1,xm;
	
	dp_compensated_horner(pol,grado,a,&fa);
	dp_compensated_horner(pol,grado,b,&fb);
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) { 
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
 		tol1 = 2.0*tol*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0)  break;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else {
			if(xm >= 0.0) b += fabs(tol1);
			else b -= fabs(tol1);
		}
		dp_compensated_horner(pol,grado,b,&fb);
	}
	*raiz = b;
	return iter;
}

/******************************************************
 Metodo combinado con Newton para encontrar el cero (corte) de un
 polinomio pol en un intervalo con signo distinto de 
 la  funcion en los extremos del intervalo. No se comprueba, 
 hay que comprobarlo antes
 *****************************************************/


int dp_rtsafe(double *pol, int grado, double x1, double x2, double tol, double *raiz)
{
	int j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;
	
	dp_ch_evaluate_poly_der(pol,grado,x1,&fl,&df);
	dp_ch_evaluate_poly_der(pol,grado,x2,&fh,&df);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		printf("Root must be bracketed in rtsafe\n");
	if (fl == 0.0) {
		*raiz = x1;
		return 0;
	}
	if (fh == 0.0) {
		*raiz = x2;
		return 0;
	}
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;	
	dp_ch_evaluate_poly_der(pol,grado,rts,&f,&df);
	for (j=1;j<=ITMAX;j++) { 
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx = (xh/xl-1)*(xl/2) ;
			rts= xl+dx ;
			if (xl == rts) {
				*raiz = rts;
				return j;
			}
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) {
				*raiz = rts;
				return j;
			}
		} 
		if (fabs(dx) < tol) {
			*raiz = rts;
			return j;
		}		
		dp_ch_evaluate_poly_der(pol,grado,rts,&f,&df);
		if (f < 0.0)  
			xl=rts; 
		else 
			xh=rts;
	}
	*raiz = rts;
	return j;
}

