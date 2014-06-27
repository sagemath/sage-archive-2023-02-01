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

#include "mpfrPOL.h"

extern int compensated_horner;
extern int BINARY_PRECISION;

/******************************************************
 Calcula la derivada del polinomio pol
 *****************************************************/
void mp_pol_derivative(mpfr_t *pol, int grado, mpfr_t *dpol)
{
	int i;
	for(i = 1; i <= grado ; i++)  
		mpfrts_mul_i (&dpol[i-1], pol[i], i);
}

/******************************************************
 Horner TIDES
 *****************************************************/

void mp_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval)
{
	if (compensated_horner) 
		mp_compensated_horner(pol, grado, t, eval);
	else
		mp_classic_horner(pol, grado, t, eval);
}

void mp_horner_der(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval)
{
	mpfr_t dpol[grado];
	int i; 
	for(i=1; i<=grado; i++) {
		mpfrts_init(&dpol[i-1]);
		mpfrts_mul_i(&dpol[i-1], pol[i], i);
	}
	mp_horner(dpol, grado-1, t, eval);
	for(i=1; i<=grado; i++) mpfr_clear(dpol[i-1]);
	mpfr_free_cache ();
}

/******************************************************
 Horner clasico
 *****************************************************/

void mp_classic_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval)
{
	int i;
	mpfr_t aux;
	mpfrts_init (&aux);
	mpfrts_set(eval, aux);
	for(i = grado; i >= 0 ; i--)  {
		mpfrts_mul (&aux, t, *eval);
		mpfrts_add (eval, aux, pol[i]);
	}
	mpfr_clear (aux);
	mpfr_free_cache ();
}

/******************************************************
 Horner compensado
 *****************************************************/

void mp_twosum(mpfr_t a, mpfr_t b, mpfr_t *x, mpfr_t *y)
{
	mpfr_t z, aux;
	mpfrts_init (&z);
	mpfrts_init (&aux);
	mpfrts_add (x, a, b);
	mpfrts_sub (&z, *x, a);
	mpfrts_sub (&aux, b, z);
	mpfrts_sub (y, *x, z);
	mpfrts_sub (y, a, *y);
	mpfrts_add (y, *y, aux);
	mpfr_clear (z);
	mpfr_clear (aux);
	mpfr_free_cache ();
}

void mp_split_real(mpfr_t a, mpfr_t *ah, mpfr_t *al)
{
	mpfr_t z, factor;
	mpfrts_init (&z);
	mpfrts_init (&factor);
	mpfrts_set_i(&factor, ceil(BINARY_PRECISION/2.));
	mpfrts_exp2(&factor, factor);
	mpfrts_mul(&z, a, factor);
	mpfrts_sub(ah, z, a);
	mpfrts_sub(ah, z, *ah);
	mpfrts_sub(al, a, *ah);
	mpfr_clear (z);
	mpfr_clear (factor);
	mpfr_free_cache ();
}

void mp_twoproduct(mpfr_t a, mpfr_t b, mpfr_t *x, mpfr_t *y)
{
	mpfr_t ah,al,bh,bl, aux;
	mpfrts_init (&ah);
	mpfrts_init (&al);
	mpfrts_init (&bh);
	mpfrts_init (&bl);
	mpfrts_init (&aux);
	mp_split_real(a, &ah, &al);
	mp_split_real(b, &bh, &bl);
	mpfrts_mul(x, a, b);
	
	mpfrts_mul(&aux, ah, bh);
	mpfrts_sub(y, *x, aux);
	mpfrts_mul(&aux, al, bh);
	mpfrts_sub(y, *y, aux);
	mpfrts_mul(&aux, ah, bl);
	mpfrts_sub(y, *y, aux);
	mpfrts_mul(&aux, al, bl);
	mpfrts_sub(y, aux, *y);
	mpfr_clear (ah);
	mpfr_clear (al);
	mpfr_clear (bh);
	mpfr_clear (bl);
	mpfr_clear (aux);
	mpfr_free_cache ();
	
	
}

void mp_compensated_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval)
{
	mpfr_t exceso, p, pe,se;
	mpfrts_init (&exceso);
	mpfrts_init (&p);
	mpfrts_init (&pe);
	mpfrts_init (&se);
	
	mpfrts_set(eval, pol[grado]);
	
	int i; 
	for(i= grado-1; i>=0; i--) {
		mp_twoproduct(*eval, t, &p, &pe);
		mp_twosum(p, pol[i], eval, &se);
		mpfrts_mul(&exceso, exceso, t);
		mpfrts_add(&exceso, exceso, pe);
		mpfrts_add(&exceso, exceso, se);
	}
	mpfrts_add(eval, *eval, exceso);
	mpfr_clear (exceso);
	mpfr_clear (p);
	mpfr_clear (pe);
	mpfr_clear (se);
	mpfr_free_cache ();
}

/******************************************************
 Evalua el polinomio pol y la derivada  en el punto t
 *****************************************************/
void mp_ch_evaluate_poly_der(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *vpol, mpfr_t *vder)
{
	mpfr_t dpol[grado];
	int i; 
	for(i=1; i<=grado; i++) {
		mpfrts_init(&dpol[i-1]);
		mpfrts_mul_i(&dpol[i-1], pol[i], i);
	}
	mp_compensated_horner(pol, grado, t, vpol);
	mp_compensated_horner(dpol, grado-1, t, vder);
	for(i=1; i<=grado; i++) mpfr_clear(dpol[i-1]);
	mpfr_free_cache ();
}

#define MPITMAX 1000


/******************************************************
 Encuentra un cero en un intervalo.
 Debe haber como máximo un unico cero en ese intervalo y ese cero
 debe ser simple (cambio de signo) o doble 
 (minimo o maximo proximo al cero).
 Devuelve 1 si lo encuentra y cero si no lo encuentra.
 *****************************************************/
int mp_one_zero(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *raiz)
{
	int dgrado = grado-1, iter, i;
	mpfr_t dpol[grado];
	for(i=0; i<grado; i++) mpfrts_init(&dpol[i]);
	mp_pol_derivative(pol, grado, dpol);
	mpfr_t fa,fb, temp, cero;
	mpfrts_init(&fa);
	mpfrts_init(&fb);
	mpfrts_init(&temp);
	mpfrts_init(&cero);
	mp_compensated_horner(pol, grado, ta, &fa);
	mp_compensated_horner(pol, grado, tb, &fb);
	mpfrts_mul(&temp, fa, fb);
 	if (mpfrts_lessequal(temp, cero)) {
		iter = mp_brent(pol,grado,ta,tb,tol, raiz);
//		iter = mp_rtsafe(pol,grado,ta,tb,tol, raiz);
		if(iter > MPITMAX) 
			mpfrts_write("Maximum number of iterations exceeded:  Bad approximation of the zero \n", *raiz);			
		for(i=0; i<grado; i++) mpfr_clear(dpol[i]);
		mpfr_clear(fa);
		mpfr_clear(fb);
		mpfr_clear(temp);
		mpfr_clear(cero);
		mpfr_free_cache ();
		return 1;
	}
	mp_compensated_horner(dpol, dgrado, ta, &fa);
	mp_compensated_horner(dpol, dgrado, tb, &fb);
	mpfrts_mul(&temp, fa, fb);
 	if (mpfrts_lessequal(temp, cero)) {
		iter = mp_brent(dpol,dgrado,ta,tb,tol, raiz); 
//		iter = mp_rtsafe(dpol,dgrado,ta,tb,tol, raiz); 
		mp_compensated_horner(pol, grado, *raiz, &fa);
		mpfrts_abs(&temp, fa);
		if(mpfrts_greater(temp,tol)) return 0;
		if(iter > MPITMAX) 
			mpfrts_write("Maximum number of iterations exceeded:  Bad approximation of the zero \n", *raiz);			
		for(i=0; i<grado; i++) mpfr_clear(dpol[i]);
		mpfr_clear(fa);
		mpfr_clear(fb);
		mpfr_clear(temp);
		mpfr_clear(cero);
		mpfr_free_cache ();
		return 1;
	}
	for(i=0; i<grado; i++) mpfr_clear(dpol[i]);
	mpfr_clear(fa);
	mpfr_clear(fb);
	mpfr_clear(temp);
	mpfr_clear(cero);
	mpfr_free_cache ();
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

int  mp_one_extremum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr)
{
	int dgrado = grado-1, ddgrado = grado-2, iter, i, sal=0;
	mpfr_t dpol[grado], ddpol[dgrado];
	for(i=0; i<grado; i++) mpfrts_init(&dpol[i]);
	for(i=0; i<dgrado; i++) mpfrts_init(&ddpol[i]);
	mp_pol_derivative(pol, grado, dpol);
	mp_pol_derivative(dpol, dgrado, ddpol);
	mpfr_t fa,fb, temp, cero;;
	mpfrts_init(&fa);
	mpfrts_init(&fb);
	mpfrts_init(&temp);
	mpfrts_init(&cero);
	mp_compensated_horner(dpol, dgrado, ta, &fa);
	mp_compensated_horner(dpol, dgrado, tb, &fb);
	mpfrts_mul(&temp, fa, fb);
 	if (mpfrts_lessequal(temp, cero)) {
		iter = mp_brent(dpol,dgrado,ta,tb,tol, extr); 
//		iter = mp_rtsafe(dpol,dgrado,ta,tb,tol, extr); 
		mp_compensated_horner(ddpol, ddgrado, *extr, &fa);
		if(mpfrts_less(fa,cero)) sal =  1;
		else if(mpfrts_greater(fa,cero)) sal = -1;
		if(iter > MPITMAX) 
			mpfrts_write("Maximum number of iterations exceeded:  Bad approximation of the extremum \n", * extr);			
	}
	for(i=0; i<grado; i++) mpfr_clear(dpol[i]);
	for(i=0; i<dgrado; i++) mpfr_clear(ddpol[i]);
	mpfr_clear(fa);
	mpfr_clear(fb);
	mpfr_clear(temp);
	mpfr_clear(cero);
	mpfr_free_cache ();
	return sal;
}

int  mp_one_maximum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr)
{
	int caso;
	caso = mp_one_extremum(pol, grado, ta, tb, tol, extr);
	if(caso == -1) caso = 0;
	return caso;
}
int  mp_one_minimum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr)
{
	int caso;
	caso = mp_one_extremum(pol, grado, ta, tb, tol, extr);
	if(caso == 1) caso = 0;
	else caso = 1;
	return caso;
}

/******************************************************
 Metodo de Brent para encontrar el cero (corte) de un
 polinomio pol en un intervalo con signo distinto de 
 la  funcion en los extremos del intervalo. 
 No se comprueba que solo hay una raiz y que es simple.
 *****************************************************/

int mp_brent(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *raiz)
{
	int iter;
	mpfr_t a,b,c,d,e,min1,min2,cero;
	mpfr_t fa,fb,fc,p,q,r,s,tol1,xm;
	mpfr_t temp, tempb, tempc;
	mpfrts_init(&a);
	mpfrts_set(&a, ta);
	mpfrts_init(&b);
	mpfrts_set(&b, tb);
	mpfrts_init(&c);
	mpfrts_set(&c, tb);
	mpfrts_init(&d);
	mpfrts_init(&e);
	mpfrts_init(&min1);
	mpfrts_init(&min2);
	mpfrts_init(&cero);
	mpfrts_init(&fa);
	mpfrts_init(&fb);
	mpfrts_init(&fc);
	mpfrts_init(&p);
	mpfrts_init(&q);
	mpfrts_init(&r);
	mpfrts_init(&s);
	mpfrts_init(&tol1);
	mpfrts_init(&xm);
	mpfrts_init(&temp);
	mpfrts_init(&tempb);
	mpfrts_init(&tempc);
	
	mp_compensated_horner(pol,grado,a,&fa);
	mp_compensated_horner(pol,grado,b,&fb);
	mpfrts_set(&fc, fb);
	for (iter=1;iter<=MPITMAX;iter++) {
		mpfrts_mul(&temp, fb, fc);
		if (mpfrts_greater(temp, cero) ) { 
			mpfrts_set(&c, a);
			mpfrts_set(&fc, fa);
			mpfrts_sub(&d, b, a);
			mpfrts_set(&e, d);
		}
		mpfrts_abs(&temp, fc);
		mpfrts_abs(&tempb, fb);
		if (mpfrts_less(temp, tempb)) {
			mpfrts_set(&a, b);
			mpfrts_set(&b, c);
			mpfrts_set(&c, a);
			mpfrts_set(&fa, fb);
			mpfrts_set(&fb, fc);
			mpfrts_set(&fc, fa);
		}
		mpfrts_div_i(&tempb, tol, 2);
		mpfrts_abs(&temp, b);
		mpfrts_mul(&temp, temp, tol);
		mpfrts_mul_i(&temp, temp, 2);
		mpfrts_add(&tol1, temp, tempb);
		mpfrts_sub(&xm, c, b);
		mpfrts_div_i(&xm, xm, 2);
		mpfrts_abs(&temp, xm);
		if(mpfrts_lessequal(temp, tol1) || mpfrts_equal(fb, cero)) break;
		mpfrts_abs(&temp, e);
		mpfrts_abs(&tempb, fa);
		mpfrts_abs(&tempc, fb);
		if (mpfrts_greaterequal(temp, tol1) && mpfrts_greater(tempb, tempc)) {
			mpfrts_div(&s, fb, fa);
			if (mpfrts_equal(a, c)) {
				mpfrts_mul(&p, xm, s);
				mpfrts_mul_i(&p, p, 2);
				mpfrts_i_sub(&q, 1, s);
			} else {
				mpfrts_div(&q, fa, fc);
				mpfrts_div(&r, fb, fc);
				mpfrts_sub_i(&temp, r, 1);
				mpfrts_sub(&tempb, b, a);
				mpfrts_sub(&tempc, q, r);
				mpfrts_mul(&tempc, tempc, q);
				mpfrts_mul(&tempc, tempc, xm);
				mpfrts_mul_i(&tempc, tempc, 2);
				mpfrts_sub(&tempc, tempc, tempb);
				mpfrts_mul(&tempc, tempc, temp);
				mpfrts_mul(&p, s, tempc);
				mpfrts_sub_i(&tempb, s, 1);
				mpfrts_sub_i(&tempc, q, 1);
				mpfrts_mul(&q, temp, tempb);
				mpfrts_mul(&q, q, tempc);
			}
			if(mpfrts_greater(p, cero)) mpfrts_neg(&q, q);
			mpfrts_abs(&p, p);
			mpfrts_mul(&temp, tol1, q);
			mpfrts_abs(&temp, temp);
			mpfrts_mul(&tempb, xm, q);
			mpfrts_mul_i(&tempb, tempb, 3);
			mpfrts_sub(&min1, tempb, temp);
			mpfrts_mul(&min2, e, q);
			mpfrts_abs(&min2, min2);
			mpfrts_mul_i(&temp, p, 2);
			if(mpfrts_less(min1, min2) ) mpfrts_set(&tempb, min1);
			else mpfrts_set(&tempb, min2);
			if (mpfrts_less(temp, tempb)) {
				mpfrts_set(&e, d);
				mpfrts_div(&d, p, q);
			} else {
				mpfrts_set(&d, xm);
				mpfrts_set(&e, d);
			}
		} else {
			mpfrts_set(&d, xm);
			mpfrts_set(&e, d);
		}
		mpfrts_set(&a, b);
		mpfrts_set(&fa, fb);
		mpfrts_abs(&temp, d);
		if (mpfrts_greater(temp, tol1))
			mpfrts_add(&b, b, d);
		else {
			mpfrts_abs(&temp, tol1);
			if(mpfrts_greaterequal(xm, cero)) mpfrts_add(&b, b, temp);
			else mpfrts_sub(&b, b, temp);
		}
		mp_compensated_horner(pol,grado,b,&fb);
	}
	mpfrts_set(raiz, b);
	
	mpfr_clear(a);
	mpfr_clear(b);
	mpfr_clear(c);
	mpfr_clear(d);
	mpfr_clear(e);
	mpfr_clear(min1);
	mpfr_clear(min2);
	mpfr_clear(cero);
	mpfr_clear(fa);
	mpfr_clear(fb);
	mpfr_clear(fc);
	mpfr_clear(p);
	mpfr_clear(q);
	mpfr_clear(r);
	mpfr_clear(s);
	mpfr_clear(tol1);
	mpfr_clear(xm);
	mpfr_clear(temp);
	mpfr_clear(tempb);
	mpfr_clear(tempc);
	mpfr_free_cache ();
	
	return iter;
}

/******************************************************
 Metodo combinado con Newton para encontrar el cero (corte) de un
 polinomio pol en un intervalo con signo distinto de 
 la  funcion en los extremos del intervalo. No se comprueba, 
 hay que comprobarlo antes
 *****************************************************/


int mp_rtsafe(mpfr_t *pol, int grado, mpfr_t x1, mpfr_t x2, mpfr_t tol, mpfr_t *raiz)
{
	int j;
	mpfr_t df,dx,dxold,f,fh,fl;
	mpfr_t xh,xl,rts, cero;
	mpfr_t temp, tempb, tempc;
	mpfrts_init(&df);
	mpfrts_init(&dx);
	mpfrts_init(&dxold);
	mpfrts_init(&f);
	mpfrts_init(&fh);
	mpfrts_init(&fl);
	mpfrts_init(&xh);
	mpfrts_init(&xl);
	mpfrts_init(&rts);
	mpfrts_init(&cero);
	mpfrts_init(&temp);
	mpfrts_init(&tempb);
	mpfrts_init(&tempc);
	
	mp_ch_evaluate_poly_der(pol,grado,x1,&fl,&df);
	mp_ch_evaluate_poly_der(pol,grado,x2,&fh,&df);
	mpfrts_mul(&temp, fl, fh);
	if (mpfrts_greater(temp, cero))
		printf("Root must be bracketed in rtsafe\n");
	if (mpfrts_equal(fl, cero)) {
		mpfrts_set(raiz, x1);
		return 0;
	}
	if (mpfrts_equal(fh, cero)) {
		mpfrts_set(raiz, x2);
		return 0;
	}
	if (mpfrts_less(fl, cero)) {
		mpfrts_set(&xl, x1);
		mpfrts_set(&xh, x2);
	} else {
		mpfrts_set(&xh, x1);
		mpfrts_set(&xl, x2);
	}
	mpfrts_add(&rts, x1, x2);
	mpfrts_div_i(&rts, rts, 2);
	mpfrts_sub(&dxold, x2, x1);
	mpfrts_abs(&dxold, dxold);
	mpfrts_set(&dx, dxold);
	mp_ch_evaluate_poly_der(pol,grado,rts,&f,&df);
	for (j=1;j<=MPITMAX;j++) { 
		mpfrts_sub(&temp, rts, xh);
		mpfrts_mul(&temp, temp, df);
		mpfrts_sub(&temp, temp, f);
		mpfrts_sub(&tempb, rts, xl);
		mpfrts_mul(&tempb, tempb, df);
		mpfrts_sub(&tempb, tempb, f);
		mpfrts_mul(&temp, temp, tempb);
		mpfrts_mul_i(&tempb, f, 2);
		mpfrts_abs(&tempb, tempb);
		mpfrts_mul(&tempc, dxold, df);
		mpfrts_abs(&tempc, tempc);
		if (mpfrts_greater(temp, cero) || mpfrts_greater(tempb, tempc)) {
			mpfrts_set(&dxold, dx);
			mpfrts_div(&temp, xh, xl);
			mpfrts_sub_i(&temp, temp, 1);
			mpfrts_div_i(&tempb, xl, 2);
			mpfrts_mul(&dx, temp, tempb);
			mpfrts_add(&rts, xl, dx);
			if(mpfrts_equal(xl, rts) ){
				mpfrts_set(raiz, rts);
				return j;
			}
		} else {
			mpfrts_set(&dxold, dx);
			mpfrts_div(&dx, f, df);
			mpfrts_set(&temp, rts);
			mpfrts_sub(&rts, rts, dx);
			if(mpfrts_equal(temp, rts) ){
				mpfrts_set(raiz, rts);
				return j;
			}
		} 
		mpfrts_abs(&temp, dx);
		if (mpfrts_less(temp, tol)) {
			mpfrts_set(raiz, rts);
			return j;
		}		
		mp_ch_evaluate_poly_der(pol,grado,rts,&f,&df);
		if (mpfrts_less(f, cero)) mpfrts_set(&xl, rts);
		else  mpfrts_set(&xh, rts);
	}
	mpfrts_set(raiz, rts);
	mpfr_clear(df);
	mpfr_clear(dx);
	mpfr_clear(dxold);
	mpfr_clear(f);
	mpfr_clear(fh);
	mpfr_clear(fl);
	mpfr_clear(xh);
	mpfr_clear(xl);
	mpfr_clear(rts);
	mpfr_clear(cero);
	mpfr_clear(temp);
	mpfr_clear(tempb);
	mpfr_clear(tempc);
	mpfr_free_cache ();
	
	return j;
}



