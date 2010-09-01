/*

   Copyright (C) 2001,2002,2003,2004 Michael Rubinstein

   This file is part of the L-function package L.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   Check the License for details. You should have received a copy of it, along
   with the package; see the file 'COPYING'. If not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

*/


#ifndef Lcommandline_numbertheory_H
#define Lcommandline_numbertheory_H

#include <time.h>
#include <stdlib.h>     //for things like srand
#include <iostream>     //for input and output
#include <iomanip>      //for manipulating output such as setprecision
#include <math.h>
//#include <ctime>
#include <limits.h>     // for INT_MAX
#include "Lglobals.h"

using namespace std;

long long nextprime(long long n);
bool issquarefree(long long n);
bool isfunddiscriminant(long long n);
long long nextfunddiscriminant(long long n);
int simplified_jacobi(int n,int m);
int simplified_jacobi(long long n,long long m);
int my_kronecker(int n,int m);
int my_kronecker(long long n, long long m);
Double L_1_chi(long long d);
int class_number(long long d);
void ramanujan_tau(Double *c0, int N_terms);
long long gcd(long long a,long long b);
void factor(long long q, long long **factors);
long long power_mod_q(long long a, long long k,long long q);
int prim_root(long long p, int alpha);
void factor(long long q, long long **factors);
int characters();

bool isprime(long long n);  //by Kyle Wichert
bool RM(long long a, long long N);  //by Kyle Wichert
long long multmodN(long long a, long long b, long long N);  //by Kyle Wichert


template <class ttype>
Complex gauss_sum(ttype *chi,long long r,bool cnj=false)
{
    Complex SUM=0;
    if(cnj==true)
        for(int n=1;n<=r; n++) SUM=SUM+conj(chi[n])*exp(n*2*I*Pi/double(r));
    else
        for(int n=1;n<=r; n++) SUM=SUM+chi[n]*exp(n*2*I*Pi/double(r));

    return SUM;
}


#endif
