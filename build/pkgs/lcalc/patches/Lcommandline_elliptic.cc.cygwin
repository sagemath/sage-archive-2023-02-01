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

#include "Lcommandline_elliptic.h"


//returns 0 if initialization goes smoothly, 1 if there's an error
// given the elliptic curve
// y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6
// i.e. computes, the sign, the conductor, and first N_terms dirichlet coefficients
// of the corresponding L-function. The nth dirichlet coefficient is normalized by
// sqrt(n)
//

#ifdef INCLUDE_PARI
void initialize_new_L(char *a1, char *a2, char *a3, char *a4, char *a6, int N_terms)
{

    // basic data for the L-function (see the class L_function for full comments)
    int what_type;
    Long Period;
    Double q;
    Complex w;
    int A;
    Double *g;
    Complex *l;
    int n_poles;
    Complex *p;
    Complex *r;

    current_L_type=2; //the normalized dirichlet coeffs are real

    what_type=2; //i.e. eliiptic curve L-functions are cusp form L-function

    Period=0;

    A=1;

    //GAMMA factor is GAMMA(s+1/2)
    g=new Double[2];
    l=new Complex[2];
    g[1]=1.;
    l[1]=.5;

    Double * coeff;
    coeff = new Double[N_terms+1];

    data_E(a1,a2,a3,a4,a6,N_terms,coeff);
    //coeff[n], if n > 1, is the nth dirichlet coefficient nomalized by sqrt(n).
    //coeff[0] comes back with the sign of the functional equation
    //coeff[1] comes back with the conductor of E. We then set it to one.


    q=sqrt(coeff[1])/(2*Pi);
    coeff[1]=1.;

    w=coeff[0];

    n_poles=0; //no poles
    p = new Complex[1];
    r = new Complex[1];

    Double_L=L_function<Double>("Elliptic curve",what_type,N_terms,coeff,Period,q,w,A,g,l,n_poles,p,r);

    delete [] g;
    delete [] l;
    delete [] p;
    delete [] r;
    delete [] coeff;


}


// SAGE -- used below -- needed for Cygwin.
#ifndef llrint
inline long long int llrint (double x)
{
    long long int llrintres;
    asm
    ("fistpll %0"
    : "=m" (llrintres) : "t" (x) : "st");
    return llrintres;
}
#endif


// given the elliptic curve
// y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6
// i.e. computes, the sign, the conductor, and first N_terms dirichlet coefficients
// of the corresponding L-function. The nth dirichlet coefficient is normalized by
// sqrt(n)
void data_E(char *a1, char *a2, char *a3, char *a4, char *a6, int N_terms,Double * coeff)
{


    int sign; //sign stores the sign of the functional equation
    Long p; //denotes a prime
    Long m,n;
    Double x,r,tmp,tmp2;
    Long conductor; //the conductor of the elliptic curve

    GEN y, F, E, C;


    y = cgeti(64);

    C = cgetg(4, t_VEC);


    //if a1 etc are integers, we can use gaffsg to
    //assign F[1] etc. However, I am treating a1 etc as character
    //strings to allow for larger integers, and therefore use gaffect


    F = cgetg(6, t_VEC);
    F[1] = lgeti(BIGDEFAULTPREC);
    F[2] = lgeti(BIGDEFAULTPREC);
    F[3] = lgeti(BIGDEFAULTPREC);
    F[4] = lgeti(BIGDEFAULTPREC);
    F[5] = lgeti(BIGDEFAULTPREC);

    //gaffsg(a1,(GEN) F[1]);
    //gaffsg(a2,(GEN) F[2]);
    //gaffsg(a3,(GEN) F[3]);
    //gaffsg(a4,(GEN) F[4]);
    //gaffsg(a6,(GEN) F[5]);

    gaffect(strtoGEN(a1), (GEN) F[1]);
    gaffect(strtoGEN(a2), (GEN) F[2]);
    gaffect(strtoGEN(a3), (GEN) F[3]);
    gaffect(strtoGEN(a4), (GEN) F[4]);
    gaffect(strtoGEN(a6), (GEN) F[5]);

    E = initell(F,BIGDEFAULTPREC);

    C=globalreduction(E);

    x=gtodouble((GEN) C[1]);


    //if(x<1e18) conductor=Long(x+.1);
    if(x<Double(1.*my_LLONG_MAX)) conductor=Long(x+.1);

    else{
        cout << "conductor equals " << x << " and is too large" << endl;
        exit(1);
    }


    gaffsg(1, (GEN) y);
    sign = ellrootno(E,y); //sign of the functional equation


    for(n=1;n<=N_terms;n++) coeff[n]=1.;

    n=2;
    do{
        if(isprime(n)){

            p=n;
            gaffsg(p,y);
            coeff[p] = Double(1.*llrint(gtodouble(apell(E,y))))/sqrt(Double(1.*p));
            //coeff[p] = Double(1.*Long(gtodouble(apell(E,y))+.1))/sqrt(Double(1.*p));

            if(gtolong(gmod((GEN) E[12],(GEN) y))==0) // if p|discriminant, i.e. bad reduction
            {

               tmp=coeff[p];
               r=tmp*tmp;
               x=Double(1.*p)*Double(1.*p);
               m=Long(x+.1);
               if(m<=N_terms)
               do{
                   coeff[m]=coeff[m]*r;
                   x=x*p;
                   m=Long(x+.1);
                   r=r*tmp;
               }while(m<=N_terms);

            }

            else{ // a(p^(j+1)) = a(p)a(p^j) - p a(p^(j-1)) so normalizng by sqrt
                  // gives coeff[p^(j+1)] = coeff[p] coeff[p^j]-coeff[p^(j-1)]

               x=Double(1.*p)*Double(1.*p);
               m=Long(x+.1);
               if(m<=N_terms)
               do{
                   tmp2=0;
                   coeff[m]=coeff[m]*(coeff[p]*coeff[m/p]-coeff[m/(p*p)]);
                   x=x*p;
                   m=Long(x+.1);
                   r=r*tmp;
               }while(m<=N_terms);
            }
        }
        else{

            p=1;
            do{
                p++;
            }while(n%p!=0);

            m=p;
            do{
                m=m*p;
            }while(n%m==0&&m<n);

            if(n%m!=0) m=m/p;

            coeff[n]=coeff[m]*coeff[n/m];
        }
//if(n%10000==1) cout << n << endl;
        n++;

    }while(n<=N_terms);

    coeff[0]=1.*sign; coeff[1]=1.*conductor; //the first two spots are used
                                       //for the sign and conductor. The rest are
                                       //used for the Dirichlet coefficients.

}

#endif //ifdef INCLUDE_PARI

void compute_rank(){
    switch(current_L_type)
    {
        case 1:
            int_L.compute_rank(true);
            break;
        case 2:
            Double_L.compute_rank(true);
            break;
        case 3:
            Complex_L.compute_rank(true);
            break;
    }
}

void verify_rank(int rank){
    switch(current_L_type)
    {
        case 1:
            int_L.verify_rank(rank);
            break;
        case 2:
            Double_L.verify_rank(rank);
            break;
        case 3:
            Complex_L.verify_rank(rank);
            break;
    }
}

