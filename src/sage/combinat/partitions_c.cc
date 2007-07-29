/*
    Program to compute the number of partitions of n.

    Author: Jonathan Bober
    Version: .1 (7/28/2007)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>

#include <cmath>

#include <iostream>
#include <iomanip>
#include <limits>

#include <mpfr.h>
#include <gmp.h>

const bool debug = false;
//const bool debug = true;

const bool debugs = false;

//const bool debugf = true;
const bool debugf = false;
//const bool debuga = true;
const bool debuga = false;
//const bool debugt = true;
const bool debugt = false;

using std::cout;
using std::endl;

const unsigned int min_prec = 53;
const unsigned int min_precision = 53;

long GCD(long a, long b);

unsigned int calc_precision(unsigned int n, unsigned int N);

//mpfr versions

mpz_t ztemp1, ztemp2, ztemp3;

mpfr_t mp_one_over_12, mp_one_over_24, mp_sqrt2, mp_sqrt3, mp_pi;
mpfr_t mp_A, mp_B, mp_C, mp_D;
mp_rnd_t round_mode = GMP_RNDN;

mpfr_t tempa1, tempa2, tempf1, tempf2, temps1, temps2, tempt1, tempt2;  // temp variables for different functions, with precision set and cleared by mp_set_precision

bool mp_vars_initialized = false;

mp_prec_t mp_precision;

void mp_init(unsigned int prec, unsigned int n);
//void mp_f_precompute(unsigned int n);
void mp_f(mpfr_t result, unsigned int k);
void mp_s(mpfr_t result, unsigned int j, unsigned int q);
void mp_a(mpfr_t result, unsigned int n, unsigned int k);
void mp_t(mpfr_t result, unsigned int n);

unsigned int compute_initial_precision(unsigned int n);     // computes the precision required to accurately compute p(n)
unsigned int compute_current_precision(unsigned int n, unsigned int N); // computed the precision required to
                                                                        // accurately compute the tail of the rademacher series
                                                                        // assuming that N terms have already been computed

double compute_remainder(unsigned int n, unsigned int N);   // Gives an upper bound on the error that occurs
                                                            // when only N terms of the Rademacher series have been
                                                            // computed. NOTE: should only be called when we already know
                                                            // that the error is small (ie, compute_current_precision returns
                                                            // a small number, eg, something < 32)

//low precision (double) versions of the functions:

double df_A, df_B, df_C, df_D;

void d_f_precompute(unsigned int n);
double d_f(unsigned int k);
double d_s(unsigned int h,unsigned int q);
double d_A(unsigned int n, unsigned int k);
double d_t(unsigned int n, unsigned int N);

unsigned int compute_initial_precision(unsigned int n) {
    // We just want to know how many bits we will need to
    // compute to get an accurate answer.

    // We know that
    //
    //          p(n) ~ exp(pi * sqrt(2n/3))/(4n sqrt(3)),
    //
    // so for now we are assuming that p(n) < exp(pi * sqrt(2n/3))/n,
    // so we need pi*sqrt(2n/3)/log(2) - log(n)/log(2) + EXTRA bits.
    //
    // EXTRA should depend on n, and should be something that ensures
    // that the TOTAL ERROR in all computations is < (something small).
    // This needs to be worked out carefully. EXTRA = log(n)/log(2) + 3
    // is probably good enough, and is convenient...
    //
    // but we really need:
    //
    //                  p(n) < something
    //
    // to be sure that we compute the correct answer

    /* I don't think that we need to be that careful about this.
     *
    mpfr_t bits;
    mpfr_t t1;                                          // we refer to t1 as temp1

    mpfr_init2(bits, 32);     // use 32 bits of precision to compute the approximation

    mpfr_init2(t1, 32);

    //    Error = exp( pi* sqrt(2n/3));

    mpfr_set_ui(t1, 2, round_mode);                     // temp1 = 2
    mpfr_mul_ui(t1, t1, n, round_mode);                 // temp1 = 2n
    mpfr_div_ui(t1, t1, 3, round_mode);                 // temp1 = 2n/3
    mpfr_sqrt(t1, t1, round_mode);                      // temp1 = sqrt(2n/3)

    mpfr_const_pi(bits, round_mode);                    // bits = pi

    mpfr_mul(bits, bits, t1, round_mode);               // bits = pi * sqrt(2n/3)

    mpfr_const_log2(t1, round_mode);                    // temp1 = log(2)
    mpfr_div(bits, bits, t1, round_mode);               // bits = pi * sqrt(2n/3)/log(2)
    */
    unsigned int result = (unsigned int)(ceil(3.1415926535897931 * sqrt(2.0 * double(n)/ 3.0) / log(2))) + 3;
    if(debug)
        cout << "Using initial precision of " << result << " bits." << endl;
    if(result > min_precision) {
        return result;
    }
    else return min_precision;

}

unsigned int compute_current_precision(unsigned int n, unsigned int N) {
    // we want to compute
    //      log(A/sqrt(N) + B*sqrt(N/(n-1))*sinh(C * sqrt(n) / N) / log(2)
    //
    //  maybe there is a better way...

    // cout << N;

    mpfr_t A, B, C;
    mpfr_init2(A, 32);
    mpfr_init2(B, 32);
    mpfr_init2(C, 32);

    mpfr_set_d(A,1.11431833485164,round_mode);
    mpfr_set_d(B,0.059238439175445,round_mode);
    mpfr_set_d(C,2.5650996603238,round_mode);

    mpfr_t error, t1, t2;
    mpfr_init2(error, 32);      // we shouldn't need much precision here since we just need the most significant bit
    mpfr_init2(t1, 32);
    mpfr_init2(t2, 32);

    mpfr_set(error, A, round_mode);      // error = A
    mpfr_sqrt_ui(t1, N, round_mode);        // t1 = sqrt(N)
    mpfr_div(error, error, t1, round_mode); // error = A/sqrt(N)


    mpfr_sqrt_ui(t1, n, round_mode);        // t1 = sqrt(n)
    mpfr_mul(t1, t1, C, round_mode);     // t1 = C * sqrt(n)
    mpfr_div_ui(t1, t1, N, round_mode);     // t1 = C * sqrt(n) / N
    mpfr_sinh(t1, t1, round_mode);          // t1 = sinh( ditto )
    mpfr_mul(t1, t1, B, round_mode);     // t1 = B * sinh( ditto )

    mpfr_set_ui(t2, N, round_mode);         // t2 = N
    mpfr_div_ui(t2, t2, n-1, round_mode);       // t2 = N/(n-1)
    mpfr_sqrt(t2, t2, round_mode);          // t2 = sqrt( ditto )

    mpfr_mul(t1, t1, t2, round_mode);       // t1 = B * sqrt(N/(n-1)) * sinh(C * sqrt(n)/N)

    mpfr_add(error, error, t1, round_mode); // error = (ERROR ESTIMATE)

    unsigned int p = mpfr_get_exp(error) + 3;   // I am almost certain that this does the right thing
                                                // (The 3 is for good luck.)

    p = p + (unsigned int)ceil(log(n)/log(2));

    if(debug) {
        cout << "Error seems to be: ";
        mpfr_out_str(stdout, 10, 0, error, round_mode);
        cout << endl;
        cout.flush();
        cout << "Switching to precision of " << p << " bits. " << endl;
    }

    mpfr_clear(error);
    mpfr_clear(t1);
    mpfr_clear(t2);
    mpfr_clear(A);
    mpfr_clear(B);
    mpfr_clear(C);




    if(p > min_precision) {
        return p;                           // don't want to return < 32
                                            // (for now, at least -- we can be more careful)
    }
    return min_precision;
}

double compute_remainder(unsigned int n, unsigned int N) {
    double A = 1.11431833485164;
    double B = 0.059238439175445;
    double C = 2.5650996603238;
    double result;
    result = A/sqrt(N) + B * sqrt(double(N)/double(n-1))*sinh(C * sqrt(double(n))/double(N));
    return result;
}


void mp_set_precision(unsigned int prec) {
    static bool init = false;
    if(init) {
        mpfr_clear(tempa1);
        mpfr_clear(tempa2);
        mpfr_clear(tempf1);
        mpfr_clear(tempf2);
        mpfr_clear(tempt1);
        mpfr_clear(tempt2);
        mpfr_clear(temps1);
        mpfr_clear(temps2);
    }
    mpfr_init2(tempa1, prec);
    mpfr_init2(tempa2, prec);
    mpfr_init2(tempf1, prec);
    mpfr_init2(tempf2, prec);
    mpfr_init2(tempt1, prec);
    mpfr_init2(tempt2, prec);
    mpfr_init2(temps1, prec);
    mpfr_init2(temps2, prec);

    init = true;
}

void mp_init(unsigned int prec, unsigned int n) {
    static bool init = false;
    mp_precision = prec;
    mp_prec_t p = mp_precision;

    if(init) {
        mpfr_clear(mp_one_over_12); mpfr_clear(mp_one_over_24); mpfr_clear(mp_sqrt2); mpfr_clear(mp_sqrt3); mpfr_clear(mp_pi);
        mpfr_clear(mp_A); mpfr_clear(mp_B); mpfr_clear(mp_C); mpfr_clear(mp_D);

        mpz_clear(ztemp1);
        mpz_clear(ztemp2);
        mpz_clear(ztemp3);
    }
    mpfr_init2(mp_one_over_12,p); mpfr_init2(mp_one_over_24,p); mpfr_init2(mp_sqrt2,p); mpfr_init2(mp_sqrt3,p); mpfr_init2(mp_pi,p);
    mpfr_init2(mp_A,p); mpfr_init2(mp_B,p); mpfr_init2(mp_C,p); mpfr_init2(mp_D,p);
    mpz_init(ztemp1);
    mpz_init(ztemp2);
    mpz_init(ztemp3);

    init = true;

    mpfr_set_ui(mp_one_over_12, 1, round_mode);                         // mp_one_over_12 = 1/12
    mpfr_div_ui(mp_one_over_12, mp_one_over_12, 12, round_mode);        //

    mpfr_set_ui(mp_one_over_24, 1, round_mode);                         // mp_one_over_24 = 1/24
    mpfr_div_ui(mp_one_over_24, mp_one_over_24, 24, round_mode);        //

    mpfr_t n_minus;                                                     //
    mpfr_init2(n_minus, p);                                             //
    mpfr_set_ui(n_minus, n, round_mode);                                // n_minus = n
    mpfr_sub(n_minus, n_minus, mp_one_over_24,round_mode);              // n_minus = n - 1/24

    mpfr_t sqrt_n_minus;                                                //
    mpfr_init2(sqrt_n_minus, p);                                        //
    mpfr_sqrt(sqrt_n_minus, n_minus, round_mode);                       // n_minus = sqrt(n-1/24)


    mpfr_sqrt_ui(mp_sqrt2, 2, round_mode);                              // mp_sqrt2 = sqrt(2)
    mpfr_sqrt_ui(mp_sqrt3, 3, round_mode);                              // mp_sqrt3 = sqrt(3)
    mpfr_const_pi(mp_pi, round_mode);                                   // mp_pi = pi

    //mp_A = sqrt(2) * 3.1415926535897931 * sqrt(n - 1.0/24.0);-----------
                                                                        //
    mpfr_set(mp_A, mp_sqrt2, round_mode);                               // mp_A = sqrt(2)
    mpfr_mul(mp_A, mp_A, mp_pi, round_mode);                            // mp_A = sqrt(2) * pi
    mpfr_mul(mp_A, mp_A, sqrt_n_minus, round_mode);                     // mp_A = sqrt(2) * pi * sqrt(n - 1/24)
    //--------------------------------------------------------------------

    //cout << "n_minus_1/24: ";
    //mpfr_out_str(stdout, 10, 20, n_minus, round_mode);
    //cout << endl;

    //mp_B = 2.0 * sqrt(3) * (n - 1.0/24.0);------------------------------
    mpfr_set_ui(mp_B, 2, round_mode);                                   // mp_A = 2
    mpfr_mul(mp_B, mp_B, mp_sqrt3, round_mode);                         // mp_A = 2*sqrt(3)
    mpfr_mul(mp_B, mp_B, n_minus, round_mode);                          // mp_A = 2*sqrt(3)*(n-1/24)
    //--------------------------------------------------------------------

    //mp_C = sqrt(2) * pi * sqrt(n - 1.0/24.0) / sqrt(3);-----------------
    mpfr_set(mp_C, mp_sqrt2, round_mode);                               // mp_C = sqrt(2)
    mpfr_mul(mp_C, mp_C, mp_pi, round_mode);                                  // mp_C = sqrt(2) * pi
    mpfr_mul(mp_C, mp_C, sqrt_n_minus, round_mode);                           // mp_C = sqrt(2) * pi * sqrt(n - 1/24)
    mpfr_div(mp_C, mp_C, mp_sqrt3, round_mode);                               // mp_C = sqrt(2) * pi * sqrt(n - 1/24) / sqrt3
    //--------------------------------------------------------------------

    //mp_D = 2.0 * (n - 1.0/24.0) * sqrt(n - 1.0/24.0);-------------------
    mpfr_set_ui(mp_D, 2, round_mode);                                   // mp_D = 2
    mpfr_mul(mp_D, mp_D, n_minus, round_mode);                          // mp_D = 2 * (n - 1/24)
    mpfr_mul(mp_D, mp_D, sqrt_n_minus, round_mode);                     // mp_D = 2 * (n - 1/24) * sqrt(n - 1/24)
    //--------------------------------------------------------------------

    mpfr_clear(n_minus);
    mpfr_clear(sqrt_n_minus);


    // some double precision versions of the above values
    df_A = sqrt(2) * 3.1415926535897931 * sqrt(n - 1.0/24.0);
    df_B = 2.0 * sqrt(3) * (n - 1.0/24.0);
    df_C = sqrt(2) * 3.1415926535897931 * sqrt(n - 1.0/24.0) / sqrt(3);
    df_D = 2.0 * (n - 1.0/24.0) * sqrt(n - 1.0/24.0);

}

void mp_f(mpfr_t result, unsigned int k) {
    //result =  pi * sqrt(2) * cosh(A/(sqrt(3)*k))/(B*k) - sinh(C/k)/D;
    mpfr_set(result, mp_pi, round_mode);                    // result = pi
    mpfr_mul(result, result, mp_sqrt2, round_mode);         // result = sqrt(2) * pi

    mpfr_div(tempf1, mp_A, mp_sqrt3, round_mode);            // temp1 = mp_A/sqrt(3)
    mpfr_div_ui(tempf1, tempf1, k, round_mode);               // temp1 = mp_A/(sqrt(3) * k)
    mpfr_cosh(tempf1, tempf1, round_mode);                    // temp1 = cosh(mp_A/(sqrt(3) * k))
    mpfr_div(tempf1, tempf1, mp_B, round_mode);               // temp1 = cosh(mp_A/(sqrt(3) * k))/mp_B
    mpfr_div_ui(tempf1, tempf1, k, round_mode);               // temp1 = cosh(mp_A/(sqrt(3) * k))/(mp_B*k)

    mpfr_mul(result, result, tempf1, round_mode);            // result = sqrt(2) * pi * cosh(mp_A/(sqrt(3) * k))/(mp_B*k)

    mpfr_div_ui(tempf1, mp_C, k, round_mode);                // temp1 = mp_C/k
    mpfr_sinh(tempf1, tempf1, round_mode);                    // temp1 = sinh(mp_C/k)
    mpfr_div(tempf1, tempf1, mp_D, round_mode);               // temp1 = sinh(mp_C/k)/D

    mpfr_sub(result, result, tempf1, round_mode);            // result = RESULT!
    //return 3.1415926535897931 * sqrt(2) * cosh(df_A/(sqrt(3)*k))/(df_B*k) - sinh(df_C/k)/df_D;
}

void mp_s(mpfr_t result, unsigned int h,unsigned int q) {
    if(q < 3) {
        mpfr_set_ui(result, 0, round_mode);
        return;
    }

    //double result, R1, R2, temp1, temp2;
    unsigned int n, r1, r2, temp3 = 0;

    if(h == 1) {
        if(q < (unsigned int)(sqrt(UINT_MAX))) {
            //elternate, assuming no overflow
            mpfr_set_ui(result, (q-1)*(q-2), round_mode);
            mpfr_div_ui(result, result, 12*q, round_mode);
        }
        else {
            mpz_set_ui(ztemp1, q-1);                                // temp = q-1
            mpz_mul_ui(ztemp1, ztemp1, q-2);                        // temp = (q-1)(q-2)

            mpfr_set_z(result, ztemp1, round_mode);                 // result = (q-1)(q-2)
            mpz_set_ui(ztemp1, q);                                  // temp = q
            mpz_mul_ui(ztemp1, ztemp1, 12);                         // temp = 12q
            mpfr_div_z(result, result, ztemp1, round_mode);         // result = (q-1)(q-2)/12q
        }
/*
        //result = (q-1)*(q-2)/(12*q);
        mpfr_set_ui(result, q-1, round_mode);
        mpfr_mul_ui(result, result, q-2, round_mode);
        mpfr_div_ui(result, result, q, round_mode);         // in this step we don't want to assume that 12*q will not overflow
        mpfr_div_ui(result, result, 12, round_mode);
*/


        return;
    }

    mpfr_set_ui(result, 0, round_mode);             // result = 0

    if(debugs)
        if(h > q) {
            cout << "oops in mp_s()" << endl;
            exit(1);
        }
    r1 = q;
    r2 = h;

    n = 0;
    while(r1 > 0 && r2 > 0) {
        if(r1 < (unsigned int)(sqrt(UINT_MAX)/2.0)) {       // if r1 is small enough we can use
                                                        // standard C integers
                                                        // NOTE: squareroot computation should be optimized by the compiler.
            mpfr_set_ui(temps1, r1*r1 + r2*r2 + 1, round_mode);
            mpfr_div_ui(temps1, temps1, r1 * r2, round_mode);

        }
        else {
            //temp1 = (R1*R1 + R2*R2 + 1.0)/(R1 * R2);      // again, we are NOT going to assume no overflow

            mpfr_set_ui(temps1, r1, round_mode);             // temp1 = r1
            mpfr_mul_ui(temps1, temps1, r1, round_mode);      // temp1 = r1 * r1

            mpfr_set_ui(temps2, r2, round_mode);             // temp2 = r2
            mpfr_mul_ui(temps2, temps2, r2, round_mode);      // temp2 = r2 * r2

            mpfr_add(temps1, temps1, temps2, round_mode);      // temp1 = r1*r1 + r2*r2
            mpfr_add_ui(temps1, temps1, 1, round_mode);       // temp1 = r1*r1 + r2*r2 + 1

            mpfr_div_ui(temps1, temps1, r1, round_mode);      // temp1 = (r1*r1 + r2*r2 + 1)/r1
            mpfr_div_ui(temps1, temps1, r2, round_mode);      // temp1 = (r1*r1 + r2*r2 + 1)/(r1 * r2)
        }
        if(n % 2 == 0){
            mpfr_add(result, result, temps1, round_mode); // result += temp1;
        }
        else {
            mpfr_sub(result, result, temps1, round_mode); // result -= temp1;
        }
        temp3 = r1 % r2;
        r1 = r2;
        r2 = temp3;
        n++;
    }

    //cout << temps1 << endl;
    //cout << result << endl;

    mpfr_mul(result, result, mp_one_over_12, round_mode); // result = result * 1.0/12.0;


    if(n % 2 == 1) {
        mpfr_set_d(temps1, .25, round_mode);
        mpfr_sub(result, result, temps1, round_mode);      // result = result - .25;
    }



    //return result;
}


void mp_a(mpfr_t result, unsigned int n, unsigned int k) {

    if (k == 1) {
        mpfr_set_ui(result, 1, round_mode);     //result = 1
        return;
    }

    mpfr_set_ui(result, 0, round_mode);

    unsigned int h = 0;
    for(h = 1; h < k+1; h++) {
        if(GCD(h,k) == 1) {
            //result += cos(pi * ( s(h,k) - (2.0 * h * n)/k) );
            mp_s(tempa1, h, k);                              // temp1 = s(h,k)

            if(debugs) {
                cout << "s(" << h << "," << k << ") = ";
                mpfr_out_str(stdout, 10, 0, tempa1, round_mode);
                cout << endl;
            }

            mpfr_set_ui(tempa2, 2, round_mode);              // temp2 = 2
            mpfr_mul_ui(tempa2, tempa2, h, round_mode);       // temp2 = 2h
            mpfr_mul_ui(tempa2, tempa2, n, round_mode);       // temp2 = 2hn
            mpfr_div_ui(tempa2, tempa2, k, round_mode);       // temp2 = 2hn/k

            mpfr_sub(tempa1, tempa1, tempa2, round_mode);      // temp1 = s(h,k) - 2hn/k
            mpfr_mul(tempa1, tempa1, mp_pi, round_mode);      // temp1 = pi * (s(h,k) - 2hn/k)
            mpfr_cos(tempa1, tempa1, round_mode);             // temp1 = cos( ditto )

            mpfr_add(result, result, tempa1, round_mode);    // result = result + temp1
        }
    }

}


void mp_t(mpfr_t result, unsigned int n) {
    // NOTE: result should NOT have been initialized when this is called,
    // as we initialize it to the proper precision in this function.

    unsigned int initial_precision = compute_initial_precision(n);
    mpfr_t t1, t2;
    mpfr_init2(t1, initial_precision);
    mpfr_init2(t2, initial_precision);
    mpfr_init2(result, initial_precision);
    mpfr_set_ui(result, 0, round_mode);

    mp_init(initial_precision, n);
    mp_set_precision(initial_precision);

    unsigned int current_precision = initial_precision;
    unsigned int new_precision;

    double remainder = 100.0;
//    d_f_precompute(n);
    unsigned int k = 1;
    for(k = 1; current_precision > min_precision; k++) {
        mpfr_sqrt_ui(t1, k, round_mode);                            // t1 = sqrt(k)

        mp_a(t2, n, k);                                             // t2 = A_k(n)

        if(debuga) {
            cout << "a(" << k <<  ") = ";
            mpfr_out_str(stdout, 10, 0, t2, round_mode);
            cout << endl;
        }

        mpfr_mul(t1, t1, t2, round_mode);                           // t1 = sqrt(k)*A_k(n)

        mp_f(t2, k);                                                // t2 = f_k(n)

        if(debugf) {
            cout << "f(" << k <<  ") = ";
            mpfr_out_str(stdout, 10, 0, t2, round_mode);
            cout << endl;
        }

        mpfr_mul(t1, t1, t2, round_mode);                           // t1 = sqrt(k)*A_k(n)*f_k(n)

        mpfr_add(result, result, t1, round_mode);                   // result += summand

        if(debugt) {
            cout << "Partial sum " << k << " = ";
            mpfr_out_str(stdout, 10, 0, result, round_mode);
            cout << endl;
        }
        //K++;
        //temp = sqrt(K) * d_A(n,k) * d_f(k);
        //result += temp;

        new_precision = compute_current_precision(n,k);
        if(new_precision != current_precision) {
            current_precision = new_precision;
            mp_set_precision(current_precision);
            mpfr_clear(t1); mpfr_clear(t2);
            mpfr_init2(t1, current_precision);
            mpfr_init2(t2, current_precision);
        }
    }

    double result2 = 0;

    for( ; remainder > .5; k++) {
        result2 += sqrt(k) * d_A(n,k) * d_f(k);
        remainder = compute_remainder(n,k);
    }

    mpfr_set_d(t1, result2, round_mode);
    mpfr_add(result, result, t1, round_mode);

//    result = result/(3.1415926535897931 * sqrt(2));
    mpfr_div(result, result, mp_pi, round_mode);
    mpfr_div(result, result, mp_sqrt2, round_mode);
}

//  double versions of the functions

void d_f_precompute(unsigned int n) {
    df_A = sqrt(2) * 3.1415926535897931 * sqrt(n - 1.0/24.0);
    df_B = 2.0 * sqrt(3) * (n - 1.0/24.0);
    df_C = sqrt(2) * 3.1415926535897931 * sqrt(n - 1.0/24.0) / sqrt(3);
    df_D = 2.0 * (n - 1.0/24.0) * sqrt(n - 1.0/24.0);
}

double d_f(unsigned int k) {
    return 3.1415926535897931 * sqrt(2) * cosh(df_A/(sqrt(3)*k))/(df_B*k) - sinh(df_C/k)/df_D;
}

double d_s(unsigned int h,unsigned int q) {
    if(q < 3) {
        return 0.0;
    }

    double result, R1, R2, temp1, temp2;
    unsigned int n, r1, r2, temp3 = 0;

    if(h == 1) {
        double Q;
        Q = q;
        result = (Q-1)*(Q-2)/(12*Q);
        return result;
    }

    result = 0;
    R1 = q;
    R2 = h;

    r1 = q;
    r2 = h;

    n = 0;
    while(r1 && r2) {
        temp1 = (R1*R1 + R2*R2 + 1.0)/(R1 * R2);
        if(n % 2 == 0){
            result += temp1;
        }
        else {
            result -= temp1;
        }
        temp3 = r1 % r2;
        r1 = r2;
        r2 = temp3;
        R1 = r1;
        R2 = r2;
        n++;
    }

    result = result * 1.0/12.0;

    if(n % 2 == 1) {
        result = result - .25;
    }

    return result;
}

double d_A(unsigned int n, unsigned int k) {
    double result;
    result = 0;

    if (k == 1) {
        return 1.0;
    }

    unsigned int h = 0;
    for(h = 1; h < k+1; h++) {
        if(GCD(h,k) == 1) {
        //    cout << "s(" << h << "," << k << ") = " << d_s(h,k) << endl;
        //    cout << "s(" << h << "," << k << ") = " << d_s(h,k) << endl;
            result += cos( 3.1415926535897931 * ( d_s(h,k) - (2.0 * h * n)/double(k)) );
        }
    }
    return result;
}

long GCD(long a, long b)
{
   long u, v, t, x;

   if (a < 0) {
      a = -a;
   }

   if (b < 0) {
      b = -b;
   }


   if (b==0)
      x = a;
   else {
      u = a;
      v = b;
      do {
         t = u % v;
         u = v;
         v = t;
      } while (v != 0);

      x = u;
   }

   return x;
}


/* answer must have already been mpz_init'd. */
int part(mpz_t answer, unsigned int n){
    mpfr_t result;

    mp_t(result, n);

    mpfr_get_z(answer, result, round_mode);

    return 0;
}
