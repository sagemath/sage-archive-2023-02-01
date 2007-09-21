/*      Author:     Jonathan Bober
 *      Version:    .4
 *
 *      This program computes p(n), the number of integer partitions of n, using Rademacher's
 *      formula. (See Hans Rademacher, On the Partition Function p(n),
 *      Proceedings of the London Mathematical Society 1938 s2-43(4):241-254; doi:10.1112/plms/s2-43.4.241,
 *      currently at
 *
 *      http://plms.oxfordjournals.org/cgi/content/citation/s2-43/4/241
 *
 *      if you have access.)
 *
 *      We use the following notation:
 *
 *      p(n) = lim_{n --> oo} t(n,N)
 *
 *      where
 *
 *      t(n,N) = sum_{k=1}^N a(n,k) f_n(k),
 *
 *      where
 *
 *      a(n,k) = sum_{h=1, (h,k) = 1}^k exp(\pi i s(h,k) - 2 \pi i h  n / k)
 *
 *      and
 *
 *      f_n(k) = \pi sqrt{2} cosh(A_n/(sqrt{3}*k))/(B_n*k) - sinh(C_n/k)/D_n;
 *
 *      where
 *
 *      s(h,k) = \sum_{j=1}^{k-1}(j/k)((hj/k))
 *
 *      A_n = sqrt{2} \pi * sqrt{n - 1/24}
 *      B_n = 2 * sqrt{3} * (n - 1/24)
 *      C_n = sqrt{2} * \pi * sqrt{n - 1.0/24.0} / sqrt{3}
 *      D_n = 2 * (n - 1/24) * sqrt{n - 1.0/24.0}
 *
 *      and, finally, where ((x)) is the sawtooth function ((x)) = {x} - 1/2 if x is not an integer, 0 otherwise.
 *
 *      Some clever tricks are used in the computation of s(h,k), and perhaps at least
 *      some of these are due to Apostol. (I don't know a reference for this.)
 *
 *      TODO:
 *
 *          -- Search source code for other TODO comments.
 *
 *      OTHER CREDITS:
 *
 *      I looked source code written by Ralf Stephan, currently available at
 *
 *              http://www.ark.in-berlin.de/part.c
 *
 *      while writing this code, but didn't really make use of it, except for the
 *      function cospi(), currently included in the source (slightly modified), but not
 *      used.
 *
 *      More useful were notes currently available at
 *
 *              http://www.ark.in-berlin.de/part.pdf
 *
 *      and others at
 *
 *              http://www.math.uwaterloo.ca/~dmjackso/CO630/ptnform1.pdf
 *
 *      Also, it is worth noting that the code for GCD(), while trivial,
 *      was directly copied from the NTL source code.
 *
 *      Also, Bill Hart made some comments about ways to speed up this computation on the SAGE
 *      mailing list.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 2 of the License, or
 *      (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      You should have received a copy of the GNU General Public License
 *      along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#if defined(__sun)
extern "C" long double fabsl (long double);
extern "C" long double sinl (long double);
extern "C" long double cosl (long double);
extern "C" long double sqrtl (long double);
extern "C" long double coshl (long double);
extern "C" long double sinhl (long double);
#endif




#include <stdio.h>
#include <cfloat>

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

const unsigned int min_precision = DBL_MANT_DIG;
const unsigned int double_precision = DBL_MANT_DIG;
const unsigned int long_double_precision = LDBL_MANT_DIG;

//const unsigned int level_two_precision = long_double_precision; // qd_real precision
const unsigned int level_two_precision = 200; // qd_real precision
const unsigned int level_three_precision = 128;
const unsigned int level_four_precision = long_double_precision;
//const unsigned int level_two_precision = double_precision;
const unsigned int level_five_precision = double_precision;

const long double ld_pi = 3.141592653589793238462643L;
const double d_pi = ld_pi;

bool test(bool longtest = false);

void cospi (mpfr_t res, mpfr_t x);




// The following function can be useful for debugging in come circumstances, but should not be used for anything else
// unless it is rewritten.
int grab_last_digits(char * output, int n, mpfr_t x) {
    // fill output with the n digits of x that occur
    // just before the decimal point
    // Note: this assumes that x has enough digits and enough
    // precision -- otherwise bad things can happen
    //
    // returns: the number of digits to the right of the decimal point

    char * temp;
    mp_exp_t e;

    temp = mpfr_get_str(NULL, &e, 10, 0, x, GMP_RNDN);

    int retval;

    if(e > 0) {
        strncpy(output, temp + e - n, n);
        retval =  strlen(temp + e);
    }
    else {
        for(int i = 0; i < n; i++)
            output[i] = '0';
        retval = strlen(temp);
    }
    output[n] = '\0';

    mpfr_free_str(temp);

    return retval;
}


long GCD(long a, long b);

unsigned int calc_precision(unsigned int n, unsigned int N);

//mpfr versions

mpz_t ztemp1, ztemp2, ztemp3;
mpq_t qtemp1, qtemp2, qtemp3;

mpfr_t mp_one_over_12, mp_one_over_24, mp_sqrt2, mp_sqrt3, mp_pi, half, fourth;
mpfr_t mp_A, mp_B, mp_C, mp_D;


mp_rnd_t round_mode = GMP_RNDN;

mpfr_t tempa1, tempa2, tempf1, tempf2, temps1, temps2, tempt1, tempt2;  // temp variables for different functions, with precision set and cleared by initialize_mpfr_variables
mpfr_t tempc1, tempc2; // temp variable used by cospi()


bool mp_vars_initialized = false;

mp_prec_t mp_precision;

void initialize_constants(unsigned int prec, unsigned int n);
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

double d_A, d_B, d_C, d_D;
long double ld_A, ld_B, ld_C, ld_D;


double d_f(unsigned int k);
double d_s(unsigned int h,unsigned int q);
double d_a(unsigned int n, unsigned int k);
double d_t(unsigned int n, unsigned int N);

long double ld_f(unsigned int k);
long double ld_s(unsigned int h,unsigned int q);
long double ld_a(unsigned int n, unsigned int k);
long double ld_t(unsigned int n, unsigned int N);

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

    unsigned int result = (unsigned int)(ceil(3.1415926535897931 * sqrt(2.0 * double(n)/ 3.0) / log(2))) + 3;
    if(debug) cout << "Using initial precision of " << result << " bits." << endl;

    if(result > min_precision) {
        return result;
    }

    else return min_precision;

}

unsigned int compute_current_precision(unsigned int n, unsigned int N) {
    // we want to compute
    //      log(A/sqrt(N) + B*sqrt(N/(n-1))*sinh(C * sqrt(n) / N) / log(2)
    //
    //  where A, B, and C are the constants listed below. These error bounds
    //  are given in the paper by Rademacher listed at the top of this file.
    //

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

    unsigned int p = mpfr_get_exp(error);   // I am almost certain that this does the right thing
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
        return p;                           // don't want to return < min_precision
                                            // Note that when we hit the minimum precision
                                            // we should switch over to using C doubles instead
                                            // of mpfr types.
    }
    return min_precision;
}

double compute_remainder(unsigned int n, unsigned int N) {
    // This computes the remainer left after N terms have been computed.
    // The formula is exactly the same as the one used to compute the required
    // precision, but once we know the necessary precision is small, we can
    // call this function to determine the actual error (rather than the precision).
    //
    // Generally, this is only called once we know that the necessary
    // precision is <= min_precision, because then the error is small
    // enough to fit into a double, and also, we know that we are
    // getting close to the correct answer.

    double A = 1.11431833485164;
    double B = 0.059238439175445;
    double C = 2.5650996603238;
    double result;
    result = A/sqrt(N) + B * sqrt(double(N)/double(n-1))*sinh(C * sqrt(double(n))/double(N));
    return result;
}


void initialize_mpfr_variables(unsigned int prec) {
    //
    // Clear and initialize some "temp" variables that are used in the computation of various functions.
    //
    // We do this in this auxilliary function (and endure the pain of
    // the extra global variables) so that we only have to initialize/clear
    // these variables once every time the precision changes.
    //
    // NAMING CONVENTIONS:
    //
    // -tempa1 and tempa2 are two variables available for use
    //      in the function mp_a()
    //
    // -tempf1 and tempf2 are two variables available for use
    //      in the function mp_f()
    //
    // -etc...
    //
    // NOTE: Calls to this function must be paired with calls to clear_mpfr_variables()

    mpfr_init2(tempa1, prec);
    mpfr_init2(tempa2, prec);
    mpfr_init2(tempf1, prec);
    mpfr_init2(tempf2, prec);
    mpfr_init2(tempt1, prec);
    mpfr_init2(tempt2, prec);
    mpfr_init2(temps1, prec);
    mpfr_init2(temps2, prec);
    mpfr_init2(tempc1, prec);
    mpfr_init2(tempc2, prec);
}

void clear_mpfr_variables() {
    mpfr_clear(tempa1);
    mpfr_clear(tempa2);
    mpfr_clear(tempf1);
    mpfr_clear(tempf2);
    mpfr_clear(tempt1);
    mpfr_clear(tempt2);
    mpfr_clear(temps1);
    mpfr_clear(temps2);

    mpfr_clear(tempc1);
    mpfr_clear(tempc2);
}


void initialize_constants(unsigned int prec, unsigned int n) {
    // The variables mp_A, mp_B, mp_C, and mp_D are used for
    // A_n, B_n, C_n, and D_n listed at the top of this file.
    //
    // They depend only on n, so we compute them just once in this function,
    // and then use them many times elsewhere.
    //
    // Also, we precompute some extra constants that we use a lot, such as
    // sqrt2, sqrt3, pi, 1/24, 1/12, etc.
    //
    // NOTE: Calls to this function must be paired with calls to clear_constants()
    static bool init = false;
    mp_precision = prec;
    mp_prec_t p = mp_precision;

    mpfr_init2(mp_one_over_12,p); mpfr_init2(mp_one_over_24,p); mpfr_init2(mp_sqrt2,p); mpfr_init2(mp_sqrt3,p); mpfr_init2(mp_pi,p);
    mpfr_init2(mp_A,p); mpfr_init2(mp_B,p); mpfr_init2(mp_C,p); mpfr_init2(mp_D,p); mpfr_init2(fourth, p); mpfr_init2(half, p);

    init = true;

    mpfr_set_ui(mp_one_over_12, 1, round_mode);                             // mp_one_over_12 = 1/12
    mpfr_div_ui(mp_one_over_12, mp_one_over_12, 12, round_mode);            //


    mpfr_set_ui(mp_one_over_24, 1, round_mode);                             // mp_one_over_24 = 1/24
    mpfr_div_ui(mp_one_over_24, mp_one_over_24, 24, round_mode);            //


    mpfr_set_ui(half, 1, round_mode);                                       //
    mpfr_div_ui(half, half, 2, round_mode);                                    // half = 1/2
    mpfr_div_ui(fourth, half, 2, round_mode);                                  // fourth = 1/4


    mpfr_t n_minus;                                                         //
    mpfr_init2(n_minus, p);                                                 //
    mpfr_set_ui(n_minus, n, round_mode);                                    // n_minus = n
    mpfr_sub(n_minus, n_minus, mp_one_over_24,round_mode);                  // n_minus = n - 1/24

    mpfr_t sqrt_n_minus;                                                    //
    mpfr_init2(sqrt_n_minus, p);                                            //
    mpfr_sqrt(sqrt_n_minus, n_minus, round_mode);                           // n_minus = sqrt(n-1/24)


    mpfr_sqrt_ui(mp_sqrt2, 2, round_mode);                                  // mp_sqrt2 = sqrt(2)
    mpfr_sqrt_ui(mp_sqrt3, 3, round_mode);                                  // mp_sqrt3 = sqrt(3)
    mpfr_const_pi(mp_pi, round_mode);                                       // mp_pi = pi


    //mp_A = sqrt(2) * 3.1415926535897931 * sqrt(n - 1.0/24.0);---------------
                                                                            //
    mpfr_set(mp_A, mp_sqrt2, round_mode);                                   // mp_A = sqrt(2)
    mpfr_mul(mp_A, mp_A, mp_pi, round_mode);                                // mp_A = sqrt(2) * pi
    mpfr_mul(mp_A, mp_A, sqrt_n_minus, round_mode);                         // mp_A = sqrt(2) * pi * sqrt(n - 1/24)
    //------------------------------------------------------------------------

    //cout << "mp_A = ";
    //mpfr_out_str(stdout, 10, 0, mp_A, round_mode);
    //cout << endl;

    //mp_B = 2.0 * sqrt(3) * (n - 1.0/24.0);----------------------------------
    mpfr_set_ui(mp_B, 2, round_mode);                                       // mp_A = 2
    mpfr_mul(mp_B, mp_B, mp_sqrt3, round_mode);                             // mp_A = 2*sqrt(3)
    mpfr_mul(mp_B, mp_B, n_minus, round_mode);                              // mp_A = 2*sqrt(3)*(n-1/24)
    //------------------------------------------------------------------------

    //cout << "mp_B = ";
    //mpfr_out_str(stdout, 10, 0, mp_B, round_mode);
    //cout << endl;

    //mp_C = sqrt(2) * pi * sqrt(n - 1.0/24.0) / sqrt(3);---------------------
    mpfr_set(mp_C, mp_sqrt2, round_mode);                                   // mp_C = sqrt(2)
    mpfr_mul(mp_C, mp_C, mp_pi, round_mode);                                // mp_C = sqrt(2) * pi
    mpfr_mul(mp_C, mp_C, sqrt_n_minus, round_mode);                         // mp_C = sqrt(2) * pi * sqrt(n - 1/24)
    mpfr_div(mp_C, mp_C, mp_sqrt3, round_mode);                             // mp_C = sqrt(2) * pi * sqrt(n - 1/24) / sqrt3
    //------------------------------------------------------------------------

    //cout << "mp_C = ";
    //mpfr_out_str(stdout, 10, 0, mp_C, round_mode);
    //cout << endl;

    //mp_D = 2.0 * (n - 1.0/24.0) * sqrt(n - 1.0/24.0);-----------------------
    mpfr_set_ui(mp_D, 2, round_mode);                                       // mp_D = 2
    mpfr_mul(mp_D, mp_D, n_minus, round_mode);                              // mp_D = 2 * (n - 1/24)
    mpfr_mul(mp_D, mp_D, sqrt_n_minus, round_mode);                         // mp_D = 2 * (n - 1/24) * sqrt(n - 1/24)
    //------------------------------------------------------------------------

    //cout << "mp_D = ";
    //mpfr_out_str(stdout, 10, 0, mp_D, round_mode);
    //cout << endl;

    mpfr_clear(n_minus);
    mpfr_clear(sqrt_n_minus);


    // some double precision versions of the above values

    d_A = sqrt(2) * d_pi * sqrt(n - 1.0/24.0);
    d_B = 2.0 * sqrt(3) * (n - 1.0/24.0);
    d_C = sqrt(2) * d_pi * sqrt(n - 1.0/24.0) / sqrt(3);
    d_D = 2.0 * (n - 1.0/24.0) * sqrt(n - 1.0/24.0);

    ld_A = sqrtl(2) * d_pi * sqrtl(n - 1.0L/24.0L);
    ld_B = 2.0L * sqrtl(3) * (n - 1.0/24.0);
    ld_C = sqrt(2) * d_pi * sqrtl(n - 1.0L/24.0L) / sqrtl(3);
    ld_D = 2.0L * (n - 1.0L/24.0L) * sqrtl(n - 1.0L/24.0L);

}

void clear_constants() {
        mpfr_clear(mp_one_over_12); mpfr_clear(mp_one_over_24); mpfr_clear(mp_sqrt2); mpfr_clear(mp_sqrt3); mpfr_clear(mp_pi);
        mpfr_clear(mp_A); mpfr_clear(mp_B); mpfr_clear(mp_C); mpfr_clear(mp_D); mpfr_clear(half); mpfr_clear(fourth);
}

void initialize_mpz_and_mpq_variables() {
    /*
     * We use a few mpz_t and mpq_t variables which need to be initialized
     * before they can be used. Initialization and clearing take some
     * time, so we initialize just once in this function, and clear in another.
     */
    mpz_init(ztemp1);
    mpz_init(ztemp2);
    mpz_init(ztemp3);

    mpq_init(qtemp1);
    mpq_init(qtemp2);
    mpq_init(qtemp3);
}

void clear_mpz_and_mpq_variables() {
    mpz_clear(ztemp1);
    mpz_clear(ztemp2);
    mpz_clear(ztemp3);

    mpq_clear(qtemp1);
    mpq_clear(qtemp2);
    mpq_clear(qtemp3);
}

void mp_f(mpfr_t result, unsigned int k) {
    // compute f_n(k) as described in the introduction
    //
    // notice that this doesn't use n - the "constants"
    // A, B, C, and D depend on n, but they are precomputed
    // once before this function is called.

    //result =  pi * sqrt(2) * cosh(A/(sqrt(3)*k))/(B*k) - sinh(C/k)/D;

                                                                    //
    mpfr_set(result, mp_pi, round_mode);                            // result = pi

    mpfr_mul(result, result, mp_sqrt2, round_mode);                 // result = sqrt(2) * pi
                                                                    //
    mpfr_div(tempf1, mp_A, mp_sqrt3, round_mode);                   // temp1 = mp_A/sqrt(3)
    mpfr_div_ui(tempf1, tempf1, k, round_mode);                     // temp1 = mp_A/(sqrt(3) * k)
    mpfr_cosh(tempf1, tempf1, round_mode);                          // temp1 = cosh(mp_A/(sqrt(3) * k))
    mpfr_div(tempf1, tempf1, mp_B, round_mode);                     // temp1 = cosh(mp_A/(sqrt(3) * k))/mp_B
    mpfr_div_ui(tempf1, tempf1, k, round_mode);                     // temp1 = cosh(mp_A/(sqrt(3) * k))/(mp_B*k)
                                                                    //
    mpfr_mul(result, result, tempf1, round_mode);                   // result = sqrt(2) * pi * cosh(mp_A/(sqrt(3) * k))/(mp_B*k)
                                                                    //
    mpfr_div_ui(tempf1, mp_C, k, round_mode);                       // temp1 = mp_C/k
    mpfr_sinh(tempf1, tempf1, round_mode);                          // temp1 = sinh(mp_C/k)
    mpfr_div(tempf1, tempf1, mp_D, round_mode);                     // temp1 = sinh(mp_C/k)/D
                                                                    //
    mpfr_sub(result, result, tempf1, round_mode);                   // result = RESULT!
                                                                    //
    //return pi * sqrt(2) * cosh(A/(sqrt(3)*k))/(B*k) - sinh(C/k)/D;
}

// call the following when 4k < sqrt(UINT_MAX)
//
// TODO: maybe a faster version of this can be written without using mpq_t,
//       or maybe this version can be smarter.
//
//       The actual value of the cosine that we compute using s(h,k)
//       only depends on {s(h,k)/2}, that is, the fractional
//       part of s(h,k)/2. It may be possible to make use of this somehow.
//
//       NOTE: when we compute p(1000000000),
//       it takes about 3m 30s (on my laptop). If we uncomment
//       the first two lines below, so that this function doesn't actually
//       compute anything, it takes about 3m 10s to run. So not that much time
//       is spent in this function when computing for large n
void q_s(mpq_t result, unsigned int h, unsigned int k) {
    //mpq_set_ui(result, 1, 1);
    //return;

    if(k < 3) {
        mpq_set_ui(result, 0, 1);
        return;
    }

    if (h == 1) {
        unsigned int d = GCD( (k-1)*(k-2), 12*k);
        if(d > 1) {
            mpq_set_ui(result, ((k-1)*(k-2))/d, (12*k)/d);
        }
        else {
            mpq_set_ui(result, (k-1)*(k-2), 12*k);
        }
        //mpq_canonicalize(result);
        return;
    }
    // TODO: In the function mp_s() there are a few special cases for special forms of h and k.
    // (And there are more special cases listed in one of the references listed in the introduction.)
    //
    // It may be advantageous to implement some here, but I'm not sure
    // if there is any real speed benefit to this.
    //
    // In the mpfr_t version of this function, the speedups didn't seem to help too much, but
    // they might make more of a difference when using mpq_t.

    // if h = 2 and k is odd, we have
    // s(h,k) = (k-1)*(k-5)/24k
    //if(h == 2 && k > 5 && k % 2 == 1) {
    //    unsigned int d = GCD( (k-1)*(k-5), 24*k);
    //    if(d > 1) {
    //        mpq_set_ui(result, ((k-1)*(k-5))/d, (24*k)/d);
    //    }
    //    else {
    //        mpq_set_ui(result, (k-1)*(k-5), 24*k);
    //    }
    //    return;
    //}

/*

    // if k % h == 1, then
    //
    //      s(h,k) = (k-1)(k - h^2 - 1)/(12hk)
    //

    // We need to be a little careful here because k - h^2 - 1 can be negative.
    if(k % h == 1) {
        int num = (k-1)*(k - h*h - 1);
        int den = 12*k*h;
        int d = GCD(num, den);
        if(d > 1) {
            mpq_set_si(result, num/d, den/d);
        }
        else {
            mpq_set_si(result, num, den);
        }
        return;
    }

    // if k % h == 2, then
    //
    //      s(h,k) = (k-2)[k - .5(h^2 + 1)]/(12hk)
    //

    //if(k % h == 2) {
    //}
*/




    mpq_set_ui(result, 0, 1);                             // result = 0

    int r1 = k;
    int r2 = h;

    int n = 0;
    int temp3;
    while(r1 > 0 && r2 > 0) {
        unsigned int d = GCD(r1 * r1 + r2 * r2 + 1, r1 * r2);
        if(d > 1) {
            mpq_set_ui(qtemp1, (r1 * r1 + r2 * r2 + 1)/d, (r1 * r2)/d);
        }
        else{
            mpq_set_ui(qtemp1, r1 * r1 + r2 * r2 + 1, r1 * r2);
        }
        //mpq_canonicalize(qtemp1);

        if(n % 2 == 0){                                             //
            mpq_add(result, result, qtemp1);                        // result += temp1;
        }                                                           //
        else {                                                      //
            mpq_sub(result, result, qtemp1);                        // result -= temp1;
        }                                                           //
        temp3 = r1 % r2;                                            //
        r1 = r2;                                                    //
        r2 = temp3;                                                 //
        n++;                                                        //
    }

    mpq_set_ui(qtemp1, 1, 12);
    mpq_mul(result, result, qtemp1);                                // result = result * 1.0/12.0;


    if(n % 2 == 1) {
        mpq_set_ui(qtemp1, 1, 4);
        mpq_sub(result, result, qtemp1);                            // result = result - .25;
    }

}


void mp_s(mpfr_t result, unsigned int h,unsigned int k) {
    // Compute s(h,k) using some clever tricks.

    // If k < 3, then we know that the result is going to be 0.
    if(k < 3) {
        mpfr_set_ui(result, 0, round_mode);
        return;
    }

    unsigned int n, r1, r2, temp3 = 0;

    // If h = 1, then the result has the form
    //
    //      s(h,k) = (k-1)(k-2)/(12k).
    //

    if(h == 1) {

        // When k is small enough, we can do some of the computation using
        // just ordinary integer arithmetic.

        if(k < (unsigned int)(sqrt(UINT_MAX))) {
            mpfr_set_ui(result, (k-1)*(k-2), round_mode);
            mpfr_div_ui(result, result, 12*k, round_mode);
        }
        else {
            // When k is very large, we might need to use this code.
            mpz_set_ui(ztemp1, k-1);                                // temp = k-1
            mpz_mul_ui(ztemp1, ztemp1, k-2);                        // temp = (k-1)(k-2)

            mpfr_set_z(result, ztemp1, round_mode);                 // result = (k-1)(k-2)
            mpz_set_ui(ztemp1, k);                                  // temp = k
            mpz_mul_ui(ztemp1, ztemp1, 12);                         // temp = 12k
            mpfr_div_z(result, result, ztemp1, round_mode);         // result = (k-1)(k-2)/12k
        }
        return;
    }

    // TODO: Below are a few special cases for special forms of h and k.
    //
    // It may be advantageous to implement a few more, but I'm not even sure
    // that there is any real speed benefit to the special forms that are here
    // now.

    // if h = 2 and k is odd, then s(h,k) is given by
    // (we need k > 5 because k is unsigned. k = 3 or 5 should
    // be special cased.)
    //
    //      s(h,k) = (k-1)(k-5)/(24k)
    //
    if(h == 2 && k > 5 && k % 2 == 1) {

        // When k is small enough, we can do some of the computation using
        // just ordinary integer arithmetic.

        if(k < (unsigned int)(sqrt(UINT_MAX))) {
            mpfr_set_ui(result, (k-1)*(k-5), round_mode);
            mpfr_div_ui(result, result, 24*k, round_mode);
        }
        else {
            // When k is very large, we might need to use this code.
            mpz_set_ui(ztemp1, k-1);                                // temp = k-1
            mpz_mul_ui(ztemp1, ztemp1, k-5);                        // temp = (k-1)(k-5)

            mpfr_set_z(result, ztemp1, round_mode);                 // result = (k-1)(k-5)
            mpz_set_ui(ztemp1, k);                                  // temp = k
            mpz_mul_ui(ztemp1, ztemp1, 24);                         // temp = 24k
            mpfr_div_z(result, result, ztemp1, round_mode);         // result = (k-1)(k-5)/24k
        }
        return;
    }

    // if k % h == 1, then
    //
    //      s(h,k) = (k-1)(k - h^2 - 1)/(12hk)
    //

    if(k % h == 1) {
        if(4*k < (unsigned int)(sqrt(UINT_MAX))) {
            mpfr_set_si(result, (k-1)*(k - h*h - 1), round_mode);
            mpfr_div_ui(result, result, 12*h*k, round_mode);
            return;
        }
    }

    // if k % h == 2, then
    //
    //      s(h,k) = (k-2)[k - .5(h^2 + 1)]/(12hk)
    //
    //
    //



    // At this point we have given up hope of using a special form and fall back on our generic
    // algorithm.


    mpfr_set_ui(result, 0, round_mode);                             // result = 0

    r1 = k;
    r2 = h;

    n = 0;
    while(r1 > 0 && r2 > 0) {
        if(r1 < (unsigned int)(sqrt(UINT_MAX)/2.0)) {               // if r1 is small enough we can use
                                                                    // standard C integers
                                                                    // NOTE: squareroot computation should be optimized by the compiler.
            mpfr_set_ui(temps1, r1*r1 + r2*r2 + 1, round_mode);
            mpfr_div_ui(temps1, temps1, r1 * r2, round_mode);

        }
        else {
            //temp1 = (R1*R1 + R2*R2 + 1.0)/(R1 * R2);              //
                                                                    //
            mpfr_set_ui(temps1, r1, round_mode);                    // temp1 = r1
            mpfr_mul_ui(temps1, temps1, r1, round_mode);            // temp1 = r1 * r1
                                                                    //
            mpfr_set_ui(temps2, r2, round_mode);                    // temp2 = r2
            mpfr_mul_ui(temps2, temps2, r2, round_mode);            // temp2 = r2 * r2
                                                                    //
            mpfr_add(temps1, temps1, temps2, round_mode);           // temp1 = r1*r1 + r2*r2
            mpfr_add_ui(temps1, temps1, 1, round_mode);             // temp1 = r1*r1 + r2*r2 + 1
                                                                    //
            mpfr_div_ui(temps1, temps1, r1, round_mode);            // temp1 = (r1*r1 + r2*r2 + 1)/r1
            mpfr_div_ui(temps1, temps1, r2, round_mode);            // temp1 = (r1*r1 + r2*r2 + 1)/(r1 * r2)
        }                                                           //
        if(n % 2 == 0){                                             //
            mpfr_add(result, result, temps1, round_mode);           // result += temp1;
        }                                                           //
        else {                                                      //
            mpfr_sub(result, result, temps1, round_mode);           // result -= temp1;
        }                                                           //
        temp3 = r1 % r2;                                            //
        r1 = r2;                                                    //
        r2 = temp3;                                                 //
        n++;                                                        //
    }


    mpfr_mul(result, result, mp_one_over_12, round_mode);           // result = result * 1.0/12.0;


    if(n % 2 == 1) {
        mpfr_set_d(temps1, .25, round_mode);
        mpfr_sub(result, result, temps1, round_mode);               // result = result - .25;
    }

}


void mp_a(mpfr_t result, unsigned int n, unsigned int k) {
    // compute a(n,k)

    if (k == 1) {
        mpfr_set_ui(result, 1, round_mode);                         //result = 1
        return;
    }

    mpfr_set_ui(result, 0, round_mode);

    unsigned int h = 0;
    for(h = 1; h < k+1; h++) {
        if(GCD(h,k) == 1) {

            // Note that we compute each term of the summand as
            //      result += cos(pi * ( s(h,k) - (2.0 * h * n)/k) );
            //
            // This is the real part of the exponential that was written
            // down in the introduction, and we don't need to compute
            // the imaginary part because we know that, in the end, the
            // imaginary part will be 0, as we are computing an integer.

            if(4*k < (unsigned int)(sqrt(UINT_MAX))) {
                q_s(qtemp2, h, k);

                //mpfr_mul_q(tempa1, mp_pi, qtemp2, round_mode);
                //mpfr_mul_ui(tempa1, tempa1, k * k, round_mode);

                //mpfr_set_q(tempa1, qtemp2, round_mode);
                unsigned int r = n % k;                                     // here we make use of the fact that the
                unsigned int d = GCD(r,k);                                  // cos() term written above only depends
                unsigned int K;                                             // on {hn/k}.
                if(d > 1) {
                    r = r/d;
                    K = k/d;
                }
                else {
                    K = k;
                }
                if(K % 2 == 0) {
                    K = K/2;
                }
                else {
                    r = r * 2;
                }
                mpq_set_ui(qtemp3, h*r, K);
                mpq_sub(qtemp2, qtemp2, qtemp3);
                /*
                mpfr_set_q(tempa2, qtemp2, round_mode);                 // This might be faster, according to
                cospi(tempa1, tempa2);                                  // the comments in Ralf Stephan's part.c, but
                                                                        // I haven't noticed a significant speed up.
                                                                        // (Perhaps a different version that takes an mpq_t
                                                                        // as an input might be faster.)
                */
                mpfr_mul_q(tempa1, mp_pi, qtemp2, round_mode);
                mpfr_cos(tempa1, tempa1, round_mode);
                mpfr_add(result, result, tempa1, round_mode);
            }
            else{
                mp_s(tempa1, h, k);                                     // temp1 = s(h,k)

                mpfr_set_ui(tempa2, 2, round_mode);                     // temp2 = 2
                mpfr_mul_ui(tempa2, tempa2, h, round_mode);             // temp2 = 2h
                mpfr_mul_ui(tempa2, tempa2, n, round_mode);             // temp2 = 2hn
                mpfr_div_ui(tempa2, tempa2, k, round_mode);             // temp2 = 2hn/k

                mpfr_sub(tempa1, tempa1, tempa2, round_mode);           // temp1 = s(h,k) - 2hn/k
                mpfr_mul(tempa1, tempa1, mp_pi, round_mode);            // temp1 = pi * (s(h,k) - 2hn/k)
                mpfr_cos(tempa1, tempa1, round_mode);                   // temp1 = cos( ditto )

                mpfr_add(result, result, tempa1, round_mode);           // result = result + temp1
            }

        }

    }

}


void mp_t(mpfr_t result, unsigned int n) {
    // This is the function that actually computes p(n).
    //
    // More specifically, it computes t(n,N) to within 1/2 of p(n), and then
    // we can find p(n) by rounding.
    //
    //
    // NOTE: result should NOT have been initialized when this is called,
    // as we initialize it to the proper precision in this function.

    unsigned int initial_precision = compute_initial_precision(n);  // We begin by computing the precision necessary to hold the final answer.
                                                                    // and then initialize both the result and some temporary variables to that
                                                                    // precision.
    mpfr_t t1, t2;                                                  //
    mpfr_init2(t1, initial_precision);                              //
    mpfr_init2(t2, initial_precision);                              //
    mpfr_init2(result, initial_precision);                          //
    mpfr_set_ui(result, 0, round_mode);                             //

    initialize_mpz_and_mpq_variables();
    initialize_constants(initial_precision, n);                     // Now that we have the precision information, we initialize some constants
                                                                    // that will be used throughout, and also
                                                                    //
    initialize_mpfr_variables(initial_precision);                            // set the precision of the "temp" variables that are used in individual functions.
//
    unsigned int current_precision = initial_precision;
    unsigned int new_precision;

    double remainder = 0.5772156649;                                // (We just need the remainder to be initialized to something. This
                                                                    // seems like as good a number as any.)


    // We start by computing with high precision arithmetic, until
    // we are sure enough that we don't need that much precision
    // anymore. Once current_precision == min_precision, we drop
    // out of this loop and switch to a computation
    // that only involves doubles.

    unsigned int k = 1;                                             // (k holds the index of the summand that we are computing.)
    for(k = 1; current_precision > level_four_precision; k++) {            //
        mpfr_sqrt_ui(t1, k, round_mode);                            // t1 = sqrt(k)
                                                                    //
        mp_a(t2, n, k);                                             // t2 = A_k(n)

        if(debuga) {
            cout << "a(" << k <<  ") = ";
            mpfr_out_str(stdout, 10, 10, t2, round_mode);
            cout << endl;
        }

        mpfr_mul(t1, t1, t2, round_mode);                           // t1 = sqrt(k)*A_k(n)
                                                                    //
        mp_f(t2, k);                                                // t2 = f_k(n)

        if(debugf) {
            cout << "f(" << n << "," << k <<  ") = ";
            mpfr_out_str(stdout, 10, 0, t2, round_mode);
            cout << endl;
        }

        mpfr_mul(t1, t1, t2, round_mode);                           // t1 = sqrt(k)*A_k(n)*f_k(n)
                                                                    //
        mpfr_add(result, result, t1, round_mode);                   // result += summand

        if(debugt) {
            //cout << "Partial sum " << k << " = ";
            //mpfr_out_str(stdout, 10, 0, result, round_mode);
            //cout << endl;
            int num_digits = 20;
            int num_extra_digits;
            char digits[num_digits + 1];
            num_extra_digits = grab_last_digits(digits, 5, t1);
            grab_last_digits(digits, num_digits, result);

            mpfr_out_str(stdout, 10, 10, t1, round_mode);
            cout << endl;

            cout << k << ": current precision:"  << current_precision << ". 20 last digits of partial result: " << digits << ". Number of extra digits: " << num_extra_digits << endl;
            cout.flush();

        }

        new_precision = compute_current_precision(n,k);             // After computing one summand, check what the new precision should be.
        if(new_precision != current_precision) {                    // If the precision changes, we need to clear
            current_precision = new_precision;                      // and reinitialize all "temp" variables to
            clear_mpfr_variables();                                 // use lower precision.
            initialize_mpfr_variables(current_precision);           //
            mpfr_clear(t1); mpfr_clear(t2);                         //
            mpfr_init2(t1, current_precision);                      //
            mpfr_init2(t2, current_precision);                      //
        }
    }


    mpfr_clear(t1); mpfr_clear(t2);

    mpfr_init2(t1, 200);
    mpfr_init2(t2, 200);


    long double tail_result_1 = 0;

    for( ; current_precision > level_five_precision; k++) {        // (don't change k -- it is already the right value)
        tail_result_1 += sqrtl(k) * ld_a(n,k) * ld_f(k);            //
        current_precision = compute_current_precision(n,k);         // The only reason that we compute the new precision
                                                                    // now is so that we know when we can change to using just doubles.
                                                                    // (There should be a 'long double' version of the compute_current_precision function.
    }


    double tail_result_2 = 0;                                       // (tail_result_2 will hold the result of the "tail end"
                                                                    // computation using doubles.)

    for( ; remainder > .5; k++) {                                   // (don't change k -- it is already the right value)
        tail_result_2 += sqrt(k) * d_a(n,k) * d_f(k);                     //
        remainder = compute_remainder(n,k);                         // Now we start computing the size of the remainder. Once
                                                                    // it is small enough, we know that we have the answer.
    }

    mpfr_set_d(t1, tail_result_2, round_mode);                      //
    mpfr_add(result, result, t1, round_mode);                       // We add together the main result and the tail ends'

    mpfr_set_ld(t1, tail_result_1, round_mode);
    mpfr_add(result, result, t1, round_mode);

    //cout << tail_result_0_str << endl;
    //mpfr_out_str(stdout, 10, 0, t1, round_mode);
    //cout << endl;

    mpfr_div(result, result, mp_pi, round_mode);                    // The actual result is the sum that we have computed
    mpfr_div(result, result, mp_sqrt2, round_mode);                 // divided by pi*sqrt(2).

    clear_constants();
    clear_mpz_and_mpq_variables();
    clear_mpfr_variables();
    mpfr_clear(t1);
    mpfr_clear(t2);                         //
}

//  Double versions of the functions, see the above functions for documentation.

double d_f(unsigned int k) {
    return  3.141592653589793238462643 * sqrt(2) * cosh(d_A/(sqrt(3)*k))/(d_B*k) - sinh(d_C/k)/d_D;
}

long double ld_f(unsigned int k) {
    return  3.141592653589793238462643 * sqrtl(2) * coshl(ld_A/(sqrtl(3)*k))/(ld_B*k) - sinhl(ld_C/k)/ld_D;
}

double d_s(unsigned int h,unsigned int k) {
    if(k < 3) {
        return 0.0;
    }

    double result, R1, R2, temp1, temp2;
    unsigned int n, r1, r2, temp3 = 0;

    if(h == 1) {
        double K;
        K = k;
        result = (K-1)*(K-2)/(12*K);
        return result;
    }

    result = 0;
    R1 = k;
    R2 = h;

    r1 = k;
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

long double ld_s(unsigned int h,unsigned int k) {
    if(k < 3) {
        return 0.0;
    }

    long double result, R1, R2, temp1, temp2;
    unsigned int n, r1, r2, temp3 = 0;

    if(h == 1) {
        long double K;
        K = k;
        result = (K-1)*(K-2)/(12*K);
        return result;
    }

    result = 0;
    R1 = k;
    R2 = h;

    r1 = k;
    r2 = h;

    n = 0;
    while(r1 && r2) {
        temp1 = (R1*R1 + R2*R2 + 1.0L)/(R1 * R2);
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

    result = result * 1.0L/12.0L;

    if(n % 2 == 1) {
        result = result - .25L;
    }

    return result;
}



double d_a(unsigned int n, unsigned int k) {
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

long double ld_a(unsigned int n, unsigned int k) {
    long double result;
    result = 0;

    if (k == 1) {
        return 1.0;
    }

    unsigned int h = 0;
    for(h = 1; h < k+1; h++) {
        if(GCD(h,k) == 1) {
        //    cout << "s(" << h << "," << k << ") = " << d_s(h,k) << endl;
        //    cout << "s(" << h << "," << k << ") = " << d_s(h,k) << endl;
            result += cosl( ld_pi * ( ld_s(h,k) - (2.0L * h * n)/(long double)(k)) );
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




void cospi (mpfr_t res,
		mpfr_t x)
{
//	mpfr_t t, tt, half, fourth;

//	mpfr_init2 (t, prec);
//	mpfr_init2 (tt, prec);
//	mpfr_init2 (half, prec);
//	mpfr_init2 (fourth, prec);

//	mpfr_set_ui (half, 1, r);
//	mpfr_div_2ui (half, half, 1, r);
//	mpfr_div_2ui (fourth, half, 1, r);

        // NOTE: switched t to tempc2
        //       and tt to tempc1


        mp_rnd_t r = round_mode;


	mpfr_div_2ui (tempc1, x, 1, r);
	mpfr_floor (tempc1, tempc1);
	mpfr_mul_2ui (tempc1, tempc1, 1, r);
	mpfr_sub (tempc2, x, tempc1, r);
	if (mpfr_cmp_ui (tempc2, 1) > 0)
		mpfr_sub_ui (tempc2, tempc2, 2, r);
	mpfr_abs (tempc1, tempc2, r);
	if (mpfr_cmp (tempc1, half) > 0)
	{
		mpfr_ui_sub (tempc2, 1, tempc1, r);
		mpfr_abs (tempc1, tempc2, r);
		if (mpfr_cmp (tempc1, fourth) > 0)
		{
			if (mpfr_sgn (tempc2) > 0)
				mpfr_sub (tempc2, half, tempc2, r);
			else
				mpfr_add (tempc2, tempc2, half, r);
			mpfr_mul (tempc2, tempc2, mp_pi, r);
			mpfr_sin (tempc2, tempc2, r);
		}
		else
		{
			mpfr_mul (tempc2, tempc2, mp_pi, r);
			mpfr_cos (tempc2, tempc2, r);
		}
		mpfr_neg (res, tempc2, r);
	}
	else
	{
		mpfr_abs (tempc1, tempc2, r);
		if (mpfr_cmp (tempc1, fourth) > 0)
		{
			if (mpfr_sgn (tempc2) > 0)
				mpfr_sub (tempc2, half, tempc2, r);
			else
				mpfr_add (tempc2, tempc2, half, r);
			mpfr_mul (tempc2, tempc2, mp_pi, r);
			mpfr_sin (res, tempc2, r);
		}
		else
		{
			mpfr_mul (tempc2, tempc2, mp_pi, r);
			mpfr_cos (res, tempc2, r);
		}
	}

//	mpfr_clear (half);
//	mpfr_clear (fourth);
//	mpfr_clear (t);
//	mpfr_clear (tt);
}



void mpz_part(mpz_t result, unsigned int n) {
    if(n == 1) {
        mpz_set_str(result, "1", 10);
        return;
    }

    mpfr_t mp_result;

    mp_t(mp_result, n);

    mpfr_get_z(result, mp_result, round_mode);

    mpfr_clear(mp_result);

    return;
}



int main(int argc, char *argv[]){
    //init();

    unsigned int n = 10;

    if(argc > 1)
        n = atoi(argv[1]);
    else {

        cout << test() << endl;
        return 0;
    }
    //mpfr_t result;

    //mp_t(result, n);

    mpz_t answer;
    mpz_init(answer);
    mpz_part(answer, n);

    //mpfr_get_z(answer, result, round_mode);

    mpz_out_str (stdout, 10, answer);

    cout << endl;

    return 0;
}


bool test(bool longtest) {
    // The values given below are confirmed by multiple sources, so are probably correct.
    // TODO: There should be some more code here to test that answers satisfy the proper congruences that the should
    //  satisfy. On the other hand, it might be better to test this file from within SAGE,
    //  since that might be easier, and would certainly be more in line with the rest of
    //  SAGE.

    mpz_t expected_value;
    mpz_t actual_value;

    mpz_init(expected_value);
    mpz_init(actual_value);

    // n = 1
    mpz_set_str(expected_value, "1", 10);
    mpz_part(actual_value, 1);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return false;

    // n = 10
    mpz_set_str(expected_value, "42", 10);
    mpz_part(actual_value, 10);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return false;

    // n = 100
    mpz_set_str(expected_value, "190569292", 10);
    mpz_part(actual_value, 100);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return false;

    // n = 1000
    mpz_set_str(expected_value, "24061467864032622473692149727991", 10);
    mpz_part(actual_value, 1000);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return false;

    // n = 10000
    mpz_set_str(expected_value, "36167251325636293988820471890953695495016030339315650422081868605887952568754066420592310556052906916435144", 10);
    mpz_part(actual_value, 10000);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return false;

    // n = 100000
    mpz_set_str(expected_value, "27493510569775696512677516320986352688173429315980054758203125984302147328114964173055050741660736621590157844774296248940493063070200461792764493033510116079342457190155718943509725312466108452006369558934464248716828789832182345009262853831404597021307130674510624419227311238999702284408609370935531629697851569569892196108480158600569421098519", 10);
    mpz_part(actual_value, 100000);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return false;

    // n = 1000000
    mpz_set_str(expected_value, "1471684986358223398631004760609895943484030484439142125334612747351666117418918618276330148873983597555842015374130600288095929387347128232270327849578001932784396072064228659048713020170971840761025676479860846908142829356706929785991290519899445490672219997823452874982974022288229850136767566294781887494687879003824699988197729200632068668735996662273816798266213482417208446631027428001918132198177180646511234542595026728424452592296781193448139994664730105742564359154794989181485285351370551399476719981691459022015599101959601417474075715430750022184895815209339012481734469448319323280150665384042994054179587751761294916248142479998802936507195257074485047571662771763903391442495113823298195263008336489826045837712202455304996382144601028531832004519046591968302787537418118486000612016852593542741980215046267245473237321845833427512524227465399130174076941280847400831542217999286071108336303316298289102444649696805395416791875480010852636774022023128467646919775022348562520747741843343657801534130704761975530375169707999287040285677841619347472368171772154046664303121315630003467104673818", 10);
    mpz_part(actual_value, 1000000);

    if(mpz_cmp(expected_value, actual_value) != 0)
        return false;




    mpz_clear(expected_value);
    mpz_clear(actual_value);

    return true;
}


/* answer must have already been mpz_init'd. */
int part(mpz_t answer, unsigned int n){
    mpfr_t result;

    mp_t(result, n);

    mpfr_get_z(answer, result, round_mode);

    mpfr_clear(result);
    return 0;
}
