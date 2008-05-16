#ifndef __GIVARO_TABLESIZE_MAX__
#define __GIVARO_TABLESIZE_MAX__
// ==========================================================================
// file: givadicqfq.h
// Time-stamp: <11 Jun 07 19:18:31 Jean-Guillaume.Dumas@imag.fr>
// (c) Givaro Team
// date: 2007
// version:
// author: Jean-Guillaume.Dumas
// Description:
//   Zech extension fitting a small enough memory space
//   t-adic max sizes for BLAS based linear algebra over extension fields
// see:
// [Dumas, Gautier, Pernet 2002] Finite field linear algebra subroutines.
// ISSAC'02: Proceedings of the 2002 International Symposium on Symbolic
// and Algebraic Computation, Lille, France pp 63--74.
// ==========================================================================


#ifndef FF_TABLE_MAX
// 2^23 ---> 2^23*4*3 = 100K
// #define FF_TABLE_MAX 8388608UL
// 2^20 ---> 2s on 735MHz
//#define FF_TABLE_MAX 1048576UL
// Now 2^21+1 seems OK
#define FF_TABLE_MAX 2097153UL
#endif

#ifndef _GIVARO_FF_MAXEXPONENT_
#define _GIVARO_FF_MAXEXPONENT_ 21
#endif


#include <iostream>
#include <vector>
#include "givaro/givprimes16.h"

#include <math.h>
/* SAGE: Added because for some reason this is
   missing from /usr/include/math.h, even though
   it *should* be there.  -- William Stein*/
extern double logb _PARAMS((double));

#include <stddef.h>

  // ---------------------------------------------  class
class AdicSize {
public:

    static size_t nmax53(const unsigned long P, const unsigned long e) {
        size_t i = 0;
        while (Primes16::ith(i) < P) ++i;
        return n_max_53[i][e-2];
    }

    static size_t qmax53(const unsigned long P, const unsigned long e) {
        size_t i = 0;
        while (Primes16::ith(i) < P) ++i;
        return qadic_53[i][e-2];
    }

    static size_t nmax64(const unsigned long P, const unsigned long e) {
        size_t i = 0;
        while (Primes16::ith(i) < P) ++i;
        return n_max_64[i][e-2];
    }

    size_t qmax64(const unsigned long P, const unsigned long e) {
        size_t i = 0;
        while (Primes16::ith(i) < P) ++i;
        return qadic_64[i][e-2];
    }

    static size_t twopmax53(const unsigned long P, const unsigned long e, const unsigned long nm) {
        double tmp = double(P-1);
        tmp *= double(P-1);
        tmp *= double(e);
        tmp *= double(nm);
        size_t k = size_t(logb(tmp));
        return ( (53/(2*e-1))>k ? ++k : 0);
    }

    static size_t twopmax53(const unsigned long P, const unsigned long e) {
        size_t k = 53/(2*e-1);
        return ( pow(2.0,double(k))>double(e*(P-1)*(P-1)) ? k: 0);
    }

private:
    static const size_t n_max_53[][_GIVARO_FF_MAXEXPONENT_];
    static const size_t n_max_64[][_GIVARO_FF_MAXEXPONENT_];
    static const size_t qadic_53[][_GIVARO_FF_MAXEXPONENT_];
    static const size_t qadic_64[][_GIVARO_FF_MAXEXPONENT_];
};





#endif
