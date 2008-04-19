#ifndef _GMPplusplus_H_
#define _GMPplusplus_H_
// ========================================================================
// Givaro version of gmp++.h
// Time-stamp: <08 Mar 07 14:07:31 Jean-Guillaume.Dumas@imag.fr>
// ========================================================================
#include <climits> // required by gcc 4.3
#include <givaro-config.h>

#ifndef __GIVARO__DONOTUSE_longlong__
#ifndef __DONOTUSE_Givaro_SIXTYFOUR__
#define __USE_64_bits__
#endif
#endif

#if !defined(GMP_VERSION_3) && !defined(GMP_NO_CXX) && !defined(__GIVARO_GMP_VERSION_3) && !defined(__GIVARO_GMP_NO_CXX)
// gmpxx.h defines __GMP_PLUSPLUS__
#include <gmpxx.h>
#endif

#ifdef __GIVARO_GMP_VERSION_3
extern "C" {
#endif

#include "gmp.h"

#ifdef __GIVARO_GMP_VERSION_3
}
#endif


#include <gmp++/gmp++_int.h>

#endif
