#ifndef _GMPplusplus_H_
#define _GMPplusplus_H_
// ========================================================================
// Givaro version of gmp++.h
// Time-stamp: <19 Dec 06 10:51:44 Jean-Guillaume.Dumas@imag.fr>
// ========================================================================
#ifndef __GIVARO__DONOTUSE_longlong__
#ifndef __DONOTUSE_Givaro_SIXTYFOUR__
#define __USE_64_bits__
#endif
#endif

#if !defined(GMP_NO_CXX) && !defined(__GIVARO_GMP_VERSION_3) && !defined(__GIVARO_GMP_NO_CXX)
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
