
#ifndef NTL_new__H
#define NTL_new__H

#include <NTL/config.h>

#if (defined(NTL_STD_CXX) || defined(NTL_PSTD_NTN))

// We use <new> and std::nothrow, even if neither NTL_STD_CXX nor
// NTL_PSTD_NHF are set.  This appears to be somewhat more compatible
// with current compilers.

#include <new>

// uncommenting std::nothrow makes this ntl work properly with Singular
#define NTL_NEW_OP new //(std::nothrow)


#else

#define NTL_NEW_OP new

#endif

#endif
