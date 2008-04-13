
#ifndef NTL_tools__H
#define NTL_tools__H

#include <NTL/ctools.h>

#if (defined(NTL_STD_CXX) || defined(NTL_PSTD_NHF))

// new header files

#include <cstdlib>
#include <cmath>
#include <iostream>

#else

// old header files

#include <stdlib.h>
#include <math.h>
#include <iostream.h>

#endif

#if (defined(NTL_STD_CXX) || defined(NTL_PSTD_NHF))

#define NTL_SNS std ::
#define NTL_USE_SNS using namespace std;

#elif (defined(NTL_PSTD_NNS))

#define NTL_SNS ::
#define NTL_USE_SNS

#else

#define NTL_SNS
#define NTL_USE_SNS

#endif

#if (defined(NTL_STD_CXX) || defined(NTL_PSTD_NNS))

#define NTL_NAMESPACE NTL
#define NTL_OPEN_NNS namespace NTL_NAMESPACE {
#define NTL_CLOSE_NNS  }
#define NTL_USE_NNS using namespace NTL_NAMESPACE;
#define NTL_NNS NTL_NAMESPACE ::

// To make things work, we have to apply using declarations of all std
// functions that are both overloaded by NTL and are used in
// the implementation of NTL.

#define NTL_START_IMPL NTL_USE_SNS NTL_OPEN_NNS \
   using NTL_SNS abs; \
   using NTL_SNS ceil; \
   using NTL_SNS exp; \
   using NTL_SNS fabs; \
   using NTL_SNS floor; \
   using NTL_SNS ldexp; \
   using NTL_SNS log; \
   using NTL_SNS sqrt;

#define NTL_END_IMPL NTL_CLOSE_NNS

#else

#define NTL_NAMESPACE
#define NTL_OPEN_NNS
#define NTL_CLOSE_NNS
#define NTL_USE_NNS
#define NTL_NNS

#define NTL_START_IMPL
#define NTL_END_IMPL

#endif

#define NTL_CLIENT NTL_USE_SNS NTL_USE_NNS



#if 0

// This is for debugging purposes only.

namespace foo_bar {

class ostream;
class istream;

typedef unsigned int size_t;

double floor(double);
float floor(float);

}
#endif



#if (!defined(NTL_CXX_ONLY))
extern "C"
#endif
double _ntl_GetTime();

typedef unsigned long _ntl_ulong;
typedef _ntl_ulong *_ntl_ulong_ptr;
// I made these have "obscure" names to avoid conflict with
// (non-standard but common) definitions in standard headers.
// Putting u_long inside namespace NTL only tends to creates ambiguities,
// for no good reason.



NTL_OPEN_NNS

struct INIT_SIZE_STRUCT { };
const INIT_SIZE_STRUCT INIT_SIZE = INIT_SIZE_STRUCT();
typedef const INIT_SIZE_STRUCT& INIT_SIZE_TYPE;

struct INIT_VAL_STRUCT { };
const INIT_VAL_STRUCT INIT_VAL = INIT_VAL_STRUCT();
typedef const INIT_VAL_STRUCT& INIT_VAL_TYPE;

struct INIT_TRANS_STRUCT { };
const INIT_TRANS_STRUCT INIT_TRANS = INIT_TRANS_STRUCT();
typedef const INIT_TRANS_STRUCT& INIT_TRANS_TYPE;


struct INIT_LOOP_HOLE_STRUCT { };
const INIT_LOOP_HOLE_STRUCT INIT_LOOP_HOLE = INIT_LOOP_HOLE_STRUCT();
typedef const INIT_LOOP_HOLE_STRUCT& INIT_LOOP_HOLE_TYPE;

struct INIT_FFT_STRUCT { };
const INIT_FFT_STRUCT INIT_FFT = INIT_FFT_STRUCT();
typedef const INIT_FFT_STRUCT& INIT_FFT_TYPE;


#ifdef NTL_NO_INIT_TRANS
#define NTL_OPT_RETURN(t, x) return x
#else
#define NTL_OPT_RETURN(t, x) return t(x, INIT_TRANS)
#endif


#ifndef NTL_NO_MIN_MAX

inline int min(int a, int b) { return (a < b) ?  a : b; }
inline int max(int a, int b) { return (a < b) ? b : a; }

inline long min(long a, long b) { return (a < b) ?  a : b; }
inline long max(long a, long b) { return (a < b) ? b : a; }

inline long min(int a, long b) { return (a < b) ?  long(a) : b; }
inline long max(int a, long b) { return (a < b) ? b : long(a); }

inline long min(long a, int b) { return (a < b) ?  a : long(b); }
inline long max(long a, int b) { return (a < b) ? long(b) : a; }

#endif


inline void swap(long& a, long& b)  {  long t;  t = a; a = b; b = t; }
inline void swap(int& a, int& b)  {  int t;  t = a; a = b; b = t; }



inline void conv(int& x, int a) { x = a; }
inline void conv(int& x, long a)
   { unsigned y = (unsigned) a;  x = NTL_UINT_TO_INT(y); }
inline void conv(int& x, float a) { x = int(NTL_SNS floor(double(a))); }
inline void conv(int& x, double a) { x = int(NTL_SNS floor(a)); }

inline void conv(int& x, unsigned a)
   { x = NTL_UINT_TO_INT(a); }

inline void conv(int& x, unsigned long a)
   { unsigned y = (unsigned) a;  x = NTL_UINT_TO_INT(y); }

inline int to_int(int a) { return a; }
inline int to_int(long a)
   { unsigned y = (unsigned) a;  return NTL_UINT_TO_INT(y); }
inline int to_int(float a) { return int(NTL_SNS floor(double(a))); }
inline int to_int(double a) { return int(NTL_SNS floor(a)); }

inline int to_int(unsigned a)
   { return NTL_UINT_TO_INT(a); }

inline int to_int(unsigned long a)
   { unsigned y = (unsigned) a;  return NTL_UINT_TO_INT(y); }


inline void conv(long& x, int a) { x = a; }
inline void conv(long& x, long a) { x = a; }
inline void conv(long& x, float a) { x = long(NTL_SNS floor(double(a))); }
inline void conv(long& x, double a) { x = long(NTL_SNS floor(a)); }

inline void conv(long& x, unsigned a)
   { unsigned long y = a;  x = NTL_ULONG_TO_LONG(y); }

inline void conv(long& x, unsigned long a)
   { x = NTL_ULONG_TO_LONG(a); }

inline long to_long(int a) { return a; }
inline long to_long(long a) { return a; }
inline long to_long(float a) { return long(NTL_SNS floor(double(a))); }
inline long to_long(double a) { return long(NTL_SNS floor(a)); }

inline long to_long(unsigned a)
   { unsigned long y = a;  return NTL_ULONG_TO_LONG(y); }

inline long to_long(unsigned long a)
   { return NTL_ULONG_TO_LONG(a); }

inline void conv(float& x, int a) { x = float(a); }
inline void conv(float& x, long a) { x = float(a); }
inline void conv(float& x, unsigned a) { x = float(a); }
inline void conv(float& x, unsigned long a) { x = float(a); }
inline void conv(float& x, float a) { x = a; }
inline void conv(float& x, double a) { x = float(a); }

inline float to_float(int a) { return float(a); }
inline float to_float(long a) { return float(a); }
inline float to_float(unsigned a) { return float(a); }
inline float to_float(unsigned long a) { return float(a); }
inline float to_float(float a) { return a; }
inline float to_float(double a) { return float(a); }

inline void conv(double& x, int a) { x = double(a); }
inline void conv(double& x, long a) { x = double(a); }
inline void conv(double& x, unsigned a) { x = double(a); }
inline void conv(double& x, unsigned long a) { x = double(a); }
inline void conv(double& x, float a) { x = double(a); }
inline void conv(double& x, double a) { x = a; }

inline double to_double(int a) { return double(a); }
inline double to_double(long a) { return double(a); }
inline double to_double(unsigned a) { return double(a); }
inline double to_double(unsigned long a) { return double(a); }
inline double to_double(float a) { return double(a); }
inline double to_double(double a) { return a; }


long SkipWhiteSpace(NTL_SNS istream& s);
long IsWhiteSpace(long c);

long CharToIntVal(long c);
char IntValToChar(long a);


/*
  This function is not present in vanilla NTL 5.4.2.
  See tools.c for documentation.
 */
void SetErrorCallbackFunction(void (*func)(const char *s, void *context), void *context);


void Error(const char *s);


inline double GetTime() { return _ntl_GetTime(); }

inline long IsFinite(double *p) { return _ntl_IsFinite(p); }


#if (NTL_EXT_DOUBLE)

inline void ForceToMem(double *p) { _ntl_ForceToMem(p); }

#else

inline void ForceToMem(double *p) { }

#endif



void PrintTime(NTL_SNS ostream& s, double t);

NTL_CLOSE_NNS


#endif

