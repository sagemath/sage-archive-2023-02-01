#include <math.h>

inline long double sage_logl(long double x)
{
#if defined(__CYGWIN__)
  return log(x);
#else
  return logl(x);
#endif
}

inline long double sage_sqrtl(long double x)
{
#if defined(__CYGWIN__)
  return sqrt(x);
#else
  return sqrtl(x);
#endif
}

inline long double sage_tgammal(long double x)
{
#if defined(__CYGWIN__)
  return tgamma(x);
#else
  return tgammal(x);
#endif
}

inline long double sage_lgammal(long double x)
{
#if defined(__CYGWIN__)
  return lgamma(x);
#else
  return lgammal(x);
#endif
}
