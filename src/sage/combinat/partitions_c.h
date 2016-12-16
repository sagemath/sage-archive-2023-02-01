#if defined(__sun)
#include <ieeefp.h>
int isinf(double x) { return !finite(x) && x==x; }
#endif

#include <gmp.h>

int part(mpz_t answer, unsigned int n);
int test(bool longtest, bool forever);

