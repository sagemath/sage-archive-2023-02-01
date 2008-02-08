// defines NAN and INFINITY on solaris, etc.


#if defined(__sun)

#define NAN (0.0/0.0)
#define INFINITY __builtin_huge_valf()

#include <ieeefp.h>  # needed to define the function "infinite"

#define roundf(a) ((fmod(a,1)<0.5)?floor(a):ceil(a))

int isinf(double x) { return !finite(x) && x==x; }

#endif
