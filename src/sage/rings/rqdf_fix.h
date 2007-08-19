// defines NAN and INFINITY on solaris, etc.


#if defined(__sun)

#define NAN (0.0/0.0)
#define INFINITY __builtin_huge_valf()

#include <ieeefp.h>  # needed to define the function "infinite"

#endif
