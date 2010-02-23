
#if defined(__sun)
typedef int int_fast32_t;
typedef long long int_fast64_t;
#else
#include <stdint.h>
#endif

#define INTEGER_MOD_INT32_LIMIT 46341          //  = ceil(sqrt(2^31-1))
#define INTEGER_MOD_INT64_LIMIT 2147483647     //  = 2^31-1 for now, should be 3037000500LL = ceil(sqrt(2^63-1))

