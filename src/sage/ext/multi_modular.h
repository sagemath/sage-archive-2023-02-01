#define mod_int size_t


/* This is the largest value we can do arithmatic on without risking overflowing.
 * It will be optimized to a constant if any optimization is turned on.
 * The expression below assumes unsigned.
 */
#define MOD_INT_MAX ((1 << (sizeof(mod_int)*4)) - 1)
#define MOD_INT_OVERFLOW ((1 << (sizeof(mod_int)*8)) - 1)

/* Don't want to be up against the limit so
 * we have room to gather.
 */
#define START_PRIME_MAX (MOD_INT_MAX/5)