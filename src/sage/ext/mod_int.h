/* See mod_int.pxd if you want to use these declarations in Cython */

/* CPython using signed long as the data type for plain integers.  We
 * use the same here so that conversion to/from Python int is fast. If
 * this were an unsigned long then it would be converted to
 * arbitrary-precision Python PyLong.
 */
#define mod_int long

/* This is the largest value we can do arithmetic on without risking overflowing.
 * It will be optimized to a constant if any optimization is turned on.
 * The expression below assumes signed.
 */
#define MOD_INT_MAX ((1 << (sizeof(mod_int)*4 - 1)) - 1)
#define MOD_INT_OVERFLOW ((1 << (sizeof(mod_int)*8 - 1)) - 1)


