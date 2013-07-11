/* See mod_int.pxd if you want to use these declarations in Cython */

/* CPython using signed long as the data type for plain integers.  We
 * use the same here so that conversion to/from Python int is fast. If
 * this were an unsigned long then it would be converted to
 * arbitrary-precision Python PyLong.
 */
#define mod_int long

/* The largest value we can do arithmetic on without risking overflowing.
 * That is, you can multiply two MOD_INT_MAX in a mod_int.
 */
#define MOD_INT_OVERFLOW ((((mod_int)1) << (sizeof(mod_int)*8 - 1)) - 1)
#define MOD_INT_MAX ((mod_int)sqrt(MOD_INT_OVERFLOW) - 1)

