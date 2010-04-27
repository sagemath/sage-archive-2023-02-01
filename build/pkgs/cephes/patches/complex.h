/* This is a C9X complex.h implementation for the GNU C compiler.
   It uses GNU extensions such as __complex__, but defines them
   so the syntax specified by C9X works.

   S. L. Moshier
   February, 1997  */


/* This may look like something else, but it is really C9X.  */

/* Use `I' as the imaginary unit.
   Initialize complex constants like this:
   double complex c = 1.0 + 2.0 * I;  */
#define _Imaginary_I (1.0fi)
#define _Complex_I (1.0fi)
#define I _Imaginary_I

/* These are not typed in gcc.  They preseve the type of the argument.  */
/* Complex conjugate function.  */
#define conj(x) (~(x))
/* Function to get imaginary part.  */
#define cimag(x) (__imag__ (x))
/* Function to get real part.  */
#define creal(x) (__real__ (x))

extern double atan2 (double, double);
#define carg(z) (atan2 ((double) cimag (z), (double) creal (z)))

/* This is the new C language key word.
   ... Actually, the key word that the compiler is supposed to reserve
   and understand as part of the language is spelled "_Complex" but
   gcc understands "__complex__."  The macro "complex" is supposed to
   exist and expand to _Complex in complex.h.  But "_Complex" will seldom
   if ever appear in user programs.  */
#define complex __complex__
/* Remove this whenever gcc implements _Complex as a key word.  */
#define _Complex __complex__

/* These pragmas can't work without a compiler modification.  */
#define CX_LIMITED_RANGE_ON
#define CX_LIMITED_RANGE_OFF
#define CX_LIMITED_RANGE_DEFAULT CX_LIMITED_RANGE_ON

/* The builtin complex absolute value in gcc currently is unreachable
   due to overzealous type checking.  */
extern double cabs ( double complex z );
extern double hypot ( double x, double y );

/* Prototypes for clog.c.
   This is how you declare complex things in the new C language.  */
extern double complex clog ( double complex z );
extern double complex cexp ( double complex z );
extern double complex csqrt ( double complex z );
extern double complex csin ( double complex z );
extern double complex ccos ( double complex z );
extern double complex ctan ( double complex z );
extern double complex ccot ( double complex z );
extern double complex casin ( double complex z );
extern double complex cacos ( double complex z );
extern double complex catan ( double complex z );
extern double complex csinh ( double complex z );
extern double complex casinh ( double complex z );
extern double complex ccosh ( double complex z );
extern double complex cacosh ( double complex z );
extern double complex ctanh ( double complex z );
extern double complex catanh ( double complex z );
extern double complex cpow (double complex a, double complex z);

/* These functions might be used if the compiler were to generate
   subroutine calls.  But their names would be spelled some other way.  */
extern double complex cadd ( double complex a, double complex b );
extern double complex csub ( double complex a, double complex b );
extern double complex cmul ( double complex a, double complex b );
extern double complex cdiv ( double complex a, double complex b );

/* There are float and long double sizes, too.  */
#define cimagf(x) ((float) __imag__ (x))
#define crealf(x) ((float) __real__ (x))
extern float atan2f (float, float);
#define cargf(z) (atan2f ((float) cimag (z), (float) creal (z)))
extern float cabsf ( float complex z );
extern float complex clogf ( float complex z );
extern float complex cexpf ( float complex z );
extern float complex csqrtf ( float complex z );
extern float complex csinf ( float complex z );
extern float complex ccosf ( float complex z );
extern float complex ctanf ( float complex z );
extern float complex ccotf ( float complex z );
extern float complex casinf ( float complex z );
extern float complex cacosf ( float complex z );
extern float complex catanf ( float complex z );
extern float complex csinhf ( float complex z );
extern float complex casinhf ( float complex z );
extern float complex ccoshf ( float complex z );
extern float complex cacoshf ( float complex z );
extern float complex ctanhf ( float complex z );
extern float complex catanhf ( float complex z );
extern float complex cpowf (float complex a, float complex z);

#define cimagl(x) ((long double) __imag__ (x))
#define creall(x) ((long double) __real__ (x))
extern long double atan2l (long double, long double);
#define cargl(z) (atan2l ((long double) cimag (z), (long double) creal (z)))
extern long double cabsl ( long double complex z );
extern long double complex clogl ( long double complex z );
extern long double complex cexpl ( long double complex z );
extern long double complex csqrtl ( long double complex z );
extern long double complex csinl ( long double complex z );
extern long double complex ccosl ( long double complex z );
extern long double complex ctanl ( long double complex z );
extern long double complex ccotl ( long double complex z );
extern long double complex casinl ( long double complex z );
extern long double complex cacosl ( long double complex z );
extern long double complex catanl ( long double complex z );
extern long double complex csinhl ( long double complex z );
extern long double complex casinhl ( long double complex z );
extern long double complex ccoshl ( long double complex z );
extern long double complex cacoshl ( long double complex z );
extern long double complex ctanhl ( long double complex z );
extern long double complex catanhl ( long double complex z );
extern long double complex cpowl (long double complex a, long double complex z);
extern float complex clgamf ( float complex z );
extern double complex clgam ( double complex z );
extern long double complex clgaml ( long double complex z );
extern float complex cgammaf ( float complex z );
extern double complex cgamma ( double complex z );
extern long double complex cgammal ( long double complex z );
