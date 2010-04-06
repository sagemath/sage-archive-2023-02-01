/*							cgammaf
 *
 *	Complex gamma function, single precision
 *
 *
 *
 * SYNOPSIS:
 *
 * #include <complex.h>
 * float complex x, y, cgammaf();
 *
 * y = cgammaf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns complex-valued gamma function of the complex argument.
 *
 * Arguments |Re(x)| <= 14 are increased by recurrence.
 * Large arguments are handled by Stirling's formula. Large negative
 * arguments are made positive using the reflection formula.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10,10      100000      5.2e-6      5.8e-7
 *    IEEE      -22,22      100000      1.1e-5      1.2e-6
 *    IEEE      -27,27      100000      4.6e-5      1.5e-6
 *
 * Error for arguments outside the test range will be larger
 * owing to error amplification by the exponential function.
 *
 */
/*							clgamf
 *
 *	Natural logarithm of complex gamma function, single precision
 *
 *
 *
 * SYNOPSIS:
 *
 * #include <complex.h>
 * float complex x, y, clgamf();
 *
 * y = clgamf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of the complex gamma
 * function of the argument.
 *
 * The logarithm of the gamma function is approximated by the
 * logarithmic version of Stirling's asymptotic formula.
 * Arguments of real part less than 8 are increased by recurrence.
 * The cosecant reflection formula is employed for arguments
 * having real part less than -8.
 *
 * Arguments greater than MAXLGM return MAXNUM and an error
 * message.  MAXLGM = 2.035093e36.
 *
 *
 *
 * ACCURACY:
 *
 *
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -20,20     9 000 000    2.6e-6       1.1e-7
 *    IEEE     -100,100      500 000                 7.5e-8
 * The error criterion was relative when the function magnitude
 * was greater than one but absolute when it was less than one.
 */

/*
Cephes Math Library Release 2.7:  March, 1998
Copyright 1984, 1998 Stephen L. Moshier
*/

#include "complex.h"
#include "mconf.h"

#define MAXGAM 171.624376956302725
static float LOGPI = 1.14472988584940017414f;

/* Stirling's formula for the gamma function */
#define NSTIR 7
static float STIR[ NSTIR ] = {
 -5.92166437353693882865E-4,
 6.97281375836585777429E-5,
 7.84039221720066627474E-4,
-2.29472093621399176955E-4,
-2.68132716049382716049E-3,
 3.47222222222222222222E-3,
 8.33333333333333333333E-2
};
#define MAXSTIR 26.0
static float SQTPI = 2.50662827463100050242f;

extern float MAXLOGF, MAXNUMF, PIF;
#ifdef ANSIPROT
extern float sinhf ( float x );
extern float coshf ( float x );
extern float sinf ( float x );
extern float cosf ( float x );
#else
float logf(), sinf(), polevlf(), p1evlf(), floorf(), fabsf();
float cosf(), sinhf(), coshf();
float complex cpowf(), cexpf(), cabsf();
#endif

/* Gamma function computed by Stirling's formula.  */

static float complex cstirf(x)
float complex x;
{
float complex y, v, w;
int i;

w = 1.0f/x;

y = STIR[0];
for (i = 1; i < NSTIR; i++)
  {
    y = y * w + STIR[i];
  }

w = 1.0f + w * y;
if( cabsf(x) > MAXSTIR )
	{ /* Avoid overflow in pow() */
	v = cpowf( x, 0.5f * x - 0.25f );
	y = v * (v * cexpf(-x));
	}
else
	{
	y = cpowf( x, x - 0.5f ) * cexpf(-x);
	}
/* y = cpowf( x, x - 0.5f ) * cexpf(-x); */

y = SQTPI * y * w;
return( y );
}



float complex cgammaf(x)
float complex x;
{
float p, q;
float complex c, u;
int k;

if( fabsf(crealf(x)) > 13.0f )
	{
	if( crealf(x) < 0.0f )
		{
		q = crealf(x);
		p = floorf(q);
		if((p == q ) && (cimag(x) == 0.0f))
			{
cgoverf:
			mtherr( "cgammaf", OVERFLOW );
			return( MAXNUMF + I * MAXNUMF);
			}
		/* Compute sin(pi x)  */
		k = q - 2.0f * floorf (0.5f * q);
		q = PIF * (q - p);
		p = PIF * cimag(x);
		c = sinf(q) * coshf(p) + cosf(q) * sinhf(p) * I;
		if (k & 1)
		  c = -c;
		/* Reflection formula.
		   ??? The denominator has been observed to overflow,
		   producing a NaN result on a SPARCstation instead of zero.
		   Example:
		   cgammaf(-19.425196335160404f - 18.188121066880587f*I) */
		c = PIF/(c * cgammaf(1.0f - x));
		}
	else
		{
		c = cstirf(x);
		}
	return( c );
	}

c = 1.0f;
p = 0.0f;
u = x;
while( crealf(u) < 13.0f )
	{
	if ((crealf(u) == 0.0f) && (cimagf(u) == 0.0f))
		goto cgoverf;
	c *= u;
	p += 1.0f;
	u = x + p;
	}
u = cstirf(u);
return( u / c );
}



/* Asymptotic expansion of log gamma  */
static float A[] = {
#if 0
-1.3924322169059011164274322169059011164274E0,
 1.7964437236883057316493849001588939669435E-1,
-2.9550653594771241830065359477124183006536E-2,
 6.4102564102564102564102564102564102564103E-3,
-1.9175269175269175269175269175269175269175E-3,
 8.4175084175084175084175084175084175084175E-4,
#endif
-5.9523809523809523809523809523809523809524E-4,
 7.9365079365079365079365079365079365079365E-4,
-2.7777777777777777777777777777777777777778E-3,
 8.3333333333333333333333333333333333333333E-2
};
/* log( sqrt( 2*pi ) ) */
static float LS2PI  =  0.91893853320467274178f;
#define MAXLGM 2.035093e36f



/* Logarithm of gamma function */

float complex clgamf(x)
float complex x;
{
float complex c, w, u, v;
float p, q, a;
int i, cj;

cj = 0;
if (cimagf(x) < 0)
  {
    cj = 1;
    x = conj(x);
  }

/* Reflection formula -z gamma(-z) gamma(z) = pi / sin(pi z)  */
if((crealf(x) < -7.0f) || (cimagf(x) < -7.0f))
	{
	q = crealf(x);
	p = floorf(q);
	if( p == q )
		goto loverf;
	if (fabsf(cimag(x)) > 18.4f)
	  {
	    /* sin z grows exponentially with Im(z).  Find ln sin(pi z)
	       from |sin z| = sqrt( sin^2 x + sinh^2 y),
               arg sin z = arctan(tanh y / tan x).  */
	    c = PIF * cimagf(x) - 0.6931471805599453094f
	      + I * PIF * (0.5f - q);
	    c = LOGPI - c - clgamf(1.0f - x);
	  }
	else
	  {
	    /* Reduce sine arg mod pi.  */
	    u = csinf( PIF * (x - p) );
	    if( u == 0.0f )
	      goto loverf;
	    w = clgamf(1.0f - x);
	    c = LOGPI - clogf( u ) - w;
	    /* Adjust for reduced sine arg.  */
	    //cimagf(c) += PIF * p;
	    (float)(__imag__ c = (cimagf(c) + PIF * p));
	  }
	goto ldone;
	}
w = 0.0f;
if( crealf(x) < 7.0f )
	{
	  /* To satisfy Im {clgam(z)} = arg cgamma(z), accumulate
	     arg u during the recurrence.  */
	  a = 0.0f;
	  w = 1.0f;
	  p = 0.0f;
	  u = x;
	  while( crealf(u) < 7.0f )
		{
		if( u == 0.0f )
			goto loverf;
		w *= u;
		a += cargf(u);
		p += 1.0f;
		u = x + p;
		}
	x = u;
	w = -logf(cabsf(w)) - I * a;
	}

if( crealf(x) > MAXLGM )
	{
loverf:
	mtherr( "clgamf", OVERFLOW );
	c = MAXNUMF + MAXNUMF * I;
	goto ldone;
	}

c = ( x - 0.5f ) * clogf(x) - x + LS2PI + w;

if( cabsf(x) > 1.0e8f )
  goto ldone;

v = 1.0f/(x*x);
u = A[0];
for (i = 1; i < 4; i++)
  {
    u = u * v + A[i];
  }
c = c + u / x;

ldone:
if (cj)
  c = conj(c);
return( c );
}
