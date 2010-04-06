/*							cgammal.c
 *
 *	Complex gamma function, long double precision
 *
 *
 *
 * SYNOPSIS:
 *
 * #include <complex.h>
 * long double complex x, y, cgammal();
 *
 * y = cgammal( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns complex-valued gamma function of the complex argument.
 *
 * Arguments |Re(x)| <= 20 are increased by recurrence.
 * Large arguments are handled by Stirling's formula. Large negative
 * arguments are made positive using the reflection formula.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 * 80-bit long double:
 *    IEEE      -10,10       40000      4.1e-18     7.0e-19
 *    IEEE      -20,20       40000      9.1e-18     1.0e-18
 *    IEEE     -100,100      40000      5.9e-17     7.4e-18
 * 128-bit long double:
 *    IEEE      -10,10       30000      4.9e-32     8.7e-33
 *    IEEE     -100,100      45000      1.2e-31     1.7e-32
 *
 * Error for arguments outside the test range will be larger
 * owing to error amplification by the exponential function.
 *
 */
/*							clgaml()
 *
 *	Natural logarithm of complex gamma function, long double precision
 *
 *
 *
 * SYNOPSIS:
 *
 * #include <complex.h>
 * long double complex x, y, clgaml();
 *
 * y = clgaml( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of the complex gamma
 * function of the complex argument.
 *
 * The logarithm of the gamma function is approximated by the
 * logarithmic version of Stirling's asymptotic formula.
 * Arguments of real part less than 16 are increased by recurrence.
 * The cosecant reflection formula is employed for arguments
 * having real part less than -16.5.
 *
 * Arguments greater than MAXLGM return MAXNUM and an error
 * message.  MAXLGM = 1.048e+4928L.
 *
 *
 *
 * ACCURACY:
 *
 *
 * arithmetic   domain     # trials      peak         rms
 * 80-bit long double:
 *    IEEE      -10,10       30000     7.9e-18      5.2e-19
 *    IEEE      -50,50       20000                  1.1e-19
 *    IEEE     -100,100      20000                  7.4e-20
 * 128-bit long double:
 *    IEEE      -10,10       21000     4.4e-32      3.6e-33
 *    IEEE     -100,100      23000                  4.4e-34
 * The error criterion was relative when the function magnitude
 * was greater than one but absolute when it was less than one.
 */

/*
Cephes Math Library Release 2.7:  April, 1998
Copyright 1998 Stephen L. Moshier
*/

#include <complex.h>
#include "mconf.h"

#ifdef LD128BITS
#define NGITER 50.0L
#define NLGITER 50.0L
#define LGMXINT 44.4L
#define GSMALL 1.e-17L
#else
#define NGITER 20.0L
#define NLGITER 16.0L
#define LGMXINT 78.3L
#define GSMALL 1.e-9L
#endif

#define MAXGAM  1755.455L
static long double LOGPIL =  1.1447298858494001741434273513530587116473L;


/* Stirling's formula for the gamma function */
#define NSTIR 18
static long double STIR[NSTIR] = {
 1.50561130400264244123842218771311273E-2L,
 1.79540117061234856107699407722226331E-1L,
-2.48174360026499773091565836874346432E-3L,
-2.95278809456991205054406510546938244E-2L,
 5.40164767892604515180467508570241736E-4L,
 6.40336283380806979482363809026579583E-3L,
-1.62516262783915816898635123980270998E-4L,
-1.91443849856547752650089885832852254E-3L,
 7.20489541602001055908571930225015052E-5L,
 8.39498720672087279993357516764983445E-4L,
-5.17179090826059219337057843002058823E-5L,
-5.92166437353693882864836225604401187E-4L,
 6.97281375836585777429398828575783308E-5L,
 7.84039221720066627474034881442288850E-4L,
-2.29472093621399176954732510288065844E-4L,
-2.68132716049382716049382716049382716E-3L,
 3.47222222222222222222222222222222222E-3L,
 8.33333333333333333333333333333333333E-2L,
};
#define MAXSTIR 1024.0L
static long double SQTPIL = 2.50662827463100050241576528481104525L;

extern long double MAXLOGL, MAXNUML, PIL;
#ifdef ANSIPROT
extern long double sinhl ( long double x );
extern long double coshl ( long double x );
extern long double sinl ( long double x );
extern long double cosl ( long double x );
#else
long double logl(), sinl(), polevll(), p1evll(), floorl(), fabsl();
long double sinhl(), coshl(), cosl();
long double complex cpowl(), cexpl(), cabsl();
#endif

/* Gamma function computed by Stirling's formula.  */

/* static double complex cstirf(x) */
long double complex cstirfl(x)
long double complex x;
{
long double complex y, w;
int i;

w = 1.0L/x;

y = STIR[0];
for (i = 1; i < NSTIR; i++)
  {
    y = y * w + STIR[i];
  }

w = 1.0L + w * y;
#if 1
y = cpowl( x, x - 0.5L ) * cexpl(-x);
#else
y = (x - 0.5L) * clogl(x) - x;
y = cexpl(y);
#endif
y = SQTPIL * y * w;
return( y );
}



long double complex cgammal(x)
long double complex x;
{
long double complex c, u;
long double p, q;
int cj, k;

cj = 0;
if (cimagl(x) < 0.0L)
  {
    cj = 1;
    x = conj(x);
  }

if( fabsl(creall(x)) > NGITER )
	{
	if( creall(x) < 0.0L )
		{
		q = creall(x);
		p = floorl(q);
		if(( p == q ) && (cimagl(x) == 0.0L))
			{
			mtherr( "cgammal", OVERFLOW );
			c = MAXNUML + I * MAXNUML;
			goto gamdone;
			}
		/*	c = csinl( PIL * x );*/
		/* Compute sin(pi x)  */
		k = q - 2.0L * floorl (0.5L * q);
		q = PIL * (q - p);
		p = PIL * cimagl(x);
		c = sinl(q) * coshl(p) + cosl(q) * sinhl(p) * I;
		if (k & 1)
		  c = -c;
		c = PIL/(c * cgammal(1.0L - x) );
		goto gamdone;
		}
	else
		{
		  c = cstirfl(x);
		  goto gamdone;
		}
	}
c = 1.0L;
p = 0.0L;
u = x;
while( creall(u) < NGITER )
	{
	if((fabsl (creall(u)) < GSMALL) && (fabsl (cimagl(u)) < GSMALL))
		goto small;
	c *= u;
	p += 1.0L;
	u = x + p;
	}
u = cstirfl(u);
c = u / c;
goto gamdone;


small:
if((creall(x) == 0.0L) && (cimagl(x) == 0.0L))
	{
	mtherr( "cgammal", SING );
	c = MAXNUML + MAXNUML * I;
	goto gamdone;
	}
else
	c = 1.0L/(((1.0L + 0.57721566490153286060651209008240243L * u) * u)*c);

gamdone:

if (cj)
  c = conj(c);
return( c );
}


/* Asymptotic expansion of log gamma  */
#define NUMA 9
static long double A[NUMA] = {
#if 0
  1.3402864044168391994478951000690131124914E1L,
 -1.3924322169059011164274322169059011164274E0L,
#endif
 1.7964437236883057316493849001588939669435E-1L,
-2.9550653594771241830065359477124183006536E-2L,
 6.4102564102564102564102564102564102564103E-3L,
-1.9175269175269175269175269175269175269175E-3L,
 8.4175084175084175084175084175084175084175E-4L,
-5.9523809523809523809523809523809523809524E-4L,
 7.9365079365079365079365079365079365079365E-4L,
-2.7777777777777777777777777777777777777778E-3L,
 8.3333333333333333333333333333333333333333E-2L
};
/* log( sqrt( 2*pi ) ) */
static long double LS2PIL  =  0.918938533204672741780329736405617639861397L;
#define MAXLGML 1.04848146839019521116e+4928L



/* Logarithm of gamma function */

long double complex clgaml(x)
long double complex x;
{
long double complex c, w, u, v;
long double p, q, a;
int i, cj;

cj = 0;
if (cimagl(x) < 0.0L)
  {
    cj = 1;
    x = conj(x);
  }

/* -z gamma(-z) gamma(z) = pi / sin(pi z) */
/* log gamma(z) = log pi - log sin(pi z) - log(-z) - log gamma(-z) */
if((creall(x) < -NLGITER) || (cimagl(x) < -NLGITER))
	{
	q = creall(x);
	p = floorl(q);
	if( p == q )
		goto loverf;
	if (fabsl(cimagl(x)) > LGMXINT)
	  {
	    /* sin z grows exponentially with Im(z).  Find ln sin(pi z)
	       from |sin z| = sqrt( sin^2 x + sinh^2 y),
               arg sin z = arctan(tanh y / tan x).  */
	    c = PIL * cimagl(x) - 0.69314718055994530941723212145817656807550L
	      + I * PIL * (0.5L - q);
	    c = LOGPIL - c - clgaml(1.0L - x);
	  }
	else
	  {
	    /* Reduce sine arg mod pi.  */
	    u = csinl( PIL * (x - p) );
	    if( u == 0.0L )
	      goto loverf;
	    w = clgaml(1.0L - x);
	    c = LOGPIL - clogl( u ) - w;
	    /* Adjust for reduced sine arg.  */
	    //cimagl(c) += PIL * p;
	    (long double)(__imag__ c = cimagl(c) + PIL * p);
	  }
	goto ldone;
	}
w = 0.0L;
if( creall(x) < NLGITER )
	{
	  /* To satisfy Im {clgam(z)} = arg cgamma(z), accumulate
	     arg u during the recurrence.  */
	  a = 0.0L;
	  w = 1.0L;
	  p = 0.0L;
	  u = x;
	  while( creall(u) < NLGITER )
		{
		if( u == 0.0L )
			goto loverf;
		w *= u;
		a += cargl(u);
		p += 1.0L;
		u = x + p;
		}
	x = u;
	/*	w = -logl(cabsl(w)) - I * a; */
	p = creall(w);
	q = cimagl(w);
	w = -0.5 * logl(p*p + q*q) - I * a;
	}

if( creal(x) > MAXLGML )
	{
loverf:
	mtherr( "clgaml", OVERFLOW );
	c = MAXNUML + MAXNUML * I;
	goto ldone;
	}

c = ( x - 0.5L ) * clogl(x) - x + LS2PIL + w;

if( cabsl(x) > 1.0e10L )
  goto ldone;

v = 1.0L/(x*x);
u = A[0];
for (i = 1; i < NUMA; i++)
  {
    u = u * v + A[i];
  }
c = c + u / x;

ldone:
if (cj)
  c = conj(c);
return( c );
}
