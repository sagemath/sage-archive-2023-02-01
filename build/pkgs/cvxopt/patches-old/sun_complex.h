/*
 * Copyright 2004 Sun Microsystems, Inc.  All rights reserved.
 * Use is subject to license terms.
 */

#ifndef _COMPLEX_H
#define	_COMPLEX_H

#pragma ident	"@(#)complex.h	1.9	04/10/23 SMI"

#if !defined(__cplusplus)

/*
 * Compilation environments for Solaris must provide the _Imaginary datatype
 * and the compiler intrinsics _Complex_I and _Imaginary_I
 */
/* ISO C99 has types `float complex' and `double complex'.  */
#define	complex		_Complex

/* Use the smaller imaginary unit `float complex' constant  1.0fi,
   instead of `double complex' constant 1.0i.  */
#define        _Complex_I      (__extension__ 1.0fi)

/* Another more descriptive name is `I'.  */
#undef	I

/* Once there is imaginary support in gcc, use it.  */
#define        I               _Complex_I

typedef        float _Complex  _Imaginary;
#define        imaginary       _Imaginary
#define        _Imaginary_I    (__extension__ 1.0fi)

extern float cabsf(float complex);
extern float cargf(float complex);
extern float cimagf(float complex);
extern float crealf(float complex);
extern float complex cacosf(float complex);
extern float complex cacoshf(float complex);
extern float complex casinf(float complex);
extern float complex casinhf(float complex);
extern float complex catanf(float complex);
extern float complex catanhf(float complex);
extern float complex ccosf(float complex);
extern float complex ccoshf(float complex);
extern float complex cexpf(float complex);
extern float complex clogf(float complex);
extern float complex conjf(float complex);
extern float complex cpowf(float complex, float complex);
extern float complex cprojf(float complex);
extern float complex csinf(float complex);
extern float complex csinhf(float complex);
extern float complex csqrtf(float complex);
extern float complex ctanf(float complex);
extern float complex ctanhf(float complex);

extern double cabs(double complex);
extern double carg(double complex);
extern double cimag(double complex);
extern double creal(double complex);
extern double complex cacos(double complex);
extern double complex cacosh(double complex);
extern double complex casin(double complex);
extern double complex casinh(double complex);
extern double complex catan(double complex);
extern double complex catanh(double complex);
extern double complex ccos(double complex);
extern double complex ccosh(double complex);
extern double complex cexp(double complex);
#if defined(__PRAGMA_REDEFINE_EXTNAME)
#pragma redefine_extname clog __clog
#else
#undef	clog
#define	clog	__clog
#endif
extern double complex clog(double complex);
extern double complex conj(double complex);
extern double complex cpow(double complex, double complex);
extern double complex cproj(double complex);
extern double complex csin(double complex);
extern double complex csinh(double complex);
extern double complex csqrt(double complex);
extern double complex ctan(double complex);
extern double complex ctanh(double complex);

extern long double cabsl(long double complex);
extern long double cargl(long double complex);
extern long double cimagl(long double complex);
extern long double creall(long double complex);
extern long double complex cacoshl(long double complex);
extern long double complex cacosl(long double complex);
extern long double complex casinhl(long double complex);
extern long double complex casinl(long double complex);
extern long double complex catanhl(long double complex);
extern long double complex catanl(long double complex);
extern long double complex ccoshl(long double complex);
extern long double complex ccosl(long double complex);
extern long double complex cexpl(long double complex);
extern long double complex clogl(long double complex);
extern long double complex conjl(long double complex);
extern long double complex cpowl(long double complex, long double complex);
extern long double complex cprojl(long double complex);
extern long double complex csinhl(long double complex);
extern long double complex csinl(long double complex);
extern long double complex csqrtl(long double complex);
extern long double complex ctanhl(long double complex);
extern long double complex ctanl(long double complex);

#endif	/* !defined(__cplusplus) */

#endif	/* _COMPLEX_H */
