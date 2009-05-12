/* emacs edit mode for this file is -*- C -*- */
/* $Id: assert.h,v 1.11 2009/01/06 13:58:20 Singular Exp $ */

/* This is for compatibility with standard assert.h */
#if defined (NDEBUG) && ! defined (NOASSERT)
#define NOASSERT
#endif

/* It should be possible to include this file multiple times for different */
/* settings of NOASSERT */

/* {{{ undefines */
#undef __ASSERT
#undef __ASSERT1
#undef STICKYASSERT
#undef STICKYASSERT1
#undef ASSERT
#undef ASSERT1

#undef __WARN
#undef STICKYWARN
#undef WARN

#undef PVIRT_VOID
#undef PVIRT_INTCF
#undef PVIRT_BOOL
#undef PVIRT_INT
#undef PVIRT_CHARCC
/* }}} */

#ifdef __cplusplus
#ifdef HAVE_CSTDIO
#include <cstdio>
#else
extern "C" {
#include <stdio.h>
#include <stdlib.h>
}
#endif
#else
#include <stdio.h>
#include <stdlib.h>
#endif

/* {{{ permanent macro definitions */
#ifndef __GNUC__
#define __ASSERT(expression, message, file, line) \
(fprintf( stderr, "error: " message "\n%s:%u: failed assertion `%s'\n", \
 file, line, expression ), abort(), 0 )
#define __ASSERT1(expression, message, parameter1, file, line)  \
(fprintf( stderr, "error: " message "\n%s:%u: failed assertion `%s'\n", \
 parameter1, file, line, expression ), abort(), 0 )

#define STICKYASSERT(expression, message) \
((void)((expression) ? 0 : __ASSERT(#expression, message, __FILE__, __LINE__)))
#define STICKYASSERT1(expression, message, parameter1) \
((void)((expression) ? 0 : __ASSERT1(#expression, message, parameter1, __FILE__, __LINE__)))

#define __WARN(expression, message, file, line)  \
(fprintf( stderr, "warning: " message "\n%s:%u: failed assertion `%s'\n", \
 file, line, expression ), 0 )
#define STICKYWARN(expression, message) \
((void)((expression) ? 0 : __WARN(#expression, message, __FILE__, __LINE__)))
#else /* __GNUCC__ */
/* use preprocessor macro __PRETTY_FUNCTION__ for more informative output */
#define __ASSERT(expression, message, file, line, function) \
(fprintf( stderr, "error: " message "\n%s:%u: In function `%s':\nfailed assertion `%s'\n", \
 file, line, function, expression ), abort(), 0 )
#define __ASSERT1(expression, message, parameter1, file, line, function)  \
(fprintf( stderr, "error: " message "\n%s:%u: In function `%s':\nfailed assertion `%s'\n", \
 parameter1, file, line, function, expression ), abort(), 0 )

#define STICKYASSERT(expression, message) \
((void)((expression) ? 0 : __ASSERT(#expression, message, __FILE__, __LINE__, __PRETTY_FUNCTION__)))
#define STICKYASSERT1(expression, message, parameter1) \
((void)((expression) ? 0 : __ASSERT1(#expression, message, parameter1, __FILE__, __LINE__, __PRETTY_FUNCTION__)))

#define __WARN(expression, message, file, line, function)  \
(fprintf( stderr, "warning: " message "\n%s:%u: In function `%s':\nfailed assertion `%s'\n", \
 file, line, function, expression ), 0 )
#define STICKYWARN(expression, message) \
((void)((expression) ? 0 : __WARN(#expression, message, __FILE__, __LINE__, __PRETTY_FUNCTION__)))
#endif /* __GNUCC__ */
/* }}} */

/* {{{ macro definitions dependent on NOASSERT */
#ifndef NOASSERT
#ifndef __GNUC__
#define ASSERT(expression, message) \
((void)((expression) ? 0 : __ASSERT(#expression, message, __FILE__, __LINE__)))
#define ASSERT1(expression, message, parameter1) \
((void)((expression) ? 0 : __ASSERT1(#expression, message, parameter1, __FILE__, __LINE__)))

#define WARN(expression, message) \
((void)((expression) ? 0 : __WARN(#expression, message, __FILE__, __LINE__)))
#else /* __GNUCC__ */
/* use preprocessor macro __PRETTY_FUNCTION__ for more informative output */
#define ASSERT(expression, message) \
((void)((expression) ? 0 : __ASSERT(#expression, message, __FILE__, __LINE__, __PRETTY_FUNCTION__)))
#define ASSERT1(expression, message, parameter1) \
((void)((expression) ? 0 : __ASSERT1(#expression, message, parameter1, __FILE__, __LINE__, __PRETTY_FUNCTION__)))

#define WARN(expression, message) \
((void)((expression) ? 0 : __WARN(#expression, message, __FILE__, __LINE__, __PRETTY_FUNCTION__)))
#endif /* __GNUCC__ */

#define PVIRT_VOID(msg) \
{ fprintf( stderr, "pure method( " msg " ) called\n" ); abort(); }
#define PVIRT_INTCF(msg) \
{ fprintf( stderr, "pure method( " msg " ) called\n" ); abort(); return 0; }
#define PVIRT_BOOL(msg) \
{ fprintf( stderr, "pure method( " msg " ) called\n" ); abort(); return false; }
#define PVIRT_INT(msg) \
{ fprintf( stderr, "pure method( " msg " ) called\n" ); abort(); return 0; }
#define PVIRT_CHARCC(msg) \
{ fprintf( stderr, "pure method( " msg " ) called\n" ); abort(); return 0; }
#else /* NOASSERT */
#define ASSERT(expression, message)
#define ASSERT1(expression, message, parameter1)

#define WARN(expression, message)

#define PVIRT_VOID(msg) = 0
#define PVIRT_INTCF(msg) = 0
#define PVIRT_BOOL(msg) = 0
#define PVIRT_INT(msg) = 0
#define PVIRT_CHARCC(msg) = 0
#endif /* NOASSERT */
/* }}} */


/* SAGE hack: This is *critically* needed to compile on some systems
   for which the logic above is broken.   Obviously it defeates the
   purpose of assert here, but at least it means that the code compiles.
*/
#ifndef assert
#define assert
#endif
