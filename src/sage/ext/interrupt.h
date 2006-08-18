/*
Interrupt handling for extension code (header file)

Copyright (c) 2006, William Stein, under "The BSD License" (see below).
All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the
following conditions are met:

  * Redistributions of source code must retain the above
    copyright notice, this list of conditions and the
    following disclaimer.
  * Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials
    provided with the distribution.
  * Neither the name of the SAGE Foundation nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

/****************************************************************
  WARNING: Only modify the version of this file in the
  sage/ext/ directory.  It overwrites all other copies
  of this file during the SAGE build.
 ****************************************************************/

#ifndef FOO_H
#define FOO_H
#include <setjmp.h>
#include <Python.h>
#include <signal.h>

sigjmp_buf env;

/* In your Pyrex code, put

    _sig_on
    [pure c code block]
    _sig_off

   to get signal handling capabilities.

   VERY VERY IMPORTANT:
       1. These *must* always come in pairs. E.g., if you have just
          a _sig_off without a corresponding _sig_on, then ctrl-c
          later in the interpreter will sigfault!

       2. Do *not* put these in the __init__ method of a Pyrex extension
          class, or you'll get crashes.

       3. Do not put these around any non-pure C code!!!!

       4. If you need to do a check of control-c inside a loop, e.g.,
          when constructing a matrix, put _sig_check instead of using
          _sig_on then _sig_off.
*/


/* Unfortunately signal handler function
   declarations vary HUGELY from platform to platform:
*/


#if defined(__CYGWIN32__)     /* Windows XP */

  _sig_func_ptr  last, last_handler;

#elif defined(__FreeBSD__)    /* FreeBSD */

  sig_t last, last_handler;

#elif defined(__APPLE__)      /* OSX */

  sig_t last, last_handler;

#elif defined (__sun__) || defined (__sun)   /* Solaris */

  typedef void (*__sighandler_t )();
  __sighandler_t last, last_handler;

#else                                   /* Other, e.g., Linux */

__sighandler_t last, last_handler;

#endif


static PyObject *_m;

void sig_handle(int n);

static int __sig__n, __sig_str__n;

#define _sig_on last = signal(SIGINT, sig_handle); \
                if (last != sig_handle) \
                   last_handler = last; \
                signal(SIGALRM, sig_handle); \
                signal(SIGSEGV, sig_handle); \
                signal(SIGABRT, sig_handle); \
                signal(SIGFPE, sig_handle); \
                if (__sig__n = sigsetjmp(env,1)) {  \
                  if (__sig__n == SIGINT)   \
                      PyErr_SetString(PyExc_KeyboardInterrupt, ""); \
                  else  \
                    if (__sig__n == SIGALRM)  \
                      PyErr_SetString(PyExc_KeyboardInterrupt, "Alarm received"); \
                  else \
                      PyErr_SetString(PyExc_RuntimeError, ""); \
                  return(0); \
	       }

#define _sig_str(s) last = signal(SIGINT, sig_handle); \
                if (last != sig_handle) \
                   last_handler = last; \
                signal(SIGSEGV, sig_handle); \
                signal(SIGABRT, sig_handle); \
                signal(SIGFPE, sig_handle); \
                signal(SIGALRM, sig_handle); \
                if (__sig_str__n = sigsetjmp(env,1)) {  \
		  if(__sig_str__n == SIGINT ) \
                      PyErr_SetString(PyExc_KeyboardInterrupt, s); \
		  else \
                     if (__sig__n == SIGALRM) \
                      PyErr_SetString(PyExc_KeyboardInterrupt, s); \
                  else  \
                      PyErr_SetString(PyExc_RuntimeError, s); \
                  return(0); \
	       }

#define _sig_off signal(SIGINT, last_handler); signal(SIGALRM, last_handler); signal(SIGFPE, last_handler); signal(SIGABRT, last_handler);

#define _sig_check _sig_on _sig_off

#endif /* FOO_H */



/*
I thought maybe the following would work nicely, instead of just
returning, but no (it crashes sometimes!).  Also the line number is
not set correctly.  The only way I can imagine setting this correctly
is if the following bit of code is generated within pyrex, since pyrex
knows the line numbers.

	             PyErr_SetString(PyExc_RuntimeError, "(complete C-level traceback not available)"); \
		     __pyx_filename = __pyx_f[1]; __pyx_lineno = 0; goto __pyx_L1; \
                }

*/
