/******************************************************************************
       Copyright (C) 2006 William Stein <wstein@gmail.com>
                     2006 Martin Albrecht <malb@informatik.uni-bremen.de>

  Distributed under the terms of the GNU General Public License (GPL), Version 2.

  The full text of the GPL is available at:
                  http://www.gnu.org/licenses/

******************************************************************************/

#ifndef FOO_H
#define FOO_H
#include <setjmp.h>
#include <Python.h>
#include <signal.h>



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

/**
 * bitmask to enable the sage signal handler
 */

#define SIG_SAGE_HANDLER 1

/**
 * bitmask to enable long jumps
 */
#define SIG_LONG_JMP 2

/**
 * bitmask that a signal has been received.
 *
 */

#define SIG_SIGNAL_RECEIVED 4

/**
 * default signal handler for SIGINT, SIGALRM, SIGSEGV, SIGABRT,
 * SIGFPE in the context of SAGE. It calls the default Python handler
 * if _signals.sage_handler == 0.
 */

void sage_signal_handler(int n);


/**
 * fallback signal handlers
 */
void sig_handle_sigsegv(int n);
void sig_handle_sigbus(int n);
void sig_handle_sigfpe(int n);

/**
 * Supposed to be called exactly once to ensure the SAGE signal
 * handler is the default signal handler. This function implements a
 * mechanism to make sure it is only called once so it is safe to call
 * it more than once.
 *
 */

void setup_signal_handler(void);

#if defined (__sun__) || defined (__sun)   /* Needed for Solaris below. */
   typedef void (*__sighandler_t )();
#endif

/**
 * All relevant data for the signal handler is bundled in this struct.
 *
 */

struct sage_signals {

  /**
   * Bitmask which determines what to do in the signal handler. Let
   * bit 0 be the least significant bit, then the bitmask is as follows:
   \verbatim
     bit | semantic
     -------------------------------------------
       0 | use sage signal handler (1) or default (0)
       1 | use longjmp (1) or return (0)
       2 | signal received (1) or not (0)
   \endverbatim
   *
   * You may want to use the provided definitions SIG_SAGE_HANDLER,
   * SIG_LONG_JMP, SIG_SIGNAL_RECEIVED.
   */

  int mpio;

  /**
   * An internal buffer holding where to jump after a signal has been
   * received, don't touch!
   */

  sigjmp_buf env;


  /**
   * An optional string may be passed to the signal handler which will
   * be printed. Use _sig_str for that.
   */
  char *s;


  /**
   * Unfortunately signal handler function declarations vary HUGELY
   * from platform to platform:
   */

#if defined(__CYGWIN32__)     /* Windows XP */

  _sig_func_ptr  python_handler;

#elif defined(__FreeBSD__)    /* FreeBSD */

  sig_t  python_handler;

#elif defined(__APPLE__)      /* OSX */

  sig_t  python_handler;

#elif defined (__sun__) || defined (__sun)   /* Solaris */

  __sighandler_t  python_handler;

#else                                   /* Other, e.g., Linux */

  __sighandler_t python_handler;

#endif

};

/**
 * The actual object.
 */

extern struct sage_signals _signals;

/**
 * Enables SAGE signal handling for the following C block. This macro
 * *MUST* be followed by _sig_off.
 *
 *
 * See also @ref _sig_str
 */

#define _sig_on if (_signals.mpio == 0) { _signals.mpio = 1+2; _signals.s = NULL;\
                 if (sigsetjmp(_signals.env,1)) { \
                  _signals.mpio = 0;   \
                  return(0); \
                } } // else { _signals.s = "Unbalanced _sig_on/_sig_off\n"; fprintf(stderr, _signals.s); sage_signal_handler(SIGABRT); }

/* /\** */
/*  * Enables SAGE signal handling for the following C block. This macro */
/*  * *MUST* be followed by _sig_off_short. */
/*  * */
/*  * If the following block takes very little time to compute this macro */
/*  * is the right choice. Otherwise _sig_on is appropriate. */
/*  * */
/*  * See also _sig_on and _sig_str. */
/*  * */
/*  *\/ */

/* #define _sig_on_short _signals.mpio = 1 */


/**
 * Same as @ref _sig_on but with string support
 *
 */

#define _sig_str(mstring) if (_signals.mpio == 0) { _signals.mpio = 1+2; \
                _signals.s = mstring; \
                if (sigsetjmp(_signals.env,1)) { \
                  _signals.mpio = 0; \
                 return(0); \
                } } //else { _signals.s = "Unbalanced _sig_str/_sig_off\n"; fprintf(stderr, _signals.s); sage_signal_handler(SIGABRT); }


/*
  This function may be called by the code that *raises* the signal (before raising
  the signal) to pass an error message back to the handler (for example to appear as the
  text accompanying the python RuntimeError exception)
 */
void set_sage_signal_handler_message(const char* s);

/*
  Above error message will get truncated at this length:
 */
#define SAGE_SIGNAL_HANDLER_MESSAGE_LEN 256


/**
 *
 *
 */

#define _sig_off _signals.mpio = 0;


/* /\** */
/*  * */
/*  * */
/*  *\/ */

/* #define _sig_off_short (_signals.mpio & 4) */


/**
 *
 *
 */

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
