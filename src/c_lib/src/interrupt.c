/******************************************************************************
       Copyright (C) 2006 William Stein <wstein@gmail.com>
                     2006 Martin Albrecht <malb@informatik.uni-bremen.de>

  Distributed under the terms of the GNU General Public License (GPL), Version 2.

  The full text of the GPL is available at:
                  http://www.gnu.org/licenses/

******************************************************************************/

#include "stdsage.h"
#include "interrupt.h"
#include <stdio.h>

struct sage_signals _signals;

void msg(char* s);

void sage_signal_handler(int sig) {

  char *s = _signals.s;
  _signals.s = NULL;

  //we override the default handler
  if ( _signals.mpio & 1 ) {

    //what to do?

    switch(sig) {

    case SIGINT:
      if( s ) {
	PyErr_SetString(PyExc_KeyboardInterrupt, s);
      } else {
	PyErr_SetString(PyExc_KeyboardInterrupt, "");
      }
      break;

    case SIGALRM:
      if( s ) {
	PyErr_SetString(PyExc_KeyboardInterrupt, s);
      } else {
	PyErr_SetString(PyExc_KeyboardInterrupt, "Alarm received");
      }
      break;

    default:
      if( s ) {
	PyErr_SetString(PyExc_RuntimeError, s);
      } else {
	PyErr_SetString(PyExc_RuntimeError, "");
      }
    }


    //notify 'calling' function

    _signals.mpio |= 4;

    signal(sig, sage_signal_handler);
    //where to go next?
    if ( _signals.mpio & 2 ) {
      siglongjmp(_signals.env, sig);
    } else {
      //this case shouldn't happen as _sig_[on|off]_short is disabled
      return;
    }

  } else {

    //we use the default handler
    _signals.mpio = 0;

    switch(sig) {
    case SIGSEGV:
      sig_handle_sigsegv(sig);
      break;
    case SIGBUS:
      sig_handle_sigbus(sig);
      break;
    case SIGFPE:
      sig_handle_sigfpe(sig);
      break;
    default:
      _signals.python_handler(sig);
      break;
    };

    signal(sig, sage_signal_handler);
  }
}

void setup_signal_handler(void) {
  void *tmp = NULL;

  //we need to make sure not to store our own signal handler as the old one
  //because that would cause infinite recursion if an error occurs and we are
  //not supposed to catch it.
  tmp = signal(SIGINT, sage_signal_handler);

  if(sage_signal_handler != tmp) {
    _signals.python_handler = tmp;
  }

  _signals.s = NULL;

/*  signal(SIGBUS, sig_handle_sigbus); */
  signal(SIGBUS, sage_signal_handler);

  signal(SIGALRM,sage_signal_handler);
  signal(SIGSEGV,sage_signal_handler);
  signal(SIGABRT,sage_signal_handler);
  signal(SIGFPE,sage_signal_handler);
}


void msg(char* s) {
  fprintf(stderr, "\n\n------------------------------------------------------------\n");
  fprintf(stderr, s);
  fprintf(stderr, "This probably occured because a *compiled* component\n");
  fprintf(stderr, "of SAGE has a bug in it (typically accessing invalid memory)\n");
  fprintf(stderr, "or is not properly wrapped with _sig_on, _sig_off.\n");
  fprintf(stderr, "You might want to run SAGE under gdb with 'sage -gdb' to debug this.\n");
  fprintf(stderr, "SAGE will now terminate (sorry).\n");
  fprintf(stderr, "------------------------------------------------------------\n\n");
}

void sig_handle_sigsegv(int n) {
  msg("Unhandled SIGSEGV: A segmentation fault occured in SAGE.\n");
  PyErr_SetString(PyExc_KeyboardInterrupt, "");
  exit(1);
}

void sig_handle_sigbus(int n) {
  msg("Unhandled SIGBUS: A bus error occured in SAGE.\n");
  PyErr_SetString(PyExc_KeyboardInterrupt, "");
  exit(1);
}

void sig_handle_sigfpe(int n) {
  msg("Unhandled SIGFPE: An unhandled floating point exception occured in SAGE.\n");
  PyErr_SetString(PyExc_KeyboardInterrupt, "");
  exit(1);
}
