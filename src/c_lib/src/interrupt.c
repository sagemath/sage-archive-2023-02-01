/*
Interrupt and signal handling for Sage

AUTHORS:

- William Stein, Martin Albrecht (2006): initial version

- Jeroen Demeyer (2010-10-03): almost complete rewrite (#9678)

*/

/*****************************************************************************
 *       Copyright (C) 2006 William Stein <wstein@gmail.com>
 *                     2006 Martin Albrecht <malb@informatik.uni-bremen.de>
 *                     2010-2013 Jeroen Demeyer <jdemeyer@cage.ugent.be>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <limits.h>
/* glibc has a backtrace() command since version 2.1 */
#ifdef __GLIBC__
#if (__GLIBC__ > 2) || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 1)
#define HAVE_BACKTRACE 1
#include <execinfo.h>
#endif
#endif
#include "stdsage.h"
#include "interrupt.h"


struct sage_signals_t _signals;

/* The default signal mask during normal operation,
 * initialized by setup_sage_signal_handler(). */
static sigset_t default_sigmask;

/* default_sigmask with SIGINT and SIGALRM added. */
static sigset_t sigmask_with_sigint;

/* Does this processor support the x86 EMMS instruction? */
#if defined(__i386__) || defined(__x86_64__)
#define CPU_ARCH_x86
static int cpu_has_emms = 0;
#endif

/* Do whatever is needed to reset the CPU to a sane state after
 * handling a signal.  In particular on x86 CPUs, we need to clear
 * the FPU (this is needed after MMX instructions have been used or
 * if an interrupt occurs during an FPU computation).
 * Linux and OS X 10.6 do this as part of their signals implementation,
 * but Solaris doesn't.  Since this code is called only when handling a
 * signal (which should be very rare), it's better to play safe and
 * always execute this instead of special-casing based on the operating
 * system.
 * See http://trac.sagemath.org/sage_trac/ticket/12873
 */
static inline void reset_CPU()
{
#ifdef CPU_ARCH_x86
    /* Clear FPU tag word */
    if (cpu_has_emms)
    {
        asm("emms");
    }
    else
    {
        asm("ffree %st(0)");
        asm("ffree %st(1)");
        asm("ffree %st(2)");
        asm("ffree %st(3)");
        asm("ffree %st(4)");
        asm("ffree %st(5)");
        asm("ffree %st(6)");
        asm("ffree %st(7)");
    }
#endif
}


/* Handler for SIGINT */
void sage_interrupt_handler(int sig)
{
#if ENABLE_DEBUG_INTERRUPT
    fprintf(stderr, "\n*** SIGINT *** %s sig_on\n", (_signals.sig_on_count > 0) ? "inside" : "outside");
    print_backtrace();
#endif

    if (_signals.sig_on_count > 0)
    {
        if (!_signals.block_sigint)
        {
            /* Actually raise an exception so Python can see it */
            sig_raise_exception(sig);

            /* Jump back to sig_on() (the first one if there is a stack) */
            reset_CPU();
            siglongjmp(_signals.env, sig);
        }
    }
    else
    {
        /* Set an internal Python flag that an interrupt has been
         * raised.  This will not immediately raise an exception, only
         * on the next call of PyErr_CheckSignals().  We cannot simply
         * raise an exception here because of Python's "global
         * interpreter lock" -- Jeroen Demeyer */
        PyErr_SetInterrupt();
    }

    /* If we are here, we cannot handle the interrupt immediately, so
     * we set interrupt_received for later use. */
    _signals.interrupt_received = 1;
}

/* Handler for SIGILL, SIGABRT, SIGFPE, SIGBUS, SIGSEGV */
void sage_signal_handler(int sig)
{
    sig_atomic_t inside = _signals.inside_signal_handler;
    _signals.inside_signal_handler = 1;

    if (inside == 0 && _signals.sig_on_count > 0)
    {
        /* We are inside sig_on(), so we can handle the signal! */
#if ENABLE_DEBUG_INTERRUPT
        fprintf(stderr, "\n*** SIG %i *** inside sig_on\n", sig);
        print_backtrace();
#endif

        /* Actually raise an exception so Python can see it */
        sig_raise_exception(sig);

        /* Jump back to sig_on() (the first one if there is a stack) */
        reset_CPU();
        siglongjmp(_signals.env, sig);
    }
    else
    {
        /* We are outside sig_on() and have no choice but to terminate Sage */

        /* Reset all signals to their default behaviour and unblock
         * them in case something goes wrong as of now. */
        signal(SIGILL, SIG_DFL);
        signal(SIGABRT, SIG_DFL);
        signal(SIGFPE, SIG_DFL);
        signal(SIGBUS, SIG_DFL);
        signal(SIGSEGV, SIG_DFL);
        sigprocmask(SIG_SETMASK, &sigmask_with_sigint, NULL);

        if (inside) sigdie(sig, "An error occured during signal handling.");

        /* Quit Sage with an appropriate message. */
        switch(sig)
        {
            case SIGILL:
                sigdie(sig, "Unhandled SIGILL: An illegal instruction occurred in Sage.");
                break;  /* This will not be reached */
            case SIGABRT:
                sigdie(sig, "Unhandled SIGABRT: An abort() occurred in Sage.");
                break;  /* This will not be reached */
            case SIGFPE:
                sigdie(sig, "Unhandled SIGFPE: An unhandled floating point exception occurred in Sage.");
                break;  /* This will not be reached */
            case SIGBUS:
                sigdie(sig, "Unhandled SIGBUS: A bus error occurred in Sage.");
                break;  /* This will not be reached */
            case SIGSEGV:
                sigdie(sig, "Unhandled SIGSEGV: A segmentation fault occurred in Sage.");
                break;  /* This will not be reached */
        };
        sigdie(sig, "Unknown signal received.\n");
    }
}


void sig_raise_exception(int sig)
{
#if ENABLE_DEBUG_INTERRUPT
    fprintf(stderr, "sig_raise_exception(sig=%i)\nPyErr_Occurred() = %p\nRaising Python exception...\n",
        sig, PyErr_Occurred());
    fflush(stderr);
#endif

    /* String to be printed in the Python exception */
    const char* msg = _signals.s;

    switch(sig)
    {
        case SIGINT:
            PyErr_SetNone(PyExc_KeyboardInterrupt);
            break;
        case SIGILL:
            if (!msg) msg = "Illegal instruction";
            PyErr_SetString(PyExc_RuntimeError, msg);
            break;
        case SIGABRT:
            if (!msg) msg = "Aborted";
            PyErr_SetString(PyExc_RuntimeError, msg);
            break;
        case SIGFPE:
            if (!msg) msg = "Floating point exception";
            PyErr_SetString(PyExc_RuntimeError, msg);
            break;
        case SIGBUS:
            if (!msg) msg = "Bus error";
            PyErr_SetString(PyExc_RuntimeError, msg);
            break;
        case SIGSEGV:
            if (!msg) msg = "Segmentation fault";
            PyErr_SetString(PyExc_RuntimeError, msg);
            break;
        default:
            PyErr_SetString(PyExc_SystemError, "Unknown signal");
    }
}


/* Handle an interrupt before sig_on(). */
int _sig_on_interrupt_received()
{
    _signals.interrupt_received = 0;
    if (PyErr_CheckSignals())
    {
        _signals.sig_on_count = 0;
        return 0;
    }
    return 1;
}

/* Recover after siglongjmp() */
void _sig_on_recover()
{
    _signals.block_sigint = 0;
    _signals.sig_on_count = 0;

    /* Reset signal mask */
    sigprocmask(SIG_SETMASK, &default_sigmask, NULL);
    _signals.inside_signal_handler = 0;
}

void _sig_off_warning(const char* file, int line)
{
    char buf[320];
    snprintf(buf, sizeof(buf), "sig_off() without sig_on() at %s:%i", file, line);
    PyErr_WarnEx(PyExc_RuntimeWarning, buf, 2);
    print_backtrace();
}


void set_sage_signal_handler_message(const char* s)
{
    _signals.s = s;
}


void setup_sage_signal_handler()
{
    /* Reset the _signals structure */
    memset(&_signals, 0, sizeof(_signals));

    /* Save the default signal mask */
    sigprocmask(SIG_BLOCK, NULL, &default_sigmask);

    /* Save the signal mask with SIGINT and SIGALRM */
    sigprocmask(SIG_BLOCK, NULL, &sigmask_with_sigint);
    sigaddset(&sigmask_with_sigint, SIGINT);
    sigaddset(&sigmask_with_sigint, SIGALRM);

    /* Install signal handlers */
    struct sigaction sa;
    memset(&sa, 0, sizeof(sa));
    /* Block SIGINT and SIGALRM during the signal handlers */
    sigemptyset(&sa.sa_mask);
    sigaddset(&sa.sa_mask, SIGINT);
    sigaddset(&sa.sa_mask, SIGALRM);

    sa.sa_handler = sage_interrupt_handler;
    if (sigaction(SIGINT, &sa, NULL)) {perror("sigaction"); exit(1);}
    sa.sa_handler = sage_signal_handler;
    /* Allow signals during signal handling, we have code to deal with
     * this case. */
    sa.sa_flags |= SA_NODEFER;
    if (sigaction(SIGILL, &sa, NULL)) {perror("sigaction"); exit(1);}
    if (sigaction(SIGABRT, &sa, NULL)) {perror("sigaction"); exit(1);}
    if (sigaction(SIGFPE, &sa, NULL)) {perror("sigaction"); exit(1);}
    if (sigaction(SIGBUS, &sa, NULL)) {perror("sigaction"); exit(1);}
    if (sigaction(SIGSEGV, &sa, NULL)) {perror("sigaction"); exit(1);}

    /* If the CPU architecture is x86, check whether the EMMS
     * instruction is supported by executing it and catching a
     * possible SIGILL (illegal instruction signal). */
#ifdef CPU_ARCH_x86
    if (!cpu_has_emms)
    {
        if (sig_on_no_except())  /* try: */
        {
            asm("emms");
            sig_off();
            cpu_has_emms = 1;
        }
        else  /* except: */
        {
            PyErr_Clear();  /* Clear Python exception */
        }
    }
#endif
}


void print_backtrace()
{
    void* backtracebuffer[1024];
    fflush(stderr);
#ifdef HAVE_BACKTRACE
    int btsize = backtrace(backtracebuffer, 1024);
    backtrace_symbols_fd(backtracebuffer, btsize, 2);
#endif
}

void print_enhanced_backtrace()
{
    /* Flush all buffers before forking */
    fflush(stdout);
    fflush(stderr);

    pid_t parent_pid = getpid();
    pid_t pid = fork();

    if (pid < 0)
    {
        /* Failed to fork: no problem, just ignore */
        perror("fork");
        return;
    }

    if (pid == 0) { /* child */
        /* We deliberately put these variables on the stack to avoid
         * malloc() calls, the heap might be messed up! */
        char path[1024];
        char pid_str[32];
        char* argv[5];

        snprintf(path, sizeof(path), "%s/bin/sage-CSI", getenv("SAGE_LOCAL"));
        snprintf(pid_str, sizeof(pid_str), "%i", parent_pid);

        argv[0] = "sage-CSI";
        argv[1] = "--no-color";
        argv[2] = "--pid";
        argv[3] = pid_str;
        argv[4] = NULL;
        execv(path, argv);
        perror("Failed to execute sage-CSI");
        exit(2);
    }
    /* Wait for sage-CSI to finish */
    waitpid(pid, NULL, 0);
}


void sigdie(int sig, const char* s)
{
    fprintf(stderr,
        "------------------------------------------------------------------------\n");
    print_backtrace();

    fprintf(stderr,
        "------------------------------------------------------------------------\n");
    print_enhanced_backtrace();

    fprintf(stderr,
        "------------------------------------------------------------------------\n"
        "%s\n"
        "This probably occurred because a *compiled* component of Sage has a bug\n"
        "in it and is not properly wrapped with sig_on(), sig_off(). You might\n"
        "want to run Sage under gdb with 'sage -gdb' to debug this.\n"
        "Sage will now terminate.\n"
        "------------------------------------------------------------------------\n",
        s);
    fflush(stderr);

    /* Suicide with signal ``sig`` */
    kill(getpid(), sig);

    /* We should be dead! */
    exit(128 + sig);
}
