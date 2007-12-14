/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/* $Id: cntrlc.cc,v 1.53 2007/03/13 18:54:15 Singular Exp $ */
/*
* ABSTRACT - interupt handling
*/

/* includes */
#ifdef DecAlpha_OSF1
#define _XOPEN_SOURCE_EXTENDED
#endif /* MP3-Y2 0.022UF */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <strings.h>
#include <signal.h>
#include "mod2.h"
#include "omalloc.h"
#include "tok.h"
#include "ipshell.h"
#include "febase.h"
#include "cntrlc.h"
#include "polys.h"
#include "feOpt.h"
#include "version.h"

/* undef, if you don't want GDB to come up on error */

#if !defined(__alpha)
#define CALL_GDB
#endif

#if defined(__OPTIMIZE__) && defined(CALL_GDB)
#undef CALL_GDB
#endif

#if defined(unix) && !defined(hpux)
 #include <unistd.h>
 #include <sys/types.h>

 #ifdef TIME_WITH_SYS_TIME
   #include <time.h>
   #ifdef HAVE_SYS_TIME_H
     #include <sys/time.h>
   #endif
 #else
   #ifdef HAVE_SYS_TIME_H
     #include <sys/time.h>
   #else
     #include <time.h>
   #endif
 #endif
 #ifdef HAVE_SYS_TIMES_H
   #include <sys/times.h>
 #endif

 #define INTERACTIVE 0
 #define STACK_TRACE 1

 #ifdef CALL_GDB
   static void debug (int);
   static void debug_stop (char **);
 #endif
 #ifndef __OPTIMIZE__
   static void stack_trace (char **);
   static void stack_trace_sigchld (int);
 #endif
#endif

/*---------------------------------------------------------------------*
 * File scope Variables (Variables share by several functions in
 *                       the same file )
 *
 *---------------------------------------------------------------------*/
/* data */
jmp_buf si_start_jmpbuf;
int siRandomStart;
short si_restart=0;
BOOLEAN siCntrlc = FALSE;

typedef void (*si_hdl_typ)(int);


/*0 implementation*/
/*---------------------------------------------------------------------*
 * Functions declarations
 *
 *---------------------------------------------------------------------*/
#ifndef MSDOS
/* signals are not implemented in DJGCC */
void sigint_handler(int sig);
#endif /* MSDOS */

si_hdl_typ si_set_signal ( int sig, si_hdl_typ signal_handler);

/*---------------------------------------------------------------------*/
/**
 * @brief meta function for binding a signal to an handler

 @param[in] sig             Signal number
 @param[in] signal_handler  Pointer to signal handler

 @return value of signal()
**/
/*---------------------------------------------------------------------*/
si_hdl_typ si_set_signal (
  int sig,
  si_hdl_typ signal_handler
  )
{
  si_hdl_typ retval=signal (sig, (si_hdl_typ)signal_handler);
  if (retval == SIG_ERR)
  {
     fprintf(stderr, "Unable to init signal %d ... exiting...\n", sig);
  }
#ifdef HAVE_SIGINTERRUPT
  siginterrupt(sig, 1);
#endif
  return retval;
}                               /* si_set_signal */


/*---------------------------------------------------------------------*/
#if defined(ix86_Linux)
  #if defined(HAVE_SIGCONTEXT) || defined(HAVE_ASM_SIGCONTEXT_H)
  #include <signal.h>
  #else
struct sigcontext_struct {
        unsigned short gs, __gsh;
        unsigned short fs, __fsh;
        unsigned short es, __esh;
        unsigned short ds, __dsh;
        unsigned long edi;
        unsigned long esi;
        unsigned long ebp;
        unsigned long esp;
        unsigned long ebx;
        unsigned long edx;
        unsigned long ecx;
        unsigned long eax;
        unsigned long trapno;
        unsigned long err;
        unsigned long eip;
        unsigned short cs, __csh;
        unsigned long eflags;
        unsigned long esp_at_signal;
        unsigned short ss, __ssh;
        unsigned long i387;
        unsigned long oldmask;
        unsigned long cr2;
};
#endif
#define HAVE_SIGSTRUCT
typedef struct sigcontext_struct sigcontext;
#endif

#if defined(x86_64_Linux)
#define HAVE_SIGSTRUCT
#endif


#if defined(HAVE_SIGSTRUCT)
/*2---------------------------------------------------------------------*/
/**
 * @brief signal handler for run time errors, linux/i386, x86_64 version

 @param[in] sig
 @param[in] s
**/
/*---------------------------------------------------------------------*/
void sigsegv_handler(int sig, sigcontext s)
{
  fprintf(stderr,"Singular : signal %d (v: %d/%u):\n",sig,SINGULAR_VERSION,feVersionId);
  if (sig!=SIGINT)
  {
    fprintf(stderr,"Segment fault/Bus error occurred at %lx because of %lx (r:%d)\n"
                   "please inform the authors\n",
		   #ifdef __i386__
                   (long)s.eip,
		   #else /* x86_64*/
                   (long)s.rip,
		   #endif
		   (long)s.cr2,siRandomStart);
  }
#ifdef __OPTIMIZE__
  if(si_restart<3)
  {
    si_restart++;
    fprintf(stderr,"trying to restart...\n");
    init_signals();
    longjmp(si_start_jmpbuf,1);
  }
#endif /* __OPTIMIZE__ */
#ifdef CALL_GDB
  if (sig!=SIGINT) debug(INTERACTIVE);
#endif /* CALL_GDB */
  exit(0);
}

/*---------------------------------------------------------------------*/
/**
 * @brief additional default signal handler

  // some newer Linux version cannot have SIG_IGN for SIGCHLD,
  // so use this nice routine here:
  //  SuSe 9.x reports -1 always
  //  Redhat 9.x/FC x reports sometimes -1
  // see also: hpux_system

 @param[in] sig
**/
/*---------------------------------------------------------------------*/
void sig_ign_hdl(int sig)
{
}

/*2
* init signal handlers, linux/i386 version
*/
void init_signals()
{
/*4 signal handler: linux*/
  if (SIG_ERR==si_set_signal(SIGSEGV,(si_hdl_typ)sigsegv_handler))
  {
    PrintS("cannot set signal handler for SEGV\n");
  }
  if (SIG_ERR==si_set_signal(SIGFPE, (si_hdl_typ)sigsegv_handler))
  {
    PrintS("cannot set signal handler for FPE\n");
  }
  if (SIG_ERR==si_set_signal(SIGILL, (si_hdl_typ)sigsegv_handler))
  {
    PrintS("cannot set signal handler for ILL\n");
  }
  if (SIG_ERR==si_set_signal(SIGIOT, (si_hdl_typ)sigsegv_handler))
  {
    PrintS("cannot set signal handler for IOT\n");
  }
  if (SIG_ERR==si_set_signal(SIGINT ,(si_hdl_typ)sigint_handler))
  {
    PrintS("cannot set signal handler for INT\n");
  }
  //si_set_signal(SIGCHLD, (void (*)(int))SIG_IGN);
  si_set_signal(SIGCHLD, (si_hdl_typ)sig_ign_hdl);
}

/*---------------------------------------------------------------------*/
#elif defined(SunOS) /*SPARC_SUNOS*/
/*2
* signal handler for run time errors, sparc sunos 4 version
*/
void sigsegv_handler(int sig, int code, struct sigcontext *scp, char *addr)
{
  fprintf(stderr,"Singular : signal %d, code %d (v: %d/%u):\n",
    sig,code,SINGULAR_VERSION,feVersionId);
  if ((sig!=SIGINT)&&(sig!=SIGABRT))
  {
    fprintf(stderr,"Segment fault/Bus error occurred at %x (r:%d)\n"
                   "please inform the authors\n",
                   (int)addr,siRandomStart);
  }
#ifdef __OPTIMIZE__
  if(si_restart<3)
  {
    si_restart++;
    fprintf(stderr,"trying to restart...\n");
    init_signals();
    longjmp(si_start_jmpbuf,1);
  }
#endif /* __OPTIMIZE__ */
#ifdef CALL_GDB
  if (sig!=SIGINT) debug(STACK_TRACE);
#endif /* CALL_GDB */
  exit(0);
}

/*2
* init signal handlers, sparc sunos 4 version
*/
void init_signals()
{
/*4 signal handler:*/
  si_set_signal(SIGSEGV,sigsegv_handler);
  si_set_signal(SIGBUS, sigsegv_handler);
  si_set_signal(SIGFPE, sigsegv_handler);
  si_set_signal(SIGILL, sigsegv_handler);
  si_set_signal(SIGIOT, sigsegv_handler);
  si_set_signal(SIGINT ,sigint_handler);
  si_set_signal(SIGCHLD, (void (*)(int))SIG_IGN);
}
#else

/*---------------------------------------------------------------------*/
/*2
* signal handler for run time errors, general version
*/
void sigsegv_handler(int sig)
{
  fprintf(stderr,"Singular : signal %d (v: %d/%u):\n",
    sig,SINGULAR_VERSION,feVersionId);
  if (sig!=SIGINT)
  {
    fprintf(stderr,"Segment fault/Bus error occurred (r:%d)\n"
                   "please inform the authors\n",
                   siRandomStart);
  }
  #ifdef __OPTIMIZE__
  if(si_restart<3)
  {
    si_restart++;
    fprintf(stderr,"trying to restart...\n");
    init_signals();
    longjmp(si_start_jmpbuf,1);
  }
  #endif /* __OPTIMIZE__ */
  #if defined(unix) && !defined(hpux)
  /* debug(..) does not work under HPUX (because ptrace does not work..) */
  #ifdef CALL_GDB
  #ifndef MSDOS
  if (sig!=SIGINT) debug(STACK_TRACE);
  #endif /* MSDOS */
  #endif /* CALL_GDB */
  #endif /* unix */
  exit(0);
}

/*2
* init signal handlers, general version
*/
void init_signals()
{
  #ifndef MSDOS
/* signals are not implemented in DJGCC */
/*4 signal handler:*/
  si_set_signal(SIGSEGV,(void (*) (int))sigsegv_handler);
  #ifdef SIGBUS
  si_set_signal(SIGBUS, sigsegv_handler);
  #endif /* SIGBUS */
  #ifdef SIGFPE
  si_set_signal(SIGFPE, sigsegv_handler);
  #endif /* SIGFPE */
  #ifdef SIGILL
  si_set_signal(SIGILL, sigsegv_handler);
  #endif /* SIGILL */
  #ifdef SIGIOT
  si_set_signal(SIGIOT, sigsegv_handler);
  #endif /* SIGIOT */
  #ifdef SIGXCPU
  si_set_signal(SIGXCPU, (void (*)(int))SIG_IGN);
  #endif /* SIGIOT */
  si_set_signal(SIGINT ,sigint_handler);
  #if defined(HPUX_9) || defined(HPUX_10)
  si_set_signal(SIGCHLD, (void (*)(int))SIG_IGN);
  #endif
  #endif /* !MSDOS */
}
#endif


#ifndef MSDOS
/*2
* signal handler for SIGINT
*/
void sigint_handler(int sig)
{
  mflush();
  #ifdef HAVE_FEREAD
  if (fe_is_raw_tty) fe_temp_reset();
  #endif /* HAVE_FEREAD */
  loop
  {
    int cnt=0;
    int c;
    fprintf(stderr,"// ** Interrupt at cmd:`%s` in line:'%s'\n",
      Tok2Cmdname(iiOp),my_yylinebuf);
    if (feGetOptValue(FE_OPT_EMACS) == NULL)
    {
      fputs("abort command(a), continue(c) or quit Singular(q) ?",stderr);fflush(stderr);
      c = fgetc(stdin);
    }
    else
    {
      c = 'a';
    }

    switch(c)
    {
      case 'q':
                m2_end(2);
      case 'r':
                longjmp(si_start_jmpbuf,1);
      case 'b':
                VoiceBackTrack();
                break;
      case 'a':
                siCntrlc++;
      case 'c':
                if (feGetOptValue(FE_OPT_EMACS) == NULL) fgetc(stdin);
                si_set_signal(SIGINT ,(si_hdl_typ)sigint_handler);
                return;
                //siCntrlc ++;
                //if (siCntrlc>2) si_set_signal(SIGINT,(si_hdl_typ) sigsegv_handler);
                //else            si_set_signal(SIGINT,(si_hdl_typ) sigint_handler);
    }
    cnt++;
    if(cnt>5) m2_end(2);
  }
}
#endif /* !MSDOS */

#ifndef MSDOS
//void test_int()
//{
//  if (siCntrlc!=0)
//  {
//    int saveecho = si_echo;
//    siCntrlc = FALSE;
//    si_set_signal(SIGINT ,sigint_handler);
//    iiDebug();
//    si_echo = saveecho;
//  }
//}
#endif /* !MSDOS */

#ifdef unix
# ifndef hpux
#  ifndef __OPTIMIZE__
#   ifndef MSDOS
int si_stop_stack_trace_x;
#    ifdef CALL_GDB
static void debug (int method)
{
  if (feOptValue(FE_OPT_NO_TTY))
  {
    dReportError("Caught Signal 11");
    return;
  }
  int pid;
  char buf[16];
  char *args[4] = { "gdb", "Singularg", NULL, NULL };

  #ifdef HAVE_FEREAD
  if (fe_is_raw_tty) fe_temp_reset();
  #endif /* HAVE_FEREAD */

  sprintf (buf, "%d", getpid ());

  args[2] = buf;

  pid = fork ();
  if (pid == 0)
  {
    switch (method)
    {
      case INTERACTIVE:
        fprintf (stderr, "debug_stop\n");
        debug_stop (args);
        break;
      case STACK_TRACE:
        fprintf (stderr, "stack_trace\n");
        stack_trace (args);
        break;
      default:
        // should not be reached:
        exit(1);
    }
  }
  else if (pid == -1)
  {
    perror ("could not fork");
    return;
  }

  si_stop_stack_trace_x = 1;
  while (si_stop_stack_trace_x) ;
}

static void debug_stop ( char **args)
{
  execvp (args[0], args);
  perror ("exec failed");
  _exit (0);
}
#    endif /* CALL_GDB */

static int stack_trace_done;

static void stack_trace (char **args)
{
  int pid;
  int in_fd[2];
  int out_fd[2];
  fd_set fdset;
  fd_set readset;
  struct timeval tv;
  int sel, index, state;
  char buffer[256];
  char c;

  stack_trace_done = 0;

  signal (SIGCHLD, stack_trace_sigchld);

  if ((pipe (in_fd) == -1) || (pipe (out_fd) == -1))
  {
    perror ("could open pipe");
    m2_end(999);
  }

  pid = fork ();
  if (pid == 0)
  {
    close (0); dup2 (in_fd[0],0);   /* set the stdin to the in pipe */
    close (1); dup2 (out_fd[1],1);  /* set the stdout to the out pipe */
    close (2); dup2 (out_fd[1],2);  /* set the stderr to the out pipe */

    execvp (args[0], args);      /* exec gdb */
    perror ("exec failed");
    m2_end(999);
  }
  else if (pid == -1)
  {
    perror ("could not fork");
    m2_end(999);
  }

  FD_ZERO (&fdset);
  FD_SET (out_fd[0], &fdset);

  write (in_fd[1], "backtrace\n", 10);
  write (in_fd[1], "p si_stop_stack_trace_x = 0\n", 28);
  write (in_fd[1], "quit\n", 5);

  index = 0;
  state = 0;

  loop
  {
    readset = fdset;
    tv.tv_sec = 1;
    tv.tv_usec = 0;

#    ifdef hpux
    sel = select (FD_SETSIZE, (int *)readset.fds_bits, NULL, NULL, &tv);
#    else /* hpux */
    sel = select (FD_SETSIZE, &readset, NULL, NULL, &tv);
#    endif /* hpux */
    if (sel == -1)
      break;

    if ((sel > 0) && (FD_ISSET (out_fd[0], &readset)))
    {
      if (read (out_fd[0], &c, 1))
      {
        switch (state)
        {
          case 0:
            if (c == '#')
            {
              state = 1;
              index = 0;
              buffer[index++] = c;
            }
            break;
          case 1:
            buffer[index++] = c;
            if ((c == '\n') || (c == '\r'))
            {
              buffer[index] = 0;
              fprintf (stderr, "%s", buffer);
              state = 0;
              index = 0;
            }
            break;
          default:
            break;
        }
      }
    }
    else if (stack_trace_done)
      break;
  }

  close (in_fd[0]);
  close (in_fd[1]);
  close (out_fd[0]);
  close (out_fd[1]);
  m2_end(0);
}

static void stack_trace_sigchld (int signum)
{
  stack_trace_done = 1;
}

#   endif /* !MSDOS */
#  endif /* !__OPTIMIZE__ */
# endif /* !hpux */
#endif /* unix */

/* Under HPUX 9, system(...) returns -1 if SIGCHLD does not equal
   SIG_DFL. However, if it stays at SIG_DFL we get zombie processes
   for terminated childs generated by fork. Therefors some special treatment
   is necessary */
#ifdef HPUX_9
# undef system
extern "C" {
  int  hpux9_system(const char* call)
  {
    int ret;
    si_set_signal(SIGCHLD, (void (*)(int))SIG_DFL);
    ret = system(call);
    si_set_signal(SIGCHLD, (void (*)(int))SIG_IGN);
    return ret;
  }
}
#endif /* HPUX_9 */
