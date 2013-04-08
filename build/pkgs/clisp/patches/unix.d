/*
 * The include file for the UNIX version of CLISP
 * Bruno Haible 1990-2006
 * Sam Steingold 1998-2005
 */

/* control character constants: */
#define BEL  7              /* ring the bell */
/* #define NL  10              new line, see <lispbibl.d> */
#define RUBOUT 127          /* Rubout = Delete */
#define CRLFstring  "\n"    /* C string containing Newline */

#define stdin_handle   0  /* the file handle for the standard input */
#define stdout_handle  1  /* the file handle for the standard output */
#define stderr_handle  2  /* the file handle for the standard error */

/* Declaration of types of I/O parameters of operating system functions */
#ifdef STDC_HEADERS
  #include <stdlib.h>
#endif
#include <sys/types.h>  /* declares pid_t, uid_t */
#if defined(TIME_WITH_SYS_TIME)
 #include <sys/time.h>
 #include <time.h>
#else
 #if defined(HAVE_SYS_TIME_H)
  #include <sys/time.h>
 #elif defined(HAVE_TIME_H)
  #include <time.h>
 #endif
#endif
#ifdef HAVE_UNISTD_H
  #include <unistd.h>
#endif
#if defined(HAVE_SYS_RESOURCE_H)
 #include <sys/resource.h>
#endif

/* the table of the system error messages */
#include <errno.h>
/* extern int errno; */ /* last error code */
/* NB: errno may be a macro which expands to a function call.
   Therefore access and assignment to errno must be wrapped in
   begin_system_call()/end_system_call() */
#define OS_errno errno
#define OS_set_errno(e) (errno=(e))
#ifdef HAVE_STRERROR
#include <string.h>
extern_C char* strerror (int errnum);
#else
/* sys_errlist[] and sys_nerr are defined in either <stdio.h> or <errno.h> */
#endif
/* perror(3)
   On UnixWare 7.0.1 some errno value is defined to an invalid negative value,
   causing an out-of-bounds array access in errunix.d. */
#if (EDQUOT < 0)
  #undef EDQUOT
#endif
/* used by ERROR, SPVW, STREAM, PATHNAME */

/* Make the main memory available */
#ifdef HAVE_GETPAGESIZE
  extern_C RETGETPAGESIZETYPE getpagesize (void); /* getpagesize(2) */
#endif
/* malloc(), free(), realloc() are defined in <stdlib.h> */
#ifdef UNIX_NEXTSTEP
  /* ignore the contents of libposix.a, since it is not documented */
  #undef HAVE_MMAP
  #undef HAVE_MUNMAP
#endif
#ifdef UNIX_RHAPSODY
/* Ignore mmap and friends, because the configure test says no working mmap. */
  #undef HAVE_MMAP
  #undef HAVE_MUNMAP
  #undef HAVE_WORKING_MPROTECT
#endif
#if defined(HAVE_MMAP) || defined(HAVE_MMAP_ANON) || defined(HAVE_MMAP_ANONYMOUS) || defined(HAVE_MMAP_DEVZERO) || defined(HAVE_MMAP_DEVZERO_SUN4_29)
  #include <sys/mman.h>
  #if defined(HAVE_MMAP_ANONYMOUS) && !defined(HAVE_MMAP_ANON)
    /* HP-UX uses MAP_ANONYMOUS instead of MAP_ANON. */
    #define MAP_ANON MAP_ANONYMOUS
    #define HAVE_MMAP_ANON
  #endif
  #if defined(UNIX_SUNOS4) || defined(UNIX_SUNOS5)
    /* for SINGLEMAP_MEMORY: */
    #if defined(HAVE_MMAP_DEVZERO_SUN4_29) && defined(SUN4_29) && !defined(HAVE_MMAP_DEVZERO)
      /* On the assumption of the SUN4_29-type code distribution
         HAVE_MMAP_DEVZERO_SUN4_29 is a sufficient replacement
         for HAVE_MMAP_DEVZERO. */
      #define HAVE_MMAP_DEVZERO
    #endif
  #endif
  #ifdef UNIX_SUNOS5
   /* NB: Under UNIX_SUNOS5, HAVE_MMAP_DEVZERO should be defined.
      There is however a limit of 25 MB mmap() memory.
      Since the shared memory facility of UNIX_SUNOS5 denies
      memory at addresses >= 0x06000000 or more than 6 times to attach,
      we must use SINGLEMAP_MEMORY */
  #endif
  #ifdef HAVE_MSYNC
    #ifdef MS_INVALIDATE
      /* tested only on UNIX_LINUX, not UNIX_SUNOS4, not UNIX_SUNOS5,
         not UNIX_FREEBSD. ?? */
      /* for MULTIMAP_MEMORY_VIA_FILE: */
      /* extern_C int msync (void* addr, size_t len, int flags); */
    #else
      /* NetBSD has a 2-argument msync(), unusable for our purposes. */
      #undef HAVE_MSYNC
    #endif
  #endif
  /* for MULTIMAP_MEMORY_VIA_FILE: */
  #if defined(UNIX_SUNOS4) || defined(UNIX_LINUX) || defined(UNIX_FREEBSD)
    #if HAVE_SYS_STATVFS_H
      #include <sys/statvfs.h>
    #elif HAVE_SYS_STATFS_H
      #include <sys/statfs.h>
    #else
      /* Old BSDs define struct statfs in <sys/mount.h>. */
      #include <sys/param.h>
      #include <sys/mount.h>
    #endif
  #endif
#endif
#ifdef HAVE_MACH_VM /* vm_allocate(), task_self(), ... available */
  /* the headers for UNIX_NEXTSTEP must look indescribable ... */
  #undef local
  #include <mach/mach_interface.h>
  #if defined(UNIX_NEXTSTEP) || defined(UNIX_RHAPSODY)
    #include <mach/mach_init.h>
  #endif
  #ifdef UNIX_OSF
    #include <mach_init.h>
  #endif
  /* #include <mach/mach.h> */
  #include <mach/mach_traps.h> /* for map_fd() */
  #include <mach/machine/vm_param.h>
  #define local static
  /* thus one can use mmap(), munmap() und mprotect(). see spvw.d. */
  #define HAVE_MMAP
  #define HAVE_MUNMAP
  #define HAVE_WORKING_MPROTECT
  /* We assume the following types are effectively the same:
       vm_address_t and void*
       vm_size_t and size_t
   */
  #define PROT_NONE  0
  #define PROT_READ  VM_PROT_READ
  #define PROT_WRITE VM_PROT_WRITE
  #define PROT_EXEC  VM_PROT_EXECUTE
#endif
#ifdef HAVE_MMAP
  /* extern_C void* mmap (void* addr, size_t len, int prot, int flags, int fd, off_t off); */ /* MMAP(2) */
#endif
#ifdef HAVE_MUNMAP
  /* extern_C int munmap (void* addr, size_t len); */ /* MUNMAP(2) */
#endif
#ifdef HAVE_WORKING_MPROTECT
  /* extern_C int mprotect (void* addr, size_t len, int prot); */ /* MPROTECT(2) */
#endif
/* Possible values of prot: PROT_NONE, PROT_READ, PROT_READ_WRITE. */
#ifndef PROT_NONE
  #define PROT_NONE  0
#endif
#define PROT_READ_WRITE  (PROT_READ | PROT_WRITE)
#ifdef HAVE_SHM
  #include <sys/ipc.h>
  #include <sys/shm.h>
  #ifdef HAVE_SYS_SYSMACROS_H
    #include <sys/sysmacros.h>
  #endif
  #ifdef UNIX_HPUX
    #include <sys/vmmac.h> /* for SHMLBA */
  #endif
  #ifdef UNIX_AUX
    #include <sys/mmu.h> /* for SHMLBA */
  #endif
  #if defined(UNIX_LINUX) && !defined(UNIX_GNU)
    #include <asm/page.h> /* for SHMLBA on Linux 2.0 */
  #endif
  #if defined(UNIX_SUNOS4) || defined(UNIX_SUNOS5)
    #define SHMMAX  0x100000 /* maximum shared memory segment size = 1 MB */
  #endif
  #ifndef SHMMAX
    #define SHMMAX  0xFFFFFFFFUL /* maximum shared memory segment size accepted to mean infinite */
  #endif
/* <sys/shm.h> declares shmget(), shmat(), shmdt(), shmctl() */
#endif
/* used by SPVW, STREAM */

/* paging control */
#ifdef HAVE_VADVISE
  #include <sys/vadvise.h> /* control codes */
  extern_C void vadvise (int param); /* paging system control, see VADVISE(2) */
#endif
/* use madvise() ?? */
/* used by SPVW */

/* make stack large enough */
#ifdef UNIX_NEXTSTEP
  extern_C int getrlimit (RLIMIT_RESOURCE_T resource, struct rlimit * rlim); /* GETRLIMIT(2) */
  extern_C int setrlimit (RLIMIT_RESOURCE_T resource, SETRLIMIT_CONST struct rlimit * rlim); /* SETRLIMIT(2) */
#endif
/* used by SPVW */

/* normal program end */
nonreturning_function(extern_C, _exit, (int status)); /* EXIT(2V) */
nonreturning_function(extern_C, exit, (int status)); /* EXIT(2V) */
/* used by SPVW, PATHNAME, STREAM */

/* Immediate abnormal termination, jump into the debugger */
extern_C ABORT_VOLATILE RETABORTTYPE abort (void); /* ABORT(3) */
/* used by SPVW, DEBUG, EVAL, IO */

/* signal handling */
#include <signal.h>
/* a signal handler is a non-returning function. */
#ifdef __cplusplus
  #ifdef SIGTYPE_DOTS
    typedef RETSIGTYPE (*signal_handler_t) (...);
  #else
    typedef RETSIGTYPE (*signal_handler_t) (int);
  #endif
#else
  typedef RETSIGTYPE (*signal_handler_t) ();
#endif
/* install a signal cleanly: */
extern_C signal_handler_t signal (int sig, signal_handler_t handler); /* SIGNAL(3V) */
#if defined(SIGNAL_NEED_UNBLOCK_OTHERS) && defined(HAVE_SIGACTION)
/* On some BSD systems (e.g. SunOS 4.1.3_U1), the call of a signal handler
   is different when the current signal is blocked.
   We therefore use sigaction() instead of signal(). */
  #define USE_SIGACTION
#endif
extern signal_handler_t install_signal_handler (int sig, signal_handler_t handler);
#define SIGNAL(sig,handler)  install_signal_handler(sig,handler)
/* a signal block and release: */
#if defined(SIGNALBLOCK_POSIX)
  /* extern_C int sigprocmask (int how, const sigset_t* set, sigset_t* oset); */ /* SIGPROCMASK(2V) */
  /* extern_C int sigemptyset (sigset_t* set); */ /* SIGSETOPS(3V) */
  /* extern_C int sigaddset (sigset_t* set, int signo); */ /* SIGSETOPS(3V) */
  #define signalblock_on(sig)  \
      { var sigset_t sigblock_mask;                                 \
        sigemptyset(&sigblock_mask); sigaddset(&sigblock_mask,sig); \
        sigprocmask(SIG_BLOCK,&sigblock_mask,NULL);
  #define signalblock_off(sig)  \
        sigprocmask(SIG_UNBLOCK,&sigblock_mask,NULL); \
      }
#elif defined(SIGNALBLOCK_SYSV)
  extern_C int sighold (int sig);
  extern_C int sigrelse (int sig);
  #define signalblock_on(sig)  sighold(sig);
  #define signalblock_off(sig)  sigrelse(sig);
#elif defined(SIGNALBLOCK_BSD)
  extern_C int sigblock (int mask); /* SIGBLOCK(2) */
  extern_C int sigsetmask (int mask); /* SIGSETMASK(2) */
  #define signalblock_on(sig)  \
      { var int old_sigblock_mask = sigblock(sigmask(sig));
  #define signalblock_off(sig)  \
        sigsetmask(old_sigblock_mask); \
      }
#else
  #error "How does one block a signal?"
#endif
/* deliver a signal some time later: */
/* extern_C {unsigned|} int alarm ({unsigned|} int seconds); / * ALARM(3V) */
#if !defined(HAVE_UALARM) && defined(HAVE_SETITIMER)
  #define NEED_OWN_UALARM /* ualarm() can be implemented with setitimer() */
  /* declares setitimer() */
  #define HAVE_UALARM
#endif
#ifdef HAVE_UALARM
  #ifdef UNIX_CYGWIN32
    /* <sys/types.h>: typedef long useconds_t; */
    extern_C useconds_t ualarm (useconds_t value, useconds_t interval);
  #else
    extern_C unsigned int ualarm (unsigned int value, unsigned int interval);
  #endif
#endif
/* acknowledge the arrival of a signal (from the signal handler): */
#ifdef USE_SIGACTION
  #ifdef SIGACTION_NEED_REINSTALL
    /* restore the handler */
    #define signal_acknowledge(sig,handler) install_signal_handler(sig,handler)
  #else /* BSD-stype signals do not need this */
    #define signal_acknowledge(sig,handler)
  #endif
#else
  #ifdef SIGNAL_NEED_REINSTALL /* UNIX_SYSV || UNIX_LINUX || ... */
    /* restore the handler */
    #define signal_acknowledge(sig,handler) install_signal_handler(sig,handler)
  #else  /* BSD-stype signals do not need this */
    #define signal_acknowledge(sig,handler)
  #endif
#endif
/* the signal one gets on termination of the child process: SIGCLD */
#if defined(SIGCHLD) && !defined(SIGCLD)
  #define SIGCLD  SIGCHLD
#endif
/* the behavior of the signals the affect system calls:
   flag=0: after the signal SIG the system call keep running.
   flag=1: after the signal SIG the system call is aborted, errno=EINTR. */
#ifdef EINTR
  extern_C int siginterrupt (int sig, int flag); /* SIGINTERRUPT(3V) */
  #ifndef HAVE_SIGINTERRUPT
    /* siginterrupt() can be implemented with sigaction() or sigvec() */
    #define NEED_OWN_SIGINTERRUPT
  #endif
#else
  #define siginterrupt(sig,flag)
#endif
/* For recovery from the SIGSEGV signal (write attempts to write
   protected ranges). See libsigsegv.
   Watch out: Hans-J. Boehm <boehm@parc.xerox.com> says that write accesses
   originating from OS calls (e.g. read()) do not trigger a signal on many
   systems, unexpectedly. (It works on Linux, though) */
#ifndef SPVW_MIXED_BLOCKS
/* We are lucky to write with read() only into the C-stack and into strings
   and not into possibly mprotect-protected ranges. */
#endif
/* raise a signal. */
#ifdef HAVE_RAISE
extern_C int raise (int sig);
#endif
/* used by SPVW */

/* check environment variables:
   getenv(), putenv(), setenv() are declared in <stdlib.h> */

/* Adjustment to locale preferences: */
#include <locale.h>
/* declares setlocale() */
/* used by SPVW */

/* get user home directory: */
#include <pwd.h>
/* declares getpwnam(), getpwuid(), getuid(), getlogin() */
/* used by PATHNAME, SPVW */

/* set working directory: */
/* chdir() - declared in <unistd.h> */
/* used by PATHNAME */

/* get working directory: */
#include <sys/param.h>
/* maximum path length (incl. terminating NULL), returned by getwd(): */
#ifndef MAXPATHLEN
  #define MAXPATHLEN  1024  /* <sys/param.h> */
#endif
#ifdef HAVE_GETCWD
extern_C char* getcwd (char* buf, GETCWD_SIZE_T bufsize);
#define getwd(buf)  getcwd(buf,MAXPATHLEN)
#else
extern_C char* getwd (char* pathname); /* GETWD(3) */
#endif
/* used by PATHNAME */

/* maximum number of symbolic links which are successively resolved: */
#ifndef MAXSYMLINKS
  #define MAXSYMLINKS  8  /* <sys/param.h> */
#endif
/* used by PATHNAME */

/* resolve symbolic links in pathname: */
#ifdef HAVE_READLINK
/* extern_C ssize_t readlink (const char* path, char* buf, size_t bufsiz); */ /* READLINK(2) */
#endif
/* used by PATHNAME */

/* get information about a file: */
#include <sys/stat.h>
#ifdef STAT_MACROS_BROKEN
  #undef S_ISDIR
  #undef S_ISLNK
  #undef S_ISREG
#endif
/* extern/extern_C int stat (const char* path, struct stat * buf); */ /* STAT(2V) */
#ifdef HAVE_LSTAT
  /* extern/extern_C int lstat (const char* path, struct stat * buf); */ /* STAT(2V) */
#else
  #define lstat stat
  #define S_ISLNK(m)  false
#endif
/* extern/extern_C int fstat (int fd, struct stat * buf); */ /* STAT(2V) */
#ifndef S_ISDIR
  #define S_ISDIR(m)  (((m)&S_IFMT) == S_IFDIR)
#endif
#ifndef S_ISLNK
  #define S_ISLNK(m)  (((m)&S_IFMT) == S_IFLNK)
#endif
#ifndef S_ISREG
  #define S_ISREG(m)  (((m)&S_IFMT) == S_IFREG)
#endif
/* used by PATHNAME, STREAM, SPVW */
struct file_id {              /* Unique ID for a file on this machine */
#if defined(SIZEOF_DEV_T)
  dev_t device;
#else
  uintL device;
#endif
#if defined(SIZEOF_INO_T)
  ino_t inode;
#else
  uintL inode;
#endif
};
/* if file NAMESTRING exists, fill file_id and call function on it,
   otherwise return NULL */
extern void* with_file_id (char *namestring, void *data,
                           void* (*func) (struct file_id *fid, void*data));
/* fill FI for an existing file handle */
typedef int errno_t;
extern errno_t handle_file_id (int fd, struct file_id *fi);
/* if the file IDs are identical, return 1, otherwise return 0 */
extern int file_id_eq (struct file_id *fi1, struct file_id *fi2);

/* unlink() - declared in <unistd.h>; used by PATHNAME, UNIXAUX */
/* rename() - declared in <stdio.h>; used by PATHNAME, UNIXAUX */

/* directory search: */
#if defined(DIRENT) || defined(_POSIX_VERSION)
  #include <dirent.h>
  #define SDIRENT  struct dirent
#else
  #ifdef SYSNDIR
    #include <sys/ndir.h>
  #else
    #ifdef SYSDIR
      #include <sys/dir.h>
    #else
      #ifdef NDIR
        #include <ndir.h>
      #else
        #include <dir.h>
      #endif
    #endif
  #endif
  #define SDIRENT  struct direct
#endif
/* declared in one of the above includes: opendir(), readdir(), closedir() */
#ifdef VOID_CLOSEDIR
  #define CLOSEDIR(dirp)  (closedir(dirp),0)
#else
  #define CLOSEDIR  closedir
#endif
/* used by PATHNAME */

/* create directory: */
/* mkdir() - declared in <sys/stat.h> or <unistd.h> */
/* used by PATHNAME */

/* remove directory: */
/* rmdir() - declared in <unistd.h> */
/* used by PATHNAME */

/* work with open files: */
#include <fcntl.h>
/* declares open() */
#if defined(ACCESS_NEEDS_SYS_FILE_H) || defined(OPEN_NEEDS_SYS_FILE_H)
  #include <sys/file.h>
#endif
/* Only a few Unices (like UNIX_CYGWIN32) have O_TEXT and O_BINARY.
   BeOS 5 has them, but they have no effect. */
#ifdef UNIX_BEOS
  #undef O_BINARY
#endif
#ifndef O_BINARY
  #define O_BINARY  0
#endif
#define my_open_mask  0644
#define Handle  int  /* the type of a file descriptor */
#define INVALID_HANDLE  -1
extern_C off_t lseek (int fd, off_t offset, int whence); /* LSEEK(2V) */
#ifndef SEEK_SET /* e.g., UNIX_NEXTSTEP */
  /* position modes, see <unistd.h> : */
  #define SEEK_SET  0
  #define SEEK_CUR  1
  #define SEEK_END  2
#endif
/* extern_C ssize_t read (int fd, void* buf, size_t nbyte); */ /* READ(2V) */
/* extern_C ssize_t write (int fd, const void* buf, size_t nbyte); */ /* WRITE(2V) */
extern_C int close (int fd); /* CLOSE(2V) */
#ifdef HAVE_FSYNC
extern_C int fsync (int fd); /* FSYNC(2) */
#endif
#if defined(HAVE_POLL)
  #include <poll.h>
  /* extern_C int poll (struct pollfd * fds, unsigned {int|long} nfds, int timeout); */
#endif
#if !defined(HAVE_SELECT) && defined(HAVE_POLL)
  #define NEED_OWN_SELECT /* select() can be implemented with poll()  */
  #ifndef _EMUL_SYS_TIME_H
    #define _EMUL_SYS_TIME_H
    struct timeval { long tv_sec; long tv_usec; };
    struct timezone { int tz_minuteswest; int tz_dsttime; };
  #endif
  #define SELECT_WIDTH_T int
  #define SELECT_SET_T fd_set
  #define SELECT_CONST
  #define HAVE_SELECT /* see unixaux.d */
#endif
#ifdef HAVE_SELECT
  #ifdef UNIX_BEOS
    #include <sys/socket.h>
  #endif
  #ifdef HAVE_SYS_SELECT_H
    #include <sys/select.h>
  #endif
  #ifndef FD_SETSIZE
    /* definition of types fd_set, err <sys/types.h> : */
    #ifdef UNIX_HPUX /* fd_set is defined, but FD_SETSIZE is not */
      #define fd_set  my_fd_set
    #endif
    #define FD_SETSIZE 256 /* maximum number of file descriptors */
    typedef int fd_mask; /* a bit group */
    #define NFDBITS (sizeof(fd_mask) * 8) /* number of bits in a bit group */
    typedef struct fd_set { fd_mask fds_bits[ceiling(FD_SETSIZE,NFDBITS)]; }
            fd_set;
    #define FD_SET(n,p)  ((p)->fds_bits[(n)/NFDBITS] |= bit((n)%NFDBITS))
    #define FD_CLR(n,p)  ((p)->fds_bits[(n)/NFDBITS] &= ~bit((n)%NFDBITS))
    #define FD_ISSET(n,p)  ((p)->fds_bits[(n)/NFDBITS] & bit((n)%NFDBITS))
    #define FD_ZERO(p)  bzero((char*)(p),sizeof(*(p)))
    #include <string.h>
    #define bzero(ptr,len)  memset(ptr,0,len)
  #endif
  extern_C int select (SELECT_WIDTH_T width, SELECT_SET_T* readfds,
                       SELECT_SET_T* writefds, SELECT_SET_T* exceptfds,
                       SELECT_CONST struct timeval * timeout); /* SELECT(2) */
#endif
#ifdef EINTR
/* wrapper around the system call, which intercepts and handles EINTR: */
extern int nonintr_open (const char* path, int flags, mode_t mode);
extern int nonintr_close (int fd);
#define OPEN nonintr_open
#define CLOSE nonintr_close
#else
#define OPEN open
#define CLOSE close
#endif
/* wrapper around the system call, get partial results and handle EINTR: */
extern ssize_t fd_read (int fd, void* buf, size_t nbyte, perseverance_t persev);
extern ssize_t fd_write (int fd, const void* buf, size_t nbyte, perseverance_t persev);
#define safe_read(fd,buf,nbyte)  fd_read(fd,buf,nbyte,persev_partial)
#define full_read(fd,buf,nbyte)  fd_read(fd,buf,nbyte,persev_full)
#define safe_write(fd,buf,nbyte)  fd_write(fd,buf,nbyte,persev_partial)
#define full_write(fd,buf,nbyte)  fd_write(fd,buf,nbyte,persev_full)
/* used by STREAM, PATHNAME, SPVW, MISC, UNIXAUX */

/* inquire the terminal, window size: */
extern_C int isatty (int fd); /* TTYNAME(3V) */
#if defined(HAVE_TERMIOS_H) && defined(HAVE_TCGETATTR) && defined(HAVE_TCSAFLUSH)
  #define UNIX_TERM_TERMIOS
  #include <termios.h> /* TERMIOS(3V) */
  /* extern_C int tcgetattr (int fd, struct termios * tp); */
  /* extern_C int tcsetattr (int fd, int optional_actions, const struct termios * tp); */
  /* extern_C int tcdrain (int fd); */ /* TERMIOS(3V) */
  /* extern_C int tcflush (int fd, int flag); */ /* TERMIOS(3V) */
  #undef TCSETATTR  /* eg. HP-UX 10 */
  #define TCSETATTR tcsetattr
  #define TCDRAIN tcdrain
  #define TCFLUSH tcflush
  #ifndef NCCS
    #define NCCS  sizeof(((struct termios *)0)->c_cc)
  #endif
  #if defined(WINSIZE_NEED_SYS_IOCTL_H) /* glibc2 needs this for "struct winsize" */
    #include <sys/ioctl.h>
  #elif defined(WINSIZE_NEED_SYS_PTEM_H) /* SCO needs this for "struct winsize" */
    #include <sys/stream.h>
    #include <sys/ptem.h>
  #endif
#elif defined(HAVE_SYS_TERMIO_H) || defined(HAVE_TERMIO_H)
  #define UNIX_TERM_TERMIO
  #if defined(HAVE_SYS_TERMIO_H)
    #include <sys/termio.h> /* TERMIO(4) */
  #elif defined(HAVE_TERMIO_H)
    #include <termio.h>
  #endif
  #ifndef NCCS
    #define NCCS  sizeof(((struct termio *)0)->c_cc)
  #endif
#elif defined(HAVE_SGTTY_H)
  /* compatibel to V7 or 4BSD, TIOC form ioctls.... */
  #define UNIX_TERM_SGTTY
  #include <sgtty.h>
  #include <sys/ioctl.h> /* TTY(4) */
#endif
#if defined(NEED_SYS_FILIO_H)
  #include <sys/filio.h>
#elif defined(NEED_SYS_IOCTL_H)
  #include <sys/ioctl.h>
#endif
#ifdef HAVE_IOCTL
  #ifdef IOCTL_DOTS
    extern_C int ioctl (int fd, IOCTL_REQUEST_T request, ...); /* IOCTL(2) */
    #define IOCTL_ARGUMENT_T  CADDR_T
  #else
    extern_C int ioctl (int fd, IOCTL_REQUEST_T request, IOCTL_ARGUMENT_T arg); /* IOCTL(2) */
    /* 3rd argument is always cast to type IOCTL_ARGUMENT_T (usually CADDR_T): */
    #define ioctl(fd,request,arg)  (ioctl)(fd,request,(IOCTL_ARGUMENT_T)(arg))
  #endif
#endif
#ifndef HAVE_SELECT
/* fcntl() will be used, declared in <fcntl.h> */
#endif
/* START_NO_BLOCK() & END_NO_BLOCK() should appear in pairs
   inside { NO_BLOCK_DECL(); ... };
 NO_BLOCK_DECL() should be before the first statement,
   but after the last declaration. */
#if defined(F_GETFL) && defined(O_NONBLOCK)
  /* non-blocking I/O a la SYSV */
  #define NO_BLOCK_DECL(handle)                                         \
    int fcntl_flags;                                                    \
    if ((fcntl_flags = fcntl(handle,F_GETFL,0))<0) { OS_error(); }
  #define START_NO_BLOCK(handle)  \
    do {                                                                  \
      if (fcntl(handle,F_SETFL,fcntl_flags|O_NONBLOCK)<0) { OS_error(); } \
    } while (0)
  #define END_NO_BLOCK(handle)  \
    do {                                                       \
      if (fcntl(handle,F_SETFL,fcntl_flags)<0) { OS_error(); } \
    } while(0)
#elif defined(F_GETFL) && defined(O_NDELAY)
  /* non-blocking I/O a la SYSV, older Unices called it O_NDELAY */
  #define NO_BLOCK_DECL(handle)                                         \
    int fcntl_flags;                                                    \
    if ((fcntl_flags = fcntl(handle,F_GETFL,0))<0) { OS_error(); }
  #define START_NO_BLOCK(handle)  \
    do {                                                                \
      if (fcntl(handle,F_SETFL,fcntl_flags|O_NDELAY)<0) { OS_error(); } \
    } while (0)
  #define END_NO_BLOCK(handle)  \
    do {                                                       \
      if (fcntl(handle,F_SETFL,fcntl_flags)<0) { OS_error(); } \
    } while(0)
#elif defined(FIONBIO)
  /* non-blocking I/O a la BSD 4.2 */
  /* semicolon in NO_BLOCK_DECL ensures no declaration after it */
  #define NO_BLOCK_DECL(handle)  \
    int non_blocking_io = 1;
  #define START_NO_BLOCK(handle)                                \
    if (ioctl(handle,FIONBIO,&non_blocking_io)) { OS_error(); }
  #define END_NO_BLOCK(handle)                                    \
    do {                                                          \
      non_blocking_io = 0;                                        \
      if (ioctl(handle,FIONBIO,&non_blocking_io)) { OS_error(); } \
    } while (0)
#endif

#if (defined(UNIX_TERM_TERMIOS) || defined(UNIX_TERM_TERMIO)) && !(defined(TCIFLUSH) && defined(TCOFLUSH))
  #define TCIFLUSH 0
  #define TCOFLUSH 1
#endif
extern_C int tgetent (const char* bp, const char* name); /* TERMCAP(3X) */
extern_C int tgetnum (const char* id); /* TERMCAP(3X) */
extern_C int tgetflag (const char* id); /* TERMCAP(3X) */
extern_C const char* tgetstr (const char* id, char** area); /* TERMCAP(3X) */
#ifdef EINTR
  /* wrapper around the system call, which intercepts and handles EINTR: */
  extern int nonintr_ioctl (int fd, IOCTL_REQUEST_T request, IOCTL_ARGUMENT_T arg);
  #undef ioctl
  #define ioctl(fd,request,arg)  nonintr_ioctl(fd,request,(IOCTL_ARGUMENT_T)(arg))
  #ifdef UNIX_TERM_TERMIOS
    extern int nonintr_tcsetattr (int fd, int optional_actions, struct termios * tp);
    extern int nonintr_tcdrain (int fd); /* TERMIOS(3V) */
    extern int nonintr_tcflush (int fd, int flag); /* TERMIOS(3V) */
    #undef TCSETATTR
    #define TCSETATTR nonintr_tcsetattr
    #undef TCDRAIN
    #define TCDRAIN nonintr_tcdrain
    #undef TCFLUSH
    #define TCFLUSH nonintr_tcflush
  #endif
#endif
/* used by SPVW, STREAM */

/* process date/time of day: */
#if defined(HAVE_GETTIMEOFDAY)
  #ifdef GETTIMEOFDAY_DOTS
    extern_C int gettimeofday (struct timeval * tp, ...); /* GETTIMEOFDAY(2) */
  #else
    extern_C int gettimeofday (struct timeval * tp, GETTIMEOFDAY_TZP_T tzp); /* GETTIMEOFDAY(2) */
  #endif
#elif defined(HAVE_FTIME)
  #include <sys/timeb.h>
  extern_C int ftime (struct timeb * tp); /* TIME(3V) */
  /* emulate gettimeofday() in unixaux.d: */
  #define NEED_OWN_GETTIMEOFDAY
  #ifndef _EMUL_SYS_TIME_H
    #define _EMUL_SYS_TIME_H
    struct timeval { long tv_sec; long tv_usec; };
    struct timezone { int tz_minuteswest; int tz_dsttime; };
  #endif
  extern int gettimeofday (struct timeval * tp, struct timezone * tzp); /* see unixaux.d */
#elif defined(HAVE_TIMES_CLOCK)
  #include <sys/times.h>
  extern_C clock_t times (struct tms * buffer); /* TIMES(3V) */
  extern_C time_t time (time_t* tloc); /* TIME(3V) */
#else
  #error "Cannot access real time with resolution finer than 1 second."
#endif
/* used by SPVW, MISC */

/* inquire used time of the process: */
#if defined(HAVE_GETRUSAGE)
  extern_C int getrusage (RUSAGE_WHO_T who, struct rusage * rusage); /* GETRUSAGE(2) */
  /* prototype useless, there 'struct rusage' /= 'struct rusage' */
#elif defined(HAVE_SYS_TIMES_H)
  #include <sys/param.h> /* define HZ, unit is 1/HZ seconds */
  #include <sys/times.h>
  extern_C clock_t times (struct tms * buffer); /* TIMES(3V) */
#endif
/* used by SPVW */

/* take a break for some time: */
extern_C unsigned int sleep (unsigned int seconds); /* SLEEP(3V) */
/* used by MISC */

/* program call: */
#define SHELL "/bin/sh"  /* the name of the shell command interpreter */
extern_C int pipe (int fd[2]); /* PIPE(2V) */
#ifdef HAVE_VFORK_H
  #include <vfork.h>
#endif
/* vfork() declared in <vfork.h> or <unistd.h> */
extern_C int dup2 (int oldfd, int newfd); /* DUP(2V) */
#if defined(HAVE_SETPGID)
  extern_C pid_t getpid (void); /* GETPID(2V) */
  extern_C int setpgid (pid_t pid, pid_t pgid); /* SETPGID(2V), SETSID(2V), TERMIO(4) */
  #define SETSID()  { register pid_t pid = getpid(); setpgid(pid,pid); }
#elif defined(HAVE_SETSID)
  extern_C pid_t setsid (void); /* SETSID(2V), TERMIO(4) */
  #define SETSID()  setsid()
#else
  #define SETSID()
#endif

/* exec*() functions are declared in <unistd.h> */

/* NB: In the period between vfork() and execv()/execl()/execlp() the child
   process may access only the data in the stack and constant data,
   because the parent process keeps running in this time already
   and can modify data in STACK, malloc() range, Lisp data range etc. */
#include <sys/wait.h>
extern_C pid_t waitpid (PID_T pid, int* statusp, int options); /* WAIT(2V) */
extern int wait2 (PID_T pid); /* see unixaux.d */
/* used by STREAM, PATHNAME, SPVW, UNIXAUX */

/* get random numbers: */
#ifndef rand /* some define rand() as macro... */
  extern_C int rand (void); /* RAND(3V) */
#endif
#if !defined(HAVE_SETPGID) /* in this case, already declared above */
  extern_C pid_t getpid (void); /* GETPID(2V) */
#endif
/* used by LISPARIT */

/* determine MACHINE-TYPE and MACHINE-VERSION and MACHINE-INSTANCE: */
#ifdef HAVE_SYS_UTSNAME_H
  #include <sys/utsname.h>
  extern_C int uname (struct utsname * buf); /* UNAME(2V) */
#endif
/* used by MISC */

/* determine MACHINE-INSTANCE: */
#ifdef HAVE_GETHOSTNAME
  /* extern_C int gethostname (char* name, size_t namelen); */ /* GETHOSTNAME(2) */
#endif
#ifdef HAVE_GETHOSTBYNAME
  #ifdef HAVE_NETDB_H
    #include <sys/socket.h>
    #include <netdb.h>
  #else
    #include <sun/netdb.h>
  #endif
/* gethostbyname() is declared in the above files */
#endif
#ifndef MAXHOSTNAMELEN
  #define MAXHOSTNAMELEN 64 /* see <sys/param.h> */
#endif
/* used by MISC */

/* work with sockets: */
#ifdef HAVE_GETHOSTBYNAME
  /* Type of a socket */
  #define SOCKET  int
  /* Error value for functions returning a socket */
  #define INVALID_SOCKET  (SOCKET)(-1)
  /* Error value for functions returning an `int' status */
  #define SOCKET_ERROR  (-1)
  /* Accessing the error code */
  #define sock_errno  errno
  #define sock_errno_is(val)  (errno == val)
  #define sock_set_errno(val)  (void)(errno = val)
  /* Signalling a socket-related error */
  #define SOCK_error()  OS_error()
  #ifdef UNIX_BEOS
    /* BeOS 5 sockets cannot be used like file descriptors.
       Reading and writing from a socket */
    extern ssize_t sock_read (int socket, void* buf, size_t size, perseverance_t persev);
    extern ssize_t sock_write (int socket, const void* buf, size_t size, perseverance_t persev);
    /* Closing a socket */
    /* extern int closesocket (int socket); */
  #else
    /* Reading and writing from a socket */
    #define sock_read(socket,buf,nbyte,persev)   fd_read(socket,buf,nbyte,persev)
    #define sock_write(socket,buf,nbyte,persev)  fd_write(socket,buf,nbyte,persev)
    /* Closing a socket */
    #define closesocket  close
  #endif
  /* Wrapping and unwrapping of a socket in a Lisp object */
  #define allocate_socket(fd)  allocate_handle(fd)
  #define TheSocket(obj)  TheHandle(obj)
#endif
/* used by SOCKET, STREAM */

/* shutdown(2) - older systems do not define SHUT_* */
#ifndef SHUT_RD
 #define SHUT_RD 0
#endif
#ifndef SHUT_WR
 #define SHUT_WR 1
#endif
#ifndef SHUT_RDWR
 #define SHUT_RDWR 2
#endif

/* Dynamic module loading:
   Even though HP-UX 10.20 and 11.00 support shl_load *and* dlopen,
   dlopen works correctly only with a patch. Because we want to avoid
   the situation where we build on a system with the patch but deploy
   on a system without, do not use dlopen on HP-UX. */
#ifdef UNIX_HPUX
  #undef HAVE_DLOPEN
#endif
#ifdef HAVE_DLOPEN
  #include <dlfcn.h>            /* declares dlopen,dlsym,dlclose,dlerror */
  #define HAVE_DYNLOAD
#endif

/* Character set conversion: */
#ifdef HAVE_ICONV
  #include <iconv.h>
  extern_C iconv_t iconv_open (const char * to_code, const char * from_code);
  extern_C size_t iconv (iconv_t cd, ICONV_CONST char * *inbuf, size_t *inbytesleft, char * *outbuf, size_t* outbytesleft);
  extern_C int iconv_close (iconv_t cd);
#endif

/* Interpretation of FILETIME structure: */
#ifdef UNIX_CYGWIN32
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
  #undef WIN32
  extern long time_t_from_filetime (const FILETIME * ptr);
  extern void time_t_to_filetime (time_t time_in, FILETIME * out);
#endif

/* close all file descriptors before exec() */
global void close_all_fd (void);

/* CLISP as a NeXTstep-GUI-Application: */
#ifdef NEXTAPP
/* Terminal-Stream, as nxterminal.m over the class LispServer implements it. */
  extern void nxterminal_send_output (void);
  extern void nxterminal_write_char (unsigned char ch);
  extern void nxterminal_write_string (unsigned char * string);
  extern unsigned char nxterminal_read_char (int* linepos);
  extern int nxterminal_unread_char (void);
  extern int nxterminal_listen (void);
  extern int nxterminal_init (void);
  extern int nxterminal_exit (void);
  extern int nxterminal_line_length;
#endif
