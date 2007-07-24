/*  This file was created by Configure. Any change made to it will be lost
 *  next time Configure is run. */
#ifndef __CONFIG_H__
#define __CONFIG_H__
#define UNIX
#define GPHELP "/Users/was/sage-v0.5.beta1-2005-08-17/source/../local/bin/gphelp"
#define GPDATADIR "/Users/was/sage-v0.5.beta1-2005-08-17/source/../local/share/pari"
#define SHELL_Q '\''

#define PARIVERSION "GP/PARI CALCULATOR Version 2.2.11 (development CHANGES-to-Aug)"
#ifdef __cplusplus
# define PARIINFO "C++ PowerPC running darwin (PPC kernel) 32-bit version"
#else
# define PARIINFO "PowerPC running darwin (PPC kernel) 32-bit version"
#endif
#define PARI_VERSION_CODE 131595
#define PARI_VERSION(a,b,c) (((a) << 16) + ((b) << 8) + (c))
#define PARI_VERSION_SHIFT 8

#define PARI_DOUBLE_FORMAT 0
#define GCC_VERSION "3.3 20030304 (Apple Computer, Inc. build 1671)"
#define ASMINLINE

/*  Location of GNU gzip program (enables reading of .Z and .gz files). */
#define GNUZCAT
#define ZCAT "/usr/bin/gzip -dc"

#define HAS_EXP2
#define HAS_LOG2
#define ULONG_NOT_DEFINED
#define USE_GETRUSAGE 1
#define USE_SIGSETMASK 1
#define HAS_DLOPEN
#define STACK_CHECK
#define HAS_VSNPRINTF
#define HAS_TIOCGWINSZ
#define HAS_STRFTIME
#define HAS_STAT
#endif
