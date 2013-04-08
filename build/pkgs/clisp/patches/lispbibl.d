/*
 * Main include-file for CLISP
 * Bruno Haible 1990-2006
 * Marcus Daniels 11.11.1994
 * Sam Steingold 1998-2006
 * German comments translated into English: Stefan Kain 2001-09-24

 Flags intended to be set through CFLAGS:
   Readline library:
     NO_READLINE
   Termcap/ncurses library:
     NO_TERMCAP_NCURSES
   Internationalization:
     NO_GETTEXT, UNICODE
   Fault handling:
     NO_SIGSEGV
   Foreign function interface:
     DYNAMIC_FFI
   Dynamic loading of modules:
     DYNAMIC_MODULES
   Safety level:
     SAFETY={0,1,2,3}
   Debugging (turned on by --with-debug configure option):
     DEBUG_GCSAFETY (requires G++)
     DEBUG_OS_ERROR
     DEBUG_SPVW
     DEBUG_BYTECODE
     DEBUG_BACKTRACE
     DEBUG_COMPILER
 Flags that may be set through CFLAGS, in order to override the defaults:
   Object representation (on 32-bit platforms only):
     TYPECODES, HEAPCODES, STANDARD_HEAPCODES, LINUX_NOEXEC_HEAPCODES, WIDE
   Advanced memory management:
     NO_SINGLEMAP, NO_TRIVIALMAP, NO_MULTIMAP_FILE, NO_MULTIMAP_SHM,
     NO_VIRTUAL_MEMORY, CONS_HEAP_GROWS_DOWN, CONS_HEAP_GROWS_UP,
     NO_MORRIS_GC, NO_GENERATIONAL_GC
   String representation:
     NO_SMALL_SSTRING


 Implementation is prepared for the following computers,
 operating systems and c-compilers
 (Only a rough listing, check the file PLATFORMS for further details.)
 Machine      Producer           Operating system              C-Compiler    recognized through
 AMIGA        Commodore          AMIGA-OS (AMIGADOS)           GNU           amiga or AMIGA, __GNUC__, maybe MC68000 or AMIGA3000
 any          any                UNIX                          GNU           unix, __GNUC__, ...
 any          any                UNIX                          CC            unix, ...
 Amiga 3000   Commodore          Amiga UNIX 2.1 SVR4.0         GNU           unix, __unix__, AMIX, __AMIX__, __svr4__, m68k, __m68k__, __motorola__, __GNUC__
 SUN-3        Sun                SUN-OS3 (UNIX BSD 4.2)        GNU           sun, unix, mc68020, __GNUC__
 SUN-3        Sun                SUN-OS4 (UNIX SUNOS 4.1)      GNU           sun, unix, mc68020, __GNUC__
 SUN-386      Sun                SUN-OS4 (UNIX SUNOS 4.0)      GNU           sun, unix, sun386, i386, __GNUC__
 SUN-386      Sun                SUN-OS4 (UNIX SUNOS 4.0)      CC            sun, unix, sun386, i386
 SUN-4        Sun                SUN-OS4 (UNIX SUNOS 4.1)      GNU           sun, unix, sparc, __GNUC__
 SUN-4        Sun                SUN-OS4 (UNIX SUNOS 4.1)      CC            sun, unix, sparc
 SUN-4        Sun                SUN-OS5 (UNIX Solaris)        GCC           sun, unix, sparc, __GNUC__
 UltraSparc   Sun                Solaris 7 (UNIX SUNOS 5.7)    CC            sun, unix, __sparc, __sparcv9
 UltraSparc   Sun                Solaris 7 (UNIX SUNOS 5.7)    GCC           sun, unix, __sparc, __arch64__, __GNUC__
 IBM-PC/386   any                SUN-OS5 (UNIX Solaris)        GCC           sun, unix, __svr4__, i386, __GNUC__
 HP9000-300   Hewlett-Packard    NetBSD 0.9 (UNIX BSD 4.3)     GNU           unix, __NetBSD__, mc68000, __GNUC__
 HP9000-300   Hewlett-Packard    HP-UX 8.0 (UNIX SYS V)        GNU           [__]hpux, [__]unix, [__]hp9000s300, mc68000, __GNUC__
 HP9000-800   Hewlett-Packard    HP-UX 8.0 (UNIX SYS V)        GNU           [__]hpux, [__]unix, [__]hp9000s800
 IRIS         Silicon Graphics   IRIX (UNIX SYS V 3.2)         GNU           unix, SVR3, mips, sgi, __GNUC__
 IRIS         Silicon Graphics   IRIX (UNIX SYS V)             cc -ansi      [__]unix, [__]SVR3, [__]mips, [__]sgi
 IRIS         Silicon Graphics   IRIX 5 (UNIX SYS V 4)         GNU           [__]unix, [__]SYSTYPE_SVR4, [__]mips, [__]host_mips, [__]MIPSEB, [__]sgi, __DSO__, [__]_MODERN_C, __GNUC__
 DECstation 5000                 RISC/OS (Ultrix V4.2A)        GNU           unix, [__]mips, [__]ultrix
 DG-UX 88k    Data General       DG/UX                         GNU           unix, m88000, DGUX
 DEC Alpha    DEC                OSF/1 1.3                     cc            [unix,] __unix__, __osf__, __alpha
 DEC Alpha    DEC                OSF/1 1.3                     GNU           unix, __unix__, __osf__, __alpha, __alpha__, _LONGLONG
 Apple MacII  Apple              A/UX (UNIX SYS V 2)           GNU           [__]unix, [__]AUX, [__]macII, [__]m68k, mc68020, mc68881, __GNUC__
 NeXT         NeXT               NeXTstep 3.1 (UNIX)           cc            NeXT, m68k; NEXTAPP for NeXTstep Application
 PowerPC      Apple              Mach 3.0 + MkLinux            GNU           unix, __powerpc__, __PPC__, _ARCH_PPC, _CALL_SYSV, __ELF__, __linux__
 PowerPC      Apple              Mach + Rhapsody               cc            __MACH__, __APPLE__, __ppc[__], __GNUC__, __APPLE_CC__
 PowerPC      Apple              Mach + MacOS X                cc            __MACH__, __APPLE__, __ppc__, __GNUC__, __APPLE_CC__
 Sequent      Sequent            PTX 3.2.0 V2.1.0 i386 (SYS V) GNU           unix, i386, _SEQUENT_, __GNUC__
 Sequent      Sequent            PTX V4.1.3                    GNU           unix, i386, _SEQUENT_, __svr4__, __GNUC__
 Convex C2    Convex             ConvexOS 10.1                 GNU           __convex__, __GNUC__
 IBM RS/6000  IBM                AIX 3.2                       GNU           _AIX, _AIX32, _IBMR2, __CHAR_UNSIGNED__, __GNUC__
 IBM-PC/386   any                LINUX (free UNIX)             GNU           unix, linux, i386, __GNUC__
 IBM-PC/386   any                LINUX (free UNIX)             Intel 5.0     __unix__, __linux__, __INTEL_COMPILER, __ICC, __USLC__
 IBM-PC/386   any                386BSD 0.1 (UNIX BSD 4.2)     GNU           unix, __386BSD__, i386, __GNUC__
 IBM-PC/386   any                NetBSD 0.9 (UNIX BSD 4.3)     GNU           unix, __NetBSD__, i386, __GNUC__
 IBM-PC/386   any                FreeBSD 4.0 (UNIX BSD 4.4)    GNU           unix, __FreeBSD__, i386, __GNUC__
 IBM-PC/386   any                EMX 0.9c (UNIXlike on OS/2)   GNU           [unix,] i386, __GNUC__, __EMX__
 IBM-PC/386   any                Cygwin32 on WinNT/Win95       GNU           _WIN32, __WINNT__, __CYGWIN32__, __POSIX__, _X86_, i386, __GNUC__
 IBM-PC/386   any                Mingw32 on WinNT/Win95        GNU           _WIN32, __WINNT__, __MINGW32__, _X86_, i386, __GNUC__
 IBM-PC/386   any                WinNT/Win95                   MSVC4.0,5.0   _WIN32, _M_IX86, _MSC_VER
 IBM-PC/386   any                WinNT/Win95                   Borland 5.0   __WIN32__, _M_IX86, __TURBOC__, __BORLANDC__
 IBM-PC/386   any                WinNT/Win95 and Cygwin32      GNU           _WIN32, __WINNT__, __CYGWIN32__, __POSIX__, __i386__, _X86_, __GNUC__
 IBM-PC/586   any                BeOS 5                        GNU           __BEOS__, __INTEL__, __i386__, _X86_, __GNUC__
 IBM-PC/586   any                HP NUE/ski, Linux             GNU           unix, linux, __ia64[__], __GNUC__, __LP64__
 RM400        Siemens-Nixdorf    SINIX-N 5.42                  c89           unix, mips, MIPSEB, host_mips, sinix, SNI, _XPG_IV
 Acorn        Risc PC            RISC OS 3.x                   GNU           [__]arm, [__]riscos, __GNUC__
 Acorn        Risc PC            RISC OS 3.x                   Norcroft      [__]arm, [__]riscos
 APPLE IIGS   Apple              ??                            ??
 For ANSI-C-Compiler: use pre-processors comment5, varbrace
   (and maybe gcc-cpp, ccpaux).
*/

/* this machine: WIN32 or GENERIC_UNIX */
#if (defined(__unix) || defined(__unix__) || defined(_AIX) || defined(sinix) || defined(__MACH__) || defined(__POSIX__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__BEOS__)) && !defined(unix)
  #define unix
#endif
#if (defined(_WIN32) && (defined(_MSC_VER) || defined(__MINGW32__))) || (defined(__WIN32__) && defined(__BORLANDC__))
  #undef WIN32                  /* because of __MINGW32__ */
  #define WIN32
#endif
#if !defined(WIN32)
  #if defined(unix)
    #define GENERIC_UNIX
  #else
    #error "Unknown machine type -- set machine again!"
  #endif
#endif
# additional specification of the machine:
#if defined(WIN32)
  # declare availability of typical PC facilities,
  # like a console with a graphics mode that differs from the text mode,
  # or a keyboard with function keys F1..F12.
  #define PC386 # IBMPC-compatible with 80386/80486-processor
#endif
#ifdef GENERIC_UNIX
  #if (defined(sun) && defined(unix) && defined(sun386))
    #define SUN386
  #endif
  #if (defined(unix) && (defined(linux) || defined(__CYGWIN32__) || defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__DragonFly__)) && (defined(i386) || defined(__i386__) || defined(__x86_64__) || defined(__amd64__)))
    #define PC386
  #endif
  #if (defined(sun) && defined(unix) && defined(mc68020))
    #define SUN3
  #endif
  #if (defined(sun) && defined(unix) && defined(sparc))
    #define SUN4
  #endif
  #if defined(sparc) || defined(__sparc__)
    # maybe SUN4_29 if only addresses <2^29 are supported
  #endif
  #if defined(hp9000s800) || defined(__hp9000s800)
    #define HP8XX
  #endif
#endif

# Determine the processor:
# MC680X0 == all processors of the Motorola 68000 series
# MC680Y0 == all processors of the Motorola 68000 series, starting at MC68020
# SPARC == the Sun SPARC processor
# HPPA == all processors of the HP Precision Architecture
# MIPS == all processors of the MIPS series
# M88000 == all processors of the Motorola 88000 series
# POWERPC == the IBM RS/6000 and PowerPC processor family
# I80386 == all processors of the Intel 8086 series, starting at 80386,
#           nowadays called IA32
# VAX == the VAX processor
# ARM == the ARM processor
# DECALPHA == the DEC Alpha superchip
# IA64 == the Intel IA-64 latecomer chip
# AMD64 == the AMD hammer chip
# S390 == the IBM S/390 processor
#if defined(__vax__)
  #define VAX
#endif
#if defined(arm) || defined(__arm) || defined(__arm__)
  #define ARM
#endif
#ifdef WIN32
  #if defined(_M_IX86) || defined(_X86_)
    #define I80386
  #endif
#endif
#ifdef GENERIC_UNIX
  #if defined(m68k) || defined(__m68k__) || defined(mc68000)
    #define MC680X0
  #endif
  #if defined(mc68020) || defined(__mc68020__) || (defined(m68k) && defined(NeXT))
    #define MC680X0
    #define MC680Y0
  #endif
  #if defined(i386) || defined(__i386) || defined(__i386__) || defined(_I386)
    #define I80386
  #endif
  #if defined(sparc) || defined(__sparc__)
    #define SPARC
    #if defined(__sparcv9) || defined(__arch64__)
      #define SPARC64
    #endif
  #endif
  #if defined(mips) || defined(__mips) || defined(__mips__)
    #define MIPS
    #if defined(_MIPS_SZLONG)
      #if (_MIPS_SZLONG == 64)
        # We should also check for (_MIPS_SZPTR == 64), but gcc keeps this at 32.
        #define MIPS64
      #endif
    #endif
  #endif
  #if defined(HP8XX) || defined(hppa) || defined(__hppa) || defined(__hppa__)
    #define HPPA
  #endif
  #if defined(m88000) || defined(__m88k__)
    #define M88000
  #endif
  #if defined(_IBMR2) || defined(__powerpc) || defined(__ppc) || defined(__ppc__) || defined(__powerpc__)
    #define POWERPC
  #endif
  #ifdef __alpha
    #define DECALPHA
  #endif
  #ifdef __ia64__
    #define IA64
  #endif
  #if defined(__x86_64__) || defined(__amd64__)
    #define AMD64
  #endif
  #ifdef __s390__
    #define S390
  #endif
#endif


# Selection of the operating system
#ifdef WIN32
  # Windows NT, Windows 95
  #define WIN32_NATIVE  # native NT API, no DOS calls
#endif
#ifdef GENERIC_UNIX
  #define UNIX
  #ifdef __linux__
    #define UNIX_LINUX  # Linux (Linus Torvalds Unix)
  #endif
  #ifdef __GNU__
    #define UNIX_HURD  # the GNU system (Hurd + glibc)
  #endif
  #ifdef __NetBSD__
    #define UNIX_NETBSD
  #endif
  #if defined(__FreeBSD__) || defined(__DragonFly__)
    # FreeBSD or its fork called DragonFly BSD.
    #define UNIX_FREEBSD
  #endif
  #ifdef __OpenBSD__
    #define UNIX_OPENBSD
  #endif
  #if defined(hpux) || defined(__hpux)
    #define UNIX_HPUX  # HP-UX
  #endif
  #if defined(SVR3) || defined(__SVR3) || defined(SVR4) || defined(__SVR4) || defined(SYSTYPE_SVR4) || defined(__SYSTYPE_SVR4) || defined(__svr4__) || defined(USG) || defined(UNIX_HPUX) # ??
    #define UNIX_SYSV  # UNIX System V
  #endif
  #if defined(UNIX_SYSV) && (defined(sgi) || defined(__sgi))
    #define UNIX_IRIX  # Irix
    #if defined(SYSTYPE_SVR4) || defined(__SYSTYPE_SVR4)
      #define UNIX_IRIX5  # Irix 5
    #endif
  #endif
  #if defined(MIPS) && (defined(ultrix) || defined(__ultrix))
    #define UNIX_DEC_ULTRIX  # DEC's (or IBM's ?) RISC/OS Ultrix on DEC MIPS
    #ifdef __GNUC__
      #define UNIX_DEC_ULTRIX_GCCBUG  # work around a bug in GCC 2.3.3
    #endif
  #endif
  #if defined(MIPS) && defined(sinix) # && defined(SNI)
    #define UNIX_SINIX # Siemens is nix
  #endif
  #if defined(USL) || (defined(__svr4__) && defined(I80386) && !defined(__sun))
    # A couple of Unices for 386s (all running under different names)
    # derive from USL SysV R 4:
    #   386 UHC UNIX System V release 4
    #   Consensys System V 4.2
    #   Onsite System V 4.2
    #   SINIX-Z
    #   DYNIX/ptx V4.1.3
    #   SunOS 5
    #define UNIX_SYSV_USL  # Unix System V R 4 by AT&T's subsidiary USL
    #define UNIX_SYSV_UHC_1 # treat like HPPA && UNIX_HPUX
    #ifdef SNI
      #define UNIX_SINIX # Siemens is nix
    #endif
  #endif
  #if defined(_SEQUENT_) && !defined(__svr4__)
    #define UNIX_SYSV_PTX  # Dynix/ptx v. 2 or 3
  #endif
  #ifdef _AIX
    #define UNIX_AIX  # IBM AIX
  #endif
  #ifdef DGUX
    #define UNIX_DGUX  # Data General DG/UX
  #endif
  #ifdef __osf__
    #define UNIX_OSF  # OSF/1
  #endif
  #ifdef AUX
    #define UNIX_AUX  # Apple A/UX, a spiced-up SVR2
  #endif
  #ifdef NeXT
    #define UNIX_NEXTSTEP  # NeXTstep
    # define NEXTAPP       # Define this to get a NeXTstep-GUI application

    #define MAYBE_NEXTAPP  # a little hack, to make the .mem files compatible
                           # between CLISP with NEXTAPP and without NEXTAPP
  #endif
  #if defined(__APPLE__) && defined(__MACH__)
    #define UNIX_MACOSX  # MacOS X
  #endif
  #ifdef AMIX
    #define UNIX_AMIX  # Amiga UNIX
  #endif
  #ifdef __CYGWIN32__
    #define UNIX_CYGWIN32  # Cygwin32 (UNIXlike on WinNT/Win95)
  #endif
  #ifdef __BEOS__
    #define UNIX_BEOS  # BeOS (UNIXlike)
  #endif
#endif
%% #ifdef WIN32_NATIVE
%%   puts("#define WIN32_NATIVE");
%% #endif
%% #ifdef UNIX
%%   puts("#define UNIX");
%% #endif


# Determine properties of compiler and environment:
#if defined(UNIX) || defined(__MINGW32__)
  #include "unixconf.h"  # configuration generated by configure
  #include "intparam.h"  # integer-type characteristics created by the machine
  #include "floatparam.h" # floating-point type characteristics
#elif defined(WIN32) && !defined(__MINGW32__)
  #include "version.h"          /* defines PACKAGE_* */
  #define char_bitsize 8
  #define short_bitsize 16
  #define int_bitsize 32
  #if defined(I80386)
    #define long_bitsize 32
  #elif defined(DECALPHA)
    #define long_bitsize 64
  #endif
  #ifdef __GNUC__
    #if (__GNUC__ >= 2) # GCC 2 got a working 'long long' type meanwhile.
      #define HAVE_LONGLONG
      #define long_long_bitsize 64
    #endif
  #endif
  #if defined(I80386)
    #define pointer_bitsize 32
  #elif defined(DECALPHA)
    #define pointer_bitsize 64
  #endif
  #define alignment_long 4
  #if defined(I80386) || defined(VAX) || defined(ARM) || defined(DECALPHA)
    #define short_little_endian
    #define long_little_endian
  #endif
  #define stack_grows_down
  #define CODE_ADDRESS_RANGE 0
  #define MALLOC_ADDRESS_RANGE 0
  #define SHLIB_ADDRESS_RANGE 0
  #define STACK_ADDRESS_RANGE ~0UL
  #define ICONV_CONST
#else
  #error "where is the configuration for your platform?"
#endif


# A more precise classification of the operating system:
#if defined(UNIX) && defined(SIGNALBLOCK_BSD) && !defined(SIGNALBLOCK_SYSV)
  #define UNIX_BSD  # BSD Unix
#endif
#if (defined(SUN3) || defined(SUN386) || defined(SUN4)) && defined(HAVE_MMAP) && defined(HAVE_VADVISE)
  #define UNIX_SUNOS4  # Sun OS Version 4
#endif
#if (defined(SUN4) || (defined(I80386) && defined(__svr4__) && defined(__sun))) && !defined(HAVE_VADVISE) # && !defined(HAVE_GETPAGESIZE)
  #define UNIX_SUNOS5  # Sun OS Version 5.[1-5] (Solaris 2)
#endif
#if defined(UNIX_MACOSX) && !defined(HAVE_MSYNC)
  #define UNIX_RHAPSODY  # MacOS X Server, a.k.a. Rhapsody
#endif
#if defined(UNIX_MACOSX) && defined(HAVE_MSYNC)
  #define UNIX_DARWIN  # MacOS X, a.k.a. Darwin
#endif


# Choose the character set:
#if defined(UNIX) || defined(WIN32)
  #define ISOLATIN_CHS  # ISO 8859-1, see isolatin.chs
  # Most Unix systems today support the ISO Latin-1 character set, in
  # particular because they have X11 and the X11 fonts are in ISO Latin-1.
  # Exceptions below.
  # On Win32, the standard character set is ISO-8859-1. Only the DOS box
  # displays CP437, but we convert from ISO-8859-1 to CP437 in the
  # low-level output routine full_write().
#endif
#ifdef UNIX_BEOS
  # The default encoding on BeOS is UTF-8, not ISO 8859-1.
  # If compiling with Unicode support, we use it. Else fall back to ASCII.
  #undef ISOLATIN_CHS
  #ifdef UNICODE
    #define UTF8_CHS  # UTF-8
  #endif
#endif
#ifdef HP8XX
  #undef ISOLATIN_CHS
  #define HPROMAN8_CHS  # HP-Roman8, see hproman8.chs
  # under X-Term however: #define ISOLATIN_CHS ??
#endif
#ifdef UNIX_NEXTSTEP
  #undef ISOLATIN_CHS
  #define NEXTSTEP_CHS  # NeXTstep, see nextstep.chs
#endif
#if !(defined(ISOLATIN_CHS) || defined(HPROMAN8_CHS) || defined(NEXTSTEP_CHS))
  #define ASCII_CHS  # Default: plain ASCII charset without special chars
#endif


# Choose the compiler:
#if defined(__GNUC__)
  #define GNU
#endif
#if defined(__STDC__) || defined(__BORLANDC__) || defined(__cplusplus)
  # ANSI C compilers define __STDC__ (but some define __STDC__=0 !).
  # Borland C has an ANSI preprocessor and compiler, but fails to define
  # __STDC__.
  # HP aCC is an example of a C++ compiler which defines __cplusplus but
  # not __STDC__.
  #define ANSI
#endif
#if defined(_MSC_VER)
  #define MICROSOFT
#endif
#if defined(__BORLANDC__)
  #define BORLAND
#endif
#if defined(__INTEL_COMPILER)
  #define INTEL
#endif


# Selection of floating-point capabilities:
# FAST_DOUBLE should be defined if there is a floating-point coprocessor
# with a 'double'-type IEEE-Floating-Points with 64 Bits.
# FAST_FLOAT should be defined if there is a floating-point co-processor
# with a 'float'-type IEEE-Floating-Points with 32 Bits,
# and a C-Compiler that generates 'float'-operations
# instead of 'double'-operations
#if (float_mant_bits == 24) && (float_rounds == rounds_to_nearest) && float_rounds_correctly && !defined(FLOAT_OVERFLOW_EXCEPTION) && !defined(FLOAT_UNDERFLOW_EXCEPTION) && !defined(FLOAT_INEXACT_EXCEPTION)
  #define FAST_FLOAT
#endif
#if (double_mant_bits == 53) && (double_rounds == rounds_to_nearest) && double_rounds_correctly && !defined(DOUBLE_OVERFLOW_EXCEPTION) && !defined(DOUBLE_UNDERFLOW_EXCEPTION) && !defined(DOUBLE_INEXACT_EXCEPTION)
  #define FAST_DOUBLE
#endif
#ifdef ARM
  # The processor is little-endian w.r.t. integer types but stores 'double'
  # floats in big-endian word order!
  #undef FAST_DOUBLE
#endif
#ifdef NO_FAST_DOUBLE
  #undef FAST_DOUBLE
#endif
#ifdef NO_FAST_FLOAT
  #undef FAST_FLOAT
#endif


# Selection of the language:
#ifdef ENGLISH
  #undef ENGLISH
  #define ENGLISH 1
  #define LANGUAGE_STATIC
#endif


# Selection of the safety-level:
# SAFETY=0 : all optimizations are turned on
# SAFETY=1 : all optimizations on, but keep STACKCHECKs
# SAFETY=2 : only simple assembler-support
# SAFETY=3 : no optimizations
#ifndef SAFETY
  #define SAFETY 0
#endif
#if SAFETY >= 3
  #define NO_ASM
  #define NO_ARI_ASM
  #define NO_FAST_DISPATCH
#endif

# We don't support pre-ANSI-C compilers any more.
#if !defined(ANSI)
  #error "An ANSI C or C++ compiler is required to compile CLISP!"
#endif

# gcc-2.7.2 has a bug: it interpretes `const' as meaning `not modified by
# other parts of the program', and thus miscompiles at least justify_empty_2
# and pr_enter_1 in io.d.
#if defined(GNU) && (__GNUC__ == 2) && (__GNUC_MINOR__ == 7)
  #undef const
  #define const
  #define __const const
  # We define __const to const, not to empty, to avoid warnings on
  # UNIX_RHAPSODY, which unconditionally defines __const to const when
  # <sys/cdefs.h> is included via <setjmp.h> below.
  #ifdef MULTITHREAD
    #warning "Multithreading will not be efficient because of a workaround to a gcc bug."
    #warning "Get a newer version of gcc."
  #endif
#endif

# A property of the processor:
# The sequence in which words/long-words are being put into bytes
#if defined(short_little_endian) || defined(int_little_endian) || defined(long_little_endian)
  # Z80, VAX, I80386, DECALPHA, MIPSEL, IA64, AMD64, ...:
  # Low Byte is the lowest, High Byte in a higher address
  #if defined(BIG_ENDIAN_P)
    #error "Bogus BIG_ENDIAN_P -- set BIG_ENDIAN_P again!"
  #endif
  #define BIG_ENDIAN_P  0
#endif
#if defined(short_big_endian) || defined(int_big_endian) || defined(long_big_endian)
  # MC680X0, SPARC, HPPA, MIPSEB, M88000, POWERPC, S390, ...:
  # High Byte is the lowest, Low Byte is a higher adress (easier to read)
  #if defined(BIG_ENDIAN_P)
    #error "Bogus BIG_ENDIAN_P -- set BIG_ENDIAN_P again"
  #endif
  #define BIG_ENDIAN_P  1
#endif
#if !defined(BIG_ENDIAN_P)
  #error "Bogus BIG_ENDIAN_P -- set BIG_ENDIAN_P again!"
#endif
%% export_def(BIG_ENDIAN_P);

# A property of the processor (and C compiler): The alignment of C functions.
# (See gcc's machine descriptions, macro FUNCTION_BOUNDARY, for information.)
#if defined(IA64)
  #define C_CODE_ALIGNMENT  16
  #define log2_C_CODE_ALIGNMENT  4
#endif
#if defined(DECALPHA)
  #define C_CODE_ALIGNMENT  8
  #define log2_C_CODE_ALIGNMENT  3
#endif
#if (defined(I80386) && defined(GNU)) || defined(SPARC) || defined(MIPS) || defined(M88000) || defined(POWERPC) || defined(ARM) || defined(AMD64) || defined(S390)
  # When using gcc on i386, this assumes that -malign-functions has not been
  # used to specify an alignment smaller than 4 bytes.
  #define C_CODE_ALIGNMENT  4
  #define log2_C_CODE_ALIGNMENT  2
#endif
#if defined(HPPA)
  # A function pointer on hppa is either
  #   - a code pointer == 0 mod 4, or
  #   - a pointer to a two-word structure (first word: a code pointer,
  #     second word: a value which will be put in register %r19),
  #     incremented by 2, hence == 2 mod 4.
  # The current compilers only emit the second kind of function pointers,
  # hence we can assume that all function pointers are == 2 mod 4.
  #define C_CODE_ALIGNMENT  2
  #define log2_C_CODE_ALIGNMENT  1
#endif
#if defined(MC680X0)
  #define C_CODE_ALIGNMENT  2
  #define log2_C_CODE_ALIGNMENT  1
#endif
#if !defined(C_CODE_ALIGNMENT) # e.g. (defined(I80386) && defined(MICROSOFT))
  #define C_CODE_ALIGNMENT  1
  #define log2_C_CODE_ALIGNMENT  0
#endif


# Flags for the system's include files.
#ifdef MULTITHREAD
  #if defined(UNIX_GNU) || defined(UNIX_SUNOS5)
    #define _REENTRANT
  #endif
#endif


# Width of object representation:
# WIDE means than an object (pointer) occupies 64 bits (instead of 32 bits).
# WIDE_HARD means on a 64-bit platform.
# WIDE_SOFT means on a 32-bit platform, each object pointer occupies 2 words.
# WIDE_AUXI means on a 32-bit platform, each object occupies 2 words, the
#           pointer and some auxiliary word.

#if defined(DECALPHA) || defined(MIPS64) || defined(SPARC64) || defined(IA64) || defined(AMD64)
  #define WIDE_HARD
#endif

#if defined(WIDE_SOFT_LARGEFIXNUM) && !defined(WIDE_SOFT)
  #define WIDE_SOFT
#endif

#if defined(WIDE) && !(defined(WIDE_HARD) || defined(WIDE_SOFT) || defined(WIDE_AUXI))
  #define WIDE_SOFT
#endif
#if (defined(WIDE_HARD) || defined(WIDE_SOFT) || defined(WIDE_AUXI)) && !defined(WIDE)
  #define WIDE
#endif
# Now: defined(WIDE) == defined(WIDE_HARD) || defined(WIDE_SOFT) || defined(WIDE_AUXI)


/* Global register declarations.
   Speed benefit: Just putting the STACK into a register, brought 5% of speed
   around 1992. Now, with an AMD Athlon CPU from 2000, with good caches, it
   still brings 4%.
   The declarations must occur before any system include files define any
   inline function, which is the case on UNIX_DGUX and UNIX_GNU.
   Only GCC supports global register variables. Not Apple's variant of GCC.
   And only the C frontend, not the C++ frontend, understands the syntax.
   And gcc-3.0 to 3.3 has severe bugs with global register variables, see
   CLISP bugs 710737 and 723097 and
   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=7871
   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=10684
   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=14937
   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=14938 */
#if defined(GNU) && !(__APPLE_CC__ > 1) && !defined(__cplusplus) && !(__GNUC__ == 3 && __GNUC_MINOR__ < 4) && !defined(MULTITHREAD) && (SAFETY < 2)
/* Overview of use of registers in gcc terminology:
 fixed: mentioned in FIXED_REGISTERS
 used:  mentioned in CALL_USED_REGISTERS but not FIXED_REGISTERS
                     (i.e. caller-saved)
 save:  otherwise (i.e. call-preserved, callee-saved)

               STACK    mv_count  value1   back_trace
 MC680X0       used
 I80386        save
 SPARC (gcc2)  fixed    fixed     fixed    used
 MIPS
 HPPA          save     save      save     save
 M88000        save     save      save
 ARM           save
 DECALPHA      save     save      save
 IA64
 AMD64
 S390          save

 Special notes:
 - gcc3/Sparc (Linux & Solaris) handles registers differently from gcc2. FIXME
 - If STACK is in a "used"/"save" register, it needs to be saved into
   saved_STACK upon begin_call(), so that asynchronous interrupts will
   be able to restore it.
 - All of the "used" registers need to be backed up upon begin_call()
   and restored during end_call().
 - All of the "save" registers need to be backed up upon begin_callback()
   and restored during end_callback().
 - When the interpreter does a longjmp(), the registers STACK, mv_count,
   value1 may need to be temporarily saved. This is highly machine
   dependent and is indicated by the NEED_temp_xxxx macros. */

  # Register for STACK.
  #if defined(MC680X0)
    #define STACK_register "a4" # highest address register after sp=A7,fp=A6/A5
  #endif
  #if defined(I80386) && !defined(UNIX_BEOS) && !defined(DYNAMIC_MODULES)
    # On BeOS, everything is compiled as PIC, hence %ebx is already booked.
    # If DYNAMIC_MODULES is defined, external modules are compiled as PIC,
    # which is why %ebx is already in use.
    #if (__GNUC__ >= 2) # The register names have changed
      #define STACK_register  "%ebx"  # one of the call-saved registers without special hardware commands
    #else
      #define STACK_register  "bx"
    #endif
  #endif
  #if defined(SPARC)
    #define STACK_register  "%g5"  # a global register
  #endif
  #if defined(HPPA) && (__GNUC__*100 + __GNUC_MINOR__ >= 2*100+7) # gcc versions earlier than 2.7 had bugs
    #define STACK_register  "%r10"  # one of the general registers %r5..%r18
  #endif
  #if defined(M88000)
    #define STACK_register  "%r14"  # one of the general registers %r14..%r25
  #endif
  #if defined(ARM)
    #define STACK_register  "%r8"   # one of the general registers %r4..%r8
  #endif
  #if defined(DECALPHA)
    #define STACK_register  "$9"    # one of the general registers $9..$14
  #endif
  #if defined(S390) && ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1)))
    # global register assignment did not work on s390 until gcc 3.1
    #define STACK_register  "9"     # one of the general registers %r8..%r9
  #endif
  # What about NEED_temp_STACK ?? Needed if STACK is in a "used" register??
  # Register for mv_count.
  #if defined(SPARC)
    #define mv_count_register  "%g6"
    #if defined(UNIX_NETBSD)
      #define NEED_temp_mv_count
    #endif
  #endif
  #if defined(HPPA)
    #define mv_count_register  "%r11"  # one of the general registers %r5..%r18
    #define NEED_temp_mv_count
  #endif
  #if defined(M88000)
    #define mv_count_register  "%r15" # one of the general registers %r14..%r25
    #define NEED_temp_mv_count
  #endif
  #if defined(DECALPHA)
    #define mv_count_register  "$10"  # one of the general registers $9..$14
    #define NEED_temp_mv_count
  #endif
  # Register for value1.
  #if !(defined(WIDE) && !defined(WIDE_HARD))
    #if defined(SPARC)
      #define value1_register  "%g7"
      #if defined(UNIX_NETBSD)
        #define NEED_temp_value1
      #endif
    #endif
    #if defined(HPPA)
      #define value1_register  "%r12"  # one of the general registers %r5..%r18
      #define NEED_temp_value1
    #endif
    #if defined(M88000)
      #define value1_register  "%r16"  # one of the general registers %r14..%r25
      #define NEED_temp_value1
    #endif
    #if defined(DECALPHA)
      #define value1_register  "$11"  # one of the general registers $9..$14
      #define NEED_temp_value1
    #endif
  #endif
  # Register for back_trace.
  #if !(defined(WIDE) && !defined(WIDE_HARD))
    #if defined(SPARC)
      #define back_trace_register  "%g4"  # a global register
      # %g4 seems to be a scratch-register as of lately with gcc 2.3
      # This causes problems with libc.so.1.6.1 (and higher) (in getwd())
      # That's why HAVE_SAVED_back_trace has been defined above.
    #endif
    #if defined(HPPA)
      #define back_trace_register  "%r13"  # one of the general registers  %r5..%r18
    #endif
  #endif
  # Declare the registers now (before any system include file which could
  # contain some inline functions).
  #ifdef STACK_register
    register long STACK_reg __asm__(STACK_register);
  #endif
  #ifdef mv_count_register
    register long mv_count_reg __asm__(mv_count_register);
  #endif
  #ifdef value1_register
    register long value1_reg __asm__(value1_register);
  #endif
  #ifdef back_trace_register
    register long back_trace_reg __asm__(back_trace_register);
  #endif
  # Saving "save" registers.
  #if (defined(I80386) || defined(HPPA) || defined(M88000) || defined(ARM) || defined(DECALPHA) || defined(S390)) && (defined(STACK_register) || defined(mv_count_register) || defined(value1_register) || defined(back_trace_register))
    #define HAVE_SAVED_REGISTERS
    struct registers {
      #ifdef STACK_register
        long STACK_register_contents;
      #endif
      #ifdef mv_count_register
        long mv_count_register_contents;
      #endif
      #ifdef value1_register
        long value1_register_contents;
      #endif
      #ifdef back_trace_register
        long back_trace_register_contents;
      #endif
    };
    #ifndef MULTITHREAD
      extern struct registers * callback_saved_registers;
    #else
      #define callback_saved_registers  (current_thread()->_callback_saved_registers)
    #endif
    #ifdef STACK_register
      #define SAVE_STACK_register(registers)     \
              registers->STACK_register_contents = STACK_reg
      #define RESTORE_STACK_register(registers)  \
              STACK_reg = registers->STACK_register_contents
    #else
      #define SAVE_STACK_register(registers)
      #define RESTORE_STACK_register(registers)
    #endif
    #ifdef mv_count_register
      #define SAVE_mv_count_register(registers)     \
              registers->mv_count_register_contents = mv_count_reg
      #define RESTORE_mv_count_register(registers)  \
              mv_count_reg = registers->mv_count_register_contents
    #else
      #define SAVE_mv_count_register(registers)
      #define RESTORE_mv_count_register(registers)
    #endif
    #ifdef value1_register
      #define SAVE_value1_register(registers)     \
              registers->value1_register_contents = value1_reg
      #define RESTORE_value1_register(registers)  \
              value1_reg = registers->value1_register_contents
    #else
      #define SAVE_value1_register(registers)
      #define RESTORE_value1_register(registers)
    #endif
    #ifdef back_trace_register
      #define SAVE_back_trace_register(registers)     \
              registers->back_trace_register_contents = back_trace_reg
      #define RESTORE_back_trace_register(registers)  \
              back_trace_reg = registers->back_trace_register_contents
    #else
      #define SAVE_back_trace_register(registers)
      #define RESTORE_back_trace_register(registers)
    #endif
    #define SAVE_REGISTERS(inner_statement)                                  \
      do {                                                                   \
        var struct registers * registers = alloca(sizeof(struct registers)); \
        SAVE_STACK_register(registers);                                      \
        SAVE_mv_count_register(registers);                                   \
        SAVE_value1_register(registers);                                     \
        SAVE_back_trace_register(registers);                                 \
        inner_statement;                                                     \
        { var gcv_object_t* top_of_frame = STACK;                            \
          pushSTACK(fake_gcv_object((aint)callback_saved_registers));        \
          finish_frame(CALLBACK);                                            \
        }                                                                    \
        callback_saved_registers = registers;                                \
      } while(0)
    #define RESTORE_REGISTERS(inner_statement)                                \
      do {                                                                    \
        var struct registers * registers = callback_saved_registers;          \
        if (!(framecode(STACK_0) == CALLBACK_frame_info)) abort();            \
        callback_saved_registers = (struct registers *)(aint)as_oint(STACK_1);\
        skipSTACK(2);                                                         \
        inner_statement;                                                      \
        RESTORE_STACK_register(registers);                                    \
        RESTORE_mv_count_register(registers);                                 \
        RESTORE_value1_register(registers);                                   \
        RESTORE_back_trace_register(registers);                               \
      } while(0)
  #endif
  # Saving the STACK (for asynchronous interrupts).
  # If STACK is a global variable or lies in a register which is left
  # untouched by operating system and library (this is the case on SUN4),
  # we don't need to worry about it.
  #if defined(STACK_register) && !defined(SUN4)
    #define HAVE_SAVED_STACK
  #endif
  # Saving "used" registers.
  #if defined(mv_count_register) && 0
    #define HAVE_SAVED_mv_count
  #endif
  #if defined(value1_register) && 0
    #define HAVE_SAVED_value1
  #endif
  #if defined(back_trace_register) && defined(SPARC)
    #define HAVE_SAVED_back_trace
  #endif
#endif
#ifndef HAVE_SAVED_REGISTERS
  #define SAVE_REGISTERS(inner_statement)
  #define RESTORE_REGISTERS(inner_statement)
#endif
%% #ifdef HAVE_SAVED_REGISTERS
%%   puts("#ifndef IN_MODULE_CC");
%%   #ifdef STACK_register
%%     printf("register long STACK_reg __asm__(\"%s\");\n",STACK_register);
%%   #endif
%%   #ifdef mv_count_register
%%     printf("register long mv_count_reg __asm__(\"%s\");\n",mv_count_register);
%%   #endif
%%   #ifdef value1_register
%%     printf("register long value1_reg __asm__(\"%s\");\n",value1_register);
%%   #endif
%%   #ifdef back_trace_register
%%     printf("register long back_trace_reg __asm__(\"%s\");\n",back_trace_register);
%%   #endif
%%   printf("struct registers { ");
%%   #ifdef STACK_register
%%     printf("long STACK_register_contents; ");
%%   #endif
%%   #ifdef mv_count_register
%%     printf("long mv_count_register_contents; ");
%%   #endif
%%   #ifdef value1_register
%%     printf("long value1_register_contents; ");
%%   #endif
%%   #ifdef back_trace_register
%%     printf("long back_trace_register_contents; ");
%%   #endif
%%   puts("};");
%%   puts("extern struct registers * callback_saved_registers;");
%%   puts("#endif");
%% #endif

#define VALUES_IF(cond)                         \
  do { value1 = (cond) ? T : NIL; mv_count = 1; } while (0)
%% export_def(VALUES_IF(C));

#define VALUES0                                 \
  do { value1 = NIL; mv_count = 0; } while (0)
%% export_def(VALUES0);

#define VALUES1(A)                               \
  do { value1 = (A); mv_count = 1; } while (0)
%% export_def(VALUES1(A));

#define VALUES2(A,B)                            \
  do { value1 = (A); value2 = (B); mv_count = 2; } while (0)
%% export_def(VALUES2(A,B));

#define VALUES3(A,B,C)                          \
  do { value1 = (A); value2 = (B); value3 = (C); mv_count = 3; } while (0)
%% export_def(VALUES3(A,B,C));

# ###################### Macros for C ##################### #

#if !defined(return_void)
  # To return a type of value void: return_void(...);
  #ifdef GNU
    #define return_void  return # 'return void;' is admissible
  #else
    # In general it is not legal to return `void' values.
    #define return_void  # Don't use 'return' for expressions of type 'void'.
  #endif
#endif
#if defined(GNU) && defined(__GNUG__)
  # Although legal, g++ warns about 'return void;'. Shut up the warning.
  #undef return_void
  #define return_void
#endif

#if !defined(GNU) && !defined(inline)
  #define inline      # inline foo() {...} --> foo() {...}
#endif
%% puts("#if !defined(__GNUC__) && !defined(inline)");
%% puts("#define inline");
%% puts("#endif");

# Definitions for C++-Compilers:
#ifdef __cplusplus
  #define BEGIN_DECLS  extern "C" {
  #define END_DECLS    }
#else
  #define BEGIN_DECLS
  #define END_DECLS
#endif
%% export_def(BEGIN_DECLS);
%% export_def(END_DECLS);

# Empty macro-arguments:
# Some compilers (ie. cc under HP-UX) seem to interpret a macro call
# foo(arg1,...,argn,) as equivalent to foo(arg1,...,argn), which will
# yield an error. _EMA_ stands for "empty macro argument".
# It will be inserted by CC_NEED_DEEMA,
# each time between comma and closing parentheses.
# It is also needed when potentially empty arguments
# are returned to other macros

#define _EMA_

# Concatenation of two macro-expanded tokens:
# Example:
#   #undef x
#   #define y 16
#   CONCAT(x,y)        ==>  'x16' (not 'xy' !)
#define CONCAT_(xxx,yyy)  xxx##yyy
#define CONCAT3_(aaa,bbb,ccc)  aaa##bbb##ccc
#define CONCAT4_(aaa,bbb,ccc,ddd)  aaa##bbb##ccc##ddd
#define CONCAT5_(aaa,bbb,ccc,ddd,eee)  aaa##bbb##ccc##ddd##eee
#define CONCAT6_(aaa,bbb,ccc,ddd,eee,fff)  aaa##bbb##ccc##ddd##eee##fff
#define CONCAT7_(aaa,bbb,ccc,ddd,eee,fff,ggg)  aaa##bbb##ccc##ddd##eee##fff##ggg
#define CONCAT(xxx,yyy)  CONCAT_(xxx,yyy)
#define CONCAT3(aaa,bbb,ccc)  CONCAT3_(aaa,bbb,ccc)
#define CONCAT4(aaa,bbb,ccc,ddd)  CONCAT4_(aaa,bbb,ccc,ddd)
#define CONCAT5(aaa,bbb,ccc,ddd,eee)  CONCAT5_(aaa,bbb,ccc,ddd,eee)
#define CONCAT6(aaa,bbb,ccc,ddd,eee,fff)  CONCAT6_(aaa,bbb,ccc,ddd,eee,fff)
#define CONCAT7(aaa,bbb,ccc,ddd,eee,fff,ggg)  CONCAT7_(aaa,bbb,ccc,ddd,eee,fff,ggg)
%% puts("#define CONCAT_(xxx,yyy)  xxx##yyy");
%% puts("#define CONCAT3_(aaa,bbb,ccc)  aaa##bbb##ccc");
%% #if notused
%% puts("#define CONCAT4_(aaa,bbb,ccc,ddd)  aaa##bbb##ccc##ddd");
%% puts("#define CONCAT5_(aaa,bbb,ccc,ddd,eee)  aaa##bbb##ccc##ddd##eee");
%% #endif
%% puts("#define CONCAT(xxx,yyy)  CONCAT_(xxx,yyy)");
%% puts("#define CONCAT3(aaa,bbb,ccc)  CONCAT3_(aaa,bbb,ccc)");
%% #if notused
%% puts("#define CONCAT4(aaa,bbb,ccc,ddd)  CONCAT4_(aaa,bbb,ccc,ddd)");
%% puts("#define CONCAT5(aaa,bbb,ccc,ddd,eee)  CONCAT5_(aaa,bbb,ccc,ddd,eee)");
%% #endif

# Generation of goto-tag macros:
# GENTAG(end)  ==>  end116
# This allows a macro defining marks to be used more than once per function
# but still only once per source-line.
#define GENTAG(xxx)  CONCAT(xxx,__LINE__)

# Converting tokens to strings:
# STRING(token)  ==>  "token"
#define STRING(token) #token
#define STRINGIFY(token) STRING(token)
%% puts("#define STRING(token) #token");
%% puts("#define STRINGIFY(token) STRING(token)");

# Storage-Class-Specifier in top-level-declarations:
# for variables:
#   global           globally visible variable
#   local            variable that is only visible in the file (local)
#   extern           pointer to a variable that's defined externally
# for functions:
#   global           globally visible function
#   local            function that is only visible in the file (local)
#   extern           pointer to a function that's defined externally
#   extern_C         pointer to a c-function that's defined externally
#   nonreturning     function that will never return
#   maygc            function that can trigger GC
#define global
#define local  static
# #define extern extern
#ifdef __cplusplus
  #define extern_C  extern "C"
#else
  #define extern_C  extern
#endif

/* Declaration of a function that will never return (nonreturning function)
 nonreturning_function(extern,abort,(void)); == extern void abort (void); */
#if defined(GNU) && !(__APPLE_CC__ > 1)
  #if (__GNUC__ >= 3) || ((__GNUC__ == 2) && (__GNUC_MINOR__ >= 7))
    /* Note:
       storclass __attribute__((__noreturn__)) void funname arguments
         works in gcc 2.95 or newer, and in g++ 2.7.2 or newer.
       storclass void __attribute__((__noreturn__)) funname arguments
         works in gcc 2.7.2 or newer and in g++ 2.7.2 or newer.
       storclass void funname arguments __attribute__((__noreturn__))
         works in gcc 2.7.2 or newer and in g++ 2.7.2 or newer, but
         only when followed by a semicolon, not in a function definition. */
    #define nonreturning_function(storclass,funname,arguments)  \
      storclass void __attribute__((__noreturn__)) funname arguments
  #else
    #define nonreturning_function(storclass,funname,arguments)  \
      storclass void funname arguments
  #endif
#elif defined(MICROSOFT)
  #define nonreturning_function(storclass,funname,arguments)  \
    __declspec(noreturn) storclass void funname arguments
#else
  #define nonreturning_function(storclass,funname,arguments)  \
    storclass void funname arguments
#endif
%% export_def(nonreturning_function(storclass,funname,arguments));

/* A function that can trigger GC is declared either as
   - maygc, if (1) all callers must assume the worst case: that it triggers GC,
            and (2) the function uses only the 'object's passed as arguments and
            on the STACK, but no objects stored in other non-GCsafe locations.
   - ／＊maygc＊／ otherwise. If (1) is not fulfilled, the functions begins
                   with an appropriate GCTRIGGER_IF statement. If (2) is not
                   fulfilled, the GCTRIGGER call needs to mention all other
                   non-GCsafe locations whose values are used by the function,
                   such as 'value1' or 'mv_space'. */
#define maygc

/* Storage-Class-Specifier in declarations at the beginning of a block:
 var                       will lead a variable declaration
 used by utils/varbrace to allow declarations mixed with other statements */
#define var

/* Generalized if-statement:
 if (cond1) ... {elif (condi) ...} [else ...]
 _deprecated_, please use standard C */
#define elif  else if

/* Infinite loop, can only be left with break;  or  return...;:
 _deprecated_, please use standard C */
#define loop  while (1)

/* Reverted stop condition in loops:
 Allows   until (expression) statement
 and      do statement until (expression);
 _deprecated_, please use standard C */
#define until(expression)  while(!(expression))

# Case-statement for a value >=0
# switchu (expression) ...
#ifdef GNU # for better optimization
  #define switchu(expression)  switch ((unsigned int)(expression))
#else
  #define switchu  switch
#endif

# Ignore C++ keyword.
#define export export_sym

# Swap the contents of two variables:  swap(register int, x1, x2);
#define swap(swap_type,swap_var1,swap_var2)  \
  do { var swap_type swap_temp;                                          \
    swap_temp = swap_var1; swap_var1 = swap_var2; swap_var2 = swap_temp; \
  } while(0)

# Marking a program line that may not be reached: NOTREACHED;
#define NOTREACHED  fehler_notreached(__FILE__,__LINE__)
%% puts("#define NOTREACHED  fehler_notreached(__FILE__,__LINE__)");

# Asserting an arithmetic expression: ASSERT(expr);
#define ASSERT(expr)  do { if (!(expr)) NOTREACHED; } while(0)
%%  puts("#define ASSERT(expr)  do { if (!(expr)) NOTREACHED; } while(0)");

# Ensure the Linux headers define nonstandard symbols like IPC_INFO.
#ifdef UNIX_LINUX
  #define _GNU_SOURCE 1
#endif

# alloca()
#ifdef GNU
  #define alloca  __builtin_alloca
#elif defined(MICROSOFT)
  #include <malloc.h>
  #define alloca _alloca
#elif defined(HAVE_ALLOCA_H)
  #include <alloca.h>
  #ifndef alloca # Manche definieren 'alloca' als Macro...
    #if !(defined(UNIX_OSF) || defined(UNIX_DEC_ULTRIX))
      # OSF/1 V3 declares `alloca' as returning char*, but in OSF/1 V4
      # it returns void*. I don't know how to distinguish the two.
      extern_C void* alloca (int size); # see MALLOC(3V)
    #endif
  #endif
#elif defined(_AIX)
  #pragma alloca /* AIX requires this to be the first thing in the file. */
#elif defined(BORLAND)
  #include <malloc.h> # defines  'alloca' as macro
#elif !defined(NO_ALLOCA)
  extern_C void* alloca (int size); # see MALLOC(3V)
#endif
%% #ifdef GNU
%%   emit_define("alloca","__builtin_alloca");
%% #elif defined(MICROSOFT)
%%   puts("#include <malloc.h>");
%%   emit_define("alloca","_alloca");
%% #elif defined(HAVE_ALLOCA_H)
%%   puts("#include <alloca.h>");
%%   #ifndef alloca
%%     #if !(defined(UNIX_OSF) || defined(UNIX_DEC_ULTRIX))
%%       puts("extern void* alloca (int size);");
%%     #endif
%%   #endif
%% #elif defined(_AIX)
%%   puts("#pragma alloca");
%% #elif !defined(NO_ALLOCA)
%%   puts("extern void* alloca (int size);");
%% #endif

#define MALLOC(size,type)   (type*)malloc((size)*sizeof(type))

# Literal constants of 64-bit integer types
# LL(nnnn)  = nnnn parsed as a sint64
# ULL(nnnn) = nnnn parsed as a uint64
#if defined(HAVE_LONGLONG)
  #define LL(nnnn) nnnn##LL
  #define ULL(nnnn) nnnn##ULL
#elif defined(MICROSOFT)
  #define LL(nnnn) nnnn##i64
  #define ULL(nnnn) nnnn##ui64
#endif
%% #if defined(HAVE_LONGLONG)
%%   puts("#define LL(nnnn) nnnn##LL");
%%   puts("#define ULL(nnnn) nnnn##ULL");
%% #elif defined(MICROSOFT)
%%   puts("#define LL(nnnn) nnnn##i64");
%%   puts("#define ULL(nnnn) nnnn##ui64");
%% #endif

# Synonyms for Byte, Word, Longword:
# SBYTE   = signed 8 bit integer
# UBYTE   = unsigned 8 bit int
# SWORD   = signed 16 bit int
# UWORD   = unsigned 16 bit int
# SLONG   = signed 32 bit int
# ULONG   = unsigned 32 bit int
# On the other hand, "char" is only used as an element of a string
# You never really compute with a "char"; it might depend on
# __CHAR_UNSIGNED___!
#if (char_bitsize==8)
  #ifdef __CHAR_UNSIGNED__
    typedef signed char  SBYTE;
  #else
    typedef char         SBYTE;
  #endif
  typedef unsigned char  UBYTE;
#else
  #error "No 8 bit integer type? -- Which Interger-type has 8 Bit?"
#endif
#if (short_bitsize==16)
  typedef short          SWORD;
  typedef unsigned short UWORD;
#else
  #error "No 16 bit integer type? -- Which Integer-type has 16 Bit?"
#endif
#if (long_bitsize==32)
  typedef long           SLONG;
  typedef unsigned long  ULONG;
#elif (int_bitsize==32)
  typedef int            SLONG;
  typedef unsigned int   ULONG;
#else
  #error "No 32 bit integer type? -- Which Integer-type has 32 Bit?"
#endif
#if (long_bitsize==64)
  typedef long           SLONGLONG;
  typedef unsigned long  ULONGLONG;
  #undef HAVE_LONGLONG
  #define HAVE_LONGLONG
#elif defined(MICROSOFT)
  typedef __int64           SLONGLONG;
  typedef unsigned __int64  ULONGLONG;
  #define HAVE_LONGLONG
#elif defined(HAVE_LONGLONG)
 #if defined(long_long_bitsize) && (long_long_bitsize==64)
  typedef long long           SLONGLONG;
  typedef unsigned long long  ULONGLONG;
 #else # useless type
  #undef HAVE_LONGLONG
 #endif
#endif
#if defined(WIDE) && !defined(HAVE_LONGLONG)
  #error "No 64 bit integer type? -- Which Integer-type has 64 Bit?"
#endif
%% #ifdef __CHAR_UNSIGNED__
%%   emit_typedef("signed char","SBYTE");
%% #else
%%   emit_typedef("char","SBYTE");
%% #endif
%% emit_typedef("unsigned char","UBYTE");
%% emit_typedef("short","SWORD");
%% emit_typedef("unsigned short","UWORD");
%% #if (long_bitsize==32)
%%   emit_typedef("long","SLONG");
%%   emit_typedef("unsigned long","ULONG");
%% #elif (int_bitsize==32)
%%   emit_typedef("int","SLONG");
%%   emit_typedef("unsigned int","ULONG");
%% #endif
%% #if (long_bitsize==64)
%%   emit_typedef("long","SLONGLONG");
%%   emit_typedef("unsigned long","ULONGLONG");
%% #elif defined(MICROSOFT)
%%   emit_typedef("__int64","SLONGLONG");
%%   emit_typedef("unsigned __int64","ULONGLONG");
%% #elif defined(HAVE_LONGLONG)
%%   emit_typedef("long long","SLONGLONG");
%%   emit_typedef("unsigned long long","ULONGLONG");
%% #endif

#include <stdbool.h>  /* boolean values */
%% #ifdef HAVE_STDBOOL_H
%%   puts("#include <stdbool.h>");
%% #else
%%   print_file("stdbool.h");
%% #endif

# Type for signed values, results of comparisons, tertiary enums
# with values +1, 0, -1
typedef signed int  signean;
#define signean_plus    1 # +1
#define signean_null    0 #  0
#define signean_minus  -1 # -1

# Null pointers
#ifdef __cplusplus
  #define NULL  0
#elif !(defined(INTEL) || defined(_AIX))
  #define NULL  ((void*) 0L)
#endif
%% puts("#undef NULL");
%% export_def(NULL);

#include <stdio.h>    /* libc i/o */

# A more precise classification of the operating system:
# (This test works only after at least one system header has been included.)
#if (__GLIBC__ >= 2)
  #define UNIX_GNU # glibc2 (may be UNIX_LINUX, UNIX_HURD or UNIX_FREEBSD)
#endif

# Determine the offset of a component 'ident' in a struct of the type 'type':
# See 0 as pointer to 'type', put a struct 'type' there and determine the
# address of its component 'ident' and return it as number:
#if defined(HAVE_OFFSETOF) || defined(__MINGW32__) || (defined(BORLAND) && defined(WIN32)) || defined(MICROSOFT)
  #include <stddef.h>
#else
  #undef offsetof
  #define offsetof(type,ident)  ((ULONG)&(((type*)0)->ident))
#endif
# Determine the offset of an array 'ident' in a struct of the type 'type':
#if defined(__cplusplus) || defined(MICROSOFT)
  #define offsetofa(type,ident)  offsetof(type,ident)
#else
  #define offsetofa(type,ident)  offsetof(type,ident[0])
#endif

# alignof(type) is a constant expression, returning the alignment of type.
#ifdef __cplusplus
  #ifdef GNU
    #define alignof(type)  __alignof__(type)
  #else
    template <class type> struct alignof_helper { char slot1; type slot2; };
    #define alignof(type)  offsetof(alignof_helper<type>, slot2)
  #endif
#else
  #define alignof(type)  offsetof(struct { char slot1; type slot2; }, slot2)
#endif

# Unspecified length of arrays in structures:
# struct { ...; ...; type x[unspecified]; }
# Instead of sizeof(..) you'll always have to use offsetof(..,x).
#if defined(GNU) # GNU-C is able to work with arrays of length 0
  #define unspecified 0
#elif 0
  # Usually one would omit the array's limit
  #define unspecified
#else
  # However, HP-UX- and IRIX-compilers will only work with this:
  #define unspecified 1
#endif
%% export_def(unspecified);

# Pointer arithmetics: add a given offset (measured in bytes)
# to a pointer.
#if defined(GNU) || (pointer_bitsize > 32)
 /* Essential for GNU-C for initialization of static-variables
   (must be a bug in 'c-typeck.c' in 'initializer_constant_valid_p'):
   The only correct way, if sizeof(ULONG) < sizeof(void*): */
  #define pointerplus(pointer,offset)  ((UBYTE*)(pointer)+(offset))
#else
  # Cheap way:
  #define pointerplus(pointer,offset)  ((void*)((ULONG)(pointer)+(offset)))
#endif
%% export_def(pointerplus(pointer,offset));

# Bit number n (0<=n<32)
#define bit(n)  (1L<<(n))
# Bit number n (0<n<=32) mod 2^32
#define bitm(n)  (2L<<((n)-1))
# Bit-test of bit n in x, n constant, x an oint:
#if !defined(SPARC)
  #define bit_test(x,n)  ((x) & bit(n))
#else
  # On SPARC-processors, long constants are slower than shifts.
  #if defined(SPARC64)
    #if !defined(GNU)
      #define bit_test(x,n)  \
        ((n)<12 ? ((x) & bit(n)) : ((sint64)((uint64)(x) << (63-(n))) < 0))
    #else # the GNU-compiler will optimize boolean expressions better this way:
      #define bit_test(x,n)  \
        (   ( ((n)<12) && ((x) & bit(n)) )                           \
         || ( ((n)>=12) && ((sint64)((uint64)(x) << (63-(n))) < 0) ) \
        )
    #endif
  #else
    #if !defined(GNU)
      #define bit_test(x,n)  \
        ((n)<12 ? ((x) & bit(n)) : ((sint32)((uint32)(x) << (31-(n))) < 0))
    #else # the GNU-compiler will optimize boolean expressions better this way:
      #define bit_test(x,n)  \
        (   ( ((n)<12) && ((x) & bit(n)) )                           \
         || ( ((n)>=12) && ((sint32)((uint32)(x) << (31-(n))) < 0) ) \
        )
    #endif
  #endif
#endif
# Minus bit number n (0<=n<32)
#define minus_bit(n)  (-1L<<(n))
# Minus bit number n (0<n<=32) mod 2^32
#define minus_bitm(n)  (-2L<<((n)-1))
%% export_def(bit(n));
%% #if notused
%% export_def(bitm(n));
%% #endif
%% export_def(bit_test(x,n));
%% export_def(minus_bit(n));
%% #if notused
%% export_def(minus_bitm(n));
%% #endif

# floor(a,b) yields for a>=0, b>0  floor(a/b).
# b should be a 'constant expression'.
#define floor(a_from_floor,b_from_floor)  ((a_from_floor) / (b_from_floor))
%% /* FIXME: Difference between lispbibl.d and clisp.h */
%% puts("#define ifloor(a_from_floor,b_from_floor)  ((a_from_floor) / (b_from_floor))");

# ceiling(a,b) yields for a>=0, b>0  ceiling(a/b) = floor((a+b-1)/b).
# b should be a 'constant expression'.
#define ceiling(a_from_ceiling,b_from_ceiling)  \
  (((a_from_ceiling) + (b_from_ceiling) - 1) / (b_from_ceiling))
%% export_def(ceiling(a_from_ceiling,b_from_ceiling));

# round_down(a,b) rounds a>=0 so that b>0 divides it.
# b should be a 'constant expression'.
#define round_down(a_from_round,b_from_round)  \
  (floor(a_from_round,b_from_round)*(b_from_round))
%% /* FIXME: Difference between lispbibl.d and clisp.h */
%% puts("#define round_down(a_from_round,b_from_round)  (ifloor(a_from_round,b_from_round)*(b_from_round))");

# round_up(a,b) rounds a>=0 so that b>0 divides it.
# b should be a 'constant expression'.
#define round_up(a_from_round,b_from_round)  \
  (ceiling(a_from_round,b_from_round)*(b_from_round))
%% export_def(round_up(a_from_round,b_from_round));

# non-local exits
#include <setjmp.h>
#if defined(UNIX) && defined(HAVE__JMP) && !defined(UNIX_LINUX) && !defined(UNIX_GNU) && !defined(UNIX_BEOS) && !defined(UNIX_CYGWIN32)
  # The following routines are more efficient (don't use with signal-masks):
  #undef setjmp
  #undef longjmp
  #define setjmp  _setjmp
  #define longjmp  _longjmp
  #ifdef LONGJMP_RETURNS
    # _longjmp(jmpbuf,value) can return if jmpbuf is invalid.
    #undef longjmp
    #define longjmp(x,y)  (_longjmp(x,y), NOTREACHED)
  #endif
#endif
# A longjmp() can only be called using an `int'.
# But if we want to use a `long' and if sizeof(int) < sizeof(long),
# we'll need a global variable:
#if (int_bitsize == long_bitsize)
  #define setjmpl(x)  setjmp(x)
  #define longjmpl(x,y)  longjmp(x,y)
#else # (int_bitsize < long_bitsize)
  extern long jmpl_value;
  #define setjmpl(x)  (setjmp(x) ? jmpl_value : 0)
  #define longjmpl(x,y)  (jmpl_value = (y), longjmp(x,1))
#endif

# An alloca() replacement, used for DYNAMIC_ARRAY and SAVE_NUM_STACK.
# See spvw_alloca.d.
#if !(defined(GNU) || (defined(UNIX) && !defined(NO_ALLOCA) && !defined(SPARC)) || defined(BORLAND) || defined(MICROSOFT))
  #define NEED_MALLOCA
  #include <stdlib.h>
  extern void* malloca (size_t size);
  extern void freea (void* ptr);
#endif

# Dynamically allocated array with dynamic extent:
# Example:
#     var DYNAMIC_ARRAY(my_array,uintL,n);
#     ...
#     FREE_DYNAMIC_ARRAY(my_array);
# Attention: depending on your implementation my_array is either the array
# itself or a pointer to the array! Always use my_array only as expression!
#if defined(GNU)
  # can deal with dynamically allocated arrays in the machine stack
  # { var uintL my_array[n]; ... }
  #define DYNAMIC_ARRAY(arrayvar,arrayeltype,arraysize)  \
    arrayeltype arrayvar[arraysize]
  #define FREE_DYNAMIC_ARRAY(arrayvar)
  #ifdef DECALPHA # GCC 2.5.5 Bug umgehen
    #undef DYNAMIC_ARRAY
    #define DYNAMIC_ARRAY(arrayvar,arrayeltype,arraysize)  \
      arrayeltype arrayvar[(arraysize)+1]
  #endif
#elif (defined(UNIX) && (defined(HAVE_ALLOCA_H) || defined(_AIX) || !defined(NO_ALLOCA))) || defined(BORLAND) || defined(MICROSOFT)
  # Allocate space in machine stack.
  # { var uintL* my_array = (uintL*)alloca(n*sizeof(uintL)); ... }
  #define DYNAMIC_ARRAY(arrayvar,arrayeltype,arraysize)  \
    arrayeltype* arrayvar = (arrayeltype*)alloca((arraysize)*sizeof(arrayeltype))
  #define FREE_DYNAMIC_ARRAY(arrayvar)
  # no error check??
#else
  # Allocate space somewhere else and then free it.
  # { var uintL* my_array = (uintL*)malloc(n*sizeof(uintL)); ... free(my_array); }
  #define DYNAMIC_ARRAY(arrayvar,arrayeltype,arraysize)  \
    arrayeltype* arrayvar = (arrayeltype*)malloca((arraysize)*sizeof(arrayeltype))
  #define FREE_DYNAMIC_ARRAY(arrayvar)  freea(arrayvar)
#endif
%% export_def(DYNAMIC_ARRAY(arrayvar,arrayeltype,arraysize));
%% export_def(FREE_DYNAMIC_ARRAY(arrayvar));

# Signed/Unsigned-Integer-types with given minumum size:
typedef UBYTE   uint1;   # unsigned 1 bit Integer
typedef SBYTE   sint1;   # signed 1 bit Integer
typedef UBYTE   uint2;   # unsigned 2 bit Integer
typedef SBYTE   sint2;   # signed 2 bit Integer
typedef UBYTE   uint3;   # unsigned 3 bit Integer
typedef SBYTE   sint3;   # signed 3 bit Integer
typedef UBYTE   uint4;   # unsigned 4 bit Integer
typedef SBYTE   sint4;   # signed 4 bit Integer
typedef UBYTE   uint5;   # unsigned 5 bit Integer
typedef SBYTE   sint5;   # signed 5 bit Integer
typedef UBYTE   uint6;   # unsigned 6 bit Integer
typedef SBYTE   sint6;   # signed 6 bit Integer
typedef UBYTE   uint7;   # unsigned 7 bit Integer
typedef SBYTE   sint7;   # signed 7 bit Integer
typedef UBYTE   uint8;   # unsigned 8 bit Integer
typedef SBYTE   sint8;   # signed 8 bit Integer
typedef UWORD   uint9;   # unsigned 9 bit Integer
typedef SWORD   sint9;   # signed 9 bit Integer
typedef UWORD   uint10;  # unsigned 10 bit Integer
typedef SWORD   sint10;  # signed 10 bit Integer
typedef UWORD   uint11;  # unsigned 11 bit Integer
typedef SWORD   sint11;  # signed 11 bit Integer
typedef UWORD   uint12;  # unsigned 12 bit Integer
typedef SWORD   sint12;  # signed 12 bit Integer
typedef UWORD   uint13;  # unsigned 13 bit Integer
typedef SWORD   sint13;  # signed 13 bit Integer
typedef UWORD   uint14;  # unsigned 14 bit Integer
typedef SWORD   sint14;  # signed 14 bit Integer
typedef UWORD   uint15;  # unsigned 15 bit Integer
typedef SWORD   sint15;  # signed 15 bit Integer
typedef UWORD   uint16;  # unsigned 16 bit Integer
typedef SWORD   sint16;  # signed 16 bit Integer
typedef ULONG   uint17;  # unsigned 17 bit Integer
typedef SLONG   sint17;  # signed 17 bit Integer
typedef ULONG   uint18;  # unsigned 18 bit Integer
typedef SLONG   sint18;  # signed 18 bit Integer
typedef ULONG   uint19;  # unsigned 19 bit Integer
typedef SLONG   sint19;  # signed 19 bit Integer
typedef ULONG   uint20;  # unsigned 20 bit Integer
typedef SLONG   sint20;  # signed 20 bit Integer
typedef ULONG   uint21;  # unsigned 21 bit Integer
typedef SLONG   sint21;  # signed 21 bit Integer
typedef ULONG   uint22;  # unsigned 22 bit Integer
typedef SLONG   sint22;  # signed 22 bit Integer
typedef ULONG   uint23;  # unsigned 23 bit Integer
typedef SLONG   sint23;  # signed 23 bit Integer
typedef ULONG   uint24;  # unsigned 24 bit Integer
typedef SLONG   sint24;  # signed 24 bit Integer
typedef ULONG   uint25;  # unsigned 25 bit Integer
typedef SLONG   sint25;  # signed 25 bit Integer
typedef ULONG   uint26;  # unsigned 26 bit Integer
typedef SLONG   sint26;  # signed 26 bit Integer
typedef ULONG   uint27;  # unsigned 27 bit Integer
typedef SLONG   sint27;  # signed 27 bit Integer
typedef ULONG   uint28;  # unsigned 28 bit Integer
typedef SLONG   sint28;  # signed 28 bit Integer
typedef ULONG   uint29;  # unsigned 29 bit Integer
typedef SLONG   sint29;  # signed 29 bit Integer
typedef ULONG   uint30;  # unsigned 30 bit Integer
typedef SLONG   sint30;  # signed 30 bit Integer
typedef ULONG   uint31;  # unsigned 31 bit Integer
typedef SLONG   sint31;  # signed 31 bit Integer
typedef ULONG   uint32;  # unsigned 32 bit Integer
typedef SLONG   sint32;  # signed 32 bit Integer
#ifdef HAVE_LONGLONG
  typedef ULONGLONG  uint33;  # unsigned 33 bit Integer
  typedef SLONGLONG  sint33;  # signed 33 bit Integer
  typedef ULONGLONG  uint48;  # unsigned 48 bit Integer
  typedef SLONGLONG  sint48;  # signed 48 bit Integer
  typedef ULONGLONG  uint64;  # unsigned 64 bit Integer
  typedef SLONGLONG  sint64;  # signed 64 bit Integer
#endif
#define exact_uint_size_p(n) (((n)==char_bitsize)||((n)==short_bitsize)||((n)==int_bitsize)||((n)==long_bitsize))
#define signed_int_with_n_bits(n) CONCAT(sint,n)
#define unsigned_int_with_n_bits(n) CONCAT(uint,n)
# Use 'uintn' and 'sintn' for Integers with exactly specified width.
# exact_uint_size_p(n) specifies, whether the uint with n Bits has really
# only n Bits.
%% { int i;
%%   for (i=1; i<=8; i++) {
%%     sprintf(buf,"uint%d",i); emit_typedef("UBYTE",buf);
%%     sprintf(buf,"sint%d",i); emit_typedef("SBYTE",buf);
%%   }
%%   for (i=9; i<=16; i++) {
%%     sprintf(buf,"uint%d",i); emit_typedef("UWORD",buf);
%%     sprintf(buf,"sint%d",i); emit_typedef("SWORD",buf);
%%   }
%%   for (i=17; i<=32; i++) {
%%     sprintf(buf,"uint%d",i); emit_typedef("ULONG",buf);
%%     sprintf(buf,"sint%d",i); emit_typedef("SLONG",buf);
%%   }
%%   #ifdef HAVE_LONGLONG
%%     for (i=33; i<=64; i++)
%%       if ((i==33) || (i==48) || (i==64)) {
%%         sprintf(buf,"uint%d",i); emit_typedef("ULONGLONG",buf);
%%         sprintf(buf,"sint%d",i); emit_typedef("SLONGLONG",buf);
%%       }
%%   #endif
%% }

# 'uintX' and 'sintX' mean unsigned bzw. signed integer - types with
# wordsize X (X=B,W,L,Q) here as well.
#define intBsize 8
  typedef signed_int_with_n_bits(intBsize)    sintB;
  typedef unsigned_int_with_n_bits(intBsize)  uintB;
#define intWsize 16
  typedef signed_int_with_n_bits(intWsize)    sintW;
  typedef unsigned_int_with_n_bits(intWsize)  uintW;
#define intLsize 32
  typedef signed_int_with_n_bits(intLsize)    sintL;
  typedef unsigned_int_with_n_bits(intLsize)  uintL;
#if defined(DECALPHA) || defined(MIPS64) || defined(SPARC64) || defined(IA64) || defined(AMD64)
  # Machine has real 64-bit integers in hardware.
  #define intQsize 64
  typedef signed_int_with_n_bits(intQsize)    sintQ;
  typedef unsigned_int_with_n_bits(intQsize)  uintQ;
  typedef sintQ  sintL2;
  typedef uintQ  uintL2;
#else
  # Emulate 64-Bit-numbers using two 32-Bit-numbers.
  typedef struct { sintL hi; uintL lo; } sintL2; # signed 64 Bit integer
  typedef struct { uintL hi; uintL lo; } uintL2; # unsigned 64 Bit integer
#endif
# Use 'uintX' and 'sintX' for Integers with approximately given width
# and a minumum of storage space.
%% sprintf(buf,"sint%d",intBsize); emit_typedef(buf,"sintB");
%% sprintf(buf,"uint%d",intBsize); emit_typedef(buf,"uintB");
%% #if notused
%% sprintf(buf,"sint%d",intWsize); emit_typedef(buf,"sintW");
%% #endif
%% sprintf(buf,"uint%d",intWsize); emit_typedef(buf,"uintW");
%% sprintf(buf,"sint%d",intLsize); emit_typedef(buf,"sintL");
%% sprintf(buf,"uint%d",intLsize); emit_typedef(buf,"uintL");
%% #if notused
%% #ifdef intQsize
%%   sprintf(buf,"sint%d",intQsize); emit_typedef(buf,"sintQ");
%%   sprintf(buf,"uint%d",intQsize); emit_typedef(buf,"uintQ");
%% #else
%%   emit_typedef("struct { sintL hi; uintL lo; }","sintL2");
%%   emit_typedef("struct { uintL hi; uintL lo; }","uintL2");
%% #endif
%% #endif

# From here on 'uintP' and 'sintP' are unsigned or signed integer types,
# which are as wide as void* - pointers
typedef signed_int_with_n_bits(pointer_bitsize)    sintP;
typedef unsigned_int_with_n_bits(pointer_bitsize)  uintP;
%% sprintf(buf,"sint%d",pointer_bitsize); emit_typedef(buf,"sintP");
%% sprintf(buf,"uint%d",pointer_bitsize); emit_typedef(buf,"uintP");

# From here on 'uintXY' and 'sintXY' mean unsigned or signed integer types,
# with word sizes X or Y (X,Y=B,W,L).
#if (defined(MC680X0) && !defined(HPUX_ASSEMBLER)) || defined(VAX)
  # The 68000 offers good processing of uintB and uintW, especially
  # DBRA-commands for uintW.
  #define intBWsize intBsize
  #define intWLsize intWsize
  #define intBWLsize intBsize
#elif (defined(MC680X0) && defined(HPUX_ASSEMBLER)) || defined(SPARC) || defined(HPPA) || defined(MIPS) || defined(M88000) || defined(POWERPC) || defined(S390)
  # The Sparc-processor computes rather badly with uintB and uintW.
  # Other 32-Bit-processoren have similar weaknesses.
  #define intBWsize intWsize
  #define intWLsize intLsize
  #define intBWLsize intLsize
#elif defined(I80386) || defined(AMD64)
  # If you compute using uintB and uintW on a 80386, there will be many
  # Zero-Extends, that will - because there aren't enough registers - load
  # other variables into memory, which is rather unnecessary.
  #define intBWsize intWsize
  #define intWLsize intLsize
  #define intBWLsize intLsize
#elif defined(ARM)
  # The ARM computes very badly when it uses uintB and uintW.
  #define intBWsize intBsize
  #define intWLsize intLsize
  #define intBWLsize intLsize
#elif defined(DECALPHA) || defined(IA64)
  # 64-bit processors also compute badly with uintB and uintW.
  #define intBWsize intWsize
  #define intWLsize intLsize
  #define intBWLsize intLsize
#else
  #error "Preferred integer sizes depend on CPU -- readjust intBWsize, intWLsize, intBWLsize!"
#endif
typedef signed_int_with_n_bits(intBWsize)     sintBW;
typedef unsigned_int_with_n_bits(intBWsize)   uintBW;
typedef signed_int_with_n_bits(intWLsize)     sintWL;
typedef unsigned_int_with_n_bits(intWLsize)   uintWL;
typedef signed_int_with_n_bits(intBWLsize)    sintBWL;
typedef unsigned_int_with_n_bits(intBWLsize)  uintBWL;
# Use 'uintXY' and 'sintXY' for integers with given minumum width,
# that allow easy computations.
%% #if notused
%% sprintf(buf,"sint%d",intBWsize); emit_typedef(buf,"sintBW");
%% sprintf(buf,"uint%d",intBWsize); emit_typedef(buf,"uintBW");
%% sprintf(buf,"sint%d",intWLsize); emit_typedef(buf,"sintWL");
%% #endif
%% sprintf(buf,"uint%d",intWLsize); emit_typedef(buf,"uintWL");
%% #if notused
%% sprintf(buf,"sint%d",intBWLsize); emit_typedef(buf,"sintBWL");
%% #endif
%% sprintf(buf,"uint%d",intBWLsize); emit_typedef(buf,"uintBWL");

# Loop that will excute as statement a certain number of times:
# dotimesW(countvar,count,statement);  if count fits into a uintW,
# dotimesL(countvar,count,statement);  if  count only fits into a uintL,
# dotimesV(countvar,count,statement);  if  count only fits into a uintV,
# dotimespW(countvar,count,statement);  if count fits into a uintW and is >0,
# dotimespL(countvar,count,statement);  if count fits only into a uintL and is >0.
# dotimespV(countvar,count,statement);  if count fits only into a uintV and is >0.
# The variable countvar has to be declared previously, be of type uintW or uintL,
# and will be changed by this expression.
# It must not be used in the statement itself!
# The expression count will only be evaluated once (at the beginning).
#if defined(GNU) && defined(MC680X0) && !defined(HPUX_ASSEMBLER)
  # GNU-C on a 680X0 can be persuaded to use the DBRA-instruction:
  #define fast_dotimesW
  # To find out, what the best was to 'persuade' GNU-C is, check the
  # code, that'll be generated for spvw.d:gc_markphase().
  # Or a small test program (dbratest.c), that is compiled with
  # "gcc -O6 -da -S dbratest.c", and take a look at dbratest.s
  # and dbratest.c.flow as well as dbratest.c.combine.
  #if (__GNUC__<2) # GNU C Version 1
    #define dotimesW_(countvar_from_dotimesW,count_from_dotimesW,statement_from_dotimesW)  \
      { countvar_from_dotimesW = (count_from_dotimesW);     \
        if (!(countvar_from_dotimesW==0))                   \
          { countvar_from_dotimesW--;                       \
            do {statement_from_dotimesW}                    \
               until ((sintW)--countvar_from_dotimesW==-1); \
      }   }
    #define dotimespW_(countvar_from_dotimespW,count_from_dotimespW,statement_from_dotimespW)  \
      { countvar_from_dotimespW = (count_from_dotimespW)-1;                         \
        do {statement_from_dotimespW} until ((sintW)--countvar_from_dotimespW==-1); \
      }
  #else
    #define dotimesW_(countvar_from_dotimesW,count_from_dotimesW,statement_from_dotimesW)  \
      { countvar_from_dotimesW = (count_from_dotimesW);        \
        if (!(countvar_from_dotimesW==0))                      \
          { countvar_from_dotimesW--;                          \
            do {statement_from_dotimesW}                       \
               until ((sintW)(--countvar_from_dotimesW)+1==0); \
      }   }
    #define dotimespW_(countvar_from_dotimespW,count_from_dotimespW,statement_from_dotimespW)  \
      { countvar_from_dotimespW = (count_from_dotimespW)-1;                            \
        do {statement_from_dotimespW} until ((sintW)(--countvar_from_dotimespW)+1==0); \
      }
  #endif
#else
  #define dotimesW_(countvar_from_dotimesW,count_from_dotimesW,statement_from_dotimesW)  \
    { countvar_from_dotimesW = (count_from_dotimesW);         \
      until (countvar_from_dotimesW==0)                       \
        {statement_from_dotimesW; countvar_from_dotimesW--; } \
    }
  #define dotimespW_(countvar_from_dotimespW,count_from_dotimespW,statement_from_dotimespW)  \
    { countvar_from_dotimespW = (count_from_dotimespW);                   \
      do {statement_from_dotimespW} until (--countvar_from_dotimespW==0); \
    }
#endif
#if defined(GNU) && defined(MC680X0) && !defined(HPUX_ASSEMBLER)
  # GNU-C on a 680X0 can be 'persuaded' to use the DBRA-instruction
  # in an intelligent manner:
  #define fast_dotimesL
  #define dotimesL_(countvar_from_dotimesL,count_from_dotimesL,statement_from_dotimesL)  \
    { countvar_from_dotimesL = (count_from_dotimesL);           \
      if (!(countvar_from_dotimesL==0))                         \
        { countvar_from_dotimesL--;                             \
          do {statement_from_dotimesL}                          \
             until ((sintL)(--countvar_from_dotimesL) == -1);   \
    }   }
  #define dotimespL_(countvar_from_dotimespL,count_from_dotimespL,statement_from_dotimespL)  \
    { countvar_from_dotimespL = (count_from_dotimespL)-1;                             \
      do {statement_from_dotimespL} until ((sintL)(--countvar_from_dotimespL) == -1); \
    }
#endif
#ifndef dotimesL_
  #define dotimesL_(countvar_from_dotimesL,count_from_dotimesL,statement_from_dotimesL)  \
    { countvar_from_dotimesL = (count_from_dotimesL);         \
      until (countvar_from_dotimesL==0)                       \
        {statement_from_dotimesL; countvar_from_dotimesL--; } \
    }
  #define dotimespL_(countvar_from_dotimespL,count_from_dotimespL,statement_from_dotimespL)  \
    { countvar_from_dotimespL = (count_from_dotimespL);                   \
      do {statement_from_dotimespL} until (--countvar_from_dotimespL==0); \
    }
#endif
#if defined(GNU) && defined(__OPTIMIZE__)
  # It happened twice to me that I used dotimesL on a
  # variable of type uintC. I check for that now, so that
  # Joerg and Marcus won't have to search for that anymore.
  # The GCC will optimize the dummy-call away, if things go by plan.
  # If not, you'll see a linker error.
  #define dotimes_check_sizeof(countvar,type)  \
    if (!(sizeof(countvar)==sizeof(type))) { dotimes_called_with_count_of_wrong_size(); }
  extern void dotimes_called_with_count_of_wrong_size (void); # non-existing function
#else
  #define dotimes_check_sizeof(countvar,type)
#endif
#define dotimesW(countvar_from_dotimesW,count_from_dotimesW,statement_from_dotimesW) \
  do { dotimes_check_sizeof(countvar_from_dotimesW,uintW); \
    dotimesW_(countvar_from_dotimesW,count_from_dotimesW,statement_from_dotimesW); \
  } while(0)
#define dotimespW(countvar_from_dotimespW,count_from_dotimespW,statement_from_dotimespW) \
  do { dotimes_check_sizeof(countvar_from_dotimespW,uintW); \
    dotimespW_(countvar_from_dotimespW,count_from_dotimespW,statement_from_dotimespW); \
  } while(0)
#define dotimesL(countvar_from_dotimesL,count_from_dotimesL,statement_from_dotimesL) \
  do { dotimes_check_sizeof(countvar_from_dotimesL,uintL); \
    dotimesL_(countvar_from_dotimesL,count_from_dotimesL,statement_from_dotimesL); \
  } while(0)
#define dotimespL(countvar_from_dotimespL,count_from_dotimespL,statement_from_dotimespL) \
  do { dotimes_check_sizeof(countvar_from_dotimespL,uintL); \
    dotimespL_(countvar_from_dotimespL,count_from_dotimespL,statement_from_dotimespL); \
  } while(0)
#define dotimesV(countvar_from_dotimesV,count_from_dotimesV,statement_from_dotimesV) \
  do { dotimes_check_sizeof(countvar_from_dotimesV,uintV); \
    dotimesL_(countvar_from_dotimesV,count_from_dotimesV,statement_from_dotimesV); \
  } while(0)
#define dotimespV(countvar_from_dotimespV,count_from_dotimespV,statement_from_dotimespV) \
  do { dotimes_check_sizeof(countvar_from_dotimespV,uintV); \
    dotimespL_(countvar_from_dotimespV,count_from_dotimespV,statement_from_dotimespV); \
  } while(0)
# doconsttimes(count,statement);
# executes a statement count times (count times the code!),
# where count is a constant-expression >=0, <=8.
#define doconsttimes(count_from_doconsttimes,statement_from_doconsttimes) \
 do { if (0 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (1 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (2 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (3 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (4 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (5 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (6 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
      if (7 < (count_from_doconsttimes)) { statement_from_doconsttimes; } \
 } while(0)
# DOCONSTTIMES(count,macroname);
# calls the macro macroname count times (count times the code!),
# where count is a constant-expression >=0, <=8.
# And macroname will get the values 0,...,count-1 in sequence.
#define DOCONSTTIMES(count_from_DOCONSTTIMES,macroname_from_DOCONSTTIMES)  \
 do { if (0 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((0 < (count_from_DOCONSTTIMES) ? 0 : 0)); } \
      if (1 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((1 < (count_from_DOCONSTTIMES) ? 1 : 0)); } \
      if (2 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((2 < (count_from_DOCONSTTIMES) ? 2 : 0)); } \
      if (3 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((3 < (count_from_DOCONSTTIMES) ? 3 : 0)); } \
      if (4 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((4 < (count_from_DOCONSTTIMES) ? 4 : 0)); } \
      if (5 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((5 < (count_from_DOCONSTTIMES) ? 5 : 0)); } \
      if (6 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((6 < (count_from_DOCONSTTIMES) ? 6 : 0)); } \
      if (7 < (count_from_DOCONSTTIMES)) { macroname_from_DOCONSTTIMES((7 < (count_from_DOCONSTTIMES) ? 7 : 0)); } \
 } while(0)

# From here on  uintC means an unsigned integer type, that'll allow
# easy counting. Subset relation: uintW <= uintC <= uintL.
#define intCsize intWLsize
#define uintC uintWL
#define sintC sintWL
#if (intCsize==intWsize)
  #define dotimesC dotimesW
  #define dotimespC dotimespW
#endif
#if (intCsize==intLsize)
  #define dotimesC dotimesL
  #define dotimespC dotimespL
#endif
# Use 'uintC' for counters, which are small most of the time.
%% export_def(uintC);
%% #if notused
%% export_def(sintC);
%% #endif

# The arithmetics use "digit sequences" of "digits".
# They are unsigned ints with intDsize bits (should be =8 or =16 or =32).
# If  HAVE_DD: "double-digits" are unsigned ints with 2*intDsize<=32 bits.
#if defined(MC680X0) && !defined(MC680Y0)
  #define intDsize 16
  #define intDDsize 32  # = 2*intDsize
  #define log2_intDsize  4  # = log2(intDsize)
#elif defined(MC680Y0) || defined(I80386) || defined(SPARC) || defined(HPPA) || defined(MIPS) || defined(M88000) || defined(POWERPC) || defined(VAX) || defined(ARM) || defined(DECALPHA) || defined(IA64) || defined(AMD64) || defined(S390)
  #define intDsize 32
  #define intDDsize 64  # = 2*intDsize
  #define log2_intDsize  5  # = log2(intDsize)
#else
  #error "Preferred digit size depends on CPU -- readjust intDsize!"
#endif
typedef unsigned_int_with_n_bits(intDsize)  uintD;
typedef signed_int_with_n_bits(intDsize)    sintD;
#if (intDDsize<=32) || ((intDDsize<=64) && (defined(DECALPHA) || defined(MIPS64) || defined(SPARC64) || defined(IA64) || defined(AMD64)))
  #define HAVE_DD 1
  typedef unsigned_int_with_n_bits(intDDsize)  uintDD;
  typedef signed_int_with_n_bits(intDDsize)    sintDD;
#else
  #define HAVE_DD 0
#endif
%% #if notused
%% sprintf(buf,"sint%d",intDsize); emit_typedef(buf,"sintD");
%% #endif
%% sprintf(buf,"uint%d",intDsize); emit_typedef(buf,"uintD");

/* Other acronyms like 'oint', 'tint', 'aint', 'cint' will be used
   for the corresponding integer types:
   Integer type      contains information equivalent to
      oint           LISP object
      tint           type code of a LISP object
      aint           address of a LISP object
      cint           LISP character

 Usually sizeof(oint) = sizeof(aint) = sizeof(uintL) = 32 Bit.
 Under the model WIDE sizeof(oint) is > sizeof(uintL).
 Model WIDE_HARD stands for sizeof(aint) > sizeof(uintL).
   This model is to be chosen if the following holds true:
   sizeof(void*) > sizeof(uintL) = 32 bit. It also requires that
   sizeof(long) = sizeof(void*) = 64 bit, because some 64-bit numbers
   appear as pre-processor constants.
 Model WIDE_SOFT stands for sizeof(oint) = 64 bit and sizeof(aint) = 32 bit.
   This model can be chosen on any 32-Bit-Machine, if the
   compiler has 64-bit numbers (in software or hardware).
   You will also need to choose it, if there would not be enough space
   for the type-bits in a 32-bit pointer.
 Model HEAPCODES stands for sizeof(oint) = sizeof(aint), and only minimal
   type information is stored in a pointer. All heap allocated objects
   (except conses) must contain the complete type and a length field in the
   first word. The heap gets somewhat bigger by this, and type tests require
   more memory accesses on average, but this model is portable even to
   systems whose memory map looks like Swiss Cheese. */

%% #if notused
%% #ifdef WIDE_HARD
%%   puts("#define WIDE_HARD");
%% #endif
%% #ifdef WIDE_SOFT
%%   puts("#define WIDE_SOFT");
%% #endif
%% #ifdef WIDE_AUXI
%%   puts("#define WIDE_AUXI");
%% #endif
%% #ifdef WIDE
%%   puts("#define WIDE");
%% #endif
%% #endif

#if defined(STANDARD_HEAPCODES) || defined(LINUX_NOEXEC_HEAPCODES)
  #define HEAPCODES
#endif

#if defined(WIDE_SOFT) && defined(HEAPCODES)
  #error "WIDE_SOFT and HEAPCODES make no sense together, no need for WIDE_SOFT"
#endif

#if !(defined(TYPECODES) || defined(HEAPCODES))
  # Choose typecodes on 64-bit machines (because there's enough room for type
  # bits), but not on 32-bit machines (because a 16 MB limit is ridiculous
  # today), except if the CPU cannot address more than 16 MB anyway.
  # HEAPCODES will normally not work if alignof(subr_t) = alignof(long) < 4,
  # but with egcs-1.1 or newer we can force alignof(subr_t) = 4.
  #if defined(WIDE_HARD) || defined(WIDE_SOFT) || defined(MC68000) || ((alignment_long < 4) && !defined(GNU))
    #define TYPECODES
  #else
    #define HEAPCODES
  #endif
#endif
%% #ifdef HEAPCODES
%%   puts("#define HEAPCODES");
%% #endif

#ifdef WIDE_SOFT
  #if defined(GNU) && !defined(WIDE_SOFT_LARGEFIXNUM)
    # Use the GNU-C extensions, to regard the wide oints as structs.
    #define WIDE_STRUCT
  #endif
  # defines the arrangement of an oint's elements:
  #define WIDE_ENDIANNESS true  # more efficient this way
#endif

#if defined(GNU) && (SAFETY >= 3)
  #if (__GNUC__ >= 2)
    #if (__GNUC__ > 2) || (__GNUC_MINOR__ >= 7) # circumvent gcc-2.6.3 bug
      # Typechecking by the C-compiler
      #define OBJECT_STRUCT
      #if !(defined(MC680X0) || defined(ARM)) && !(defined(__GNUG__) && (__GNUC__ == 3) && (__GNUC_MINOR__ == 3)) # only if struct_alignment==1, and not with g++ 3.3
        #define CHART_STRUCT
      #endif
    #endif
  #endif
#endif


# ###################### OS-related routines  ##################### #

# general standard constants for control chars:
#define BS    8  #  #\Backspace     Backspace
#define TAB   9  #  #\Tab           Tabulator
#define LF   10  #  #\Linefeed      linefeed
#define CR   13  #  #\Return        carriage return
#define PG   12  #  #\Page          form feed, new page

# Desired reaction when an I/O operation cannot be completed immediately.
typedef enum {
  persev_full,      /* Continue the I/O operation until the whole buffer is
                       handled or EOF or an error occurred. May hang. */
  persev_partial,   /* Continue the I/O operation until some (non-empty) part
                       of the buffer is handled or EOF or an error occurred.
                       May hang. */
  persev_immediate, /* Act immediately. Perform I/O only if we know in advance
                       that it will not block. In case of doubt, perform it
                       anyway. May return with 0 bytes handled. Does usually
                       not hang. */
  persev_bonus      /* Act immediately. Perform I/O only if we know in advance
                       that it will not block. In case of doubt, don't perform
                       it. May return with 0 bytes handled. Does not hang. */
} perseverance_t;
%% printf("typedef enum { persev_full=%d, persev_partial=%d, persev_immediate=%d, persev_bonus=%d } perseverance_t;\n",persev_full,persev_partial,persev_immediate,persev_bonus);

#if defined(UNIX) || defined(WIN32)

#ifdef UNIX
  #include "unix.c"
#endif
#ifdef WIN32_NATIVE
  #include "win32.c"
#endif
%% puts("#include <stdlib.h>");
%% puts("#include <sys/types.h>");
%% #if defined(UNIX)
%%   emit_typedef("int","Handle");
%%   emit_typedef("int","SOCKET");
%%   #ifdef UNIX_CYGWIN32
%%     puts("#include <windows.h>");
%%     puts("#undef WIN32");
%%     puts("extern long time_t_from_filetime (const FILETIME * ptr);");
%%     puts("extern void time_t_to_filetime (time_t time_in, FILETIME * out);");
%%   #endif
%% #elif defined(WIN32_NATIVE)
%%   puts("#include <windows.h>");
%%   export_def(Handle);
%%   puts("#include <winsock2.h>"); /* defines SOCKET */
%% #else
%%   puts("#error \"what is Handle on your platform?!\"");
%% #endif
%% #if defined(UNIX)
%%   puts("extern ssize_t fd_read (int fd, void* buf, size_t nbyte, perseverance_t persev);");
%%   puts("extern ssize_t fd_write (int fd, const void* buf, size_t nbyte, perseverance_t persev);");
%% #elif defined(WIN32_NATIVE)
%%   puts("extern ssize_t fd_read (Handle fd, void* buf, size_t nbyte, perseverance_t persev);");
%%   puts("extern ssize_t fd_write (Handle fd, const void* buf, size_t nbyte, perseverance_t persev);");
%% #endif

# execute statement on interrupt:
# interruptp(statement);
#if defined(UNIX) || defined(WIN32_NATIVE)
  # A keyboard interrupt (signal SIGINT, generated by Ctrl-C)
  # is pending for one second. It can be treated with 'interruptp' in
  # a continuing manner in that time. After this time has passed, the
  # program will be interrupted and can't be continued..
  #define PENDING_INTERRUPTS
  extern uintB interrupt_pending;
  #define interruptp(statement)  if (interrupt_pending) { statement; }
#endif
# used by EVAL, IO, SPVW, STREAM

#endif # UNIX || WIN32

#if (defined(UNIX) || defined(WIN32_NATIVE)) && !defined(NO_SIGSEGV)
  /* Support for fault handling. */
  #include <sigsegv.h>
  #if defined(UNIX_CYGWIN32)
    /* <sigsegv.h> includes <windows.h> */
    #undef WIN32
  #endif
#endif

/* Ignoring of a value (instead of assigning it to a variable)
 unused ...
 <sigsegv.h> includes <windows.h> which uses unused! */
#ifndef unused                  /* win32.d defines unused */
 #ifdef GNU     /* to prevent a gcc-warning "statement with no effect" */
  #define unused  (void)
 #else
  #define unused
 #endif
#endif
%% export_def(unused);

# Consensys and Solaris: "#define DS 3", "#define SP ESP", "#define EAX 11".
# Grr...
#undef DS
#undef SP
#undef EAX
# 386BSD does "#define CBLOCK 64". Grr...
#undef CBLOCK
# AIX 3.2.5 does "#define hz 100". Grr...
#undef hz
# MacOS X does "#define TIME_ABSOLUTE 0x00" and "#define TIME_RELATIVE 0x01".
# Grr...
#undef TIME_ABSOLUTE
#undef TIME_RELATIVE

#ifdef UNIX
  # Handling of UNIX errors
  # OS_error();
  # > int errno: error code
    nonreturning_function(extern, OS_error, (void));
  # used by SPVW, STREAM, PATHNAME, GRAPH
#endif
#if defined(WIN32_NATIVE)
  # Handling of Win32 errors
  # OS_error();
  # > GetLastError(): error code
    nonreturning_function(extern, OS_error, (void));
  # Handling of Winsock errors
  # SOCK_error();
  # > WSAGetLastError(): error code
    nonreturning_function(extern, SOCK_error, (void));
#endif
#if defined(DEBUG_OS_ERROR)
  # Show the file and line number of the caller of OS_error(). For debugging.
  #define OS_error()  \
    (fprintf(stderr,"\n[%s:%d] ",__FILE__,__LINE__), (OS_error)())
#endif
%% puts("nonreturning_function(extern, OS_error, (void));");

# Handling of ANSI C errors
# ANSIC_error();
# > int errno: error code
#ifdef UNIX
  #define ANSIC_error OS_error
#else
  nonreturning_function(extern, ANSIC_error, (void));
#endif
# used by SPVW, STREAM

#ifdef MULTITHREAD

#include "xthread.c"

#if !(defined(HAVE_MMAP_ANON) || defined(HAVE_MMAP_DEVZERO) || defined(HAVE_MACH_VM) || defined(HAVE_WIN32_VM))
  #error "Multithreading requires memory mapping facilities!"
#endif

#endif

# ##################### Further system-dependencies ##################### #

# At first dependencies that are visible to the LISP-level:

# setting of the table of character-names:
#ifdef WIN32
  #define WIN32_CHARNAMES
#endif
#ifdef UNIX
  #define UNIX_CHARNAMES
#endif
# When changed: extend CONSTOBJ, CHARSTRG, FORMAT.LISP.

# Whether to use the GNU gettext library for internationalization:
#if defined(ENABLE_NLS) && !defined(NO_GETTEXT)
  #define GNU_GETTEXT
#endif

# Whether to create a stream *KEYBOARD-INPUT*
# and whether it will be used for the stream *TERMINAL-IO*:
#if ((defined(UNIX) && !defined(NEXTAPP) || defined(MAYBE_NEXTAPP)) && !defined(NO_TERMCAP_NCURSES)) || defined(WIN32_NATIVE)
  #define KEYBOARD
  #if 0
    #define TERMINAL_USES_KEYBOARD
  #endif
#endif
# When changed: extend stream.d, keyboard.lisp

# Whether to use the GNU readline library for *TERMINAL-IO*:
#if defined(HAVE_READLINE) && !defined(NO_READLINE)
  #define GNU_READLINE
#endif

# Whether there are Window-streams and a package SCREEN:
#if defined(WIN32_NATIVE) || ((defined(UNIX) && !defined(NEXTAPP) || defined(MAYBE_NEXTAPP)) && !defined(NO_TERMCAP_NCURSES))
  #define SCREEN
#endif
# When changed: extend stream.d (loads of work!).

# Whether there are Pipe-streams:
#if defined(UNIX) || defined(WIN32_NATIVE)
  #define PIPES
  #if defined(UNIX) || defined(WIN32_NATIVE)
    #define PIPES2  # bidirectional pipes
  #endif
#endif
# When changed: extend stream.d and runprog.lisp.

# If the system has sockets, we support socket streams:
# We assume that if we have gethostbyname(), we have a networking OS
# (Unix or Win32). Then we decide independently about UNIX domain connections
# and TCP/IP connections.
#if defined(HAVE_GETHOSTBYNAME) # ==> defined(UNIX) || defined(WIN32_NATIVE)
  #ifdef HAVE_SYS_UN_H  # have <sys/un.h> and Unix domain sockets?
    #define UNIXCONN  # use Unix domain sockets
  #endif
  #if defined(HAVE_NETINET_IN_H) || defined(WIN32_NATIVE)  # have <netinet/in.h> ?
    #define TCPCONN  # use TCP/IP sockets
  #endif
  # Now, which kinds of socket streams:
  #define X11SOCKETS  # works even without TCPCONN (very young Linux)
  #ifdef TCPCONN
    #define SOCKET_STREAMS
  #endif
#endif
# When changed: extend stream.d, socket.d

# Whether there are generic streams:
#if 1
  #define GENERIC_STREAMS
#endif
# When changed: do nothing

# Whether the OS provides the required information for the
# functions  MACHINE-TYPE, MACHINE-VERSION, MACHINE-INSTANCE
#if defined(UNIX) || defined(WIN32_NATIVE)
  #define MACHINE_KNOWN
#endif
# When changed: extend misc.d, socket.d

# Whether there are LOGICAL-PATHNAMEs:
#if 1
  #define LOGICAL_PATHNAMES
#endif
# When changed: do nothing

# Whether the function USER-HOMEDIR-PATHNAME exists:
#if defined(UNIX) || defined(WIN32)
  #define USER_HOMEDIR
#endif
# When changed: extend pathname.d

# Whether the operating system manages an environment that associates Strings
# with Strings
#if defined(UNIX) || defined(WIN32)
  #define HAVE_ENVIRONMENT
#endif
# When changed: do nothing

# Whether the operating system has a preferred command-interpreter:
#if defined(UNIX) || defined(WIN32_NATIVE)
  #define HAVE_SHELL
#endif
# When changed: extend pathname.d

# Whether a foreign function interface is provided:
#if (defined(UNIX) && !defined(UNIX_BINARY_DISTRIB)) || defined(DYNAMIC_FFI)
  #define HAVE_FFI
#endif
#if 0
  #define HAVE_AFFI # Amiga-minded FFI
#endif
# When changed: ??

# Now the ones that are only relevant internally:

# Whether the GC closes files that aren't referenced any longer:
#if defined(UNIX) || defined(WIN32)
  #define GC_CLOSES_FILES
#endif
# When changed: do nothing

# How time is measured:
#ifdef UNIX
  #if defined(HAVE_GETTIMEOFDAY) || defined(HAVE_FTIME)
    #define TIME_UNIX
  #elif defined(HAVE_TIMES_CLOCK)
    #define TIME_UNIX_TIMES
  #endif
#endif
#ifdef WIN32_NATIVE
  #define TIME_WIN32
#endif
#if defined(TIME_UNIX_TIMES)
  # There's only a medium time resolution, so you can use 32-bit numbers
  # to store the time-differences without any problems.
  #define TIME_1
  # We fetch the time once on system sart. All further times are taken
  # relatively to that one.
  #define TIME_RELATIVE
#endif
#if defined(TIME_UNIX) || defined(TIME_WIN32)
  # The time resolution is so high that you need two 32-bit numbers to
  # measure time differences: seconds and and fractions of seconds.
  #define TIME_2
  # In this case we can use absolute and relative times for measurements.
  #define TIME_ABSOLUTE
#endif
# When changed: extend time.d

# Whether the operating system can give us the run-time, or whether we'll have
# to accumulate it ourselves:
#if defined(UNIX) || defined(WIN32_NATIVE)
  #define HAVE_RUN_TIME
#endif
# When changed: extend time.d

# Whether the operating system provides virtual memory.
#if (defined(UNIX) || defined(WIN32)) && !defined(NO_VIRTUAL_MEMORY)
  #define VIRTUAL_MEMORY
#endif
# When changed: do nothing

# Whether the operating system is capable of sending interruptions
# (Ctrl-C and others) as signal:
#if defined(UNIX)
  #define HAVE_SIGNALS
#endif
/* pass on for clx/new-clx */
%% #ifdef HAVE_SIGNALS
%%   puts("#define HAVE_SIGNALS");
%% #endif
# Whether we can even react to asynchronous signals:
# (If WIDE && !WIDE_HARD, writing a pointer is usually no elementary
# operation anymore!)
#if (defined(WIDE) && !defined(WIDE_HARD)) && !(defined(GNU) && defined(SPARC))
  #define NO_ASYNC_INTERRUPTS
#endif
#if defined(NO_ASYNC_INTERRUPTS) && defined(MULTITHREAD)
  #error "No multithreading possible with this memory model!"
#endif
# When changed: extend SPVW, write a interruptp().

# Flavors of Pathname-management:
#ifdef UNIX
  #define PATHNAME_UNIX
#endif
#ifdef WIN32
  #define PATHNAME_WIN32
#endif
# Components of pathnames:
#ifdef PATHNAME_WIN32
  #define HAS_HOST      1
  #define HAS_DEVICE    1
#endif
#ifdef PATHNAME_UNIX
  #define HAS_HOST      0
  #define HAS_DEVICE    0
#endif
# Handling of the file "extension" (pathname-type):
#if 0
  #define PATHNAME_EXT  # Name and Type are separated, so no limitation of the length
#endif
#if defined(PATHNAME_UNIX) || defined(PATHNAME_WIN32)
  #define PATHNAME_NOEXT  # no explicit extension.
#endif
# Whether "//" at the beginning of a pathname has to remain (and not to be shortened to "/"):
#ifdef UNIX_CYGWIN32
  #define PATHNAME_UNIX_UNC
#endif
# When changed: extend pathname.d

# Whether there is a type FOREIGN (a wrapper for various kinds of pointers):
#if defined(UNIX) || defined(DYNAMIC_FFI) || defined(WIN32_NATIVE)
  # (Used by FFI and by CLX.)
  #define FOREIGN  void*
#endif
# When changed: do nothing
%% #ifdef FOREIGN
%%   export_def(FOREIGN);
%% #endif

# Whether the STACK is checked at certain key points:
#define STACKCHECKS  (SAFETY >= 1) # when SUBRs and FSUBRs are called
#define STACKCHECKC  (SAFETY >= 1) # when compiled closures are interpreted
#define STACKCHECKR  (SAFETY >= 1) # in the reader
#define STACKCHECKP  (SAFETY >= 1) # in the printer
#define STACKCHECKB  (SAFETY >= 1) # in the bindings
# When changed: do nothing


# Feature dependent include files.

#ifdef HAVE_ICONV
  #include <stdlib.h>
  #include <iconv.h>
  #if _LIBICONV_VERSION
    # We use GNU libiconv.
    #define GNU_LIBICONV
    #define HAVE_GOOD_ICONV
  #elif (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 2))
    # glibc-2.2 iconv is also very reliable. Use it.
    #define HAVE_GOOD_ICONV
  #else
    # Other iconv implementations are too unreliable.
    # Don't define HAVE_GOOD_ICONV.
  #endif
#endif


# ############### List of implemented CLtL2-features ################ #

#define X3J13_005  # 18.5.1993
#define X3J13_014  # 22.1.1995
#define X3J13_149  # 22.7.1993
#define X3J13_175  # 25.7.1993


# ##################### Memory representation of objects #################### #

/*

Memory Representation and the Type Code of the various data types
=================================================================

1. The type code
----------------

An object consists of - in the same word - some type information and, for
immediate types, a couple of data bits, or, for heap allocated types,
a pointer to memory. There are many models of mixing type and pointer.
In the standard model, 6 to 8 bits (the word's high bits) are used for the
type. In the WIDE_HARD and WIDE_SOFT models, type and pointer are each 32
bits. In the HEAPCODES model, there are only 2 to 6 bits.

One bit (normally bit 31) is used as mark bit by the garbage collector.
Outside of GC, it is always cleared. (Except for the get_circularities and
subst_circ routines, and in the STACK, the GC bit is used for marking frames.)

2. Memory formats
-----------------

2.1. Immediate objects

2.1.1. Machine pointers

Machine pointers are immediate objects. They may point to the code area
(.text segment), to data areas (.bss, .data segments, malloc'ed areas).
Other values (e.g. pointers to text/data in shared libraries) are not
allowed, because they may contain bits which are interpreted as a type code.
To use such machine addresses, you must wrap them in foreign-pointers or
simple-bit-vectors.

2.1.2. Other immediate objects

Character, Fixnum, Short-Float, and, if IMMEDIATE_FFLOAT, Single-Float.
Furthermore: Frame-Pointer, Small-Read-Label, System. (System means some
finite number of special values, such as #<UNBOUND>.)

2.2. SUBRs

They are immediate in the sense that they do not move (they do not need to,
because they are allocated statically), but they have to be traversed by GC.

2.3. Pairs

These are heap objects containing just two pointers: Cons and, if SPVW_PURE,
Ratio and Complex.

2.4. Varobjects

These are heap objects of varying size. GC needs a header word at the
beginning of the object.

2.4.1. Records

These are varobjects which have additional type information and flags
in the second header word. Closure, Structure, Stream, Instance are always
records. Depending on the memory model, arrays, symbols etc. may or may
not be records.

2.4.2. Arrays

Simple-Bit-Vector, Simple-String, Simple-Vector are the "simple" arrays.
The non-simple ones are represented by a Iarray, yet the type code gives
some information about the rank, the representation and the element type:

                                |    "simple"     |  "not simple"  |
                                |    Sarray       |     Iarray     |
  ------------------------------+-----------------+----------------+
   (vector bit)                 | sbvector_type   | bvector_type   |
  ------------------------------+-----------------+----------------+
   (vector (unsigned-byte 2))   | sb2vector_type  | b2vector_type  |
  ------------------------------+-----------------+----------------+
   (vector (unsigned-byte 4))   | sb4vector_type  | b4vector_type  |
  ------------------------------+-----------------+----------------+
   (vector (unsigned-byte 8))   | sb8vector_type  | b8vector_type  |
  ------------------------------+-----------------+----------------+
   (vector (unsigned-byte 16))  | sb16vector_type | b16vector_type |
  ------------------------------+-----------------+----------------+
   (vector (unsigned-byte 32))  | sb32vector_type | b32vector_type |
  ------------------------------+-----------------+----------------+
   (vector character)           | sstring_type    | string_type    |
  ------------------------------+-----------------+----------------+
   (vector t)                   | svector_type    | vector_type    |
  ------------------------------+-----------------+----------------+
   array of dimension /= 1      |       --        |  mdarray_type  |
  ------------------------------+-----------------+----------------+

2.4.3. Other varobjects

Symbol has some special flags (keyword, constant, special) in the header,
if possible.

FSUBR, Bignum, Single-Float (unless IMMEDIATE_FFLOAT), Double-Float,
Long-Float, Ratio and Complex (only if SPVW_MIXED).

*/

# ######################## LISP-objects in general ######################### #

#if defined(DEBUG_GCSAFETY)
  #ifndef __cplusplus
    #error "DEBUG_GCSAFETY works only with a C++ compiler! Reconfigure with CC=g++."
  #endif
  #if defined(WIDE_SOFT) || defined(WIDE_AUXI)
    #error "DEBUG_GCSAFETY cannot be used together with WIDE_SOFT or WIDE_AUXI (not yet implemented)!"
  #endif
  # The 'gcv_object_t' and 'object' types share the major part of their innards.
  #ifndef OBJECT_STRUCT
    #define OBJECT_STRUCT
  #endif
#endif

/* The type 'object' denotes an object in registers or in memory that is
 not seen by the GC.

 The type `gcv_object_t' denotes a GC visible object, i.e. a slot inside
 a heap-allocated object or a STACK slot. If its value is not an immediate
 object, any call that can trigger GC can modify the pointer value.
 NEVER write "var gcv_object_t foo;" - this is forbidden!
 You can write "var gcunsafe_object_t foo;" instead - but then you must not
 trigger GC during the entire lifetime of the variable 'foo'!
*/

%% #if (defined(WIDE_AUXI) || defined(OBJECT_STRUCT) || defined(WIDE_STRUCT)) && defined(WIDE) && !defined(WIDE_HARD) && defined(GENERATIONAL_GC)
%%   #define attribute_aligned_object " __attribute__ ((aligned(8)))"
%% #else
%%   #define attribute_aligned_object ""
%% #endif

#if !defined(WIDE_SOFT)

  # An object pointer is an empty pointer to begin with (so you cannot do
  # anything unwanted with it in C):
  #if defined(WIDE_AUXI)
    # Make room for an auxiliary word in every object.
    # The struct around the union is needed to work around a gcc-2.95 bug.
    #if BIG_ENDIAN_P
      #define INNARDS_OF_GCV_OBJECT  \
        union {                                         \
          struct { uintP auxi_ob; uintP one_ob; } both; \
          oint align_o _attribute_aligned_object_;      \
        } u _attribute_aligned_object_;
    #else
      #define INNARDS_OF_GCV_OBJECT  \
        union {                                         \
          struct { uintP one_ob; uintP auxi_ob; } both; \
          oint align_o _attribute_aligned_object_;      \
        } u _attribute_aligned_object_;
    #endif
    #define one_o  u.both.one_ob
    #define auxi_o  u.both.auxi_ob
  #elif defined(OBJECT_STRUCT)
    #define INNARDS_OF_GCV_OBJECT  \
      uintP one_o;
  #else
    typedef  void *  gcv_object_t;
  #endif
  # But there is an address and type bits in the representation.

  # An (unsigned) Integer of the object's size:
  #ifdef WIDE_AUXI
    typedef  uint64  oint;
    typedef  sint64  soint;
  #else
    typedef  uintP  oint;
    typedef  sintP  soint;
  #endif

#else # defined(WIDE_SOFT)

  # An object consists of a separated 32 bit address and a 32 bit type info.
  typedef  uint64  oint;
  typedef  sint64  soint;
  #ifdef WIDE_STRUCT
    # The struct around the union is needed to work around a gcc-2.95 bug.
    #if BIG_ENDIAN_P==WIDE_ENDIANNESS
      #define INNARDS_OF_GCV_OBJECT                                      \
        union {                                                          \
          struct { /* tint */ uintL type; /* aint */ uintL addr; } both; \
          oint one_u _attribute_aligned_object_;                         \
        } u _attribute_aligned_object_;
    #else
      #define INNARDS_OF_GCV_OBJECT                                      \
        union {                                                          \
          struct { /* aint */ uintL addr; /* tint */ uintL type; } both; \
          oint one_u _attribute_aligned_object_;                         \
        } u _attribute_aligned_object_;
    #endif
    #define one_o  u.one_u
  #else
    typedef  oint  gcv_object_t;
  #endif

#endif

# sizeof(gcv_object_t) = sizeof(oint) must hold true!

%% #if !defined(WIDE_SOFT)
%%   #if defined(WIDE_AUXI)
%%     strcpy(buf,"struct { union { struct { ");
%%     #if BIG_ENDIAN_P
%%       strcat(buf,"uintP auxi_ob; uintP one_ob;");
%%     #else
%%       strcat(buf,"uintP one_ob; uintP auxi_ob;");
%%     #endif
%%     strcat(buf," } both; oint align_o");
%%     strcat(buf,attribute_aligned_object);
%%     strcat(buf,"; } u");
%%     strcat(buf,attribute_aligned_object);
%%     strcat(buf,"; }");
%%     emit_typedef(buf,"gcv_object_t");
%%     emit_define("one_o","u.both.one_ob");
%%     emit_define("auxi_o","u.both.auxi_ob");
%%   #elif defined(OBJECT_STRUCT)
%%     #ifdef DEBUG_GCSAFETY
%%       puts("struct object { uintP one_o; uintL allocstamp; };");
%%       puts("struct gcv_object_t { uintP one_o; operator object () const; gcv_object_t (object obj); gcv_object_t (struct fake_gcv_object obj); gcv_object_t (); };");
%%     #else
%%       emit_typedef("struct { uintP one_o; }","gcv_object_t");
%%     #endif
%%   #else
%%     emit_typedef("void *","gcv_object_t");
%%   #endif
%%   #ifdef WIDE_AUXI
%%     emit_typedef("uint64","oint");
%%     emit_typedef("sint64","soint");
%%   #else
%%     emit_typedef("uintP","oint");
%%     emit_typedef("sintP","soint");
%%   #endif
%% #else
%%   emit_typedef("uint64","oint");
%%   emit_typedef("sint64","soint");
%%   #ifdef WIDE_STRUCT
%%     strcpy(buf,"struct { union {\n");
%%     #if BIG_ENDIAN_P==WIDE_ENDIANNESS
%%       strcat(buf,"  struct { /*tint*/ uintL type; /*aint*/ uintL addr; } both;\n");
%%     #else
%%       strcat(buf,"  struct { /*aint*/ uintL addr; /*tint*/ uintL type; } both;\n");
%%     #endif
%%     strcat(buf,"  oint one_u");
%%     strcat(buf,attribute_aligned_object);
%%     strcat(buf,"; } u");
%%     strcat(buf,attribute_aligned_object);
%%     strcat(buf,"; }");
%%     emit_typedef(buf,"gcv_object_t");
%%     emit_define("one_o","u.one_u");
%%   #else
%%     emit_typedef("oint","gcv_object_t");
%%   #endif
%% #endif

# conversion between gcv_object_t/object and oint:
# as_oint(expr)   gcv_object_t/object --> oint
# as_object(x)    oint --> gcv_object_t
# The conversion  gcv_object_t --> object
# is implicit.
#if defined(WIDE_STRUCT) || defined(OBJECT_STRUCT)
  #define as_oint(expr)  ((expr).one_o)
  #if defined(WIDE_STRUCT)
    #define as_object(o)  ((object){u:{one_u:(o)}INIT_ALLOCSTAMP})
  #elif defined(OBJECT_STRUCT)
    #define as_object(o)  ((object){one_o:(o)INIT_ALLOCSTAMP})
  #else
    extern __inline__ object as_object (register oint o)
      { register object obj; obj.one_o = o; return obj; }
  #endif
#elif defined(WIDE_AUXI)
  #define as_oint(expr)  ((expr).u.align_o)
  # These could store arbitrary information in auxi_o.
  #define as_object_with_auxi(o)  ((object){u:{both:{ one_ob: (o), auxi_ob: 0 }} INIT_ALLOCSTAMP })
  #define as_object(o)  ((object){u:{align_o:(o)}INIT_ALLOCSTAMP})
#else
  #define as_oint(expr)  (oint)(expr)
  #define as_object(o)  (gcv_object_t)(o)
#endif
%% export_def(as_oint(expr));
%% export_def(as_object(o));

# Separation of an oint in type bits and address:
# oint_type_mask  is always subset  (2^oint_type_len-1)<<oint_type_shift
# and  oint_addr_mask superset (2^oint_addr_len-1)<<oint_addr_shift .
#if !defined(TYPECODES)
  # HEAPCODES model:
  # For pointers, the address takes the full word (with type info in the
  # lowest two bits). For immediate objects, we use 24 bits for the data
  # (but exclude the highest available bit, which is the garcol_bit).
  #if !(defined(STANDARD_HEAPCODES) || defined(LINUX_NOEXEC_HEAPCODES))
    # Choose the appropriate HEAPCODES variant for the machine.
    # On most systems, one of the high bits is suitable as GC bit; here we
    # choose STANDARD_HEAPCODES.
    # On some Linux/x86 systems, starting in 2004, a "no-exec" kernel patch
    # is used that distributes virtual addresses over the entire address
    # space from 0x00000000 to 0xBFFFFFFF (as a function of its access
    # permissions); here we use LINUX_NOEXEC_HEAPCODES.
    #if defined(I80386) && defined(UNIX_LINUX)
      #define LINUX_NOEXEC_HEAPCODES
    #else
      #define STANDARD_HEAPCODES
    #endif
  #endif
  #ifdef STANDARD_HEAPCODES
    # The portable case. Assumes only that the GC bit can be chosen.
    #if defined(SPARC) && defined(UNIX_LINUX) && (__GLIBC__ < 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ < 2))
      #define LINUX_SPARC_OLD_GLIBC
    #endif
    #if defined(WIDE_HARD)
      #define oint_type_shift 0
      #define oint_type_len 8
      #define oint_type_mask 0x000000000000007FUL
      #define oint_data_shift 7
      #define oint_data_len 48
      #define oint_data_mask 0x007FFFFFFFFFFF80UL
      #define garcol_bit_o 63
    #elif !((defined(MC680X0) && defined(UNIX_LINUX)) || (defined(I80386) && defined(UNIX_BEOS)) || defined(LINUX_SPARC_OLD_GLIBC))
      #define oint_type_shift 0
      #define oint_type_len 8
      #define oint_type_mask 0x0000007FUL
      #define oint_data_shift 7
      #define oint_data_len 24
      #define oint_data_mask 0x7FFFFF80UL
      #define garcol_bit_o 31
    #elif defined(I80386) && defined(UNIX_BEOS)
      # On BeOS 5, malloc()ed addresses are of the form 0x80...... Bit 31
      # is therefore part of an address and cannot be used as garcol_bit.
      #define oint_type_shift 0
      #define oint_type_len 8
      #define oint_type_mask 0x0000003FUL
      #define oint_data_shift 6
      #define oint_data_len 24
      #define oint_data_mask 0x3FFFFFC0UL
      #define garcol_bit_o 30
    #elif (defined(MC680X0) && defined(UNIX_LINUX)) || defined(LINUX_SPARC_OLD_GLIBC)
      # On Sparc-Linux with glibc 2.1 and older:
      # malloc()ed addresses are of the form 0x0....... or 0xe........
      # Bits 31..29 are therefore part of an address and cannot
      # be used as garcol_bit. We therefore choose bit 28 as garcol_bit.
      # Now, the 24 data bits of an immediate value must not intersect the
      # garcol_bit, so we use bits 27..4 for that (we could use bits 26..3
      # as well).
      # On m68k-Linux, malloc()ed addresses are of the form 0x80...... or
      # 0xc0....... Bits 31..30 are therefore part of an address and cannot
      # be used as garcol_bit. We therefore have three choices:
      #   data bits: bits 26..3, garcol_bit_o = 28/27
      #   data bits: bits 27..4, garcol_bit_o = 28/3
      #   data bits: bits 28..5, garcol_bit_o = 4/3
      #define oint_type_shift 0
      #define oint_type_len 32
      #define oint_type_mask 0xE000000FUL
      #define oint_data_shift 4
      #define oint_data_len 24
      #define oint_data_mask 0x0FFFFFF0UL
      #define garcol_bit_o 28
    #endif
  #endif # STANDARD_HEAPCODES
  #ifdef LINUX_NOEXEC_HEAPCODES
    # The Linux/32-bit case. Assumes 1. that the virtual memory addresses end
    # at 0xC0000000, or at least that we can put a black hole on the range
    # 0xC0000000..0xDFFFFFFF, 2. that the compiler and linker can enforce an
    # 8-byte alignment of symbol_tab and subr_tab.
    # Only bit 0 or 1 can be used as GC-bit.
    #define oint_type_shift 0
    #define oint_type_len 32
    #define oint_type_mask 0xE000001FUL
    #define oint_data_shift 5
    #define oint_data_len 24
    #define oint_data_mask 0x1FFFFFE0UL
    #define garcol_bit_o 0
  #endif # LINUX_NOEXEC_HEAPCODES
  #if defined(WIDE_HARD)
    #define oint_addr_shift 0
    #define oint_addr_len 64
    #define oint_addr_mask 0xFFFFFFFFFFFFFFFFUL
  #else
    #define oint_addr_shift 0
    #define oint_addr_len 32
    #define oint_addr_mask 0xFFFFFFFFUL
  #endif
/* Now come the platforms with TYPECODES. oint_type_len should be >= 8,
 and oint_type_mask should have at least 8 bits set and at most one bit in
 common with oint_addr_mask. */
#elif defined(WIDE_HARD)
  #if defined(DECALPHA) && (defined(UNIX_OSF) || defined(UNIX_LINUX) || defined(UNIX_FREEBSD) || defined(UNIX_NETBSD))
  /* UNIX_OSF:
     Ordinary pointers are in the range 1*2^32..2*2^32.
     Code address range:    0x000000012xxxxxxx
     Malloc address range:  0x000000014xxxxxxx
     Shared libraries:      0x000003FFCxxxxxxx
   UNIX_LINUX:
     Code address range:    0x000000012xxxxxxx
     Malloc address range:  0x000000012xxxxxxx
                      and:  0x0000015555xxxxxx
     Shared libraries:      0x0000015555xxxxxx
     Virtual address limit: 0x0000040000000000
   UNIX_FREEBSD
     Code address range:    0x0000000120000000
     Malloc address range:  0x0000000120000000
     Shared libraries:      0x0000000160000000
     Stack address range:   0x0000000011000000
   UNIX_NETBSD
     Code address range:    0x0000000120000000
     Malloc address range:  0x0000000120000000
     Shared libraries:      0x0000000160000000
     Stack address range:   0x00000001FF000000 */
    # This is the safest.
    # Bits 63..48 = type code, Bits 47..0 = address
    #define oint_type_shift 48
    #define oint_type_len 16
    #define oint_type_mask 0xFFFF000000000000UL
    #define oint_addr_shift 0
    #define oint_addr_len 48
    #define oint_addr_mask 0x0000FFFFFFFFFFFFUL
    #define oint_data_shift oint_addr_shift
    #define oint_data_len oint_addr_len
    #define oint_data_mask oint_addr_mask
  #endif
  #if defined(MIPS64)
    # Bits 63..48 = type code, bits 31..0 = address
    #define oint_type_shift 48
    #define oint_type_len 16
    #define oint_type_mask 0xFFFF000000000000UL
    #define oint_addr_shift 0
    #define oint_addr_len 64
    #define oint_addr_mask 0x00000000FFFFFFFFUL
    #define oint_data_shift 0
    #define oint_data_len 48
    #define oint_data_mask 0x0000FFFFFFFFFFFFUL
  #endif
  #if defined(SPARC64)
    # Virtual address limit on some systems: -2^43..2^43.
    # This is the safest.
    # Bits 63..48 = type code, bits 47..0 = address
    #define oint_type_shift 48
    #define oint_type_len 16
    #define oint_type_mask 0xFFFF000000000000UL
    #define oint_addr_shift 0
    #define oint_addr_len 48
    #define oint_addr_mask 0x0000FFFFFFFFFFFFUL
    #define oint_data_shift oint_addr_shift
    #define oint_data_len oint_addr_len
    #define oint_data_mask oint_addr_mask
  #endif
  #if defined(IA64) && defined(UNIX_LINUX)
    # Bits 63..61 = region code,
    # bits 60..39 all zero or all one,
    # virtual address limit: R*2^61..R*2^61+2^39, (R+1)*2^61-2^39..(R+1)*2^61.
    # SHLIB_ADDRESS_RANGE  = 0x2000000000000000UL (region 1)
    # CODE_ADDRESS_RANGE   = 0x4000000000000000UL (region 2)
    # MALLOC_ADDRESS_RANGE = 0x6000000000000000UL (region 3)
    # STACK_ADDRESS_RANGE  = 0x9FFFFFFFFF000000UL (region 4)
    # This is the safest.
    # Bits 63..48 = Typcode, Bits 47..0 = address
    #define oint_type_shift 48
    #define oint_type_len 16
    #define oint_type_mask 0x1FFF000000000000UL
    #define oint_addr_shift 0
    #define oint_addr_len 64
    #define oint_addr_mask 0xE000FFFFFFFFFFFFUL
    #define oint_data_shift 0
    #define oint_data_len 48
    #define oint_data_mask 0x0000FFFFFFFFFFFFUL
  #endif
  #if defined(AMD64)
   /* UNIX_LINUX:
     CODE_ADDRESS_RANGE     0x0000000000000000UL
     MALLOC_ADDRESS_RANGE   0x0000000000000000UL
     SHLIB_ADDRESS_RANGE    0x00000034F5000000UL
     STACK_ADDRESS_RANGE    0x0000007FBF000000UL
   UNIX_FREEBSD
     CODE_ADDRESS_RANGE     0x0000000000000000UL
     MALLOC_ADDRESS_RANGE   0x0000000000000000UL
     SHLIB_ADDRESS_RANGE    0x0000000800000000UL
     STACK_ADDRESS_RANGE    0x00007FFFFF000000UL
     Bits 63..48 = type code, Bits 47..0 = address */
    #define oint_type_shift 48
    #define oint_type_len 16
    #define oint_type_mask 0xFFFF000000000000UL
    #define oint_addr_shift 0
    #define oint_addr_len 48
    #define oint_addr_mask 0x0000FFFFFFFFFFFFUL
    #define oint_data_shift oint_addr_shift
    #define oint_data_len oint_addr_len
    #define oint_data_mask oint_addr_mask
  #endif
#elif defined(WIDE_SOFT)
  # separate one 32-bit word for typcode and address.
  #if defined(WIDE_SOFT_LARGEFIXNUM)
    # Used to test large fixnums on 32-bit platforms.
    # Bits 63..48 = Typcode, Bits 47..32 = zero, Bits 31..0 = address
    #define oint_type_shift 48
    #define oint_type_len 16
    #define oint_type_mask ULL(0xFFFF000000000000)
    #define oint_addr_shift 0
    #define oint_addr_len 48
    #define oint_addr_mask ULL(0x0000FFFFFFFFFFFF)
  #elif WIDE_ENDIANNESS
    # Bits 63..32 = Typcode, Bits 31..0 = address
    #define oint_type_shift 32
    #define oint_type_len 32
    #define oint_type_mask ULL(0xFFFFFFFF00000000)
    #define oint_addr_shift 0
    #define oint_addr_len 32
    #define oint_addr_mask ULL(0x00000000FFFFFFFF)
  #else # conversely it is a little slower:
    # Bits 63..32 = Adress, Bits 31..0 = Typcode
    #define oint_type_shift 0
    #define oint_type_len 32
    #define oint_type_mask ULL(0x00000000FFFFFFFF)
    #define oint_addr_shift 32
    #define oint_addr_len 32
    #define oint_addr_mask ULL(0xFFFFFFFF00000000)
  #endif
# Now come the 32-bit platforms with TYPECODES. We need to support it only on
# MC680X0 platforms without new gcc.
# It worked on the following platforms in the past, and may still work on:
#   (defined(MC680X0) && !defined(UNIX_AMIX) && !defined(UNIX_NEXTSTEP) && !(defined(UNIX_LINUX) && CODE_ADDRESS_RANGE))
#   (defined(I80386) && !(defined(UNIX_LINUX) && (CODE_ADDRESS_RANGE != 0)) && !defined(UNIX_HURD) && !defined(UNIX_SYSV_UHC_1) && !defined(UNIX_NEXTSTEP) && !defined(UNIX_SYSV_PTX) && !defined(UNIX_SUNOS5) && !defined(UNIX_CYGWIN32) && !defined(WIN32_NATIVE))
#   (defined(SPARC) && !defined(SUN4_29))
#   (defined(MIPS) && !defined(UNIX_IRIX) && !defined(UNIX_DEC_ULTRIX))
#   defined(M88000)
#   (defined(POWERPC) && !defined(UNIX_AIX) && !defined(UNIX_LINUX))
#   defined(VAX)
#elif (defined(I80386) && ((defined(UNIX_LINUX) && (CODE_ADDRESS_RANGE != 0)) || (defined(UNIX_FREEBSD) && !defined(UNIX_GNU)))) || (defined(POWERPC) && defined(UNIX_DARWIN)) || defined(TRY_TYPECODES_1)
  # You can add more platforms here provided that
  # 1. you need it,
  # 2. CODE_ADDRESS_RANGE | MALLOC_ADDRESS_RANGE has at most one bit set,
  # 3. it works.
  #define oint_type_shift 24
  #define oint_type_len 8
  #define oint_type_mask (0xFF000000UL & ~(CODE_ADDRESS_RANGE | MALLOC_ADDRESS_RANGE))
  #define oint_addr_shift 0
  #define oint_addr_len 24
  #define oint_addr_mask (0x00FFFFFFUL | CODE_ADDRESS_RANGE | MALLOC_ADDRESS_RANGE)
  #define oint_data_shift 0
  #define oint_data_len 24
  #define oint_data_mask 0x00FFFFFFUL
#elif 0 || defined(TRY_TYPECODES_2)
  # You can add more platforms here provided that
  # 1. you need it,
  # 2. it works.
  # Bits 31..24 = Typcode, Bits 23..0 = Adress
  #define oint_type_shift 24
  #define oint_type_len 8
  #define oint_type_mask 0xFF000000UL
  #define oint_addr_shift 0
  #define oint_addr_len 24
  #define oint_addr_mask 0x00FFFFFFUL
#else
  #error "TYPECODES maybe not supported any more on this platform. Try defining TRY_TYPECODES_1 or TRY_TYPECODES_2, or use -DHEAPCODES."
#endif
%% #if notused
%%  export_def(oint_type_shift);
%%  export_def(oint_type_len);
%%  export_def(oint_type_mask);
%%  export_def(oint_addr_shift);
%%  export_def(oint_addr_len);
%%  export_def(oint_addr_mask);
%% #endif
#ifndef oint_type_len
#error "CLISP has not been ported to this platform - oint_type_len undefined"
#endif

# Generally we use all of the space of an address for the data of Fixnums etc.
# Always     [oint_data_shift..oint_data_shift+oint_data_len-1] is subset of
#            [oint_addr_shift..oint_addr_shift+oint_addr_len-1],
# thus       oint_data_len <= oint_addr_len.
#ifndef oint_data_len
  #define oint_data_shift oint_addr_shift
  #define oint_data_len oint_addr_len
  #define oint_data_mask oint_addr_mask
#endif
%% #if notused
%%  export_def(oint_data_shift);
%%  export_def(oint_data_len);
%%  export_def(oint_data_mask);
%% #endif

# Integer type for typebits:
typedef unsigned_int_with_n_bits(oint_type_len)  tint;
%% sprintf(buf,"uint%d",oint_type_len); emit_typedef(buf,"tint");

# Integer type for addresses:
typedef unsigned_int_with_n_bits(oint_addr_len)  aint;
typedef signed_int_with_n_bits(oint_addr_len)  saint;
%% sprintf(buf,"uint%d",oint_addr_len); emit_typedef(buf,"aint");
%% #if notused
%% sprintf(buf,"sint%d",oint_addr_len); emit_typedef(buf,"saint");
%% #endif

# Integer type for immediate values:
# Always 32 = intLsize <= intVsize <= 64.
#if (oint_data_len <= 32)
  #define intVsize 32
#else
  #define intVsize 64
#endif
typedef unsigned_int_with_n_bits(intVsize)  uintV;
typedef signed_int_with_n_bits(intVsize)  sintV;
%% sprintf(buf,"uint%d",intVsize); emit_typedef(buf,"uintV");
%% sprintf(buf,"sint%d",intVsize); emit_typedef(buf,"sintV");

# Integer type used to represent an amount of memory:
# (This may be larger than size_t or ptrdiff_t: size_t is required by ISO C to
# be enough for the size of a single memory block; ptrdiff_t is required by
# ISO C to be enough for the size of a single memory block plus room for a sign
# bit. But on segmented architectures which allow many medium-sized memory
# blocks, like the 80286 was, the total available memory size may be bigger.
# Also, we avoid size_t because it's likely to be wrong on 64-bit Woe32.)
#if !defined(WIDE_HARD) || ((oint_addr_mask & ~0xFFFFFFFFUL) == 0)
  # A 32-bit integer is sufficient.
  #define intMsize  intLsize
  typedef uintL uintM;
  typedef sintL sintM;
#else
  # An integer as wide as a pointer may be required.
  #define intMsize  pointer_bitsize
  typedef uintP uintM;
  typedef sintP sintM;
#endif

# Number of bits by which an address is finally being shifted:
#ifndef addr_shift
  #define addr_shift 0
#endif
%% #if notused
%%  export_def(addr_shift);
%% #endif

# Verify the values w.r.t. the autoconfigured CODE_ADDRESS_RANGE and
# MALLOC_ADDRESS_RANGE values.
#if !defined(WIDE_SOFT)
  #if (CODE_ADDRESS_RANGE >> addr_shift) & ~(oint_addr_mask >> oint_addr_shift)
     #error "oint_addr_mask doesn't cover CODE_ADDRESS_RANGE !!"
  #endif
  #if (MALLOC_ADDRESS_RANGE >> addr_shift) & ~(oint_addr_mask >> oint_addr_shift)
     #error "oint_addr_mask doesn't cover MALLOC_ADDRESS_RANGE !!"
  #endif
#endif


#if (oint_addr_shift == 0) && (addr_shift == 0) && defined(TYPECODES) && !defined(WIDE_SOFT) && !(defined(SUN3) && !defined(UNIX_SUNOS4) && !defined(WIDE_SOFT)) && !(defined(AMD64) && defined(UNIX_LINUX))
# If the address bits are the lower ones and not WIDE_SOFT,
# memory mapping may be possible.

  #if (defined(HAVE_MMAP_ANON) || defined(HAVE_MMAP_DEVZERO) || defined(HAVE_MACH_VM) || defined(HAVE_WIN32_VM)) && !defined(MULTIMAP_MEMORY) && !(defined(UNIX_SINIX) || defined(UNIX_AIX)) && !defined(NO_SINGLEMAP)
    # Access to LISP-objects is made easier by putting each LISP-object
    # to an address that already contains its type information.
    # But this does not work on SINIX and AIX.
      #define SINGLEMAP_MEMORY
  #endif

  #if defined(UNIX_SUNOS4) && !defined(MULTIMAP_MEMORY) && !defined(SINGLEMAP_MEMORY) && !defined(NO_MULTIMAP_FILE)
    # Access to Lisp-objects is done through memory-mapping: Each
    # memory page can be accessed at several addresses.
      #define MULTIMAP_MEMORY
      #define MULTIMAP_MEMORY_VIA_FILE
  #endif

  #if defined(HAVE_SHM) && !defined(MULTIMAP_MEMORY) && !defined(SINGLEMAP_MEMORY) && !defined(NO_MULTIMAP_SHM)
    # Access to Lisp-objects is done through memory-mapping: Each
    # memory page can be accessed at several addresses.
      #define MULTIMAP_MEMORY
      #define MULTIMAP_MEMORY_VIA_SHM
  #endif

  #if (defined(UNIX_LINUX) || defined(UNIX_FREEBSD)) && !defined(MULTIMAP_MEMORY) && !defined(SINGLEMAP_MEMORY) && !defined(NO_MULTIMAP_FILE)
     # Access to Lisp-objects is done through memory-mapping: Each
     # memory page can be accessed at several addresses.
      #define MULTIMAP_MEMORY
      #define MULTIMAP_MEMORY_VIA_FILE
  #endif

#endif

#if defined(MULTIMAP_MEMORY) || defined(SINGLEMAP_MEMORY)
  #define MAP_MEMORY
#endif

#if (defined(HAVE_MMAP_ANON) || defined(HAVE_MMAP_DEVZERO) || defined(HAVE_MACH_VM) || defined(HAVE_WIN32_VM)) && !defined(MAP_MEMORY) && !(defined(UNIX_HPUX) || defined(UNIX_AIX)) && !defined(NO_TRIVIALMAP)
  # mmap() allows for a more flexible way of memory management than malloc().
  # It's not really memory-mapping, but a more comfortable way to
  # manage two large memory blocks.
  # But is doesn't work on HP-UX 9 and AIX.
  #define TRIVIALMAP_MEMORY
#endif


# Flavor of the garbage collection: normal or generational.
#if defined(VIRTUAL_MEMORY) && (defined(SINGLEMAP_MEMORY) || defined(TRIVIALMAP_MEMORY) || (defined(MULTIMAP_MEMORY) && (defined(UNIX_LINUX) || defined(UNIX_FREEBSD)))) && defined(HAVE_WORKING_MPROTECT) && defined(HAVE_SIGSEGV_RECOVERY) && !defined(UNIX_IRIX) && !defined(WIDE_SOFT_LARGEFIXNUM) && (SAFETY < 3) && !defined(NO_GENERATIONAL_GC)
  # "generational garbage collection" has some requirements.
  # With Linux, it will only work with 1.1.52, and higher, which will be checked in makemake.
  # On IRIX 6, it worked in the past, but leads to core dumps now. Reason unknown. FIXME!
  #define GENERATIONAL_GC
#endif


#ifdef MAP_MEMORY
  # Some type-bit combinations might not be allowed
  #ifdef vm_addr_mask
    #define tint_allowed_type_mask  ((oint_type_mask & vm_addr_mask) >> oint_type_shift)
  #endif
#endif


# Complete the definition of the type 'gcv_object_t'.
#if defined(WIDE_AUXI) || defined(OBJECT_STRUCT) || defined(WIDE_STRUCT)
  #if defined(WIDE) && !defined(WIDE_HARD)
    #ifdef GENERATIONAL_GC
      # The generational GC can't deal with an object-pointer that points
      # towards two memory pages.
      # Thus we enforce alignof(gcv_object_t) = sizeof(gcv_object_t).
      #define _attribute_aligned_object_  __attribute__ ((aligned(8)))
    #else
      #define _attribute_aligned_object_
    #endif
  #endif
  #ifdef DEBUG_GCSAFETY
    struct object;
    struct gcv_object_t {
      INNARDS_OF_GCV_OBJECT
      # Conversion to object.
      operator object () const;
      # Conversion from object.
      gcv_object_t (object obj);
      # Conversion from fake_gcv_object.
      gcv_object_t (struct fake_gcv_object obj);
      # Uninitialized object.
      gcv_object_t ();
    };
  #else
    typedef struct { INNARDS_OF_GCV_OBJECT } gcv_object_t;
  #endif
#endif
#ifndef _attribute_aligned_object_
  #define _attribute_aligned_object_
#endif


# Define the type 'object'.
#ifdef DEBUG_GCSAFETY
  # A counter that is incremented each time an allocation occurs that could
  # trigger GC.
  extern uintL alloccount;

  # A register-allocated object contains, if not GC-invariant, the timestamp
  # of when it was fetched from a GC-visible location.
  struct object {
    INNARDS_OF_GCV_OBJECT
    uintL allocstamp;
  };
  # Always initialize allocstamp with the current(!) value of alloccount.
  #define INIT_ALLOCSTAMP  , allocstamp: alloccount
#else
  typedef gcv_object_t object;
  #define INIT_ALLOCSTAMP
#endif
%% #ifdef DEBUG_GCSAFETY
%%   puts("extern uintL alloccount;");
%% #else
%%   emit_typedef("gcv_object_t","object");
%% #endif

# fake_gcv_object(value)
# creates a gcv_object that is actually not seen by GC,
# for use as second word in SKIP2 frames.
#ifdef DEBUG_GCSAFETY
  struct fake_gcv_object {
    oint fake_value;
    fake_gcv_object (oint value) : fake_value (value) {}
  };
#else
  #define fake_gcv_object(value)  as_object((oint)(value))
#endif

# Hack for use only in areas where no GC can be triggered.
#ifdef DEBUG_GCSAFETY
  struct gcunsafe_object_t : gcv_object_t {
    uintL allocstamp;
    # Conversion from object.
    gcunsafe_object_t (object obj);
    # Conversion from gcv_object_t.
    gcunsafe_object_t (gcv_object_t obj);
    # Verification that no GC has been triggered.
    ~gcunsafe_object_t ();
  };
#else
  typedef gcv_object_t gcunsafe_object_t;
#endif


# mask of those bits of a tint, which really belong to the type:
# tint_type_mask = oint_type_mask >> oint_type_shift
# (a constant expression, without any 'long long's in it!)
#ifdef WIDE_SOFT
  #define tint_type_mask  (bitm(oint_type_len)-1)
#else
  #define tint_type_mask  (oint_type_mask >> oint_type_shift)
#endif
%% #if notused
%% export_def(tint_type_mask);
%% #endif

# To add something to an object/oint:
# objectplus(obj,offset)
#if !(defined(WIDE_SOFT) || defined(WIDE_AUXI) || defined(OBJECT_STRUCT))
  #define objectplus(obj,offset)  ((object)pointerplus(obj,offset))
#elif defined(WIDE_AUXI)
  static inline object objectplus (object obj, saint offset) {
    return (object){u:{both:{ one_ob: obj.one_o+offset, auxi_ob: obj.auxi_o }}};
  }
#else # defined(WIDE_SOFT) || defined(OBJECT_STRUCT)
  #define objectplus(obj,offset)  as_object(as_oint(obj)+(soint)(offset))
#endif
%% #if !(defined(WIDE_SOFT) || defined(WIDE_AUXI) || defined(OBJECT_STRUCT))
%%   emit_define("objectplus(obj,offset)","((object)pointerplus(obj,offset))");
%% #elif defined(WIDE_AUXI)
%%   puts("static inline object objectplus (object obj, saint offset) { return (object){u:{both:{ one_ob: obj.one_o+offset, auxi_ob: obj.auxi_o }}}; }");
%% #else
%%   emit_define("objectplus(obj,offset)","as_object(as_oint(obj)+(soint)(offset))");
%% #endif

# Bit operations on entities of type uintV:
# ...vbit... instead of ...bit..., "v" = "value".
#if (intVsize > 32)
  #define vbit(n)  (LL(1)<<(n))
  #define vbitm(n)  (LL(2)<<((n)-1))
  #define vbit_test(x,n)  ((x) & vbit(n))
  #define minus_vbit(n)  (-LL(1)<<(n))
#else
  #define vbit  bit
  #define vbitm  bitm
  #define vbit_test  bit_test
  #define minus_vbit  minus_bit
#endif
%% #if notused
%%  export_def(vbit(n));
%%  export_def(vbitm(n));
%%  export_def(vbit_test(x,n));
%%  export_def(minus_vbit(n));
%% #endif

# Bit operations on entities of type oint:
# ...wbit... instead of ...bit..., "w" = "wide".
#if defined(WIDE_SOFT) || defined(WIDE_AUXI)
  #define wbit(n)  (LL(1)<<(n))
  #define wbitm(n)  (LL(2)<<((n)-1))
  #define wbit_test(x,n)  ((x) & wbit(n))
  #define minus_wbit(n)  (-LL(1)<<(n))
#else
  #define wbit  bit
  #define wbitm  bitm
  #define wbit_test  bit_test
  #define minus_wbit  minus_bit
#endif
%% export_def(wbit);
%% #if notused
%%  export_def(wbitm);
%% #endif
%% export_def(wbit_test);
%% export_def(minus_wbit);

#ifdef TYPECODES

  # Type info:
  # typecode(object) and mtypecode(object) yield the type code of
  # an object obj. For mtypecode it has to be in memory.
  #if !(exact_uint_size_p(oint_type_len) && (tint_type_mask == bit(oint_type_len)-1))
    #define typecode(expr)  \
      ((tint)(as_oint(expr) >> oint_type_shift) & (oint_type_mask >> oint_type_shift))
    #define mtypecode(expr)  typecode(expr)
  #else
    # The type 'tint' has exactly oint_type_len bits,
    # and tint_type_mask = 2^oint_type_len-1.
    # So it's not necessary for you to AND.
    # On the other hand on a 68000 a ROL.L #8 is faster,
    # as is a shift on a SPARC.
    #define typecode(expr)  ((tint)(as_oint(expr) >> oint_type_shift))
    #if defined(MC68000) && defined(GNU) && !defined(NO_ASM) && (oint_type_shift==24) && (oint_type_len==8)
      # GNU C on a 68000, replace LSR.L #24 with ROL.L #8 :
      #undef typecode
      #define typecode(expr)  \
        ({var tint __typecode;                                               \
          __asm__ ("roll #8,%0" : "=d" (__typecode) : "0" (as_oint(expr)) ); \
          __typecode;                                                        \
         })
      #elif defined(SPARC) && !defined(WIDE)
      #undef typecode
      #define typecode(expr)  \
        ((as_oint(expr) << (32-oint_type_len-oint_type_shift)) >> (32-oint_type_len))
    #elif defined(WIDE) && defined(WIDE_STRUCT)
      #undef typecode
      #define typecode(expr)  ((expr).u.both.type)
    #endif
    # Furthermore you can do accesses in memory without shift:
    #if !defined(WIDE) && (((oint_type_shift==24) && BIG_ENDIAN_P) || ((oint_type_shift==0) && !BIG_ENDIAN_P))
      #define mtypecode(expr)  (*(tint*)&(expr))
      #define fast_mtypecode
    #elif !defined(WIDE) && (((oint_type_shift==24) && !BIG_ENDIAN_P) || ((oint_type_shift==0) && BIG_ENDIAN_P))
      #define mtypecode(expr)  (*((tint*)&(expr)+3))
      #define fast_mtypecode
    #elif defined(WIDE)
      #ifdef WIDE_STRUCT
        #define mtypecode(expr)  ((expr).u.both.type)
      #elif (oint_type_len==16)
        #if (oint_type_shift==0) == BIG_ENDIAN_P
          #define mtypecode(expr)  (*((tint*)&(expr)+3))
        #else # (oint_type_shift==48) == BIG_ENDIAN_P
          #define mtypecode(expr)  (*(tint*)&(expr))
        #endif
      #elif (oint_type_len==32)
        #if (oint_type_shift==0) == BIG_ENDIAN_P
          #define mtypecode(expr)  (*((tint*)&(expr)+1))
        #else # (oint_type_shift==32) == BIG_ENDIAN_P
          #define mtypecode(expr)  (*(tint*)&(expr))
        #endif
      #endif
      #define fast_mtypecode
    #else # no optimization is possible
      #define mtypecode(expr)  typecode(expr)
    #endif
  #endif

  # Extraction of the address field without type info.
  # untype(obj)
  #if defined(WIDE) && defined(WIDE_STRUCT)
    #define untype(expr)  ((expr).u.both.addr)
  #elif !(defined(SPARC) && (oint_addr_len+oint_addr_shift<32))
    #define untype(expr)    \
      ((aint)(as_oint(expr) >> oint_addr_shift) & (aint)(oint_addr_mask >> oint_addr_shift))
  #else
    # On a SPARC processor long constants are slower than shifts:
    # Possibly, one does not need to use AND here.
    #define untype(expr)  \
      ((aint)((as_oint(expr) << (32-oint_addr_len-oint_addr_shift)) >> (32-oint_addr_len)))
  #endif

  # Object from type info and address field:
  # type_untype_object(type,address)
  #if defined(WIDE) && defined(WIDE_STRUCT)
    #if BIG_ENDIAN_P==WIDE_ENDIANNESS
      #define type_untype_object(type,address)  ((object){{(tint)(type),(aint)(address)}INIT_ALLOCSTAMP})
    #else
      #define type_untype_object(type,address)  ((object){{(aint)(address),(tint)(type)}INIT_ALLOCSTAMP})
    #endif
  #elif !(oint_addr_shift==0)
    #define type_untype_object(type,address)  \
      (as_object(  ((oint)(tint)(type) << oint_type_shift) + \
                   ((oint)(aint)(address) << oint_addr_shift) ))
  #else # you don't have to shift if oint_addr_shift=0:
    #if defined(WIDE_SOFT)
      # Beware: Conversion of  address  to oint by Zero-Extend!
      #define type_untype_object(type,address)              \
        objectplus((oint)(aint)(address),(oint)(tint)(type)<<oint_type_shift)
    #elif defined(WIDE_AUXI)
      #define type_untype_object(type,address)              \
        as_object_with_auxi((aint)pointerplus((address),(aint)(tint)(type)<<oint_type_shift))
    #elif defined(OBJECT_STRUCT)
      #define type_untype_object(type,address)              \
        as_object((oint)pointerplus((address),(oint)(tint)(type)<<oint_type_shift))
    #else # normal case
      # In order for this (NIL_IS_CONSTANT) to be a valid initializer
      # under gcc-2.5.8, you must not cast from pointer to oint and then
      # back to pointer, but you'll have to stay in the pointer's range..
      #define type_untype_object(type,address)              \
        as_object(pointerplus((address),(oint)(tint)(type)<<oint_type_shift))
    #endif
  #endif

  # Object from type info and direct data (as "address"):
  # type_data_object(type,data)
  #if defined(WIDE) && defined(WIDE_STRUCT)
    #if BIG_ENDIAN_P==WIDE_ENDIANNESS
      #define type_data_object(type,data)  ((object){{(tint)(type),(aint)(data)}INIT_ALLOCSTAMP})
    #else
      #define type_data_object(type,data)  ((object){{(aint)(data),(tint)(type)}INIT_ALLOCSTAMP})
    #endif
  #elif !(oint_addr_shift==0)
    #define type_data_object(type,data)  \
      (as_object(  ((oint)(tint)(type) << oint_type_shift) + \
                   ((oint)(aint)(data) << oint_addr_shift) ))
  #else # if oint_addr_shift=0, you don't have to shift:
    #define type_data_object(type,data)  \
      (as_object( ((oint)(tint)(type) << oint_type_shift) + (oint)(aint)(data) ))
  #endif

  # Extraction of the address without type info:
  # upointer(obj)
  # (upointer means "untyped pointer".)
  #if (addr_shift==0)
    #define upointer  untype
  #else
    #define optimized_upointer(obj)  \
      ((aint)((as_oint(obj) << (32-oint_addr_len-oint_addr_shift)) >> (32-oint_addr_len-addr_shift)))
    #define upointer(obj)  (untype(obj)<<addr_shift)
  #endif

  # Object from type info and address:
  # type_pointer_object(type,address)
  #if defined(WIDE_SOFT) && !defined(WIDE_STRUCT)
    # Cast to uintP, so that conversion of  address  to aint is done by Zero-Extend!
    #define type_pointer_object(type,address)  \
      type_untype_object(type,(aint)(uintP)(address)>>addr_shift)
  #elif (addr_shift==0)
    # (No cast to aint, so NIL can be used to initialize.)
    #define type_pointer_object(type,address)  \
      type_untype_object(type,address)
  #else # more efficient,
    # but this requires address to be divisible by 2^addr_shift:
    #define type_pointer_object(type,address)  \
      (as_object(((oint)(tint)(type) << oint_type_shift) + \
                 ((oint)(aint)(address) << (oint_addr_shift-addr_shift))))
  #endif

  # Object from constant type info and constant address:
  # type_constpointer_object(type,address)
  #define type_constpointer_object(type,address)  type_pointer_object(type,address)

  # oint from constant type info and address = 0:
  # type_zero_oint(type)
  #if defined(WIDE_SOFT) && defined(WIDE_STRUCT)
    #define type_zero_oint(type)  as_oint(type_untype_object(type,0))
  #else
    #define type_zero_oint(type)  ((oint)(tint)(type) << oint_type_shift)
  #endif

#else # HEAPCODES

  #ifdef STANDARD_HEAPCODES

    # We can assume a general alignment of 4 bytes, and thus have the low 2
    # bits for encoding type. Here's how we divide the address space:
    #   machine, frame_pointer  1/4
    #   subr                    1/4
    #   cons                    1/8
    #   varobject               1/4 (not 1/8 because symbol_tab is not 8-aligned)
    #   immediate               > 0 (anything >= 7/256 does it).
    # Note that cons and varobject cannot have the same encoding mod 8
    # (otherwise gc_mark:up wouldn't work).
    # So, here are the encodings.
    #   machine             ... .00   encodes pointers, offset 0
    #   subr                ... .10   encodes pointers, offset 2
    #   varobject           ... .01   offset 1, the pointers are == 0 mod 4
    #   cons                ... 011   offset 3, the pointers are == 0 mod 8
    #   immediate           ... 111
    #     fixnum            00s 111   s = sign bit
    #     sfloat            01s 111   s = sign bit
    #     char              100 111
    #     small-read-label  110 111
    #     system            111 111
    # Varobjects all start with a word containing the type (1 byte) and a
    # length field (up to 24 bits).

    /* These are the biases, mod 8. */
      #define machine_bias    0UL  /* mod 4 */
      #define subr_bias       2UL  /* mod 4 */
      #define varobject_bias  1UL  /* mod 4 */
      #define cons_bias       3UL  /* mod 8 */
      #define immediate_bias  7UL  /* mod 8 */

    # Immediate objects have a second type field.
      #if defined(LINUX_SPARC_OLD_GLIBC)
        #define imm_type_shift  29
      #else
        #define imm_type_shift  3
      #endif

    # The types of immediate objects.
      #define fixnum_type            ((0 << imm_type_shift) + immediate_bias)
      #define sfloat_type            ((2 << imm_type_shift) + immediate_bias)
      #define char_type              ((4 << imm_type_shift) + immediate_bias)
      #define small_read_label_type  ((6 << imm_type_shift) + immediate_bias)
      #define system_type            ((7 << imm_type_shift) + immediate_bias)

    # The sign bit, for immediate numbers only.
      #define sign_bit_t  (0 + imm_type_shift)
      #define sign_bit_o  (sign_bit_t+oint_type_shift)
    # Distinction between fixnums and bignums.
      #define bignum_bit_o  1
      #define NUMBER_BITS_INVERTED
    # Distinction between fixnums, short-floats and other kinds of numbers.
    # (NB: IMMEDIATE_FFLOAT is not defined for HEAPCODES.)
      #define number_immediatep(obj)  ((as_oint(obj) & wbit(1)) != 0)

    # For masking out the nonimmediate biases.
    # This must be 3, not 7, otherwise gc_mark won't work.
      #define nonimmediate_bias_mask  3
      #define nonimmediate_heapcode_mask  3

    # Combine an object from type info and immediate data.
    # type_data_object(type,data)
      #define type_data_object(type,data)  \
          (as_object(  ((oint)(tint)(type) << oint_type_shift) + \
                       ((oint)(aint)(data) << oint_data_shift) ))

    # An oint made up with a given type info, and address = 0.
    # type_zero_oint(type)
      #define type_zero_oint(type)  ((oint)(tint)(type) << oint_type_shift)

    # The GC bit. Addresses may not have this bit set.
      # define garcol_bit_o  (already defined above)  # only set during garbage collection

    # Test for immediate object.
    # immediate_object_p(obj)
      #define immediate_object_p(obj)  ((7 & ~as_oint(obj)) == 0)

    # Test for gc-invariant object. (This includes immediate, machine, subr.)
    # gcinvariant_object_p(obj)
      #define gcinvariant_object_p(obj)  \
        (((as_oint(obj) & 1) == 0) || immediate_object_p(obj))
      #define gcinvariant_oint_p(obj_o)  \
        ((((obj_o) & 1) == 0) || ((7 & ~(obj_o)) == 0))

    # Test for gc-invariant object, given only the bias.
      #define gcinvariant_bias_p(bias)  \
        ((((bias) & 1) == 0) || ((7 & ~(bias)) == 0))

    # The heap of a heap allocated object. 0 for varobjects, 1 for conses.
      #define nonimmediate_heapnr(obj)  \
        ((as_oint(obj) >> 1) & 1)

  #endif # STANDARD_HEAPCODES

  #ifdef LINUX_NOEXEC_HEAPCODES

    # We must assume a general alignment of 4 bytes and an enforced alignment
    # of 8 bytes for Lisp objects, and thus have the low 2 to 3 bits for
    # encoding heap and the garcol_bit. Here's how we divide the address space:
    #   machine, frame_pointer  1/4 * 3/4
    #   immediate               1/4 * 1/4
    #   cons                    1/8
    #   varobject               1/8
    # Note that cons and varobject cannot have the same encoding mod 8
    # (otherwise gc_mark:up wouldn't work).
    # Immediates look like pointers in the range 0xC0000000..0xFFFFFFFF.
    # We know that the Linux kernel never assigns virtual memory in this area.
    # So, here are the encodings. Bit 0 is used as the garcol_bit.
    #   machine                ... ... .00   encodes pointers, offset 0
    #   cons                   ... ... 010   offset 2, the pointers are == 0 mod 8
    #   varobject              ... ... 110   offset 6, the pointers are == 4 mod 8
    #   immediate           11 ... ...  00
    #     fixnum            11 ... 00s  00   s = sign bit
    #     sfloat            11 ... 01s  00   s = sign bit
    #     char              11 ... 100  00
    #     small-read-label  11 ... 110  00
    #     system            11 ... 111  00
    # Varobjects all start with a word containing the type (1 byte) and a
    # length field (up to 24 bits).

    # These are the biases.
      #define machine_bias    0  # + 0 mod 4
      #define varobject_bias  2  # + 4 mod 8
      #define cons_bias       2  # + 0 mod 8
      #define immediate_bias  0xC0000000  # + 0 mod 4
      #define subr_bias       varobject_bias

    # Immediate objects have a second type field.
      #define imm_type_shift  2  # could also be 3, if oint_data_shift == 6

    # The types of immediate objects.
      #define fixnum_type            ((0 << imm_type_shift) + immediate_bias)
      #define sfloat_type            ((2 << imm_type_shift) + immediate_bias)
      #define char_type              ((4 << imm_type_shift) + immediate_bias)
      #define small_read_label_type  ((6 << imm_type_shift) + immediate_bias)
      #define system_type            ((7 << imm_type_shift) + immediate_bias)

    # The sign bit, for immediate numbers only.
      #define sign_bit_t  (0 + imm_type_shift)
      #define sign_bit_o  (sign_bit_t+oint_type_shift)
    # Distinction between fixnums and bignums.
      #define bignum_bit_o  1
    # Distinction between fixnums, short-floats and other kinds of numbers.
    # (NB: IMMEDIATE_FFLOAT is not defined for HEAPCODES.)
      #define number_immediatep(obj)  ((as_oint(obj) & wbit(1)) == 0)

    # The misalignment of varobjects, modulo varobject_alignment.
      #define varobjects_misaligned  4

    # For masking out the nonimmediate biases.
      #define nonimmediate_bias_mask  3
      #define nonimmediate_heapcode_mask  7

    # Combine an object from type info and immediate data.
    # type_data_object(type,data)
      #define type_data_object(type,data)  \
          (as_object(  ((oint)(tint)(type) << oint_type_shift) + \
                       ((oint)(aint)(data) << oint_data_shift) ))

    # An oint made up with a given type info, and address = 0.
    # type_zero_oint(type)
      #define type_zero_oint(type)  ((oint)(tint)(type) << oint_type_shift)

    # The GC bit. Addresses may not have this bit set.
      # define garcol_bit_o  (already defined above)  # only set during garbage collection

    # Test for immediate object.
    # immediate_object_p(obj)
      #define immediate_object_p(obj)  \
        ((as_oint(obj) & 0xE0000003) == (immediate_bias & 0xE0000003))

    # Test for gc-invariant object. (This includes immediate, machine.)
    # gcinvariant_object_p(obj)
      #define gcinvariant_object_p(obj)  \
        ((as_oint(obj) & bit(1)) == 0)
      #define gcinvariant_oint_p(obj_o)  \
        (((obj_o) & bit(1)) == 0)
    # NB: Subrs are not included in this test, because subrp(obj) require a
    # memory access.

    # Test for gc-invariant object, given only the bias.
      #define gcinvariant_bias_p(bias)  \
        (((bias) & 2) == 0)

    # The heap of a heap allocated object. 0 for varobjects, 1 for conses.
      #define nonimmediate_heapnr(obj)  \
        (1 & ~(as_oint(obj) >> 2))

  #endif # LINUX_NOEXEC_HEAPCODES

#endif # TYPECODES
%% #ifdef TYPECODES
%%  export_def(typecode(expr));
%%  export_def(mtypecode(expr));
%%  export_def(type_untype_object(type,address));
%%  export_def(upointer(obj));
%%  export_def(type_pointer_object(type,address));
%%  export_def(type_constpointer_object(type,address));
%% #else
%%  export_def(number_immediatep(obj));
%%  export_def(immediate_object_p(obj));
%%  export_def(gcinvariant_object_p(obj));
%%  export_def(gcinvariant_oint_p(obj_o));
%%  export_def(gcinvariant_bias_p(bias));
%% #endif
%% export_def(type_data_object(type,address));
%% export_def(type_zero_oint(obj));

# The misalignment of varobjects, modulo varobject_alignment.
#ifndef varobjects_misaligned
  #define varobjects_misaligned  0
#endif
#if varobjects_misaligned
  #define VAROBJECTS_ALIGNMENT_DUMMY_DECL  char alignment_dummy[varobjects_misaligned];
#else
  #define VAROBJECTS_ALIGNMENT_DUMMY_DECL
#endif
%% export_def(varobjects_misaligned);
%% export_def(VAROBJECTS_ALIGNMENT_DUMMY_DECL);

# The misalignment of conses, modulo 2*sizeof(gcv_object_t).
#ifndef conses_misaligned
  #define conses_misaligned  0
#endif


# Objects with variable length must reside at addresses that are divisable by 2
#if defined(VAX) # ?? gcc/config/vax/vax.h sagt: Alignment = 4
  #define varobject_alignment  1
#endif
#if defined(MC680X0)
  #if addr_shift!=0
    #define varobject_alignment  bit(addr_shift)  # because of the condensed distribution of typecodes
  #else
    #define varobject_alignment  2
  #endif
#endif
#if defined(I80386) || defined(POWERPC) || defined(ARM) || defined(S390)
  #define varobject_alignment  4
#endif
#if defined(SPARC) || defined(HPPA) || defined(MIPS) || defined(M88000) || defined(DECALPHA) || defined(IA64) || defined(AMD64)
  #define varobject_alignment  8
#endif
#if (!defined(TYPECODES) || defined(GENERATIONAL_GC)) && (varobject_alignment < 4)
  #undef varobject_alignment
  #define varobject_alignment  4
#endif
#if ((defined(GENERATIONAL_GC) && defined(WIDE)) || defined(LINUX_NOEXEC_HEAPCODES)) && (varobject_alignment < 8)
  #undef varobject_alignment
  #define varobject_alignment  8
#endif
# varobject_alignment should be defined:
#ifndef varobject_alignment
  #error "varobject_alignment depends on CPU -- readjust varobject_alignment!!"
#endif
# varobject_alignment should be a power of 2:
#if !((varobject_alignment & (varobject_alignment-1)) ==0)
  #error "Bogus varobject_alignment -- readjust varobject_alignment!!"
#endif
# varobject_alignment should be a multiple of 2^addr_shift :
#if (varobject_alignment % bit(addr_shift))
  #error "Bogus varobject_alignment -- readjust varobject_alignment!!"
#endif
%% export_def(varobject_alignment);


#ifdef TYPECODES

# Now we'll define the various type bits and type codes.

# Single-floats can be immediate objects, like short-floats, if there are
# enough bits in a 'gcv_object_t'.
#if defined(WIDE_HARD) || defined(WIDE_SOFT)
  #define IMMEDIATE_FFLOAT
#endif

# Determine whether a type isn't changed by the GC
# (ie. if it's not a pointer):
  #if 0 && (defined(GNU) || defined(INTEL))
    #define gcinvariant_type_p(type)  \
      ({var bool _erg;                         \
        switch (type)                          \
          { case_machine:                      \
            case_char: case_subr: case_system: \
            case_fixnum: case_sfloat:          \
            /* with WIDE also: case_ffloat: */ \
              _erg = true; break;              \
            default: _erg = false; break;      \
          }                                    \
        _erg;                                  \
       })
  #endif

#ifndef tint_allowed_type_mask
  #define tint_allowed_type_mask  tint_type_mask
#endif

# There are 7 to 8 type bits available: TB7, [TB6,] TB5, TB4, ..., TB0.
# All of them have to be set in tint_allowed_type_mask and thus in tint_type_mask as well
# We distribute them under the assumption that only one bit is missing in tint_type_mask.
# TB6 will be set to -1, if it can't be used.
#if ((0xFF & ~tint_allowed_type_mask) == 0)
  #define TB7 7
  #define TB6 6
  #define TB5 5
  #define TB4 4
  #define TB3 3
  #define TB2 2
  #define TB1 1
  #define TB0 0
#elif (oint_type_len==7)
  #define TB7 6
  #define TB6 -1
  #define TB5 5
  #define TB4 4
  #define TB3 3
  #define TB2 2
  #define TB1 1
  #define TB0 0
#else
  # Some bits have to be avoided
  #define tint_avoid  ((bitm(oint_type_len)-1) & ~tint_allowed_type_mask)
  # tint_avoid must only contain one bit:
  #if (tint_avoid & (tint_avoid-1))
    #error "Bogus oint_type_mask -- oint_type_mask has more than one extraneous bit!!"
  #endif
  # tint_avoid consists of exactly one bit that has to be avoided.
  #if (tint_avoid > bit(0))
    #define TB0 0
  #else
    #define TB0 1
  #endif
  #if (tint_avoid > bit(1))
    #define TB1 1
  #else
    #define TB1 2
  #endif
  #if (tint_avoid > bit(2))
    #define TB2 2
  #else
    #define TB2 3
  #endif
  #if (tint_avoid > bit(3))
    #define TB3 3
  #else
    #define TB3 4
  #endif
  #if (tint_avoid > bit(4))
    #define TB4 4
  #else
    #define TB4 5
  #endif
  #if (tint_avoid > bit(5))
    #define TB5 5
  #else
    #define TB5 6
  #endif
  #if ((tint_allowed_type_mask & ~0xFF) == 0)
    #define TB6 -1
    #if (tint_avoid > bit(6))
      #define TB7 6
    #else
      #define TB7 7
    #endif
  #else
    #if (tint_avoid > bit(6))
      #define TB6 6
    #else
      #define TB6 7
    #endif
    #if (tint_avoid > bit(7))
      #define TB7 7
    #else
      #define TB7 8
    #endif
  #endif
#endif

# bit masks for the type bits:
  #define BTB0  bit(TB0)
  #define BTB1  bit(TB1)
  #define BTB2  bit(TB2)
  #define BTB3  bit(TB3)
  #define BTB4  bit(TB4)
  #define BTB5  bit(TB5)
  #define BTB6  bit(TB6)
  #define BTB7  bit(TB7)

#define STANDARD_8BIT_TYPECODES

#ifdef STANDARD_8BIT_TYPECODES

#if defined(I80386) && defined(UNIX_LINUX) && (CODE_ADDRESS_RANGE == 0)
  # At 0x60000000 there are the shared-libraries.
  # At 0x50000000 (Linux 1.2) resp. 0x40000000 (Linux 2.0) there are several
  # mmap-pages,for example ones allocated  by setlocale() or gettext().
  # Therefore we only have to do a few changes to the distribution of the type codes.
#endif

#if defined(I80386) && defined(UNIX_LINUX) && (CODE_ADDRESS_RANGE != 0)
  # Code and malloc memory is at 0x08000000.
  # Therefore avoid allocating typecode 0x08 for the moment.
#endif

#if (defined(MC680X0) || (defined(SPARC) && !defined(SUN4_29))) && defined(UNIX_LINUX)
  # At 0x50000000 there are shared libraries located.
  # But this doesn't mean we have to change the type code distribution.
#endif

#if (defined(MIPS) || defined(POWERPC)) && defined(UNIX_LINUX)
  # At 0x2AAAB000 there are shared libraries located.
  # But this doesn't mean we have to change the type code distribution.
#endif

#if defined(DECALPHA) && defined(UNIX_OSF) && !(defined(NO_SINGLEMAP) || defined(NO_TRIVIALMAP))
# mmap() only works with addresses >=0, <2^38, but since ordinary pointers are in the range
# 1*2^32..2*2^32, only the Bits 37..33 remain as type-bits.
#endif

#if defined(SPARC64) && defined(UNIX_LINUX)
  # At 0x70000000 there are shared libraries located.
  # But this doesn't mean we have to change the type code distribution.
#endif

# Type bits:
# in Typcodes (tint):
  #define garcol_bit_t     TB7  # only set during GC
  #if (TB6 >= 0)
    #define cons_bit_t     TB6  # only set for CONS
  #endif
  #define number_bit_t     TB5  # only set for numbers
  #define notsimple_bit_t  TB3  # for arrays: deleted for simple arrays
  #define sign_bit_t       TB0  # Sign for real numbers (set <==> number <0)
  #define float_bit_t      TB1
  #define float1_bit_t     TB3
  #define float2_bit_t     TB2
  #define ratio_bit_t      TB3
  #define bignum_bit_t     TB2
# in Objects (oint):
  #define garcol_bit_o     (garcol_bit_t+oint_type_shift)    # only set during the garbage collection!
  #if (TB6 >= 0)
    #define cons_bit_o     (cons_bit_t+oint_type_shift)      # only set for cons CONS
  #endif
  #define number_bit_o     (number_bit_t+oint_type_shift)    # only set for numbers
  #define notsimple_bit_o  (notsimple_bit_t+oint_type_shift) # for arrays: deleted for simple arrays
  #define sign_bit_o       (sign_bit_t+oint_type_shift)      # Sign for real numbers
  #define float_bit_o      (float_bit_t+oint_type_shift)
  #define float1_bit_o     (float1_bit_t+oint_type_shift)
  #define float2_bit_o     (float2_bit_t+oint_type_shift)
  #define ratio_bit_o      (ratio_bit_t+oint_type_shift)
  #define bignum_bit_o     (bignum_bit_t+oint_type_shift)

# constant type codes:
  #define machine_type    (0)                                  # 0x00  # %00000000  ; machine pointer
  #define subr_type       (                              BTB0) # 0x01  # %00000001  ; SUBR
  #define char_type       (                         BTB1     ) # 0x02  # %00000010  ; character
  #define system_type     (                         BTB1|BTB0) # 0x03  # %00000011  ; frame-pointer, small-read-label, system
  #define symbol_type     (                    BTB2          ) # 0x04  # %000001xx  ; symbol
          # bits for symbols in the GCself pointer:
          #define var_bit0_t  TB0  # set if the symbol is proclaimed SPECIAL or constant
          #define var_bit1_t  TB1  # set if the symbol is a symbol-macro or constant
  #if (TB6 < 0)
  #define cons_type       (               BTB3               ) # 0x08  # %00001000  ; cons
  #endif
  #define closure_type    (               BTB3          |BTB0) # 0x09  # %00001001  ; closure
  #define structure_type  (               BTB3     |BTB1     ) # 0x0A  # %00001010  ; structure
  #define stream_type     (               BTB3     |BTB1|BTB0) # 0x0B  # %00001011  ; stream
  #define orecord_type    (               BTB3|BTB2          ) # 0x0C  # %00001100  ; OtherRecord (Package, Byte, ...)
  #define instance_type   (               BTB3|BTB2     |BTB0) # 0x0D  # %00001101  ; CLOS instance
  #define lrecord_type    (               BTB3|BTB2|BTB1     ) # 0x0E  # %00001110  ; LongRecord (WeakList, WeakAlist, ...)
  #define mdarray_type    (               BTB3|BTB2|BTB1|BTB0) # 0x0F  # %00001111  ; other array (rank/=1 or other eltype)
  #define sbvector_type   (          BTB4                    ) # 0x10  # %00010000  ; simple-bit-vector
  #define sb2vector_type  (          BTB4               |BTB0) # 0x11  # %00010001  ; simple (VECTOR (UNSIGNED-BYTE 2))
  #define sb4vector_type  (          BTB4          |BTB1     ) # 0x12  # %00010010  ; simple (VECTOR (UNSIGNED-BYTE 4))
  #define sb8vector_type  (          BTB4          |BTB1|BTB0) # 0x13  # %00010011  ; simple (VECTOR (UNSIGNED-BYTE 8))
  #define sb16vector_type (          BTB4     |BTB2          ) # 0x14  # %00010100  ; simple (VECTOR (UNSIGNED-BYTE 16))
  #define sb32vector_type (          BTB4     |BTB2     |BTB0) # 0x15  # %00010101  ; simple (VECTOR (UNSIGNED-BYTE 32))
  #define sstring_type    (          BTB4     |BTB2|BTB1     ) # 0x16  # %00010110  ; simple-string
  #define svector_type    (          BTB4     |BTB2|BTB1|BTB0) # 0x17  # %00010111  ; simple-vector
  #define bvector_type    (          BTB4|BTB3               ) # 0x18  # %00011000  ; non-simple bit-vector
  #define b2vector_type   (          BTB4|BTB3          |BTB0) # 0x19  # %00011001  ; non-simple (VECTOR (UNSIGNED-BYTE 2))
  #define b4vector_type   (          BTB4|BTB3     |BTB1     ) # 0x1A  # %00011010  ; non-simple (VECTOR (UNSIGNED-BYTE 4))
  #define b8vector_type   (          BTB4|BTB3     |BTB1|BTB0) # 0x1B  # %00011011  ; non-simple (VECTOR (UNSIGNED-BYTE 8))
  #define b16vector_type  (          BTB4|BTB3|BTB2          ) # 0x1C  # %00011100  ; non-simple (VECTOR (UNSIGNED-BYTE 16))
  #define b32vector_type  (          BTB4|BTB3|BTB2     |BTB0) # 0x1D  # %00011101  ; non-simple (VECTOR (UNSIGNED-BYTE 32))
  #define string_type     (          BTB4|BTB3|BTB2|BTB1     ) # 0x1E  # %00011110  ; non-simple string
  #define vector_type     (          BTB4|BTB3|BTB2|BTB1|BTB0) # 0x1F  # %00011111  ; non-simple (VECTOR T)
  #define fixnum_type     (     BTB5                         ) # 0x20  # %00100000  ; fixnum
  #define sfloat_type     (     BTB5               |BTB1     ) # 0x22  # %00100010  ; short-float
  #define bignum_type     (     BTB5          |BTB2          ) # 0x24  # %00100100  ; bignum
  #define ffloat_type     (     BTB5          |BTB2|BTB1     ) # 0x26  # %00100110  ; single-float
  #define ratio_type      (     BTB5     |BTB3               ) # 0x28  # %00101000  ; ratio
  #define dfloat_type     (     BTB5     |BTB3     |BTB1     ) # 0x2A  # %00101010  ; double-float
  #define complex_type    (     BTB5     |BTB3|BTB2          ) # 0x2C  # %00101100  ; complex
  #define lfloat_type     (     BTB5     |BTB3|BTB2|BTB1     ) # 0x2E  # %00101110  ; long-float
  #if (TB6 >= 0)
  #define cons_type       (BTB6                              ) # 0x40  # %01000000  ; cons
  #endif

# Bits for symbols in VAR/FUN-Frames (in LISP-Stack):
# aren't in the oint_type-part, but in the oint_addr-part.
  #define active_bit  0  # set: binding is active
  #define dynam_bit   1  # set: binding is dynamic
  #define svar_bit    2  # set: next parameter is supplied-p-parameter for this
#if (varobject_alignment >= bit(3))
  #define oint_symbolflags_shift  oint_addr_shift
#else
  #define NO_symbolflags # there's no space in the symbol for active_bit, dynam_bit, svar_bit
#endif

#ifndef IMMEDIATE_FFLOAT
  # type is GC-invariant, if
  # type-info-byte >=0, <= system_type or >= fixnum_type, < bignum_type.
    #define gcinvariant_type_p(type)  \
      (((type) & ~(BTB5|BTB1|BTB0)) == 0)
#else
  # type is GC-invariant, if
  # type-info-byte is one of 0x00..0x03,0x20..0x23,0x26..0x27 ist.
  #if (TB1==TB0+1) && (TB2==TB0+2) && (TB3==TB0+3) && (TB4==TB0+4) && (TB5==TB0+5)
    #define gcinvariant_type_p(type)  \
      ((((type)>>(TB0+1))<0x14) && ((bit((type)>>(TB0+1)) & 0xFFF4FFFCUL) == 0))
  #else
    # Test whether ((type)>>TB1) is one of
    #   0, 1, bit(TB5-TB1), bit(TB5-TB1) | 1, bit(TB5-TB1) | bit(TB2-TB1) | 1.
    #define gcinvariant_type_p(type)  gcinvariant_type_aux((type)>>TB1)
    #define gcinvariant_type_sum(type)  \
      (((type) | ((type)>>(TB5-(TB2+1)))) & (((BTB2<<1)+BTB2+BTB1)>>TB1))
    #define gcinvariant_type_aux(type)                         \
      (((type) < ((BTB5+(BTB2<<1))>>TB1))                      \
       && ((type) & ~((BTB5|BTB2|BTB1)>>TB1)) == 0             \
       && (bit(gcinvariant_type_sum(type))                     \
           & (  bit(0)                                         \
              | bit(1)                                         \
              | bit(bit(TB2+1-TB1))                            \
              | bit(bit(TB2+1-TB1) | 1)                        \
              | bit(bit(TB2+1-TB1) | bit(TB2-TB1) | 1))) != 0)
  #endif
#endif

#endif # STANDARD_8BIT_TYPECODES

#if !(gcinvariant_type_p(ffloat_type) == defined(IMMEDIATE_FFLOAT))
  #error "gcinvariant_type_p() incorrectly implemented!"
#endif

# Test for gc-invariant object. (This includes immediate, machine, subr.)
# gcinvariant_object_p(obj)
  #define gcinvariant_object_p(obj)  \
    gcinvariant_type_p(typecode(obj))
  #define gcinvariant_oint_p(obj_o)  \
    gcinvariant_type_p((tint)((obj_o) >> oint_type_shift) & (oint_type_mask >> oint_type_shift))

#else # no TYPECODES

# Bits for symbols in VAR/FUN-Frames (on LISP-Stack):
# are not located in the oint_type-part, but in the oint_data-part.
  #define active_bit  0  # set: binding is active
  #define dynam_bit   1  # set: binding is dynamic
  #define svar_bit    2  # set: next parameter is supplied-p-parameter for this one
#define NO_symbolflags # there's no space in the symbol for active_bit, dynam_bit, svar_bit

# Bits for symbols in the flags:
  #define var_bit0_f  0  # set if the symbol is proclaimed SPECIAL or constant
  #define var_bit1_f  1  # set if the symbol is a symbol-macro or constant

#endif # TYPECODES
%% #ifdef TYPECODES
%%  export_def(gcinvariant_type_p(type));
%%  export_def(gcinvariant_type_sum(type));
%%  export_def(gcinvariant_type_aux(type));
%%  export_def(gcinvariant_object_p(obj));
%%  export_def(gcinvariant_oint_p(obj_o));
%% #endif


# What's really being sent from an address to the address-bus
#if defined(MC68000)
  #define hardware_  0x00FFFFFFUL  # 68000 drops 8
#elif defined(SUN3) && !defined(UNIX_SUNOS4)
  #define hardware_addressbus_mask  0x0FFFFFFFUL  # SUN3 unter SunOS 3.5 wirft 4 Bits weg
#else
  #define hardware_addressbus_mask  ~0UL  # Default: nothing is dropped
#endif
# Clever memory-mapping spares us from masking out of certain
# bits before one accesses the address
#define addressbus_mask  hardware_addressbus_mask
#ifdef MAP_MEMORY
  #if defined(SUN4_29)
    # Memory-mapping makes the bits 28..24 of an address redundant now.
    #undef addressbus_mask
    #define addressbus_mask  0xE0FFFFFFUL
  #elif defined(DECALPHA) && defined(UNIX_OSF)
    # Memory-mapping makes the bits 39..33 of an address redundant now.
    #undef addressbus_mask
    #define addressbus_mask  0xFFFFFF01FFFFFFFFUL
  #elif !defined(WIDE_SOFT)
    # Memory-mapping makes the bits 31..24 of an address redundant now.
    #undef addressbus_mask
    #define addressbus_mask  oint_addr_mask  # most of the time it's = 0x00FFFFFFUL
  #endif
#endif
%% #if notused
%% export_def(addressbus_mask);
%% #endif


#ifdef TYPECODES

# You have to remove the typebits in order to access the components
# of an object.
# pointable(obj) does this, returning a void*.
# pointable_unchecked(obj) likewise, but without the DEBUG_GCSAFETY check.
# pointable_address_unchecked(obj_o) likewise, but takes an oint as argument
#                                    and returns an aint.
  #if !((oint_addr_shift==0) && (addr_shift==0))
    #define pointable_unchecked(obj)  ((void*)upointer(obj))
    #define pointable_address_unchecked(obj_o)  \
      (((aint)((obj_o) >> oint_addr_shift) & (aint)(oint_addr_mask >> oint_addr_shift)) << addr_shift)
  #else
    #define pointable_unchecked(obj)  \
      ((void*)pointable_address_unchecked(as_oint(obj)))
    # If oint_addr_shift=0 and addr_shift=0, you don't have to shift.
    #if !((tint_type_mask & (addressbus_mask>>oint_type_shift)) == 0)
      #define pointable_address_unchecked(obj_o)  \
        ((aint)(obj_o) & ((aint)oint_addr_mask | ~addressbus_mask))
    #else
      # Moreover if oint_type_mask and addressbus_mask are disjoint
      # (((tint_type_mask<<oint_type_shift) & addressbus_mask) == 0),
      # no typebits are being sent to the address bus anyway.
      # So there's nothing to be done:
      #define pointable_address_unchecked(obj_o)  (aint)(obj_o)
    #endif
  #endif
  #ifdef DEBUG_GCSAFETY
    # Check that obj has not been held in a GC-unsafe variable while a
    # memory allocation was made.
    static inline void* pointable (gcv_object_t obj) {
      return pointable_unchecked(obj);
    }
    static inline void* pointable (object obj) {
      return pointable_unchecked((gcv_object_t)obj); # The cast does the check.
    }
  #else
    #define pointable(obj)  pointable_unchecked(obj)
  #endif

# If you want to access an object with a known type-info whose
# set typebits are being swallowed by the address bus (the
# typebits, that are =0 don't matter), you can do without 'untype':
  #if defined(DEBUG_GCSAFETY)
    #define type_pointable(type,obj)  pointable(obj)
  #elif defined(WIDE_STRUCT)
    #define type_pointable(type,obj)  ((void*)((obj).u.both.addr))
  #elif !((oint_addr_shift==0) && (addr_shift==0) && (((tint_type_mask<<oint_type_shift) & addressbus_mask) == 0))
    #if (addr_shift==0)
      #define type_pointable(type,obj)  \
        ((oint_addr_shift==0) && ((type_zero_oint(type) & addressbus_mask) == 0) \
         ? (void*)(aint)as_oint(obj)                                             \
         : (void*)(aint)pointable(obj)                                           \
        )
    #elif !(addr_shift==0)
      # Analogous, but here the macro 'optimized_upointer'
      # assumes the role of the address bus:
      #define type_pointable(type,obj)  \
        ((optimized_upointer(type_data_object(type,0)) == 0) \
         ? (void*)(aint)optimized_upointer(obj)              \
         : (void*)(aint)pointable(obj)                       \
        )
    #endif
  #else
    # If pointable(obj) = obj, type_pointable() doesn't do anything as well:
    #define type_pointable(type,obj)  ((void*)(aint)as_oint(obj))
  #endif

# If you want to access an object that has one of several known
# type infos, you can probably omit the 'untype'.
# The  OR of the type infos is more authoritative.
  #define types_pointable(ORed_types,obj)  type_pointable(ORed_types,obj)

#else # HEAPCODES

  #define pointable_address_unchecked(obj_o)  \
    (((aint)((obj_o) >> oint_addr_shift) & (aint)(oint_addr_mask >> oint_addr_shift)) << addr_shift)

#endif
%% #ifdef TYPECODES
%%  export_def(pointable_unchecked(obj));
%%  export_def(pointable_address_unchecked(obj_o));
%%  #ifdef DEBUG_GCSAFETY
%%   puts("static inline void* pointable (gcv_object_t obj) { return pointable_unchecked(obj); }");
%%   puts("static inline void* pointable (object obj) { return pointable_unchecked((gcv_object_t)obj); }");
%%  #else
%%   emit_define("pointable(obj)","pointable_unchecked(obj)");
%%  #endif
%% #endif


#if defined(SINGLEMAP_MEMORY) && (((system_type*1UL << oint_type_shift) & addressbus_mask) == 0)
  # The STACK resides in a singlemap-area as well, Typinfo system_type.
  #define SINGLEMAP_MEMORY_STACK
#endif


#ifdef oint_symbolflags_shift
  #if defined(SINGLEMAP_MEMORY) && (oint_symbolflags_shift==oint_type_shift)
    # Since we can't multimap the symbol_tab, we can't use extrabits in
    # a symbol's typecode.
    #undef oint_symbolflags_shift
    #define NO_symbolflags
  #endif
#endif
#ifdef NO_symbolflags
  #define oint_symbolflags_shift  -1 # invalid value
#endif


# Whether we try to initialize subr_tab statically.
# (g++ 3.3 doesn't accept compound expressions as initializers: PR#12615.
# g++ 3.4 similarly: PR#15180.)
#if !(defined(WIDE_SOFT) && !defined(WIDE_STRUCT)) && !(defined(__GNUG__) && (__GNUC__ == 3) && (__GNUC_MINOR__ == 3 || __GNUC_MINOR__ == 4) && defined(OBJECT_STRUCT))
  #define INIT_SUBR_TAB
#endif
# NB: This has to be defined so external modules can work.
# When changed: do nothing

# Whether we try to initialize symbol_tab statically.
# (Make initialization easier, but there is not enough space for the
# compilation of SPVWTABS on some systems.
# EMX 0.9c (gcc-2.7.2.1) says "Virtual memory exhausted".
# g++ 3.3 doesn't accept compound expressions as initializers: PR#12615.
# g++ 3.4 similarly: PR#15180.)
#if !(defined(WIDE_SOFT) && !defined(WIDE_STRUCT)) && !(defined(__GNUG__) && (__GNUC__ == 3) && (__GNUC_MINOR__ == 3 || __GNUC_MINOR__ == 4) && defined(OBJECT_STRUCT))
  #define INIT_SYMBOL_TAB
#endif
# When changed: nothing to do

# Whether we try to initialize object_tab statically.
# (g++ 3.3 doesn't accept compound expressions as initializers: PR#12615.
# g++ 3.4 similarly: PR#15180.)
#if !(defined(WIDE_SOFT) && !defined(WIDE_STRUCT)) && !(defined(__GNUG__) && (__GNUC__ == 3) && (__GNUC_MINOR__ == 3 || __GNUC_MINOR__ == 4) && defined(OBJECT_STRUCT))
  #define INIT_OBJECT_TAB
#endif
# When changed: do nothing


# Set during the core of GC.
# When this is set, unexpected handle_fault() calls that can
# - if defined(MORRIS_GC) && defined(GENERATIONAL_GC) - copy
# Morris-chain backpointers from a cons cell to an old_new_pointer_t with set
# garcol_bit(!) into the heap, where they are guaranteed to lead to a crash
# later. So, uncontrolled memory accesses are forbidden while inside_gc.
extern bool inside_gc;
%% puts("extern bool inside_gc;");


#ifdef DEBUG_GCSAFETY

  # Forward declarations.
  static inline bool gcinvariant_symbol_p (object obj);
  #ifdef LINUX_NOEXEC_HEAPCODES
  static inline bool nonimmsubrp (object obj);
  #else
  #define nonimmsubrp(obj)  false
  #endif

  # Force a crash if a memory pointer points to nonexistent memory.
  #define nonimmprobe(obj_o)  \
    do {                                                                       \
      /* Don't do probes inside GC. It leads to unexpected handle_fault()      \
         calls that can - if defined(MORRIS_GC) && defined(GENERATIONAL_GC) -  \
         copy Morris-chain backpointers from a cons cell to an old_new_pointer_t \
         with set garcol_bit(!) into the heap, where they are guaranteed to    \
         lead to a crash later. */                                             \
      if (!inside_gc)                                                          \
        if (((obj_o) & wbit(garcol_bit_o)) == 0) /* exclude frame words from the STACK */ \
          if (!gcinvariant_oint_p(obj_o)) /* exclude immediate objects */      \
            /* Access a single char, without needing to subtract the bias. */  \
            *(volatile char *)pointable_address_unchecked(obj_o);              \
    } while (0)

  # When a gcv_object_t is fetched from a GC visible location (in the heap or
  # on the STACK) we can assume that GC has updated it.
  inline gcv_object_t::operator object () const {
    nonimmprobe(one_o);
    return (object){ one_o: one_o INIT_ALLOCSTAMP };
  }

  # When an object is put into a GC visible location (in the heap or
  # on the STACK) we check that it has not been held in a GC-unsafe variable
  # while a memory allocation was made.
  inline gcv_object_t::gcv_object_t (object obj) {
    if (!(gcinvariant_object_p(obj) || gcinvariant_symbol_p(obj)
          || obj.allocstamp == alloccount || nonimmsubrp(obj)))
      abort();
    one_o = as_oint(obj);
    nonimmprobe(one_o);
  }
  # The only exception are fake gcv_objects.
  inline gcv_object_t::gcv_object_t (fake_gcv_object obj) {
    one_o = obj.fake_value;
  }

  # Uninitialized.
  inline gcv_object_t::gcv_object_t () {
  }

  # Start of an area where no GC can be triggered.
  inline gcunsafe_object_t::gcunsafe_object_t (object obj)
    : gcv_object_t (obj), allocstamp (alloccount) {}
  inline gcunsafe_object_t::gcunsafe_object_t (gcv_object_t obj)
    : gcv_object_t (obj), allocstamp (alloccount) {}
  # End of an area where no GC can be triggered.
  inline gcunsafe_object_t::~gcunsafe_object_t () {
    if (!(allocstamp == alloccount))
      abort();
  }

#endif
%% #ifdef DEBUG_GCSAFETY
%%   puts("static inline bool gcinvariant_symbol_p (object obj);");
%%   #ifdef LINUX_NOEXEC_HEAPCODES
%%     puts("static inline bool nonimmsubrp (object obj);");
%%   #else
%%     emit_define("nonimmsubrp(obj)","false");
%%   #endif
%%   export_def(nonimmprobe(obj_o));
%%   puts("inline gcv_object_t::operator object () const { nonimmprobe(one_o); return (object){ one_o: one_o, allocstamp: alloccount }; }");
%%   puts("inline gcv_object_t::gcv_object_t (object obj) { if (!(gcinvariant_object_p(obj) || gcinvariant_symbol_p(obj) || obj.allocstamp == alloccount || nonimmsubrp(obj))) abort(); one_o = as_oint(obj); nonimmprobe(one_o); }");
%%   puts("inline gcv_object_t::gcv_object_t () {}");
%% #endif


/* Force a memory allocation for all functions which can trigger GC but
   sometimes do and sometimes don't. This makes DEBUG_GCSAFETY more efficient.
   GCTRIGGER()             does a no-op memory allocation
   GCTRIGGER1(obj1)        likewise, but saves obj1 temporarily
   GCTRIGGER2(obj1,obj2)   likewise, but saves obj1, obj2 temporarily
   ...
   GCTRIGGER_IF(condition, statement)
                           likewise, but only if the condition is fulfilled */
#ifdef DEBUG_GCSAFETY
  # When these macros are used in C macros, obj1, obj2 etc. can sometimes be
  # expressions of type 'object' and sometimes of type 'gcv_object_t'. Need
  # two implementations of inc_allocstamp.
  inline void inc_allocstamp (object& obj) { obj.allocstamp++; }
  inline void inc_allocstamp (gcv_object_t& obj) { }
  #define GCTRIGGER()  \
    (void)(alloccount++)
  #define GCTRIGGER1(obj1)  \
    (void)(inc_allocstamp(obj1), alloccount++)
  #define GCTRIGGER2(obj1,obj2)  \
    (void)(inc_allocstamp(obj1), inc_allocstamp(obj2), alloccount++)
  #define GCTRIGGER3(obj1,obj2,obj3)  \
    (void)(inc_allocstamp(obj1), inc_allocstamp(obj2), inc_allocstamp(obj3), alloccount++)
  #define GCTRIGGER4(obj1,obj2,obj3,obj4)  \
    (void)(inc_allocstamp(obj1), inc_allocstamp(obj2), inc_allocstamp(obj3), inc_allocstamp(obj4), alloccount++)
  #define GCTRIGGER5(obj1,obj2,obj3,obj4,obj5)  \
    (void)(inc_allocstamp(obj1), inc_allocstamp(obj2), inc_allocstamp(obj3), inc_allocstamp(obj4), inc_allocstamp(obj5), alloccount++)
  #define GCTRIGGER6(obj1,obj2,obj3,obj4,obj5,obj6)  \
    (void)(inc_allocstamp(obj1), inc_allocstamp(obj2), inc_allocstamp(obj3), inc_allocstamp(obj4), inc_allocstamp(obj5), inc_allocstamp(obj6), alloccount++)
  #define GCTRIGGER_IF(condition,statement)  \
    if (condition) statement
#else
  #define GCTRIGGER()  (void)0
  #define GCTRIGGER1(obj1)  (void)0
  #define GCTRIGGER2(obj1,obj2)  (void)0
  #define GCTRIGGER3(obj1,obj2,obj3)  (void)0
  #define GCTRIGGER4(obj1,obj2,obj3,obj4)  (void)0
  #define GCTRIGGER5(obj1,obj2,obj3,obj4,obj5)  (void)0
  #define GCTRIGGER6(obj1,obj2,obj3,obj4,obj5,obj6)  (void)0
  #define GCTRIGGER_IF(condition,statement)  (void)0
#endif


# ################### Methods for memory management ###################### #

# SPVW_BLOCKS : Memory management with few memory blocks
# SPVW_PAGES  : Memory management with many memory blocks
# SPVW_MIXED  : Objects of mixed types are possible on the same page or block
# SPVW_PURE   : Every memory block/every memory page contains only objects
#               of exactly one type
#if defined(MAP_MEMORY) || defined(TRIVIALMAP_MEMORY)
  # Multimapping of single pages isn't implemented yet.??
  # Singlemapping of single pages isn't implemented yet.??
  # If you use mmap() as malloc()-replacement, single pages aren't needed.
  #define SPVW_BLOCKS
#elif defined(VIRTUAL_MEMORY)
  # On Unix-systems you can still fetch more memory afterwards,
  # but you should concentrate the data - if possible - on few pages.
  #define SPVW_PAGES
#else
  #define SPVW_BLOCKS
#endif
#if defined(MULTIMAP_MEMORY)
  # MULTIMAP_MEMORY -> Mixed pages allow a better usage of memory.
  #define SPVW_MIXED
#elif defined(SINGLEMAP_MEMORY)
  # SINGLEMAP_MEMORY -> Ony pure pages/blocks make sense, since
  # the address of a page determines the type of the objects it contains.
  #define SPVW_PURE
#elif !defined(TYPECODES) || defined(MC68000) || defined(SUN3) || defined(SPVW_BLOCKS) || defined(TRIVIALMAP_MEMORY)
  # !TYPECODES -> there aren't real typecodes, only Cons and Varobject.
  # MC68000 or SUN3 -> type_pointable(...) costs little or nothing.
  # SPVW_BLOCKS -> SPVW_PURE_BLOCKS is only implemented for SINGLEMAP_MEMORY.
  # TRIVIALMAP_MEMORY -> not many blocks available, small adress space.
  #define SPVW_MIXED
#elif 1 # provisionally!??
  #define SPVW_MIXED
#endif
#if !(defined(SPVW_BLOCKS) || defined(SPVW_PAGES))
  #error "readjust SPVW_BLOCKS/SPVW_PAGES!"
#endif
#if !(defined(SPVW_MIXED) || defined(SPVW_PURE))
  #error "readjust SPVW_MIXED/SPVW_PURE!"
#endif
#if (defined(SPVW_BLOCKS) && defined(SPVW_PURE)) != defined(SINGLEMAP_MEMORY)
  #error "SINGLEMAP_MEMORY <==> SPVW_PURE_BLOCKS!"
#endif
#if (defined(SPVW_BLOCKS) && defined(SPVW_MIXED)) < defined(TRIVIALMAP_MEMORY)
  #error "TRIVIALMAP_MEMORY ==> SPVW_MIXED_BLOCKS!"
#endif
#if defined(SPVW_PURE) && !defined(TYPECODES)
  #error "SPVW_PURE ==> TYPECODES!"
#endif
#if (defined(SPVW_BLOCKS) && (defined(SPVW_PURE) || defined(SPVW_MIXED))) < defined(GENERATIONAL_GC)
  #error "GENERATIONAL_GC ==> SPVW_PURE_BLOCKS or SPVW_MIXED_BLOCKS_STAGGERED or SPVW_MIXED_BLOCKS_OPPOSITE!"
#endif

# Algorithm by Morris, that compacts Conses without mixing them up:
#if defined(SPVW_BLOCKS) && defined(VIRTUAL_MEMORY) && !defined(NO_MORRIS_GC)
  # Morris-GC is recommended, as it preserves the locality.
  #define MORRIS_GC
#endif

# Put subr_tab and symbol_tab to given addresses through memory-mapping.
# (The Morris-GC uses the macro upointer() for MULTIMAP_MEMORY. For
# &symbol_tab = 0x20000000 it'd be upointer(NIL)=0. Darn!)
#if defined(MAP_MEMORY) && !defined(WIDE_SOFT) && !(defined(MULTIMAP_MEMORY) && defined(MORRIS_GC))
  #define MAP_MEMORY_TABLES
#endif


# ################# definitions by cases with respect to type codes ################## #

#ifdef TYPECODES

# Has to start with switch (typecode(obj)), after that it's like a
# switch-statement with arbitrarily many case-labels.
# Example:  switch (typecode(arg)) { case_string: ...; break; ... }
  #define case_machine    case machine_type   # machine-pointer
  #define case_sstring    case sstring_type   # Simple-String
  #define case_ostring    case string_type    # Other String
  #define case_sbvector   case sbvector_type   # Simple-Bit-Vector
  #define case_obvector   case bvector_type    # Other Bit-Vector
  #define case_sb2vector  case sb2vector_type  # Simple-2Bit-Vector
  #define case_ob2vector  case b2vector_type   # Other 2Bit-Vector
  #define case_sb4vector  case sb4vector_type  # Simple-4Bit-Vector
  #define case_ob4vector  case b4vector_type   # Other 4Bit-Vector
  #define case_sb8vector  case sb8vector_type  # Simple-8Bit-Vector
  #define case_ob8vector  case b8vector_type   # Other 8Bit-Vector
  #define case_sb16vector case sb16vector_type # Simple-16Bit-Vector
  #define case_ob16vector case b16vector_type  # Other 16Bit-Vector
  #define case_sb32vector case sb32vector_type # Simple-32Bit-Vector
  #define case_ob32vector case b32vector_type  # Other 32Bit-Vector
  #define case_svector    case svector_type    # Simple-(General-)Vector
  #define case_ovector    case vector_type    # Other (General-)Vector
  #define case_mdarray    case mdarray_type   # other Array
  #define case_string     case_sstring: case_ostring # general string
  #define case_bvector    case_sbvector: case_obvector # general bit vector
  #define case_b2vector   case_sb2vector: case_ob2vector # general 2bit vector
  #define case_b4vector   case_sb4vector: case_ob4vector # general 4bit vector
  #define case_b8vector   case_sb8vector: case_ob8vector # general 8bit vector
  #define case_b16vector  case_sb16vector: case_ob16vector # general 16bit vector
  #define case_b32vector  case_sb32vector: case_ob32vector # general 32bit vector
  #define case_vector     case_svector: case_ovector # general vector
  #define case_array      case_string: case_bvector: case_b2vector: case_b4vector: case_b8vector: case_b16vector: case_b32vector: case_vector: case_mdarray # general Array
  #define case_closure    case closure_type   # Closure
  #ifdef structure_type
  #define case_structure  case structure_type # Structure
  #define _case_structure case_structure:
  #else
  #define structure_type  orecord_type        # Structures are OtherRecords
  #define _case_structure
  #endif
  #ifdef stream_type
  #define case_stream     case stream_type    # Stream
  #define _case_stream    case_stream:
  #else
  #define stream_type     orecord_type        # Streams are OtherRecords
  #define _case_stream
  #endif
  #define case_orecord    case orecord_type   # Other Record
  #define case_instance   case instance_type  # CLOS-Instance
  #define case_lrecord    case lrecord_type   # Long Record
  #define case_char       case char_type      # Character
  #define case_subr       case subr_type      # SUBR
  #define case_system     case system_type    # Frame-Pointer, Small-Read-Label, System
  #define case_posfixnum  case fixnum_type    # Fixnum >=0
  #define case_negfixnum  case fixnum_type|bit(sign_bit_t) # Fixnum <0
  #define case_fixnum     case_posfixnum: case_negfixnum # Fixnum
  #define case_posbignum  case bignum_type    # Bignum >0
  #define case_negbignum  case bignum_type|bit(sign_bit_t) # Bignum <0
  #define case_bignum     case_posbignum: case_negbignum # Bignum
  #define case_integer    case_fixnum: case_bignum # Integer
  #define case_ratio      case ratio_type: case ratio_type|bit(sign_bit_t) # Ratio
  #ifdef SPVW_MIXED
  #define _case_ratio     case_ratio:
  #else
  #define _case_ratio
  #endif
  #define case_rational   case_integer: case_ratio # Rational
  #define case_sfloat     case sfloat_type: case sfloat_type|bit(sign_bit_t) # Short-Float
  #define case_ffloat     case ffloat_type: case ffloat_type|bit(sign_bit_t) # Single-Float
  #define case_dfloat     case dfloat_type: case dfloat_type|bit(sign_bit_t) # Double-Float
  #define case_lfloat     case lfloat_type: case lfloat_type|bit(sign_bit_t) # Long-Float
  #define case_float      case_sfloat: case_ffloat: case_dfloat: case_lfloat # Float
  #define case_real       case_rational: case_float # Real
  #define case_complex    case complex_type # Complex
  #ifdef SPVW_MIXED
  #define _case_complex   case_complex:
  #else
  #define _case_complex
  #endif
  #define case_number     case_real: case_complex # Number
  #define case_symbol     case symbol_type # Symbol
  #define case_sxrecord   case_closure: _case_structure _case_stream _case_ratio _case_complex case_orecord: case_instance # Srecord/Xrecord general
  #define case_record     case_sxrecord: case_lrecord # Lrecord/Srecord/Xrecord general
  #if /* !defined(NO_symbolflags) && */ (oint_symbolflags_shift==oint_type_shift)
  #define case_symbolflagged  # Symbol with Flags                       \
          case symbol_type:                                             \
          case symbol_type|bit(active_bit):                             \
          case symbol_type|bit(dynam_bit):                              \
          case symbol_type|bit(dynam_bit)|bit(active_bit):              \
          case symbol_type|bit(svar_bit):                               \
          case symbol_type|bit(svar_bit)|bit(active_bit):               \
          case symbol_type|bit(svar_bit)|bit(dynam_bit):                \
          case symbol_type|bit(svar_bit)|bit(dynam_bit)|bit(active_bit)
  #else
  #define case_symbolflagged  case_symbol # Symbol with flags
  #endif
  #define case_cons       case cons_type # Cons

#else

  #define _case_structure
  #define _case_stream

#endif


# ################## Structure of memory of LISP objects #################### #

# uintWC is the Integer type for the lengths of Bignum, Lfloat, Iarray.
# Subset relation: uintW <= uintWC <= uintC.
#ifdef TYPECODES
  #define intWCsize intCsize
  typedef uintC uintWC;
  typedef sintC sintWC;
#else
  # Type and sign are stored in the heap - only 16 bits for the length.
  #define intWCsize intWsize
  typedef uintW uintWC;
  typedef sintW sintWC;
#endif
# uintWCoverflow(x) checks, if there has been an overflow after the execution
# of an x++.
#define uintWCoverflow(x)  ((intWCsize<intLsize) && ((uintWC)(x)==0))

# ---------------------- Objects with two pointers ---------------------- #
# They contain just the two pointers, no header. The type must already be
# known when the object is accessed.

# Normally, Cons, Ratio, Complex can all be considered as pairs. But if
# SPVW_MIXED, the heap statistics are a little unspecific if we mix the
# three types; therefore in that case we let Ratio and Complex be Varobjects.
#ifdef SPVW_MIXED
  #define case_pair  case_cons
#else
  #define case_pair  case_cons: case_ratio: case_complex
#endif

# ---------------------- Objects of varying length ---------------------- #
# The first word is reserved for garbage collection. Outside of garbage
# collection, it contains a pointer to the object itself. Note that the
# GC, when it moves an object, takes care not to modify the typecode of
# this first word (except the GC bit, which it temporarily uses).

# Type of the header flags:
#if (oint_type_len<=8) && !defined(ARM) && !defined(DECALPHA) && !defined(IA64) && !defined(DEBUG_GCSAFETY)
  # Access to an individual byte is possible
  #define hfintsize  intBsize
  typedef uintB  hfint;
#else
  # access to a full word
  #define hfintsize  pointer_bitsize
  typedef uintP  hfint;
#endif
%% #if (oint_type_len<=8) && !defined(ARM) && !defined(DECALPHA) && !defined(IA64) && !defined(DEBUG_GCSAFETY)
%%   emit_typedef("uintB","hfint");
%% #else
%%   emit_typedef("uintP","hfint");
%% #endif

# Objecs with variable length
#ifdef TYPECODES
  #ifdef DEBUG_GCSAFETY
    #define VAROBJECT_HEADER  \
               gcv_object_t _GCself;  # Self pointer for GC, contains flags
  #else
    #define VAROBJECT_HEADER  \
               union {                                                    \
                 gcv_object_t _GCself;  # Self pointer for GC             \
                 hfint flags[sizeof(gcv_object_t)/sizeof(hfint)]; # Flags \
               } header;
  #endif
#else
  #define VAROBJECT_HEADER  \
               gcv_object_t GCself;  # Self pointer for GC \
               uintL tfl;            # type, flags, length
#endif
typedef struct {
  VAROBJECT_HEADER
} varobject_;
typedef varobject_ *  Varobject;
#ifdef TYPECODES
  #ifdef DEBUG_GCSAFETY
    #define GCself  _GCself
    #define header_flags  _GCself.one_o
  #else
    #define GCself  header._GCself
    # The typecode can be found in the byte ((Varobject)p)->header_flags.
    #if !(oint_type_len>=hfintsize ? oint_type_shift%hfintsize==0 : floor(oint_type_shift,hfintsize)==floor(oint_type_shift+oint_type_len-1,hfintsize))
      #error "Bogus header_flags -- redefine header_flags!"
    #endif
    #if BIG_ENDIAN_P
      #define header_flags  header.flags[sizeof(gcv_object_t)/sizeof(hfint)-1-floor(oint_type_shift,hfintsize)]
    #else
      #define header_flags  header.flags[floor(oint_type_shift,hfintsize)]
    #endif
  #endif
  # it applies  mtypecode(((Varobject)p)->GCself) =
  # (((Varobject)p)->header_flags >> (oint_type_shift%hfintsize)) & tint_type_mask
  # Bits for Symbols in the self pointer (see above):
  # define var_bit0_t  ...  # set if the symbol is proclaimed SPECIAL or constant
  # define var_bit1_t  ...  # set if the symbol is a symbol-macro or constant
  #define var_bit0_hf  (var_bit0_t+(oint_type_shift%hfintsize))
  #define var_bit1_hf  (var_bit1_t+(oint_type_shift%hfintsize))
#else
  # Three possible layouts of type, flags, length:
  #   8 bits type, 24 bits length [Vrecord, Lrecord]
  #   8 bits type, 8 bits flags, 16 bits length [Srecord]
  #   8 bits type, 8 bits flags, 8 bits length, 8 bits xlength [Xrecord]
  #define vrecord_tfl(type,length)  \
    ((uintL)(uintB)(type)+((uintL)(length)<<8))
  #define lrecord_tfl(type,length)  \
    ((uintL)(uintB)(type)+((uintL)(length)<<8))
  #define srecord_tfl(type,flags,length)  \
    ((uintL)(uintB)(type)+((uintL)(uintB)(flags)<<8)+((uintL)(length)<<16))
  #define xrecord_tfl(type,flags,length,xlength)  \
    ((uintL)(uintB)(type)+((uintL)(uintB)(flags)<<8)+((uintL)(uintB)(length)<<16)+((uintL)(uintB)(xlength)<<24))
  #define varobject_type(ptr) ((sintB)((ptr)->tfl & 0xFF))
  #if defined(__GNUC__) && (__GNUC__ == 2) && ((__GNUC_MINOR__ == 8) || (__GNUC_MINOR__ == 90))
    # Work around for a gcc bug present (at least) in gcc-2.8.1 on hppa and
    # egcs-1.0.3a on i386. It miscompiles xpathnamep.
    #undef varobject_type
    #define varobject_type(ptr) ((sintB)((sintL)((ptr)->tfl) & 0xFF))
  #endif
    # Bits for symbols in the flags:
    #define header_flags  tfl
    #define var_bit0_hf  (var_bit0_f+8)
    #define var_bit1_hf  (var_bit1_f+8)
#endif
%% export_def(VAROBJECT_HEADER);
%% #ifndef TYPECODES
%%  export_def(GCself);
%%  export_def(varobject_type(ptr));
%% #endif

# Records
# These are varobjects with a one-byte type field in memory.
# There are three types of records:
#   Vector-Records can have up to 16777215 elements, but have no flags and
#   if TYPECODES also no type (because the type info is in the pointer).
#   Long-Records can have up to 16777215 elements, but have no flags.
#   Simple-Records can have up to 65535 elements,
#   Extended-Records have room for up to 255 elements and 255 extra (non-Lisp)
#   elements.
# Vector-Records are recognized by their type field:
#   rectype == Rectype_Sbvector, Rectype_Sb[2|4|8|16|32]vector,
#              Rectype_S[8|16|32]string, Rectype_Imm_S[8|16|32]string,
#              Rectype_Svector.
# Long-Records are recognized by rectype >= rectype_longlimit, or if TYPECODES
# equivalently by their typecode lrecord_type.
# The others are partitioned into:
#   - Simple-Records, if rectype < rectype_limit.
#   - Extended-Records, if rectype >= rectype_limit.

#ifdef TYPECODES
  #define RECORD_HEADER                                               \
    VAROBJECT_HEADER   /* self-pointer GC */                          \
    sintB rectype;     /* for OtherRecord and LongRecord: sub-type */ \
    uintB recflags;    /* for OtherRecord: flags */                   \
    uintW reclength;   /* lengths and others */
#else
  #define RECORD_HEADER  \
    VAROBJECT_HEADER   /* self-pointer for GC, tfl */
#endif
typedef struct {
  RECORD_HEADER
  gcv_object_t recdata[unspecified] _attribute_aligned_object_; # elements
} record_;
typedef record_ *  Record;
# access to type, flags:
#ifdef TYPECODES
  #define record_type(ptr)  ((ptr)->rectype)
#else
  #define record_type(ptr)  varobject_type(ptr)
#endif
#define Record_type(obj)  record_type(TheRecord(obj))
#ifdef TYPECODES
  #define record_flags(ptr)  ((ptr)->recflags)
#else
  #define record_flags(ptr)  (((ptr)->tfl >> 8) & 0xFF)
#endif
#define Record_flags(obj)  record_flags(TheRecord(obj))
#ifdef TYPECODES
  #define record_flags_clr(ptr,bits)  ((ptr)->recflags &= ~(bits))
  #define record_flags_set(ptr,bits)  ((ptr)->recflags |= (bits))
  #define record_flags_replace(ptr,newflags)  ((ptr)->recflags = (newflags))
#else
  #define record_flags_clr(ptr,bits)  ((ptr)->tfl &= ~((uintL)(bits) << 8))
  #define record_flags_set(ptr,bits)  ((ptr)->tfl |= ((uintL)(bits) << 8))
  #define record_flags_replace(ptr,newflags)  \
    ((ptr)->tfl ^= (((ptr)->tfl ^ (uintL)(newflags)<<8) & 0xFF00))
#endif
%% export_def(RECORD_HEADER);
%% sprintf(buf,"struct { RECORD_HEADER gcv_object_t recdata[unspecified]%s; }",attribute_aligned_object);
%% emit_typedef(buf,"record_");
%% emit_typedef("record_ *","Record");
%% export_def(record_type(ptr));
%% export_def(Record_type(obj));
%% export_def(record_flags(ptr));
%% export_def(record_flags_set(ptr,bits));
%% export_def(Record_flags(obj));

#ifdef TYPECODES
  #define VRECORD_HEADER  \
                 VAROBJECT_HEADER # self-pointer for GC \
                 uintL length;    # length
#else
  #define VRECORD_HEADER  \
                 VAROBJECT_HEADER # self-pointer for GC, tfl
#endif
typedef struct {
  VRECORD_HEADER
} vrecord_;
typedef vrecord_ *  Vrecord;
#ifdef TYPECODES
  #define vrecord_length(ptr)  ((ptr)->length)
#else
  #define vrecord_length(ptr)  ((ptr)->tfl >> 8)
#endif
%% export_def(VRECORD_HEADER);
%% emit_typedef("struct { VRECORD_HEADER }","vrecord_");
%% emit_typedef("vrecord_ *","Vrecord");
%% export_def(vrecord_length(ptr));

#ifdef TYPECODES
  #define LRECORD_HEADER  \
                 VAROBJECT_HEADER # self-pointer for GC \
                 uintL tfl;       # subtype (1 byte), then length (3 bytes)
#else
  #define LRECORD_HEADER  \
                 VAROBJECT_HEADER # self-pointer for GC, tfl
#endif
typedef struct {
  LRECORD_HEADER
  gcv_object_t recdata[unspecified] _attribute_aligned_object_; # reclength elements
} lrecord_;
typedef lrecord_ *  Lrecord;
#ifdef TYPECODES
  #if BIG_ENDIAN_P
    #define lrecord_tfl(type,length)  \
      (((uintL)(uintB)(type)<<24)+(uintL)(length))
    #define lrecord_length(ptr)  ((ptr)->tfl & 0xFFFFFF)
  #else
    #define lrecord_tfl(type,length)  \
      ((uintL)(uintB)(type)+((uintL)(length)<<8))
    #define lrecord_length(ptr)  ((ptr)->tfl >> 8)
  #endif
#else
  #define lrecord_length(ptr)  ((ptr)->tfl >> 8)
#endif
#define Lrecord_length(obj)  lrecord_length(TheLrecord(obj))

#ifdef TYPECODES
  #define SRECORD_HEADER                                        \
                 VAROBJECT_HEADER # self-pointer GC             \
                 sintB rectype;   # subtype, < rectype_limit    \
                 uintB recflags;  # flags                       \
                 uintW reclength; # lengths in objects
#else
  #define SRECORD_HEADER  \
                 VAROBJECT_HEADER # self-pointer for GC, tfl
#endif
typedef struct {
  SRECORD_HEADER
  gcv_object_t recdata[unspecified] _attribute_aligned_object_; # reclength elements
} srecord_;
typedef srecord_ *  Srecord;
#ifdef TYPECODES
  #define srecord_length(ptr)  ((ptr)->reclength)
#else
  #define srecord_length(ptr)  ((ptr)->tfl >> 16)
#endif
#define Srecord_length(obj)  srecord_length(TheSrecord(obj))
%% export_def(SRECORD_HEADER);
%% sprintf(buf,"struct { SRECORD_HEADER gcv_object_t recdata[unspecified]%s; }",attribute_aligned_object);
%% emit_typedef(buf,"srecord_");
%% emit_typedef("srecord_ *","Srecord");
%% export_def(srecord_length(ptr));

#ifdef TYPECODES
  #define XRECORD_HEADER                                                \
                 VAROBJECT_HEADER  # self-pointer for GC                \
                 sintB rectype;    # subtype, >= rectype_limit          \
                 uintB recflags;   # flags                              \
                 uintB reclength;  # lengths in objects                 \
                 uintB recxlength; # lengths of the extra objects
#else
  #define XRECORD_HEADER  \
                 VAROBJECT_HEADER  # self-pointer for GC, tfl
#endif
typedef struct {
  XRECORD_HEADER
  gcv_object_t recdata[unspecified] _attribute_aligned_object_; # reclength elements
  # uintB      recxdata[unspecified]; # recxlength extra elements
} xrecord_;
typedef xrecord_ *  Xrecord;
#ifdef TYPECODES
  #define xrecord_length(ptr)  ((ptr)->reclength)
  #define xrecord_xlength(ptr)  ((ptr)->recxlength)
#else
  #define xrecord_length(ptr)  (((ptr)->tfl >> 16) & 0xFF)
  #define xrecord_xlength(ptr)  ((ptr)->tfl >> 24)
#endif
#define Xrecord_length(obj)  xrecord_length(TheXrecord(obj))
#define Xrecord_xlength(obj)  xrecord_xlength(TheXrecord(obj))
%% export_def(XRECORD_HEADER);
%% #if notused
%% sprintf(buf,"struct { XRECORD_HEADER gcv_object_t recdata[unspecified]%s; }",attribute_aligned_object);
%% emit_typedef(buf,"xrecord_");
%% emit_typedef("xrecord_ *","Xrecord");
%% #endif

/* *** Possible rectype values for records. *** */
enum {
  enum_rectype_first = -4, /* Try to keep rectype_limit = 0. */
  Rectype_Closure,
%% printf("#define Rectype_Closure %d\n",Rectype_Closure);
  Rectype_Structure,     /* only used #ifndef case_structure */
%% printf("#define Rectype_Structure %d\n",Rectype_Structure);
  Rectype_Instance,
%% printf("#define Rectype_Instance %d\n",Rectype_Instance);
     rectype_limit, /* Here is the limit between Srecord and Xrecord. */
  Rectype_Hashtable = rectype_limit,
%% printf("#define Rectype_Hashtable %d\n",Rectype_Hashtable);
#ifndef TYPECODES
%% #ifndef TYPECODES
  /* Rectype_vector is the bottom ARRAY & VECTOR */
  Rectype_vector,       /* 1 -- Iarray, not Srecord/Xrecord */
%% printf("#define Rectype_vector %d\n",Rectype_vector);
  Rectype_bvector,      /* 2 -- Iarray, not Srecord/Xrecord */
%% printf("#define Rectype_bvector %d\n",Rectype_bvector);
  Rectype_b2vector,     /* 3 -- Iarray, not Srecord/Xrecord */
%% printf("#define Rectype_b2vector %d\n",Rectype_b2vector);
  Rectype_b4vector,     /* 4 -- Iarray, not Srecord/Xrecord */
%% printf("#define Rectype_b4vector %d\n",Rectype_b4vector);
  Rectype_b8vector,     /* 5 -- Iarray, not Srecord/Xrecord */
%% printf("#define Rectype_b8vector %d\n",Rectype_b8vector);
  Rectype_b16vector,    /* 6 -- Iarray, not Srecord/Xrecord */
%% printf("#define Rectype_b16vector %d\n",Rectype_b16vector);
  Rectype_b32vector,    /* 7 -- Iarray, not Srecord/Xrecord */
%% printf("#define Rectype_b32vector %d\n",Rectype_b32vector);
  Rectype_unused1,           /* 8 */
  /* Rectype_Svector is the bottom SIMPLE VECTOR */
  Rectype_Svector,     /* 9 -- Svector, not Srecord/Xrecord */
%% printf("#define Rectype_Svector %d\n",Rectype_Svector);
  Rectype_Sbvector,   /* 10 -- Sbvector, not Srecord/Xrecord */
%% printf("#define Rectype_Sbvector %d\n",Rectype_Sbvector);
  Rectype_Sb2vector,  /* 11 -- Sbvector, not Srecord/Xrecord */
%% printf("#define Rectype_Sb2vector %d\n",Rectype_Sb2vector);
  Rectype_Sb4vector,  /* 12 -- Sbvector, not Srecord/Xrecord */
%% printf("#define Rectype_Sb4vector %d\n",Rectype_Sb4vector);
  Rectype_Sb8vector,  /* 13 -- Sbvector, not Srecord/Xrecord */
%% printf("#define Rectype_Sb8vector %d\n",Rectype_Sb8vector);
  Rectype_Sb16vector, /* 14 -- Sbvector, not Srecord/Xrecord */
%% printf("#define Rectype_Sb16vector %d\n",Rectype_Sb16vector);
  Rectype_Sb32vector, /* 15 -- Sbvector, not Srecord/Xrecord */
%% printf("#define Rectype_Sb32vector %d\n",Rectype_Sb32vector);
  Rectype_unused2,          /* 16 */
  /* Rectype_S8string is the bottom STRING */
  Rectype_S8string,   /* 17 -- S8string, not Srecord/Xrecord */
%% printf("#define Rectype_S8string %d\n",Rectype_S8string);
  Rectype_Imm_S8string, /* 18 -- immutable S8string, not Srecord/Xrecord */
%% printf("#define Rectype_Imm_S8string %d\n",Rectype_Imm_S8string);
  Rectype_S16string, /* 19 -- S16string, not Srecord/Xrecord */
%% printf("#define Rectype_S16string %d\n",Rectype_S16string);
  Rectype_Imm_S16string, /* 20 -- immutable S16string, not Srecord/Xrecord */
%% printf("#define Rectype_Imm_S16string %d\n",Rectype_Imm_S16string);
  Rectype_S32string, /* 21 -- S32string, not Srecord/Xrecord */
%% printf("#define Rectype_S32string %d\n",Rectype_S32string);
  Rectype_Imm_S32string, /* 22 -- immutable S32string, not Srecord/Xrecord */
%% printf("#define Rectype_Imm_S32string %d\n",Rectype_Imm_S32string);
  Rectype_reallocstring, /* 23 -- reallocated simple string, an Sistring, only used #ifdef HAVE_SMALL_SSTRING */
%% printf("#define Rectype_reallocstring %d\n",Rectype_reallocstring);
  /* Rectype_reallocstring is the top SIMPLE-STRING & SIMPLE-VECTOR */
  Rectype_string,      /* 24 -- Iarray, not Srecord/Xrecord */
%% printf("#define Rectype_string %d\n",Rectype_string);
  /* Rectype_string is the top STRING */
  Rectype_mdarray,     /* 25 -- Iarray, not Srecord/Xrecord */
%% printf("#define Rectype_mdarray %d\n",Rectype_mdarray);
  /* Rectype_mdarray is the top ARRAY */
  /* Rectype_Bignum is the bottom NUMBER */
  Rectype_Bignum,        /* Bignum, not Srecord/Xrecord */
%% printf("#define Rectype_Bignum %d\n",Rectype_Bignum);
  Rectype_Lfloat,        /* Lfloat, not Srecord/Xrecord */
%% printf("#define Rectype_Lfloat %d\n",Rectype_Lfloat);
  Rectype_Dfloat,
%% printf("#define Rectype_Dfloat %d\n",Rectype_Dfloat);
  Rectype_Ffloat,
%% printf("#define Rectype_Ffloat %d\n",Rectype_Ffloat);
%% #endif
#endif  /* TYPECODES */
#ifdef SPVW_MIXED
%% #ifdef SPVW_MIXED
  Rectype_Ratio,
%% printf("#define Rectype_Ratio %d\n",Rectype_Ratio);
  Rectype_Complex,
%% printf("#define Rectype_Complex %d\n",Rectype_Complex);
%% #endif
#endif /* SPVW_MIXED */
  /* *** Here the numbers end. *** */
#ifndef TYPECODES
%% #ifndef TYPECODES
  Rectype_Symbol,
%% printf("#define Rectype_Symbol %d\n",Rectype_Symbol);
%% #endif
#endif  /* TYPECODES */
  Rectype_Package,
%% printf("#define Rectype_Package %d\n",Rectype_Package);
  Rectype_Readtable,
%% printf("#define Rectype_Readtable %d\n",Rectype_Readtable);
  Rectype_Pathname,
%% printf("#define Rectype_Pathname %d\n",Rectype_Pathname);
#ifdef LOGICAL_PATHNAMES
%% #ifdef LOGICAL_PATHNAMES
  Rectype_Logpathname,
%% printf("#define Rectype_Logpathname %d\n",Rectype_Logpathname);
%% #endif
#endif  /* LOGICAL_PATHNAMES */
  Rectype_Random_State,
%% printf("#define Rectype_Random_State %d\n",Rectype_Random_State);
#ifndef case_stream
%% #ifndef case_stream
  Rectype_Stream,
%% printf("#define Rectype_Stream %d\n",Rectype_Stream);
%% #endif
#endif  /* case_stream */
  Rectype_Byte,
%% printf("#define Rectype_Byte %d\n",Rectype_Byte);
  Rectype_Subr,
%% printf("#define Rectype_Subr %d\n",Rectype_Subr);
  Rectype_Fsubr,
%% printf("#define Rectype_Fsubr %d\n",Rectype_Fsubr);
  Rectype_Loadtimeeval,
%% printf("#define Rectype_Loadtimeeval %d\n",Rectype_Loadtimeeval);
  Rectype_Symbolmacro,
%% printf("#define Rectype_Symbolmacro %d\n",Rectype_Symbolmacro);
  Rectype_GlobalSymbolmacro,
%% printf("#define Rectype_GlobalSymbolmacro %d\n",Rectype_GlobalSymbolmacro);
  Rectype_Macro,
%% printf("#define Rectype_Macro %d\n",Rectype_Macro);
  Rectype_FunctionMacro,
%% printf("#define Rectype_FunctionMacro %d\n",Rectype_FunctionMacro);
  Rectype_BigReadLabel,
%% printf("#define Rectype_BigReadLabel %d\n",Rectype_BigReadLabel);
  Rectype_Encoding,
%% printf("#define Rectype_Encoding %d\n",Rectype_Encoding);
  Rectype_Fpointer,      /* only used #ifdef FOREIGN */
%% printf("#define Rectype_Fpointer %d\n",Rectype_Fpointer);
#ifdef DYNAMIC_FFI
%% #ifdef DYNAMIC_FFI
  Rectype_Faddress,
%% printf("#define Rectype_Faddress %d\n",Rectype_Faddress);
  Rectype_Fvariable,
%% printf("#define Rectype_Fvariable %d\n",Rectype_Fvariable);
  Rectype_Ffunction,
%% printf("#define Rectype_Ffunction %d\n",Rectype_Ffunction);
%% #endif
#endif  /* DYNAMIC_FFI */
  Rectype_Weakpointer,
%% printf("#define Rectype_Weakpointer %d\n",Rectype_Weakpointer);
  Rectype_MutableWeakList,
%% printf("#define Rectype_MutableWeakList %d\n",Rectype_MutableWeakList);
  Rectype_MutableWeakAlist,
%% printf("#define Rectype_MutableWeakAlist %d\n",Rectype_MutableWeakAlist);
  Rectype_Weakmapping,
%% printf("#define Rectype_Weakmapping %d\n",Rectype_Weakmapping);
  Rectype_Finalizer,
%% printf("#define Rectype_Finalizer %d\n",Rectype_Finalizer);
#ifdef SOCKET_STREAMS
%% #ifdef SOCKET_STREAMS
  Rectype_Socket_Server,
%% printf("#define Rectype_Socket_Server %d\n",Rectype_Socket_Server);
%% #endif
#endif  /* SOCKET_STREAMS */
#ifdef YET_ANOTHER_RECORD
%% #ifdef YET_ANOTHER_RECORD
  Rectype_Yetanother,
%% printf("#define Rectype_Yetanother %d\n",Rectype_Yetanother);
%% #endif
#endif /* YET_ANOTHER_RECORD */
  rectype_longlimit, /* the boundary between Srecord/Xrecord and Lrecord. */
  Rectype_WeakList,
%% printf("#define Rectype_WeakList %d\n",Rectype_WeakList);
  Rectype_WeakAnd,
%% printf("#define Rectype_WeakAnd %d\n",Rectype_WeakAnd);
  Rectype_WeakOr,
%% printf("#define Rectype_WeakOr %d\n",Rectype_WeakOr);
  Rectype_WeakAndMapping,
%% printf("#define Rectype_WeakAndMapping %d\n",Rectype_WeakAndMapping);
  Rectype_WeakOrMapping,
%% printf("#define Rectype_WeakOrMapping %d\n",Rectype_WeakOrMapping);
  Rectype_WeakAlist_Key,
%% printf("#define Rectype_WeakAlist_Key %d\n",Rectype_WeakAlist_Key);
  Rectype_WeakAlist_Value,
%% printf("#define Rectype_WeakAlist_Value %d\n",Rectype_WeakAlist_Value);
  Rectype_WeakAlist_Either,
%% printf("#define Rectype_WeakAlist_Either %d\n",Rectype_WeakAlist_Either);
  Rectype_WeakAlist_Both,
%% printf("#define Rectype_WeakAlist_Both %d\n",Rectype_WeakAlist_Both);
  Rectype_WeakHashedAlist_Key,
%% printf("#define Rectype_WeakHashedAlist_Key %d\n",Rectype_WeakHashedAlist_Key);
  Rectype_WeakHashedAlist_Value,
%% printf("#define Rectype_WeakHashedAlist_Value %d\n",Rectype_WeakHashedAlist_Value);
  Rectype_WeakHashedAlist_Either,
%% printf("#define Rectype_WeakHashedAlist_Either %d\n",Rectype_WeakHashedAlist_Either);
  Rectype_WeakHashedAlist_Both,
%% printf("#define Rectype_WeakHashedAlist_Both %d\n",Rectype_WeakHashedAlist_Both);
  rectype_for_broken_compilers_that_dont_like_trailing_commas
};

/* -------------------------- the various types -------------------------- */

/* Cons */
typedef struct {
  gcv_object_t cdr _attribute_aligned_object_; /* CDR */
  gcv_object_t car _attribute_aligned_object_; /* CAR */
} cons_;
typedef cons_ *  Cons;
%% sprintf(buf,"struct { gcv_object_t cdr%s; gcv_object_t car%s; }",attribute_aligned_object,attribute_aligned_object);
%% emit_typedef(buf,"cons_");
%% emit_typedef("cons_ *","Cons");

# Ratio
typedef struct {
  #ifdef SPVW_MIXED
  XRECORD_HEADER
  #endif
  gcv_object_t rt_num _attribute_aligned_object_; # numerator, Integer
  gcv_object_t rt_den _attribute_aligned_object_; # denominator, Integer >0
} ratio_;
typedef ratio_ *  Ratio;
%% #if notused
%% #ifdef SPVW_MIXED
%%   sprintf(buf,"struct { XRECORD_HEADER gcv_object_t rt_num%s; gcv_object_t rt_den%s; }",attribute_aligned_object,attribute_aligned_object);
%% #else
%%   sprintf(buf,"struct { gcv_object_t rt_num%s; gcv_object_t rt_den%s; }",attribute_aligned_object,attribute_aligned_object);
%% #endif
%% emit_typedef(buf,"ratio_");
%% emit_typedef("ratio_ *","Ratio");
%% #endif

# Complex
typedef struct {
  #ifdef SPVW_MIXED
  XRECORD_HEADER
  #endif
  gcv_object_t c_real _attribute_aligned_object_; # real part, real number
  gcv_object_t c_imag _attribute_aligned_object_; # imaginary part, real number
} complex_;
typedef complex_ *  Complex;
%% #if notused
%% #ifdef SPVW_MIXED
%%   sprintf(buf,"struct { XRECORD_HEADER gcv_object_t c_real%s; gcv_object_t c_imag%s; }",attribute_aligned_object,attribute_aligned_object);
%% #else
%%   sprintf(buf,"struct { gcv_object_t c_real%s; gcv_object_t c_imag%s; }",attribute_aligned_object,attribute_aligned_object);
%% #endif
%% emit_typedef(buf,"complex_");
%% emit_typedef("complex_ *","Complex");
%% #endif

# Symbol
typedef struct {
  VAROBJECT_HEADER
  gcv_object_t symvalue    _attribute_aligned_object_; # value cell
  gcv_object_t symfunction _attribute_aligned_object_; # function definition cell
  gcv_object_t hashcode    _attribute_aligned_object_; # hash code
  gcv_object_t proplist    _attribute_aligned_object_; # property list
  gcv_object_t pname       _attribute_aligned_object_; # Printname
  gcv_object_t homepackage _attribute_aligned_object_; # Home-Package or NIL
  # If necessary, add fillers here to ensure sizeof(subr_t) is a multiple of
  # varobject_alignment.
  #if defined(LINUX_NOEXEC_HEAPCODES) && 0
  gcv_object_t filler      _attribute_aligned_object_;
  #endif
} symbol_;
typedef symbol_ *  Symbol;
#ifdef LINUX_NOEXEC_HEAPCODES
  # Compile-time check: sizeof(symbol_) is a multiple of varobject_alignment.
  typedef int symbol_size_check[1 - 2 * (int)(sizeof(symbol_) % varobject_alignment)];
#endif
#define symbol_objects_offset  offsetof(symbol_,symvalue)
#define symbol_length  6
%% #if defined(LINUX_NOEXEC_HEAPCODES) && 0
%%   sprintf(buf,"struct { VAROBJECT_HEADER gcv_object_t symvalue%s; gcv_object_t symfunction%s; gcv_object_t hashcode%s; gcv_object_t proplist%s; gcv_object_t pname%s; gcv_object_t homepackage%s; gcv_object_t filler%s; }",attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object);
%% #else
%%   sprintf(buf,"struct { VAROBJECT_HEADER gcv_object_t symvalue%s; gcv_object_t symfunction%s; gcv_object_t hashcode%s; gcv_object_t proplist%s; gcv_object_t pname%s; gcv_object_t homepackage%s; }",attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object);
%% #endif
%% emit_typedef(buf,"symbol_");
%% emit_typedef("symbol_ *","Symbol");

# Every keyword is a constant.

# Tests whether a symbol is a keyword:
  #define keywordp(sym)  \
    (eq(TheSymbol(sym)->homepackage,O(keyword_package)))

# For constants, the special-bit is meaningless (since constants
# can't be bound lexically nor dynamically).

# Tests whether a symbol is a constant:
  #define constant_var_p(sym)  \
    (((bit(var_bit0_hf)|bit(var_bit1_hf)) & ~((sym)->header_flags)) == 0)

# Tests whether a symbol is a SPECIAL-proclaimed variable or a constant:
  #define special_var_p(sym)  (((sym)->header_flags) & bit(var_bit0_hf))

# Tests whether a symbol is a symbol-macro:
  #define symmacro_var_p(sym)  \
    ((((sym)->header_flags) & bit(var_bit1_hf))            \
     && ((((sym)->header_flags) & bit(var_bit0_hf)) == 0))

# Set the constant-flag of a non-symbol-macro symbol:
  #define set_const_flag(sym)  \
    (((sym)->header_flags) |= bit(var_bit0_hf)|bit(var_bit1_hf))

# Delete the constant-flag of a symbol that is a constant:
# (Symbol must not be a Keyword, comp. spvw.d:case_symbolwithflags)
  #define clear_const_flag(sym)  \
    (((sym)->header_flags) &= ~(bit(var_bit0_hf)|bit(var_bit1_hf)))

# Set the special-flag of a non-symbol-macro symbol:
  #define set_special_flag(sym)  \
    (((sym)->header_flags) |= bit(var_bit0_hf))

# Delete the special-flag of a symbol that is special and non-constant:
  #define clear_special_flag(sym)  \
    (((sym)->header_flags) &= ~bit(var_bit0_hf))

# Set the symbol-macro-flag of a non-special/constant symbol:
  #define set_symmacro_flag(sym)  \
    (((sym)->header_flags) |= bit(var_bit1_hf))

# Delete the symbol-macro-flag of a symbol that is a symbol-macro:
  #define clear_symmacro_flag(sym)  \
    (((sym)->header_flags) &= ~bit(var_bit1_hf))

# Define symbol as constant with given value val.
# val must not trigger the GC!
  #define define_constant(sym,val)                              \
    do { var Symbol sym_from_define_constant = TheSymbol(sym);  \
         set_const_flag(sym_from_define_constant);              \
         sym_from_define_constant->symvalue = (val);            \
    } while(0)

# Define symbol as variable and initialize it with a given value val.
# val must not trigger the GC!
  #define define_variable(sym,val)                              \
    do { var Symbol sym_from_define_variable = TheSymbol(sym);  \
         set_special_flag(sym_from_define_variable);            \
         sym_from_define_variable->symvalue = (val);            \
    } while(0)

# Remove flag-bits of a symbol:
#if defined(NO_symbolflags)
  #define symbol_without_flags(symbol)  symbol
#elif (oint_symbolflags_shift==oint_type_shift)
  #define symbol_without_flags(symbol)  \
    as_object(as_oint(symbol) & (type_zero_oint(symbol_type) | oint_addr_mask))
#else
  #define symbol_without_flags(symbol)  \
    as_object(as_oint(symbol) & ~((wbit(active_bit)|wbit(dynam_bit)|wbit(svar_bit))<<oint_symbolflags_shift))
#endif
/* add a flag to the object */
#define SET_BIT(o,b)  as_object(as_oint(o) | wbit(b));
/* remove a flag from the object */
#define CLR_BIT(o,b)  as_object(as_oint(o) & ~wbit(b));

# Characters

# Integer type holding the data of a character:
#ifdef UNICODE
  #define char_int_len 24  /* anything between 21 and 32 will do */
  #define char_int_limit 0x110000UL
#else
  #define char_int_len 8
  #define char_int_limit 0x100UL
#endif
typedef unsigned_int_with_n_bits(char_int_len)  cint;
#define char_code_limit  char_int_limit
# Converting an integral code to a character:
#define int_char(int_from_int_char)  \
  type_data_object(char_type,(aint)(cint)(int_from_int_char))
# Converting a character to an integral code:
#if !((oint_data_shift==0) && (char_int_len<=oint_data_len) && (exact_uint_size_p(char_int_len)))
  #ifdef TYPECODES
    #define char_int(char_from_char_int)  \
      ((cint)(untype(char_from_char_int)))
  #else
    #if (char_type>>oint_data_shift)==0 || (char_int_len<=16)
      #define char_int(char_from_char_int)  \
        ((cint)(as_oint(char_from_char_int)>>oint_data_shift))
    #else
      #define char_int(char_from_char_int)  \
        ((cint)((as_oint(char_from_char_int)>>oint_data_shift)&(bitm(oint_data_len)-1)))
    #endif
  #endif
#else
  # If oint_data_shift=0, untype does not need to shift. If also
  # char_int_len<=oint_data_len, and if a cint has exactly char_int_len
  # bits, untype does not need to AND.
  #define char_int(char_from_char_int)  \
    ((cint)as_oint(char_from_char_int))
#endif
# Characters can therefore be compared for equality using EQ, this is an
# oint comparison, among the characters a comparison of their integral code.
%% sprintf(buf,"uint%d",char_int_len); emit_typedef(buf,"cint");
%% export_def(int_char(int_from_int_char));
%% export_def(char_int(char_from_char_int));

# A standalone character. Prefer `chart' to `cint' wherever possible because
# it is typesafe. sizeof(chart) = sizeof(cint).
#ifdef CHART_STRUCT
  #ifdef __cplusplus
    struct chart { chart() {} chart(int c) : one_c(c) {} cint one_c; };
  #else
    typedef struct { cint one_c; } chart;
  #endif
#else
  typedef cint chart;
#endif
# Conversions between both:
# as_cint(ch)   chart --> cint
# as_chart(c)   cint --> chart
#ifdef CHART_STRUCT
  #define as_cint(ch)  ((ch).one_c)
  #if 1
    #ifdef __cplusplus
      inline chart as_chart(int c) { return chart(c); }
    #else
      #define as_chart(c)  ((chart){one_c:(c)})
    #endif
  #else
    extern __inline__ chart as_chart (register cint c)
      { register chart ch; ch.one_c = c; return ch; }
  #endif
#else
  #define as_cint(ch)  (ch)
  #define as_chart(c)  (c)
#endif
# Conversion chart --> object.
#define code_char(ch)  int_char(as_cint(ch))
# Conversion object --> chart.
#define char_code(obj)  as_chart(char_int(obj))
# Comparison operations.
#define chareq(ch1,ch2)  (as_cint(ch1) == as_cint(ch2))
#define charlt(ch1,ch2)  (as_cint(ch1) < as_cint(ch2))
#define chargt(ch1,ch2)  (as_cint(ch1) > as_cint(ch2))
%% #ifdef CHART_STRUCT
%%   emit_typedef("struct { cint one_c; }","chart");
%% #else
%%   emit_typedef("cint","chart");
%% #endif
%% export_def(as_cint(ch));
%% export_def(as_chart(c));
%% export_def(code_char(ch));
%% export_def(char_code(obj));

# Conversion standard char (in ASCII encoding) --> chart.
#define ascii(x)  as_chart((uintB)(x))
# Conversion standard char (in ASCII encoding) --> object.
#define ascii_char(x)  code_char(ascii(x))

# Test for STANDARD-CHAR.
#define standard_cint_p(x)  ((('~' >= (x)) && ((x) >= ' ')) || ((x) == NL))

# Whether to use three different kinds of string representations.
#if defined(UNICODE) && (defined(GNU) || (defined(UNIX) && !defined(NO_ALLOCA) && !defined(SPARC)) || defined(BORLAND) || defined(MICROSOFT)) && !defined(NO_SMALL_SSTRING)
#define HAVE_SMALL_SSTRING
#endif

#ifdef HAVE_SMALL_SSTRING
  #define if_HAVE_SMALL_SSTRING(statement)  statement
  # We have three kinds of simple strings, with 8-bit codes (ISO-8859-1
  # strings), with 16-bit codes (UCS-2 strings) and with 32-bit codes
  # (UCS-4/UTF-32 strings).
  typedef uint8 cint8;
  #define cint8_limit (1UL<<8)
  typedef uint16 cint16;
  #define cint16_limit (1UL<<16)
  typedef uint32 cint32;
  #define cint32_limit 0x110000UL
#else
  #define if_HAVE_SMALL_SSTRING(statement)  /*nop*/
  # Only one kind of simple strings.
  typedef cint cint8;
  #define cint8_limit char_int_limit
  typedef cint cint16;
  #define cint16_limit char_int_limit
  typedef cint cint32;
  #define cint32_limit char_int_limit
#endif
%% #ifdef HAVE_SMALL_SSTRING
%%   emit_typedef("uint8","cint8");
%%   emit_typedef("uint16","cint16");
%%   emit_typedef("uint32","cint32");
%% #endif

# Base characters.
#define base_char_int_len char_int_len
#define base_char_code_limit  char_code_limit
# The BASE-CHAR type is defined as
#     (upgraded-array-element-type 'standard-char),
# i.e. the element-type of arrays created with (make-array 'standard-char ...).
# Since it defeats the purpose of UNICODE to have different 8-bit, 16-bit
# and 24-bit character types, we define BASE-CHAR=CHARACTER.

# Fixnums

# fixnum(x) is a fixnum with value x>=0.
# x is an expression with 0 <= x < 2^oint_data_len.
# (Should really be called posfixnum(x).)
#define fixnum(x)  type_data_object(fixnum_type,x)
%% export_def(fixnum(x));

# Fixnum_0 is the number 0, Fixnum_1 is the number 1,
# Fixnum_minus1 is the number -1
#define Fixnum_0  fixnum(0)
#define Fixnum_1  fixnum(1)
#define Fixnum_minus1  type_data_object( fixnum_type | bit(sign_bit_t), vbitm(oint_data_len)-1 )
%% export_def(Fixnum_0);
%% export_def(Fixnum_1);
%% export_def(Fixnum_minus1);

# Value of a non-negative fixnum:
# posfixnum_to_V(obj)
# result is >= 0, < 2^oint_data_len.
#if !(defined(SPARC) && (oint_data_len+oint_data_shift<32))
  #define posfixnum_to_V(obj)  \
    ((uintV)((as_oint(obj)&((oint)wbitm(oint_data_len+oint_data_shift)-1))>>oint_data_shift))
#else
  # Long constants are slower than shifts on a SPARC-processor:
  #define posfixnum_to_V(obj)  \
    ((uintV)((as_oint(obj) << (32-oint_data_len-oint_data_shift)) >> (32-oint_data_len)))
#endif
%% export_def(posfixnum_to_V(obj));

# Value of a negative fixnum:
# negfixnum_to_V(obj)
# Result is >= - 2^oint_data_len, < 0.
#define negfixnum_to_V(obj)  (posfixnum_to_V(obj) | (-vbitm(oint_data_len)))
%% #if notused
%% export_def(negfixnum_to_V(obj));
%% #endif

# Absolute value of a negative fixnum:
# negfixnum_abs_V(obj)
# Result is > 0, <= 2^oint_data_len.
# Beware: Possible wraparound at oint_data_len=intVsize!
#define negfixnum_abs_V(obj)  \
  ((uintV)((as_oint(fixnum_inc(Fixnum_minus1,1))-as_oint(obj))>>oint_data_shift))

# Value of a fixnum, obj should be a variable:
# fixnum_to_V(obj)
# Result is >= - 2^oint_data_len, < 2^oint_data_len and of Type sintV.
# This macro should only be used for oint_data_len+1 <= intLsize!
#if (oint_data_len>=intVsize)
  # No space left for the sign-bit, thus fixnum_to_V = posfixnum_to_V = negfixnum_to_V !
  #define fixnum_to_V(obj)  (sintV)posfixnum_to_V(obj)
#elif (sign_bit_o == oint_data_len+oint_data_shift)
  #define fixnum_to_V(obj)  \
    (((sintV)as_oint(obj) << (intVsize-1-sign_bit_o)) >> (intVsize-1-sign_bit_o+oint_data_shift))
#else
  #if !defined(SPARC)
    #define fixnum_to_V(obj)  \
      (sintV)( ((((sintV)as_oint(obj) >> sign_bit_o) << (intVsize-1)) >> (intVsize-1-oint_data_len)) \
              |((uintV)((as_oint(obj) & ((oint)wbitm(oint_data_len+oint_data_shift)-1)) >> oint_data_shift)) \
             )
  #else
    # Long constants are slower than shifts on a SPARC-processor:
    #define fixnum_to_V(obj)  \
      (sintV)( ((((sintV)as_oint(obj) >> sign_bit_o) << (intVsize-1)) >> (intVsize-1-oint_data_len)) \
              |(((uintV)as_oint(obj) << (intVsize-oint_data_len-oint_data_shift)) >> (intVsize-oint_data_len)) \
             )
  #endif
#endif
%% export_def(fixnum_to_V(obj));

#ifdef intQsize
  # Value of a fixnum, obj should be a variable:
  # fixnum_to_Q(obj)
  # Result is >= - 2^oint_data_len, < 2^oint_data_len.
  #if (sign_bit_o == oint_data_len+oint_data_shift)
    #define fixnum_to_Q(obj)  \
      (((sintQ)as_oint(obj) << (intQsize-1-sign_bit_o)) >> (intQsize-1-sign_bit_o+oint_data_shift))
  #else
    #define fixnum_to_Q(obj)  \
      (sintQ)( ((((sintQ)as_oint(obj) >> sign_bit_o) << (intQsize-1)) >> (intQsize-1-oint_data_len)) \
              |((uintQ)((as_oint(obj) & (wbitm(oint_data_len+oint_data_shift)-1)) >> oint_data_shift)) \
             )
  #endif
#endif

# Add a constant to a non-negative fixnum, given that
# the result is a non-negative fixnum as well:
# fixnum_inc(obj,delta)
# > obj: a fixnum
# > delta: a constant
# < result: incremented fixnum
#define fixnum_inc(obj,delta)                                           \
    objectplus(obj, (soint)(delta) << oint_data_shift)
%% export_def(fixnum_inc(obj,delta));

# posfixnum(x) is a fixnum with value x>=0.
#define posfixnum(x)  fixnum_inc(Fixnum_0,x)
%% export_def(posfixnum(x));

# negfixnum(x) is a fixnum with value x<0.
# (Beware if x is unsigned!)
#define negfixnum(x)  fixnum_inc(fixnum_inc(Fixnum_minus1,1),x)
%% export_def(negfixnum(x));

# sfixnum(x) is a fixnum with value x,
# x is a constant-expression with -2^oint_data_len <= x < 2^oint_data_len.
#define sfixnum(x) ((x)>=0 ? posfixnum(x) : negfixnum(x))
%% export_def(sfixnum(x));

# Convert a character into a fixnum >=0 (the same as for char-int):
#ifdef WIDE_STRUCT
  #define char_to_fixnum(obj)  \
    type_data_object(fixnum_type,untype(obj))
#else
  #define char_to_fixnum(obj)  \
    objectplus(obj,type_zero_oint(fixnum_type)-type_zero_oint(char_type))
#endif

# Make a character from a fitting fixnum >=0 (the same as for int-char):
#ifdef WIDE_STRUCT
  #define fixnum_to_char(obj)  \
    type_data_object(char_type,untype(obj))
#else
  #define fixnum_to_char(obj)  \
    objectplus(obj,type_zero_oint(char_type)-type_zero_oint(fixnum_type))
#endif

# Bignums
typedef struct {
  VAROBJECT_HEADER  # self-pointer for GC
  #ifdef TYPECODES
  uintC length;     # length in digits
  #endif
  uintD data[unspecified]; # number as its two's complement representation
} bignum_;
typedef bignum_ *  Bignum;
# The length is actually an uintWC.
#ifdef TYPECODES
  #define bignum_length(ptr)  ((ptr)->length)
#else
  #define bignum_length(ptr)  srecord_length(ptr)
#endif
#define Bignum_length(obj)  bignum_length(TheBignum(obj))
%% #ifdef TYPECODES
%%   emit_typedef("struct { VAROBJECT_HEADER uintC length; uintD data[unspecified]; }","bignum_");
%% #else
%%   emit_typedef("struct { VAROBJECT_HEADER uintD data[unspecified]; }","bignum_");
%% #endif
%% emit_typedef("bignum_ *","Bignum");
%% export_def(bignum_length(ptr));
%% export_def(Bignum_length(obj));

# Single-Floats
typedef uint32 ffloat; # 32-Bit-Float in IEEE-format
typedef union {
  ffloat eksplicit;    # Value, explicit
  #ifdef FAST_FLOAT
  float machine_float; # Value, as C-'float'
  #endif
} ffloatjanus;
#if !defined(IMMEDIATE_FFLOAT)
typedef struct {
  VAROBJECT_HEADER            # self-pointer for GC
  ffloatjanus representation; # Value
} ffloat_;
typedef ffloat_ *  Ffloat;
#define ffloat_value(obj)  (TheFfloat(obj)->float_value)
#else
# The float-value is stored in the pointer itself, like short-floats.
#define ffloat_value(obj)  ((ffloat)untype(obj))
#endif
%% emit_typedef("uint32","ffloat");
%% emit_typedef("union { ffloat eksplicit; }","ffloatjanus");

# Double-Floats
typedef # 64-Bit-Float in IEEE-format:
        #ifdef intQsize
          # Sign/Exponent/Mantissa
          uint64
        #else
          # Sign/Exponent/MantissaHigh and MantissaLow
          #if BIG_ENDIAN_P || defined(ARM)
            struct {uint32 semhi,mlo;}
          #else
            struct {uint32 mlo,semhi;}
          #endif
        #endif
  dfloat;
typedef union {
  dfloat eksplicit;      # Value, explicit
  #ifdef FAST_DOUBLE
  double machine_double; # Value, as C-'double'
  #endif
} dfloatjanus;
typedef struct {
  VAROBJECT_HEADER            # self-pointer for GC
  dfloatjanus representation; # value
} dfloat_;
typedef dfloat_ *  Dfloat;
%% #ifdef intQsize
%%   emit_typedef("uint64","dfloat");
%% #else
%%   #if BIG_ENDIAN_P
%%     emit_typedef("struct {uint32 semhi,mlo;}","dfloat");
%%   #else
%%     emit_typedef("struct {uint32 mlo,semhi;}","dfloat");
%%   #endif
%% #endif
%% emit_typedef("union { dfloat eksplicit; }","dfloatjanus");

# Single- and Double-Floats
  #define float_value  representation.eksplicit

# Long-Floats
typedef struct {
  VAROBJECT_HEADER   # Self-pointer for GC
  #ifdef TYPECODES
  uintC  len;        # length of the mantissa in digits
  #endif
  uint32 expo;       # exponent
  uintD  data[unspecified]; # mantissa
} lfloat_;
typedef lfloat_ *  Lfloat;
# The length is actually an uintWC.
#ifdef TYPECODES
  #define lfloat_length(ptr)  ((ptr)->len)
#else
  #define lfloat_length(ptr)  srecord_length(ptr)
#endif
#define Lfloat_length(obj)  lfloat_length(TheLfloat(obj))

# simple array (cover simple linear arrays: simple bit vector, simple vector)
typedef struct {
  VRECORD_HEADER # Self-pointer for GC, length in elements
} sarray_;
typedef sarray_ *  Sarray;
#define sarray_length(ptr)  vrecord_length(ptr)
#define Sarray_length(obj)  sarray_length(TheSarray(obj))
%% #if notused
%% emit_typedef("struct { VRECORD_HEADER }","sarray_");
%% emit_typedef("sarray_ *","Sarray");
%% #endif
%% export_def(sarray_length(ptr));
%% export_def(Sarray_length(obj));

# simple bit vector
typedef struct {
  VRECORD_HEADER # self-pointer for GC, length in bits
  uint8  data[unspecified]; # Bits, divided into bytes
} sbvector_;
typedef sbvector_ *  Sbvector;
#define sbvector_length(ptr)  sarray_length(ptr)
#define Sbvector_length(obj)  sbvector_length(TheSbvector(obj))
%% emit_typedef("struct { VRECORD_HEADER uint8  data[unspecified]; }","sbvector_");
%% emit_typedef("sbvector_ *","Sbvector");
%% export_def(sbvector_length(ptr));
%% export_def(Sbvector_length(obj));

# simple string template
#ifdef TYPECODES
  #define SSTRING_HEADER  \
                 VAROBJECT_HEADER # self-pointer for GC \
                 uintL tfl;       # type, flags, length
#else
  #define SSTRING_HEADER  \
                 VAROBJECT_HEADER # self-pointer for GC, tfl
#endif
typedef struct {
  SSTRING_HEADER
} sstring_;
typedef sstring_ *  Sstring;
#define STRUCT_SSTRING(cint_type) \
  struct {                                                                     \
    SSTRING_HEADER /* self-pointer for GC, type+flags, length in characters */ \
    cint_type  data[unspecified];  /* characters */                            \
  }
#ifdef HAVE_SMALL_SSTRING
  typedef STRUCT_SSTRING(cint8)  s8string_;
  typedef s8string_ *  S8string;
  typedef STRUCT_SSTRING(cint16)  s16string_;
  typedef s16string_ *  S16string;
  typedef STRUCT_SSTRING(cint32)  s32string_;
  typedef s32string_ *  S32string;
#else
  # Only one kind of simple strings.
  #ifdef UNICODE
    typedef STRUCT_SSTRING(cint32)  s32string_;
    typedef s32string_ *  S32string;
    # Aliases.
    typedef s32string_  s16string_;
    typedef S32string  S16string;
    typedef s32string_  s8string_;
    typedef S32string  S8string;
  #else
    typedef STRUCT_SSTRING(cint8)  s8string_;
    typedef s8string_ *  S8string;
    # Aliases.
    typedef s8string_  s16string_;
    typedef S8string  S16string;
    typedef s8string_  s32string_;
    typedef S8string  S32string;
  #endif
#endif
# A "normal simple string" is one of maximum-width element type.
# It cannot be reallocated. Only strings with smaller element type
# (called "small simple strings") can be reallocated.
  typedef STRUCT_SSTRING(chart)  snstring_;
  typedef snstring_ *  Snstring;
# These accessors work on any simple string, except reallocated simple-strings.
#ifdef TYPECODES
  #define sstring_length(ptr)  ((ptr)->tfl >> 6)
#else
  #define sstring_length(ptr)  ((ptr)->tfl >> 10)
#endif
#define Sstring_length(obj)  sstring_length(TheSstring(obj))
# Maximum allowed simple-string length:
#ifdef TYPECODES
  #define stringsize_limit_1  ((uintL)(bit(intLsize-6)-1))
#else
  #define stringsize_limit_1  ((uintL)(bit(intLsize-10)-1))
#endif
# Constructing the tfl word:
#ifdef TYPECODES
  #define sstring_tfl(eltype,imm,flags,length)  \
    (((length) << 6) + ((eltype) << 4) + ((imm) << 3) + (flags))
#else
  # This must be consistent with vrecord_tfl.
  #define sstringrecord_tfl(rectype,flags,length)  \
    (((length) << 10) + ((flags) << 8) + (rectype))
  #define sstring_tfl(eltype,imm,flags,length)  \
    sstringrecord_tfl(Rectype_S8string + ((eltype) << 1) + (imm),flags,length)
#endif
# Test whether a simple string is reallocated:
#ifdef HAVE_SMALL_SSTRING
  #ifdef TYPECODES
    #define sstringflags_forwarded_B  bit(2)
    #define sstring_reallocatedp(ptr)  ((ptr)->tfl & sstringflags_forwarded_B)
  #else
    #define sstring_reallocatedp(ptr)  (record_type(ptr) == Rectype_reallocstring)
  #endif
#else
  #define sstring_reallocatedp(ptr)  0
#endif
# Extract the element type of a not-reallocated simple string:
#ifdef TYPECODES
  #define sstring_eltype(ptr)  (((ptr)->tfl >> 4) & 3)
#else
  #define sstring_eltype(ptr)  ((record_type(ptr) - Rectype_S8string) >> 1)
#endif
# Possible values of sstring_eltype:
  #define Sstringtype_8Bit   0
  #define Sstringtype_16Bit  1
  #define Sstringtype_32Bit  2
# Extract the immutable bit of a simple string (reallocated or not):
#ifdef TYPECODES
  #define sstring_immutable(ptr)  (((ptr)->tfl >> 3) & 1)
#else
  #define sstring_immutable(ptr)  ((record_type(ptr) - Rectype_S8string) & 1)
#endif
# Extract the flags of a simple string (reallocated or not):
#ifdef TYPECODES
  # Three bits, containing also sstringflags_forwarded_B.
  #define sstring_flags(ptr)  ((ptr)->tfl & 7)
  #define sstring_flags_clr(ptr,bits)  ((ptr)->tfl &= ~(uintL)(bits))
  #define sstring_flags_set(ptr,bits)  ((ptr)->tfl |= (uintL)(bits))
#else
  #define sstring_flags(ptr)  (((ptr)->tfl >> 8) & 3)
  #define sstring_flags_clr(ptr,bits)  ((ptr)->tfl &= ~((uintL)(bits) << 8))
  #define sstring_flags_set(ptr,bits)  ((ptr)->tfl |= ((uintL)(bits) << 8))
#endif
# Bit masks in the flags. Only used during garbage collection.
  #define sstringflags_backpointer_B  bit(0)
  #define sstringflags_relocated_B    bit(1)
  #define mark_sstring_clean(ptr)  \
    sstring_flags_clr(ptr,sstringflags_backpointer_B|sstringflags_relocated_B)
%% export_def(SSTRING_HEADER);
%% emit_typedef("struct { SSTRING_HEADER }","sstring_");
%% emit_typedef("sstring_ *","Sstring");
%% #ifdef HAVE_SMALL_SSTRING
%%   export_def(STRUCT_SSTRING(cint_type));
%%   emit_typedef("STRUCT_SSTRING(cint8)","s8string_");
%%   emit_typedef("s8string_ *","S8string");
%%   emit_typedef("STRUCT_SSTRING(cint16)","s16string_");
%%   emit_typedef("s16string_ *","S16string");
%%   emit_typedef("STRUCT_SSTRING(cint32)","s32string_");
%%   emit_typedef("s32string_ *","S32string");
%% #endif
%% emit_typedef("struct { SSTRING_HEADER chart data[unspecified]; }","snstring_");
%% emit_typedef("snstring_*","Snstring");
%% export_def(sstring_length(ptr));
%% export_def(Sstring_length(obj));
%% export_def(sstring_eltype(ptr));

# simple vector
typedef struct {
  VRECORD_HEADER # self-pointer for GC, length in objects
  gcv_object_t data[unspecified] _attribute_aligned_object_; # elements
} svector_;
typedef svector_ *  Svector;
#define svector_length(ptr)  sarray_length(ptr)
#define Svector_length(obj)  svector_length(TheSvector(obj))
%% sprintf(buf,"struct { VRECORD_HEADER gcv_object_t data[unspecified]%s; }",attribute_aligned_object);
%% emit_typedef(buf,"svector_");
%% emit_typedef("svector_ *","Svector");

# simple indirect string
typedef struct {
  SSTRING_HEADER   # self-pointer for GC, tfl
  gcv_object_t data _attribute_aligned_object_; # data vector
} sistring_;
typedef sistring_ *  Sistring;
#define sistring_data_offset  offsetof(sistring_,data)

# non-simple indirect Array
typedef struct {
  VAROBJECT_HEADER   # self-pointer for GC
  #ifdef TYPECODES
  uintB flags;       # flags
  uintC rank;        # rank n
  #endif
  gcv_object_t data _attribute_aligned_object_; # data vector
  uintL totalsize;   # totalsize = product of the n dimensions
  uintL dims[unspecified]; # poss. displaced-offset,
                           # n dimensions,
                           # poss. fill-pointer
} iarray_;
typedef iarray_ *  Iarray;
#define iarray_data_offset  offsetof(iarray_,data)
# The rank is actually an uintWC.
# access Rang, Flags:
#ifdef TYPECODES
  #define iarray_rank(ptr)  ((ptr)->rank)
#else
  #define iarray_rank(ptr)  srecord_length(ptr)
#endif
#define Iarray_rank(obj)  iarray_rank(TheIarray(obj))
#ifdef TYPECODES
  #define iarray_flags(ptr)  ((ptr)->flags)
#else
  #define iarray_flags(ptr)  record_flags(ptr)
#endif
#define Iarray_flags(obj)  iarray_flags(TheIarray(obj))
#ifdef TYPECODES
  #define iarray_flags_clr(ptr,bits)  ((ptr)->flags &= ~(bits))
  #define iarray_flags_set(ptr,bits)  ((ptr)->flags |= (bits))
  #define iarray_flags_replace(ptr,newflags)  ((ptr)->flags = (newflags))
#else
  #define iarray_flags_clr(ptr,bits)  record_flags_clr(ptr,bits)
  #define iarray_flags_set(ptr,bits)  record_flags_set(ptr,bits)
  #define iarray_flags_replace(ptr,newflags) record_flags_replace(ptr,newflags)
#endif
# Bits in the Flags:
#define arrayflags_adjustable_bit  7 # set, if array is adjustable
#define arrayflags_fillp_bit       6 # set, if a fill-pointer exists (only possible for n=1)
#define arrayflags_displaced_bit   5 # set, if array is displaced
#define arrayflags_dispoffset_bit  4 # set, if there is space for the
                                     # displaced-offset
                                     # (<==> array adjustable or displaced)
#define arrayflags_atype_mask  0x0F  # mask for the element-type
# Element-types of arrays in Bits 3..0 of its flags:
# The first ones are chosen, so that 2^Atype_nBit = n.
#define Atype_Bit    0  # storage vector is of type sbvector_type
#define Atype_2Bit   1  # storage vector is of type sb2vector_type
#define Atype_4Bit   2  # storage vector is of type sb4vector_type
#define Atype_8Bit   3  # storage vector is of type sb8vector_type
#define Atype_16Bit  4  # storage vector is of type sb16vector_type
#define Atype_32Bit  5  # storage vector is of type sb32vector_type
#define Atype_T      6  # storage vector is of type svector_type
#define Atype_Char   7  # storage vector is of type sstring_type
#define Atype_NIL    8  # storage vector is NIL
%% export_def(Atype_Bit);
%% export_def(Atype_8Bit);
%% export_def(Atype_32Bit);
%% export_def(Atype_T);

# array-types
#ifdef TYPECODES
  #define Array_type(obj)  typecode(obj)
  #define Array_type_bvector     bvector_type      # Iarray
  #define Array_type_b2vector    b2vector_type     # Iarray
  #define Array_type_b4vector    b4vector_type     # Iarray
  #define Array_type_b8vector    b8vector_type     # Iarray
  #define Array_type_b16vector   b16vector_type    # Iarray
  #define Array_type_b32vector   b32vector_type    # Iarray
  #define Array_type_string      string_type       # Iarray
  #define Array_type_vector      vector_type       # Iarray
  #define Array_type_mdarray     mdarray_type      # Iarray
  #define Array_type_sbvector    sbvector_type     # Sbvector
  #define Array_type_sb2vector   sb2vector_type    # Sbvector
  #define Array_type_sb4vector   sb4vector_type    # Sbvector
  #define Array_type_sb8vector   sb8vector_type    # Sbvector
  #define Array_type_sb16vector  sb16vector_type   # Sbvector
  #define Array_type_sb32vector  sb32vector_type   # Sbvector
  #define Array_type_sstring     sstring_type      # Sstring
  #define Array_type_svector     svector_type      # Svector
  #define Array_type_snilvector  symbol_type       # Symbol NIL
  # Array_type_simple_bit_vector(atype)
  # maps Atype_[n]Bit to Array_type_sb[n]vector. Depends on TB0, TB1, TB2.
  # The formula works because there are only 4 possible cases:
  #  (TB0,TB1,TB2)   formula
  #    (0, 1, 2)      atype
  #    (0, 1, 3)      atype + (atype & -4)
  #    (0, 2, 3)      atype + (atype & -2)
  #    (1, 2, 3)      atype + (atype & -1) = atype << 1
  #define Array_type_simple_bit_vector(atype)  \
    (Array_type_sbvector + ((atype)<<TB0) + ((atype)&(bit(TB0+1)-bit(TB1))) + ((atype)&(bit(TB1+1)-bit(TB2))))
#else
  #define Array_type(obj)  Record_type(obj)
  #define Array_type_bvector     Rectype_bvector     # Iarray
  #define Array_type_b2vector    Rectype_b2vector    # Iarray
  #define Array_type_b4vector    Rectype_b4vector    # Iarray
  #define Array_type_b8vector    Rectype_b8vector    # Iarray
  #define Array_type_b16vector   Rectype_b16vector   # Iarray
  #define Array_type_b32vector   Rectype_b32vector   # Iarray
  #define Array_type_string      Rectype_string      # Iarray
  #define Array_type_vector      Rectype_vector      # Iarray
  #define Array_type_mdarray     Rectype_mdarray     # Iarray
  #define Array_type_sbvector    Rectype_Sbvector    # Sbvector
  #define Array_type_sb2vector   Rectype_Sb2vector   # Sbvector
  #define Array_type_sb4vector   Rectype_Sb4vector   # Sbvector
  #define Array_type_sb8vector   Rectype_Sb8vector   # Sbvector
  #define Array_type_sb16vector  Rectype_Sb16vector  # Sbvector
  #define Array_type_sb32vector  Rectype_Sb32vector  # Sbvector
  #define Array_type_sstring     Rectype_S8string: case Rectype_Imm_S8string: case Rectype_S16string: case Rectype_Imm_S16string: case Rectype_S32string: case Rectype_Imm_S32string: case Rectype_reallocstring   # S[8|16|32]string, reallocated string
  #define Array_type_svector     Rectype_Svector     # Svector
  #define Array_type_snilvector  Rectype_Symbol      # Symbol NIL
#endif
# Determining the atype of a [simple-]bit-array:
#define sbNvector_atype(obj)  \
  type_bits_to_atype(Array_type(obj) - Array_type_sbvector)
#define bNvector_atype(obj)  \
  type_bits_to_atype(Array_type(obj) - Array_type_bvector)
#ifdef TYPECODES
  # There are only 4 cases:
  #  (TB0,TB1,TB2)   formula
  #    (0, 1, 2)      type
  #    (0, 1, 3)      (type + (type & 3)) >> 1 = type - ((type & -8) >> 1)
  #    (0, 2, 3)      (type + (type & 1)) >> 1 = type - ((type & -4) >> 1)
  #    (1, 2, 3)      type >> 1                = type - ((type & -2) >> 1)
  #if TB2 > 2
    #define type_bits_to_atype(type)  \
      (((type) + ((type)&(bit(6-TB0-TB1-TB2)-1))) >> 1)
  #else
    #define type_bits_to_atype(type)  (type)
  #endif
#else
  #define type_bits_to_atype(type)  (type)
#endif
%% #ifdef TYPECODES
%%  export_def(Array_type_simple_bit_vector(atype));
%% #endif

/* Packages */
typedef struct {
  XRECORD_HEADER
  gcv_object_t pack_external_symbols  _attribute_aligned_object_;
  gcv_object_t pack_internal_symbols  _attribute_aligned_object_;
  gcv_object_t pack_shadowing_symbols _attribute_aligned_object_;
  gcv_object_t pack_use_list          _attribute_aligned_object_;
  gcv_object_t pack_used_by_list      _attribute_aligned_object_;
  gcv_object_t pack_name              _attribute_aligned_object_;
  gcv_object_t pack_nicknames         _attribute_aligned_object_;
  gcv_object_t pack_docstring         _attribute_aligned_object_;
  gcv_object_t pack_shortest_name     _attribute_aligned_object_;
} *  Package;
#define package_length  ((sizeof(*(Package)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))
/* Some packages are case-sensitive. */
#define mark_pack_casesensitive(obj)  record_flags_set(ThePackage(obj),bit(0))
#define mark_pack_caseinsensitive(obj) record_flags_clr(ThePackage(obj),bit(0))
#define pack_casesensitivep(obj)      (record_flags(ThePackage(obj)) & bit(0))
/* Some packages are case-inverted. */
#define mark_pack_caseinverted(obj)  record_flags_set(ThePackage(obj),bit(1))
#define mark_pack_casepreserved(obj) record_flags_clr(ThePackage(obj),bit(1))
#define pack_caseinvertedp(obj)      (record_flags(ThePackage(obj)) & bit(1))
/* Some packages, such as COMMON-LISP, are locked. */
#define mark_pack_locked(obj)    record_flags_set(ThePackage(obj),bit(2))
#define mark_pack_unlocked(obj)  record_flags_clr(ThePackage(obj),bit(2))
#define pack_locked_p(obj)       (record_flags(ThePackage(obj)) & bit(2))
/* Do not do anything with deleted packages. */
#define mark_pack_deleted(obj)  record_flags_set(ThePackage(obj),bit(7))
#define pack_deletedp(obj)      (record_flags(ThePackage(obj)) & bit(7))
%% #if notused
%% sprintf(buf,"struct { XRECORD_HEADER gcv_object_t pack_external_symbols%s; gcv_object_t pack_internal_symbols%s; gcv_object_t pack_shadowing_symbols%s; gcv_object_t pack_use_list%s; gcv_object_t pack_used_by_list%s; gcv_object_t pack_name%s; gcv_object_t pack_nicknames%s; gcv_object_t pack_docstring%s; gcv_object_t pack_shortest_name%s; } *",attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object);
%% emit_typedef(buf,"Package");
%% #endif

# Hash-Tables
typedef struct {
  XRECORD_HEADER
  #ifdef GENERATIONAL_GC
  gcv_object_t ht_lastrehash         _attribute_aligned_object_;
  #endif
  gcv_object_t ht_maxcount           _attribute_aligned_object_;
  gcv_object_t ht_kvtable            _attribute_aligned_object_;
  gcv_object_t ht_lookupfn           _attribute_aligned_object_;
  gcv_object_t ht_hashcodefn         _attribute_aligned_object_;
  gcv_object_t ht_testfn             _attribute_aligned_object_;
  gcv_object_t ht_gcinvariantfn      _attribute_aligned_object_;
  gcv_object_t ht_rehash_size        _attribute_aligned_object_;
  gcv_object_t ht_mincount_threshold _attribute_aligned_object_;
  gcv_object_t ht_mincount           _attribute_aligned_object_;
  gcv_object_t ht_test               _attribute_aligned_object_; /* hash-table-test - for define-hash-table-test */
  gcv_object_t ht_hash               _attribute_aligned_object_; /* hash function */
  uintL ht_size;
} *  Hashtable;
#ifdef GENERATIONAL_GC
  #define hashtable_length  12
#else
  #define hashtable_length  11
#endif
#define hashtable_xlength  (sizeof(*(Hashtable)0)-offsetofa(record_,recdata)-hashtable_length*sizeof(gcv_object_t))
# Mark a Hash Table as now to reorganize
# set_ht_invalid(TheHashtable(ht));
# mark_ht_invalid(TheHashtable(ht));
# A bit that is set when the list structure is invalid and a rehash is needed.
#define htflags_invalid_B  bit(7)
# A bit that is set if the table has a key whose hash code is not GC-invariant.
#define htflags_gc_rehash_B  bit(6)
#ifdef GENERATIONAL_GC
  #define mark_ht_invalid(ptr)  \
    (record_flags_set(ptr,htflags_invalid_B), \
     (ptr)->ht_lastrehash = unbound)
  #define mark_ht_valid(ptr)  \
    (record_flags_clr(ptr,htflags_invalid_B), \
     (ptr)->ht_lastrehash = O(gc_count))
  #define ht_validp(ptr)  \
    ((record_flags(ptr) & htflags_invalid_B) == 0       \
     && ((record_flags(ptr) & htflags_gc_rehash_B) == 0 \
         || eq((ptr)->ht_lastrehash,O(gc_count))))
#else
  #define mark_ht_invalid(ptr)  record_flags_set(ptr,htflags_invalid_B)
  #define mark_ht_valid(ptr)  record_flags_clr(ptr,htflags_invalid_B)
  #define ht_validp(ptr)  ((record_flags(ptr) & htflags_invalid_B) == 0)
#endif
#ifdef GENERATIONAL_GC
  #define set_ht_invalid(ptr)  mark_ht_invalid(ptr)
  #define set_ht_valid(ptr)  mark_ht_valid(ptr)
#else
  extern bool hash_lookup_builtin (object ht, object obj, bool allowgc, gcv_object_t** KVptr_, gcv_object_t** Iptr_);
  extern bool hash_lookup_builtin_with_rehash (object ht, object obj, bool allowgc, gcv_object_t** KVptr_, gcv_object_t** Iptr_);
  #define set_ht_invalid(ptr)  \
    (mark_ht_invalid(ptr),                                               \
     eq((ptr)->ht_lookupfn,P(hash_lookup_builtin))                       \
     ? ((ptr)->ht_lookupfn = P(hash_lookup_builtin_with_rehash), 0) : 0)
  #define set_ht_valid(ptr)  \
    (mark_ht_valid(ptr),                                       \
     eq((ptr)->ht_lookupfn,P(hash_lookup_builtin_with_rehash)) \
     ? ((ptr)->ht_lookupfn = P(hash_lookup_builtin), 0) : 0)
#endif
#define set_ht_invalid_if_needed(ptr)  \
  if (record_flags(ptr) & htflags_gc_rehash_B) \
    set_ht_invalid(ptr)/*;*/
# A bit that indicates whether to warn about this situation.
#define htflags_warn_gc_rehash_B  bit(5)
# Extract the part of the flags that encodes the test.
#define ht_test_code(flags)  \
  (flags & (bit(0) | bit(1) | bit(2) | bit(3)))
# Tests whether a test code indicates a user-defined test function.
#define ht_test_code_user_p(test_code)  \
  (((test_code) & bit(2)) != 0)
# Test whether a hash table is weak.
#define ht_weak_p(ht)  \
  !simple_vector_p(TheHashtable(ht)->ht_kvtable)
# The kvtable array is either a HashedAlist or a WeakHashedAlist.
# Both share the same layout, i.e.
#   &((HashedAlist)0)->hal_data == &((WeakHashedAlist)0)->whal_data.
typedef struct {
  VRECORD_HEADER # self-pointer for GC, length in objects
  gcv_object_t hal_filler            _attribute_aligned_object_; # for consistency with WeakHashedAlist
  gcv_object_t hal_itable            _attribute_aligned_object_; # index-vector
  gcv_object_t hal_count             _attribute_aligned_object_; # remaining pairs
  gcv_object_t hal_freelist          _attribute_aligned_object_; # start index of freelist
  gcv_object_t hal_data[unspecified] _attribute_aligned_object_; # (key, value, next) triples
} * HashedAlist;
# TheHashedAlist is used to access both HashedAlist and WeakHashedAlist.
#define TheHashedAlist(obj)  ((HashedAlist)TheVarobject(obj))

# Readtables
typedef struct {
  XRECORD_HEADER
  gcv_object_t readtable_syntax_table _attribute_aligned_object_;
  gcv_object_t readtable_macro_table  _attribute_aligned_object_;
  gcv_object_t readtable_case         _attribute_aligned_object_;
} *  Readtable;
#define readtable_length  ((sizeof(*(Readtable)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# Pathnames
typedef struct {
  XRECORD_HEADER
  #if HAS_HOST
    gcv_object_t pathname_host      _attribute_aligned_object_;
  #endif
  #if HAS_DEVICE
    gcv_object_t pathname_device    _attribute_aligned_object_;
  #endif
  #if 1
    gcv_object_t pathname_directory _attribute_aligned_object_;
    gcv_object_t pathname_name      _attribute_aligned_object_;
    gcv_object_t pathname_type      _attribute_aligned_object_;
    gcv_object_t pathname_version   _attribute_aligned_object_;
  #endif
} *  Pathname;
#define pathname_length  ((sizeof(*(Pathname)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

#ifdef LOGICAL_PATHNAMES
# Logical Pathnames
typedef struct {
  XRECORD_HEADER
  gcv_object_t pathname_host      _attribute_aligned_object_;
  gcv_object_t pathname_directory _attribute_aligned_object_;
  gcv_object_t pathname_name      _attribute_aligned_object_;
  gcv_object_t pathname_type      _attribute_aligned_object_;
  gcv_object_t pathname_version   _attribute_aligned_object_;
} *  Logpathname;
#define logpathname_length  ((sizeof(*(Logpathname)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))
#endif

# Random-States
typedef struct {
  XRECORD_HEADER
  gcv_object_t random_state_seed _attribute_aligned_object_;
} *  Random_state;
#define random_state_length  ((sizeof(*(Random_state)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# Bytes
typedef struct {
  XRECORD_HEADER
  gcv_object_t byte_size     _attribute_aligned_object_;
  gcv_object_t byte_position _attribute_aligned_object_;
} *  Byte;
#define byte_length  ((sizeof(*(Byte)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# Fsubrs
typedef struct {
  XRECORD_HEADER
  gcv_object_t name    _attribute_aligned_object_;
  gcv_object_t argtype _attribute_aligned_object_;
  void* function; # actually a fsubr_function_t*
} *  Fsubr;
#define fsubr_length  2
#define fsubr_xlength  (sizeof(*(Fsubr)0)-offsetofa(record_,recdata)-fsubr_length*sizeof(gcv_object_t))

# Load-time-evals
typedef struct {
  XRECORD_HEADER
  gcv_object_t loadtimeeval_form _attribute_aligned_object_;
} *  Loadtimeeval;
#define loadtimeeval_length  ((sizeof(*(Loadtimeeval)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# Symbol-macros
typedef struct {
  XRECORD_HEADER
  gcv_object_t symbolmacro_expansion _attribute_aligned_object_;
} *  Symbolmacro;
#define symbolmacro_length  ((sizeof(*(Symbolmacro)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# Global-Symbol-macros
typedef struct {
  XRECORD_HEADER
  gcv_object_t globalsymbolmacro_definition _attribute_aligned_object_;
} *  GlobalSymbolmacro;
#define globalsymbolmacro_length  ((sizeof(*(GlobalSymbolmacro)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# Macros
typedef struct {
  XRECORD_HEADER
  gcv_object_t macro_expander _attribute_aligned_object_;
} *  Macro;
#define macro_length  ((sizeof(*(Macro)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# FunctionMacros
typedef struct {
  XRECORD_HEADER
  gcv_object_t functionmacro_macro_expander _attribute_aligned_object_;
  gcv_object_t functionmacro_function       _attribute_aligned_object_;
} *  FunctionMacro;
#define functionmacro_length  ((sizeof(*(FunctionMacro)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# BigReadLabel
typedef struct {
  XRECORD_HEADER
  gcv_object_t brl_value _attribute_aligned_object_;
} *  BigReadLabel;
#define bigreadlabel_length  ((sizeof(*(BigReadLabel)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# Encoding
typedef struct {
  XRECORD_HEADER
  gcv_object_t enc_eol         _attribute_aligned_object_; # line termination, a keyword (:UNIX, :MAC, :DOS)
  gcv_object_t enc_towcs_error _attribute_aligned_object_; # input error action, :ERROR or :IGNORE or a character
  gcv_object_t enc_tombs_error _attribute_aligned_object_; # output error action, :ERROR or :IGNORE or a character or an uint8
  #ifdef UNICODE
  gcv_object_t enc_charset     _attribute_aligned_object_; # character set, a symbol in the CHARSET package
                                                           # or a simple-string
  # Functions to convert bytes to characters.
    gcv_object_t enc_mblen     _attribute_aligned_object_; # uintL (*) (object encoding, const uintB* src, const uintB* srcend);
    gcv_object_t enc_mbstowcs  _attribute_aligned_object_; # void (*) (object encoding, object stream, const uintB* *srcp, const uintB* srcend, chart* *destp, chart* destend);
  # Functions to convert characters to bytes.
    gcv_object_t enc_wcslen    _attribute_aligned_object_; # uintL (*) (object encoding, const chart* src, const chart* srcend);
    gcv_object_t enc_wcstombs  _attribute_aligned_object_; # void (*) (object encoding, object stream, const chart* *srcp, const chart* srcend, uintB* *destp, uintB* destend);
  # Function to return the set of defined characters in the range [start,end],
  # as a simple-string of intervals #(start1 end1 ... startm endm).
    gcv_object_t enc_range     _attribute_aligned_object_; # object (*) (object encoding, uintL start, uintL end, uintL maxintervals);
  # An auxiliary pointer.
  gcv_object_t enc_table       _attribute_aligned_object_;
  # Minimum number of bytes needed to represent a character
  # caveat: correct only for some encodings, defaults to 1
  uintL min_bytes_per_char;
  # Maximum number of bytes needed to represent a character
  # caveat: correct only for some encodings, defaults to 8
  uintL max_bytes_per_char;
  #endif
} *  Encoding;
#ifdef UNICODE
  #define encoding_length  10
#else
  #define encoding_length  3
#endif
#define encoding_xlength  (sizeof(*(Encoding)0)-offsetofa(record_,recdata)-encoding_length*sizeof(gcv_object_t))
#ifdef UNICODE
  #define Encoding_mblen(encoding)  ((uintL (*) (object, const uintB*, const uintB*)) ThePseudofun(TheEncoding(encoding)->enc_mblen))
  #define Encoding_mbstowcs(encoding)  ((void (*) (object, object, const uintB**, const uintB*, chart**, chart*)) ThePseudofun(TheEncoding(encoding)->enc_mbstowcs))
  #define Encoding_wcslen(encoding)  ((uintL (*) (object, const chart*, const chart*)) ThePseudofun(TheEncoding(encoding)->enc_wcslen))
  #define Encoding_wcstombs(encoding)  ((void (*) (object, object, const chart**, const chart*, uintB**, uintB*)) ThePseudofun(TheEncoding(encoding)->enc_wcstombs))
  #define Encoding_range(encoding)  ((object (*) (object, uintL, uintL, uintL)) ThePseudofun(TheEncoding(encoding)->enc_range))
#endif
#ifdef UNICODE
  #define cslen(encoding,src,srclen)  \
    Encoding_wcslen(encoding)(encoding,src,(src)+(srclen))
  #define cstombs_help_(encoding,src,srclen,dest,destlen,A)         \
    do { var const chart* _srcptr = (src);                          \
      var const chart* _srcendptr = _srcptr+(srclen);               \
      var uintB* _destptr = (dest);                                 \
      var uintB* _destendptr = _destptr+(destlen);                  \
      Encoding_wcstombs(encoding)(encoding,nullobj,&_srcptr,_srcendptr,&_destptr,_destendptr); \
      A((_srcptr == _srcendptr) && (_destptr == _destendptr)); \
    } while(0)
#else
  #define cslen(encoding,src,srclen)  (srclen)
  #define cstombs_help_(encoding,src,srclen,dest,destlen,A)           \
    do { A((srclen) == (destlen));                                   \
         begin_system_call(); memcpy(dest,src,srclen); end_system_call(); \
    } while(0)
#endif
#define cstombs(encoding,src,srclen,dest,destlen)  cstombs_help_(encoding,src,srclen,dest,destlen,ASSERT)
%% sprintf(buf,"struct { XRECORD_HEADER gcv_object_t enc_eol%s; gcv_object_t enc_towcs_error%s; gcv_object_t enc_tombs_error%s;",attribute_aligned_object,attribute_aligned_object,attribute_aligned_object);
%% #ifdef UNICODE
%%   sprintf(buf+strlen(buf)," gcv_object_t enc_charset%s; gcv_object_t enc_mblen%s; gcv_object_t enc_mbstowcs%s; gcv_object_t enc_wcslen%s; gcv_object_t enc_wcstombs%s; gcv_object_t enc_range%s; gcv_object_t enc_table%s; uintL min_bytes_per_char; uintL max_bytes_per_char;",attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object,attribute_aligned_object);
%% #endif
%% strcat(buf," } *");
%% emit_typedef(buf,"Encoding");
%% #ifdef UNICODE
%%  export_def(Encoding_wcslen(encoding));
%%  export_def(Encoding_wcstombs(encoding));
%% #endif
%% export_def(cslen(encoding,src,srclen));
%% export_def(cstombs_help_(encoding,src,srclen,dest,destlen,A));
%% puts("#define cstombs(encoding,src,srclen,dest,destlen)  cstombs_help_(encoding,src,srclen,dest,destlen,ASSERT)");

#ifdef FOREIGN
# foreign pointer wrap
typedef struct {
  XRECORD_HEADER
  void* fp_pointer;
} *  Fpointer;
#define fpointer_length  0
#define fpointer_xlength  (sizeof(*(Fpointer)0)-offsetofa(record_,recdata)-fpointer_length*sizeof(gcv_object_t))
#define mark_fp_invalid(ptr)  record_flags_set(ptr,bit(7))
#define mark_fp_valid(ptr)  record_flags_clr(ptr,bit(7))
#define fp_validp(ptr)  ((record_flags(ptr) & bit(7)) == 0)
#else
#define mark_fp_invalid(ptr)
#endif
%% #ifdef FOREIGN
%%   emit_typedef("struct { XRECORD_HEADER void* fp_pointer;} *","Fpointer");
%%   export_def(fp_validp(ptr));
%%   export_def(mark_fp_invalid(ptr));
%% #endif

#ifdef DYNAMIC_FFI

# foreign adresses
typedef struct {
  XRECORD_HEADER
  gcv_object_t fa_base _attribute_aligned_object_;
  sintP fa_offset;
} * Faddress;
#define faddress_length  1
#define faddress_xlength  (sizeof(*(Faddress)0)-offsetofa(record_,recdata)-faddress_length*sizeof(gcv_object_t))

# foreign variables
typedef struct {
  XRECORD_HEADER
  gcv_object_t fv_name    _attribute_aligned_object_;
  gcv_object_t fv_address _attribute_aligned_object_;
  gcv_object_t fv_size    _attribute_aligned_object_;
  gcv_object_t fv_type    _attribute_aligned_object_;
} * Fvariable;
#define fvariable_length  ((sizeof(*(Fvariable)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# foreign functions
typedef struct {
  XRECORD_HEADER
  gcv_object_t ff_name       _attribute_aligned_object_;
  gcv_object_t ff_address    _attribute_aligned_object_;
  gcv_object_t ff_resulttype _attribute_aligned_object_;
  gcv_object_t ff_argtypes   _attribute_aligned_object_;
  gcv_object_t ff_flags      _attribute_aligned_object_;
  gcv_object_t ff_properties _attribute_aligned_object_;
} * Ffunction;
#define ffunction_length  ((sizeof(*(Ffunction)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

#endif

# weak pointer
typedef struct {
  XRECORD_HEADER
  gcv_object_t wp_cdr   _attribute_aligned_object_; # active weak-pointers form a chained list
  gcv_object_t wp_value _attribute_aligned_object_; # the referenced object
} * Weakpointer;
/* Both wp_cdr and wp_value are invisible to gc_mark routines.
 When the weak-pointer becomes inactive, both fields are turned to unbound.
 When wp_value is GC-invariant, WP does not have to be on the
 O(all_weakpointers) list!  WP is on the list <==> ( wp_cdr != unbound ) */
#define weakpointer_length  ((sizeof(*(Weakpointer)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))
#define weakpointer_broken_p(wp) (!boundp(TheWeakpointer(wp)->wp_value))

# weak list
typedef struct {
  LRECORD_HEADER
  gcv_object_t wp_cdr                   _attribute_aligned_object_; # active weak-pointers form a chained list
  gcv_object_t wl_count                 _attribute_aligned_object_; # remaining objects
  gcv_object_t wl_elements[unspecified] _attribute_aligned_object_; # the referenced objects
} * WeakList;

# mutable weak list
typedef struct {
  XRECORD_HEADER
  gcv_object_t mwl_list _attribute_aligned_object_;
} * MutableWeakList;
#define mutableweaklist_length  ((sizeof(*(MutableWeakList)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# weak "and" relation
typedef struct {
  LRECORD_HEADER
  gcv_object_t wp_cdr                _attribute_aligned_object_; # active weak-pointers form a chained list
  gcv_object_t war_keys_list         _attribute_aligned_object_; # list to copy the keys into
  gcv_object_t war_keys[unspecified] _attribute_aligned_object_; # the referenced objects
} * WeakAnd;

# weak "or" relation
typedef struct {
  LRECORD_HEADER
  gcv_object_t wp_cdr                _attribute_aligned_object_; # active weak-pointers form a chained list
  gcv_object_t wor_keys_list         _attribute_aligned_object_; # list to copy the keys into
  gcv_object_t wor_keys[unspecified] _attribute_aligned_object_; # the referenced objects
} * WeakOr;

# weak mapping
typedef struct {
  XRECORD_HEADER
  gcv_object_t wp_cdr   _attribute_aligned_object_; # active weak-pointers form a chained list
  gcv_object_t wm_value _attribute_aligned_object_; # the dependent referenced object
  gcv_object_t wm_key   _attribute_aligned_object_; # the weak referenced object
} * Weakmapping;
#define weakmapping_length  ((sizeof(*(Weakmapping)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# weak "and" mapping
typedef struct {
  LRECORD_HEADER
  gcv_object_t wp_cdr                _attribute_aligned_object_; # active weak-pointers form a chained list
  gcv_object_t wam_value             _attribute_aligned_object_; # the dependent referenced object
  gcv_object_t wam_keys_list         _attribute_aligned_object_; # list to copy the keys into
  gcv_object_t wam_keys[unspecified] _attribute_aligned_object_; # the referenced objects
} * WeakAndMapping;

# weak "or" mapping
typedef struct {
  LRECORD_HEADER
  gcv_object_t wp_cdr                _attribute_aligned_object_; # active weak-pointers form a chained list
  gcv_object_t wom_value             _attribute_aligned_object_; # the dependent referenced object
  gcv_object_t wom_keys_list         _attribute_aligned_object_; # list to copy the keys into
  gcv_object_t wom_keys[unspecified] _attribute_aligned_object_; # the referenced objects
} * WeakOrMapping;

# weak alist (rectype = Rectype_WeakAlist_{Key,Value,Either,Both})
typedef struct {
  LRECORD_HEADER
  gcv_object_t wp_cdr                _attribute_aligned_object_; # active weak-pointers form a chained list
  gcv_object_t wal_count             _attribute_aligned_object_; # remaining pairs
  gcv_object_t wal_data[unspecified] _attribute_aligned_object_; # key, value alternating
} * WeakAlist;

# mutable weak alist
typedef struct {
  XRECORD_HEADER
  gcv_object_t mwal_list _attribute_aligned_object_;
} * MutableWeakAlist;
#define mutableweakalist_length  ((sizeof(*(MutableWeakAlist)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# weak hashed alist (rectype = Rectype_WeakHashedAlist_{Key,Value,Either,Both})
typedef struct {
  LRECORD_HEADER
  gcv_object_t wp_cdr                 _attribute_aligned_object_; # active weak-pointers form a chained list
  gcv_object_t whal_itable            _attribute_aligned_object_; # index-vector
  gcv_object_t whal_count             _attribute_aligned_object_; # remaining pairs
  gcv_object_t whal_freelist          _attribute_aligned_object_; # start index of freelist
  gcv_object_t whal_data[unspecified] _attribute_aligned_object_; # (key, value, next) triples
} * WeakHashedAlist;

# Finalizer
typedef struct {
  XRECORD_HEADER
  gcv_object_t fin_alive    _attribute_aligned_object_; # only if this object is alive
  gcv_object_t fin_trigger  _attribute_aligned_object_; # wait for the death of this object
  gcv_object_t fin_function _attribute_aligned_object_; # then this function is called
  gcv_object_t fin_cdr      _attribute_aligned_object_;
} * Finalizer;
#define finalizer_length  ((sizeof(*(Finalizer)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

#ifdef SOCKET_STREAMS
# Socket-Server
typedef struct {
  XRECORD_HEADER
  gcv_object_t socket_handle _attribute_aligned_object_; # socket handle
  gcv_object_t host          _attribute_aligned_object_; # host string
  gcv_object_t port          _attribute_aligned_object_; # port number
} * Socket_server;
#define socket_server_length  ((sizeof(*(Socket_server)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

# Information about any of the two ends of a socket connection.
#ifndef MAXHOSTNAMELEN
  #define MAXHOSTNAMELEN 64
#endif
typedef struct host_data_t {
  char hostname[45+1];   # IP address in presentable, printable format
                         # (IPv4 max. 15 characters, IPv6 max. 45 characters)
  char truename[MAXHOSTNAMELEN+1]; # hostname, with or without domain name
  unsigned int port;
} host_data_t;
#endif

#ifdef YET_ANOTHER_RECORD

# Yet another record
typedef struct {
  XRECORD_HEADER
  gcv_object_t yetanother_x _attribute_aligned_object_;
  gcv_object_t yetanother_y _attribute_aligned_object_;
  gcv_object_t yetanother_z _attribute_aligned_object_;
} * Yetanother;
#define yetanother_length  ((sizeof(*(Yetanother)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))

#endif

# Streams with metaclass BUILT-IN-CLASS
typedef struct {
  #ifdef case_stream
    VAROBJECT_HEADER # self-pointer for GC
    uintB strmtype;  # subtype (as sintB >=0 !)
    uintB strmflags; # flags
    uintB reclength; # length in object
    uintB recxlength; # lengths of the extra-elements
  #else
    # Because of space requirements, I have to put strmflags and strmtype
    # into a fixnum in recdata[0].
    #if !((oint_addr_len+oint_addr_shift>=24) && (8>=oint_addr_shift))
      #error "No room for stream flags -- re-accommodate Stream-Flags!!"
    #endif
    XRECORD_HEADER
    #if defined(WIDE) && BIG_ENDIAN_P
      uintL strmfiller0;
    #endif
    uintB strmfiller1;
    uintB strmflags; # Flags
    uintB strmtype;  # Subtype
    uintB strmfiller2;
    #if defined(WIDE) && !BIG_ENDIAN_P
      uintL strmfiller0;
    #endif
  #endif
  gcv_object_t strm_rd_by            _attribute_aligned_object_;
  gcv_object_t strm_rd_by_array      _attribute_aligned_object_;
  gcv_object_t strm_wr_by            _attribute_aligned_object_;
  gcv_object_t strm_wr_by_array      _attribute_aligned_object_;
  gcv_object_t strm_rd_ch            _attribute_aligned_object_;
  gcv_object_t strm_pk_ch            _attribute_aligned_object_;
  gcv_object_t strm_rd_ch_array      _attribute_aligned_object_;
  gcv_object_t strm_rd_ch_last       _attribute_aligned_object_;
  gcv_object_t strm_wr_ch            _attribute_aligned_object_;
  gcv_object_t strm_wr_ch_array      _attribute_aligned_object_;
  gcv_object_t strm_wr_ch_npnl       _attribute_aligned_object_;
  gcv_object_t strm_wr_ch_array_npnl _attribute_aligned_object_;
  gcv_object_t strm_wr_ch_lpos       _attribute_aligned_object_;
  gcv_object_t strm_other[unspecified] _attribute_aligned_object_; # type-specific components
} *  Stream;
# The macro TheStream actually means TheBuiltinStream.
#define strm_len  ((sizeof(*(Stream)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t)-unspecified)
#define stream_length(ptr)  xrecord_length(ptr)
#define stream_xlength(ptr)  xrecord_xlength(ptr)
#define Stream_length(obj)  stream_length(TheStream(obj))
#define Stream_xlength(obj)  stream_xlength(TheStream(obj))
# Bit-masks in the Flags:
  #define strmflags_open_bit_B   0  # set, if the Stream is open
  #define strmflags_immut_bit_B  1  # set if read literals are immutable
  #define strmflags_reval_bit_B  2  # set, if Read-Eval is permitted
  #define strmflags_rd_by_bit_B  4  # set, if READ-BYTE is possible
  #define strmflags_wr_by_bit_B  5  # set, if WRITE-BYTE is possible
  #define strmflags_rd_ch_bit_B  6  # set, if READ-CHAR is possible
  #define strmflags_wr_ch_bit_B  7  # set, if WRITE-CHAR is possible
  #define strmflags_open_B   bit(strmflags_open_bit_B)
  #define strmflags_rd_by_B  bit(strmflags_rd_by_bit_B)
  #define strmflags_wr_by_B  bit(strmflags_wr_by_bit_B)
  #define strmflags_rd_ch_B  bit(strmflags_rd_ch_bit_B)
  #define strmflags_wr_ch_B  bit(strmflags_wr_ch_bit_B)
  #define strmflags_rd_B  (strmflags_rd_by_B | strmflags_rd_ch_B)
  #define strmflags_wr_B  (strmflags_wr_by_B | strmflags_wr_ch_B)
# approach Typinfo:
  enum { # The ordered values of this enumeration are 0,1,2,...
  # First the OS independent streams.
                              enum_strmtype_synonym,
  #define strmtype_synonym    (uintB)enum_strmtype_synonym
                              enum_strmtype_broad,
  #define strmtype_broad      (uintB)enum_strmtype_broad
                              enum_strmtype_concat,
  #define strmtype_concat     (uintB)enum_strmtype_concat
                              enum_strmtype_twoway,
  #define strmtype_twoway     (uintB)enum_strmtype_twoway
                              enum_strmtype_echo,
  #define strmtype_echo       (uintB)enum_strmtype_echo
                              enum_strmtype_str_in,
  #define strmtype_str_in     (uintB)enum_strmtype_str_in
                              enum_strmtype_str_out,
  #define strmtype_str_out    (uintB)enum_strmtype_str_out
                              enum_strmtype_str_push,
  #define strmtype_str_push   (uintB)enum_strmtype_str_push
                              enum_strmtype_pphelp,
  #define strmtype_pphelp     (uintB)enum_strmtype_pphelp
                              enum_strmtype_buff_in,
  #define strmtype_buff_in    (uintB)enum_strmtype_buff_in
                              enum_strmtype_buff_out,
  #define strmtype_buff_out   (uintB)enum_strmtype_buff_out
  #ifdef GENERIC_STREAMS
                              enum_strmtype_generic,
  #define strmtype_generic    (uintB)enum_strmtype_generic
  #endif
  # Then the OS dependent streams.
                              enum_strmtype_file,
  #define strmtype_file       (uintB)enum_strmtype_file
  #ifdef KEYBOARD
                              enum_strmtype_keyboard,
  #define strmtype_keyboard   (uintB)enum_strmtype_keyboard
  #endif
                              enum_strmtype_terminal,
  #define strmtype_terminal   (uintB)enum_strmtype_terminal
  #ifdef SCREEN
                              enum_strmtype_window,
  #define strmtype_window     (uintB)enum_strmtype_window
  #endif
  #ifdef PRINTER
                              enum_strmtype_printer,
  #define strmtype_printer    (uintB)enum_strmtype_printer
  #endif
  #ifdef PIPES
                              enum_strmtype_pipe_in,
  #define strmtype_pipe_in    (uintB)enum_strmtype_pipe_in
                              enum_strmtype_pipe_out,
  #define strmtype_pipe_out   (uintB)enum_strmtype_pipe_out
  #endif
  #ifdef X11SOCKETS
                              enum_strmtype_x11socket,
  #define strmtype_x11socket  (uintB)enum_strmtype_x11socket
  #endif
  #ifdef SOCKET_STREAMS
                              enum_strmtype_socket,
  #define strmtype_socket     (uintB)enum_strmtype_socket
                              enum_strmtype_twoway_socket,
  #define strmtype_twoway_socket (uintB)enum_strmtype_twoway_socket
  #endif
                              enum_strmtype_dummy
  };
  # When this table is changed, also adapt
  # - the 12 jumptables for STREAM-ELEMENT-TYPE, SET-STREAM-ELEMENT-TYPE,
  #   STREAM-EXTERNAL-FORMAT, SET-STREAM-EXTERNAL-FORMAT, INTERACTIVE-STREAM-P,
  #   CLOSE, LISTEN-CHAR, CLEAR_INPUT, LISTEN-BYTE, FINISH_OUTPUT,
  #   FORCE_OUTPUT, CLEAR_OUTPUT in STREAM.D and
  # - the name-table in CONSTOBJ.D and
  # - the jumptable for PR_STREAM in IO.D and
  # - the pseudo-function-table in PSEUDOFUN.D
  #
# more type-specific components:
  #define strm_eltype          strm_other[0] # CHARACTER or ([UN]SIGNED-BYTE n)
  #define strm_encoding        strm_other[1] # an encoding
  #define strm_file_name       strm_other[6] # filename, a pathname or NIL
  #define strm_file_truename   strm_other[7] # truename, a non-logical pathname or NIL
  #define strm_buffered_channel  strm_other[5] # packed Handle
  #define strm_synonym_symbol  strm_other[0]
  #define strm_broad_list      strm_other[0] # list of Streams
  #define strm_concat_list     strm_other[0] # list of Streams
  #define strm_pphelp_lpos     strm_wr_ch_lpos # Line Position (Fixnum>=0)
  #define strm_pphelp_strings  strm_other[0]   # Semi-Simple-Strings for Output
  #define strm_pphelp_modus    strm_other[1]   # Mode (NIL=Single line, T=multiple lines)
  #define strm_pphelp_miserp   strm_other[2] # miser mode indicator
  #define strm_pphelp_offset   strm_other[3] # initial line offset (indent)
  #define strm_buff_in_fun     strm_other[0] # read function
  #define strm_buff_out_fun    strm_other[0] # output function
  #define strm_twoway_input    strm_other[0] # stream for input
  #define strm_twoway_output   strm_other[1] # stream for output
  #ifdef PIPES
  #define strm_pipe_pid        strm_other[6] # process-Id, a Fixnum >=0
  #endif
  #ifdef X11SOCKETS
  #define strm_x11socket_connect  strm_other[6] # List (host display)
  #endif
  #ifdef SOCKET_STREAMS
  #define strm_socket_port     strm_other[6] # port, a fixnum >=0
  #define strm_socket_host     strm_other[7] # host, NIL or a string
  #define strm_twoway_socket_input  strm_other[0] # input side, a socket stream
  #endif
  #ifdef GENERIC_STREAMS
  #define strm_controller_object strm_other[0] # Controller (usually a CLOS-instance)
  #endif
  #define strm_buffered_bufflen 4096   # buffer length, a power of 2, <2^16
# is used by stream.d, pathname.d, io.d
%% export_def(strm_buffered_bufflen);

# Structures
typedef Srecord  Structure;
  #define structure_types   recdata[0]
#define structure_length(ptr)  srecord_length(ptr)
#define Structure_length(obj)  structure_length(TheStructure(obj))
%% emit_typedef("Srecord","Structure");
%% export_def(structure_types);

# CLOS class-versions, see clos.lisp
typedef struct {
  VRECORD_HEADER
  gcv_object_t cv_newest_class             _attribute_aligned_object_; # the CLASS object describing the newest available version
  gcv_object_t cv_class                    _attribute_aligned_object_; # the CLASS object describing the slots
  gcv_object_t cv_shared_slots             _attribute_aligned_object_; # simple-vector with the values of all shared slots, or nil
  gcv_object_t cv_serial                   _attribute_aligned_object_; # serial number of this class version
  gcv_object_t cv_next                     _attribute_aligned_object_; # next class-version, or nil
  gcv_object_t cv_slotlists_valid_p        _attribute_aligned_object_; # true if the following fields are already computed
  gcv_object_t cv_kept_slot_locations      _attribute_aligned_object_; # plist of old and new slot locations of those slots that remain local or were shared and become local
  gcv_object_t cv_added_slots              _attribute_aligned_object_; # list of local slots that are added in the next version
  gcv_object_t cv_discarded_slots          _attribute_aligned_object_; # list of local slots that are removed or become shared in the next version
  gcv_object_t cv_discarded_slot_locations _attribute_aligned_object_; # plist of local slots and their old slot locations that are removed or become shared in the next version
} *  ClassVersion;
#define classversion_length  ((sizeof(*(ClassVersion)0)-offsetofa(svector_,data))/sizeof(gcv_object_t))

# CLOS-instances
typedef struct {
  SRECORD_HEADER
  gcv_object_t inst_class_version _attribute_aligned_object_; # indirect pointer to a CLOS-class
  gcv_object_t other[unspecified] _attribute_aligned_object_;
} *  Instance;
# Bit masks in the flags:
  #define instflags_forwarded_B    bit(0)
  #define instflags_beingupdated_B bit(3)
  # The following are only used during garbage collection.
  #define instflags_backpointer_B  bit(1)
  #define instflags_relocated_B    bit(2)
  #define mark_inst_clean(ptr)  \
    record_flags_clr(ptr,instflags_backpointer_B|instflags_relocated_B)
%% sprintf(buf,"struct { SRECORD_HEADER gcv_object_t inst_class_version%s; gcv_object_t other[unspecified]%s; } *",attribute_aligned_object,attribute_aligned_object);
%% emit_typedef(buf,"Instance");

# Structures that inherit from <structure-stablehash>
typedef struct {
  SRECORD_HEADER
  gcv_object_t _structure_types   _attribute_aligned_object_;
  gcv_object_t stablehashcode     _attribute_aligned_object_;
  gcv_object_t other[unspecified] _attribute_aligned_object_;
} *  StablehashStructure;

# CLOS instances that inherit from <standard-stablehash>
typedef struct {
  SRECORD_HEADER
  gcv_object_t inst_class_version _attribute_aligned_object_; # indirect pointer to a CLOS-class
  gcv_object_t stablehashcode     _attribute_aligned_object_;
  gcv_object_t other[unspecified] _attribute_aligned_object_;
} *  StablehashInstance;

# Slot definitions (= instances of <slot-definition>, see clos-slotdef1.lisp
typedef struct {
  SRECORD_HEADER
  gcv_object_t inst_class_version         _attribute_aligned_object_;
  gcv_object_t slotdef_name               _attribute_aligned_object_;
  gcv_object_t slotdef_initargs           _attribute_aligned_object_;
  gcv_object_t slotdef_type               _attribute_aligned_object_;
  gcv_object_t slotdef_allocation         _attribute_aligned_object_;
  gcv_object_t slotdef_inheritable_initer _attribute_aligned_object_;
  gcv_object_t slotdef_inheritable_doc    _attribute_aligned_object_;
  # from here on only for class ⊆ <effective-slot-definition>
  gcv_object_t slotdef_location           _attribute_aligned_object_;
  gcv_object_t slotdef_efm_svuc           _attribute_aligned_object_;
  gcv_object_t slotdef_efm_ssvuc          _attribute_aligned_object_;
  gcv_object_t slotdef_efm_sbuc           _attribute_aligned_object_;
  gcv_object_t slotdef_efm_smuc           _attribute_aligned_object_;
} *  SlotDefinition;

# CLOS-Classes (= instances of <class>), see clos-class1.lisp
typedef struct {
  SRECORD_HEADER
  gcv_object_t inst_class_version       _attribute_aligned_object_; # indirect pointer to a CLOS-class
  gcv_object_t hashcode                 _attribute_aligned_object_; # GC invariant hash code
  gcv_object_t direct_methods           _attribute_aligned_object_; # set of methods that use this specializer
  gcv_object_t classname                _attribute_aligned_object_; # a symbol
  gcv_object_t direct_subclasses        _attribute_aligned_object_; # weak-list or weak-hash-table of all direct subclasses
  # from here on only for metaclass ⊆ <defined-class>
  gcv_object_t direct_superclasses      _attribute_aligned_object_; # direct superclasses
  gcv_object_t all_superclasses         _attribute_aligned_object_; # all superclasses, including itself
  gcv_object_t precedence_list          _attribute_aligned_object_; # ordered list of all superclasses
  gcv_object_t direct_slots             _attribute_aligned_object_;
  gcv_object_t slots                    _attribute_aligned_object_;
  gcv_object_t slot_location_table      _attribute_aligned_object_; # hashtable slotname -> where the slot is located
  gcv_object_t direct_default_initargs  _attribute_aligned_object_;
  gcv_object_t default_initargs         _attribute_aligned_object_;
  gcv_object_t documentation            _attribute_aligned_object_; # string or NIL
  gcv_object_t listeners                _attribute_aligned_object_; # list of objects to be notified upon a change
  gcv_object_t initialized              _attribute_aligned_object_; # describes which parts of the class are initialized
  # from here on only for metaclass ⊆ <standard-class> or metaclass ⊆ <funcallable-standard-class> or metaclass ⊆ <structure-class>
  gcv_object_t subclass_of_stablehash_p _attribute_aligned_object_; /* true if <standard-stablehash> or <structure-stablehash> is among the superclasses */
  gcv_object_t generic_accessors        _attribute_aligned_object_;
  gcv_object_t direct_accessors         _attribute_aligned_object_;
  gcv_object_t valid_initargs_from_slots _attribute_aligned_object_;
  gcv_object_t instance_size            _attribute_aligned_object_;
  # from here on only for metaclass ⊆ <standard-class> or metaclass ⊆ <funcallable-standard-class>
  gcv_object_t current_version          _attribute_aligned_object_; /* most recent class-version, points back to this class */
  gcv_object_t funcallablep             _attribute_aligned_object_;
  gcv_object_t fixed_slot_locations     _attribute_aligned_object_;
  gcv_object_t instantiated             _attribute_aligned_object_;
  gcv_object_t direct_instance_specializers _attribute_aligned_object_;
  gcv_object_t finalized_direct_subclasses _attribute_aligned_object_; /* weak-list or weak-hash-table of all finalized direct subclasses */
  gcv_object_t prototype                _attribute_aligned_object_; /* class prototype - an instance or NIL */
  # from here on only for metaclass ⊆ <standard-class>
  gcv_object_t other[unspecified]       _attribute_aligned_object_;
} *  Class;

# Length of a <defined-class>.
#define defined_class_length ((((aint)&((Class)0)->initialized-offsetofa(record_,recdata))/sizeof(gcv_object_t))+1)
# Length of a <built-in-class>.
#define built_in_class_length  (defined_class_length+1) # = clos::*<built-in-class>-instance-size*

# Closures
typedef struct {
  SRECORD_HEADER
  gcv_object_t clos_name_or_class_version _attribute_aligned_object_;
  gcv_object_t clos_codevec               _attribute_aligned_object_;
  gcv_object_t other[unspecified]         _attribute_aligned_object_;
} *  Closure;
# interpreted Closure:
typedef struct {
  SRECORD_HEADER
  gcv_object_t clos_name       _attribute_aligned_object_;
  gcv_object_t clos_form       _attribute_aligned_object_;
  gcv_object_t clos_docstring  _attribute_aligned_object_;
  gcv_object_t clos_body       _attribute_aligned_object_;
  gcv_object_t clos_var_env    _attribute_aligned_object_;
  gcv_object_t clos_fun_env    _attribute_aligned_object_;
  gcv_object_t clos_block_env  _attribute_aligned_object_;
  gcv_object_t clos_go_env     _attribute_aligned_object_;
  gcv_object_t clos_decl_env   _attribute_aligned_object_;
  gcv_object_t clos_vars       _attribute_aligned_object_;
  gcv_object_t clos_varflags   _attribute_aligned_object_;
  gcv_object_t clos_spec_anz   _attribute_aligned_object_;
  gcv_object_t clos_req_anz    _attribute_aligned_object_;
  gcv_object_t clos_opt_anz    _attribute_aligned_object_;
  gcv_object_t clos_opt_inits  _attribute_aligned_object_;
  gcv_object_t clos_key_anz    _attribute_aligned_object_;
  gcv_object_t clos_keywords   _attribute_aligned_object_;
  gcv_object_t clos_key_inits  _attribute_aligned_object_;
  gcv_object_t clos_allow_flag _attribute_aligned_object_;
  gcv_object_t clos_rest_flag  _attribute_aligned_object_;
  gcv_object_t clos_aux_anz    _attribute_aligned_object_;
  gcv_object_t clos_aux_inits  _attribute_aligned_object_;
} *  Iclosure;
#define iclos_length  ((sizeof(*(Iclosure)0)-offsetofa(record_,recdata))/sizeof(gcv_object_t))
# compiled Closure:
typedef struct {
  SRECORD_HEADER
  gcv_object_t clos_name_or_class_version _attribute_aligned_object_;
  gcv_object_t clos_codevec               _attribute_aligned_object_;
  gcv_object_t clos_consts[unspecified]   _attribute_aligned_object_; # Closure-constants
} *  Cclosure;
#define cclosure_length(ptr)  srecord_length(ptr)
#define Cclosure_length(obj)  cclosure_length(TheCclosure(obj))
# Flags in a closure. They must be disjoint from the instflags_* bits.
#ifdef TYPECODES
  #define closure_flags(ptr)  ((ptr)->recflags)
#else
  #define closure_flags(ptr)  record_flags(ptr)
#endif
#define Closure_flags(obj)  closure_flags(TheClosure(obj))
#define Cclosure_seclass(obj)  ((Closure_flags(obj) >> 4) & 0x07)
#define Cclosure_set_seclass(obj,se)  \
  (record_flags_clr(TheCclosure(obj),0x07<<4), \
   record_flags_set(TheCclosure(obj),(se)<<4))
#define closflags_instance_B  bit(7)
#define closure_instancep(ptr)  (closure_flags(ptr) & closflags_instance_B)
#define Closure_instancep(obj)  closure_instancep(TheClosure(obj))
# Closed-over environment, as a set of nested simple-vectors.
#define clos_venv  clos_consts[0]
# The function's name. Depends on whether instancep or not.
#define Closure_name(obj)  \
  (Closure_instancep(obj)             \
   ? TheCclosure(obj)->clos_consts[1] \
   : TheClosure(obj)->clos_name_or_class_version)
typedef struct {
  VRECORD_HEADER               /* self-pointer for GC, length in bits */
  /* Here: Content of the Bitvector. */
  uintW  ccv_spdepth_1;          /* maximal SP-depth, 1-part */
  uintW  ccv_spdepth_jmpbufsize; /* maximal SP-depth, jmpbufsize-part */
  uintW  ccv_numreq;             /* number of required parameters */
  uintW  ccv_numopt;             /* number of optional parameters */
  uintB  ccv_flags; /* Flags: Bit 0: &REST - parameter given?
                              Bit 1: full lambda list at the end of const vec
                              Bit 2: docstring at the end of const vec
                              Bit 3: generic function with call-inhibition?
                              Bit 4: generic function?
                              Bit 5: not used
                              Bit 6: &ALLOW-OTHER-KEYS-Flag
                              Bit 7: keyword-parameter given? */
  uintB  ccv_signature; /* abbreviated argument type, for faster FUNCALL */
  /* If keyword-parameters are given: */
  uintW  ccv_numkey;    /* Number of keyword-parameters */
  uintW  ccv_keyconsts; /* Offset in FUNC of the keywords */
} *  Codevec;
#define CCV_SPDEPTH_1           0
#define CCV_SPDEPTH_JMPBUFSIZE  2
#define CCV_NUMREQ              4
#define CCV_NUMOPT              6
#define CCV_FLAGS               8
#define CCV_SIGNATURE           9
#define CCV_NUMKEY             10
#define CCV_KEYCONSTS          12
#define CCV_START_NONKEY       10
#define CCV_START_KEY          14
/* Compiled closures, where Bit 4 has been set in the flags of clos_codevec
   are generic functions. */
%% export_def(closure_flags(ptr));
%% export_def(closure_instancep(ptr));
%% export_def(Closure_instancep(obj));

/* the position of the last const (or doc or lalist!) */
#define Cclosure_last_const(obj)  (Cclosure_length(obj) - 1 -           \
   (sizeof(*(Cclosure)0) - offsetofa(srecord_,recdata))/sizeof(gcv_object_t))
#define ccv_flags_lambda_list_p(ccv_flags)    (((ccv_flags) & bit(1)) != 0)
#define ccv_flags_documentation_p(ccv_flags)  (((ccv_flags) & bit(2)) != 0)

# A compiled LISP-function gets its arguments on the STACK
# and returns its values in MULTIPLE_VALUE_SPACE.
# It does not return a value as a C-function.
  # Return of multiple values is completely done through
  # MULTIPLE_VALUE_SPACE. As C-function: result-type Values.
  #ifndef Values
    typedef void Values;
  #endif
  # To pass a type of the value Values: return_Values(...);
  #define return_Values  return_void
  # A Lisp-function is a pointer to a C-function without returned value.
  typedef Values (*lisp_function_t)();
# If this is changed, every call of a C-function with the result type
# 'Values' (especially 'funcall', 'apply', 'eval') is to be checked.
%% puts("typedef void Values;"); /* emit_typedef useless: no sizeof(void) */
%% emit_typedef_f("Values (*%s)()","lisp_function_t");

# FSUBRs
# As C-functions: of type fsubr_function_t (no arguments, no value):
  typedef Values fsubr_function_t (void);
# The addesses of these C-functions are jumped to directly
# For SAVEMEM/LOADMEM there is a table containing all FSUBRs.
  typedef fsubr_function_t * fsubr_t;
# Signature of FSUBRs in the Lisp-way:
#         argtype          short for the argument type     fsubr_argtype_t
#         req_anz          number of required parameters   uintW
#         opt_anz          number of optional parameters   uintW
#         body_flag        Body-Flag                       fsubr_body_t
# The component body_flag contains one uintW, but we mean:
  typedef enum {
    fsubr_nobody,
    fsubr_body
  } fsubr_body_t;
# The component argtype contains a Fixnum, but it's supposed to be:
  typedef enum {
    fsubr_argtype_1_0_nobody,
    fsubr_argtype_2_0_nobody,
    fsubr_argtype_1_1_nobody,
    fsubr_argtype_2_1_nobody,
    fsubr_argtype_0_body,
    fsubr_argtype_1_body,
    fsubr_argtype_2_body
  } fsubr_argtype_t;
# conversion: see SPVW:
# extern fsubr_argtype_t fsubr_argtype (uintW req_anz, uintW opt_anz, fsubr_body_t body_flag);

# SUBRs
# SUBR table entry:
  typedef struct {
    XRECORD_HEADER
    gcv_object_t name     _attribute_aligned_object_; # name
    gcv_object_t keywords _attribute_aligned_object_; # NIL or vector with the keywords
    lisp_function_t function; # function
    uintW argtype;          # short for the argument-type
    uintW req_anz;          # number of required parameters
    uintW opt_anz;          # number of optional parameters
    uintB rest_flag;        # flag for arbitrary number of arguments
    uintB key_flag;         # flag for keywords
    uintW key_anz;          # number of keyword parameter
    uintW seclass;          /* side-effect class */
    # If necessary, add fillers here to ensure sizeof(subr_t) is a multiple of
    # varobject_alignment.
  } subr_t
    #if defined(HEAPCODES) && (alignment_long < 4) && defined(GNU)
      # Force all Subrs to be allocated with a 4-byte alignment. GC needs this.
      __attribute__ ((aligned (4)))
    #endif
    ;
  typedef subr_t *  Subr;
  # Compile-time check: sizeof(subr_t) is a multiple of varobject_alignment.
  typedef int subr_size_check[1 - 2 * (int)(sizeof(subr_t) % varobject_alignment)];
# GC needs information where objects are in here:
  #define subr_length  2
  #define subr_xlength  (sizeof(*(Subr)0)-offsetofa(record_,recdata)-subr_length*sizeof(gcv_object_t))
# the rest_flag component is a uintB, while we really mean:
  typedef enum {
    subr_norest,
    subr_rest
  } subr_rest_t;
# the key_flag component is a uintB, while we really mean:
  typedef enum {
    subr_nokey,
    subr_key,
    subr_key_allow
  } subr_key_t;
# the argtype component is a uintW, while we really mean:
  typedef enum {
    subr_argtype_0_0,
    subr_argtype_1_0,
    subr_argtype_2_0,
    subr_argtype_3_0,
    subr_argtype_4_0,
    subr_argtype_5_0,
    subr_argtype_6_0,
    subr_argtype_0_1,
    subr_argtype_1_1,
    subr_argtype_2_1,
    subr_argtype_3_1,
    subr_argtype_4_1,
    subr_argtype_0_2,
    subr_argtype_1_2,
    subr_argtype_2_2,
    subr_argtype_3_2,
    subr_argtype_0_3,
    subr_argtype_1_3,
    subr_argtype_2_3,
    subr_argtype_0_4,
    subr_argtype_0_5,
    subr_argtype_0_0_rest,
    subr_argtype_1_0_rest,
    subr_argtype_2_0_rest,
    subr_argtype_3_0_rest,
    subr_argtype_0_0_key,
    subr_argtype_1_0_key,
    subr_argtype_2_0_key,
    subr_argtype_3_0_key,
    subr_argtype_4_0_key,
    subr_argtype_0_1_key,
    subr_argtype_1_1_key,
    subr_argtype_1_2_key
  } subr_argtype_t;
# Conversion: see SPVW:
# extern subr_argtype_t subr_argtype (uintW req_anz, uintW opt_anz, subr_rest_t rest_flag, subr_key_t key_flag);
%% sprintf(buf,"struct { XRECORD_HEADER gcv_object_t name%s; gcv_object_t keywords%s; lisp_function_t function; uintW argtype; uintW req_anz; uintW opt_anz; uintB rest_flag; uintB key_flag; uintW key_anz; uintW seclass; } %%s",attribute_aligned_object,attribute_aligned_object);
%% #if defined(HEAPCODES) && (alignment_long < 4) && defined(GNU)
%%   strcat(buf," __attribute__ ((aligned (4)))");
%% #endif
%% emit_typedef_f(buf,"subr_t");
%% emit_typedef("subr_t *","Subr");
%% emit_typedef("enum { subr_norest, subr_rest }","subr_rest_t");
%% emit_typedef("enum { subr_nokey, subr_key, subr_key_allow }","subr_key_t");

/* side-effect class is really seclass_t: */
typedef enum {
  seclass_foldable, /* the function allows Constant-Folding:
     two calls with identical arguments give the same result,
     and calls with constant arguments can be evaluated at compile time.
     In particular, no side effects, do not depend on global variables or such,
     do not even look "inside" their arguments */
  seclass_no_se, /* no side effects, do not depend on global variables or such,
     do not even look "inside" their arguments, but not "foldable". */
  seclass_read, /* no side effects, but depend on global variables
     or look "inside" their arguments. */
  seclass_write, /* only side effects: does not read anything,
     just sets some global variables. */
  seclass_default /* may do side effects */
} seclass_t;
%% puts("enum { seclass_foldable, seclass_no_se, seclass_read, seclass_write, seclass_default};");

# Small-Read-Label
#ifdef TYPECODES
  #define make_small_read_label(n)  \
    type_data_object(system_type, ((uintV)(n)<<1) + bit(0))
  #define small_read_label_integer_p(obj)  \
    (posfixnump(obj) && (posfixnum_to_V(obj) < vbit(oint_data_len-2)))
  #define small_read_label_value(obj)  \
    fixnum((as_oint(obj) >> (oint_data_shift+1)) & (vbit(oint_data_len-2)-1))
#else
  #define make_small_read_label(n)  \
    type_data_object(small_read_label_type, (uintV)(n))
  #define small_read_label_integer_p(obj)  posfixnump(obj)
  #define small_read_label_value(obj)  \
    fixnum((as_oint(obj) >> oint_data_shift) & (vbit(oint_data_len)-1))
#endif

# Machine pointers:
# make_machine(ptr)
#ifdef TYPECODES
  #define make_machine(ptr)  type_pointer_object(machine_type,ptr)
#else
  #if defined(WIDE_AUXI)
    #define make_machine(ptr)  as_object_with_auxi((aint)(ptr)+machine_bias)
  #else
    #define make_machine(ptr)  as_object((oint)(ptr)+machine_bias)
  #endif
#endif
%% export_def(make_machine(ptr));

# Pointer to machine code
# make_machine_code(ptr)
#if defined(TYPECODES) || (log2_C_CODE_ALIGNMENT >= 2)
  #define make_machine_code(ptr)  make_machine(ptr)
#elif defined(HPPA)
  #define make_machine_code(ptr)  make_machine((uintP)(ptr)&~(uintP)3)
#else
  #define make_machine_code(ptr)  make_machine((uintP)(ptr)<<(2-log2_C_CODE_ALIGNMENT))
#endif

# System-Pointer
#define make_system(data)  \
  type_data_object(system_type, vbit(oint_data_len-1) | bit(0) | ((vbitm(oint_data_len)-1) & (data)))
# all such go into the special print routine io.d:pr_system()
%% export_def(make_system(data));

# missing value
#define unbound  make_system(0xFFFFFFUL)
%% export_def(unbound);

# missing object (internal use only):
#define nullobj  make_machine(0)  # = as_object((oint)0)
#ifdef DEBUG_GCSAFETY
  #define gcv_nullobj  (gcv_object_t)nullobj
#else
  #define gcv_nullobj  nullobj
#endif
%% export_def(nullobj);
%% export_def(gcv_nullobj);


# cgci_pointable(obj)  converts a certainly GC-invariant object to an aint.
# pgci_pointable(obj)  converts a possibly GC-invariant object to an aint.
# ngci_pointable(obj)  converts a not GC-invariant object to an aint.
#if defined(DEBUG_GCSAFETY)
  static inline aint cgci_pointable (object obj) {
    return obj.one_o;
  }
  static inline aint cgci_pointable (gcv_object_t obj) {
    return obj.one_o;
  }
  static inline aint pgci_pointable (object obj) {
    if (!(gcinvariant_object_p(obj) || gcinvariant_symbol_p(obj)
          || obj.allocstamp == alloccount || nonimmsubrp(obj)))
      abort();
    nonimmprobe(obj.one_o);
    return obj.one_o;
  }
  static inline aint pgci_pointable (gcv_object_t obj) {
    nonimmprobe(obj.one_o);
    return obj.one_o;
  }
  static inline aint ngci_pointable (object obj) {
    if (!(gcinvariant_symbol_p(obj)
          || obj.allocstamp == alloccount || nonimmsubrp(obj)))
      abort();
    nonimmprobe(obj.one_o);
    return obj.one_o;
  }
  static inline aint ngci_pointable (gcv_object_t obj) {
    nonimmprobe(obj.one_o);
    return obj.one_o;
  }
#elif defined(WIDE_AUXI)
  #define cgci_pointable(obj)  (obj).one_o
  #define pgci_pointable(obj)  (obj).one_o
  #define ngci_pointable(obj)  (obj).one_o
#else
  #define cgci_pointable(obj)  as_oint(obj)
  #define pgci_pointable(obj)  as_oint(obj)
  #define ngci_pointable(obj)  as_oint(obj)
#endif
%% #if defined(DEBUG_GCSAFETY)
%%   puts("static inline aint cgci_pointable (object obj) { return obj.one_o; }");
%%   puts("static inline aint cgci_pointable (gcv_object_t obj) { return obj.one_o; }");
%%   puts("static inline aint pgci_pointable (object obj) { if (!(gcinvariant_object_p(obj) || gcinvariant_symbol_p(obj) || obj.allocstamp == alloccount || nonimmsubrp(obj))) abort(); nonimmprobe(obj.one_o); return obj.one_o; }");
%%   puts("static inline aint pgci_pointable (gcv_object_t obj) { nonimmprobe(obj.one_o); return obj.one_o; }");
%%   puts("static inline aint ngci_pointable (object obj) { if (!(gcinvariant_symbol_p(obj) || obj.allocstamp == alloccount || nonimmsubrp(obj))) abort(); nonimmprobe(obj.one_o); return obj.one_o; }");
%%   puts("static inline aint ngci_pointable (gcv_object_t obj) { nonimmprobe(obj.one_o); return obj.one_o; }");
%% #else
%%   export_def(cgci_pointable(obj));
%%   export_def(pgci_pointable(obj));
%%   export_def(ngci_pointable(obj));
%% #endif

# TheCons(object) yields the Cons that's equivalent to object.
# The information that it is a Cons has to be put into it.
# The other type conversions are similar.
#ifdef TYPECODES
  #ifdef DEBUG_GCSAFETY
    #define cgci_types_pointable(ORed_types,obj)  cgci_pointable(obj)
    #define pgci_types_pointable(ORed_types,obj)  pgci_pointable(obj)
    #define ngci_types_pointable(ORed_types,obj)  ngci_pointable(obj)
  #else
    #define cgci_types_pointable(ORed_types,obj)  types_pointable(ORed_types,obj)
    #define pgci_types_pointable(ORed_types,obj)  types_pointable(ORed_types,obj)
    #define ngci_types_pointable(ORed_types,obj)  types_pointable(ORed_types,obj)
  #endif
  #define TheCons(obj)  ((Cons)(ngci_types_pointable(cons_type,obj)))
  #define TheRatio(obj)  ((Ratio)(ngci_types_pointable(ratio_type|bit(sign_bit_t),obj)))
  #define TheComplex(obj)  ((Complex)(ngci_types_pointable(complex_type,obj)))
  #define TheSymbol(obj)  ((Symbol)(ngci_types_pointable(symbol_type,obj)))
  #if (oint_symbolflags_shift==oint_type_shift)
  #define TheSymbolflagged(obj)  ((Symbol)(ngci_types_pointable(symbol_type|bit(active_bit)|bit(dynam_bit)|bit(svar_bit),obj)))
  #else
  #define TheSymbolflagged(obj)  TheSymbol(symbol_without_flags(obj))
  #endif
  #define TheBignum(obj)  ((Bignum)(ngci_types_pointable(bignum_type|bit(sign_bit_t),obj)))
  #ifndef IMMEDIATE_FFLOAT
  #define TheFfloat(obj)  ((Ffloat)(ngci_types_pointable(ffloat_type|bit(sign_bit_t),obj)))
  #endif
  #define TheDfloat(obj)  ((Dfloat)(ngci_types_pointable(dfloat_type|bit(sign_bit_t),obj)))
  #define TheLfloat(obj)  ((Lfloat)(ngci_types_pointable(lfloat_type|bit(sign_bit_t),obj)))
  #define TheSarray(obj)  ((Sarray)(ngci_types_pointable(sbvector_type|sb2vector_type|sb4vector_type|sb8vector_type|sb16vector_type|sb32vector_type|sstring_type|svector_type,obj)))
  #define TheSbvector(obj)  ((Sbvector)(ngci_types_pointable(sbvector_type|sb2vector_type|sb4vector_type|sb8vector_type|sb16vector_type|sb32vector_type,obj)))
  #define TheCodevec(obj)  ((Codevec)(ngci_types_pointable(sb8vector_type,obj)))
  #define TheS8string(obj)  ((S8string)(ngci_types_pointable(sstring_type,obj)))
  #define TheS16string(obj)  ((S16string)(ngci_types_pointable(sstring_type,obj)))
  #define TheS32string(obj)  ((S32string)(ngci_types_pointable(sstring_type,obj)))
  #define TheSnstring(obj)  ((Snstring)(ngci_types_pointable(sstring_type,obj)))
  #define TheSistring(obj)  ((Sistring)(ngci_types_pointable(sstring_type,obj)))
  #define TheSstring(obj)  ((Sstring)(ngci_types_pointable(sstring_type,obj)))
  #define TheSvector(obj)  ((Svector)(ngci_types_pointable(svector_type,obj)))
  #define TheIarray(obj)  ((Iarray)(ngci_types_pointable(mdarray_type|bvector_type|b2vector_type|b4vector_type|b8vector_type|b16vector_type|b32vector_type|string_type|vector_type,obj)))
  #define TheRecord(obj)  ((Record)(ngci_types_pointable(closure_type|structure_type|stream_type|orecord_type|instance_type|lrecord_type,obj)))
  #define TheLrecord(obj)  ((Lrecord)(ngci_types_pointable(lrecord_type,obj)))
  #define TheSrecord(obj)  ((Srecord)(ngci_types_pointable(closure_type|structure_type|orecord_type|instance_type,obj)))
  #define TheXrecord(obj)  ((Xrecord)(ngci_types_pointable(stream_type|orecord_type,obj)))
  #define ThePackage(obj)  ((Package)(ngci_types_pointable(orecord_type,obj)))
  #define TheHashtable(obj)  ((Hashtable)(ngci_types_pointable(orecord_type,obj)))
  #define TheReadtable(obj)  ((Readtable)(ngci_types_pointable(orecord_type,obj)))
  #define ThePathname(obj)  ((Pathname)(ngci_types_pointable(orecord_type,obj)))
  #ifdef LOGICAL_PATHNAMES
  #define TheLogpathname(obj)  ((Logpathname)(ngci_types_pointable(orecord_type,obj)))
  #endif
  #define The_Random_state(obj)  ((Random_state)(ngci_types_pointable(orecord_type,obj)))
  #define TheByte(obj)  ((Byte)(ngci_types_pointable(orecord_type,obj)))
  #define TheFsubr(obj)  ((Fsubr)(ngci_types_pointable(orecord_type,obj)))
  #define TheLoadtimeeval(obj)  ((Loadtimeeval)(ngci_types_pointable(orecord_type,obj)))
  #define TheSymbolmacro(obj)  ((Symbolmacro)(ngci_types_pointable(orecord_type,obj)))
  #define TheGlobalSymbolmacro(obj)  ((GlobalSymbolmacro)(ngci_types_pointable(orecord_type,obj)))
  #define TheMacro(obj)  ((Macro)(ngci_types_pointable(orecord_type,obj)))
  #define TheFunctionMacro(obj)  ((FunctionMacro)(ngci_types_pointable(orecord_type,obj)))
  #define TheBigReadLabel(obj)  ((BigReadLabel)(ngci_types_pointable(orecord_type,obj)))
  #define TheEncoding(obj)  ((Encoding)(ngci_types_pointable(orecord_type,obj)))
  #ifdef FOREIGN
  #define TheFpointer(obj)  ((Fpointer)(ngci_types_pointable(orecord_type,obj)))
  #endif
  #ifdef DYNAMIC_FFI
  #define TheFaddress(obj)  ((Faddress)(ngci_types_pointable(orecord_type,obj)))
  #define TheFvariable(obj)  ((Fvariable)(ngci_types_pointable(orecord_type,obj)))
  #define TheFfunction(obj)  ((Ffunction)(ngci_types_pointable(orecord_type,obj)))
  #endif
  #define TheWeakpointer(obj)  ((Weakpointer)(ngci_types_pointable(orecord_type,obj)))
  #define TheMutableWeakList(obj)  ((MutableWeakList)(ngci_types_pointable(orecord_type,obj)))
  #define TheWeakList(obj)  ((WeakList)(ngci_types_pointable(lrecord_type,obj)))
  #define TheWeakAnd(obj)  ((WeakAnd)(ngci_types_pointable(lrecord_type,obj)))
  #define TheWeakOr(obj)  ((WeakOr)(ngci_types_pointable(lrecord_type,obj)))
  #define TheWeakmapping(obj)  ((Weakmapping)(ngci_types_pointable(orecord_type,obj)))
  #define TheWeakAndMapping(obj)  ((WeakAndMapping)(ngci_types_pointable(lrecord_type,obj)))
  #define TheWeakOrMapping(obj)  ((WeakOrMapping)(ngci_types_pointable(lrecord_type,obj)))
  #define TheMutableWeakAlist(obj)  ((MutableWeakAlist)(ngci_types_pointable(orecord_type,obj)))
  #define TheWeakAlist(obj)  ((WeakAlist)(ngci_types_pointable(lrecord_type,obj)))
  #define TheWeakHashedAlist(obj)  ((WeakHashedAlist)(ngci_types_pointable(lrecord_type,obj)))
  #define TheFinalizer(obj)  ((Finalizer)(ngci_types_pointable(orecord_type,obj)))
  #ifdef SOCKET_STREAMS
  #define TheSocketServer(obj) ((Socket_server)(ngci_types_pointable(orecord_type,obj)))
  #endif
  #ifdef YET_ANOTHER_RECORD
  #define TheYetanother(obj)  ((Yetanother)(ngci_types_pointable(orecord_type,obj)))
  #endif
  #define TheStream(obj)  ((Stream)(ngci_types_pointable(stream_type,obj)))
  #define TheStructure(obj)  ((Structure)(ngci_types_pointable(structure_type,obj)))
  #define TheClosure(obj)  ((Closure)(ngci_types_pointable(closure_type,obj)))
  #define TheIclosure(obj)  ((Iclosure)(ngci_types_pointable(closure_type,obj)))
  #define TheCclosure(obj)  ((Cclosure)(ngci_types_pointable(closure_type,obj)))
  #define TheInstance(obj)  ((Instance)(ngci_types_pointable(instance_type|closure_type,obj)))
  #define TheSubr(obj)  ((Subr)(cgci_types_pointable(subr_type,obj)))
  #define TheFramepointer(obj)  ((gcv_object_t*)(cgci_types_pointable(system_type,obj)))
  #define TheMachine(obj)  ((void*)(cgci_types_pointable(machine_type,obj)))
  #define TheMachineCode(obj)  TheMachine(obj)
  #define ThePseudofun(obj)  ((Pseudofun)TheMachineCode(obj))
  #ifdef FOREIGN_HANDLE
  # pack Handle in Sbvector
  #define TheHandle(obj)  (*(Handle*)(&TheSbvector(obj)->data[0]))
  #else
  # pack Handle in Fixnum>=0
  #define TheHandle(obj)  ((Handle)posfixnum_to_V(obj))
  #endif
  # variable length object:
  #define TheVarobject(obj)                                              \
    ((Varobject)                                                         \
     (ngci_types_pointable                                               \
      (sbvector_type|sb2vector_type|sb4vector_type|sb8vector_type        \
         |sb16vector_type|sb32vector_type                                \
       |sstring_type|svector_type                                        \
       |mdarray_type                                                     \
       |bvector_type|b2vector_type|b4vector_type|b8vector_type           \
         |b16vector_type|b32vector_type                                  \
       |string_type|vector_type                                          \
       |closure_type|structure_type|stream_type|orecord_type             \
       |instance_type|lrecord_type|symbol_type                           \
       |bignum_type|ffloat_type|dfloat_type|lfloat_type|bit(sign_bit_t), \
       obj                                                               \
    )))
  # Object that represents a pointer into the memory:
  #define ThePointer(obj)                                               \
    (pgci_types_pointable                                               \
     (sbvector_type|sb2vector_type|sb4vector_type|sb8vector_type        \
        |sb16vector_type|sb32vector_type                                \
      |sstring_type|svector_type                                        \
      |mdarray_type                                                     \
      |bvector_type|b2vector_type|b4vector_type|b8vector_type           \
        |b16vector_type|b32vector_type                                  \
      |string_type|vector_type                                          \
      |closure_type|structure_type|stream_type|orecord_type             \
      |instance_type|lrecord_type|symbol_type                           \
      |cons_type                                                        \
      |bignum_type|ffloat_type|dfloat_type|lfloat_type                  \
      |ratio_type|complex_type|bit(sign_bit_t),                         \
      obj                                                               \
    ))
#else # no TYPECODES
  #define TheCons(obj)  ((Cons)(ngci_pointable(obj)-cons_bias))
  #define TheRatio(obj)  ((Ratio)(ngci_pointable(obj)-varobject_bias))
  #define TheComplex(obj)  ((Complex)(ngci_pointable(obj)-varobject_bias))
  #define TheSymbol(obj)  ((Symbol)(ngci_pointable(obj)-varobject_bias))
  #define TheSymbolflagged(obj)  TheSymbol(symbol_without_flags(obj))
  #define TheBignum(obj)  ((Bignum)(ngci_pointable(obj)-varobject_bias))
  #define TheFfloat(obj)  ((Ffloat)(ngci_pointable(obj)-varobject_bias))
  #define TheDfloat(obj)  ((Dfloat)(ngci_pointable(obj)-varobject_bias))
  #define TheLfloat(obj)  ((Lfloat)(ngci_pointable(obj)-varobject_bias))
  #define TheSarray(obj)  ((Sarray)(ngci_pointable(obj)-varobject_bias))
  #define TheSbvector(obj)  ((Sbvector)(ngci_pointable(obj)-varobject_bias))
  #define TheCodevec(obj)  ((Codevec)TheSbvector(obj))
  #define TheS8string(obj)  ((S8string)(ngci_pointable(obj)-varobject_bias))
  #define TheS16string(obj)  ((S16string)(ngci_pointable(obj)-varobject_bias))
  #define TheS32string(obj)  ((S32string)(ngci_pointable(obj)-varobject_bias))
  #define TheSnstring(obj)  ((Snstring)(ngci_pointable(obj)-varobject_bias))
  #define TheSistring(obj)  ((Sistring)(ngci_pointable(obj)-varobject_bias))
  #define TheSstring(obj)  ((Sstring)(ngci_pointable(obj)-varobject_bias))
  #define TheSvector(obj)  ((Svector)(ngci_pointable(obj)-varobject_bias))
  #define TheIarray(obj)  ((Iarray)(ngci_pointable(obj)-varobject_bias))
  #define TheRecord(obj)  ((Record)(ngci_pointable(obj)-varobject_bias))
  #define TheLrecord(obj)  ((Lrecord)(ngci_pointable(obj)-varobject_bias))
  #define TheSrecord(obj)  ((Srecord)(ngci_pointable(obj)-varobject_bias))
  #define TheXrecord(obj)  ((Xrecord)(ngci_pointable(obj)-varobject_bias))
  #define ThePackage(obj)  ((Package)(ngci_pointable(obj)-varobject_bias))
  #define TheHashtable(obj)  ((Hashtable)(ngci_pointable(obj)-varobject_bias))
  #define TheReadtable(obj)  ((Readtable)(ngci_pointable(obj)-varobject_bias))
  #define ThePathname(obj)  ((Pathname)(ngci_pointable(obj)-varobject_bias))
  #ifdef LOGICAL_PATHNAMES
  #define TheLogpathname(obj)  ((Logpathname)(ngci_pointable(obj)-varobject_bias))
  #endif
  #define The_Random_state(obj)  ((Random_state)(ngci_pointable(obj)-varobject_bias))
  #define TheByte(obj)  ((Byte)(ngci_pointable(obj)-varobject_bias))
  #define TheFsubr(obj)  ((Fsubr)(ngci_pointable(obj)-varobject_bias))
  #define TheLoadtimeeval(obj)  ((Loadtimeeval)(ngci_pointable(obj)-varobject_bias))
  #define TheSymbolmacro(obj)  ((Symbolmacro)(ngci_pointable(obj)-varobject_bias))
  #define TheGlobalSymbolmacro(obj)  ((GlobalSymbolmacro)(ngci_pointable(obj)-varobject_bias))
  #define TheMacro(obj)  ((Macro)(ngci_pointable(obj)-varobject_bias))
  #define TheFunctionMacro(obj)  ((FunctionMacro)(ngci_pointable(obj)-varobject_bias))
  #define TheBigReadLabel(obj)  ((BigReadLabel)(ngci_pointable(obj)-varobject_bias))
  #define TheEncoding(obj)  ((Encoding)(ngci_pointable(obj)-varobject_bias))
  #ifdef FOREIGN
  #define TheFpointer(obj)  ((Fpointer)(ngci_pointable(obj)-varobject_bias))
  #endif
  #ifdef DYNAMIC_FFI
  #define TheFaddress(obj)  ((Faddress)(ngci_pointable(obj)-varobject_bias))
  #define TheFvariable(obj)  ((Fvariable)(ngci_pointable(obj)-varobject_bias))
  #define TheFfunction(obj)  ((Ffunction)(ngci_pointable(obj)-varobject_bias))
  #endif
  #define TheWeakpointer(obj)  ((Weakpointer)(ngci_pointable(obj)-varobject_bias))
  #define TheMutableWeakList(obj)  ((MutableWeakList)(ngci_pointable(obj)-varobject_bias))
  #define TheWeakList(obj)  ((WeakList)(ngci_pointable(obj)-varobject_bias))
  #define TheWeakAnd(obj)  ((WeakAnd)(ngci_pointable(obj)-varobject_bias))
  #define TheWeakOr(obj)  ((WeakOr)(ngci_pointable(obj)-varobject_bias))
  #define TheWeakmapping(obj)  ((Weakmapping)(ngci_pointable(obj)-varobject_bias))
  #define TheWeakAndMapping(obj)  ((WeakAndMapping)(ngci_pointable(obj)-varobject_bias))
  #define TheWeakOrMapping(obj)  ((WeakOrMapping)(ngci_pointable(obj)-varobject_bias))
  #define TheMutableWeakAlist(obj)  ((MutableWeakAlist)(ngci_pointable(obj)-varobject_bias))
  #define TheWeakAlist(obj)  ((WeakAlist)(ngci_pointable(obj)-varobject_bias))
  #define TheWeakHashedAlist(obj)  ((WeakHashedAlist)(ngci_pointable(obj)-varobject_bias))
  #define TheFinalizer(obj)  ((Finalizer)(ngci_pointable(obj)-varobject_bias))
  #ifdef SOCKET_STREAMS
  #define TheSocketServer(obj) ((Socket_server)(ngci_pointable(obj)-varobject_bias))
  #endif
  #ifdef YET_ANOTHER_RECORD
  #define TheYetanother(obj)  ((Yetanother)(ngci_pointable(obj)-varobject_bias))
  #endif
  #define TheStream(obj)  ((Stream)(ngci_pointable(obj)-varobject_bias))
  #define TheStructure(obj)  ((Structure)(ngci_pointable(obj)-varobject_bias))
  #define TheClosure(obj)  ((Closure)(ngci_pointable(obj)-varobject_bias))
  #define TheIclosure(obj)  ((Iclosure)(ngci_pointable(obj)-varobject_bias))
  #define TheCclosure(obj)  ((Cclosure)(ngci_pointable(obj)-varobject_bias))
  #define TheInstance(obj)  ((Instance)(ngci_pointable(obj)-varobject_bias))
  #define TheSubr(obj)  ((Subr)(cgci_pointable(obj)-subr_bias))
  #define TheFramepointer(obj)  ((gcv_object_t*)(cgci_pointable(obj)-machine_bias))
  #define TheMachine(obj)  ((void*)(cgci_pointable(obj)-machine_bias))
  #if (log2_C_CODE_ALIGNMENT >= 2)
    #define TheMachineCode(obj)  TheMachine(obj)
  #elif defined(HPPA)
    #define TheMachineCode(obj)  ((void*)((uintP)TheMachine(obj)+2))
  #else
    #define TheMachineCode(obj)  ((void*)(((uintP)TheMachine(obj)>>(2-log2_C_CODE_ALIGNMENT))|(CODE_ADDRESS_RANGE&~((~(uintP)0)>>(2-log2_C_CODE_ALIGNMENT)))))
  #endif
  #define ThePseudofun(obj)  ((Pseudofun)TheMachineCode(obj))
  #ifdef FOREIGN_HANDLE
  # pack Handle in Sbvector
  #define TheHandle(obj)  (*(Handle*)(&TheSbvector(obj)->data[0]))
  #else
  # pack Handle in Fixnum>=0
  #define TheHandle(obj)  ((Handle)posfixnum_to_V(obj))
  #endif
  # Object of variable length:
  #define TheVarobject(obj)  ((Varobject)(ngci_pointable(obj)-varobject_bias))
  # Object, represents a pointer into the memory:
  #define ThePointer(obj)  ((void*)(pgci_pointable(obj) & ~(aint)nonimmediate_bias_mask))
#endif
#define TheClassVersion(obj)  ((ClassVersion)TheSvector(obj))
#define TheSlotDefinition(obj)  ((SlotDefinition)TheInstance(obj))
#define TheClass(obj)  ((Class)TheInstance(obj))
%% export_def(TheCons(obj));
%% #if notused
%%   export_def(TheRatio(obj));
%%   export_def(TheComplex(obj));
%% #endif
%% export_def(TheSymbol(obj));
%% export_def(TheBignum(obj));
%% #if notused
%%   export_def(TheSarray(obj));
%% #endif
%% export_def(TheSbvector(obj));
%% #ifdef HAVE_SMALL_SSTRING
%%   export_def(TheS8string(obj));
%%   export_def(TheS16string(obj));
%%   export_def(TheS32string(obj));
%% #endif
%% export_def(TheSstring(obj));
%% export_def(TheSvector(obj));
%% export_def(TheRecord(obj));
%% export_def(TheSrecord(obj));
%% #if notused
%%   export_def(TheXrecord(obj));
%%   export_def(ThePackage(obj));
%% #endif
%% export_def(TheEncoding(obj));
%% #ifdef FOREIGN
%%   export_def(TheFpointer(obj));
%% #endif
%% export_def(TheStructure(obj));
%% export_def(TheClosure(obj));
%% export_def(TheInstance(obj));
%% export_def(TheSubr(obj));
%% export_def(TheMachine(obj));
%% export_def(TheMachineCode(obj));
%% export_def(ThePseudofun(obj));

# Some acronyms
# Access to objects that are conses:
#define Car(obj)  (TheCons(obj)->car)
#define Cdr(obj)  (TheCons(obj)->cdr)
# Access to objects that are symbols:
#define Symbol_value(obj)  (TheSymbol(obj)->symvalue)
#define Symbol_function(obj)  (TheSymbol(obj)->symfunction)
#define Symbol_plist(obj)  (TheSymbol(obj)->proplist)
#define Symbol_name(obj)  (TheSymbol(obj)->pname)
#define Symbol_package(obj)  (TheSymbol(obj)->homepackage)
# Length (number of objects) of a record, obj has to be a Srecord/Xrecord:
#define SXrecord_length(obj)  \
  (Record_type(obj) < rectype_limit ? Srecord_length(obj) : Xrecord_length(obj))
# Likewise, but ignoring weak pointers:
#define SXrecord_nonweak_length(obj)  \
  (Record_type(obj) < rectype_limit              \
   ? Srecord_length(obj)                         \
   : ((Record_type(obj)==Rectype_Weakpointer     \
       || Record_type(obj)==Rectype_Weakmapping) \
      ? 0                                        \
      : Xrecord_length(obj)))
# Length of an Lrecord, ignoring weak pointers:
#define Lrecord_nonweak_length(obj)  \
  ((Record_type(obj) >= Rectype_WeakList                 \
    && Record_type(obj) <= Rectype_WeakHashedAlist_Both) \
   ? 0                                                   \
   : Lrecord_length(obj))
# Length (number of objects) of a record, obj has to be a Record:
#define Record_length(obj)  \
  (Record_type(obj) >= rectype_longlimit \
   ? Lrecord_length(obj)                 \
   : SXrecord_length(obj))
# Likewise, but ignoring weak pointers:
#define Record_nonweak_length(obj)  \
  (Record_type(obj) >= rectype_longlimit \
   ? Lrecord_nonweak_length(obj)         \
   : SXrecord_nonweak_length(obj))
%% export_def(Car(obj));
%% export_def(Cdr(obj));
%% export_def(Symbol_value(obj));
%% export_def(Symbol_function(obj));
%% export_def(Symbol_plist(obj));
%% export_def(Symbol_name(obj));
%% export_def(Symbol_package(obj));


# ####################### type test predicates ############################### #
# There are two kinds of predicates:
# 1.  ???p, query with 'if':  if ???p(object)
# 2.  if_???p, called as
#         if_???p(object, statement1, statement2)
#       instead of
#         if ???p(object) statement1 else statement2

# UP: tests for equality of pointers EQ
# eq(obj1,obj2)
# > obj1,obj2: Lisp-objects
# < result: true, if objects are equal
#if defined(DEBUG_GCSAFETY)
  #define eq(obj1,obj2)  (pgci_pointable(obj1) == pgci_pointable(obj2))
#elif defined(WIDE_STRUCT) || defined(OBJECT_STRUCT)
  #define eq(obj1,obj2)  (as_oint(obj1) == as_oint(obj2))
#elif defined(WIDE_AUXI)
  #define eq(obj1,obj2)  ((obj1).one_o == (obj2).one_o)
#else
  #define eq(obj1,obj2)  ((obj1) == (obj2))
#endif
%% export_def(eq(obj1,obj2));

# Test for NIL
#define nullp(obj)  (eq(obj,NIL))
%% export_def(nullp(obj));

# Shorthand: Test a fixed symbol's value for NIL
#define nullpSv(sym) ( nullp(Symbol_value(S(sym))))

# Test for an argument's value, whether the argument was provided.
#define boundp(obj)  (!eq(obj,unbound))
%% export_def(boundp(obj));

# Test for an argument's value, whether the argument was not provided or NIL.
#define missingp(obj)  (!boundp(obj) || nullp(obj))
%% export_def(missingp(obj));

# Test for Cons
#ifdef TYPECODES
  #if defined(cons_bit_o)
    # define consp(obj)  (as_oint(obj) & wbit(cons_bit_o))
    #define consp(obj)  (wbit_test(as_oint(obj),cons_bit_o))
    #ifdef fast_mtypecode
      #ifdef WIDE_STRUCT
        #undef consp
        #define consp(obj)  (typecode(obj) & bit(cons_bit_t))
      #endif
      #define mconsp(obj)  (mtypecode(obj) & bit(cons_bit_t))
    #else
      #define mconsp(obj)  consp(obj)
    #endif
  #else
    #define consp(obj)  (typecode(obj) == cons_type)
    #define mconsp(obj)  (mtypecode(obj) == cons_type)
  #endif
#else
  #define consp(obj)  ((as_oint(obj) & 7) == (cons_bias+conses_misaligned))
  #define mconsp(obj)  consp(obj)
#endif
%% export_def(consp(obj));
%% export_def(mconsp(obj));

# Test for Atom
#ifdef TYPECODES
  #if defined(cons_bit_o)
    # define atomp(obj)  ((as_oint(obj) & wbit(cons_bit_o))==0)
    #define atomp(obj)  (!wbit_test(as_oint(obj),cons_bit_o))
    #ifdef fast_mtypecode
      #ifdef WIDE_STRUCT
        #undef atomp
        #define atomp(obj)  ((typecode(obj) & bit(cons_bit_t))==0)
      #endif
      #define matomp(obj)  ((mtypecode(obj) & bit(cons_bit_t))==0)
    #else
      #define matomp(obj)  atomp(obj)
    #endif
  #else
    #define atomp(obj)  (!(typecode(obj) == cons_type))
    #define matomp(obj)  (!(mtypecode(obj) == cons_type))
  #endif
#else
  #define atomp(obj)  (!consp(obj))
  #define matomp(obj)  atomp(obj)
#endif
%% export_def(atomp(obj));
%% export_def(matomp(obj));

# For all type tests below this line, the argument must be side-effect-free.
# Ideally a variable, but a STACK_(n) reference works as well.

# Test for List
#define listp(obj)  (nullp(obj) || consp(obj))
%% export_def(listp(obj));

#ifndef TYPECODES
  # Test for Object with variable length
  #define varobjectp(obj)  ((as_oint(obj) & nonimmediate_heapcode_mask) == (varobject_bias+varobjects_misaligned))
#endif
%% #ifndef TYPECODES
%%   export_def(varobjectp(obj));
%% #endif

# Test for Symbol
#ifdef TYPECODES
  #if defined(symbol_bit_o)
    # define symbolp(obj)  (as_oint(obj) & wbit(symbol_bit_o))
    #define symbolp(obj)  (wbit_test(as_oint(obj),symbol_bit_o))
    #ifdef WIDE_STRUCT
      #undef symbolp
      #define symbolp(obj)  (typecode(obj) & bit(symbol_bit_t))
    #endif
  #else
    #define symbolp(obj)  (typecode(obj) == symbol_type)
  #endif
#else
  #define symbolp(obj)  \
    (varobjectp(obj) && (Record_type(obj) == Rectype_Symbol))
#endif
%% export_def(symbolp(obj));

# Test for number
#ifdef TYPECODES
  # define numberp(obj)  (as_oint(obj) & wbit(number_bit_o))
  #define numberp(obj)  (wbit_test(as_oint(obj),number_bit_o))
  #ifdef WIDE_STRUCT
    #undef numberp
    #define numberp(obj)  (typecode(obj) & bit(number_bit_t))
  #endif
#else
  #define immediate_number_p(obj)  \
    ((as_oint(obj) & ((4 << imm_type_shift) | immediate_bias)) == (fixnum_type&sfloat_type))
  #define numberp(obj)  \
    (immediate_number_p(obj) \
     || (varobjectp(obj)     \
         && ((uintB)(Record_type(obj)-Rectype_Bignum) <= Rectype_Complex-Rectype_Bignum)))
#endif
%% #if notused
%% #ifdef TYPECODES
%%  export_def(numberp(obj));
%% #else
%%  export_def(immediate_number_p(obj));
%% #endif
%% #endif

# Test for Vector (typebytes %001,%010,%011,%101,%110,%111)
#ifdef TYPECODES
  #define vectorp(obj)  \
    ((tint)(typecode(obj) - sbvector_type) <= (tint)(vector_type - sbvector_type))
#else
  # cases: Rectype_Sbvector, Rectype_Sb[2|4|8|16|32]vector, Rectype_Svector, Rectype_[Imm_]S[8|16|32]string,
  #        Rectype_bvector, Rectype_b[2|4|8|16|32]vector, Rectype_vector, Rectype_reallocstring, Rectype_string
  #define vectorp(obj)  \
    (varobjectp(obj) && ((uintB)(Record_type(obj) - Rectype_vector) \
                         <= Rectype_string - Rectype_vector))
#endif
%% export_def(vectorp(obj));

# Test for simple-vector or simple-bit-vector or simple-string
#ifdef TYPECODES
  #define simplep(obj)  \
    ((tint)(typecode(obj) - sbvector_type) <= (tint)(svector_type - sbvector_type))
#else
  # cases: Rectype_Sbvector, Rectype_Sb[2|4|8|16|32]vector, Rectype_Svector, Rectype_[Imm_]S[8|16|32]string,
  #        Rectype_reallocstring
  #define simplep(obj)  \
    (varobjectp(obj) && ((uintB)(Record_type(obj) - Rectype_Svector) \
                         <= Rectype_reallocstring - Rectype_Svector))
#endif

# Tests an Array for simple-vector or simple-bit-vector or simple-string
#ifdef TYPECODES
  #define array_simplep(obj)  \
    ((typecode(obj) & bit(notsimple_bit_t)) == 0)
#else
  # cases: Rectype_Sbvector, Rectype_Sb[2|4|8|16|32]vector, Rectype_Svector, Rectype_[Imm_]S[8|16|32]string,
  #        Rectype_reallocstring
  #define array_simplep(obj)  \
    ((uintB)(Record_type(obj) - Rectype_Svector) \
     <= Rectype_reallocstring - Rectype_Svector)
#endif

# Test for simple-vector
#ifdef TYPECODES
  #define simple_vector_p(obj)  \
    (typecode(obj) == svector_type)
#else
  # cases: Rectype_Svector
  #define simple_vector_p(obj)  \
    (varobjectp(obj) && (Record_type(obj) == Rectype_Svector))
#endif
%% export_def(simple_vector_p(obj));

# Test for general-vector=(vector t)
#ifdef TYPECODES
  #define general_vector_p(obj)  \
    ((typecode(obj) & ~bit(notsimple_bit_t)) == svector_type)
#else
  # cases: Rectype_Svector, Rectype_vector
  #define general_vector_p(obj)  \
    (varobjectp(obj) \
     && ((Record_type(obj) & ~(Rectype_Svector ^ Rectype_vector)) == (Rectype_Svector & Rectype_vector)) \
    )
#endif
%% export_def(general_vector_p(obj));

# Test for simple-string
#ifdef TYPECODES
  #define simple_string_p(obj)  \
    (typecode(obj) == sstring_type)
#else
  # cases: Rectype_[Imm_]S[8|16|32]string, Rectype_reallocstring
  #define simple_string_p(obj)  \
    (varobjectp(obj) && ((uintB)(Record_type(obj) - Rectype_S8string) \
                         <= Rectype_reallocstring - Rectype_S8string))
#endif
%% export_def(simple_string_p(obj));

# Test for string
#ifdef TYPECODES
  #define stringp(obj)  \
    ((typecode(obj) & ~bit(notsimple_bit_t)) == sstring_type)
#else
  # cases: Rectype_[Imm_]S[8|16|32]string, Rectype_reallocstring, Rectype_string
  #define stringp(obj)  \
    (varobjectp(obj) && ((uintB)(Record_type(obj) - Rectype_S8string) \
                         <= Rectype_string - Rectype_S8string))
#endif
%% export_def(stringp(obj));

/* test for (VECTOR NIL) */
#ifdef TYPECODES
  #define nil_vector_p(obj)  \
    (typecode(obj) == vector_type \
     && (Iarray_flags(obj) & arrayflags_atype_mask) == Atype_NIL \
    )
#else
  # cases: Rectype_Svector, Rectype_vector
  #define nil_vector_p(obj)  \
    (varobjectp(obj) \
     && (Record_type(obj) == Rectype_vector \
         && (Iarray_flags(obj) & arrayflags_atype_mask) == Atype_NIL \
    )   )
#endif

# Test for simple-bit[n]-vector
#ifdef TYPECODES
  #define simple_bit_vector_p(atype,obj)  \
    (typecode(obj) == Array_type_simple_bit_vector(atype))
#else
  # cases: Rectype_Sb[2^n]vector
  #define simple_bit_vector_p(atype,obj)  \
    (varobjectp(obj) && (Record_type(obj) == Rectype_Sbvector+(atype)))
#endif
%% export_def(simple_bit_vector_p(atype,obj));

# Test for bit[n]-vector
#ifdef TYPECODES
  #define bit_vector_p(atype,obj)  \
    ((typecode(obj) & ~bit(notsimple_bit_t)) == Array_type_simple_bit_vector(atype))
#else
  # cases: Rectype_Sb[2^n]vector, Rectype_b[2^n]vector
  #define bit_vector_p(atype,obj)  \
    (varobjectp(obj) \
     && ((Record_type(obj) & ~(Rectype_Sbvector ^ Rectype_bvector)) == (Rectype_Sbvector & Rectype_bvector) + (atype)) \
    )
#endif
%% export_def(bit_vector_p(atype,obj));

# Test for Array (general)
#ifdef TYPECODES
  #define arrayp(obj)  \
    ((tint)(typecode(obj) - mdarray_type) <= (tint)(vector_type - mdarray_type))
#else
  # cases: Rectype_Sbvector, Rectype_Sb[2|4|8|16|32]vector, Rectype_Svector, Rectype_[Imm_]S[8|16|32]string,
  #        Rectype_bvector, Rectype_b[2|4|8|16|32]vector, Rectype_vector, Rectype_reallocstring, Rectype_string,
  #        Rectype_mdarray
  #define arrayp(obj)  \
    (varobjectp(obj) && ((uintB)(Record_type(obj) - Rectype_vector) \
                         <= Rectype_mdarray - Rectype_vector))
#endif
%% export_def(arrayp(obj));

# Test for Array, that isn't a Vector (type byte %100)
#ifdef TYPECODES
  #define mdarrayp(obj)  \
    (typecode(obj) == mdarray_type)
#else
  # cases: Rectype_mdarray
  #define mdarrayp(obj)  \
    (varobjectp(obj) && (Record_type(obj) == Rectype_mdarray))
#endif

#ifdef TYPECODES
  # Test for Closure/Structure/Stream/Instance/OtherRecord/LongRecord
    #define if_recordp(obj,statement1,statement2)  \
      switch (typecode(obj)) {          \
        case_record: statement1; break; \
        default: statement2; break;     \
      }
#else
  # Test for Lrecord/Srecord/Xrecord
    #define if_recordp(obj,statement1,statement2)  \
      if (orecordp(obj))                                                     \
        switch (Record_type(obj)) {                                          \
          case Rectype_Sbvector:                                             \
          case Rectype_S8string: case Rectype_Imm_S8string:                  \
          case Rectype_S16string: case Rectype_Imm_S16string:                \
          case Rectype_S32string: case Rectype_Imm_S32string:                \
          case Rectype_Svector:                                              \
          case Rectype_mdarray:                                              \
          case Rectype_bvector: case Rectype_string: case Rectype_vector:    \
          case Rectype_reallocstring:                                        \
          case Rectype_Bignum: case Rectype_Lfloat:                          \
            goto not_record;                                                 \
          default: { statement1 } break;                                     \
        }                                                                    \
      else                                                                   \
        not_record: { statement2 }
#endif

# Test for Closure
#ifdef TYPECODES
  #define closurep(obj)  (typecode(obj)==closure_type)
#else
  #define closurep(obj)  \
    (varobjectp(obj) && (Record_type(obj) == Rectype_Closure))
#endif

# Test for compiled Closure
# The second component of a closure is either a list
# (the Lambdabody for interpreted Closures)
# or a Simple-Bit-Vector (the code vector for compiled Closures).
#define cclosurep(obj)  \
  (closurep(obj)        \
   && simple_bit_vector_p(Atype_8Bit,TheClosure(obj)->clos_codevec))

# Test for a function with a code vector produced by %GENERIC-FUNCTION-LAMBDA.
#define genericlambda_function_p(obj)  \
  (cclosurep(obj) \
   && (TheCodevec(TheClosure(obj)->clos_codevec)->ccv_flags & bit(4)))

# Test for CLOS-Instance
#ifdef TYPECODES
  #define instancep(obj)  \
    (typecode(obj)==instance_type                                \
     || (typecode(obj)==closure_type && Closure_instancep(obj)))
#else
  #define instancep(obj)  \
    (varobjectp(obj)                                                  \
     && (Record_type(obj) == Rectype_Instance                         \
         || (Record_type(obj) == Rectype_Closure && Closure_instancep(obj))))
#endif
# Test for non-funcallable CLOS-Instance
#ifdef TYPECODES
  #define regular_instance_p(obj)  (typecode(obj)==instance_type)
#else
  #define regular_instance_p(obj)  \
    (varobjectp(obj) && (Record_type(obj) == Rectype_Instance))
#endif
# Test for funcallable CLOS-Instance
#define funcallable_instance_p(obj)  \
  (closurep(obj) && Closure_instancep(obj))
%% export_def(instancep(obj));

# Test for CLOS-class or forward-reference.
# Our CLOS implements all classes as instances of a
# (not necessarily direct) subclass of <class>.
#define if_potential_class_p(obj,statement1,statement2)  \
  if (instancep(obj)) {                                               \
    {                                                                 \
      var object obj_forwarded = obj;                                 \
      instance_un_realloc(obj_forwarded);                             \
      /*instance_update(obj,obj_forwarded); - not needed since we don't access a slot */ \
     {var object cv = TheInstance(obj_forwarded)->inst_class_version; \
      /* Treat the most frequent cases first, for speed. */           \
      if (eq(cv,O(class_version_standard_class))) # direct instance of STANDARD-CLASS? \
        goto obj##_classp_yes;                                        \
      if (eq(cv,O(class_version_structure_class))) # direct instance of STRUCTURE-CLASS? \
        goto obj##_classp_yes;                                        \
      if (eq(cv,O(class_version_built_in_class))) # direct instance of BUILT-IN-CLASS? \
        goto obj##_classp_yes;                                        \
      /* Now a slow, but general instanceof test. */                  \
      {var object objclas = TheClassVersion(cv)->cv_newest_class;     \
       if (eq(gethash(O(class_potential_class),TheClass(objclas)->all_superclasses,false),nullobj)) \
         goto obj##_classp_no;                                        \
    }}}                                                               \
   obj##_classp_yes: statement1;                                      \
  } else {                                                            \
   obj##_classp_no: statement2;                                       \
  }

# Test for CLOS-class.
# Our CLOS implements all classes as instances of a
# (not necessarily direct) subclass of <defined-class>.
#define if_defined_class_p(obj,statement1,statement2)  \
  if (instancep(obj)) {                                               \
    {                                                                 \
      var object obj_forwarded = obj;                                 \
      instance_un_realloc(obj_forwarded);                             \
      /*instance_update(obj,obj_forwarded); - not needed since we don't access a slot */ \
     {var object cv = TheInstance(obj_forwarded)->inst_class_version; \
      /* Treat the most frequent cases first, for speed. */           \
      if (eq(cv,O(class_version_standard_class))) # direct instance of STANDARD-CLASS? \
        goto obj##_classp_yes;                                        \
      if (eq(cv,O(class_version_structure_class))) # direct instance of STRUCTURE-CLASS? \
        goto obj##_classp_yes;                                        \
      if (eq(cv,O(class_version_built_in_class))) # direct instance of BUILT-IN-CLASS? \
        goto obj##_classp_yes;                                        \
      /* Now a slow, but general instanceof test. */                  \
      {var object objclas = TheClassVersion(cv)->cv_newest_class;     \
       if (eq(gethash(O(class_defined_class),TheClass(objclas)->all_superclasses,false),nullobj)) \
         goto obj##_classp_no;                                        \
    }}}                                                               \
   obj##_classp_yes: statement1;                                      \
  } else {                                                            \
   obj##_classp_no: statement2;                                       \
  }

# Test for Other-Record
# This is not really a type test (because there is no well-defined type
# Other-Record). It's just a precondition for calling Record_type(obj).
#ifdef TYPECODES
  #define orecordp(obj)  (typecode(obj)==orecord_type)
#else
  #define orecordp(obj)  varobjectp(obj)
#endif
%% export_def(orecordp(obj));

# Test for Long-Record
# This is not really a type test (because there is no well-defined type
# Long-Record). It's just a precondition for calling Record_type(obj).
#ifdef TYPECODES
  #define lrecordp(obj)  (typecode(obj)==lrecord_type)
#else
  #define lrecordp(obj)  varobjectp(obj)
#endif

# Test for Structure
#ifdef case_structure
  #define structurep(obj)  (typecode(obj)==structure_type)
#else
  #define structurep(obj)  \
    (orecordp(obj) && (Record_type(obj) == Rectype_Structure))
#endif
%% export_def(structurep(obj));

# Test for Builtin-Stream
#ifdef case_stream
  #define builtin_stream_p(obj)  (typecode(obj)==stream_type)
#else
  #define builtin_stream_p(obj)  \
    (orecordp(obj) && (Record_type(obj) == Rectype_Stream))
#endif
%% export_def(builtin_stream_p(obj));

# Test for Stream
#define streamp(obj)  \
  (builtin_stream_p(obj) || instanceof(obj,O(class_fundamental_stream)))

# Test for Package
#define packagep(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Package))
%% #if notused
%% export_def(packagep(obj));
%% #endif

# Test for Hash-Table
#define hash_table_p(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Hashtable))

# Test for Readtable
#define readtablep(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Readtable))

# Test for Pathname
#define pathnamep(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Pathname))

# Test for Logical Pathname
#ifdef LOGICAL_PATHNAMES
  #define logpathnamep(obj)  \
    (orecordp(obj) && (Record_type(obj) == Rectype_Logpathname))
#else
  #define logpathnamep(obj)  false
#endif

# Test for Extended Pathname (i.e., Pathname or Logical Pathname)
# define xpathnamep(obj)  (pathnamep(obj) || logpathnamep(obj))
#ifdef LOGICAL_PATHNAMES
  #define xpathnamep(obj)  \
    (orecordp(obj)                                    \
     && ((Record_type(obj) == Rectype_Pathname)       \
         || (Record_type(obj) == Rectype_Logpathname)))
#else
  #define xpathnamep(obj)  pathnamep(obj)
#endif

# Test for Random-State
#define random_state_p(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Random_State))

# Test for Byte
#define bytep(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Byte))

# Test for Fsubr
#define fsubrp(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Fsubr))

# Test for Loadtimeeval
#define loadtimeevalp(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Loadtimeeval))

# Test for Symbolmacro
#define symbolmacrop(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Symbolmacro))

# Test for GlobalSymbolmacro
#define globalsymbolmacrop(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_GlobalSymbolmacro))

# Test for Macro
#define macrop(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Macro))

# Test for FunctionMacro
#define functionmacrop(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_FunctionMacro))

# Test for BigReadLabel
#define big_read_label_p(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_BigReadLabel))

# Test for Encoding
#define encodingp(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Encoding))

# Test for Fpointer
#define fpointerp(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Fpointer))
%% #ifdef FOREIGN
%%   export_def(fpointerp(obj));
%% #endif

# Test for Faddress
#define faddressp(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Faddress))

# Test for Fvariable
#define fvariablep(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Fvariable))

# Test for Ffunction
#ifdef DYNAMIC_FFI
  #define ffunctionp(obj)  \
    (orecordp(obj) && (Record_type(obj) == Rectype_Ffunction))
#else
  #define ffunctionp(obj)  ((void)(obj), 0)
#endif

# Test for Function
#define functionp(obj) (subrp(obj) || closurep(obj) || ffunctionp(obj))

# Test for Weakpointer
#define weakpointerp(obj)  \
  (orecordp(obj) && (Record_type(obj) == Rectype_Weakpointer))

# test for socket-server and for socket-stream
#ifdef SOCKET_STREAMS
  #define socket_server_p(obj)  \
    (orecordp(obj) && (Record_type(obj) == Rectype_Socket_Server))
  #define socket_stream_p(obj)  \
    (builtin_stream_p(obj) && (TheStream(obj)->strmtype==strmtype_socket))
#endif

#ifdef YET_ANOTHER_RECORD
  # Test for Yetanother
  #define yetanotherp(obj)  \
    (orecordp(obj) && (Record_type(obj) == Rectype_Yetanother))
#endif

# Test for Character
#ifdef TYPECODES
  #define charp(obj)  (typecode(obj)==char_type)
#else
  #define charp(obj)  ((as_oint(obj) & ((7 << imm_type_shift) | immediate_bias)) == char_type)
#endif
%% export_def(charp(obj));

#if (base_char_code_limit < char_code_limit)
# Test for base character
  #define base_char_p(obj)  \
    ((as_oint(obj) & ~((oint)(bit(base_char_int_len)-1)<<oint_data_shift)) == type_zero_oint(char_type))
#endif

# Test for SUBR (compiled functional object)
#ifdef TYPECODES
  #define subrp(obj)  (typecode(obj)==subr_type)
#else
  #ifdef STANDARD_HEAPCODES
    #define subrp(obj)  ((as_oint(obj) & 3) == subr_bias)
    #define immsubrp(obj)  subrp(obj)
  #endif
  #ifdef LINUX_NOEXEC_HEAPCODES
    #define subrp(obj)  (orecordp(obj) && (Record_type(obj) == Rectype_Subr))
    #define immsubrp(obj)  false
    #ifdef DEBUG_GCSAFETY
      # This is used by pgci_pointable, so it cannot use pgci_pointable itself.
      static inline bool nonimmsubrp (object obj) {
        return (varobjectp(obj)
                && (inside_gc /* Avoid doing memory accesses during GC. */
                    || (varobject_type((Record)(cgci_pointable(obj)-varobject_bias)) == Rectype_Subr)));
      }
    #endif
  #endif
#endif
%% #ifndef TYPECODES
%%   #ifdef LINUX_NOEXEC_HEAPCODES
%%     #ifdef DEBUG_GCSAFETY
%%       printf2("static inline bool nonimmsubrp (object obj) { return (varobjectp(obj) && (varobject_type((Record)(cgci_pointable(obj)-%d)) == %d)); }\n",varobject_bias,Rectype_Subr);
%%     #endif
%%   #endif
%% #endif

# Test for pointer into the STACK (usually at a frame)
#ifdef TYPECODES
  #define framepointerp(obj)  (typecode(obj)==system_type) # other cases??
#else
  #define framepointerp(obj)  ((as_oint(obj) & 3) == machine_bias) # other cases??
#endif

#ifndef TYPECODES

  # Test for Machine-Pointer
  #ifdef STANDARD_HEAPCODES
    #define machinep(obj)  ((as_oint(obj) & 3) == machine_bias)
  #endif
  #ifdef LINUX_NOEXEC_HEAPCODES
    #define machinep(obj)  \
      ((as_oint(obj) & 3) == machine_bias            \
       && (as_oint(obj) & 0xE0000000) != 0xC0000000)
  #endif

  # Test for Small-Read-Label
  #define small_read_label_p(obj)  ((as_oint(obj) & ((7 << imm_type_shift) | immediate_bias)) == small_read_label_type)

  # Test for System-Pointer
  #define systemp(obj)  ((as_oint(obj) & ((7 << imm_type_shift) | immediate_bias)) == system_type)

#endif

# Test for real number
#ifdef TYPECODES
  #define if_realp(obj,statement1,statement2)                           \
    do {                                                                \
      var object obj_from_if_realp = (obj);                             \
      var tint type_from_if_realp = typecode(obj_from_if_realp);        \
      if ( (type_from_if_realp & bit(number_bit_t))                     \
           && !(type_from_if_realp==complex_type) )                     \
        { statement1 } else { statement2 }                              \
    } while(0)
#else
  #define if_realp(obj,statement1,statement2)                           \
    do { if (((as_oint(obj) & ((4 << imm_type_shift) | immediate_bias)) \
              == fixnum_type)                                           \
             || (varobjectp(obj)                                        \
                 && ((uintB)(Record_type(obj)-Rectype_Bignum) <=        \
                     Rectype_Ratio-Rectype_Bignum)))                    \
           { statement1 } else { statement2 }                           \
    } while(0)
#endif

# Test for rational number
#ifdef TYPECODES
  #define if_rationalp(obj,statement1,statement2)                        \
    do {                                                                 \
      var object obj_from_if_rationalp = (obj);                          \
      var tint type_from_if_rationalp = typecode(obj_from_if_rationalp); \
      if ((type_from_if_rationalp != complex_type)                       \
           && ((type_from_if_rationalp &                                 \
                ~((fixnum_type|bignum_type|ratio_type|bit(sign_bit_t))   \
                  & ~(fixnum_type&bignum_type&ratio_type)))              \
               == (fixnum_type&bignum_type&ratio_type)))                 \
        { statement1 } else { statement2 }                               \
    } while(0)
#else
  #define if_rationalp(obj,statement1,statement2)                        \
    do { if (((as_oint(obj) & ((6 << imm_type_shift) | immediate_bias))  \
              == fixnum_type)                                            \
             || (varobjectp(obj)                                         \
                 && ((Record_type(obj) == Rectype_Bignum)                \
                     || (Record_type(obj) == Rectype_Ratio))))           \
           { statement1 } else { statement2 }                            \
    } while(0)

#endif

# Test for Integer
#ifdef TYPECODES
  #define integerp(obj)  \
    ((typecode(obj) &                                                        \
      ~((fixnum_type|bignum_type|bit(sign_bit_t)) & ~(fixnum_type&bignum_type)) \
     ) == (fixnum_type&bignum_type))
#else
  #define integerp(obj)  \
   (((as_oint(obj) & ((6 << imm_type_shift) | immediate_bias)) == fixnum_type) \
    || (varobjectp(obj) && (Record_type(obj) == Rectype_Bignum)))
#endif
%% export_def(integerp(obj));

# Test for Fixnum
#ifdef TYPECODES
  #define fixnump(obj)  ((typecode(obj) & ~bit(sign_bit_t)) == fixnum_type)
#else
  #define fixnump(obj)  ((as_oint(obj) & ((6 << imm_type_shift) | immediate_bias)) == fixnum_type)
#endif
%% export_def(fixnump(obj));

# Test for Fixnum >=0
#ifdef TYPECODES
  #define posfixnump(obj)  (typecode(obj) == fixnum_type)
#else
  #define posfixnump(obj)  ((as_oint(obj) & ((7 << imm_type_shift) | immediate_bias)) == fixnum_type)
#endif
%% export_def(posfixnump(obj));

# Test for Bignum
#ifdef TYPECODES
  #define bignump(obj)  ((typecode(obj) & ~bit(sign_bit_t)) == bignum_type)
#else
  #define bignump(obj)  \
    (varobjectp(obj) && (Record_type(obj) == Rectype_Bignum))
#endif
%% export_def(bignump(obj));

# Test for Bignum >=0
#ifdef TYPECODES
  #define posbignump(obj)  (typecode(obj) == bignum_type)
#else
  #define posbignump(obj)  \
    (varobjectp(obj)                         \
     && (Record_type(obj) == Rectype_Bignum) \
     && ((Record_flags(obj) & bit(7)) == 0))
#endif
%% export_def(posbignump(obj));

# Test for Ratio
#ifdef TYPECODES
  #define ratiop(obj)  ((typecode(obj) & ~bit(sign_bit_t)) == ratio_type)
#else
  #define ratiop(obj)  (varobjectp(obj) && (Record_type(obj) == Rectype_Ratio))
#endif
%% #if notused
%%   export_def(ratiop(obj));
%% #endif

# Test for Float
#ifdef TYPECODES
  #define floatp(obj)  \
    ((typecode(obj) &  \
     ~((sfloat_type|ffloat_type|dfloat_type|lfloat_type|bit(sign_bit_t)) & ~(sfloat_type&ffloat_type&dfloat_type&lfloat_type)) \
     ) == (sfloat_type&ffloat_type&dfloat_type&lfloat_type))
#else
  #define floatp(obj)  \
    (((as_oint(obj) & ((6 << imm_type_shift) | immediate_bias)) == sfloat_type) \
     || (varobjectp(obj)                    \
         && ((uintB)(Record_type(obj)-Rectype_Lfloat) <= Rectype_Ffloat-Rectype_Lfloat)))
#endif
%% #if notused
%%   export_def(floatp(obj));
%% #endif

# Test for Short-Float
#ifdef TYPECODES
  #define short_float_p(obj)  ((typecode(obj) & ~bit(sign_bit_t)) == sfloat_type)
#else
  #define short_float_p(obj)  ((as_oint(obj) & ((6 << imm_type_shift) | immediate_bias)) == sfloat_type)
#endif
%% #if notused
%%   export_def(short_float_p(obj));
%% #endif

# Test for Single-Float
#ifdef TYPECODES
  #define single_float_p(obj)  ((typecode(obj) & ~bit(sign_bit_t)) == ffloat_type)
#else
  #define single_float_p(obj)  (varobjectp(obj) && (Record_type(obj) == Rectype_Ffloat))
#endif
%% export_def(single_float_p(obj));

# Test for Double-Float
#ifdef TYPECODES
  #define double_float_p(obj)  ((typecode(obj) & ~bit(sign_bit_t)) == dfloat_type)
#else
  #define double_float_p(obj)  (varobjectp(obj) && (Record_type(obj) == Rectype_Dfloat))
#endif
%% export_def(double_float_p(obj));

# Test for Long-Float
#ifdef TYPECODES
  #define long_float_p(obj)  ((typecode(obj) & ~bit(sign_bit_t)) == lfloat_type)
#else
  #define long_float_p(obj)  (varobjectp(obj) && (Record_type(obj) == Rectype_Lfloat))
#endif
%% #if notused
%%   export_def(long_float_p(obj));
%% #endif

# Test for Complex
#ifdef TYPECODES
  #define complexp(obj)  (typecode(obj) == complex_type)
#else
  #define complexp(obj)  (varobjectp(obj) && (Record_type(obj) == Rectype_Complex))
#endif
%% #if notused
%%   export_def(complexp(obj));
%% #endif

# Test if a real number is >=0:
#ifdef TYPECODES
  # define positivep(obj)  ((as_oint(obj) & wbit(sign_bit_o)) == 0)
  #define positivep(obj)  (!wbit_test(as_oint(obj),sign_bit_o))
  #ifdef WIDE_STRUCT
    #undef positivep
    #define positivep(obj)  ((typecode(obj) & bit(sign_bit_t)) == 0)
  #endif
#else
  #define positivep(obj)  \
    (number_immediatep(obj)                                        \
     ? /* fixnum, sfloat */ (as_oint(obj) & wbit(sign_bit_o)) == 0 \
     : /* bignum, [fdl]float */ (Record_flags(obj) & bit(7)) == 0)
#endif
%% export_def(positivep(obj));


# switch with typcodes:
# example:
#   switch (typecode(obj)) {
#     case_symbol: ....
#     case_orecord:
#       switch (Record_type(obj)) {
#         case_Rectype_Symbol_above;
#         ...
#       }
#   }

#ifdef case_structure
  #define case_Rectype_Structure_above
#else
  #define case_Rectype_Structure_above  \
    case Rectype_Structure: goto case_structure;
#endif

#ifdef case_stream
  #define case_Rectype_Stream_above
#else
  #define case_Rectype_Stream_above  \
    case Rectype_Stream: goto case_stream;
#endif

#ifdef TYPECODES
  #define case_Rectype_Closure_above
  #define case_Rectype_Instance_above
  #define case_Rectype_Sbvector_above
  #define case_Rectype_Sb2vector_above
  #define case_Rectype_Sb4vector_above
  #define case_Rectype_Sb8vector_above
  #define case_Rectype_Sb16vector_above
  #define case_Rectype_Sb32vector_above
  #define case_Rectype_Sstring_above
  #define case_Rectype_Svector_above
  #define case_Rectype_mdarray_above
  #define case_Rectype_obvector_above
  #define case_Rectype_ob2vector_above
  #define case_Rectype_ob4vector_above
  #define case_Rectype_ob8vector_above
  #define case_Rectype_ob16vector_above
  #define case_Rectype_ob32vector_above
  #define case_Rectype_ostring_above
  #define case_Rectype_ovector_above
  #define case_Rectype_Bignum_above
  #define case_Rectype_Lfloat_above
  #define case_Rectype_Dfloat_above
  #define case_Rectype_Ffloat_above
  #define case_Rectype_Ratio_above
  #define case_Rectype_Complex_above
  #define case_Rectype_Symbol_above
  # Composite cases:
  #define case_Rectype_string_above
  #define case_Rectype_bvector_above
  #define case_Rectype_b2vector_above
  #define case_Rectype_b4vector_above
  #define case_Rectype_b8vector_above
  #define case_Rectype_b16vector_above
  #define case_Rectype_b32vector_above
  #define case_Rectype_vector_above
  #define case_Rectype_array_above
  #define case_Rectype_number_above
  #define case_Rectype_float_above
  #define case_Rectype_integer_above
#else
  #define case_Rectype_Closure_above  \
    case Rectype_Closure: goto case_closure;
  #define case_Rectype_Instance_above  \
    case Rectype_Instance: goto case_instance;
  #define case_Rectype_Sbvector_above  \
    case Rectype_Sbvector: goto case_sbvector;
  #define case_Rectype_Sb2vector_above  \
    case Rectype_Sb2vector: goto case_sb2vector;
  #define case_Rectype_Sb4vector_above  \
    case Rectype_Sb4vector: goto case_sb4vector;
  #define case_Rectype_Sb8vector_above  \
    case Rectype_Sb8vector: goto case_sb8vector;
  #define case_Rectype_Sb16vector_above  \
    case Rectype_Sb16vector: goto case_sb16vector;
  #define case_Rectype_Sb32vector_above  \
    case Rectype_Sb32vector: goto case_sb32vector;
  #define case_Rectype_Sstring_above  \
    case Rectype_S8string: case Rectype_Imm_S8string: case Rectype_S16string: case Rectype_Imm_S16string: case Rectype_S32string: case Rectype_Imm_S32string: case Rectype_reallocstring: goto case_sstring;
  #define case_Rectype_Svector_above  \
    case Rectype_Svector: goto case_svector;
  #define case_Rectype_mdarray_above  \
    case Rectype_mdarray: goto case_mdarray;
  #define case_Rectype_obvector_above  \
    case Rectype_bvector: goto case_obvector;
  #define case_Rectype_ob2vector_above  \
    case Rectype_b2vector: goto case_ob2vector;
  #define case_Rectype_ob4vector_above  \
    case Rectype_b4vector: goto case_ob4vector;
  #define case_Rectype_ob8vector_above  \
    case Rectype_b8vector: goto case_ob8vector;
  #define case_Rectype_ob16vector_above  \
    case Rectype_b16vector: goto case_ob16vector;
  #define case_Rectype_ob32vector_above  \
    case Rectype_b32vector: goto case_ob32vector;
  #define case_Rectype_ostring_above  \
    case Rectype_string: goto case_ostring;
  #define case_Rectype_ovector_above  \
    case Rectype_vector: goto case_ovector;
  #define case_Rectype_Bignum_above  \
    case Rectype_Bignum: goto case_bignum;
  #define case_Rectype_Lfloat_above  \
    case Rectype_Lfloat: goto case_lfloat;
  #define case_Rectype_Dfloat_above  \
    case Rectype_Dfloat: goto case_dfloat;
  #define case_Rectype_Ffloat_above  \
    case Rectype_Ffloat: goto case_ffloat;
  #define case_Rectype_Ratio_above  \
    case Rectype_Ratio: goto case_ratio;
  #define case_Rectype_Complex_above  \
    case Rectype_Complex: goto case_complex;
  #define case_Rectype_Symbol_above  \
    case Rectype_Symbol: goto case_symbol;
  # Composite cases:
  #define case_Rectype_string_above  \
    case Rectype_S8string: case Rectype_Imm_S8string: case Rectype_S16string: case Rectype_Imm_S16string: case Rectype_S32string: case Rectype_Imm_S32string: case Rectype_reallocstring: case Rectype_string: goto case_string;
  #define case_Rectype_bvector_above  \
    case Rectype_Sbvector: case Rectype_bvector: goto case_bvector;
  #define case_Rectype_b2vector_above  \
    case Rectype_Sb2vector: case Rectype_b2vector: goto case_b2vector;
  #define case_Rectype_b4vector_above  \
    case Rectype_Sb4vector: case Rectype_b4vector: goto case_b4vector;
  #define case_Rectype_b8vector_above  \
    case Rectype_Sb8vector: case Rectype_b8vector: goto case_b8vector;
  #define case_Rectype_b16vector_above  \
    case Rectype_Sb16vector: case Rectype_b16vector: goto case_b16vector;
  #define case_Rectype_b32vector_above  \
    case Rectype_Sb32vector: case Rectype_b32vector: goto case_b32vector;
  #define case_Rectype_vector_above  \
    case Rectype_Svector: case Rectype_vector: goto case_vector;
  #define case_Rectype_array_above                            \
    case Rectype_S8string: case Rectype_Imm_S8string:         \
    case Rectype_S16string: case Rectype_Imm_S16string:       \
    case Rectype_S32string: case Rectype_Imm_S32string:       \
    case Rectype_reallocstring: case Rectype_string:          \
    case Rectype_Sbvector: case Rectype_bvector:              \
    case Rectype_Sb2vector: case Rectype_b2vector:            \
    case Rectype_Sb4vector: case Rectype_b4vector:            \
    case Rectype_Sb8vector: case Rectype_b8vector:            \
    case Rectype_Sb16vector: case Rectype_b16vector:          \
    case Rectype_Sb32vector: case Rectype_b32vector:          \
    case Rectype_Svector: case Rectype_vector:                \
    case Rectype_mdarray:                                     \
      goto case_array;
  #define case_Rectype_number_above  /* don't forget immediate_number_p */ \
    case Rectype_Complex: case Rectype_Ratio:                      \
    case Rectype_Ffloat: case Rectype_Dfloat: case Rectype_Lfloat: \
    case Rectype_Bignum:                                           \
      goto case_number;
  #define case_Rectype_float_above  /* don't forget short_float_p */ \
    case Rectype_Ffloat: case Rectype_Dfloat: case Rectype_Lfloat: \
      goto case_float;
  #define case_Rectype_integer_above  /* don't forget fixnump */ \
    case Rectype_Bignum: goto case_integer;
#endif

#if defined(TYPECODES) || defined(STANDARD_HEAPCODES)
  #define case_Rectype_Subr_above
#else # LINUX_NOEXEC_HEAPCODES
  #define case_Rectype_Subr_above  \
    case Rectype_Subr: goto case_subr;
#endif


# ################# Declarations for the arithmetics ######################## #


# Type hierachy :
# Number (N) =
#    Real (R) =
#       Float (F) =
#          Short float (SF)
#          Single float (FF)
#          Double float (DF)
#          Long float (LF)
#       Rational (RA) =
#          Integer (I) =
#             Fixnum (FN)
#             Bignum (BN)
#          Ratio (RT)
#    Complex (C)


# Type field:
# Bytes for testing whether it's that type (Bit set, is yes).
# _bit_t to test in the type byte (tint)
# _bit_o to test in the object (oint)

#ifndef NUMBER_BITS_INVERTED
  #define number_wbit_test  wbit_test
#else
  #define number_wbit_test  !wbit_test
#endif

#ifdef TYPECODES

# see above:
# #define number_bit_t     4  # set only for numbers
# #define number_bit_o     (number_bit_t+oint_type_shift)    # set only for numbers

# float_bit:
# in a number : Bit set, if it's a Float.
#                Bit unset, if it's a rational or complex number.
# (For NUMBER_BITS_INVERTED it's exactly the other way around.)
# #define float_bit_t      1
# #define float_bit_o      (float_bit_t+oint_type_shift)

# float1_bit:
# In a floating-point: discriminates further:
#ifndef NUMBER_BITS_INVERTED
# Float-Bit   1 2
#             0 0    Short Float (SF)
#             0 1    Single Float (FF)
#             1 0    Double Float (DF)
#             1 1    Long Float (LF)
#else
# Float-Bit   1 2
#             0 0    Long Float (LF)
#             0 1    Double Float (DF)
#             1 0    Single Float (FF)
#             1 1    Short Float (SF)
#endif
# #define float1_bit_t     3
# #define float1_bit_o     (float1_bit_t+oint_type_shift)
# #define float2_bit_t     2
# #define float2_bit_o     (float2_bit_t+oint_type_shift)

# ratio_bit:
# For rational numbers: Bit set , if it's a real fraction.
#                       Bit unset, if it's an Integer.
# (For NUMBER_BITS_INVERTED it's exactly the other way around..)
# #define ratio_bit_t      3
# #define ratio_bit_o      (ratio_bit_t+oint_type_shift)

# bignum_bit:
# For Integers:     Bit set, if it's a Bignum.
#                   Bit unset, if it's a Fixnum.
# (For NUMBER_BITS_INVERTED it's exactly the other way around..)
# #define bignum_bit_t     2
# #define bignum_bit_o     (bignum_bit_t+oint_type_shift)

# vorz_bit: (sign bit)
# For Reals:
# returns the sign of the number.
# Bit set, if number < 0,
# Bit unset, if number >=0.
  #define vorz_bit_t       sign_bit_t
                           # should be = 0, so the sign-extend
                           # is easier for Fixnums.
  #define vorz_bit_o       (vorz_bit_t+oint_type_shift)

#endif

# return the sign of a real number (0 if >=0, -1 if <0)
#ifdef TYPECODES
  #if (vorz_bit_o<32) && !defined(WIDE_STRUCT)
    #define R_sign(obj)  ((signean)sign_of_sint32( (sint32)((uint32)as_oint(obj) << (31-vorz_bit_o)) ))
  #else
    # define R_sign(obj)  ((signean)sign_of_sint32( (sint32)(uint32)(as_oint(obj) >> (vorz_bit_o-31)) ))
    #define R_sign(obj)  ((signean)sign_of_sint32( (sint32)((uint32)typecode(obj) << (31-vorz_bit_t)) ))
  #endif
#else
  #define R_sign(obj)  ((signean)sign_of_sint32(_R_sign(obj)))
  #define _R_sign(obj)  \
    (number_immediatep(obj)                                         \
     ? /* fixnum, sfloat */ (sint32)as_oint(obj) << (31-sign_bit_o) \
     : /* [fdl]float */ (sint32)(sintB)Record_flags(obj))
#endif

# Gives the sign of a Fixnum/Bignum/Ratio/
# Short-/Single-/Double-/Long-Float.
#ifdef TYPECODES
  #define FN_sign(obj)  R_sign(obj)
  #define BN_sign(obj)  R_sign(obj)
  #define RT_sign(obj)  R_sign(obj)
  #define SF_sign(obj)  R_sign(obj)
  #define FF_sign(obj)  R_sign(obj)
  #define DF_sign(obj)  R_sign(obj)
  #define LF_sign(obj)  R_sign(obj)
#else
  #define FN_sign(obj)  \
    ((signean)sign_of_sint32((sint32)as_oint(obj) << (31-sign_bit_o)))
  #define BN_sign(obj)  \
    ((signean)sign_of_sint32((sint32)(sintB)Record_flags(obj)))
  #define RT_sign(obj)  \
    ((signean)sign_of_sint32((sint32)(sintB)Record_flags(obj)))
  #define SF_sign(obj)  \
    ((signean)sign_of_sint32((sint32)as_oint(obj) << (31-sign_bit_o)))
  #define FF_sign(obj)  \
    ((signean)sign_of_sint32((sint32)(sintB)Record_flags(obj)))
  #define DF_sign(obj)  \
    ((signean)sign_of_sint32((sint32)(sintB)Record_flags(obj)))
  #define LF_sign(obj)  \
    ((signean)sign_of_sint32((sint32)(sintB)Record_flags(obj)))
#endif

# Checks whether two real numbers have the same sign:
#ifdef TYPECODES
  #define same_sign_p(obj1,obj2)  \
    (wbit_test(as_oint(obj1)^as_oint(obj2),vorz_bit_o)==0)
#else
  #define same_sign_p(obj1,obj2)  \
    ((sint32)(_R_sign(obj1) ^ _R_sign(obj2)) >= 0)
#endif


# Type test macros:
# (Return /=0, if satisfied. Prefix 'm', if argument is in memory)

# Tests an objects whether it's a number: (see above)
# define numberp(obj)  ...

# Tests a number whether it's a Float.
#ifdef TYPECODES
  #ifndef NUMBER_BITS_INVERTED
    # define N_floatp(obj)  ( as_oint(obj) & wbit(float_bit_o) )
    #define N_floatp(obj)  (wbit_test(as_oint(obj),float_bit_o))
  #else
    #define N_floatp(obj)  (!wbit_test(as_oint(obj),float_bit_o))
  #endif
#else
  #define N_floatp(obj)  floatp(obj)
#endif

# Tests a number whether it's an Integer.
#ifdef TYPECODES
  #ifndef NUMBER_BITS_INVERTED
    #define N_integerp(obj)  (!( as_oint(obj) & (wbit(float_bit_o)|wbit(ratio_bit_o)) ))
  #else
    #define N_integerp(obj)  (!( (wbit(float_bit_o)|wbit(ratio_bit_o)) & ~as_oint(obj) ))
  #endif
#else
  #define N_integerp(obj)  integerp(obj)
#endif

# Tests a real number whether it's rational.
#ifdef TYPECODES
  #ifndef NUMBER_BITS_INVERTED
    # define R_rationalp(obj)  (!( as_oint(obj) & wbit(float_bit_o) ))
    #define R_rationalp(obj)  (!wbit_test(as_oint(obj),float_bit_o))
  #else
    #define R_rationalp(obj)  (wbit_test(as_oint(obj),float_bit_o))
  #endif
#else
  #define R_rationalp(obj)  (!floatp(obj))
#endif

# Tests a real number whether it's a Float.
#ifdef TYPECODES
  #ifndef NUMBER_BITS_INVERTED
    # define R_floatp(obj)  ( as_oint(obj) & wbit(float_bit_o) )
    #define R_floatp(obj)  (wbit_test(as_oint(obj),float_bit_o))
  #else
    #define R_floatp(obj)  (!wbit_test(as_oint(obj),float_bit_o))
  #endif
#else
  #define R_floatp(obj)  floatp(obj)
#endif

# Tests a real number whether it's <0.
#ifdef TYPECODES
  # define R_minusp(obj)  ( as_oint(obj) & wbit(vorz_bit_o) )
  #define R_minusp(obj)  (wbit_test(as_oint(obj),vorz_bit_o))
#else
  #define R_minusp(obj)  (!positivep(obj))
#endif
%% export_def(R_minusp(obj));

# Tests a rational number whether it's an Integer.
#ifdef TYPECODES
  #ifndef NUMBER_BITS_INVERTED
    # define RA_integerp(obj)  (!( as_oint(obj) & wbit(ratio_bit_o) ))
    #define RA_integerp(obj)  (!wbit_test(as_oint(obj),ratio_bit_o))
  #else
    #define RA_integerp(obj)  (wbit_test(as_oint(obj),ratio_bit_o))
  #endif
#else
  #define RA_integerp(obj)  (!ratiop(obj))
#endif

# Tests a rational number whether it's a fraction.
#ifdef TYPECODES
  #ifndef NUMBER_BITS_INVERTED
    # define RA_ratiop(obj)  ( as_oint(obj) & wbit(ratio_bit_o) )
    #define RA_ratiop(obj)  (wbit_test(as_oint(obj),ratio_bit_o))
  #else
    #define RA_ratiop(obj)  (!wbit_test(as_oint(obj),ratio_bit_o))
  #endif
#else
  #define RA_ratiop(obj)  ratiop(obj)
#endif

# Tests an Integer whether it's a Bignum.
#ifndef NUMBER_BITS_INVERTED
  # define I_bignump(obj)  ( as_oint(obj) & wbit(bignum_bit_o) )
  #define I_bignump(obj)  (wbit_test(as_oint(obj),bignum_bit_o))
#else
  #define I_bignump(obj)  (!wbit_test(as_oint(obj),bignum_bit_o))
#endif

# Tests an Integer whether it's a Fixnum.
#ifndef NUMBER_BITS_INVERTED
  # define I_fixnump(obj)  (!( as_oint(obj) & wbit(bignum_bit_o) ))
  #define I_fixnump(obj)  (!wbit_test(as_oint(obj),bignum_bit_o))
#else
  #define I_fixnump(obj)  (wbit_test(as_oint(obj),bignum_bit_o))
#endif

# Tests a Fixnum whether it is >=0.
#ifdef TYPECODES
  #define FN_positivep(obj)  positivep(obj)
#else
  #define FN_positivep(obj)  ((as_oint(obj) & wbit(sign_bit_o)) == 0)
#endif
%% export_def(FN_positivep(obj));

# Tests a Bignum whether it is >=0.
#ifdef TYPECODES
  #define BN_positivep(obj)  positivep(obj)
#else
  #define BN_positivep(obj)  ((Record_flags(obj) & bit(7)) == 0)
#endif
%% export_def(BN_positivep(obj));

# Tests a number whether it's a real number
#define N_realp(obj)  (!complexp(obj))

# Tests a number whether it's a complex number
#define N_complexp(obj)  complexp(obj)

# Tests two Integers whether both are Bignum.
#ifndef NUMBER_BITS_INVERTED
  #define I_I_bignums_p(obj1,obj2)  \
    (wbit_test(as_oint(obj1)&as_oint(obj2),bignum_bit_o))
#else
  #define I_I_bignums_p(obj1,obj2)  \
    (!wbit_test(as_oint(obj1)|as_oint(obj2),bignum_bit_o))
#endif

# Tests for an Integer from a given range.
# obj should be a variable
#define uint1_p(obj)  \
  ((as_oint(obj) & ~((oint)0x01 << oint_data_shift)) == as_oint(Fixnum_0))
#define uint2_p(obj)  \
  ((as_oint(obj) & ~((oint)0x03 << oint_data_shift)) == as_oint(Fixnum_0))
#define uint4_p(obj)  \
  ((as_oint(obj) & ~((oint)0x0F << oint_data_shift)) == as_oint(Fixnum_0))
#define uint8_p(obj)  \
  ((as_oint(obj) & ~((oint)0xFF << oint_data_shift)) == as_oint(Fixnum_0))
#define sint8_p(obj)  \
  (((as_oint(obj) ^ (FN_positivep(obj) ? 0 : as_oint(Fixnum_minus1)^as_oint(Fixnum_0))) & ~((oint)0x7F << oint_data_shift)) == as_oint(Fixnum_0))
#define uint16_p(obj)  \
  ((as_oint(obj) & ~((oint)0xFFFF << oint_data_shift)) == as_oint(Fixnum_0))
#define sint16_p(obj)  \
  (((as_oint(obj) ^ (FN_positivep(obj) ? 0 : as_oint(Fixnum_minus1)^as_oint(Fixnum_0))) & ~((oint)0x7FFF << oint_data_shift)) == as_oint(Fixnum_0))
#if (oint_data_len>=32)
  #define uint32_p(obj)  \
    ((as_oint(obj) & ~((oint)0xFFFFFFFFUL << oint_data_shift)) == as_oint(Fixnum_0))
#else
  #define uint32_p(obj)  \
    (posfixnump(obj) \
     || (posbignump(obj) \
         && (Bignum_length(obj) <= ceiling(33,intDsize)) \
         && ((Bignum_length(obj) < ceiling(33,intDsize)) \
             || (TheBignum(obj)->data[0] < (uintD)bit(32%intDsize)))))
#endif
#if (oint_data_len>=31)
  #define sint32_p(obj)  \
    (((as_oint(obj) ^ (FN_positivep(obj) ? 0 : as_oint(Fixnum_minus1)^as_oint(Fixnum_0))) & ~((oint)0x7FFFFFFFUL << oint_data_shift)) == as_oint(Fixnum_0))
#else
  #define sint32_p(obj)  \
    (fixnump(obj) \
     || (bignump(obj) \
         && (Bignum_length(obj) <= ceiling(32,intDsize)) \
         && ((Bignum_length(obj) < ceiling(32,intDsize)) \
             || ((TheBignum(obj)->data[0] ^ (BN_positivep(obj) ? (uintD)0 : ~(uintD)0)) < (uintD)bit(31%intDsize)))))
#endif
#define uint64_p(obj)  \
  (posfixnump(obj) \
   || (posbignump(obj) \
       && (Bignum_length(obj) <= ceiling(65,intDsize)) \
       && ((Bignum_length(obj) < ceiling(65,intDsize)) \
           || (TheBignum(obj)->data[0] < (uintD)bit(64%intDsize)))))
#define sint64_p(obj)  \
  (fixnump(obj) \
   || (bignump(obj) \
       && (Bignum_length(obj) <= ceiling(64,intDsize)) \
       && ((Bignum_length(obj) < ceiling(64,intDsize)) \
           || ((TheBignum(obj)->data[0] ^ (BN_positivep(obj) ? (uintD)0 : ~(uintD)0)) < (uintD)bit(63%intDsize)))))
#if (int_bitsize==16)
  #define uint_p  uint16_p
  #define sint_p  sint16_p
#else # (int_bitsize==32)
  #define uint_p  uint32_p
  #define sint_p  sint32_p
#endif
#if (long_bitsize==32)
  #define ulong_p  uint32_p
  #define slong_p  sint32_p
#else # (long_bitsize==64)
  #define ulong_p  uint64_p
  #define slong_p  sint64_p
#endif
%% export_def(uint8_p(obj));
%% export_def(sint8_p(obj));
%% export_def(uint16_p(obj));
%% export_def(sint16_p(obj));
%% export_def(uint32_p(obj));
%% export_def(sint32_p(obj));
%% export_def(uint64_p(obj));
%% export_def(sint64_p(obj));
%% export_def(uint_p);
%% export_def(sint_p);
%% export_def(ulong_p);
%% export_def(slong_p);


# ####################### TIMEBIBL in TIME.D ############################## #

# (* 25567 24 60 60) => 2208988800
# the number of seconds from 1900-01-01 to 1970-01-01
#define UNIX_LISP_TIME_DIFF 2208988800UL
%% export_def(UNIX_LISP_TIME_DIFF);

# Type which is used for 'Internal Time':
#ifdef TIME_1
  typedef uintL internal_time_t;      # measured value of the ticking counter
  #if defined(TIME_UNIX_TIMES)
    #define ticks_per_second  CLK_TCK
  #endif
  #define sub_internal_time(x,y, z)  z = (x) - (y)
  #define add_internal_time(x,y, z)  z = (x) + (y)
#endif
#ifdef TIME_2
  #ifdef TIME_UNIX
    typedef struct {
      uintL tv_sec;    # number of seconds since 1.1.1970 00:00 GMT,
                       # 'uintL' for tv_sec is good for 136 years.
      uintL tv_usec;   # additional microseconds
    } internal_time_t;
    #define ticks_per_second  1000000UL  # 1 Tick = 1 mu-sec
    #define sub_internal_time(x,y, z)  # z:=x-y  \
      do { (z).tv_sec = (x).tv_sec - (y).tv_sec;                \
        if ((x).tv_usec < (y).tv_usec)                          \
          { (x).tv_usec += ticks_per_second; (z).tv_sec -= 1; } \
        (z).tv_usec = (x).tv_usec - (y).tv_usec;                \
      } while(0)
    #define add_internal_time(x,y, z)  # z:=x+y  \
      do { (z).tv_sec = (x).tv_sec + (y).tv_sec;                \
        (z).tv_usec = (x).tv_usec + (y).tv_usec;                \
        if ((z).tv_usec >= ticks_per_second)                    \
          { (z).tv_usec -= ticks_per_second; (z).tv_sec += 1; } \
      } while(0)
  #endif
  #ifdef TIME_WIN32
    typedef # struct _FILETIME { DWORD dwLowDateTime; DWORD dwHighDateTime; }
            FILETIME  # number of 0.1 mu-sec since 1.1.1601 00:00 GMT.
            internal_time_t;
    #define ticks_per_second  10000000UL  # 1 Tick = 0.1 mu-sec
    #define sub_internal_time(x,y, z)  # z:=x-y  \
      do { (z).dwHighDateTime = (x).dwHighDateTime - (y).dwHighDateTime;      \
        if ((x).dwLowDateTime < (y).dwLowDateTime) { (z).dwHighDateTime -= 1;}\
        (z).dwLowDateTime = (x).dwLowDateTime - (y).dwLowDateTime;            \
      } while(0)
    #define add_internal_time(x,y, z)  # z:=x+y  \
      do { (z).dwHighDateTime = (x).dwHighDateTime + (y).dwHighDateTime;      \
        (z).dwLowDateTime = (x).dwLowDateTime + (y).dwLowDateTime;            \
        if ((z).dwLowDateTime < (x).dwLowDateTime) { (z).dwHighDateTime += 1;}\
      } while(0)
  #endif
#endif

#ifndef HAVE_RUN_TIME
# UP: Stops the run-time timer
  # run_time_stop();
  extern void run_time_stop (void);
  # is used by STREAM

  # UP: restarts the run-time timer
  # run_time_restart();
  extern void run_time_restart (void);
  # is used by STREAM
#else
  # You don't need a run-time timer
  #define run_time_stop()
  #define run_time_restart()
#endif

#ifdef TIME_1

# UP: Yields the real-time
# get_real_time()
# < uintL result: time since LISP-system-start (in 1/200 sec resp. in 1/50 sec resp. in 1/100 sec resp. in 1/CLK_TCK sec)
  extern uintL get_real_time (void);
# is used by STREAM, LISPARIT

#endif

#ifdef TIME_2

# UP: yields the real-time
# get_real_time()
# < internal_time_t* result: absolute time
  extern void get_real_time (internal_time_t*);
# is used by LISPARIT

#endif

# UP: Yields the run-time
# get_running_times(&timescore);
# < timescore.runtime:  Run-time since LISP-system-start (in Ticks)
# < timescore.realtime: Real-time since LISP-system-start (in Ticks)
# < timescore.gctime:   GC-Time since LISP-system-start (in Ticks)
# < timescore.gccount:  Number of GC's since LISP-system-start
# < timescore.gcfreed:  Size of the space reclaimed by the GC's so far
  typedef struct {
    internal_time_t runtime;
    internal_time_t realtime;
    internal_time_t gctime;
    uintL gccount;
    uintL2 gcfreed;
  } timescore_t;
  extern void get_running_times (timescore_t*);
# is used by

# UP: yields the run-time
# get_running_time(runtime);
# < runtime: Run-time (in Ticks)
  #ifndef HAVE_RUN_TIME
    #define get_running_time(runtime)  runtime = get_time()
    extern uintL get_time (void);
  #endif
  #if defined(TIME_UNIX) || defined(TIME_WIN32) || defined(TIME_UNIX_TIMES)
    #define get_running_time(runtime)  get_run_time(&runtime)
    #if defined(TIME_UNIX) || defined(TIME_WIN32)
      extern void get_run_time (internal_time_t* runtime);
    #endif
    #ifdef TIME_UNIX_TIMES
      extern uintL get_run_time (internal_time_t* runtime);
    #endif
  #endif
# is used by SPVW

# Time in decoded-time:
  typedef struct {
    object Sekunden;
    object Minuten;
    object Stunden;
    object Tag;
    object Monat;
    object Jahr;
  } decoded_time_t;

#ifdef UNIX
# UP: Converts the system-time-format into Decoded-Time.
# convert_time(&time,&timepoint);
# > time_t time: time in the system-time-format
# < timepoint.Sekunden, timepoint.Minuten, timepoint.Stunden,
#   timepoint.Tag, timepoint.Monat, timepoint.Jahr, each a Fixnum
  extern void convert_time (const time_t* time, decoded_time_t* timepoint);
# is used by PATHNAME
#endif
#ifdef WIN32_NATIVE
# UP: Converts the system-time-format into Decoded-Time.
# convert_time(&time,&timepoint);
# > FILETIME time: time in the system-time-format
# < timepoint.Sekunden, timepoint.Minuten, timepoint.Stunden,
#   timepoint.Tag, timepoint.Monat, timepoint.Jahr, each a Fixnum
  extern void convert_time (const FILETIME* time, decoded_time_t* timepoint);
# is used by PATHNAME
#endif

#ifdef UNIX
# UP: Converts the system time-format into Universal-Time.
# convert_time_to_universal(&time)
# > time_t time: time in the system time-format
# < result: integer denoting the seconds since 1900-01-01 00:00 GMT
# can trigger GC
  extern maygc object convert_time_to_universal (const time_t* time);
# is used by PATHNAME
#endif
#ifdef WIN32_NATIVE
# UP: converts the system time-format into Universal-Time.
# convert_time_to_universal(&time)
# > FILETIME time: Time in the system-time-format
# < result: integer denoting the seconds since 1900-01-01 00:00 GMT
# can trigger GC
  extern maygc object convert_time_to_universal (const FILETIME* time);
# is used by PATHNAME
#endif
%% #ifdef UNIX
%%   puts("extern object convert_time_to_universal (const time_t* time);");
%% #endif
%% #ifdef WIN32_NATIVE
%%   puts("extern object convert_time_to_universal (const FILETIME* time);");
%% #endif

#ifdef UNIX
/* the inverse of convert_time_to_universal() */
extern void convert_time_from_universal (object universal, time_t* time);
#endif
#ifdef WIN32_NATIVE
/* the inverse of convert_time_to_universal() */
extern void convert_time_from_universal (object universal, FILETIME* time);
#endif
%% #ifdef UNIX
%%   puts("extern void convert_time_from_universal (object universal, time_t* time);");
%% #endif
%% #ifdef WIN32_NATIVE
%%   puts("extern void convert_time_from_universal (object universal, FILETIME* time);");
%% #endif

# UP: Initializes the time variables upon the LISP-System-Start.
# init_time();
  extern void init_time (void);
# is used by SPVW


# ####################### SPVWBIBL for SPVW.D ############################## #

/*
                          The Stacks
                          ==========

Two Stacks are being used :
  - the C-program stack (Stackpointer SP = Register A7),
  - the LISP-Stack (Stackpointer STACK).
All calls of sub-programs are done through BSR/JSR via the program stack;
it's also used to temporarily store data, that is not a LISP-object.
The LISP-Stack is used to store frames and for the temporary storage
of LISP-objects.
For both stacks the limits of growth are controlled by the memory management
and the following macros:
  check_SP();             tests the program stack for overflow
  check_STACK();          tests the LISP-Stack for overflow
  get_space_on_STACK(n);  tests, whether there are still D0.L
                          Bytes free on the LISP-Stack
Basically only long words may be stored on the LISP-Stack.
If FRAME_BIT is set, it's the lower end of a frame;
this long word is a pointer above the Frame, together with a
Frame-type-Byte; if SKIP2_BIT is unset in it, the longword above
it is not a LISP-object.
All other long words on the LISP-Stack are LISP-objects.
*/

# machine stack: SP
# SP() returns the current value of the  SP.
# setSP(adresse); sets the SP to a given value. Extremely dangerous!
# FAST_SP defined, if SP-accesses are fast.
#if defined(GNU) && !(__APPLE_CC__ > 1)
  # definition of the register, in which the SP resides.
  #ifdef MC680X0
    #define SP_register "sp"  # %sp = %a7
  #endif
  #ifdef SPARC
    #define SP_register "%sp"  # %sp = %o6
  #endif
  #ifdef HPPA
    #define SP_register "%r30"  # %sp = %r30
  #endif
  #ifdef MIPS
    #define SP_register "$sp"  # $sp = $29
  #endif
  #ifdef M88000
    #define SP_register "%r31"  # %sp = %r31
  #endif
  #ifdef POWERPC
    #define SP_register "r1"
  #endif
  #ifdef ARM
    #define SP_register "%sp"  # %sp = %r13
  #endif
  #ifdef DECALPHA
    #define SP_register "$30"  # $sp = $30
  #endif
  #ifdef I80386
    #define SP_register "%esp"
  #endif
  #ifdef VAX
    #define SP_register "sp"
  #endif
  #ifdef IA64
    #define SP_register "r12"
  #endif
  #ifdef AMD64
    #define SP_register "%rsp"
  #endif
  #ifdef S390
    #define SP_register "15"
  #endif
#endif
#if (defined(GNU) || defined(INTEL)) && !defined(NO_ASM)
  # Assembler-instruction that copies the SP-register into a variable.
  #ifdef MC680X0
    #ifdef __REGISTER_PREFIX__ # GNU C Version >= 2.4 has %/ and __REGISTER_PREFIX__
      # But the value of __REGISTER_PREFIX__ is useless, because we might be
      # cross-compiling.
      #define REGISTER_PREFIX  "%/"
    #else
      #define REGISTER_PREFIX  "" # or "%%", depends on the assembler that's being used
    #endif
    #define ASM_get_SP_register(resultvar)  ("movel "REGISTER_PREFIX"sp,%0" : "=g" (resultvar) : )
  #endif
  #ifdef SPARC
    #ifdef SPARC64
      #define ASM_get_SP_register(resultvar)  ("add %%sp,2048,%0" : "=r" (resultvar) : )
    #else
      #define ASM_get_SP_register(resultvar)  ("mov %%sp,%0" : "=r" (resultvar) : )
    #endif
  #endif
  #ifdef HPPA
    #define ASM_get_SP_register(resultvar)  ("copy %%r30,%0" : "=r" (resultvar) : )
  #endif
  #ifdef MIPS
    #define ASM_get_SP_register(resultvar)  ("move\t%0,$sp" : "=r" (resultvar) : )
  #endif
  #ifdef M88000
    #define ASM_get_SP_register(resultvar)  ("or %0,#r0,#r31" : "=r" (resultvar) : )
  #endif
  #ifdef POWERPC
    #define ASM_get_SP_register(resultvar)  ("mr %0,1" : "=r" (resultvar) : )
  #endif
  #ifdef ARM
    #define ASM_get_SP_register(resultvar)  ("mov\t%0, sp" : "=r" (resultvar) : )
  #endif
  #ifdef DECALPHA
    #define ASM_get_SP_register(resultvar)  ("bis $30,$30,%0" : "=r" (resultvar) : )
  #endif
  #ifdef I80386
    #define ASM_get_SP_register(resultvar)  ("movl %%esp,%0" : "=g" (resultvar) : )
  #endif
  #ifdef IA64
    #define ASM_get_SP_register(resultvar)  ("mov %0 = r12" : "=r" (resultvar) : )
  #endif
  #ifdef AMD64
    #define ASM_get_SP_register(resultvar)  ("movq %%rsp,%0" : "=g" (resultvar) : )
  #endif
  #ifdef S390
    #define ASM_get_SP_register(resultvar)  ("lr %0,%%r15" : "=r" (resultvar) : )
  #endif
#endif
#if defined(GNU) && defined(MC680X0) && !defined(NO_ASM)
  # Access to a global register-"variable" SP
  #define SP()  \
    ({var aint __SP;                                                          \
      __asm__ __volatile__ ("movel "REGISTER_PREFIX"sp,%0" : "=g" (__SP) : ); \
      __SP;                                                                   \
     })
  #define setSP(adresse)  \
    ({ __asm__ __volatile__ ("movel %0,"REGISTER_PREFIX"sp" : : "g" ((aint)(adresse)) : "sp" ); })
  #define FAST_SP
#elif (defined(GNU) || defined(INTEL)) && defined(I80386) && !defined(NO_ASM)
  # Access to a register-"variable" %esp
  #define SP()  \
    ({var aint __SP;                                           \
      __asm__ __volatile__ ("movl %%esp,%0" : "=g" (__SP) : ); \
      __SP;                                                    \
     })
  # Doesn't work with gcc 3.1 any more.
  #if (__GNUC__ < 3) || (__GNUC__ == 3 && __GNUC_MINOR__ < 1)
    #define setSP(adresse)  \
      ({ __asm__ __volatile__ ("movl %0,%%esp" : : "g" ((aint)(adresse)) : "sp" ); })
    #define FAST_SP
  #endif
#elif defined(GNU) && defined(SP_register)
  register __volatile__ aint __SP __asm__(SP_register);
  #ifdef SPARC64
    #define SP()  (__SP+2048)
  #else
    #define SP()  __SP
  #endif
  #if defined(SPARC)
    # We must not do a setSP() here without taking care that
    # 1. %sp has to pay attention to an alignment of 8 Bytes,
    # 2. above %sp 92 Bytes have to be kept free (that's where the
    #    register contents are saved, if a 'register window overflow trap'
    #    is triggered by a 'save' in a sub-program).
  #endif
#elif defined(MICROSOFT) && defined(I80386) && !defined(NO_ASM)
  # access the register %esp
  #define SP  getSP
  static __inline aint getSP () { __asm mov eax,esp }
  static __inline aint setSP (aint address) { __asm mov esp,address }
#elif defined(MC680X0) || defined(SPARC) || defined(MIPS) || defined(I80386)
  # access functions extern, in assembler
  #define SP  getSP
  extern_C void* SP (void);
  extern_C void setSP (void* adresse);
#else
  # access function portable in C
  #define SP()  getSP()
  extern void* getSP (void);
  #define NEED_OWN_GETSP
#endif
#if defined(stack_grows_down) # defined(MC680X0) || defined(I80386) || defined(SPARC) || defined(MIPS) || defined(M88000) || defined(DECALPHA) || defined(IA64) || defined(AMD64) || defined(S390) || ...
  #define SP_DOWN # SP grows downward
  #define SPoffset 0 # top-of-SP ist *(SP+SPoffset)
#endif
#if defined(stack_grows_up) # defined(HPPA) || ...
  #define SP_UP # SP grows upward
  #define SPoffset -1 # top-of-SP ist *(SP+SPoffset)
#endif
#if (defined(SP_DOWN) && defined(SP_UP)) || (!defined(SP_DOWN) && !defined(SP_UP))
  #error "Unknown SP direction -- readjust SP_DOWN/SP_UP!"
#endif
# Derived from that:
# SPint  is the type of the elements on the SP, an Integer type at least as
#        wide as uintL and at least as wide as aint resp. void*.
# SP_(n) = (n+1)th longword on the SP.
# _SP_(n) = &SP_(n).
# pushSP(item)  puts a longword on the SP. Synonym: -(SP).
# popSP(item=)  returns item=SP_(0) and takes it off the SP.
# skipSP(n);  takes n long words of the SP.
#if (oint_addr_len <= intLsize)
  typedef uintL  SPint;
#else
  typedef aint  SPint;
#endif
#ifdef SP_DOWN
  #define skipSPop  +=
  #define SPop      +
#endif
#ifdef SP_UP
  #define skipSPop  -=
  #define SPop      -
#endif
#define _SP_(n)  (((SPint*)SP()) + SPoffset SPop (uintP)(n))
#if !(defined(GNU) && (defined(MC680X0)) && !defined(NO_ASM)) # generally
  #define SP_(n)  (((SPint*)SP())[SPoffset SPop (uintP)(n)])
  #define skipSP(n)                             \
    do { var register SPint* sp = (SPint*)SP(); \
         sp skipSPop (uintP)(n);                \
         setSP(sp);                             \
    } while(0)
  #define pushSP(item)                                                     \
    do { var register SPint* sp = (SPint*)SP();                            \
         sp skipSPop -1;                                                   \
         setSP(sp); # First decrease SP (because of a possible interrupt!) \
         sp[SPoffset] = (item); /* then insert item as top-of-SP */        \
    } while(0)
  #define popSP(item_assignment)                                        \
    do { var register SPint* sp = (SPint*)SP();                         \
         item_assignment sp[SPoffset]; # First fetch top-of-SP          \
         sp skipSPop 1;                                                 \
         setSP(sp); # then (danger of interrupt!) increase SP           \
    } while(0)
#endif
#if defined(GNU) && defined(MC680X0) && !defined(NO_ASM)
  # With GNU on as 680X0 SP is in a register. Thus access and
  # modification of SP are a unit that cannot be interrupted.
  # And SP_DOWN as well as SPoffset=0 hold.
  #define SP_(n)  \
    ({var register uintL __n = sizeof(SPint) * (n); \
      var register SPint __item;                    \
      __asm__ __volatile__ ("movel "REGISTER_PREFIX"sp@(%1:l),%0" : "=g" (__item) : "r" (__n) ); \
      __item;                                       \
     })
  #define skipSP(n)  \
    do { var register uintL __n = sizeof(SPint) * (n);  \
     __asm__ __volatile__ ("addl %0,"REGISTER_PREFIX"sp" : : "g" (__n) : "sp" ); \
    } while(0)
  #define pushSP(item)  \
    do { var register SPint __item = (item); \
     __asm__ __volatile__ ("movel %0,"REGISTER_PREFIX"sp@-" : : "g" (__item) : "sp" ); \
    } while(0)
  #define popSP(item_assignment)  \
    do {  var register SPint __item; \
     __asm__ __volatile__ ("movel "REGISTER_PREFIX"sp@+,%0" : "=r" (__item) : : "sp" ); \
     item_assignment __item;                                                 \
    } while(0)
#endif
# An sp_jmp_buf is exactly the same as a jmp_buf,
# except that on Irix 6.5 in 32-bit mode, a jmp_buf has alignment 8,
# whereas an SPint only has alignment 4.
# Need to add some padding.
# Then jmpbufsize = sizeof(sp_jmp_buf)/sizeof(SPint).
#define sp_jmp_buf_incr  (alignof(jmp_buf)>alignof(SPint)?alignof(jmp_buf)-alignof(SPint):0)
#define sp_jmp_buf_to_jmp_buf(x)  (*(jmp_buf*)(((long)&(x)+(long)sp_jmp_buf_incr)&-(long)(alignof(jmp_buf)>alignof(SPint)?alignof(jmp_buf):1)))
#define setjmpspl(x)  setjmpl(sp_jmp_buf_to_jmp_buf(x))
#define longjmpspl(x,y)  longjmpl(sp_jmp_buf_to_jmp_buf(x),y)
#define jmpbufsize  ceiling(sizeof(jmp_buf)+sp_jmp_buf_incr,sizeof(SPint))
typedef SPint sp_jmp_buf[jmpbufsize];
# The initial value of SP() during main().
extern void* SP_anchor;
%% #if (defined(GNU) || defined(INTEL)) && defined(I80386) && !defined(NO_ASM)
%%   printf("%s\n","#define SP()  ({aint __SP; __asm__ __volatile__ (\"movl %%esp,%0\" : \"=g\" (__SP) : ); __SP; })");
%% #endif

# LISP-Stack: STACK
#if !defined(STACK_register)
  # a global variable
  #ifndef MULTITHREAD
    extern gcv_object_t* STACK;
  #else
    #define STACK  (current_thread()->_STACK)
  #endif
#else
  # a global register variable
  register gcv_object_t* STACK __asm__(STACK_register);
#endif
#if defined(SPARC) && !defined(GNU) && !defined(__SUNPRO_C) && !defined(MULTITHREAD) && (SAFETY < 2)
  # a global register variable, but access functions externally in assembler
  #define STACK  _getSTACK()
  extern_C gcv_object_t* _getSTACK (void);
  #define setSTACK(allocation)  /* hem, yuck! */ \
    do { var gcv_object_t* tempSTACK; _setSTACK(temp##allocation); } while(0)
  extern_C void _setSTACK (void* new_STACK);
#else
  #define setSTACK(allocation)  allocation
#endif
#if defined(UNIX) || defined(WIN32) || defined(HYPERSTONE)
  #define STACK_UP # STACK grows upward
#endif
#if (defined(STACK_DOWN) && defined(STACK_UP)) || (!defined(STACK_DOWN) && !defined(STACK_UP))
  #error "Unknown STACK direction -- readjust STACK_DOWN/STACK_UP!"
#endif
%% #if !defined(STACK_register)
%%   puts("extern gcv_object_t* STACK;");
%% #else
%%   puts("#ifndef IN_MODULE_CC");
%%   printf("register gcv_object_t* STACK __asm__(\"%s\");\n",STACK_register);
%%   puts("#endif");
%% #endif

/* A singly-linked list of all currently active function calls.
   Resides in the C stack. */
struct backtrace_t {
  const struct backtrace_t* bt_next; /* Link to the caller */
  gcv_object_t bt_function;          /* Function or FSUBR being called */
  gcv_object_t *bt_stack;            /* STACK value where the frame area begins */
  int bt_num_arg;                    /* Number of arguments, if known, or -1 */
};
extern void back_trace_check (const struct backtrace_t *bt,
                              const char* label, const char* file, int line);
#ifdef DEBUG_BACKTRACE
#define BT_CHECK(b,l) back_trace_check(b,l,__FILE__,__LINE__)
#else
#define BT_CHECK(b,l)
#endif
#define BT_CHECK1(l)  BT_CHECK(back_trace,l)
%% puts("struct backtrace_t {\n  struct backtrace_t* bt_next;\n  gcv_object_t bt_function;\n  gcv_object_t *bt_stack;\n  int bt_num_arg;\n};");

#if defined(DEBUG_BACKTRACE) && defined(__cplusplus)
struct p_backtrace_t {
  const struct backtrace_t * ba_tr_p;
  p_backtrace_t (void* bt) { ba_tr_p = (struct backtrace_t*)bt; }
  /* assignment should check for circularities */
  p_backtrace_t& operator= (const struct backtrace_t *bt) {
    if (this->ba_tr_p != bt) {
      BT_CHECK(bt,"=: new value");
      BT_CHECK(ba_tr_p,"=: current value");
      this->ba_tr_p = bt;
    }
    return *this;
  };
  /* back_trace->foo means back_trace.ba_tr_p->foo */
  const struct backtrace_t* operator-> () {
    BT_CHECK(ba_tr_p,"->");
    return this->ba_tr_p;
  };
  /* cast p_backtrace_t to struct backtrace_t* */
  operator const struct backtrace_t* () const {
    BT_CHECK(ba_tr_p,"(struct backtrace_t*)");
    return ba_tr_p;
  }
};
#else
typedef const struct backtrace_t* p_backtrace_t;
#endif
%% emit_typedef("struct backtrace_t *","p_backtrace_t");

/* Returns the top-of-frame of a back_trace element. */
extern gcv_object_t* top_of_back_trace_frame (const struct backtrace_t *bt);

#define bt_beyond_stack_p(bt,st) \
  ((bt) != NULL && !((aint)(st) cmpSTACKop (aint)top_of_back_trace_frame(bt)))
/* unwind backtrace to the stack location */
#define unwind_back_trace(bt,st)                                        \
  do { BT_CHECK(bt,"unwind_back_trace");                                \
    while (bt_beyond_stack_p(bt,st))                                    \
      bt = bt->bt_next;                                                 \
  } while(0)

/* Evaluate statement, augmenting back_trace with an activation record for
   the given function.
   stack permits to locate the top-of-frame, namely
     - for FSUBRs:
         stack = top-of-frame - (req + opt + (body-flag ? 1 : 0))
     - for SUBRs:
         stack = top-of-frame - (req + opt + length(keyword-list))
     - for compiled closures:
         stack = top-of-frame - (req + opt + (rest-flag ? 1 : 0) + length(keyword-list))
     - for interpreted closures:
         stack = top-of-frame
 */
#if STACKCHECKS || STACKCHECKC
#define with_saved_back_trace(fun,stack,num_arg,statement)              \
  do {                                                                  \
    p_backtrace_t bt_save = back_trace;                                 \
    struct backtrace_t bt_here;                                         \
    bt_here.bt_next = back_trace;                                       \
    bt_here.bt_function = (fun);                                        \
    bt_here.bt_stack = (stack);                                         \
    bt_here.bt_num_arg = (num_arg);                                     \
    BT_CHECK1("w/s/b/t: before");                                       \
    back_trace = &bt_here;                                              \
    statement;                                                          \
    if (back_trace != &bt_here) abort();                                \
    if (back_trace->bt_next != bt_save) abort();                        \
    BT_CHECK1("w/s/b/t: after");                                        \
    back_trace = back_trace->bt_next;                                   \
  } while(0)
#else
#define with_saved_back_trace(fun,stack,num_arg,statement)              \
  do {                                                                  \
    struct backtrace_t bt_here;                                         \
    bt_here.bt_next = back_trace;                                       \
    bt_here.bt_function = (fun);                                        \
    bt_here.bt_stack = (stack);                                         \
    bt_here.bt_num_arg = (num_arg);                                     \
    back_trace = &bt_here;                                              \
    statement;                                                          \
    back_trace = back_trace->bt_next;                                   \
  } while(0)
#endif
#define with_saved_back_trace_fsubr(fun,statement)  \
  with_saved_back_trace(fun,STACK,-1,statement)
#define with_saved_back_trace_subr(fun,stack,num_arg,statement)  \
  with_saved_back_trace(fun,stack,num_arg,statement)
#define with_saved_back_trace_cclosure(fun,statement)  \
  with_saved_back_trace(fun,STACK,-1,statement)
#define with_saved_back_trace_iclosure(fun,stack,num_arg,statement)  \
  with_saved_back_trace(fun,stack,num_arg,statement)

# Every call of an external function (or a sequence of those) has to be framed
# with
#   begin_call();
# and
#   end_call();
# Purpose: The stack, if it resides in a register,
# should be brought to a halfway recent value
# in case of an interrupt during the corresponding timespan.
#
# If you want to access the STACK while an external function run,
# you have to frame the corresponding code with
#   begin_callback();
# and
#   end_callback();
#ifdef HAVE_SAVED_mv_count
  #ifndef MULTITHREAD
    extern uintC saved_mv_count;
  #else
    #define saved_mv_count  (current_thread()->_saved_mv_count)
  #endif
  #define SAVE_mv_count()     saved_mv_count = mv_count
  #define RESTORE_mv_count()  mv_count = saved_mv_count
#else
  #define SAVE_mv_count()
  #define RESTORE_mv_count()
#endif
#ifdef HAVE_SAVED_value1
  #ifndef MULTITHREAD
    extern object saved_value1;
  #else
    #define saved_value1  (current_thread()->_saved_value1)
  #endif
  #define SAVE_value1()     saved_value1 = value1
  #define RESTORE_value1()  value1 = saved_value1
#else
  #define SAVE_value1()
  #define RESTORE_value1()
#endif
#ifdef HAVE_SAVED_back_trace
  #ifndef MULTITHREAD
    extern p_backtrace_t saved_back_trace;
  #else
    #define saved_back_trace  (current_thread()->_saved_back_trace)
  #endif
  #define SAVE_back_trace()     saved_back_trace = back_trace
  #define RESTORE_back_trace()  back_trace = saved_back_trace
#else
  #define SAVE_back_trace()
  #define RESTORE_back_trace()
#endif
#define SAVE_GLOBALS()     SAVE_mv_count(); SAVE_value1(); SAVE_back_trace();
#define RESTORE_GLOBALS()  RESTORE_mv_count(); RESTORE_value1(); RESTORE_back_trace();
#if defined(HAVE_SAVED_STACK)
  #ifndef MULTITHREAD
    extern gcv_object_t* saved_STACK;
  #else
    #define saved_STACK  (current_thread()->_saved_STACK)
  #endif
  #define begin_call()  SAVE_GLOBALS(); saved_STACK = STACK
  #define end_call()  RESTORE_GLOBALS(); saved_STACK = (gcv_object_t*)NULL
  #define begin_callback()  SAVE_REGISTERS( STACK = saved_STACK; ); end_call()
  #define end_callback()  SAVE_GLOBALS(); RESTORE_REGISTERS( saved_STACK = STACK; )
#else
  #define begin_call()  SAVE_GLOBALS()
  #define end_call()  RESTORE_GLOBALS()
  #define begin_callback()  SAVE_REGISTERS(;); end_call()
  #define end_callback()  SAVE_GLOBALS(); RESTORE_REGISTERS(;)
#endif
%% #ifdef HAVE_SAVED_mv_count
%%   puts("extern uintC saved_mv_count;");
%% #endif
%% #ifdef HAVE_SAVED_value1
%%   puts("extern object saved_value1;");
%% #endif
%% #ifdef HAVE_SAVED_back_trace
%%   puts("extern p_backtrace_t saved_back_trace;");
%% #endif
%% #if defined(HAVE_SAVED_STACK)
%%   puts("extern gcv_object_t* saved_STACK;");
%% #endif
%% export_def(begin_call());
%% export_def(end_call());
%% export_def(begin_callback());
%% export_def(end_callback());

# Every OS-call (or a sequence thereof) has to be framed with
#   begin_system_call();
# and
#   end_system_call();
# Purpose: The STACK - if it resides in a register -
# should be brought to a halfway recent value,
# if an interrupt happens during the corresponding timespan.
#
# While a break-semaphore has been set, you don't have to use the macros
# because of that.
#ifdef NO_ASYNC_INTERRUPTS
  # NO_ASYNC_INTERRUPTS: if we don't react to asynchronous Interrupts,
  #   the program can't be interruped..
  #define begin_system_call()
  #define end_system_call()
#else
  #define begin_system_call()  begin_call()
  #define end_system_call()  end_call()
#endif
# The same holds for setjmp()/longjmp(). Here we avoid an unneeded overhead
# if at all possible.
# You don't have to use these macros when a break-semaphore has been
# set.
#if 0
  # Disassembly of setjmp() and longjmp() shows, that the STACK-register
  # isn't used arbitrarily.
  #define begin_setjmp_call()
  #define end_setjmp_call()
  #define begin_longjmp_call()
  #define end_longjmp_call()
#elif defined(I80386) && (defined(UNIX_LINUX) || defined(UNIX_GNU))
  # Disassembly of setjmp() shows, that the STACK-register %ebx
  # isn't used arbitrarily.
  #define begin_setjmp_call()
  #define end_setjmp_call()
  #define begin_longjmp_call()  begin_system_call()
  #define end_longjmp_call()  end_system_call()
#else
  #define begin_setjmp_call()  begin_system_call()
  #define end_setjmp_call()  end_system_call()
  #define begin_longjmp_call()  begin_system_call()
  #define end_longjmp_call()  end_system_call()
#endif
# The same holds for arithmetics-functions that use the STACK_registers.
# On I80386 (%ebx) these are SHIFT_LOOPS, MUL_LOOPS, DIV_LOOPS.
#if defined(I80386) && !defined(NO_ARI_ASM) && defined(HAVE_SAVED_STACK)
  #define begin_arith_call()  begin_system_call()
  #define end_arith_call()  end_system_call()
#else
  #define begin_arith_call()
  #define end_arith_call()
#endif
%% export_def(begin_system_call());
%% export_def(end_system_call());

#if defined(HAVE_STACK_OVERFLOW_RECOVERY)
  # Detection of SP-overflow through a Guard-Page or other mechanisms.
  #define NOCOST_SP_CHECK
#else
  # The OS is responsible for the SP.
  # From where should we get a reasonable value for SP_bound?
  #define NO_SP_CHECK
#endif

# Tests for SP-overflow.
# check_SP();            tests for overflow
# check_SP_notUNIX();    dito, except when a temporary overflow doesn't matter
#define check_SP()  if (SP_overflow()) SP_ueber()
#if !(defined(NO_SP_CHECK) || defined(NOCOST_SP_CHECK))
  #ifdef SP_DOWN
    #define SP_overflow()  ( (aint)SP() < (aint)SP_bound )
  #endif
  #ifdef SP_UP
    #define SP_overflow()  ( (aint)SP() > (aint)SP_bound )
  #endif
#else # NO_SP_CHECK || NOCOST_SP_CHECK
  #define SP_overflow()  false
  #ifdef NOCOST_SP_CHECK
    #ifdef WIN32_NATIVE
      #ifdef SP_DOWN
        #define near_SP_overflow()  ( (aint)SP() < (aint)SP_bound+0x1000 )
      #endif
      #ifdef SP_UP
        #define near_SP_overflow()  ( (aint)SP() > (aint)SP_bound-0x1000 )
      #endif
    #else
      extern bool near_SP_overflow (void);
    #endif
  #endif
#endif
#ifndef MULTITHREAD
  extern void* SP_bound;
#else
  #define SP_bound  (current_thread()->_SP_bound)
#endif
nonreturning_function(extern, SP_ueber, (void));
#ifdef UNIX
  #define check_SP_notUNIX()
#else
  #define check_SP_notUNIX()  check_SP()
#endif

# Tests for STACK-overflow.
# check_STACK();
#define check_STACK()  if (STACK_overflow()) STACK_ueber()
#ifdef STACK_DOWN
  #define STACK_overflow()  ( (aint)STACK < (aint)STACK_bound )
#endif
#ifdef STACK_UP
  #define STACK_overflow()  ( (aint)STACK > (aint)STACK_bound )
#endif
#ifndef MULTITHREAD
  extern void* STACK_bound;
#else
  #define STACK_bound  (current_thread()->_STACK_bound)
#endif
nonreturning_function(extern, STACK_ueber, (void));
%% #if notused
%% export_def(check_STACK());
%% export_def(STACK_overflow());
%% export_def(get_space_on_STACK(n));
%% puts("extern void* STACK_bound;");
%% puts("nonreturning_function(extern, STACK_ueber, (void));");
%% #endif

# Tests, if there are still n Bytes free on the STACK.
# get_space_on_STACK(n);
#ifdef STACK_DOWN
  #define get_space_on_STACK(n)  \
    if ( (aint)STACK < (aint)STACK_bound + (aint)(n) ) STACK_ueber()
#else
  #define get_space_on_STACK(n)  \
    if ( (aint)STACK + (aint)(n) > (aint)STACK_bound ) STACK_ueber()
#endif

# Exit the LISP-Interpreter
# quit();
# > final_exitcode: 0 for a normal end, >0 for failure
nonreturning_function(extern, quit, (void));
extern int final_exitcode;
# is used by CONTROL

# Error message if an unreachable program part has been reached.
# Does not return.
# fehler_notreached(file,line);
# > file: Filename (with quotation marks) as constant ASCIZ-String
# > line: line number
nonreturning_function(extern, fehler_notreached, (const char * file, uintL line));
# used by all modules
%% puts("nonreturning_function(extern, fehler_notreached, (const char * file, uintL line));");

# Language that's used to communicate with the user:
#ifdef LANGUAGE_STATIC
  #if ENGLISH
    #define GETTEXT(english)   english
    #define GETTEXTL(english)  english
  #endif
#else
  #define language_english   0
  #ifndef GNU_GETTEXT
    # Language is determined at runtime by the variable language.
    extern uintL language;
    #define ENGLISH  (language==language_english)
    #define GETTEXT(english)  english
    #define GETTEXTL(english)  english
  #else # GNU_GETTEXT
    #ifndef COMPILE_STANDALONE
      #include <libintl.h>
    #endif
    # Fetch the message translations from a message catalog.
    #ifndef gettext  # Sometimes `gettext' is a macro...
      extern char* gettext (const char * msgid);
    #endif
    extern const char * clgettext (const char * msgid);
    extern const char * clgettextl (const char * msgid);
    # GETTEXT(english_message) fetches the translation of english_message
    # and returns it in UTF-8 (if UNICODE is defined).
    # GETTEXTL(english_message) fetches the translation of english_message
    # and returns it in the locale encoding.
    # GETTEXT and GETTEXTL are special tags recognized by clisp-xgettext. We
    # choose English because it's the only language understood by all CLISP
    # developers.
    #define GETTEXT  clgettext
    #define GETTEXTL clgettextl
  #endif

  # init the language and the locale
  extern void init_language (const char*, const char*);
#endif
%% #if !defined(LANGUAGE_STATIC) && defined(GNU_GETTEXT)
%%   puts("#define GNU_GETTEXT");
%%   puts("#ifndef COMPILE_STANDALONE");
%%   puts("#include <libintl.h>");
%%   puts("#endif");
%%   puts("extern const char * clgettext (const char * msgid);");
%%   export_def(GETTEXT);
%% #else
%%   export_def(GETTEXT(english));
%% #endif

# Fetch the message translations of a string: "CL String getTEXT"
# CLSTEXT(string)
# > obj: C string
# < result: String
# can trigger GC
extern maygc object CLSTEXT (const char*);
%% #ifndef LANGUAGE_STATIC
%%   #ifndef GNU_GETTEXT
%%     emit_define("CLSTEXT","ascii_to_string");
%%   #else
%%     puts("extern object CLSTEXT (const char* asciz);");
%%   #endif
%% #endif

# Fetch the "translation" of a Lisp object: "CL Object getTEXT"
# CLOTEXT(string)
# > obj: String
# can trigger GC
extern maygc object CLOTEXT (const char*);

# Print a Lisp object in Lisp notation relatively directly
# through the operating system:
# object_out(obj);
# can trigger GC
extern maygc object object_out (object obj);
# can trigger GC
# print the object with label, file name and line number
# this can trigger GC, but will save and restore OBJ
#define OBJECT_OUT(obj,label)                                           \
  (printf("[%s:%d] %s: %s:\n",__FILE__,__LINE__,STRING(obj),label),     \
   obj=object_out(obj))
/* print the object to a C stream - not all objects can be handled yet!
 non-consing, STACK non-modifying */
extern maygc object nobject_out (FILE* out, object obj);
# used for debugging purposes
%% puts("extern object object_out (object obj);");
%% puts("#define OBJECT_OUT(obj,label)  (printf(\"[%s:%d] %s: %s:\\n\",__FILE__,__LINE__,STRING(obj),label),obj=object_out(obj))");

# After allocating memory for an object, add the type infos.
#ifdef TYPECODES
  #define bias_type_pointer_object(bias,type,ptr) type_pointer_object(type,ptr)
#else
  #ifdef WIDE_AUXI
    #define bias_type_pointer_object(bias,type,ptr) as_object_with_auxi((aint)(ptr)+(bias))
  #else
    #define bias_type_pointer_object(bias,type,ptr) as_object((oint)(ptr)+(bias))
  #endif
#endif
# used by SPVW, macros SP_allocate_bit_vector, SP_allocate_string

# UP: executes a Garbage Collection
# gar_col();
# can trigger GC
extern maygc void gar_col(void);
# is used by DEBUG

# GC-statistics
extern uintL gc_count;
extern uintL2 gc_space;
extern internal_time_t gc_time;
# is used by TIME

# UP:  allocates a Cons
# allocate_cons()
# < result: pointer to a new CONS, with CAR and CDR =NIL
# can trigger GC
extern maygc object allocate_cons (void);
# is used by LIST, SEQUENCE, PACKAGE, EVAL, CONTROL, RECORD,
#            PREDTYPE, IO, STREAM, PATHNAME, SYMBOL, ARRAY, LISPARIT
%% puts("extern object allocate_cons (void);");

# UP: Returns a newly created uninterned symbol with a given Printname.
# make_symbol(string)
# > string: immutable Simple-String
# < result: new symbol with this name, with Home-Package=NIL.
# can trigger GC
extern maygc object make_symbol (object string);
# is used by PACKAGE, IO, SYMBOL
%% #if notused
%% puts("extern object make_symbol (object string);");
%% #endif

# UP: allocates a general vector
# allocate_vector(len)
# > len: length of the vector
# < result: fresh simple general vector (elements are initialized with NIL)
# can trigger GC
extern maygc object allocate_vector (uintL len);
# is used by ARRAY, IO, EVAL, PACKAGE, CONTROL, HASHTABL
%% puts("extern object allocate_vector (uintL len);");

# Function: Allocates a bit/byte vector.
# allocate_bit_vector(atype,len)
# > uintB atype: Atype_nBit
# > uintL len: length (number of n-bit blocks)
# < result: fresh simple bit/byte-vector of the given length
# can trigger GC
extern maygc object allocate_bit_vector (uintB atype, uintL len);
# is used by ARRAY, IO, RECORD, LISPARIT, STREAM, CLX
%% puts("extern object allocate_bit_vector (uintB atype, uintL len);");

# Macro: Allocates a 8bit-vector on the stack, with dynamic extent.
#   { var DYNAMIC_8BIT_VECTOR(obj,len);
#     ...
#     FREE_DYNAMIC_8BIT_VECTOR(obj);
#   }
# > uintL len: length (number of bytes)
# < object obj: simple-8bit-vector with dynamic extent
#   (may or may not be heap-allocated, therefore not GC-invariant)
# can trigger GC
#if defined(SPVW_PURE) || ((((STACK_ADDRESS_RANGE << addr_shift) >> garcol_bit_o) & 1) != 0)
  # No way to allocate a Lisp object on the stack.
  #define DYNAMIC_8BIT_VECTOR(objvar,len)  \
    var uintL objvar##_len = (len);               \
    var object objvar = O(dynamic_8bit_vector);   \
    O(dynamic_8bit_vector) = NIL;                 \
    if (!(simple_bit_vector_p(Atype_8Bit,objvar) && (Sbvector_length(objvar) >= objvar##_len))) \
      objvar = allocate_bit_vector(Atype_8Bit,objvar##_len); \
    GCTRIGGER1(objvar)
  #define FREE_DYNAMIC_8BIT_VECTOR(objvar)  \
    O(dynamic_8bit_vector) = objvar
#else
  # Careful: Fill GCself with pointers to itself, so that GC will leave
  # pointers to this object untouched.
  #ifdef TYPECODES
    #define DYNAMIC_8BIT_VECTOR(objvar,len)  \
      DYNAMIC_ARRAY(objvar##_storage,object,ceiling((uintL)(len)+offsetofa(sbvector_,data),sizeof(gcv_object_t))); \
      var object objvar = ((Sbvector)objvar##_storage)->GCself = bias_type_pointer_object(varobject_bias,sb8vector_type,(Sbvector)objvar##_storage); \
      ((Sbvector)objvar##_storage)->length = (len); \
      GCTRIGGER1(objvar)
  #else
    #define DYNAMIC_8BIT_VECTOR(objvar,len)  \
      DYNAMIC_ARRAY(objvar##_storage,object,ceiling((uintL)(len)+offsetofa(sbvector_,data)+varobjects_misaligned,sizeof(gcv_object_t))); \
      var object* objvar##_address = (object*)((uintP)objvar##_storage | varobjects_misaligned); \
      var object objvar = ((Sbvector)objvar##_address)->GCself = bias_type_pointer_object(varobject_bias,sb8vector_type,(Sbvector)objvar##_address); \
      ((Sbvector)objvar##_address)->tfl = vrecord_tfl(Rectype_Sb8vector,len); \
      GCTRIGGER1(objvar)
  #endif
  #define FREE_DYNAMIC_8BIT_VECTOR(objvar)  \
    FREE_DYNAMIC_ARRAY(objvar##_storage)
#endif
# used by STREAM, PATHNAME

# Macro: Wraps a GC-invariant uintB* pointer in a fake simple-8bit-vector.
# FAKE_8BIT_VECTOR(ptr)
# > uintB* ptr: pointer to GC-invariant data
# < gcv_object_t obj: a fake simple-8bit-vector,
#                     with TheSbvector(obj)->data == ptr,
#                     that must *not* be stored in GC-visible locations
#ifdef TYPECODES
  #define FAKE_8BIT_VECTOR(ptr)  \
    type_pointer_object(0, (const char*)(ptr) - offsetofa(sbvector_,data))
#else
  #define FAKE_8BIT_VECTOR(ptr)  \
    fake_gcv_object((aint)((const char*)(ptr) - offsetofa(sbvector_,data)) + varobject_bias)
#endif

#if !defined(UNICODE) || defined(HAVE_SMALL_SSTRING)
# UP, provides 8-bit character string
# allocate_s8string(len)
# > len: length of the string (in characters), must be <= stringsize_limit_1
# < result: new 8-bit character simple-string (LISP-object)
# can trigger GC
extern maygc object allocate_s8string (uintL len);
# used by
#endif
#if defined(UNICODE) && !defined(HAVE_SMALL_SSTRING)
#define allocate_s8string(len)  allocate_s32string(len)
#endif
%% #if !defined(UNICODE)
%%   puts("extern object allocate_s8string (uintL len);");
%% #endif

#if !defined(UNICODE) || defined(HAVE_SMALL_SSTRING)
# UP, provides immutable 8-bit character string
# allocate_imm_s8string(len)
# > len: length of the string (in characters), must be <= stringsize_limit_1
# < result: new immutable 8-bit character simple-string (LISP-object)
# can trigger GC
extern maygc object allocate_imm_s8string (uintL len);
# used by
#endif

#ifdef HAVE_SMALL_SSTRING
# UP, provides 16-bit character string
# allocate_s16string(len)
# > len: length of the string (in characters), must be <= stringsize_limit_1
# < result: new 16-bit character simple-string (LISP-object)
# can trigger GC
extern maygc object allocate_s16string (uintL len);
# used by
#endif
#if defined(UNICODE) && !defined(HAVE_SMALL_SSTRING)
#define allocate_s16string(len)  allocate_s32string(len)
#endif

#ifdef HAVE_SMALL_SSTRING
# UP, provides immutable 16-bit character string
# allocate_imm_s16string(len)
# > len: length of the string (in characters), must be <= stringsize_limit_1
# < result: new immutable 16-bit character simple-string (LISP-object)
# can trigger GC
extern maygc object allocate_imm_s16string (uintL len);
# used by
#endif

#ifdef UNICODE
# UP, provides 32-bit character string
# allocate_s32string(len)
# > len: length of the string (in characters), must be <= stringsize_limit_1
# < result: new 32-bit character simple-string (LISP-object)
# can trigger GC
extern maygc object allocate_s32string (uintL len);
#endif
%% #ifdef UNICODE
%%   puts("extern object allocate_s32string (uintL len);");
%% #endif

#ifdef UNICODE
# UP, provides immutable 32-bit character string
# allocate_imm_s32string(len)
# > len: length of the string (in characters), must be <= stringsize_limit_1
# < result: new immutable 32-bit character simple-string (LISP-object)
# can trigger GC
extern maygc object allocate_imm_s32string (uintL len);
#endif

# UP: allocates String
# allocate_string(len)
# > len: length of the Strings (in Characters), must be <= stringsize_limit_1
# < result: new Normal-Simple-String (LISP-object)
# can trigger GC
#ifdef UNICODE
  #define allocate_string(len)  allocate_s32string(len)
#else
  #define allocate_string(len)  allocate_s8string(len)
#endif
# is used by ARRAY, CHARSTRG, STREAM, PATHNAME
%% export_def(allocate_string(len));

# Macro: Allocates a normal string on the stack, with dynamic extent.
#   { var DYNAMIC_STRING(obj,len);
#     ...
#     FREE_DYNAMIC_STRING(obj);
#   }
# > uintL len: length (number of characters)
# < object obj: normal-simple-string with dynamic extent
#   (may or may not be heap-allocated, therefore not GC-invariant)
# can trigger GC
#if defined(SPVW_PURE) || ((((STACK_ADDRESS_RANGE << addr_shift) >> garcol_bit_o) & 1) != 0)
  # No way to allocate a Lisp object on the stack.
  #define DYNAMIC_STRING(objvar,len)  \
    var uintL objvar##_len = (len);           \
    var object objvar = O(dynamic_string);    \
    O(dynamic_string) = NIL;                  \
    if (!(simple_string_p(objvar) && (Sstring_length(objvar) >= objvar##_len))) { \
      if (objvar##_len > stringsize_limit_1)  \
        fehler_stringsize(objvar##_len);      \
      objvar = allocate_string(objvar##_len); \
    }                                         \
    GCTRIGGER1(objvar)
  #define FREE_DYNAMIC_STRING(objvar)  \
    O(dynamic_string) = objvar;
#else
  # Careful: Fill GCself with pointers to itself, so that GC will leave
  # pointers to this object untouched.
  #ifdef UNICODE
    #define DYNAMIC_STRING(objvar,len)  \
      DYNAMIC_ARRAY(objvar##_storage,object,ceiling((uintL)(len)*sizeof(chart)+offsetofa(s32string_,data)+varobjects_misaligned,sizeof(gcv_object_t))); \
      var object* objvar##_address = (object*)((uintP)objvar##_storage | varobjects_misaligned); \
      var object objvar = ((Sstring)objvar##_address)->GCself = bias_type_pointer_object(varobject_bias,sstring_type,(Sstring)objvar##_address); \
      ((Sstring)objvar##_address)->tfl = sstring_tfl(Sstringtype_32Bit,0,0,len); \
      GCTRIGGER1(objvar)
  #else
    #define DYNAMIC_STRING(objvar,len)  \
      DYNAMIC_ARRAY(objvar##_storage,object,ceiling((uintL)(len)*sizeof(chart)+offsetofa(s8string_,data)+varobjects_misaligned,sizeof(gcv_object_t))); \
      var object* objvar##_address = (object*)((uintP)objvar##_storage | varobjects_misaligned); \
      var object objvar = ((Sstring)objvar##_address)->GCself = bias_type_pointer_object(varobject_bias,sstring_type,(Sstring)objvar##_address); \
      ((Sstring)objvar##_address)->tfl = sstring_tfl(Sstringtype_8Bit,0,0,len); \
      GCTRIGGER1(objvar)
  #endif
  #define FREE_DYNAMIC_STRING(objvar)  \
    FREE_DYNAMIC_ARRAY(objvar##_storage)
#endif
# used by LISPARIT

# UP: allocates an immutable String
# allocate_imm_string(len)
# > len: length of the String (in Characters)
# < result: new immutable Normal-Simple-String (LISP-object)
# can trigger GC
#ifdef UNICODE
  #define allocate_imm_string(len)  allocate_imm_s32string(len)
#else
  #define allocate_imm_string(len)  allocate_imm_s8string(len)
#endif
# is used by CHARSTRG

#ifdef HAVE_SMALL_SSTRING
# UP: Changes the allocation of a Small-String to an Sistring, while
# copying the contents to a fresh normal string.
# reallocate_small_string(string)
# > string: a nonempty Small-String
# > newtype: new wider string type, Sstringtype_16Bit or Sstringtype_32Bit
# < result: an Sistring pointing to a wider String
# can trigger GC
  extern maygc object reallocate_small_string (object string, uintB newtype);
# is used by ARRAY
#endif

# Attempts to reallocate a simple-string, for debugging purposes.
# DBGREALLOC(string);
#if defined(DEBUG_SMALL_SSTRING) && defined(HAVE_SMALL_SSTRING)
  #define DBGREALLOC(string)  \
    if (simple_string_p(string) && !sstring_reallocatedp(TheSstring(string)) \
        && !sstring_immutable(TheSstring(string))                            \
        && sstring_eltype(TheSstring(string)) != Sstringtype_32Bit           \
        && sstring_length(TheSstring(string)) > 0)                           \
      string = reallocate_small_string(string,sstring_eltype(TheSstring(string))+1)/*;*/
#else
  #define DBGREALLOC(string)  (void)0 /*nop*/
#endif

# UP: allocates indirect array
# allocate_iarray(flags,rank,type)
# > uintB flags: Flags
# > uintC (actually uintWC) rank: rank
# > tint type: Typinfo
# < result: LISP-object Array
# can trigger GC
extern maygc object allocate_iarray (uintB flags, uintC rank, tint type);
# is used by ARRAY, IO

# UP: allocates Long-Record
# allocate_lrecord(rectype,reclen,type)
# > sintB rectype: further type-info
# > uintL reclen: length
# > tint type: type-info
# < result: LISP-object Record (elements are initialized with NIL)
# can trigger GC
#ifdef TYPECODES
  extern maygc object allocate_lrecord (uintB rectype, uintL reclen, tint type);
#else
  #define allocate_lrecord(rectype,reclen,type)  /* ignore type */ \
    allocate_lrecord_(rectype,reclen)
  extern object allocate_lrecord_ (uintB rectype, uintL reclen);
#endif
# is used by WEAK

# UP: allocates Simple-Record
# allocate_srecord(flags,rectype,reclen,type)
# > uintB flags: Flags
# > sintB rectype: further type-info
# > uintC (actually uintW) reclen: length
# > tint type: type-info
# < result: LISP-object Record (elements are initialized with NIL)
# can trigger GC
#ifdef TYPECODES
  #define allocate_srecord(flags,rectype,reclen,type)  \
    allocate_srecord_(                                                     \
       (BIG_ENDIAN_P ? (uintW)(flags)+((uintW)(uintB)(rectype)<<intBsize)  \
                     : ((uintW)(flags)<<intBsize)+(uintW)(uintB)(rectype)),\
       reclen,                                                             \
       type)
  extern maygc object allocate_srecord_ (uintW flags_rectype, uintC reclen, tint type);
#else
  #define allocate_srecord(flags,rectype,reclen,type)  /* ignore type */ \
    allocate_srecord_(((uintW)(flags)<<8)+(uintW)(uintB)(rectype),reclen)
  extern maygc object allocate_srecord_ (uintW flags_rectype, uintC reclen);
#endif
# is used by RECORD, EVAL

# UP: allocates Extended-Record
# allocate_xrecord(flags,rectype,reclen,recxlen,type)
# > uintB flags: Flags
# > sintB rectype: further type-info
# > uintC (actually uintB) reclen: length
# > uintC (actually uintB) recxlen: extra-length
# > tint type: Typinfo
# < result: LISP-object Record (elements are initialized with NIL resp. 0)
# can trigger GC
#ifdef TYPECODES
  #define allocate_xrecord(flags,rectype,reclen,recxlen,type)  \
    allocate_xrecord_(                                                     \
       (BIG_ENDIAN_P ? (uintW)(flags)+((uintW)(uintB)(rectype)<<intBsize)  \
                     : ((uintW)(flags)<<intBsize)+(uintW)(uintB)(rectype)),\
       reclen,                                                             \
       recxlen,                                                            \
       type)
  extern maygc object allocate_xrecord_ (uintW flags_rectype, uintC reclen, uintC recxlen, tint type);
#else
  #define allocate_xrecord(flags,rectype,reclen,recxlen,type)  \
    allocate_xrecord_((((uintW)(flags)<<8)+(uintW)(uintB)(rectype)),reclen,recxlen)
  extern maygc object allocate_xrecord_ (uintW flags_rectype, uintC reclen, uintC recxlen);
#endif
# is used by

# UP: allocates Closure
# allocate_closure(reclen)
# > uintC reclen: length
# < result: LISP-object Closure (elements are initialized with NIL)
#define allocate_closure(reclen,flags)                               \
  allocate_srecord(flags,Rectype_Closure,reclen,closure_type)
# is used by EVAL, RECORD

/* copy a section of memory */
#define copy_mem_b(dest,orig,len) /* bytes */                   \
  do { var char* newptr = (char*)(dest);                        \
       var const char* oldptr = (const char*)(orig);            \
       var uintL count;                                         \
       var uintL leng = (len);                                  \
       dotimespL(count,leng,{ *newptr++ = *oldptr++; });        \
  } while(0)
#define copy_mem_o(dest,orig,len) /* objects */                  \
  do { var gcv_object_t* newptr = (dest);                \
       var const gcv_object_t* oldptr = (orig);          \
       var uintC count;                                  \
       var uintC leng = (len);                           \
       dotimespC(count,leng,{ *newptr++ = *oldptr++; }); \
  } while(0)
#if 0 /* the libc alternative turns out to be ~3-5% slower */
#define copy_mem_b(dest,orig,len)                       \
    do { begin_system_call(); memcpy(dest,orig,len);    \
         end_system_call(); } while(0)
#define copy_mem_o(dest,orig,len)                                           \
    do { begin_system_call(); memcpy(dest,orig,(len)*sizeof(gcv_object_t)); \
         end_system_call(); } while(0)
#endif

# Copying a compiled closure:
# newclos = allocate_cclosure_copy(oldclos);
# can trigger GC
#define allocate_cclosure_copy(oldclos)  \
  allocate_closure(Cclosure_length(oldclos),Closure_flags(oldclos))
# do_cclosure_copy(newclos,oldclos);
#define do_cclosure_copy(newclos,oldclos)               \
  copy_mem_o(((Srecord)TheCclosure(newclos))->recdata,  \
             ((Srecord)TheCclosure(oldclos))->recdata,  \
             Cclosure_length(oldclos))
# is used by EVAL, IO, RECORD

# UP: allocates Structure
# allocate_structure(reclen)
# > uintC reclen: length
# < result: LISP-Object Structure (Elements are initialized with NIL)
# can trigger GC
#ifdef case_structure
  #define allocate_structure(reclen)  \
    allocate_srecord(0,Rectype_Structure,reclen,structure_type)
#else
  #define allocate_structure(reclen)  \
    allocate_srecord(0,Rectype_Structure,reclen,orecord_type)
#endif
# is used by RECORD

# UP: allocates Stream
# allocate_stream(strmflags,strmtype,reclen,recxlen)
# > uintB strmflags: Flags
# > uintB strmtype: further type-info
# > uintC reclen: length in objects
# > uintC recxlen: extra-length in bytes
# < result: LISP-object Stream (elements are initialized with NIL)
# can trigger GC
#ifdef case_stream
  #define allocate_stream(strmflags,strmtype,reclen,recxlen)  \
    allocate_xrecord(strmflags | strmflags_open_B,strmtype,reclen,recxlen,stream_type)
#else
  extern maygc object allocate_stream (uintB strmflags, uintB strmtype, uintC reclen, uintC recxlen);
#endif
# is used by STREAM

# UP: allocates Package
# allocate_package()
# < result: LISP-object Package
# can trigger GC
#define allocate_package()  \
  allocate_xrecord(0,Rectype_Package,package_length,0,orecord_type)
# is used by PACKAGE

# UP: allocates Hash-Table
# allocate_hash_table()
# < result: LISP-object Hash-Table
# can trigger GC
#define allocate_hash_table()  \
  allocate_xrecord(0,Rectype_Hashtable,hashtable_length,hashtable_xlength, \
                   orecord_type)
# is used by

# UP: allocates  Readtable
# allocate_readtable()
# < result: LISP-object Readtable
# can trigger GC
#define allocate_readtable()  \
  allocate_xrecord(0,Rectype_Readtable,readtable_length,0,orecord_type)
# is used by IO

# UP: allocates Pathname
# allocate_pathname()
# < result: LISP-object Pathname
# can trigger GC
#define allocate_pathname()  \
  allocate_xrecord(0,Rectype_Pathname,pathname_length,0,orecord_type)
# is used by PATHNAME

#ifdef LOGICAL_PATHNAMES
# UP: allocates Logical Pathname
# allocate_logpathname()
# < result: LISP-object Logical Pathname
# can trigger GC
#define allocate_logpathname()  \
  allocate_xrecord(0,Rectype_Logpathname,logpathname_length,0,orecord_type)
# is used by PATHNAME
#endif

# UP: allocates Random-State
# allocate_random_state()
# < result: LISP-object Random-State
# can trigger GC
#define allocate_random_state()  \
  allocate_xrecord(0,Rectype_Random_State,random_state_length,0,orecord_type)
# is used by IO, LISPARIT

# UP: allocates Byte
# allocate_byte()
# < result: LISP-object Byte
# can trigger GC
#define allocate_byte()  \
  allocate_xrecord(0,Rectype_Byte,byte_length,0,orecord_type)
# is used by LISPARIT

# UP: allocates Fsubr
# allocate_fsubr()
# < result: LISP-object Fsubr
# can trigger GC
#define allocate_fsubr()  \
  allocate_xrecord(0,Rectype_Fsubr,fsubr_length,fsubr_xlength,orecord_type)
# is used by SPVW

# UP: allocates Load-time-Eval
# allocate_loadtimeeval()
# < result: LISP-object Load-time-Eval
# can trigger GC
#define allocate_loadtimeeval()  \
  allocate_xrecord(0,Rectype_Loadtimeeval,loadtimeeval_length,0,orecord_type)
# is used by IO, RECORD

# UP: allocates Symbol-Macro
# allocate_symbolmacro()
# < result: LISP-object Symbol-Macro
# can trigger GC
#define allocate_symbolmacro()  \
  allocate_xrecord(0,Rectype_Symbolmacro,symbolmacro_length,0,orecord_type)
# is used by CONTROL, RECORD

# UP: allocates Global-Symbol-Macro
# allocate_globalsymbolmacro()
# < result: LISP-object Global-Symbol-Macro
# can trigger GC
#define allocate_globalsymbolmacro()  \
  allocate_xrecord(0,Rectype_GlobalSymbolmacro,globalsymbolmacro_length,0,orecord_type)
# is used by RECORD

# UP: allocates a Macro
# allocate_macro()
# < result: a fresh Macro
# can trigger GC
#define allocate_macro()  \
  allocate_xrecord(0,Rectype_Macro,macro_length,0,orecord_type)
# is used by RECORD

# UP: allocates a FunctionMacro
# allocate_functionmacro()
# < result: a fresh FunctionMacro
# can trigger GC
#define allocate_functionmacro()  \
  allocate_xrecord(0,Rectype_FunctionMacro,functionmacro_length,0,orecord_type)
# is used by RECORD

# UP: allocates a BigReadLabel
# allocate_big_read_label()
# < result: a fresh BigReadLabel
# can trigger GC
#define allocate_big_read_label()  \
  allocate_xrecord(0,Rectype_BigReadLabel,bigreadlabel_length,0,orecord_type)
# is used by IO

# UP: allocates an Encoding
# allocate_encoding()
# < result: a fresh Encoding
# can trigger GC
#define allocate_encoding()  \
  allocate_xrecord(0,Rectype_Encoding,encoding_length,encoding_xlength,orecord_type)
# is used by ENCODING

#ifdef FOREIGN
# UP: allocates a foreign-pointer packing
# allocate_fpointer(foreign)
# > foreign: of Type FOREIGN
# < result: LISP-object, contains the foreign pointer
# can trigger GC
  extern maygc object allocate_fpointer (FOREIGN foreign);
/* used by FFI & modules */
#endif
%% #ifdef FOREIGN
%%   puts("extern object allocate_fpointer (FOREIGN foreign);");
%% #endif

# UP: allocates foreign address
# allocate_faddress()
# < result: LISP-object foreign address
# can trigger GC
#define allocate_faddress()  \
  allocate_xrecord(0,Rectype_Faddress,faddress_length,faddress_xlength,orecord_type)
# is used by FOREIGN

# UP: allocates foreign variable
# allocate_fvariable()
# < result: LISP-object foreign variable
# can trigger GC
#define allocate_fvariable()  \
  allocate_xrecord(0,Rectype_Fvariable,fvariable_length,0,orecord_type)
# is used by FOREIGN

# UP: allocates foreign function
# allocate_ffunction()
# < result: LISP-object foreign function
# can trigger GC
#define allocate_ffunction()  \
  allocate_xrecord(0,Rectype_Ffunction,ffunction_length,0,orecord_type)
# is used by FOREIGN

# UP: allocates finalizer
# allocate_finalizer()
# < result: LISP-object finalizer
# can trigger GC
#define allocate_finalizer()  \
  allocate_xrecord(0,Rectype_Finalizer,finalizer_length,0,orecord_type)
# is used by RECORD

# UP: allocates Socket-Server
# allocate_socket_server()
# < result: LISP-object Socket-Server
#ifdef SOCKET_STREAMS
  #define allocate_socket_server() \
    allocate_xrecord(0,Rectype_Socket_Server,socket_server_length,0,orecord_type)
#endif

#ifdef YET_ANOTHER_RECORD
# UP: allocates Yetanother
# allocate_yetanother()
# < result: LISP-object Yetanother
# can trigger GC
  #define allocate_yetanother()  \
    allocate_xrecord(0,Rectype_Yetanother,yetanother_length,0,orecord_type)
# is used by
#endif

# UP: allocates handle
# allocate_handle(handle)
# < result: LISP-object, that contains handle
# can trigger GC
#ifdef FOREIGN_HANDLE
  # can trigger GC
  extern maygc object allocate_handle (Handle handle);
#else
  #define allocate_handle(handle)  fixnum((uintL)(handle))
#endif
%% #if defined(FOREIGN_HANDLE)
%%   puts("extern object allocate_handle (Handle handle);");
%% #else
%%   export_def(allocate_handle(handle));
%% #endif

# UP: allocates Bignum
# allocate_bignum(len,sign)
# > uintC (actually uintWC) len: length of the number (in Digits)
# > sintB sign: flag for sign (0 = +, -1 = -)
# < result: new Bignum (LISP-object)
# can trigger GC
extern maygc object allocate_bignum (uintC len, sintB sign);
# is used by LISPARIT, STREAM

# UP: allocates Single-Float
# allocate_ffloat(value)
# > ffloat value: value (Bit 31 = sign)
# < result: new Single-Float (LISP-object)
# can trigger GC
extern maygc object allocate_ffloat (ffloat value);
# is used by LISPARIT

# UP: allocates Double-Float
#ifdef intQsize
# allocate_dfloat(value)
# > dfloat value: value (Bit 63 = sign)
# < result: new Double-Float (LISP-object)
# can trigger GC
  extern maygc object allocate_dfloat (dfloat value);
#else
# allocate_dfloat(semhi,mlo)
# > semhi,mlo: value (Bit 31 of semhi = sign )
# < result: new Double-Float (LISP-object)
# can trigger GC
  extern maygc object allocate_dfloat (uint32 semhi, uint32 mlo);
#endif
# is used by LISPARIT

# UP: allocates Long-Float
# allocate_lfloat(len,expo,sign)
# > uintC (actually uintWC) len: length of the mantissa (in Digits)
# > uintL expo: exponent
# > signean sign: sign (0 = +, -1 = -)
# < result: new Long-Float, without mantissa
# It will only be a LISP-object when the mantissa has been entered!
# can trigger GC
extern maygc object allocate_lfloat (uintC len, uintL expo, signean sign);
# is used by LISPARIT

# UP: makes a rational number
# make_ratio(num,den)
# > object num: numerator (has to be an integer /= 0, relatively prime to den)
# > object den: denominator (has to be an Integer > 1)
# < result: rational number
# can trigger GC
extern maygc object make_ratio (object num, object den);
# is used by LISPARIT

# UP: makes a complex number
# make_complex(real,imag)
# > real: real part (has to be a real number)
# > imag: imaginary part (has to be a real number /= Fixnum 0)
# < result: complex number
# can trigger GC
extern maygc object make_complex (object real, object imag);
# is used by LISPARIT

/* Adds a freshly allocated object to the list of weak pointers.
 activate_weak(obj);
 > obj: A fresh but filled object of type Rectype_Weak* */
extern void activate_weak (object obj);
# is used by WEAK

# UP: return the length of the ASCIZ-String
# asciz_length(asciz)
# > char* asciz: ASCIZ-String
#       (added with a NULL byte determines the end of string)
# < result: Length of the character sequence (without the NULL byte)
extern uintL asciz_length (const char * asciz);
#if defined(GNU) && (SAFETY < 2)
  #ifdef HAVE_BUILTIN_STRLEN
    #define asciz_length(a)  ((uintL)__builtin_strlen(a))
  #endif
#endif
#ifndef asciz_length
  #ifdef HAVE_SAVED_STACK
    # can not use strlen() instead of asciz_length() , because this would
    # require a begin_system_call()/end_system_call() .
  #else
    # let us presume, that strlen() is implemented efficiently.
    #ifdef STDC_HEADERS
      #include <string.h> # declares strlen()
    #endif
    #define asciz_length(a)  ((uintL)strlen(a))
  #endif
#endif
# is used by SPVW
%% #ifdef asciz_length
%%   export_def(asciz_length(a));
%% #else
%%   puts("extern uintL asciz_length (const char * asciz);");
%% #endif

# UP: Compares two ASCIZ-Strings.
# asciz_equal(asciz1,asciz2)
# > char* asciz1: first ASCIZ-String
# > char* asciz2: second ASCIZ-String
# < result: true if the number-sequences are equal
extern bool asciz_equal (const char * asciz1, const char * asciz2);
# is used by STREAM
%% #if notused
%% #ifdef asciz_length
%%   export_def(asciz_equal(a1,a2));
%% #else
%%   puts("extern bool asciz_equal (const char * asciz1, const char * asciz2);");
%% #endif
%% #endif

/* allocate memory and check for success */
extern void* my_malloc (size_t size);
/* used by FOREIGN and modules */
%% puts("extern void* my_malloc (size_t size);");

/* reallocate memory and check for success */
extern void* my_realloc (void* ptr, size_t size);
/* used by modules */
%% puts("extern void* my_realloc (void *ptr, size_t size);");

# UP: Returns a Table of all circularities within an Object.
# (A circularity is a Sub-Object contained within this Object,
# which has more than one access-path to it.)
# get_circularities(obj,pr_array,pr_closure)
# > object obj: Object
# > bool pr_array: Flag, if Array-Elements recursively count as Sub-Objects
# > bool pr_closure: Flag, if Closure-Components recursively count as Sub-Objects
# < result: T if Stack-Overflow occurred,
#             NIL if no circularities available,
#             #(0 ...) an (n+1)-element Vector, that contains the number 0 and the n
#                      circularities as Elements, n>0.
# can trigger GC
extern maygc object get_circularities (object obj, bool pr_array, bool pr_closure);
# is used by IO

# UP: unentangles #n# - References in Object *ptr with help from Aliste alist.
# > *ptr : Object
# > alist : Alist (Read-Label --> Object, to be substituted)
# < *ptr : Object with unentangled References
# < result : erroneous Reference or nullobj if everything is OK
extern object subst_circ (gcv_object_t* ptr, object alist);
# is used by IO

# UP: Runs through the whole memory, and calls for each
# Object obj: fun(arg,obj,bytelen) .
# map_heap_objects(fun,arg);
# > fun: C-Function
# > arg: arbitrary given Argument
typedef void map_heap_function_t (void* arg, object obj, uintM bytelen);
extern void map_heap_objects (map_heap_function_t* fun, void* arg);
# is used by PREDTYPE

# UP: returns the size (in Bytes) of an object.
# varobject_bytelength(obj)
# > obj: Heap-object with variable length
# < result; the number of bytes occupied by it (header included)
extern uintM varobject_bytelength (object obj);
# is used by PREDTYPE

# Break-Semaphores
# As long as a Break-Semaphore is set, the Lisp-Program can not
# be interrupted. Purpose:
# - backup of Consistencies,
# - Non-reentrant Data-Structures (like e.g. DTA_buffer) can not
#   be used recursively.
typedef union {uintB einzeln[8]; uintL gesamt[2]; } break_sems_;
extern break_sems_ break_sems;
#define break_sem_0  break_sems.einzeln[0]
#define break_sem_1  break_sems.einzeln[1]
#define break_sem_2  break_sems.einzeln[2]
#define break_sem_3  break_sems.einzeln[3]
#define break_sem_4  break_sems.einzeln[4]
#define break_sem_5  break_sems.einzeln[5]
#define break_sem_6  break_sems.einzeln[6]
#define break_sem_7  break_sems.einzeln[7]
# is used by SPVW, Macros set/clr_break_sem_0/1/2/3/4/5/6/7

# Tests whether all break-semaphores have been cleared.
#define break_sems_cleared()  \
  (break_sems.gesamt[0] == 0 && break_sems.gesamt[1] == 0)
# is used by SPVW, WIN32AUX

# clears all break-semaphores. Very dangerous!
#define clear_break_sems()  \
  (break_sems.gesamt[0] = 0, break_sems.gesamt[1] = 0)
# is used by SPVW

# sets break-semaphore 0 and thus protects against interrupts
# set_break_sem_0();
#define set_break_sem_0()  (break_sem_0 = 1)
# is used by SPVW

# clears the break-semaphore 0 and thus releases the interrupts
# clr_break_sem_0();
#define clr_break_sem_0()  (break_sem_0 = 0)
# is used by SPVW

# sets break-semaphore 1 and thus protects against interrupts
# set_break_sem_1();
#define set_break_sem_1()  (break_sem_1 = 1)
# is used by SPVW, ARRAY

# clears the break-semaphore 1 and thus releases the interrupts
# clr_break_sem_1();
#define clr_break_sem_1()  (break_sem_1 = 0)
# is used by SPVW, ARRAY

# sets break-semaphore 2 and thus protects against interrupts
# set_break_sem_2();
#define set_break_sem_2()  (break_sem_2 = 1)
# is used by PACKAGE, HASHTABL

# clears the break-semaphore 2 and thus releases the interrupts
# clr_break_sem_2();
#define clr_break_sem_2()  (break_sem_2 = 0)
# is used by PACKAGE, HASHTABL

# sets break-semaphore 3 and thus protects against interrupts
# set_break_sem_3();
#define set_break_sem_3()  (break_sem_3 = 1)
# is used by PACKAGE

# clears the break-semaphore 3 and thus releases the interrupts
# clr_break_sem_3();
#define clr_break_sem_3()  (break_sem_3 = 0)
# is used by PACKAGE

# sets break-semaphore 4 and thus protects against interrupts
# set_break_sem_4();
#define set_break_sem_4()  (break_sem_4 = 1)
# is used by STREAM, PATHNAME

# clears the break-semaphore 4 and thus releases the interrupts
# clr_break_sem_4();
#define clr_break_sem_4()  (break_sem_4 = 0)
# is used by STREAM, PATHNAME

# increments break-semaphore 5 and thus protects against interrupts
# inc_break_sem_5();
#define inc_break_sem_5()  (break_sem_5++)
# is used by SPVW

# decrements break-semaphore 5 and thus releases interrupts
# dec_break_sem_5();
#define dec_break_sem_5()  (break_sem_5--)
# is used by SPVW

# clears the break-semaphore 5 and thus releases the interrupts
# clr_break_sem_5();
#define clr_break_sem_5()  (break_sem_5 = 0)
# is used by SPVW

# Flag, whether SYS::READ-FORM should behave compatible to ILISP
extern bool ilisp_mode;

# returns the amount of space occupied by static LISP-objects
extern uintM static_space (void);
# is used by DEBUG

# returns the amount of space occupied by LISP-objects
extern uintM used_space (void);
# is used by TIME, DEBUG

# returns the amount of space still available for LISP-objects
extern uintM free_space (void);
# is used by DEBUG

/* UP: saves memory image to disc
 savemem(stream);
 > object stream: open File-Output-Stream, will be closed
 > bool exec_p: should the result include runtime?
 can trigger GC */
extern maygc void savemem (object stream, bool exec_p);
/* used by PATHNAME */

#ifdef HAVE_SIGNALS
# Temporarily do not ignore the status of subprocesses.
  extern void begin_want_sigcld (void);
  extern void end_want_sigcld (void);
# is used by PATHNAME
#endif

#if defined(HAVE_SIGNALS) && defined(SIGPIPE)
  # Set ONLY during write() calls to pipes directed to subprocesses.
extern bool writing_to_subprocess;
#endif


# Declaration of the FSUBRs.
# As C-functions: C_name, of the type fsubr_function_t (no arguments, no value)

# make C-functions visible:
#define LISPSPECFORM  LISPSPECFORM_A
#include "fsubr.c"
#undef LISPSPECFORM
# is used by

# make Fsubr-table visible:
#define LISPSPECFORM  LISPSPECFORM_C
struct fsubr_tab_ {
  #include "fsubr.c"
};
#undef LISPSPECFORM
extern const struct fsubr_tab_ fsubr_tab;
# is used by CONTROL, SPVW


# Declaration of the SUBR-table:
# As C-functions: C_name
# of the type subr_norest_function_t (no arguments, no value)
# resp. subr_rest_function_t (two arguments, no value):
typedef Values subr_norest_function_t (void);
typedef Values subr_rest_function_t (uintC argcount, gcv_object_t* rest_args_pointer);
%% #if notused
%% emit_typedef_f("Values %s(void)","subr_norest_function_t");
%% emit_typedef_f("Values %s(uintC argcount, object* rest_args_pointer)","subr_rest_function_t");
%% #endif

# As LISP-Subr:    L(name)

# Make C-functions visible:
#define LISPFUN  LISPFUN_A
#include "subr.c"
#undef LISPFUN
# is used by

# Make Subr-tables visible:
#define LISPFUN  LISPFUN_C
extern struct subr_tab_ {
  VAROBJECTS_ALIGNMENT_DUMMY_DECL
  #include "subr.c"
} subr_tab_data;
#undef LISPFUN
# is used by Macro L
%% puts("extern struct subr_tab_ {");
%% puts("  VAROBJECTS_ALIGNMENT_DUMMY_DECL");
%% #undef LISPFUN
%% #define LISPFUN(name,sec,req_anz,opt_anz,rest_flag,key_flag,key_anz,keywords) \
%%   printf("  subr_t %s;\n",STRING(D_##name));
%% #include "subr.c"
%% #undef LISPFUN
%% puts("} subr_tab_data;");

# Abbreviation for LISP-Subr with a given name: L(name)
#if !defined(MAP_MEMORY_TABLES)
  #define subr_tab  subr_tab_data
  #ifdef TYPECODES
    #define subr_tab_ptr_as_object(subr_addr)  (type_constpointer_object(subr_type,subr_addr))
  #else
    #if defined(WIDE_AUXI)
      #define subr_tab_ptr_as_object(subr_addr)  as_object_with_auxi((aint)(subr_addr)+subr_bias)
    #elif defined(OBJECT_STRUCT)
      #define subr_tab_ptr_as_object(subr_addr)  as_object((oint)(subr_addr)+subr_bias)
    #else
      #define subr_tab_ptr_as_object(subr_addr)  objectplus(subr_addr,subr_bias)
    #endif
  #endif
  #define L_help_(name)  subr_tab_ptr_as_object(&subr_tab.name)
#else
  # define subr_tab_addr  ((struct subr_tab_ *)type_constpointer_object(subr_type,0))
  #define subr_tab_addr  ((struct subr_tab_ *)type_zero_oint(subr_type))
  #define subr_tab  (*subr_tab_addr)
  #define subr_tab_ptr_as_object(subr_addr)  (as_object((oint)(subr_addr)))
  #define L_help_(name)  subr_tab_ptr_as_object(&subr_tab_addr->name)
#endif
#define L(name)  L_help_(D_##name)
# is used by all modules
%% #if defined(MAP_MEMORY_TABLES)
%%   export_def(subr_tab_addr);
%% #endif
%% export_def(subr_tab);
%% export_def(subr_tab_ptr_as_object(subr_addr));
%% export_def(L_help_(name));
%% emit_define("L(name)","L_help_(D_##name)");


# Pseudofunctions are addresses of C functions (to be called directly, not via
# FUNCALL) or constant C data.
# For SAVEMEM/LOADMEM we have a table of all such pseudofunctions.
typedef const void *  Pseudofun; # assume function pointers fit in a void*
%% puts("typedef const void *  Pseudofun;");

# Declaration of the tables of relocatable pointers:
#define PSEUDO  PSEUDO_A
extern struct pseudocode_tab_ {
  #include "pseudofun.c"
} pseudocode_tab;
#undef PSEUDO
#define PSEUDO  PSEUDO_B
extern struct pseudodata_tab_ {
  #include "pseudofun.c"
  #if defined(MICROSOFT) && !defined(UNICODE)
  Pseudofun dummy_pseudofun;
  #endif
} pseudodata_tab;
#undef PSEUDO
# is used by STREAM, SPVW

# Declaration of the functions that can be stored in Lisp objects.
#define PSEUDO  PSEUDO_C
#include "pseudofun.c"
#undef PSEUDO
# is used by STREAM, and to avoid gcc -Wmissing-declarations warnings

# Return an ADDRESS object encapsulating a pseudofunction.
#ifdef TYPECODES
  #define P(fun)  type_constpointer_object(machine_type,(Pseudofun)&(fun))
#else
  #define P(fun)  make_machine_code((Pseudofun)&(fun))
#endif
# is used by STREAM, ENCODING


# Declaration if the Symbol-table:
#define LISPSYM  LISPSYM_A
extern struct symbol_tab_ {
  VAROBJECTS_ALIGNMENT_DUMMY_DECL
  #include "constsym.c"
} symbol_tab_data;
#undef LISPSYM
# is used by Macro S, gcinvariant_symbol_p
%% puts("extern struct symbol_tab_ {");
%% puts("  VAROBJECTS_ALIGNMENT_DUMMY_DECL");
%% #define LISPSYM(name,printname,package)  \
%%   printf("  symbol_ %s;\n",STRING(S_##name));
%% #include "constsym.c"
%% #undef LISPSYM
%% puts("} symbol_tab_data;");

# Abbreviation for LISP-Symbol with a given name: S(name)
#define S(name)  S_help_(S_##name)
#if !defined(MAP_MEMORY_TABLES)
  #define symbol_tab  symbol_tab_data
  #ifdef TYPECODES
    #define S_help_(name)  (type_constpointer_object(symbol_type,&symbol_tab.name))
  #else
    #if defined(WIDE_AUXI)
      #define S_help_(name)  as_object_with_auxi((aint)&symbol_tab.name+varobject_bias)
    #elif defined(OBJECT_STRUCT)
      #define S_help_(name)  as_object((oint)&symbol_tab.name+varobject_bias)
    #else
      #define S_help_(name)  objectplus(&symbol_tab.name,varobject_bias)
    #endif
  #endif
#else
  # define symbol_tab_addr ((struct symbol_tab_ *)type_constpointer_object(symbol_type,0))
  #define symbol_tab_addr ((struct symbol_tab_ *)type_zero_oint(symbol_type))
  #define symbol_tab  (*symbol_tab_addr)
  #define S_help_(name)  (as_object((oint)(&symbol_tab_addr->name)))
  #if 0 # Some compilers do not allow the above expression
        # - even though it's a 'constant expression' -
        # as initializer of static variables.
        # We have to assist:
    #undef S_help_
    #define S_help_(name)  (as_object( (char*)(&((struct symbol_tab_ *)0)->name) + (uintP)symbol_tab_addr ))
  #endif
#endif
# is used by all modules
%% emit_define("S(name)","S_help_(S_##name)");
%% #if defined(MAP_MEMORY_TABLES)
%%   export_def(symbol_tab_addr);
%% #endif
%% export_def(symbol_tab);
%% export_def(S_help_(name));

#define NIL  S(nil)
#define T    S(t)
%% export_def(NIL);
%% export_def(T);

#if defined(DEBUG_GCSAFETY)
# gcinvariant_symbol_p(obj)
# > obj: an object
# < result: true if obj is a symbol in symbol_tab
static inline bool gcinvariant_symbol_p (object obj) {
  if (
      #ifdef TYPECODES
        symbolp(obj)
      #else
        varobjectp(obj)
      #endif
      &&
      (
       #if !defined(MAP_MEMORY_TABLES)
         #ifdef TYPECODES
           (as_oint(obj) >> (oint_addr_shift-addr_shift)) - ((aint)(tint)symbol_type<<oint_type_shift)
         #else
           as_oint(obj) - varobject_bias
         #endif
       #else
         # FIXME: MAP_MEMORY_TABLES possibly uses MULTIMAP_MEMORY_SYMBOL_TAB.
         as_oint(obj)
       #endif
       - (aint)&symbol_tab < sizeof(symbol_tab))
     )
    return true;
  else
    return false;
}
#endif
%% #if defined(DEBUG_GCSAFETY)
%%   printf("static inline bool gcinvariant_symbol_p (object obj) { if (");
%%   #ifdef TYPECODES
%%     printf("symbolp(obj)");
%%   #else
%%     printf("varobjectp(obj)");
%%   #endif
%%   printf(" && (");
%%   #if !defined(MAP_MEMORY_TABLES)
%%     #ifdef TYPECODES
%%       printf2("(as_oint(obj) >> %d) - %d", oint_addr_shift-addr_shift, (aint)(tint)symbol_type<<oint_type_shift);
%%     #else
%%       printf1("as_oint(obj) - %d", varobject_bias);
%%     #endif
%%   #else
%%     printf("as_oint(obj)");
%%   #endif
%%   puts(" - (aint)&symbol_tab < sizeof(symbol_tab))) return true; else return false; }");
%% #endif

# The macro NIL_IS_CONSTANT tells , whether NIL is recognized
# as 'constant expression' by the C-Compiler. If so, tables can
# already be initialized largely by the C-Compiler.
#if (oint_addr_shift==0)
  #define NIL_IS_CONSTANT  true
#else
  #define NIL_IS_CONSTANT  false
#endif

# Declaration of the table with the remaining constant objects:
#define LISPOBJ  LISPOBJ_A
extern struct object_tab_ {
  #include "constobj.c"
} object_tab;
#undef LISPOBJ
# is used by Macro O
%% puts("extern struct object_tab_ {");
%% #define LISPOBJ(name,init)  printf("  gcv_object_t %s;\n",STRING(name));
%% #include "constobj.c"
%% #undef LISPOBJ
%% puts("} object_tab;");

# Abbreviation for other LISP-object with a given Name:
#define O(name)  (object_tab.name)
%% /* FIXME: Difference between lispbibl.d and clisp.h */
%% puts("#define GLO(name)  (object_tab.name)");

#if defined(GENERATIONAL_GC) && defined(SPVW_MIXED)
# handle_fault_range(PROT_READ,start,end) makes an address range readable.
# handle_fault_range(PROT_READ_WRITE,start,end) makes an address range writable.
  extern bool handle_fault_range (int prot, aint start_address, aint end_address);
#endif
%% export_def(PROT_READ);
%% export_def(PROT_READ_WRITE);
%% #if defined(GENERATIONAL_GC) && defined(SPVW_MIXED)
%% puts("extern bool handle_fault_range (int prot, aint start_address, aint end_address);");
%% #else
%% puts("#define handle_fault_range(p,s,e)");
%% #endif


# ###################### MODBIBL for MODULES.D ############################ #

#if defined(DYNAMIC_MODULES) && !defined(HAVE_DYNLOAD) && !defined(WIN32_NATIVE)
  /* if you want DYNAMIC_MODULES to work on a non-WIN32_NATIVE platform
     which does not HAVE_DYNLOAD (e.g., via ltdl), you will need to
     implement libopen() and find_name() in spvw.d for your platform */
  #error "Dynamic modules require dynamic loading!"
#endif

# Number of external modules:
extern uintC module_count;
%% puts("extern uintC module_count;");

# Data for initialization of a module's subr_tab:
typedef struct {
  const char* packname; # Name of the Home-Package of the Symbol or NULL
  const char* symname; # Name of the Symbol
} subr_initdata_t;
%% emit_typedef("struct { const char* packname; const char* symname; }","subr_initdata_t");

# Data for initialization of a module's object_tab:
typedef struct {
  const char* initstring; # Initialization-String
} object_initdata_t;
%% emit_typedef("struct { const char* initstring; }","object_initdata_t");

# Table resp. List of Modules:
typedef struct module_t {
  const char* name; # Name
  subr_t* stab; const uintC* stab_size; # a separate subr_tab
  gcv_object_t* otab; const uintC* otab_size; # a separate object_tab
  bool initialized;
  # Data for Initialization:
  const subr_initdata_t* stab_initdata;
  const object_initdata_t* otab_initdata;
  # Functions for Initialization
  void (*initfunction1) (struct module_t *); # only once
  void (*initfunction2) (struct module_t *); # always at start up
  void (*finifunction) (struct module_t *);  /* before termination */
  #ifdef DYNAMIC_MODULES
    struct module_t * next; # linked List
  #endif
} module_t;
#ifdef DYNAMIC_MODULES
  extern module_t modules[]; # List-Start
  BEGIN_DECLS
  extern void add_module (module_t * new_module);
  END_DECLS
#else
  extern module_t modules[]; # 1+module_count entries, then an empty entry
#endif
%% strcpy(buf,"struct module_t { const char* name; subr_t* stab; const uintC* stab_size; gcv_object_t* otab; const uintC* otab_size; bool initialized; const subr_initdata_t* stab_initdata; const object_initdata_t* otab_initdata; void (*initfunction1) (struct module_t *); void (*initfunction2) (struct module_t *); void (*finifunction) (struct module_t *);");
%% #ifdef DYNAMIC_MODULES
%%   strcat(buf," struct module_t * next;");
%% #endif
%% strcat(buf," }"); emit_typedef(buf,"module_t");
%% #ifdef DYNAMIC_MODULES
%%   puts("BEGIN_DECLS");
%%   puts("extern void add_module (module_t * new_module);");
%%   puts("END_DECLS");
%% #else
%%   puts("extern module_t modules[];");
%% #endif

#if defined(HAVE_DYNLOAD) || defined(WIN32_NATIVE)
/* open the dynamic library
 libname is the name of the library
 returns a handle suitable for find_name()
 calls dlopen() or LoadLibrary() */
extern void * libopen (const char* libname);
/* used by FOREIGN and spvw.d:dynload_modules() */

/* find the name in the dynamic library handle
 calls dlsym() or GetProcAddress()
 handle is an object returned by libopen()
        or NULL, which means emulate RTLD_DEFAULT on UNIX_FREEBSD
        and WIN32_NATIVE by searching through all libraries
 name is the name of the function (or variable) in the library */
extern void* find_name (void *handle, const char *name);
/* used by FOREIGN and spvw.d:dynload_modules() */
#endif

#if defined(DYNAMIC_MODULES)
/* Attaches a shared library to this process' memory, and attempts to load
   a number of clisp modules from it. */
extern void dynload_modules (const char * library, uintC modcount,
                             const char * const * modnames);
#endif

/* find the module with the given name */
extern module_t* find_module (const char *name);
/* push all module names to STACK and return the number of modules pushed
 can trigger GC */
extern maygc uintC modules_names_to_stack (void);

# ####################### EVALBIBL for EVAL.D ############################## #

/*

Specifications for the Evaluator
################################

SUBRs and FSUBRs
================

They're constructed through
  LISPFUN             for general LISP-functions,
  LISPFUNN            for normal  LISP-functions (only required parameters),
  LISPSPECFORM        for special forms (FSUBRs).
Note that SUBRs with KEY_ANZ=0 will be seen as SUBRs without keyword-
parameters by the evaluator (which in consequence means that in this case the
ALLOW_FLAG is meaningless and no keyword, not even :ALLOW-OTHER-KEYS,
will be accepted)!

Values
======

The following format is used for the passing of multiple values:
value1 contains the first value (NIL if there aren't values).
mv_count contains the number of values..
If there is at least one value       : value1 = first value.
If there are at least two values     : value2 = second value.
If there are at least three values   : value3 = third value .
All values are in mv_space .
Recommended commands for returning of values to the caller:
  0 values:   VALUES0;
  1 value :     VALUES1(...);
  2 values:   value1=...; value2=...; mv_count=2;
  3 values:   value1=...; value2=...; value3=...; mv_count=3;
  more than 3 values:
              if (number of values >= mv_limit) goto error_too_many_values;
              Put the values one after another onto the STACK
              STACK_to_mv(number of values);

Passing of parameters to SUBRs
==============================

The arguments are passed on the LISP-stack, with the first one being on the
top. The required arguments come first, then the optional ones
(each #UNBOUND, if not specified), then come the
keyword-arguments (again, each #UNBOUND, if not specified).
The SUBR-object can be found in back_trace.
This is all if no &REST-argument is planned. But if a &REST-argument
is planned, all further arguments follow (the optional ones) on the stack
one by one, and this will be passed: the number of these arguments and a pointer
above the first of these arguments. (This means that the number of LISP-objects on
the stack is not always the same!)
All arguments have to be removed from the LISP-stack at the return jump.
(for example. for SUBRs with &REST: the stackpointer STACK has to have the value
args_pointer = rest_args_pointer STACKop (fixed number of arguments)
= pointer above the very first argument), and mv_count/mv_space
has to hold the values.

Passing of parameters to FSUBRs
===============================

The parameters are passed on the LISP-stack with the first one being on top.
At first there are the required parametes, followed by the optional ones
(#UNBOUND, if not specifired), then - if body-flag true -
the whole rest of the body (most of the time a list).
So the number of objects on the LISP-stack is always the same, namely
numReqParameter + numOptParameter + (0 or 1 if body-flag).
At the call, back_trace holds the FSUBR-object, and the whole form is
in the EVAL-frame, directly above the parameters.
All parameters have to be removed from the LISP-stack at the return jump
(ie. the stackpointer STACK has to be incemented by the number of objects),
and mv_count/mv_space has to hold the values.

Environments
============

General
-------
The lexical environment is separated into 5 components:
  - the variables-environment (VAR_ENV),
  - the functions- and macro-environment (FUN_ENV),
  - the block-environment (BLOCK_ENV),
  - the tagbody-environment (GO_ENV),
  - the declarations-environment (DECL_ENV).
The environment is kept in 5 "global variables". They are dynamically bound
with special frames on change.
A single functions- and macro environment is passed to SYM_FUNCTION,
MACROEXP, MACROEXP0, PARSE_DD.
GET_CLOSURE expects a pointer to all environments en bloc: A3 with
VAR_(A3)=VAR_ENV, FUN_(A3)=FUN_ENV, BLOCK_(A3)=BLOCK_ENV, GO_(A3)=GO_ENV,
DECL_(A3)=DECL_ENV.

The variables-environment
-------------------------
It contains the local variable-bindings.
A variables-enviroment is given through a pointer to a
variable-binding frame, or NIL  (which means an empty lexical
environment) or a vector that is built as follows:
The vector contains n bindings and has the length 2n+1. The elements are
n-times each variable (a symbol) and  the value that belongs to it ("value" can
be #<SPECDECL> as well, and then the variable has to be referenced dynamically)
and as last element the predecessor environment.

The functions- and macro-environment
------------------------------------
It contains the local function- and macro-definitions.
A functions- and macro-environment is given through a pointer to
a functions- or macrobindings-frame or NIL (which means an empty
lexical environment) or through a vector that is built as follows:
The vector contains n bindings and has length 2n+1. The elements are
n-time each function-name (a symbol) and the definiton that belongs to it (a
closure or NIL or a SYS::MACRO object) and as last element
the predecessor environment.

The block-environment
---------------------
It contains the lexically visible block-exitpoints.
A block-environment is given through a pointer to a block-frame
or through an association-list, whose elements each have the block-name (a symbol)
as CAR and as CDR either the pointer to the appropriate
frame or #DISABLED, if the block has already been left.

The tagbody-environment
-----------------------
It contains the lexically visible Go-labels of the tagbodies.
A tagbody-environment is given through a pointer to a
tagbody-frame or an associations-list, whose elements have a vector (with the
Go-tags as elements) as CAR and as CDR either the pointer to the
related frame or #DISABLED, if the tagbody has already
been left.

The declarations-environment
----------------------------
It contains the lexically visible declarations.
A declarations-environment is given through a list of declaration-
specifiers, whose CAR is each either OPTIMIZE or DECLARATION or
a user-specified declaration-type.

Passing of environtments to LISP-functions
------------------------------------------
There are two data structures for this:
When it is passed as second argument to macro-expander-functions (CLTL p.
145-146) and when it is receipted by MACROEXPAND and MACROEXPAND-1 (CLTL p. 151)
it is simply a Simple-Vector with 2 elements, consisting of a nested
variable-environment and a nested functions- and macro-environment.
The same for passing to  SYSTEM::%EXPAND-LAMBDABODY-MAIN and the like.
If it is passed as second argument to the value of *EVALHOOK* or as third one
to the value of *APPLYHOOK* (CLTL p. 322) and on reception by
EVALHOOK and APPLYHOOK (CLTL p. 323) it is a Simple-Vector with
five elements with all five components nested.

Frames
======
Frames are not used to call SUBRs, FSUBRs and compiled closures.

There are the following 14 kinds of frames:
  - Environmentbinding-Frame (ENV_FRAME),
  - APPLY-frame (APPLY_FRAME),
  - EVAL-frame (EVAL_FRAME),
  - dynamic variable-bindings-frame (DYNBIND_FRAME),
  - Variable-bindings-frame (VAR_FRAME),
  - Function- or Macrobindings-Frame (FUN_FRAME),
  - interpreted block-frame (IBLOCK_FRAME),
  - compiled block-frame (CBLOCK_FRAME),
  - interpreted tagbody-frame (ITAGBODY_FRAME),
  - compiled tagbody-frame (CTAGBODY_FRAME),
  - Catch-Frame (CATCH_FRAME),
  - Unwind-Protect-frame (UNWIND_PROTECT_FRAME),
  - Handler-frame (HANDLER_FRAME),
  - Driver-frame (DRIVER_FRAME).
Right at the bottom of a frame there is a long-word, that contains the
frame-type information and a pointer above the frame (= the value of the
STACK before and after the frame has been built).
In the frame-info there are the bits
  SKIP2_BIT      deleted, if another long-word comes above it,
                   that is not a LISP-object and thus has to be skipped
                   by the GC,
  EXITPOINT_BIT  set for all but VAR and FUN,
  NESTED_BIT     set for IBLOCK and ITAGBODY, if the exitpoint or
                   the Go-label has already been put into an Alist.
The default-values for  the frame-type info-bytes are ENVxx_FRAME_INFO,
APPLY_FRAME_INFO, EVAL_FRAME_INFO, VAR_FRAME_INFO, FUN_FRAME_INFO,
IBLOCK_FRAME_INFO, CBLOCK_FRAME_INFO, ITAGBODY_FRAME_INFO, CTAGBODY_FRAME_INFO,
CATCH_FRAME_INFO, UNWIND_PROTECT_FRAME_INFO, DRIVER_FRAME_INFO.
The routine that is in (SP).L with SP=SP_(STACK) (for IBLOCK-, CBLOCK-,
ITAGBODY-, CTAGBODY-, CATCH-, UNWIND-PROTECT-frames), is being
jumped to by MOVE.L SP_(STACK),SP ! RTS  .
For DRIVER-frames by MOVE.L SP_(STACK),SP ! MOVE.L (SP),-(SP) ! RTS  .
In the portable C-version in SP_(STACK) there is a pointer to a
setjmp/longjmp-buffer.

Environmentbindings-frames
--------------------------
They contain dynamic bindings of a maximum of 5 environments.
ENVxx_FRAME_INFO  is frame-info (xx depending on the environment that is
bound here). Structure:
    Offset        Stack-Contents
  20/16/12/8/4  [old value ofDECL_ENV]
  16/12/8/4     [old value ofGO_ENV]
  12/8/4        [old value ofBLOCK_ENV]
  8/4           [old value ofFUN_ENV]
  4             [old value ofVAR_ENV]
  0             Frame-Info; pointer above frame

ENV1V_frame    for 1 VAR_ENV
ENV1F_frame    for 1 FUN_ENV
ENV1B_frame    for 1 BLOCK_ENV
ENV1G_frame    for 1 GO_ENV
ENV1D_frame    for 1 DECL_ENV
ENV2VD_frame   for 1 VAR_ENV and 1 DECL_ENV
ENV5_frame     for all 5 environments

APPLY-frames
------------
They are created at every call (APPLY or FUNCALL) of an interpreted
closure.
Structure:
  Offset     Stack-contents
  4n+12
  4n+8      Argument 1
  ...
  12        Argument n
  8         Function that is being called
  4         SP
  0         Frame-info; pointer above frame
SP is a pointer into the program-stack. Jumping back to (SP).L after dissolving
the APPLY-fame returns the contents of A0/... as values of the form.
The frame-info has the value APPLY_FRAME_INFO or TRAPPED_APPLY_FRAME_INFO.

EVAL-frames
-----------
They are created for every call of the EVAL-procedure.
Layout:
  Offset     Stack-content
  8         Form that is being evaluated
  4         SP
  0         Frame-info; pointer above frame
SP is a pointer into the program stack. Jumping back to (SP).L after dissolving
the EVAL-frame returns the contents of A0/... as values of the form.
The frame-info has the value EVAL_FRAME_INFO or TRAPPED_EVAL_FRAME_INFO.

Dynamic variable-bindings frames
-----------------------------------
They bind symbols to values dynamically.
The structure of such a frame with n bindings is as follows::
  Offset  stack contents
  8n+4
  8n      value 1
  8n-4    symbol 1
  ...     ...
  8       value n
  4       symbol n
  0       frame-info; pointer above frame
The content of the frameinfo-byte is DYNBIND_FRAME_INFO.

Variable-bindings-frames
------------------------
They are created when interpreted closures are being used (for the variable
bindings specified in the Lambda-list and in the dynamic references that might
be specified in the declarations) and by LET and LET*, as well as by all
constructs, that use LET or LET* implicitly (such as DO, DO*, PROG, PROG*,
DOLIST, DOTIMES, ...).
The structure of a variable-bindings-frame with n bindings is as follows:
#ifndef NO_symbolflags
  Offset  stack contents
  12+8n
  8+8n    value 1
  4+8n    symbol 1
  ...     ...
  16      value n
  12      symbol n
  8       NEXT_ENV
  4       m
  0       frame-info; pointer above frame
#else
  Offset  stack contents
  12+12n
  8+12n   value 1
  4+12n   symbol 1
  12n     marker bits 1
  ...     ...
  20      value n
  16      symbol n
  12      marker bits n
  8       NEXT_ENV
  4       m
  0       frame-info; pointer above frame
#endif
The symbol/value-pairs are numbered and stored in the order in which the
bindings become active (i.e. for interpreted closures: at first the dynamic
references (SPECIAL-declarations), then the required-parameters, then the
optional parameters, then the remaining parameters, then the keyword
parameters, then the AUX-variables).
The symbols contain the following marker bits on the stack: ACTIVE_BIT, is
set, if the binding is active, DYNAM_BIT is set, if the binding is
dynamic. (Dynamic references are marked as lexical with
the special value #SPECDECL!).
NEXT_ENV is next upper variables-environment.
m is a long-word, 0 <= m <= n, and stands for the number of bindings that
have not yet been put into a vector by NEST-operations. Thus
the symbol/value-pairs 1,...,n-m have been active but been nested meanwhile
and thus inactive again on the stack (if the bindings were static).
Only some of the pairs n-m+1,...,n can be static and active.
The frameinfo-byte contains VAR_FRAME_INFO.

Function- and Macrobindings-Frames
-----------------------------------
They are created by FLET and MACROLET.
The structure of a variable-bindings-frame with n bindings is as follows:
  Offset  stack contents
  12+8n
  8+8n    value 1
  4+8n    symbol 1
  ...     ...
  16      value n
  12      symbol n
  8       NEXT_ENV
  4       m
  0       Frame-Info; pointer above frame
NEXT_ENV is the next higher function-environment.
m is a long word, 0 <= m <= n, and stands for the number of bindings, that
have not yet been put into a vector by NEST-operations. So the
symbol/value pais 1,...,n-m have been active, but nested meanwhile and thus
inactive on the stack again. Only the pairs n-m+1,...,n are active.
Marker bits are not needed here, as opposed to the variable-bindings frames

All values are closures or SYS::MACRO objects.
The content of the Frameinfo-bytes is FUN_FRAME_INFO.

Interpreted Block-Frames
------------------------
They are created by BLOCK and all constructs that contain an implicit
BLOCK (e.g. DO, DO*, LOOP, PROG, PROG*, ...). The structure is as follows:
  Offset  stack contents
  16
  12       NAME
  8        NEXT_ENV
  4        SP
  0        Frame-Info; pointer above frame
NAME is the name of the block. NEXT_ENV is the next higher Block-Environment.
SP is a pointer into the program stack, (SP).L is a routine, that unwinds the
Block-Frame and leaves the block with the values A0-A2/...
Frame-Info is IBLOCK_FRAME_INFO, possibly with set NESTED_BIT (then NEXT_ENV
points to an Alist, whose first element is the pair (NAME . <Framepointer>),
because the block is not DISABLED yet).

Compiled Block-Frames
---------------------
Structure:
  Offset  stack contents
   12
   8        Cons (NAME . <Framepointer>)
   4        SP
   0        Frame-Info; pointer above frame
NAME is the name of the block.
SP is a pointer into the program stack, (SP).L is a routine, that
unwinds the Block-Frame and leaves the block with the values A0-A2/...
Frame-Info is CBLOCK_FRAME_INFO.

Interpreted Tagbody-Frames
--------------------------
They are created by TAGBODY and all constructs that contain an implicit
TAGBODY (e.g. DO, DO*, PROG, PROG*, ...).
The structure of a Tagbody-Frames with n tags is as follows:
  Offset  stack contents
  12+8n
  8+8n     BODY 1
  4+8n     TAG 1
  ...      ...
  16       BODY n
  12       TAG n
  8        NEXT_ENV
  4        SP
  0        Frame-Info; pointer above frame
The tags are the jump destinations ; they are symbols and Integers, that are in
the Body. The corresponding "value" BODY i contains the part of the body
that follows TAG i. NEXT_ENV is the next higher Tagbody-Environment.
SP is a pointer into the program stack, (SP).L is a routine, that executes
the action (GO TAGi), if it is jumped to with BODYi in A0.
Frame-Info is ITAGBODY_FRAME_INFO, poss. with set NESTED_BIT (then
NEXT_ENV points to an Alist, whose first element has the form
(#(TAG1 ... TAGn) . <Framepointer>), because the Tagbody is not
DISABLED yet).

Compiled Tagbody-Frames
-----------------------
Structure:
  Offset  stack contents
   12
   8        Cons (#(TAG1 ... TAGn) . <Framepointer>)
   4        SP
   0        Frame-Info; above frame
TAG1, ..., TAGn are the names of the tags (actually only contained in
the compiled code to create error messages).
SP is a pointer into the program stack, (SP).L is a routine, that executes
the action (GO TAGi), if it has been jumped at with value1 = i (1 <= i <= n)

Frame-Info is CTAGBODY_FRAME_INFO.

Catch-Frames
------------
They are created by the  Special-Form CATCH. Its structure is as follows:
  Offset  stack contents
   12
   8        TAG
   4        SP
   0        Frame-Info; pointer above frame
TAG is the tag of the catcher.
SP is a pointer into the program stack, (SP).L is a routine, that unwinds
the Frame and returns the values A0-A2/...
Frame-Info is CATCH_FRAME_INFO.

Unwind-Protect-Frames
---------------------
They are created by the Special-Form UNWIND-PROTECT and all constructs
that contain an implicit UNWIND-PROTECT (like WITH-OPEN-STREAM or
WITH-OPEN-FILE). Their structure is as follows:
  Offset  Stack-contents
   8
   4        SP
   0        Frame-Info; pointer above frame
SP is a pointer into the program stack. (SP).L a routine, that unwinds the
Frame,saves the current values A0-A2/...  executes the cleanup,
writes the saved values back and finally jumps to the address
(with RTS), that has been entered into the program stack in place of their own
and leaves D6 unchanged.

Handler-Frames
--------------
They are created by the macro HANDLER-BIND. Their structure is as follows:
  Offset  Stack-contens
   16
   12       Cons (#(type1 label1 ... typem labelm) . SPdepth)
   8        Closure
   4        SP
   0        Frame-Info; pointer above frame
SP is a pointer into the program stack.
If there is a condition of the type typei
the closure starting at Bte labeli is interpreted as Handler, where at first
a piece of the program stack with the length SPdepth  is duplicated.
One variant of Handler-Frames calls a C-Handler:
  Offset  Stack-contents
   16
   12       Cons (#(type1 label1 ... typem labelm))
   8        Handler-function
   4        SP
   0        Frame-Info; pointer above frame
SP is a pointer into the program stack.
If there is a condition of the type typei
the handler-function is called with the arguments SP
(arbitrary pointer into the C-Stack), frame (pointer above the frame),
labeli (arbitrary Lisp-object), condition.
If the Handler wants to yield control via unwind_upto(FRAME) by itself,
the Frame has to be created with finish_entry_frame.

Driver-Frames
-------------
They are created upon entry into a top-level loop (most of the time
a READ-EVAL-PRINT-loop) and are used to continue the previous top-level
loop after an error message. The structure is simple
  Offset  Stack-contens
   8
   4        SP
   0        Frame-Info; pointer above Frame
SP is a pointer into the program stack. (SP).L is a routine, that
re-enters the corresponding top-level loop.

*/

# STACK:
# STACK is the LISP-Stack.
# STACK_0 is the first object on the STACK.
# STACK_1 is the second object on the STACK.
# etc., generally STACK_(n) = (n+1)th object on the STACK.
# pushSTACK(object)  puts an object onto the Stack. Synonym: -(STACK).
# popSTACK()  returns STACK_0 and removes it from the stack.
# skipSTACK(n);  removes n objects from the STACK.
# If you want to save the value of the stack, you do this:
#   var gcv_object_t* temp = STACK; ... (no access through temp !) ... setSTACK(STACK = temp);
#   but: access through STACKpointable(temp)  is possible.
# If you want a pointer that can traverse through the Stack, you do this:
#   var gcv_object_t* ptr = &STACK_0;  or = STACKpointable(STACK);
#   assert( *(ptr STACKop 0) == STACK_0 );
#   assert( *(ptr STACKop 1) == STACK_1 );
#   ...
#   ptr skipSTACKop n;
#   assert( *(ptr STACKop 0) == STACK_(n) );
#   ...
#   This pointer must not be assigned to the STACK again!
# If you store blocks of objects on the STACK and want to get the (n+1)-th block,
#   you do this:  STACKblock_(type,n). type should be a
#   struct-type with sizeof(type) a multiple of sizeof(gcv_object_t).

#ifdef STACK_DOWN
  #define STACK_(n)  (STACK[(sintP)(n)])
  #define STACKpointable(STACKvar)  ((gcv_object_t*)(STACKvar))
  #define skipSTACKop  +=
  #define STACKop      +
  #define cmpSTACKop   <
  #define STACKblock_(type,n)  (((type*)STACK)[(sintP)(n)])
#endif
#ifdef STACK_UP
  #define STACK_(n)  (STACK[-1-(sintP)(n)])
  #define STACKpointable(STACKvar)  ((gcv_object_t*)(STACKvar)-1)
  #define skipSTACKop  -=
  #define STACKop      -
  #define cmpSTACKop   >
  #define STACKblock_(type,n)  (((type*)STACK)[-1-(sintP)(n)])
#endif
#define pushSTACK(obj)  (STACK_(-1) = (obj), STACK skipSTACKop -1)
  # Almost equivalent with *--STACK = obj  resp.  *STACK++ = obj  , but
  # Careful: first enter the object into STACK_(-1), THEN modify the STACK!
#define popSTACK()  (STACK skipSTACKop 1, STACK_(-1))
#define skipSTACK(n)  (STACK skipSTACKop (sintP)(n))

#if defined(GNU) && defined(MC680X0) && !defined(NO_ASM) && !defined(WIDE) && defined(STACK_register)
  # With GNU and a 680X0 STACK is in a register. Access and
  # modification of the STACK are an atomic unit that cannot be interrupted.
  #undef pushSTACK
  #undef popSTACK
  #ifdef STACK_DOWN
    # define pushSTACK(obj)  (*--STACK = (obj))
    #define pushSTACK(obj)  \
      ({ __asm__ __volatile__ ("movel %0,"REGISTER_PREFIX""STACK_register"@-" : : "g" ((object)(obj)) : STACK_register ); })
    # define popSTACK()  (*STACK++)
    #define popSTACK()  \
      ({var object __result;                                                                                         \
        __asm__ __volatile__ ("movel "REGISTER_PREFIX""STACK_register"@+,%0" : "=g" (__result) : : STACK_register ); \
        __result;                                                                                                    \
       })
  #endif
  #ifdef STACK_UP
    # define pushSTACK(obj)  (*STACK++ = (obj))
    #define pushSTACK(obj)  \
      ({ __asm__ __volatile__ ("movel %0,"REGISTER_PREFIX""STACK_register"@+" : : "g" ((object)(obj)) : STACK_register ); })
    # define popSTACK()  (*--STACK)
    #define popSTACK()  \
      ({var object __result;                                                                                         \
        __asm__ __volatile__ ("movel "REGISTER_PREFIX""STACK_register"@-,%0" : "=g" (__result) : : STACK_register ); \
        __result;                                                                                                    \
       })
  #endif
#endif
#if defined(SPARC) && !defined(GNU) && !defined(__SUNPRO_C) && !defined(MULTITHREAD) && (SAFETY < 2)
  #undef pushSTACK
  #undef popSTACK
  #undef skipSTACK
  #define pushSTACK(obj)  (STACK_(-1) = (obj), _setSTACK(STACK STACKop -1))
  #define popSTACK()  (_setSTACK(STACK STACKop 1), STACK_(-1))
  #define skipSTACK(n)  (_setSTACK(STACK STACKop (sintP)(n)))
#endif
%% export_def(STACK_(n));
%% export_def(skipSTACKop);
%% export_def(STACKop);
%% export_def(pushSTACK(obj));
%% export_def(popSTACK());
%% export_def(skipSTACK(n));

#define STACK_0  (STACK_(0))
#define STACK_1  (STACK_(1))
#define STACK_2  (STACK_(2))
#define STACK_3  (STACK_(3))
#define STACK_4  (STACK_(4))
#define STACK_5  (STACK_(5))
#define STACK_6  (STACK_(6))
#define STACK_7  (STACK_(7))
#define STACK_8  (STACK_(8))
#define STACK_9  (STACK_(9))
#define STACK_10  (STACK_(10))
# etc.
%% { int i;
%%   for (i=0; i<=10; i++)
%%     printf("#define STACK_%d  (STACK_(%d))\n",i,i);
%% }

# Values:

# Highest number of multiple values + 1
#define mv_limit  128
# Values are always passed in the MULTIPLE_VALUE_SPACE mv_space:
# uintC mv_count : number of values, >=0, <mv_limit
# object mv_space [mv_limit-1] : the values.
#   For mv_count>0 the first mv_count elements are occupied.
#   For mv_count=0 the first value = NIL.
#   The values in mv_space are not subject to the Garbage Collection!
#if !defined(mv_count_register)
  # a global Variable
  #ifndef MULTITHREAD
    extern uintC mv_count;
  #else
    #define mv_count  (current_thread()->_mv_count)
  #endif
#else
  # a global register
  register uintC mv_count __asm__(mv_count_register);
#endif
#ifndef MULTITHREAD
  extern object mv_space [mv_limit-1];
#else
  #define mv_space  (current_thread()->_mv_space)
#endif
# Synonyms:
#if !defined(value1_register)
  #ifndef MULTITHREAD
    #define value1  mv_space[0]
  #else
    /* The first value mv_space[0] is moved to the beginning of struct
       clisp_thread_t: */
    #define value1  (current_thread()->_value1)
    #define VALUE1_EXTRA # and thus has to be treated extra every time...
  #endif
#else
  # The first value mv_space[0] is stored permanently in a register:
  register object value1 __asm__(value1_register);
  #define VALUE1_EXTRA # and thus has to be treated extra every time...
#endif
#define value2  mv_space[1]
#define value3  mv_space[2]
#define value4  mv_space[3]
#define value5  mv_space[4]
#define value6  mv_space[5]
#define value7  mv_space[6]
#define value8  mv_space[7]
#define value9  mv_space[8]
# You might need global variables to pass with setjmp/longjmp:
#ifdef NEED_temp_mv_count
  #ifndef MULTITHREAD
    extern uintC temp_mv_count;
  #else
    #define temp_mv_count  (current_thread()->_temp_mv_count)
  #endif
  #define LONGJMP_SAVE_mv_count()  temp_mv_count = mv_count
  #define LONGJMP_RESTORE_mv_count()  mv_count = temp_mv_count
#else
  #define LONGJMP_SAVE_mv_count()
  #define LONGJMP_RESTORE_mv_count()
#endif
#ifdef NEED_temp_value1
  #ifndef MULTITHREAD
    extern object temp_value1;
  #else
    #define temp_value1  (current_thread()->_temp_value1)
  #endif
  #define LONGJMP_SAVE_value1()  temp_value1 = value1
  #define LONGJMP_RESTORE_value1()  value1 = temp_value1
#else
  #define LONGJMP_SAVE_value1()
  #define LONGJMP_RESTORE_value1()
#endif
# is used by EVAL, CONTROL,
#                    Macros LIST_TO_MV, MV_TO_LIST, STACK_TO_MV, MV_TO_STACK
%% #if notused
%% export_def(mv_limit);
%% #endif
%% #if !defined(mv_count_register)
%%   puts("extern uintC mv_count;");
%% #else
%%   puts("#ifndef IN_MODULE_CC");
%%   printf("register uintC mv_count __asm__(\"%s\");\n",mv_count_register);
%%   puts("#endif");
%% #endif
%% printf("extern object mv_space [%d];\n",mv_limit-1);
%% #if !defined(value1_register)
%%   emit_define("value1","mv_space[0]");
%% #else
%%   puts("#ifndef IN_MODULE_CC");
%%   printf("register object value1 __asm__(\"%s\");\n",value1_register);
%%   puts("#endif");
%% #endif
%% { int i = 2;
%%   for (; i <=9 ; i++)
%%     printf("#define value%d  mv_space[%d]\n",i,i-1);
%% }

#ifdef DEBUG_GCSAFETY
  # Add support for the 'mv_space' expression to the GCTRIGGER1/2/... macros.
  inline void inc_allocstamp (object (&mvsp)[mv_limit-1]) {
    inc_allocstamp(value1);
    var uintC count = mv_count;
    if (count > 1) {
      var object* mvp = &mv_space[1];
      dotimespC(count,count-1, { inc_allocstamp(*mvp++); });
    }
  }
#endif

# Returns the bottom objects from the STACK as multiple values.
# STACK_to_mv(count)
# count: number of objects, < mv_limit.
#if !defined(VALUE1_EXTRA)
  #define STACK_to_mv(countx)                                   \
    do { var uintC count = (countx);                            \
      mv_count = count;                                         \
      if (count == 0) value1 = NIL;                             \
      else { # pointer behind space for last value              \
       object* mvp = &mv_space[count];                          \
       dotimespC(count,count, { *--mvp = popSTACK(); } );       \
    }  } while(0)
#else
  #define STACK_to_mv(countx)                                   \
    do { var uintC count = (countx);                            \
      mv_count = count;                                         \
      if (count == 0) value1 = NIL;                             \
      else {                                                    \
        count--;                                                \
        if (count > 0) { # pointer behind space for last value  \
          object* mvp = &mv_space[1+count];                     \
          dotimespC(count,count, { *--mvp = popSTACK(); } );    \
        }                                                       \
        value1 = popSTACK();                                    \
    }  } while(0)
#endif
# is used by EVAL, CONTROL
%% export_def(STACK_to_mv(countx));

# Puts all values onto the STACK.
# mv_to_STACK()
# > mv_count/mv_space : values
# < values on the Stack (first value on top)
# STACK-Overflow is checked.
# modifies STACK
#if !defined(VALUE1_EXTRA)
  #define mv_to_STACK()                                         \
    do { var uintC count = mv_count;                            \
         if (count!=0) { # no values-> nothing onto the STACK   \
           var object* mvp = &mv_space[0];                      \
           get_space_on_STACK(count);                           \
           dotimespC(count,count, { pushSTACK(*mvp++); } );     \
    }  } while(0)
#else
  #define mv_to_STACK()                                         \
    do { var uintC count = mv_count;                            \
         if (count!=0) { # no values -> nothing onto the STACK  \
           get_space_on_STACK(count);                           \
           pushSTACK(value1);                                   \
           count--;                                             \
           if (count > 0) {                                     \
             var object* mvp = &mv_space[1];                    \
             dotimespC(count,count, { pushSTACK(*mvp++); } );   \
           }                                                    \
    }  } while(0)
#endif
# is used by EVAL, CONTROL

# Returns the elements of a list as multiple values.
# list_to_mv(list,fehler_statement)
# fehler_statement: if there's an error (too many values).
#define NEXT_MV  *mvp++ = Car(l); l = Cdr(l); count++
#if !defined(VALUE1_EXTRA)
  #define list_to_mv(lst,fehler_statement)                              \
    do { var object l = (lst);                                          \
     var uintC count = 0;                                               \
     if (atomp(l)) value1 = NIL;                                        \
     else {                                                             \
       var object* mvp = &mv_space[0];                                  \
       NEXT_MV; if (atomp(l)) goto mv_fertig;                           \
       NEXT_MV; if (atomp(l)) goto mv_fertig;                           \
       NEXT_MV; if (atomp(l)) goto mv_fertig;                           \
       do { if (count==mv_limit-1) { fehler_statement; } NEXT_MV;       \
       } while (consp(l));                                              \
     }                                                                  \
     mv_fertig:                                                         \
     if (!nullp(l)) fehler_proper_list_dotted(S(values_list),l);        \
     mv_count = count;                                                  \
    } while(0)
#else
  #define list_to_mv(lst,fehler_statement)                              \
    do { var object l = (lst);                                          \
     var uintC count = 0;                                               \
     if (atomp(l)) value1 = NIL;                                        \
     else {                                                             \
       value1 = Car(l); l = Cdr(l); count++; if (atomp(l)) goto mv_fertig; \
       {var object* mvp = &mv_space[1];                                 \
        NEXT_MV; if (atomp(l)) goto mv_fertig;                          \
        NEXT_MV; if (atomp(l)) goto mv_fertig;                          \
        do { if (count==mv_limit-1) { fehler_statement; } NEXT_MV;      \
        } while (consp(l));                                             \
     }}                                                                 \
     mv_fertig:                                                         \
     if (!nullp(l)) fehler_proper_list_dotted(S(values_list),l);        \
     mv_count = count;                                                  \
    } while(0)
#endif
# is used by EVAL, CONTROL

# Gives the list of the multiple values on -(STACK).
# mv_to_list()
# can trigger GC
#define mv_to_list()                                                  \
  do {                                                                \
    mv_to_STACK(); # at first all values onto the stack               \
    GCTRIGGER();                                                      \
    pushSTACK(NIL); # head of the list                                \
    { var uintC count;                                                \
      dotimesC(count,mv_count, { # until all values have been used:   \
        var object l = allocate_cons(); # new cell                    \
        Cdr(l) = popSTACK(); # list so far                            \
        Car(l) = STACK_0; # next value                                \
        STACK_0 = l; # save new cons                                  \
      });                                                             \
  } } while(0)
# is used by EVAL, CONTROL, DEBUG

# Error message if there are too many values
# fehler_mv_zuviel(caller);
# > caller: caller, a Symbol
nonreturning_function(extern, fehler_mv_zuviel, (object caller));
# is used by EVAL, CONTROL, LISPARIT
%% #if notused
%% puts("nonreturning_function(extern, fehler_mv_zuviel, (object caller));");
%% #endif

#if !defined(back_trace_register)
  #ifndef MULTITHREAD
    extern p_backtrace_t back_trace;
  #else
    #define back_trace  (current_thread()->_back_trace)
  #endif
#else
  register p_backtrace_t back_trace __asm__(back_trace_register);
#endif
#define subr_self  back_trace->bt_function
%% #if !defined(back_trace_register)
%%   puts("extern p_backtrace_t back_trace;");
%% #else
%%   puts("#ifndef IN_MODULE_CC");
%%   printf("register p_backtrace_t back_trace __asm__(\"%s\");\n",back_trace_register);
%%   puts("#endif");
%% #endif
%% export_def(subr_self);

# Within the body of a SUBR: Access to the arguments.
# A SUBR with a fixed number of arguments can access them through the STACK:
#   STACK_0 = last argument, STACK_1 = second to last argument etc.
#   Clean STACK: with skipSTACK(number of arguments) .
# A SUBR with arbitrarily many arguments (&REST-Parameter) gets passed:
#     uintC argcount                    the number of the remaining arguments
#     gcv_object_t* rest_args_pointer   Pointer above the remaining arguments
#   Additionally:
#     gcv_object_t* args_end_pointer    Pointer below all arguments, depends on the STACK
#   Additionally possible:
#     gcv_object_t* args_pointer = rest_args_pointer STACKop (fixed number of arguments);
#                                       Pointer above the first argument
#   Typical Loop-Processing:
#     from the front:
#       until (argcount==0) {
#         var object arg = NEXT(rest_args_pointer); ...; argcount--;
#       }
#       until (rest_args_pointer==args_end_pointer) {
#         var object arg = NEXT(rest_args_pointer); ...;
#       }
#     from the back:
#       until (argcount==0) {
#         var object arg = BEFORE(args_end_pointer); ...; argcount--;
#       }
#       until (rest_args_pointer==args_end_pointer) {
#         var object arg = BEFORE(args_end_pointer); ...;
#       }
#   The macros NEXT and BEFORE modify their arguments!
#   Clean STACK: with set_args_end_pointer(args_pointer)
#     or skipSTACK((fixed number of arguments) + (uintL) (number of remainung arguments)) .
#define args_end_pointer  STACK
#define set_args_end_pointer(new_args_end_pointer)  \
  setSTACK(STACK = (new_args_end_pointer))
#ifdef STACK_DOWN
  #define NEXT(argpointer)  (*(--(argpointer)))
  #define BEFORE(argpointer)  (*((argpointer)++))
#endif
#ifdef STACK_UP
  #define NEXT(argpointer)  (*((argpointer)++))
  #define BEFORE(argpointer)  (*(--(argpointer)))
#endif
# Next(pointer) yields the same value as NEXT(pointer),
# but without changing the value of pointer.
# Before(pointer) yields the same value as BEFORE(pointer),
# but without changing the value of pointer.
#define Next(pointer)  (*(STACKpointable(pointer) STACKop -1))
#define Before(pointer)  (*(STACKpointable(pointer) STACKop 0))
%% emit_define("args_end_pointer","STACK");
%% #if notused
%% emit_define("set_args_end_pointer(new_args_end_pointer)","STACK = (new_args_end_pointer)");
%% export_def(NEXT(argpointer));
%% export_def(BEFORE(argpointer));
%% emit_define("Next(pointer)","(*(STACKpointable(pointer) STACKop -1))");
%% emit_define("Before(pointer)","(*(STACKpointable(pointer) STACKop 0))");
%% #endif

# Environments:

typedef struct {
  object var_env;   # Variable-Bindings-Environment
  object fun_env;   # Function-Bindings-Environment
  object block_env; # Block-Environment
  object go_env;    # Tagbody/Go-Environment
  object decl_env;  # Declarations-Environment
} environment_t;
typedef struct {
  gcv_object_t var_env;   # Variable-Bindings-Environment
  gcv_object_t fun_env;   # Function-Bindings-Environment
  gcv_object_t block_env; # Block-Environment
  gcv_object_t go_env;    # Tagbody/Go-Environment
  gcv_object_t decl_env;  # Declarations-Environment
} gcv_environment_t;

# The current Environment:
#ifndef MULTITHREAD
  extern gcv_environment_t aktenv;
#else
  #define aktenv  (current_thread()->_aktenv)
#endif

# Macro: Puts five single Environments on the STACK
# and makes a single Environment out of them.
# make_STACK_env(venv,fenv,benv,genv,denv, env5 = );
# > object venv,fenv,benv,genv,denv: 5 single Environments
# < gcv_environment_t* env5: pointer to the Environment on the Stack
#ifdef STACK_UP
  #define make_STACK_env(venv,fenv,benv,genv,denv,env5_assignment)      \
    do { pushSTACK(venv); pushSTACK(fenv); pushSTACK(benv);             \
         pushSTACK(genv); pushSTACK(denv);                              \
         env5_assignment &STACKblock_(gcv_environment_t,0); } while(0)
#endif
#ifdef STACK_DOWN
  #define make_STACK_env(venv,fenv,benv,genv,denv,env5_assignment)      \
    do { pushSTACK(denv); pushSTACK(genv); pushSTACK(benv);             \
         pushSTACK(fenv); pushSTACK(venv);                              \
         env5_assignment &STACKblock_(gcv_environment_t,0); } while(0)
#endif

# Frameinfobits in Frames:
# in the Frame-Info-Byte (tint):
#if (oint_type_len>=7) && 0 # provisionally??
  # Bit numbers in the Frame-Info-Byte:
  # occupy Bits 6..0 (resp. Bits 7,5..0 if garcol_bit_t=7).
  #ifdef TYPECODES
    #define FB7  garcol_bit_t
    #define FB6  (garcol_bit_t>TB5 ? TB5 : TB6)
    #define FB5  (garcol_bit_t>TB4 ? TB4 : TB5)
    #define FB4  (garcol_bit_t>TB3 ? TB3 : TB4)
    #define FB3  (garcol_bit_t>TB2 ? TB2 : TB3)
    #define FB2  (garcol_bit_t>TB1 ? TB1 : TB2)
    #define FB1  (garcol_bit_t>TB0 ? TB0 : TB1)
  #else
    #define FB7  garcol_bit_o
    #define FB6  30
    #define FB5  29
    #define FB4  28
    #define FB3  27
    #define FB2  26
    #define FB1  25
  #endif
  # depending on it:
  #define frame_bit_t    FB7  # garcol_bit as FRAME-identifier
  #define skip2_bit_t    FB6  # unset if GC has to skip two longwords
  #define unwind_bit_t   FB5  # set if there's something to do while
                              # unwinding the frame
  # skip2-Bit=1 ==> unwind-Bit=1.
  # for further Information within the Frames with skip2-Bit=1:
  #define envbind_bit_t  FB4  # Bit set for ENV-Frames.
                              # Bit is unset for DYNBIND-Frames.
  # for further identification of the ENV-Frames:
  #define envbind_case_mask_t  (bit(FB3)|bit(FB2)|bit(FB1))
  # for further discrimination within the Frames with skip2-Bit=0:
  #define entrypoint_bit_t  FB4  # Bit is set, if FRAME contains
  # a non-local entrypoint, with Offset SP_, SP is on the STACK.
  # Bit is unset for VAR/FUN-Frame and CALLBACK-Frame.
  # for further discrimination in BLOCK/TAGBODY/APPLY/EVAL/CATCH/UNWIND_PROTECT/HANDLER/DRIVER:
  #define blockgo_bit_t    FB3  # Bit set for BLOCK- and TAGBODY-FRAME
  # for further discrimination in BLOCK/TAGBODY:
  #define cframe_bit_t     FB1  # set for compiled BLOCK/TAGBODY-Frames,
                                # unset for interpreted BLOCK/TAGBODY-Frames
  #define nested_bit_t unwind_bit_t # set for IBLOCK and ITAGBODY,
                                    # if Exitpoint resp. Tags were nested
  # for further discrimination in APPLY/EVAL/CATCH/UNWIND_PROTECT/HANDLER/DRIVER:
  #define dynjump_bit_t  FB2    # unset for APPLY and EVAL, set
                                # for CATCH/UNWIND_PROTECT/DRIVER-Frames
  #define trapped_bit_t unwind_bit_t # set for APPLY and EVAL, if
                                # interrupted while unwinding the Frame
  # unwind-Bit set for UNWIND_PROTECT/DRIVER/TRAPPED_APPLY/TRAPPED_EVAL,
  # else unset.
  #define eval_bit_t     FB1    # set for EVAL-Frames,
                                # unset for APPLY-Frames
  #define driver_bit_t   FB1    # set for DRIVER-Frames,
                                # unset for UNWIND_PROTECT-Frames
  #define handler_bit_t  FB1    # set for HANDLER-Frames,
                                # unset for CATCH-Frames
  # for further discrimination in VAR/FUN/CALLBACK:
  #define callback_bit_t   FB3  # Bit is unset for CALLBACK-Frames.
                                # Bit is set for VAR/FUN-Frames.
  # for further discrimination in VAR/FUN:
  #define fun_bit_t      FB2    # set for FUN-Frame, unset for VAR-Frame
  # on objects on the STACK (oint):
  #define      frame_bit_o      (frame_bit_t+oint_type_shift)
  #define      skip2_bit_o      (skip2_bit_t+oint_type_shift)
  #define     unwind_bit_o     (unwind_bit_t+oint_type_shift)
  #define    envbind_bit_o    (envbind_bit_t+oint_type_shift)
  #define   callback_bit_o   (callback_bit_t+oint_type_shift)
  #define entrypoint_bit_o (entrypoint_bit_t+oint_type_shift)
  #define    blockgo_bit_o    (blockgo_bit_t+oint_type_shift)
  #define     cframe_bit_o     (cframe_bit_t+oint_type_shift)
  #define     nested_bit_o     (nested_bit_t+oint_type_shift)
  #define    dynjump_bit_o    (dynjump_bit_t+oint_type_shift)
  #define    trapped_bit_o    (trapped_bit_t+oint_type_shift)
  #define       eval_bit_o       (eval_bit_t+oint_type_shift)
  #define     driver_bit_o     (driver_bit_t+oint_type_shift)
  #define    handler_bit_o    (handler_bit_t+oint_type_shift)
  #define        fun_bit_o        (fun_bit_t+oint_type_shift)
  # single Frame-Info-Bytes:
  #define DYNBIND_frame_info          /* %1110... */ (bit(FB7)|bit(FB6)|bit(FB5))
  #define ENV1V_frame_info            /* %1111000 */ (bit(FB7)|bit(FB6)|bit(FB5)|bit(FB4))
  #define ENV1F_frame_info            /* %1111001 */ (bit(FB7)|bit(FB6)|bit(FB5)|bit(FB4)|bit(FB1))
  #define ENV1B_frame_info            /* %1111010 */ (bit(FB7)|bit(FB6)|bit(FB5)|bit(FB4)|bit(FB2))
  #define ENV1G_frame_info            /* %1111011 */ (bit(FB7)|bit(FB6)|bit(FB5)|bit(FB4)|bit(FB2)|bit(FB1))
  #define ENV1D_frame_info            /* %1111100 */ (bit(FB7)|bit(FB6)|bit(FB5)|bit(FB4)|bit(FB3))
  #define ENV2VD_frame_info           /* %1111101 */ (bit(FB7)|bit(FB6)|bit(FB5)|bit(FB4)|bit(FB3)|bit(FB1))
  #define ENV5_frame_info             /* %1111110 */ (bit(FB7)|bit(FB6)|bit(FB5)|bit(FB4)|bit(FB3)|bit(FB2))
  #ifdef HAVE_SAVED_REGISTERS
    #define CALLBACK_frame_info         /* %10100.. */ (bit(FB7)|bit(FB5))
  #endif
  #define VAR_frame_info              /* %101010. */ (bit(FB7)|bit(FB5)|bit(FB3))
  #define FUN_frame_info              /* %101011. */ (bit(FB7)|bit(FB5)|bit(FB3)|bit(FB2))
  #define IBLOCK_frame_info           /* %1001100 */ (bit(FB7)|bit(FB4)|bit(FB3))
  #define NESTED_IBLOCK_frame_info    /* %1011100 */ (bit(FB7)|bit(FB5)|bit(FB4)|bit(FB3))
  #define ITAGBODY_frame_info         /* %1001110 */ (bit(FB7)|bit(FB4)|bit(FB3)|bit(FB2))
  #define NESTED_ITAGBODY_frame_info  /* %1011110 */ (bit(FB7)|bit(FB5)|bit(FB4)|bit(FB3)|bit(FB2))
  #define CBLOCK_CTAGBODY_frame_info  /* %1011101 */ (bit(FB7)|bit(FB5)|bit(FB4)|bit(FB3)|bit(FB1))
  #define APPLY_frame_info            /* %1001000 */ (bit(FB7)|bit(FB4))
  #define TRAPPED_APPLY_frame_info    /* %1011000 */ (bit(FB7)|bit(FB5)|bit(FB4))
  #define EVAL_frame_info             /* %1001001 */ (bit(FB7)|bit(FB4)|bit(FB1))
  #define TRAPPED_EVAL_frame_info     /* %1011001 */ (bit(FB7)|bit(FB5)|bit(FB4)|bit(FB1))
  #define CATCH_frame_info            /* %1001010 */ (bit(FB7)|bit(FB4)|bit(FB2))
  #define HANDLER_frame_info          /* %1001011 */ (bit(FB7)|bit(FB4)|bit(FB2)|bit(FB1))
  #define UNWIND_PROTECT_frame_info   /* %1011010 */ (bit(FB7)|bit(FB5)|bit(FB4)|bit(FB2))
  #define DRIVER_frame_info           /* %1011011 */ (bit(FB7)|bit(FB5)|bit(FB4)|bit(FB2)|bit(FB1))
#endif
#if (oint_type_len==6) || 1 # provisionally??
  # bit numbers in Frame-Info-Byte:
  # occupy Bits 5..0 (resp. Bits 7,4..0 if garcol_bit_t=7).
  #ifdef TYPECODES
    #define FB6  garcol_bit_t
    #define FB5  (garcol_bit_t>TB4 ? TB4 : TB5)
    #define FB4  (garcol_bit_t>TB3 ? TB3 : TB4)
    #define FB3  (garcol_bit_t>TB2 ? TB2 : TB3)
    #define FB2  (garcol_bit_t>TB1 ? TB1 : TB2)
    #define FB1  (garcol_bit_t>TB0 ? TB0 : TB1)
  #else # HEAPCODES
    #define FB6  garcol_bit_o
    #ifdef STANDARD_HEAPCODES
      #define FB5  (garcol_bit_o-1)
      #define FB4  (garcol_bit_o-2)
      #define FB3  (garcol_bit_o-3)
      #define FB2  (garcol_bit_o-4)
      #define FB1  (garcol_bit_o-5)
    #endif
    #ifdef LINUX_NOEXEC_HEAPCODES
      #define FB5  5
      #define FB4  4
      #define FB3  3
      #define FB2  2
      #define FB1  1
    #endif
  #endif
  # depending on it:
  #define frame_bit_t    FB6  # garcol_bit as FRAME-indicator
  #define skip2_bit_t    FB5  # unset if the GC has to skip two long words
  # define unwind_limit_t  ...  # above:
  # if there's something to be done while to unwind a Frame
  # skip2-Bit=1 ==> >= unwind-limit.
  # for further information within the Frames with skip2-Bit=1:
  #define envbind_bit_t  FB4  # Bit is set for ENV-Frames.
                              # Bit unset for DYNBIND-Frames.
  # for further identification within the ENV-Frames:
  #define envbind_case_mask_t  (bit(FB3)|bit(FB2)|bit(FB1))
  # for further discrimination with the Frames with skip2-Bit=0:
  # define entrypoint_limit_t  ...  # below:
  # if FRAME contains a non-local entry point
  # with Offset SP_ SP is on the STACK.
  # above: for VAR/FUN-Frame and CALLBACK-Frame.
  # for further discrimination in BLOCK/TAGBODY/APPLY/EVAL/CATCH/UNWIND_PROTECT/HANDLER/DRIVER:
  #define blockgo_bit_t    FB3  # Bit set for BLOCK- and TAGBODY-FRAME
  # for further discrimination in BLOCK/TAGBODY:
  #define cframe_bit_t   FB4  # set for compiled, unset for
                              # interpreted BLOCK/TAGBODY-Frames
  #define nested_bit_t   FB2  # set for IBLOCK and ITAGBODY,
                              # if exit point or Tags have been nested
  # for further discrimination in APPLY/EVAL/CATCH/UNWIND_PROTECT/HANDLER/DRIVER:
  #define dynjump_bit_t  FB2  # unset for APPLY and EVAL, set
                              # for CATCH/UNWIND_PROTECT/HANDLER/DRIVER-Frames
  #define trapped_bit_t  FB4  # set for APPLY and EVAL, if interruped while
                              # unwinding the Frames
  # >= unwind_limit_t for UNWIND_PROTECT/DRIVER/TRAPPED_APPLY/TRAPPED_EVAL,
  # < unwind_limit_t else.
  #define eval_bit_t     FB1  # set for EVAL-Frames,
                              # unset for APPLY-Frames
  #define driver_bit_t   FB1  # set for DRIVER-Frames,
                              # unset for UNWIND_PROTECT-Frames
  #define handler_bit_t  FB1  # set for HANDLER-Frames,
                              # unset for CATCH-Frames
  # for further discrimination in VAR/FUN/CALLBACK:
  #define callback_bit_t   FB2  # Bit is unset for CALLBACK-Frames.
                                # Bit is set for VAR/FUN-Frames.
  # for further discrimination in VAR/FUN:
  #define fun_bit_t      FB1  # set for FUN-Frame, unset for VAR-Frame
  # in Objects on the STACK (oint):
  #define    frame_bit_o    (frame_bit_t+oint_type_shift)
  #define    skip2_bit_o    (skip2_bit_t+oint_type_shift)
  #define  envbind_bit_o  (envbind_bit_t+oint_type_shift)
  #define callback_bit_o (callback_bit_t+oint_type_shift)
  #define  blockgo_bit_o  (blockgo_bit_t+oint_type_shift)
  #define   cframe_bit_o   (cframe_bit_t+oint_type_shift)
  #define   nested_bit_o   (nested_bit_t+oint_type_shift)
  #define  dynjump_bit_o  (dynjump_bit_t+oint_type_shift)
  #define  trapped_bit_o  (trapped_bit_t+oint_type_shift)
  #define     eval_bit_o     (eval_bit_t+oint_type_shift)
  #define   driver_bit_o   (driver_bit_t+oint_type_shift)
  #define  handler_bit_o  (handler_bit_t+oint_type_shift)
  #define      fun_bit_o      (fun_bit_t+oint_type_shift)
  # single Frame-Info-Bytes:
  #define APPLY_frame_info            /* %100000 */ (bit(FB6))
  #define EVAL_frame_info             /* %100001 */ (bit(FB6)|bit(FB1))
  #define CATCH_frame_info            /* %100010 */ (bit(FB6)|bit(FB2))
  #define HANDLER_frame_info          /* %100011 */ (bit(FB6)|bit(FB2)|bit(FB1))
  #define IBLOCK_frame_info           /* %100100 */ (bit(FB6)|bit(FB3))
  #define ITAGBODY_frame_info         /* %100101 */ (bit(FB6)|bit(FB3)|bit(FB1))
  #define unwind_limit_t                            (bit(FB6)|bit(FB3)|bit(FB2))
  #define NESTED_IBLOCK_frame_info    /* %100110 */ (bit(FB6)|bit(FB3)|bit(FB2))
  #define NESTED_ITAGBODY_frame_info  /* %100111 */ (bit(FB6)|bit(FB3)|bit(FB2)|bit(FB1))
  #define TRAPPED_APPLY_frame_info    /* %101000 */ (bit(FB6)|bit(FB4))
  #define TRAPPED_EVAL_frame_info     /* %101001 */ (bit(FB6)|bit(FB4)|bit(FB1))
  #define UNWIND_PROTECT_frame_info   /* %101010 */ (bit(FB6)|bit(FB4)|bit(FB2))
  #define DRIVER_frame_info           /* %101011 */ (bit(FB6)|bit(FB4)|bit(FB2)|bit(FB1))
  #define CBLOCK_CTAGBODY_frame_info  /* %101100 */ (bit(FB6)|bit(FB4)|bit(FB3))
  #define entrypoint_limit_t                        (bit(FB6)|bit(FB4)|bit(FB3)|bit(FB1))
  #ifdef HAVE_SAVED_REGISTERS
  #define CALLBACK_frame_info         /* %101101 */ (bit(FB6)|bit(FB4)|bit(FB3)|bit(FB1))
  #endif
  #define VAR_frame_info              /* %101110 */ (bit(FB6)|bit(FB4)|bit(FB3)|bit(FB2))
  #define FUN_frame_info              /* %101111 */ (bit(FB6)|bit(FB4)|bit(FB3)|bit(FB2)|bit(FB1))
  #define DYNBIND_frame_info          /* %110... */ (bit(FB6)|bit(FB5))
  #define ENV1V_frame_info            /* %111000 */ (bit(FB6)|bit(FB5)|bit(FB4))
  #define ENV1F_frame_info            /* %111001 */ (bit(FB6)|bit(FB5)|bit(FB4)|bit(FB1))
  #define ENV1B_frame_info            /* %111010 */ (bit(FB6)|bit(FB5)|bit(FB4)|bit(FB2))
  #define ENV1G_frame_info            /* %111011 */ (bit(FB6)|bit(FB5)|bit(FB4)|bit(FB2)|bit(FB1))
  #define ENV1D_frame_info            /* %111100 */ (bit(FB6)|bit(FB5)|bit(FB4)|bit(FB3))
  #define ENV2VD_frame_info           /* %111101 */ (bit(FB6)|bit(FB5)|bit(FB4)|bit(FB3)|bit(FB1))
  #define ENV5_frame_info             /* %111110 */ (bit(FB6)|bit(FB5)|bit(FB4)|bit(FB3)|bit(FB2))
#endif
#define CBLOCK_frame_info  CBLOCK_CTAGBODY_frame_info
#define CTAGBODY_frame_info  CBLOCK_CTAGBODY_frame_info
%% #ifdef HAVE_SAVED_REGISTERS
%%   export_def(CALLBACK_frame_info);
%% #endif

# Bits for Symbols in VAR-Frames:
# bit(active_bit),bit(dynam_bit),bit(svar_bit) must fit into one uintB:
#if !((active_bit<intBsize) && (dynam_bit<intBsize) && (svar_bit<intBsize))
  #error "Symbol bits don't fit in a single byte -- Symbol-Bits passen nicht in ein Byte!"
#endif
#ifdef NO_symbolflags
  # Bits are separatly stored on the Stack as Fixnums.
  #undef oint_symbolflags_shift
  #define oint_symbolflags_shift  oint_data_shift
#else
  #if (oint_symbolflags_shift==oint_addr_shift)
    # bit(active_bit),bit(dynam_bit),bit(svar_bit) must be true divisors
    # of varobject_alignment:
    #if (varobject_alignment % bit(active_bit+1)) || (varobject_alignment % bit(dynam_bit+1)) || (varobject_alignment % bit(svar_bit+1))
      #error "No more room for three bits in a symbol -- Kein Platz fuer drei Bits in der Adresse eines Symbols!"
    #endif
  #endif
#endif
#define active_bit_o  (active_bit+oint_symbolflags_shift)  # set: binding is active
#define dynam_bit_o   (dynam_bit+oint_symbolflags_shift)   # set: binding is dynamic
#define svar_bit_o    (svar_bit+oint_symbolflags_shift)    # set: next parameter is supplied-p-parameter for this

# Offsets for data in Frames, to be addressed via STACK_(Offset)
#define frame_form      2  # EVAL
#define frame_closure   2  # APPLY, HANDLER
#define frame_anz       1  # VAR, FUN
#define frame_SP        1  # IBLOCK, CBLOCK, ITAGBODY, CTAGBODY,
                           # EVAL, CATCH, UNWIND-PROTECT, HANDLER, DRIVER
#define frame_next_env  2  # VAR, FUN, IBLOCK, ITAGBODY
#define frame_ctag      2  # CBLOCK, CTAGBODY
#define frame_tag       2  # CATCH
#define frame_handlers  3  # HANDLER
#define frame_name      3  # IBLOCK
#define frame_args      3  # APPLY
#define frame_bindings  3  # VAR, FUN, ITAGBODY
# Structure of the different bindings in VAR-Frames:
#ifdef NO_symbolflags
  #define varframe_binding_size  3
  #define varframe_binding_mark   0
  #define varframe_binding_sym    1
  #define varframe_binding_value  2
  #define pushSTACK_symbolwithflags(symbol,flags)  \
    pushSTACK(symbol); pushSTACK(as_object(as_oint(Fixnum_0) | (oint)(flags)))
#else
  #define varframe_binding_size  2
  #define varframe_binding_mark   0
  #define varframe_binding_sym    0
  #define varframe_binding_value  1
  #define pushSTACK_symbolwithflags(symbol,flags)  \
    pushSTACK(as_object(as_oint(symbol) | (oint)(flags)))
#endif

/* Special value to mark BLOCK- and TAGBODY-references that are not 'live'
   anymore (replaces the Frame-Pointer in the CDR of the corresponding Cons) */
#define disabled  make_system(0xDDDDDDUL)

# Value to mark specially declared references
#define specdecl  make_system(0xECDECDUL)

# Handling Frames:
# A local variable FRAME contains the value of STACK after
# creating a Frame. Then you can access with FRAME_(n) just like
# with likeSTACK_(n):
#ifdef STACK_DOWN
  #define FRAME_(n)  (FRAME[(sintP)(n)])
#endif
#ifdef STACK_UP
  #define FRAME_(n)  (FRAME[-1-(sintP)(n)])
#endif
# make_framepointer(FRAME) is the Frame-Pointer as LISP-object.
# framecode(FRAME_(0)) is the Frame-Info-Byte (of Type fcint),
# topofframe(FRAME_(0)) is a Pointer above the Frame.
# FRAME = uTheFramepointer(obj) is a Frame-Pointer as pointer into the Stack.
#         [uTheFramepointer is the exact opposite of make_framepointer.]
# FRAME = TheFramepointer(obj) as well, but possibly still with type info!
#         [An attenuation of uTheFramepointer, that is enough for access.]
#ifdef TYPECODES
  #if !defined(SINGLEMAP_MEMORY_STACK)
    #define make_framepointer(stack_ptr)  type_pointer_object(system_type,stack_ptr)
    #define topofframe(bottomword)  (gcv_object_t*)upointer(bottomword)
    #define uTheFramepointer(obj)  (gcv_object_t*)upointer(obj)
  #else
    #define make_framepointer(stack_ptr)  (as_object((oint)(stack_ptr)))
    #define topofframe(bottomword)  (gcv_object_t*)as_oint(type_pointer_object(system_type,upointer(bottomword)))
    #define uTheFramepointer(obj)  TheFramepointer(obj) # = (gcv_object_t*)(obj)
  #endif
  #define framecode(bottomword)  mtypecode(bottomword)
  typedef tint fcint;
#else
  # Here the bottomword consists of the frame size, not the top of frame itself.
  # This leaves room for the frame info byte.
  #define make_framepointer(stack_ptr)  make_machine(stack_ptr)
  #ifdef STANDARD_HEAPCODES
    #define makebottomword(type,size)  as_object((oint)(type)+(oint)(size))
    #define framecode(bottomword)  (as_oint(bottomword) & minus_wbit(FB1))
    #define framesize(bottomword)  (as_oint(bottomword)&(wbit(FB1)-1))
  #endif
  #ifdef LINUX_NOEXEC_HEAPCODES
    #define makebottomword(type,size)  as_object((oint)(type)+((oint)(size)<<6))
    #define framecode(bottomword)  (as_oint(bottomword) & 0x3F)
    #define framesize(bottomword)  (as_oint(bottomword) >> 6)
  #endif
  #ifdef STACK_UP
    #define topofframe(bottomword)  \
      (gcv_object_t*)((uintP)(&(bottomword))-(uintP)framesize(bottomword)+sizeof(gcv_object_t))
  #endif
  #ifdef STACK_DOWN
    #define topofframe(bottomword)  \
      (gcv_object_t*)((uintP)(&(bottomword))+(uintP)framesize(bottomword))
  #endif
  #define uTheFramepointer(obj)  TheFramepointer(obj) # = (gcv_object_t*)(obj)
  typedef oint fcint;
#endif
# is used by EVAL, CONTROL, DEBUG
%% #ifdef HEAPCODES
%%  export_def(makebottomword(type,size));
%% #endif
%% export_def(framecode(bottomword));

# To determine the size of a frame:
# STACK_item_count(new_STACK_ptr,old_STACK_ptr)
# calculates the number of STACK-elements between an older stack pointer
# old_STACK_ptr and a new one new_STACK_ptr.
# (That's count with old_STACK_ptr = new_STACK_ptr STACKop count .)
#ifdef STACK_DOWN
  #define STACK_item_count(new_STACK_ptr,old_STACK_ptr)  \
    (uintL)((old_STACK_ptr) - (new_STACK_ptr))
#endif
#ifdef STACK_UP
  #define STACK_item_count(new_STACK_ptr,old_STACK_ptr)  \
    (uintL)((new_STACK_ptr) - (old_STACK_ptr))
#endif

# Finishes a frame.
# finish_frame(frametype);
# > gcv_object_t* top_of_frame: pointer to the top of the frame
# decreases STACK by 1
#ifdef TYPECODES
  #if !defined(SINGLEMAP_MEMORY_STACK)
    #define framebottomword(type,top_of_frame,bot_of_frame)  \
      type_pointer_object(type,top_of_frame)
  #else # top_of_frame has already Typinfo system_type
    #define framebottomword(type,top_of_frame,bot_of_frame)  \
      as_object(type_zero_oint(type)-type_zero_oint(system_type)+(oint)(top_of_frame))
  #endif
  #define finish_frame(frametype)  \
    pushSTACK(framebottomword(frametype##_frame_info,top_of_frame,bot_of_frame_ignored))
#else
  #ifdef STACK_UP
    #define framebottomword(type,top_of_frame,bot_of_frame)  \
      makebottomword(type,(uintP)(bot_of_frame)-(uintP)(top_of_frame))
  #endif
  #ifdef STACK_DOWN
    #define framebottomword(type,top_of_frame,bot_of_frame)  \
      makebottomword(type,(uintP)(top_of_frame)-(uintP)(bot_of_frame))
  #endif
  #define finish_frame(frametype)  \
    (STACK_(-1) = framebottomword(frametype##_frame_info,top_of_frame,STACK STACKop -1), skipSTACK(-1))
#endif
# is used by EVAL, CONTROL
%% export_def(framebottomword(type,top_of_frame,bot_of_frame));
%% export_def(finish_frame(frametype));

# Makes a Frame for all 5 Environments
# make_ENV5_frame();
# decreases STACK by 5
#define make_ENV5_frame()                       \
  do { var gcv_object_t* top_of_frame = STACK;  \
       pushSTACK(aktenv.decl_env);              \
       pushSTACK(aktenv.go_env);                \
       pushSTACK(aktenv.block_env);             \
       pushSTACK(aktenv.fun_env);               \
       pushSTACK(aktenv.var_env);               \
       finish_frame(ENV5);                      \
  } while(0)
# is used by EVAL, CONTROL, DEBUG

# Finishes a Frame with entry point and places jump-point here.
# finish_entry_frame(frametype,returner,retval_assignment,reentry_statement);
# > gcv_object_t* top_of_frame: pointer to the top of the frame
# > sp_jmp_buf* returner: longjmp-Buffer for re-entry
# > retval_assignment: allocated of the setjmp()-value to a variable
# > reentry_statement: what is to be done immediately after re-entry.
# decreases STACK by 1
#define finish_entry_frame(frametype,returner,retval_assignment,reentry_statement)  \
  do { pushSTACK(fake_gcv_object((aint)(returner))); # SP onto the Stack      \
    pushSTACK(nullobj); # dummy onto the Stack, until re-entry is permitted   \
    begin_setjmp_call();                                                      \
    if ((retval_assignment setjmpspl(returner))!=0) # set point for returner  \
      { end_longjmp_call(); LONGJMP_RESTORE_mv_count(); LONGJMP_RESTORE_value1(); reentry_statement } # after re-entry \
    else                                                                      \
      { end_setjmp_call(); STACK_0 = framebottomword(frametype##_frame_info,top_of_frame,STACK); } \
  } while(0)
# is used by EVAL, CONTROL, DEBUG

# Jumps to a Frame with entry point that starts at STACK.
# (Important: The STACK has to have the same values it had when the
# frame was created, since the STACK might not be saved at setjmp/longjmp)
# Never returns and cleans the SP!!
# The multiple values are passed.
# enter_frame_at_STACK();
#define enter_frame_at_STACK()                                do {      \
  /* the returner of finish_entry_frame: */                             \
  var sp_jmp_buf* returner = (sp_jmp_buf*)(aint)as_oint(STACK_(frame_SP)); \
  unwind_back_trace(back_trace,STACK);                                  \
  LONGJMP_SAVE_value1(); LONGJMP_SAVE_mv_count();                       \
  begin_longjmp_call();                                                 \
  longjmpspl(*returner,(aint)returner);/* jump there, pass own addess (/=0) */\
  NOTREACHED;                                                           \
 } while(0)
# is used by EVAL

# Makes a HANDLER-Frame with C-Handler.
# make_HANDLER_frame(types_labels_vector_list,handler,sp_arg);
# make_HANDLER_entry_frame(types_labels_vector_list,handler,returner,reentry_statement);
# > object types_labels_vector_list: a list containing a simple-vector: (#(type1 label1 ... typem labelm))
# > handler: void (*) (void* sp, gcv_object_t* frame, object label, object condition)
# > sp_arg: any void*
# > sp_jmp_buf* returner: longjmp-Buffer for re-entry
# > reentry_statement: what is to be done right after the re-entry.
#define make_HANDLER_frame(types_labels_vector_list,handler,sp_arg)  \
  do { var gcv_object_t* top_of_frame = STACK;     \
       pushSTACK(types_labels_vector_list);        \
       pushSTACK(make_machine_code(handler));      \
       pushSTACK(fake_gcv_object((aint)(sp_arg))); \
       finish_frame(HANDLER);                      \
  } while(0)
#define make_HANDLER_entry_frame(types_labels_vector_list,handler,returner,reentry_statement)  \
  do { var gcv_object_t* top_of_frame = STACK;                  \
       pushSTACK(types_labels_vector_list);                     \
       pushSTACK(make_machine_code(handler));                   \
       finish_entry_frame(HANDLER,returner,,reentry_statement); \
  } while(0)
#define unwind_HANDLER_frame()  skipSTACK(4)

# UP: Applies a function to its arguments.
# apply(function,args_on_stack,other_args);
# > function: function
# > Arguments: args_on_stack arguments on the STACK,
#              remaining argument-list in other_args
# < STACK: cleaned (ie. STACK is increased by args_on_stack)
# < mv_count/mv_space: values
# modifies STACK, can trigger GC
extern maygc Values apply (object fun, uintC args_on_stack, object other_args);
# is used by EVAL, CONTROL, IO, PATHNAME, ERROR
%% #if notused
%% puts("extern Values apply (object fun, uintC args_on_stack, object other_args);");
%% #endif

# UP: Applies a function to its arguments.
# funcall(function,argcount);
# > function: function
# > Arguments: argcount arguments on the STACK
# < STACK: cleaned (ie. STACK is increased by argcount)
# < mv_count/mv_space: values
# modifies STACK, can trigger GC
extern maygc Values funcall (object fun, uintC argcount);
# is used by all Modules
%% puts("extern Values funcall (object fun, uintC argcount);");

# UP: Evaluates a Form in the current Environment.
# eval(form);
# > form: Form
# < mv_count/mv_space: values
# can trigger GC
extern maygc Values eval (object form);
# is used by CONTROL, DEBUG
%% #if notused
%% puts("extern Values eval (object form);");
%% #endif

# UP: Evaluates a Form in a given Environment.
# eval_5env(form,var,fun,block,go,decl);
# > var_env: Value for VAR_ENV
# > fun_env: Value for FUN_ENV
# > block_env: Value for BLOCK_ENV
# > go_env: Value for GO_ENV
# > decl_env: Value for DECL_ENV
# > form: Form
# < mv_count/mv_space: Values
# can trigger GC
extern maygc Values eval_5env (object form, object var_env, object fun_env, object block_env, object go_env, object decl_env);
# is used by

# UP: Evaluates a Form in an empty Environment.
# eval_noenv(form);
# > form: Form
# < mv_count/mv_space: Values
# can trigger GC
extern maygc Values eval_noenv (object form);
# is used by CONTROL, IO, DEBUG, SPVW

# UP: Evaluates a Form in the current Environment. Doesn't care about
# *EVALHOOK* and *APPLYHOOK*.
# eval_no_hooks(form);
# > form: Form
# < mv_count/mv_space: Values
# can trigger GC
extern maygc Values eval_no_hooks (object form);
# is used by CONTROL

# UP: binds *EVALHOOK* and *APPLYHOOK* dynamically to the given values.
# bindhooks(evalhook_value,applyhook_value);
# > evalhook_value: Value for *EVALHOOK*
# > applyhook_value: Value for *APPLYHOOK*
# modifies STACK
extern void bindhooks (object evalhook_value, object applyhook_value);
# is used by CONTROL

# UP: Unwinds a Frame, to which STACK points.
# unwind();
# The values mv_count/mv_space aren't changed.
# If it is no Unwind-Protect-Frame: returns normally.
# If it is an Unwind-Protect-Frame:
#   saves the values, climbs up the STACK and SP
#   and then jumps to unwind_protect_to_save.fun.
# modifies STACK
# can trigger GC
#ifdef __cplusplus
  /* g++-3.4 doesn't like nonreturning in a typedef */
  typedef /* nonreturning */ /*maygc*/ void (*restartf_t)(gcv_object_t* upto_frame);
#else
  nonreturning_function(typedef /*maygc*/, (*restartf_t), (gcv_object_t* upto_frame));
#endif
typedef struct {
  restartf_t fun;
  gcv_object_t* upto_frame;
} unwind_protect_caller_t;
#ifndef MULTITHREAD
  extern unwind_protect_caller_t unwind_protect_to_save;
#else
  #define unwind_protect_to_save  (current_thread()->_unwind_protect_to_save)
#endif
extern /*maygc*/ void unwind (void);
# is used by CONTROL, DEBUG, SPVW

/* UP: "unwinds" the STACK to the next DRIVER_FRAME and
 jumps to the corresponding Top-Level-loop
 if count=0, unwind to TOP; otherwise reset that many times */
nonreturning_function(extern, reset, (uintL count));
/* is used by SPVW, CONTROL */

/* UP: binds the symbols of the list symlist dynamically
 to the values of the list vallist.
 progv(symlist,vallist);
 > symlist, vallist: two lists
 Exactly one variable-bindings-frame is created.
 modifies STACK
 can trigger GC */
extern maygc void progv (object symlist, object vallist);
/* used by CONTROL, EVAL */

# UP: Unwinds the dynamic nesting on the STACK until the frame
# (exclusively), to which upto points, and jumps to it.
# unwind_upto(upto);
# > upto: Pointer to a Frame (into the Stack, without type-info).
# Saves the values mv_count/mv_space.
# modifies STACK,SP
# can trigger GC
# Jumps to the found Frame.
nonreturning_function(extern /*maygc*/, unwind_upto, (gcv_object_t* upto_frame));
# is used by CONTROL, DEBUG

# UP: throws to the Tag tag and passes the values mv_count/mv_space.
# Only returns, if there is no CATCH-Frame of this tag.
# throw_to(tag);
extern void throw_to (object tag);
# is used by CONTROL

# UP: Invokes all handlers for the condition cond. Only returns, if none
# of the handlers feels responsible (ie. if every handler returns).
# invoke_handlers(cond);
# can trigger GC
extern maygc void invoke_handlers (object cond);
typedef struct {
  object condition;
  gcv_object_t* stack;
  SPint* sp;
  object spdepth;
} handler_args_t;
#ifndef MULTITHREAD
  extern handler_args_t handler_args;
#else
  #define handler_args  (current_thread()->_handler_args)
#endif
typedef struct stack_range_t {
  struct stack_range_t * next;
  gcv_object_t* low_limit;
  gcv_object_t* high_limit;
} stack_range_t;
#ifndef MULTITHREAD
  extern stack_range_t* inactive_handlers;
#else
  #define inactive_handlers  (current_thread()->_inactive_handlers)
#endif
# is used by ERROR

# UP: Determines, whether an Object is a function name, ie. a Symbol or
# a list of the form (SETF symbol).
# funnamep(obj)
# > obj: Objekt
# < result: true if function name
extern bool funnamep (object obj);
# is used by CONTROL

# Gives the block-name that belongs to the function name.
# funname_blockname(obj)
# > obj: a Symbol or (SETF symbol)
# < result: Block-name, a Symbol
#define funname_blockname(obj)  \
  (atomp(obj) ? (object)obj : (object)Car(Cdr(obj)))

# UP: Determines, whether a Symbol is a Macro in the current Environment.
# sym_macrop(symbol)
# > symbol: Symbol
# < result: true if sym is a Symbol-Macro
extern bool sym_macrop (object sym);
# is used by CONTROL

/* UP: Sets the value of a Symbol in the current Environment.
 setq(symbol,value);
 > symbol: Symbol, not a constant
 > value: desired value of the Symbol in the current Environment
 < result: value
 can trigger GC */
extern maygc object setq (object sym, object value);
/* used by CONTROL */

# UP: Gives the definition of the function for a Symbol in an Environment
# sym_function(sym,fenv)
# > sym: name of the function (a Symbol for example)
# > fenv: a Functions- and Macrobindings-Environment
# < result: Definition of the function, either unbound (if function is undefined)
#             or Closure/SUBR/FSUBR/Macro/FunctionMacro.
extern object sym_function (object sym, object fenv);
# is used by CONTROL

# UP: "nests" an FUN-Environment, ie. writes all active bindings
# from the Stack to freshly allocated Vectors..
# nest_fun(env)
# > env: FUN-Env
# < result: same Environment, no pointer into the Stack
# can trigger GC
extern maygc object nest_fun (object env);
# is used by CONTROL

# UP: Nests the Environments in *env (ie. write all information
# to Stack-independent structures) and pushes it onto the STACK.
# nest_env(env)
# > gcv_environment* env: Pointer to five single Environments
# < gcv_environment* result: Pointer to the Environments on the STACK
# modifies STACK, can trigger GC
extern maygc gcv_environment_t* nest_env (gcv_environment_t* env);
# is used by Macro nest_aktenv

# UP: Nests the current environments (ie. writes all information
# to Stack-independent structures) and pushes them onto the STACK.
# (The values VAR_ENV, FUN_ENV, BLOCK_ENV, GO_ENV, DECL_ENV aren't
# modified, since there might be inactive bindings in frames that cannot
# be activated without modifying VAR_ENV .)
# nest_aktenv()
# < gcv_environment* result: Pointer to the Environments on the STACK
# modifies STACK, can trigger GC
# extern gcv_environment* nest_aktenv (void);
#define nest_aktenv()  nest_env(&aktenv)
# is used by CONTROL

# UP: Augments a Declarations-Environment with one decl-spec.
# augment_decl_env(declspec,env)
# > declspec: Declarations-Specifier, a Cons
# > env: Declarations-Environment
# < result: new (possibly augmented) Declarations-Environment
# can trigger GC
extern maygc object augment_decl_env (object new_declspec, object env);
# is used by CONTROL

# UP: expands a Form, if possible, (but not, if FSUBR-call
# or Symbol or FunctionMacro-call) in an Environment
# macroexp(form,venv,fenv);
# > form: Form
# > venv: a Variable- and Symbolmacro-Environment
# > fenv: a Function- and Macrobindings-Environment
# < value1: the expansion
# < value2: NIL, if not expanded,
#           T, if expanded
# can trigger GC
extern maygc void macroexp (object form, object venv, object fenv);
# is used by CONTROL

# UP: expands a Form if possible, (also, if FSUBR-call or
# Symbol, but not if FunctionMacro-call) in an environment
# macroexp0(form,env);
# > form: Form
# > env: a macro-expansion environment
# < value1: the expansion
# < value2: NIL, if not expanded,
#           T, if expanded
# can trigger GC
extern maygc void macroexp0 (object form, object env);
# is used by CONTROL

# UP: Parse-Declarations-Docstring. Detaches from a Form-list those,
# that must be viewed as Declarations resp. Documentation-string.
# parse_dd(formlist)
# > formlist: ( {decl|doc-string} . body )
# < value1: body
# < value2: list of decl-specs
# < value3: Doc-String or NIL
# < result: true if an (COMPILE)-Declaration occurred, else false
# can trigger GC
extern maygc bool parse_dd (object formlist);
# is used by CONTROL

# UP: Creates a corresponding Closure for a Lambda-body by disassembling
# the Lambda-list and possibly macro-expanding of all forms.
# get_closure(lambdabody,name,blockp,env)
# > lambdabody: (lambda-list {decl|doc} {form})
# > name: Name, a Symbol or (SETF symbol)
# > blockp: whether an implicit BLOCK is to be added
# > env: Pointer to the five individual environments:
#        env->var_env = VENV, env->fun_env = FENV,
#        env->block_env = BENV, env->go_env = GENV,
#        env->decl_env = DENV.
# < result: Closure
# can trigger GC
extern maygc object get_closure (object lambdabody, object name, bool blockp, gcv_environment_t* env);
# is used by CONTROL, SYMBOL, PREDTYPE

# UP: Converts an argument to a function.
# coerce_function(obj)
# > obj: Object
# < result: Object as function (SUBR or Closure)
# can trigger GC
extern maygc object coerce_function (object obj);
# is used by IO, FOREIGN

# Binds a Symbol dynamically to a value.
# Creates a dynamic variable-bindings frame for 1 variable.
# dynamic_bind(var,val)
# > var: a Symbol
# > val: the new value
# decreases STACK by 3 entries
# modifies STACK
#define dynamic_bind(variable,val_to_use)      \
  do { var gcv_object_t* top_of_frame = STACK; \
    var object sym_to_bind = (variable);       \
    # Create frame :                           \
    pushSTACK(Symbol_value(sym_to_bind));      \
    pushSTACK(sym_to_bind);                    \
    finish_frame(DYNBIND);                     \
    # modify value                             \
    Symbol_value(sym_to_bind) = (val_to_use);  \
  } while(0)
# is used by IO, EVAL, DEBUG, ERROR

# Unbinds a dynamic variable-bindings frame for one variable.
# dynamic_unbind_g()  - generic
# dynamic_unbind(sym) - with an additional check when STACKCHECKB
# increases STACK by 3 entries
# modifies STACK
#if STACKCHECKB
  #define CHECK_DYNBIND                                                 \
    if (!((as_oint(STACK_0) & wbit(frame_bit_o)) && framecode(STACK_0))) \
      abort()
#else
  #define CHECK_DYNBIND
#endif
#define dynamic_unbind_g()    do {                                      \
  CHECK_DYNBIND;                                                        \
  Symbol_value(STACK_(1)) = STACK_(2); /* restore the value */          \
  skipSTACK(3);    /* dismantle Frame */                                \
} while(0)
#if STACKCHECKB
  #define dynamic_unbind(sym)   do {              \
    if (!eq(sym,STACK_1)) abort();                \
    dynamic_unbind_g();                           \
  } while(0)
#else
  #define dynamic_unbind(sym)   dynamic_unbind_g()
#endif
# is used by IO, DEBUG, ERROR, EVAL, PREDTYPE, SPVW

# Executes "implicit PROGN" .
# implicit_progn(body,default)
# Executes body as implicit PROGN.
#  If the body is empty, the value is the default one.
# can trigger GC
#define implicit_progn(body,default)                                   \
  do {                                                                 \
    var object rest = (body);                                          \
    if (atomp(rest)) {                                                 \
      VALUES1(default); # default as value                             \
    } else                                                             \
      do { pushSTACK(Cdr(rest)); eval(Car(rest)); rest = popSTACK(); } \
      while (consp(rest));                                             \
  } while(0)
# is used by EVAL, CONTROL

# Highest number of parameters in a lambda-list (< bitm(intCsize))
# (= value of LAMBDA-PARAMETERS-LIMIT - 1)
#define lp_limit_1  ((uintL)(bitm(12)-1))

# Highest number of arguments for a function call
# (= value of CALL-ARGUMENTS-LIMIT - 1)
#define ca_limit_1  lp_limit_1

# The macro LISPSPECFORM initiates the declaration of a LISP-Special-Form.
# LISPSPECFORM(name,req_anz,opt_anz,body_flag)
# > name: C-name of the function and the Symbol
# > req_anz: number of required parameters
# > opt_anz: number of optional parameters
# > body_flag: body or nobody, depending on whether &BODY exists or not
# See FSUBR.D
#define LISPSPECFORM  LISPSPECFORM_B
# is used by CONTROL

/* The macro LISPFUN initiates a declaration of a LISP functions.
 LISPFUN(name,seclass,req_anz,opt_anz,rest_flag,key_flag,key_anz,keywords)
 > name: the name of the function (a C-Identifier)
 > seclass: the side-effect class (seclass_t, see above)
 > req_anz: number of required parameters (a number)
 > opt_anz: number of optional parameters (a number)
 > rest_flag: either norest or rest, depending on whether &REST exists or not
 > key_flag: either nokey or key or key_allow, depending on whether &KEY
             exists or not and whether &ALLOW-OTHER-KEYS is present
 > key_anz: number of keyword-parameters, a number (0 if nokey)
 > keywords: either NIL or an expression of the form
             v(kw(keyword1),...,kw(keywordn))   (NIL if nokey)
 See SUBR.D */
#define LISPFUN  LISPFUN_B
/* is used by all modules */

/* The macro LISPFUNN initiates a simple declaration of a LISP-function.
 LISPFUNN(name,req_anz)
 > name: the function-name (a C-Identifier)
 > req_anz: the (fixed) number of arguments (a number)
 LISPFUNNF - ditto, but seclass_foldable instead of seclass_default
 LISPFUNNR - ditto, but seclass_read instead of seclass_default
 See SUBR.D
 is used by all modules */

/* UP: initialize hand-made compiled closures
 init_cclosures();
 can trigger GC */
extern maygc void init_cclosures (void);


# ##################### CTRLBIBL for CONTROL.D ############################# #

/* the variables declared special appear on the stack twice:
   with binding SPECDECL (added when processing declarations)
   and the actual value (added when processing bindings).
 here we activate the SPECDECL bindings */
#define specdecled_p(sym,ptr,nn) (nn>0 ? specdecled_(sym,ptr,nn) : NULL)
/* Find the SPECDECL binding for the symbol
 > spec_pointer & spec_anz are returned by make_variable_frame()
 < return the pointer to the flags (or symbol+flags)
 i.e., something suitable to SET_BIT,
 or NULL if no such binding is found */
extern gcv_object_t* specdecled_ (object symbol, gcv_object_t* spec_pointer,
                                  uintL spec_anz);
/* used by CONTROL, EVAL */

/* activate the SPECDECL binding if found */
#define activate_specdecl(sym,ptr,nn) do {                      \
  var gcv_object_t *spec = specdecled_p(sym,ptr,nn);            \
  if (spec)                                                     \
    *spec = SET_BIT(*spec,active_bit_o); /* activate binding */ \
 } while(0)

/* activate all SPECDECL declarations */
extern void activate_specdecls (gcv_object_t* spec_ptr, uintC spec_count);
/* used by CONTROL, EVAL */

# Error if a block has already been left.
# fehler_block_left(name);
# > name: Block-name
nonreturning_function(extern, fehler_block_left, (object name));
# is used by EVAL

/* convert the numeric side-effect class as stored in subr_t or cclosure_t
 to the object - CONS or NIL - as used in compiler.lisp and for #Y i/o */
#ifndef COMPILE_STANDALONE
static inline object seclass_object (seclass_t sec) {
  switch (sec) {
    case seclass_foldable: return NIL;
    case seclass_no_se:    return O(seclass_no_se);
    case seclass_read:     return O(seclass_read);
    case seclass_write:    return O(seclass_write);
    case seclass_default:  return O(seclass_default);
    default: NOTREACHED;
  }
}
#endif
/* used by IO and CONTROL */

# ########################## for ENCODING.D ################################ #

# Initialize the encodings.
# init_encodings();
extern void init_encodings_1 (void);
extern void init_encodings_2 (void);
# is used by SPVW

# Initialize the encodings which depend on environment variables.
# init_dependent_encodings();
extern void init_dependent_encodings (void);
# is used by SPVW

# Maximum number of bytes needed to form a character, over all encodings.
# max_bytes_per_chart
#ifdef UNICODE
  #define max_bytes_per_chart  8  # 6 for JAVA, 7 for ISO-2022-KR, 8 for ISO-2022-CN[-EXT]
#else
  #define max_bytes_per_chart  1
#endif
# is used by STREAM

# UP: Creates a LISP-String with given contents.
# n_char_to_string(charptr,len,encoding)
# > char* charptr: address of a character sequence
# > uintL len: length of the sequence
# > object encoding: Encoding
# < result: Normal-Simple-String with len characters starting from charptr as contents
# can trigger GC
#ifdef UNICODE
  extern maygc object n_char_to_string (const char* charptr, uintL len, object encoding);
#else
  #define n_char_to_string(charptr,len,encoding)  n_char_to_string_(charptr,len)
  extern maygc object n_char_to_string_ (const char* charptr, uintL len);
#endif
# is used by PATHNAME
%% #ifdef UNICODE
%%   puts("extern object n_char_to_string (const char* charptr, uintL len, object encoding);");
%% #else
%%   emit_define("n_char_to_string(charptr,len,encoding)","n_char_to_string_(charptr,len)");
%%   puts("extern object n_char_to_string_ (const char* charptr, uintL len);");
%% #endif

# UP: Converts an ASCIZ-String to a LISP-String.
# asciz_to_string(asciz,encoding)
# ascii_to_string(asciz)
# > char* asciz: ASCIZ-String
#       (address of a null-terminated character-sequence)
# > object encoding: Encoding
# < result: Normal-Simple-String with the character sequence (without null-byte) as contents.
# can trigger GC
#ifdef UNICODE
  extern maygc object asciz_to_string (const char * asciz, object encoding);
#else
  #define asciz_to_string(asciz,encoding)  asciz_to_string_(asciz)
  extern maygc object asciz_to_string_ (const char * asciz);
#endif
extern maygc object ascii_to_string (const char * asciz);
# is used by SPVW/CONSTSYM, STREAM, PATHNAME, PACKAGE, GRAPH
%% #ifdef UNICODE
%%   puts("extern object asciz_to_string (const char * asciz, object encoding);");
%% #else
%%   emit_define("asciz_to_string(asciz,encoding)","asciz_to_string_(asciz)");
%%   puts("extern object asciz_to_string_ (const char * asciz);");
%% #endif
%% puts("extern object ascii_to_string (const char * asciz);");

# UP: Converts a String to an ASCIZ-String.
# string_to_asciz(obj,encoding)
# > object obj: String
# > object encoding: Encoding
# < result: Simple-Bit-Vector with the same characters as bytes and one
#             additional null-byte at the end
# < TheAsciz(result): address of the byte-sequence contained in there
# can trigger GC
#ifdef UNICODE
  extern maygc object string_to_asciz (object obj, object encoding);
#else
  #define string_to_asciz(obj,encoding)  string_to_asciz_(obj)
  extern maygc object string_to_asciz_ (object obj);
#endif
#define TheAsciz(obj)  ((char*)(&TheSbvector(obj)->data[0]))
# is used by STREAM, PATHNAME
%% #ifdef UNICODE
%%   puts("extern object string_to_asciz (object obj, object encoding);");
%% #else
%%   export_def(string_to_asciz(obj,encoding));
%%   puts("extern object string_to_asciz_ (object obj);");
%% #endif
%% export_def(TheAsciz(obj));

# Converts a String to an ASCIZ-String on the C-Stack.
# with_string_0(string,encoding,asciz,statement);
# with_sstring_0(simple_string,encoding,asciz,statement);
# copies the contents of string (which should be a Lisp string) to a safe area
# (zero-terminating it), binds the variable asciz pointing to it, and
# executes the statement.
#if 0
  #define with_string_0(string,encoding,ascizvar,statement)  \
    do { var char* ascizvar = TheAsciz(string_to_asciz(string,encoding)); \
         statement                                                        \
    } while(0)
  #define with_sstring_0  with_string_0
#else
  #define with_string_0_help_(string,encoding,ascizvar,statement,ascizvar_len,ascizvar_offset,ascizvar_string,ascizvar_bytelen,ascizvar_data,A,NR) \
    do { var uintL ascizvar_len;                                        \
      var uintL ascizvar_offset;                                        \
      var object ascizvar_string = unpack_string_ro(string,&ascizvar_len,&ascizvar_offset); \
      var const chart* ptr1;                                            \
      unpack_sstring_alloca_help_(ascizvar_string,ascizvar_len,ascizvar_offset, ptr1=,NR); \
     {var uintL ascizvar_bytelen = cslen(encoding,ptr1,ascizvar_len);   \
      var DYNAMIC_ARRAY(ascizvar_data,uintB,ascizvar_bytelen+1);        \
      cstombs_help_(encoding,ptr1,ascizvar_len,&ascizvar_data[0],ascizvar_bytelen,A); \
      ascizvar_data[ascizvar_bytelen] = '\0';                           \
     {var char* ascizvar = (char*) &ascizvar_data[0];                   \
      statement}                                                        \
      FREE_DYNAMIC_ARRAY(ascizvar_data);                                \
    }} while(0)
  #define with_string_0(string,encoding,ascizvar,statement) \
    with_string_0_help_(string,encoding,ascizvar,statement,ascizvar##_len,ascizvar##_offset,ascizvar##_string,ascizvar##_bytelen,ascizvar##_data,ASSERT,NOTREACHED)
  #define with_sstring_0_help_(string,encoding,ascizvar,statement,ascizvar_len,ascizvar_string,ascizvar_bytelen,ascizvar_data,A,NR) \
    do { var object ascizvar_string = (string);                         \
      sstring_un_realloc(ascizvar_string);                              \
     {var uintL ascizvar_len = Sstring_length(ascizvar_string);         \
      var const chart* ptr1;                                            \
      unpack_sstring_alloca_help_(ascizvar_string,ascizvar_len,0, ptr1=,NR); \
     {var uintL ascizvar_bytelen = cslen(encoding,ptr1,ascizvar_len);   \
      var DYNAMIC_ARRAY(ascizvar_data,uintB,ascizvar_bytelen+1);        \
      cstombs_help_(encoding,ptr1,ascizvar_len,&ascizvar_data[0],ascizvar_bytelen,A); \
      ascizvar_data[ascizvar_bytelen] = '\0';                           \
     {var char* ascizvar = (char*) &ascizvar_data[0];                   \
      statement}                                                        \
      FREE_DYNAMIC_ARRAY(ascizvar_data);                                \
    }}} while(0)
  #define with_sstring_0(string,encoding,ascizvar,statement) \
    with_sstring_0_help_(string,encoding,ascizvar,statement,ascizvar##_len,ascizvar##_string,ascizvar##_bytelen,ascizvar##_data,ASSERT,NOTREACHED)
#endif
# is used by PATHNAME, MISC, FOREIGN
%% export_def(with_string_0_help_(string,encoding,ascizvar,statement,ascizvar_len,ascizvar_offset,ascizvar_string,ascizvar_bytelen,ascizvar_data,A,NR));
%% export_def(with_sstring_0_help_(string,encoding,ascizvar,statement,ascizvar_len,ascizvar_string,ascizvar_bytelen,ascizvar_data,A,NR));
%% /* cannot use emit_define because Rectype_* is not a define in lispbibl.d */
%% puts("#define with_string_0(string,encoding,ascizvar,statement) with_string_0_help_(string,encoding,ascizvar,statement,ascizvar##_len,ascizvar##_offset,ascizvar##_string,ascizvar##_bytelen,ascizvar##_data,ASSERT,NOTREACHED)");
%% puts("#define with_sstring_0(string,encoding,ascizvar,statement) with_sstring_0_help_(string,encoding,ascizvar,statement,ascizvar##_len,ascizvar##_string,ascizvar##_bytelen,ascizvar##_data,ASSERT,NOTREACHED)");

# In some foreign modules, we call library functions that can do callbacks.
# When we pass a parameter to such a library function, maybe it first does a
# callback - which may involve garbage collection - and only then looks at
# the parameter. Therefore all the parameters, especially strings, must be
# located in areas that are not moved by garbage collection. The following
# macro helps achieving this.

# Converts a String to a String on the C-Stack.
# with_string(string,encoding,charptr,len,statement);
# with_sstring(simple_string,encoding,charptr,len,statement);
# copies the contents of string (which should be a Lisp string) to a safe area,
# binds the variable charptr pointing to it and the variable len to its length,
# and executes the statement.
#define with_string(string,encoding,charptrvar,lenvar,statement)  \
  do { var uintL charptrvar##_len;                                        \
    var uintL charptrvar##_offset;                                        \
    var object charptrvar##_string = unpack_string_ro(string,&charptrvar##_len,&charptrvar##_offset); \
    var const chart* ptr1;                                                \
    unpack_sstring_alloca(charptrvar##_string,charptrvar##_len,charptrvar##_offset, ptr1=); \
   {var uintL lenvar = cslen(encoding,ptr1,charptrvar##_len);             \
    var DYNAMIC_ARRAY(charptrvar##_data,uintB,lenvar);                    \
    cstombs(encoding,ptr1,charptrvar##_len,&charptrvar##_data[0],lenvar); \
    {var char* charptrvar = (char*) &charptrvar##_data[0];                \
     statement                                                            \
    }                                                                     \
    FREE_DYNAMIC_ARRAY(charptrvar##_data);                                \
  }} while(0)
#define with_sstring(string,encoding,charptrvar,lenvar,statement)  \
  do { var object charptrvar##_string = (string);                         \
    sstring_un_realloc(charptrvar##_string);                              \
   {var uintL charptrvar##_len = Sstring_length(charptrvar##_string);     \
    var const chart* ptr1;                                                \
    unpack_sstring_alloca(charptrvar##_string,charptrvar##_len,0, ptr1=); \
   {var uintL lenvar = cslen(encoding,ptr1,charptrvar##_len);             \
    var DYNAMIC_ARRAY(charptrvar##_data,uintB,lenvar);                    \
    cstombs(encoding,ptr1,charptrvar##_len,&charptrvar##_data[0],lenvar); \
    {var char* charptrvar = (char*) &charptrvar##_data[0];                \
     statement                                                            \
    }                                                                     \
    FREE_DYNAMIC_ARRAY(charptrvar##_data);                                \
  }}} while(0)
# is used by PATHNAME

/* Error, when a character cannot be converted to an encoding.
 fehler_unencodable(encoding,ch); */
nonreturning_function(extern, fehler_unencodable, (object encoding, chart ch));
# is used by STREAM

# ####################### ARRBIBL for ARRAY.D ############################## #

# ARRAY-TOTAL-SIZE-LIMIT is chosen as large as possible, respecting the
# constraint that the total-size of any array is a fixnum and (from ANSI CL)
# that ARRAY-TOTAL-SIZE-LIMIT itself is also a fixnum.
# (>=0, <2^oint_data_len):
#if (oint_data_len<=intLsize)
  #define arraysize_limit_1  ((uintV)(vbitm(oint_data_len)-2))
#else
  # Respect the constraint that the total-size of any array is an uintL.
  #define arraysize_limit_1  ((uintV)(vbitm(intLsize)-1))
#endif

/* ARRAY-RANK-LIMIT is chosen as large as possible, respecting the constraint
 that the rank of any array is an uintWC:
  #define arrayrank_limit_1  ((uintL)(bitm(intWCsize)-1))
 array dimensions are pushed on STACK in array_dimensions()
 so we are limited like with LAMBDA-PARAMETERS-LIMIT */
#define arrayrank_limit_1  lp_limit_1

# Macro: Follows the Sistring chain, to get from a simple array (actually,
# a string) to its storage vector.
# sstring_un_realloc(array);
# sstring_un_realloc1(array); [when at most one Sistring is involved]
# > array: a simple array
# < array: its storage vector
#ifdef HAVE_SMALL_SSTRING
  #ifdef TYPECODES
    #define sarray_reallocstringp(array)  \
      (typecode(array) == sstring_type                                  \
       && (sstring_flags(TheSstring(array)) & sstringflags_forwarded_B) \
      )
  #else
    #define sarray_reallocstringp(array)  \
      (Record_type(array) == Rectype_reallocstring)
  #endif
  #define sstring_un_realloc(array)  \
    while (sarray_reallocstringp(array)) \
      (array) = TheSistring(array)->data/*;*/
  #define sstring_un_realloc1(array)  \
    if (sarray_reallocstringp(array)) \
      (array) = TheSistring(array)->data/*;*/
#else
  #define sstring_un_realloc(array)  (void)0 /*nop*/
  #define sstring_un_realloc1(array)  (void)0 /*nop*/
#endif

# Function: Copies a simple-vector.
# copy_svector(vector)
# > vector: simple-vector
# < result: fresh simple-vector with the same contents
# can trigger GC
extern maygc object copy_svector (object vector);
# used by IO

# Function: Copies a simple-bit/byte-vector.
# copy_sbvector(vector)
# > vector: simple-bit/byte-vector
# < result: fresh simple-bit/byte-vector with the same contents
# can trigger GC
extern maygc object copy_sbvector (object vector);
# used by RECORD

# Function: Returns the active length of a vector (same as LENGTH).
# vector_length(vector)
# > vector: a vector
# < result: its length
extern uintL vector_length (object vector);
# used by many modules
%% puts("extern uintL vector_length (object vector);");

# Function: Canonicalizes an array element-type and returns its
# element type code.
# eltype_code(element_type)
# > element_type: type specifier
# < result: element type code Atype_xxx
# The canonicalized types are the possible results of ARRAY-ELEMENT-TYPE
# (symbols T, BIT, CHARACTER and lists (UNSIGNED-BYTE n)).
# The result type is a supertype of element_type.
# can trigger GC
extern maygc uintB eltype_code (object element_type);
# is used by SEQUENCE

# Function: Creates a simple-vector with given elements.
# vectorof(len)
# > uintC len: desired vector length
# > STACK_(len-1), ..., STACK_(0): len objects
# < result: simple-vector containing these objects
# Pops n objects off STACK.
# can trigger GC
extern maygc object vectorof (uintC len);
# used by PREDTYPE
%% puts("extern object vectorof (uintC len);");

# Function: For an indirect array, returns the storage vector and the offset.
# Also verifies that all elements of the array are physically present.
# iarray_displace_check(array,size,&index)
# > object array: indirect array
# > uintL size: size
# < result: storage vector
# < index: is incremented by the offset into the storage vector
extern object iarray_displace_check (object array, uintL size, uintL* index);
# used by IO, CHARSTRG, HASHTABL, PREDTYPE, STREAM, SEQUENCE, AFFI

# Function: For an array, returns the storage vector and the offset.
# Also verifies that all elements of the array are physically present.
# array_displace_check(array,size,&index)
# > object array: array
# > uintV size: size
# < result: storage vector
# < index: is incremented by the offset into the storage vector
extern object array_displace_check (object array, uintV size, uintL* index);
# used by HASHTABL, PREDTYPE, IO, FOREIGN
%% puts("extern object array_displace_check (object array, uintV size, uintL* index);");

# Tests for the storage vector of an array of element type NIL.
# simple_nilarray_p(obj)
#define simple_nilarray_p(obj)  nullp(obj)
%% export_def(simple_nilarray_p(obj));

# error-message
# > array: array (usually a Vector)
# > STACK_0: (erroneous) Index
nonreturning_function(extern, fehler_index_range, (object array, uintL bound));
# used by SEQUENCE

# error message: attempt to retrieve a value from (ARRAY NIL)
nonreturning_function(extern, fehler_nilarray_retrieve, (void));
# used by PREDTYPE
%% puts("nonreturning_function(extern, fehler_nilarray_retrieve, (void));");

# error message: attempt to store a value in (ARRAY NIL)
nonreturning_function(extern, fehler_nilarray_store, (void));

# error message: attempt to access a value from (ARRAY NIL)
nonreturning_function(extern, fehler_nilarray_access, (void));

# Function: Performs an AREF access.
# storagevector_aref(storagevector,index)
# > storagevector: a storage vector (simple vector or semi-simple byte vector)
# > index: (already checked) index into the storage vector
# < result: (AREF storagevector index)
# can trigger GC - if the element type is (UNSIGNED-BYTE 32)
extern /*maygc*/ object storagevector_aref (object storagevector, uintL index);
# used by IO

# Error when attempting to store an invalid value in an array.
# fehler_store(array,value);
nonreturning_function(extern, fehler_store, (object array, object value));
# used by SEQUENCE

# Macro: Tests a bit in a simple-bit-vector.
# if (sbvector_btst(sbvector,index)) ...
# > sbvector: a simple-bit-vector
# > index: index (a variable, must be < (length sbvector))
#define sbvector_btst(sbvector,index)  \
  ( # in byte number (index div 8), the bit number 7 - (index mod 8) : \
   TheSbvector(sbvector)->data[(uintL)(index)/8]                       \
     & bit((~(uintL)(index)) % 8)                                      \
  )
# used by SEQUENCE, IO

# Macro: Clears a bit in a simple-bit-vector.
# sbvector_bclr(sbvector,index);
# > sbvector: a simple-bit-vector
# > index: index (a variable, must be < (length sbvector))
#define sbvector_bclr(sbvector,index)  \
  ( # in byte number (index div 8), the bit number 7 - (index mod 8) : \
    TheSbvector(sbvector)->data[(uintL)(index)/8]                      \
      &= ~bit((~(uintL)(index)) % 8)                                   \
  )
# used by IO

# Macro: Sets a bit in a simple-bit-vector.
# sbvector_bset(sbvector,index);
# > sbvector: a simple-bit-vector
# > index: index (a variable, must be < (length sbvector))
#define sbvector_bset(sbvector,index)  \
  ( # in byte number (index div 8), the bit number 7 - (index mod 8) : \
    TheSbvector(sbvector)->data[(uintL)(index)/8]                      \
      |= bit((~(uintL)(index)) % 8)                                    \
  )
# used by SEQUENCE, IO

/* return Atype for the given array */
global uintBWL array_atype (object array);
/* used by socket.d and modules */
%% puts("extern uintBWL array_atype (object array);");

# Function: Returns the element-type of an array.
# array_element_type(array)
# > array: an array
# < result: element-type, one of the symbols T, BIT, CHARACTER, or a list
# can trigger GC
extern maygc object array_element_type (object array);
# used by PREDTYPE, IO

/* Returns the rank of an array.
 array_rank(array)
 > array: an array
 < uintL result: its rank = number of dimensions */
extern uintL array_rank (object array);
# used by modules
%% puts("extern uintL array_rank (object array);");

/* Returns the dimensions of an array.
 get_array_dimensions(array,rank,&dimensions[]);
 > array: an array
 > uintL rank: = array_rank(array)
 > uintL dimensions[0..rank-1]: room for rank dimensions
 < uintL dimensions[0..rank-1]: the array's dimensions */
extern void get_array_dimensions (object array, uintL rank, uintL* dimensions);
# used by modules
%% puts("extern void get_array_dimensions (object array, uintL rank, uintL* dimensions);");

# Function: Returns the list of dimensions of an array.
# array_dimensions(array)
# > array: an array
# < result: list of its dimensions
# can trigger GC
extern maygc object array_dimensions (object array);
# used by PREDTYPE, IO

# Function: Returns the dimensions of an array and their partial products.
# iarray_dims_sizes(array,&dims_sizes);
# > array: indirect array of rank r
# > struct { uintL dim; uintL dimprod; } dims_sizes[r]: room for the result
# < for i=1,...r:  dims_sizes[r-i] = { Dim_i, Dim_i * ... * Dim_r }
typedef struct { uintL dim; uintL dimprod; }  array_dim_size_t;
extern void iarray_dims_sizes (object array, array_dim_size_t* dims_sizes);
# used by IO

# Function: Returns the total-size of an array.
# array_total_size(array)
# > array: an array (a variable)
# < uintL result: its total-size
#ifndef COMPILE_STANDALONE
static inline uintL array_total_size (object array) {
  if (array_simplep(array)) {
    sstring_un_realloc(array);
    if (simple_string_p(array))
      return Sstring_length(array); # simple string: total length
    else
      return Sarray_length(array); # simple vector: total length
  } else
    return TheIarray(array)->totalsize; # indirect array: contains totalsize
}
#endif
# used by ARRAY, SEQUENCE, FOREIGN

# Function: Compares two slices of simple-bit-vectors.
# bit_compare(array1,index1,array2,index2,count)
# > array1: first simple-bit-vector
# > index1: absolute index into array1
# > array2: second simple-bit-vector
# > index2: absolute index into array2
# > count: number of bits to be compared, > 0
# < result: true, if both slices are the same, bit for bit, else false.
extern bool bit_compare (object array1, uintL index1,
                         object array2, uintL index2,
                         uintL bitcount);
# used by PREDTYPE

# Function: Copies a slice of an array array1 into another array array2.
# elt_copy(dv1,index1,dv2,index2,count);
# > dv1: source storage-vector
# > index1: start index in dv1
# > dv2: destination storage-vector
# > index2: start index in dv2
# > count: number of elements to be copied, > 0
# can trigger GC - if dv1 and dv2 have different element types or
#                  if both are strings and dv1 is wider than dv2
extern /*maygc*/ void elt_copy (object dv1, uintL index1, object dv2, uintL index2, uintL count);
# used by SEQUENCE, STREAM

# Function: Copies a slice of an array array1 into another array array2 of
# the same element type. Handles overlapping arrays correctly.
# elt_move(dv1,index1,dv2,index2,count);
# > dv1: source storage-vector
# > index1: start index in dv1
# > dv2: destination storage-vector
# > index2: start index in dv2
# > count: number of elements to be copied, > 0
# can trigger GC - if both are strings and dv1 is wider than dv2
extern /*maygc*/ void elt_move (object dv1, uintL index1, object dv2, uintL index2, uintL count);
# used by SEQUENCE

# Function: Fills a slice of an array with an element.
# elt_fill(dv,index,count,element)
# > dv: destination storage-vector
# > index: start index in dv
# > count: number of elements to be filled
# < result: true if element does not fit, false when done
# can trigger GC
extern maygc bool elt_fill (object dv, uintL index, uintL count, object element);
# used by SEQUENCE

# Function: Reverses a slice of an array, copying it into another array
# of the same element type.
# elt_reverse(dv1,index1,dv2,index2,count);
# > dv1: source storage-vector
# > index1: start index in dv1
# > dv2: destination storage-vector
# > index2: start index in dv2
# > count: number of elements to be copied, > 0
# can trigger GC
extern maygc void elt_reverse (object dv1, uintL index1, object dv2, uintL index2, uintL count);
# used by SEQUENCE

# Function: Reverses a slice of an array destructively.
# elt_nreverse(dv,index,count);
# > dv: storage-vector
# > index: start index in dv
# > count: number of elements to be reversed, > 0
extern void elt_nreverse (object dv, uintL index, uintL count);
# used by SEQUENCE

# Function: Tests whether an array has a fill-pointer.
# array_has_fill_pointer_p(array)
# > array: ein Array
# < result: true, if it has a fill-pointer, else false.
extern bool array_has_fill_pointer_p (object array);
# used by SEQUENCE, STREAM, IO

# Function: Allocates a new simple-bit-vector, filled with zeroes.
# allocate_bit_vector_0(len)
# > uintL len: length of the desired bit-vector (number of bits)
# < result: fresh simple-bit-vector, filled with zeroes
# can trigger GC
extern maygc object allocate_bit_vector_0 (uintL len);
# used by SEQUENCE
%% #if notused
%% puts("extern object allocate_bit_vector_0 (uintL len);");
%% #endif

# The following functions work on "semi-simple string"s.
# That are CHARACTER arrays with FILL-POINTER, (pro forma) not adjustable and
# not displaced, whose storagevector is a normal-simple-string. When their
# length is exceeded, the length is doubled (so that the resizing effort
# becomes unimportant: adding a character is still O(1) on average.)

# Function: Returns a fresh semi-simple-string of given length, with
# fill-pointer = 0.
# make_ssstring(len)
# > uintL len: desired length, must be >0
# < result: fresh semi-simple-string of the given length
# can trigger GC
extern maygc object make_ssstring (uintL len);
#define SEMI_SIMPLE_DEFAULT_SIZE 50
# used by STREAM, IO

# Function: Adds a character to a semi-simple-string, thereby possibly
# extending it.
# ssstring_push_extend(ssstring,ch)
# > ssstring: a semi-simple-string
# > ch: a character
# < result: the same semi-simple-string
# can trigger GC
extern maygc object ssstring_push_extend (object ssstring, chart ch);
# used by STREAM, IO

# Function: Ensures that a semi-simple-string has at least a given length,
# possibly extending it.
# ssstring_extend(ssstring,size)
# > ssstring: a semi-simple-string
# > size: desired minimum length
# < result: the same semi-simple-string
# can trigger GC
extern maygc object ssstring_extend (object ssstring, uintL needed_len);
# used by STREAM

# Function: Adds a substring to a semi-simple-string, thereby possibly
# extending it.
# ssstring_append_extend(ssstring,srcstring,start,len)
# > ssstring: a semi-simple-string
# > srcstring: a simple-string
# > start: the start index into the sstring
# > len: the number of characters to be pushed, starting from start
# < result: the same semi-simple-string
# can trigger GC
extern maygc object ssstring_append_extend (object ssstring, object srcstring, uintL start, uintL len);
# used by STREAM

# The following functions work on "semi-simple byte-vector"s.
# That are bit vectors with FILL-POINTER, (pro forma) not adjustable and
# not displaced, whose storagevector is a simple-bit-vector. When their
# length is exceeded, the length is doubled (so that the resizing effort
# becomes unimportant: adding a character is still O(1) on average.)

# Function: Returns a fresh semi-simple byte-vector of given length, with
# fill-pointer = 0.
# make_ssbvector(len)
# > uintL len: length (number of bytes!), must be >0
# < fresh: fresh semi-simple byte-vector of the given length
# can trigger GC
extern maygc object make_ssbvector (uintL len);
# used by IO

# Function: Adds a byte to a semi-simple byte vector, thereby possibly
# extending it.
# ssbvector_push_extend(ssbvector,b)
# > ssbvector: a semi-simple byte-vector
# > b: byte
# < result: the same semi-simple byte-vector
# can trigger GC
extern maygc object ssbvector_push_extend (object ssbvector, uintB b);
# used by IO

# ##################### CHARBIBL for CHARSTRG.D ############################ #

# Special Characters: (refer to above)
# #define BEL   7  #  #\Bell
# #define BS    8  #  #\Backspace
# #define TAB   9  #  #\Tab
# #define LF   10  #  #\Linefeed
# #define CR   13  #  #\Return
# #define PG   12  #  #\Page
#define NL   10  #  #\Newline
#define NLstring  "\n"  # C-String, that contains #\Newline
#define ESC  27  #  #\Escape
#define ESCstring  "\033"  # C-String, that contains #\Escape

# Converts Byte ch to upcase
# up_case(ch)
extern chart up_case (chart ch);
# is used by IO, PREDTYPE, PATHNAME
%% #if notused
%% puts("extern chart up_case (chart ch);");
%% #endif

# Converts Byte ch to downcase
# down_case(ch)
extern chart down_case (chart ch);
# is used by IO, PATHNAME
%% #if notused
%% puts("extern chart down_case (chart ch);");
%% #endif

# Checks whether a Character is alphanumeric.
# alphanumericp(ch)
# > ch: Character-Code
# < result: true if alphanumeric, else false.
extern bool alphanumericp (chart ch);
# is used by IO, PATHNAME

# Checks, whether a Character is a Graphic-Character ("printing").
# graphic_char_p(ch)
# > ch: Character-Code
# < result: true if printing, else false.
extern bool graphic_char_p (chart ch);
# is used by STREAM, PATHNAME

# Returns the screen display width of a character.
# char_width(ch)
# > ch: character code
# < result: number of output columns occupied by ch
extern uintL char_width (chart ch);
# is used by IO, STREAM

#if !defined(UNICODE) || defined(HAVE_SMALL_SSTRING)
# Copies an array of uint8 to an array of uint8.
# copy_8bit_8bit(src,dest,len);
# > uint8* src: source
# > uint8* dest: destination
# > uintL len: number of elements to be copied, > 0
extern void copy_8bit_8bit (const uint8* src, uint8* dest, uintL len);
#endif

#if defined(HAVE_SMALL_SSTRING)
# Copies an array of uint8 to an array of uint16.
# copy_8bit_16bit(src,dest,len);
# > uint8* src: source
# > uint16* dest: destination
# > uintL len: number of elements to be copied, > 0
extern void copy_8bit_16bit (const uint8* src, uint16* dest, uintL len);
#endif
%% #ifdef HAVE_SMALL_SSTRING
%%   puts("extern void copy_8bit_16bit (const uint8* src, uint16* dest, uintL len);");
%% #endif

#if defined(HAVE_SMALL_SSTRING)
# Copies an array of uint8 to an array of uint32.
# copy_8bit_32bit(src,dest,len);
# > uint8* src: source
# > uint32* dest: destination
# > uintL len: number of elements to be copied, > 0
extern void copy_8bit_32bit (const uint8* src, uint32* dest, uintL len);
#endif
%% #ifdef HAVE_SMALL_SSTRING
%%   puts("extern void copy_8bit_32bit (const uint8* src, uint32* dest, uintL len);");
%% #endif

#if defined(HAVE_SMALL_SSTRING)
# Copies an array of uint16 to an array of uint8.
# All source elements must fit into uint8.
# copy_16bit_8bit(src,dest,len);
# > uint16* src: source
# > uint8* dest: destination
# > uintL len: number of elements to be copied, > 0
extern void copy_16bit_8bit (const uint16* src, uint8* dest, uintL len);
#endif
%% #ifdef HAVE_SMALL_SSTRING
%%   puts("extern void copy_16bit_8bit (const uint16* src, uint8* dest, uintL len);");
%% #endif

#if defined(HAVE_SMALL_SSTRING)
# Copies an array of uint16 to an array of uint16.
# copy_16bit_16bit(src,dest,len);
# > uint16* src: source
# > uint16* dest: destination
# > uintL len: number of elements to be copied, > 0
extern void copy_16bit_16bit (const uint16* src, uint16* dest, uintL len);
#endif
%% #ifdef HAVE_SMALL_SSTRING
%%   puts("extern void copy_16bit_16bit (const uint16* src, uint16* dest, uintL len);");
%% #endif

#if defined(HAVE_SMALL_SSTRING)
# Copies an array of uint16 to an array of uint32.
# copy_16bit_32bit(src,dest,len);
# > uint16* src: source
# > uint32* dest: destination
# > uintL len: number of elements to be copied, > 0
extern void copy_16bit_32bit (const uint16* src, uint32* dest, uintL len);
#endif
%% #ifdef HAVE_SMALL_SSTRING
%%   puts("extern void copy_16bit_32bit (const uint16* src, uint32* dest, uintL len);");
%% #endif

#if defined(HAVE_SMALL_SSTRING)
# Copies an array of uint32 to an array of uint8.
# All source elements must fit into uint8.
# copy_32bit_8bit(src,dest,len);
# > uint32* src: source
# > uint8* dest: destination
# > uintL len: number of elements to be copied, > 0
extern void copy_32bit_8bit (const uint32* src, uint8* dest, uintL len);
#endif
%% #ifdef HAVE_SMALL_SSTRING
%%   puts("extern void copy_32bit_8bit (const uint32* src, uint8* dest, uintL len);");
%% #endif

#if defined(HAVE_SMALL_SSTRING)
# Copies an array of uint32 to an array of uint16.
# All source elements must fit into uint16.
# copy_32bit_16bit(src,dest,len);
# > uint32* src: source
# > uint16* dest: destination
# > uintL len: number of elements to be copied, > 0
extern void copy_32bit_16bit (const uint32* src, uint16* dest, uintL len);
#endif
%% #ifdef HAVE_SMALL_SSTRING
%%   puts("extern void copy_32bit_16bit (const uint32* src, uint16* dest, uintL len);");
%% #endif

#if defined(UNICODE)
# Copies an array of uint32 to an array of uint32.
# copy_32bit_32bit(src,dest,len);
# > uint32* src: source
# > uint32* dest: destination
# > uintL len: number of elements to be copied, > 0
extern void copy_32bit_32bit (const uint32* src, uint32* dest, uintL len);
#endif

#if defined(HAVE_SMALL_SSTRING)

# Determines the smallest string element type capable of holding a
# set of 8-bit characters.
#define smallest_string_flavour8(src,len)  \
  (unused (src), unused(len), Sstringtype_8Bit)

# Determines the smallest string element type capable of holding a
# set of 16-bit characters.
# smallest_string_flavour16(src,len)
# > uint16* src: source
# > uintL len: number of characters at src
# < result: Sstringtype_8Bit or Sstringtype_16Bit
extern uintBWL smallest_string_flavour16 (const uint16* src, uintL len);

# Determines the smallest string element type capable of holding a
# set of 32-bit characters.
# smallest_string_flavour32(src,len)
# > uint32* src: source
# > uintL len: number of characters at src
# < result: Sstringtype_8Bit or Sstringtype_16Bit or Sstringtype_32Bit
extern uintBWL smallest_string_flavour32 (const uint32* src, uintL len);

# Determines the smallest string element type capable of holding a
# set of characters.
# smallest_string_flavour(src,len)
# > chart* src: source
# > uintL len: number of characters at src
# < result: Sstringtype_8Bit or Sstringtype_16Bit or Sstringtype_32Bit
#ifndef COMPILE_STANDALONE
static inline uintBWL smallest_string_flavour (const chart* src, uintL len) {
  return smallest_string_flavour32((const uint32*)src,len);
}
#endif

#endif

# Dispatches among S8string, S16string, S32string and nilvector.
# SstringCase(string,s8string_statement,s16string_statement,s32string_statement,nilvector_statement);
# > string: a not-reallocated simple-string or simple-nilvector (i.e. NIL)
# Executes one of the three statement, depending on the element size of string.
#ifdef UNICODE
  #ifdef HAVE_SMALL_SSTRING
    #define SstringCase(string,s8string_statement,s16string_statement,s32string_statement,nilvector_statement)  \
      if (Array_type(string) == Array_type_snilvector) { nilvector_statement } else \
      if (sstring_eltype(TheSstring(string)) == Sstringtype_8Bit) { s8string_statement } else \
      if (sstring_eltype(TheSstring(string)) == Sstringtype_16Bit) { s16string_statement } else \
      if (sstring_eltype(TheSstring(string)) == Sstringtype_32Bit) { s32string_statement } else \
      NOTREACHED;
  #else
    #define SstringCase(string,s8string_statement,s16string_statement,s32string_statement,nilvector_statement)  \
      if (Array_type(string) == Array_type_snilvector) { nilvector_statement } else \
      { s32string_statement }
  #endif
#else
  # In this case we take the s32string_statement, not the s8string_statement,
  # because the s32string_statement is the right one for normal simple strings.
  #define SstringCase(string,s8string_statement,s16string_statement,s32string_statement,nilvector_statement)  \
    if (Array_type(string) == Array_type_snilvector) { nilvector_statement } else \
    { /*s8string_statement*/ s32string_statement }
#endif
# is used by CHARSTRG, ARRAY, HASHTABL, PACKAGE, PATHNAME, PREDTYPE, STREAM

# Dispatches among S8string, S16string, S32string and nilvector.
# SstringDispatch(string,suffix,statement)
# > string: a not-reallocated simple-string or simple-nilvector (i.e. NIL)
# Executes the statement with cint##suffix being bound to the appropriate
# integer type (cint8, cint16 or cint32) and with Sstring being bound to the
# appropriate struct pointer type (S8string, S16string or S32string).
# Gives an error for simple-nilvector; must therefore only be called if the
# contents of the string is really to be accessed.
#define SstringDispatch(string,suffix,statement)  \
  SstringCase(string,                                                 \
    { typedef cint8 cint##suffix; typedef S8string Sstring##suffix;   \
      statement                                                       \
    },                                                                \
    { typedef cint16 cint##suffix; typedef S16string Sstring##suffix; \
      statement                                                       \
    },                                                                \
    { typedef cint32 cint##suffix; typedef S32string Sstring##suffix; \
      statement                                                       \
    },                                                                \
    { fehler_nilarray_access();                                       \
    })
# is used by CHARSTRG, ARRAY, HASHTABL, PACKAGE, PATHNAME, PREDTYPE, STREAM

# Tests whether a simple-string is a normal-simple-string.
# sstring_normal_p(string)
# > string: a not-reallocated simple-string
#ifdef HAVE_SMALL_SSTRING
  #define sstring_normal_p(string)  \
    (sstring_eltype(TheSstring(string)) == Sstringtype_32Bit)
#else
  #define sstring_normal_p(string)  1
#endif

# Makes a string contents available.
# unpack_sstring_alloca(string,len,offset, charptr = );
# > object string: a not-reallocated simple-string
# > uintL len: the number of characters to be accessed
# > uintL offset: where the characters to be accessed start
# < const chart* charptr: pointer to the characters
#   (may be in string, may be on the stack)
#ifdef HAVE_SMALL_SSTRING
  #define unpack_sstring_alloca_help_(string,len,offset,charptr_assignment,u) \
    if (simple_nilarray_p(string)) {                                           \
      if ((len) > 0) fehler_nilarray_retrieve();                               \
      charptr_assignment NULL;                                                 \
    } else if (sstring_eltype(TheSstring(string)) == Sstringtype_32Bit) {      \
      charptr_assignment (const chart*) &TheS32string(string)->data[offset];   \
    } else {                                                                   \
      var chart* _unpacked_ = (chart*)alloca((len)*sizeof(chart));             \
      if ((len) > 0) {                                                         \
        if (sstring_eltype(TheSstring(string)) == Sstringtype_16Bit)           \
          copy_16bit_32bit(&TheS16string(string)->data[offset],(cint32*)_unpacked_,len);\
        else if (sstring_eltype(TheSstring(string)) == Sstringtype_8Bit) \
          copy_8bit_32bit(&TheS8string(string)->data[offset],(cint32*)_unpacked_,len);\
        else                                                                   \
          u;                                                            \
      }                                                                        \
      charptr_assignment (const chart*) _unpacked_;                            \
    }
#else
  #define unpack_sstring_alloca_help_(string,len,offset,charptr_assignment,u) \
    if (simple_nilarray_p(string)) {                                           \
      if ((len) > 0) fehler_nilarray_retrieve();                               \
      charptr_assignment NULL;                                                 \
    } else {                                                                   \
      charptr_assignment (const chart*) &TheSnstring(string)->data[offset];    \
    }
#endif
#define unpack_sstring_alloca(string,len,offset,charptr_assignment)     \
  unpack_sstring_alloca_help_(string,len,offset,charptr_assignment,NOTREACHED)
# is used by
%% export_def(unpack_sstring_alloca_help_(string,len,offset,charptr_assignment,u));
%% puts("#define unpack_sstring_alloca(s,l,o,c) unpack_sstring_alloca_help_(s,l,o,c,NOTREACHED)");

# UP: Fetches a character from a simple string.
# schar(string,index)
# > object string: a not-reallocated simple-string or simple-nilvector (i.e. NIL)
# > uintL index: >= 0, < length of string
# < chart result: character at the given position
#ifndef COMPILE_STANDALONE
static inline chart schar (object string, uintL index) {
  SstringDispatch(string,X, {
    return as_chart(((SstringX)TheVarobject(string))->data[index]);
  });
  return as_chart(0); /* not reached - just pacify the compiler */
}
#endif
# is used by PATHNAME, STREAM

# UP: unpacks a String.
# unpack_string_rw(string,&len,&offset)  [for read-write access]
# > object string: a String.
# < uintL len: number of characters of the String.
# < uintL offset: offset in the datastorage vector
# < object result: datastorage vector, a simple-string or NIL
extern object unpack_string_rw (object string, uintL* len, uintL* offset);
# is used by AFFI

# UP: unpacks a String.
# unpack_string_ro(string,&len,&offset)  [for read-only access]
# > object string: a String.
# < uintL len: number of characters of the String.
# < uintL offset: offset into the datastorage vector
# < object result: datastorage vector, a simple-string or NIL
extern object unpack_string_ro (object string, uintL* len, uintL* offset);
# is used by STREAM, HASHTABL, PACKAGE, SEQUENCE, ENCODING
%% puts("extern object unpack_string_ro (object string, uintL* len, uintL* offset);");

# UP: tests two Strings for equality
# string_gleich(string1,string2)
# > string1: String
# > string2: simple-string
# < result: /=0, if equal
extern bool string_gleich (object string1, object string2);
# is used by PACKAGE, STREAM, IO

# UP: tests two Strings for equality, case-insensitive
# string_equal(string1,string2)
# > string1: String
# > string2: simple-string
# < result: /=0, if equal
extern bool string_equal (object string1, object string2);
# is used by IO, PATHNAME
%% puts("extern bool string_equal (object string1, object string2);");

# UP: Stores a character in a string.
# > string: a mutable string that is or was simple
# > index: (already checked) index into the string
# > element: a character
# < result: the possibly reallocated string
# can trigger GC
  extern maygc object sstring_store (object string, uintL index, chart element);
# is used by STREAM

# UP: Stores an array of characters in a string.
# > string: a mutable string that is or was simple
# > offset: (already checked) offset into the string
# > charptr[0..len-1]: a character array, not GC affected
# < result: the possibly reallocated string
# can trigger GC
  extern maygc object sstring_store_array (object string, uintL offset,
                                           const chart *charptr, uintL len);
# is used by AFFI

#ifdef UNICODE
# UP: Creates a Simple-String with given elements.
# stringof(len)
# > uintL len: desired length of vector
# > on STACK: len Characters, first one on top
# < result: Simple-String with these objects
# increases STACK
# modifies STACK, can trigger GC
  extern maygc object stringof (uintL len);
# is used by ENCODING, STREAM
#endif

# UP: copies a String and turns it into a Simple-String.
# copy_string_normal(string)
# > string: String
# < result: mutable Normal-Simple-String with the same characters
# can trigger GC
  extern maygc object copy_string_normal (object string);
# is used by IO, PATHNAME

# UP: copies a String and turns it into a Simple-String.
# copy_string(string)
# > string: String
# < result: mutable Simple-String with the same characters
# can trigger GC
#ifdef HAVE_SMALL_SSTRING
  extern maygc object copy_string (object string);
#else
  #define copy_string(string)  copy_string_normal(string)
#endif
# is used by IO, PATHNAME

# UP: Converts a String to a Simple-String.
# coerce_ss(obj)
# > obj: Lisp-Object, should be a String.
# < result: Simple-String with the same characters.
# can trigger GC
extern maygc object coerce_ss (object obj);
# is used by STREAM, PATHNAME

# UP: Converts a String to a immutable Simple-String.
# coerce_imm_ss(obj)
# > obj: Lisp-Object, should be a String.
# < result: immutable Simple-String with the same characters.
# can trigger GC
extern maygc object coerce_imm_ss (object obj);
# is used by PACKAGE

# UP: Converts a String to a Normal-Simple-String.
# coerce_normal_ss(obj)
# > obj: Lisp-Object, should be a String.
# < result: Normal-Simple-String with the same characters.
# can trigger GC
#ifndef HAVE_SMALL_SSTRING
  #define coerce_normal_ss coerce_ss
#else
  extern maygc object coerce_normal_ss (object obj);
#endif
# is used by PATHNAME

#if 0 # unused
# UP: converts a String to an immutable Normal-Simple-String.
# coerce_imm_normal_ss(obj)
# > obj: Lisp-Object, should be a String.
# < result: immutable Normal-Simple-String with the same characters
# can trigger GC
  #ifndef HAVE_SMALL_SSTRING
    #define coerce_imm_normal_ss coerce_imm_ss
  #else
    extern maygc object coerce_imm_normal_ss (object obj);
  #endif
# is used by
#endif

# UP: Converts Object to a Character
# coerce_char(obj)
# > obj: Lisp-Object
# < result: Character or NIL
extern object coerce_char (object obj);
# is used by PREDTYPE

# UP: Returns the name of a character
# char_name(code)
# > chart code: Code of a character
# < result: Simple-String (this char's name) or NIL
# can trigger GC
extern maygc object char_name (chart code);
# is used by IO

# UP: Determines the Character with a given Name
# name_char(string)
# > string: String
# < result: Character with this Name, or NIL if none exists
extern object name_char (object string);
# is used by IO

/* Converts a character to opposite case.
 invert_case(ch)
 > ch: a character
 < result: a character, either ch or up_case(ch) or down_case(ch)
 Note that always invert_case(invert_case(ch)) == ch. */
extern chart invert_case (chart ch);
/* is used by PACKAGE */

/* UP: compares two strings for equality modulo case-invert
 string_gleich_inverted(string1,string2)
 > string1: string
 > string2: simple-string
 < result: /=0, if equal modulo case-invert */
extern bool string_gleich_inverted (object string1, object string2);
/* is used by PACKAGE */

/* UP: converts a string to opposite case
 string_invertcase(string)
 > string: string
 < result: new normal-simple-string
 can trigger GC */
extern object string_invertcase (object string);
/* is used by SYMBOL, PACKAGE */

/* UP: tests the limits for a String argument
 test_string_limits_ro(&arg)  [for read-only access]
 > STACK_2: String-Argument
 > STACK_1: optional :start-Argument
 > STACK_0: optional :end-Argument
 < stringarg arg: description of the argument
 < result: String-Argument
 increases STACK by 3
 can trigger GC */
typedef struct stringarg {
  object string; /* data vector - not-reallocated simple-string or -array */
  uintL offset;                 /* offset into this string */
  uintL index;                  /* :start index */
  uintL len;                    /* :end - :start */
} stringarg;
%% puts("typedef struct stringarg { object string; uintL offset; uintL index; uintL len; } stringarg;");
extern maygc object test_string_limits_ro (stringarg* arg);
/* used by STREAM, PATHNAME, IO, ENCODING */

/* UP: checks :START and :END limits for a vector argument
 > STACK_1: optional :start-argument
 > STACK_0: optional :end-argument
 > stringarg arg: arg.string its data vector,
                  [arg.offset .. arg.offset+arg.len-1] the range within the
                  data vector corresponding to the entire vector-argument
 < stringarg arg: arg.string and arg.offset unchanged,
                  [arg.offset+arg.index .. arg.offset+arg.index+arg.len-1] the
                  range within the data vector corresponding to the selected
                  vector slice
 removes 2 elements from STACK */
extern void test_vector_limits (stringarg* arg);
/* used by ENCODING */
%% puts("extern void test_vector_limits (stringarg* arg);");
/* used by RAWSOCK */

/* UP: checks a string/symbol/character-argument
 test_stringsymchar_arg(obj,invert)
 > obj: argument
 > invert: whether to implicitly case-invert a symbol's printname
 < result: argument as string
 can trigger GC */
extern maygc object test_stringsymchar_arg (object obj, bool invert);
/* used by IO, PACKAGE */

# UP: tests two equally long strings for equality
# > string1,offset1: Chars in String1 start from here
# > string2,offset2: Chars in String2 start from here
# > len: number of chars in String1 and String2, > 0
# < result: true if equal, else false.
extern bool string_eqcomp (object string1, uintL offset1, object string2, uintL offset2, uintL len);
# is used by PREDTYPE

# UP: compares two equally long strings, case-insensitive
# > string1,offset1: Chars in String1 start from here
# > string2,offset2: Chars in String2 start from here
# > len: number of chars in String1 and String2, > 0
# < result: true if equal, else false.
extern bool string_eqcomp_ci (object string1, uintL offset1, object string2, uintL offset2, uintL len);
# is used by PREDTYPE

# UP: converts the Characters of a partial string to upcase
# nstring_upcase(dv,offset,len);
# > object dv: the character storage vector
# > uintL offset: index of first affected character
# > uintL len: number of affected characters
# can trigger GC
extern maygc void nstring_upcase (object dv, uintL offset, uintL len);
# is not used at this time (except in CHARSTRG, of course)

# UP: converts the Characters of a partial string to downcase
# nstring_downcase(dv,offset,len);
# > object dv: the character storage vector
# > uintL offset: index of first affected character
# > uintL len: number of affected characters
# can trigger GC
extern maygc void nstring_downcase (object dv, uintL offset, uintL len);
# is used by PATHNAME

# UP: changes the words of a part of a string so they start
# with capital letters and continue with lowercase ones
# nstring_capitalize(dv,offset,len);
# > object dv: the character storage vector
# > uintL offset: index of first affected character
# > uintL len: number of affected characters
# can trigger GC
extern maygc void nstring_capitalize (object dv, uintL offset, uintL len);
# is used by PATHNAME

# UP: converts a String to upcase
# string_upcase(string)
# > string: String
# < result: new Normal-Simple-String, in upcase
# can trigger GC
extern maygc object string_upcase (object string);
# is used by MISC, PATHNAME

# UP: converts a String to downcase
# string_downcase(string)
# > string: String
# < result: new Normal-Simple-String, in downcase
# can trigger GC
extern maygc object string_downcase (object string);
# is used by PATHNAME

# Returns a substring of a simple-string.
# subsstring(string,start,end)
# > object string: a simple-string
# > uintL start: start index
# > uintL end: end index
# with 0 <= start <= end <= Sstring_length(string)
# < object result: (subseq string start end),
#                  a freshly created normal-simple-string
# can trigger GC
extern maygc object subsstring (object string, uintL start, uintL end);
# is used by CHARSTRG, PATHNAME

# UP: Concatenates several Strings to one String.
# string_concat(argcount)
# > uintC argcount: Number of arguments
# > on STACK: the argument (should be Strings)
# < result: newly created string
# < STACK: cleaned
# can trigger GC
extern maygc object string_concat (uintC argcount);
# is used by PACKAGE, PATHNAME, DEBUG, SYMBOL
%% puts("extern object string_concat (uintC argcount);");

# ###################### DEBUGBIB for DEBUG.D ############################ #

# Starts the normal driver (Read-Eval-Print-Loop)
# driver();
extern void driver (void);
# is used by SPVW

# Starts a secondary driver (Read-Eval-Print-Loop)
# break_driver(continuable_p);
# > continuable_p == can be continued after the driver finishes
# can trigger GC
extern maygc void break_driver (bool continuable_p);
# is used by ERROR, EVAL

# ##################### HASHBIBL for HASHTABL.D ########################## #

# UP: Gets the hash of an object from a hash-table.
# gethash(obj,ht,allowgc)
# > obj: Object, as key
# > ht: hash-table
# > allowgc: whether GC is allowed during hash lookup
#            (should be true if the hash-table has a user-defined test)
# < result: corresponding value, if found, else nullobj
# can trigger GC - if allowgc is true
extern /*maygc*/ object gethash (object obj, object ht, bool allowgc);
# is used by EVAL, RECORD, PATHNAME, FOREIGN
%% puts("extern object gethash (object obj, object ht, bool allowgc);");

# UP: Locates a key in a hash-table and gives the older value.
# shifthash(ht,obj,value) == (SHIFTF (GETHASH obj ht) value)
# > ht: hash-table
# > obj: object
# > value: new value
# > allowgc: whether GC is allowed during hash lookup
#            (should be true if the hash-table has a user-defined test or
#             if the hash-table is not known to already contain a value for obj)
# < result: old value
# can trigger GC - if allowgc is true
extern /*maygc*/ object shifthash (object ht, object obj, object value, bool allowgc);
# is used by SEQUENCE, PATHNAME, FOREIGN

/* hash_table_weak_type(ht)
 > ht: hash-table
 < result: symbol NIL/:KEY/:VALUE/:KEY-AND-VALUE/:KEY-OR-VALUE */
extern object hash_table_weak_type (object ht);
/* used by PREDTYPE */

/* HASH-TABLE-TEST (EQ/EQL/EQUAL/EQUALP)
 > ht: hash-table
 < result: symbol EQ/EQL/EQUAL/EQUALP or cons (TEST . HASH)
 can trigger GC - for user-defined ht_test */
extern maygc object hash_table_test (object ht);
/* used by HASHTABL, IO */

# Macro: Runs through a Hash-Tabelle.
# map_hashtable(ht,key,value,statement)
# map_hashtable_nogc(ht,key,value,statement)
# > ht: Hash-Tabelle
# Calls 'statement', where key and value are a pair from the table.
# The first form is necessary, if the statement can trigger GC.
#define map_hashtable(ht,key,value,statement)                           \
  do {                                                                  \
    var object ht_from_map_hashtable = (ht);                            \
    var uintL index_from_map_hashtable =                                \
      3*posfixnum_to_V(TheHashtable(ht_from_map_hashtable)->ht_maxcount); \
    pushSTACK(TheHashtable(ht_from_map_hashtable)->ht_kvtable);         \
    loop {                                                              \
      if (index_from_map_hashtable==0) break;                           \
      index_from_map_hashtable -= 3;                                    \
      {var gcv_object_t* KVptr_from_map_hashtable =                     \
         &TheHashedAlist(STACK_0)->hal_data[index_from_map_hashtable];  \
       var object key = KVptr_from_map_hashtable[0];                    \
       if (boundp(key)) {                                               \
         var object value = KVptr_from_map_hashtable[1];                \
              statement;                                                \
      } }  }                                                            \
    skipSTACK(1);                                                       \
  } while(0)
#define map_hashtable_nogc(ht,key,value,statement)                      \
  do {                                                                  \
    var object ht_from_map_hashtable = (ht);                            \
    var uintL index_from_map_hashtable =                                \
      posfixnum_to_V(TheHashtable(ht_from_map_hashtable)->ht_maxcount); \
    var gcv_object_t* KVptr_from_map_hashtable =                        \
      &TheHashedAlist(TheHashtable(ht_from_map_hashtable)->ht_kvtable)->hal_data[3*index_from_map_hashtable]; \
    loop {                                                              \
      if (index_from_map_hashtable==0) break;                           \
      index_from_map_hashtable--; KVptr_from_map_hashtable -= 3;        \
      { var object key = KVptr_from_map_hashtable[0];                   \
        if (boundp(key)) {                                              \
          var object value = KVptr_from_map_hashtable[1];               \
          statement;                                                    \
    } } }                                                               \
  } while(0)
# is used by IO

# ######################### IOBIBL for IO.D ############################## #

# check a cint for being a whitespace
#define cint_white_p(c)   \
  ((c)==' ' || (c)=='\n' || (c)=='\r' || (c)=='\t' || (c)=='\v' || (c)=='\f')

# special Object, that indicates EOF
#define eof_value  make_system(0xE0FE0FUL)
# is used by IO, STREAM, DEBUG, SPVW

# aux. value to recognize certain Dots
#define dot_value  make_system(0xD0DD0DUL)
# is used by IO, SPVW

# UP: Initializes the reader.
# init_reader();
# can trigger GC
extern maygc void init_reader (void);
# is used by SPVW

# UP:
# (setf (strm-pphelp-strings *stream_)
#    (list* (make-Semi-Simple-String 50)
#           (cons nl_type *PRIN-INDENTATION*)
#           (strm-pphelp-strings *stream_)))
extern object cons_ssstring (const gcv_object_t* stream_, object nl_type);
# used by io.d and stream.d

# UP: Reads an object.
# stream_read(&stream,recursive-p,whitespace-p)
# > recursive-p: tells whether there's a recursive call of READ,
#                with error at EOF
# > whitespace-p: tells, whether whitespace is to be consumed afterwards
# > stream: Stream
# < stream: Stream
# < result: read object (eof_value at EOF, dot_value for single dot)
# can trigger GC
extern maygc object stream_read (const gcv_object_t* stream_, object recursive_p, object whitespace_p);
# is used by SPVW, DEBUG

# UP: Write a Simple-String to a Stream element by element.
# write_sstring(&stream,string);
# > string: Simple-String
# > stream: Stream
# < stream: Stream
# can trigger GC
extern maygc void write_sstring (const gcv_object_t* stream_, object string);
# is used by EVAL, DEBUG, ERROR, PACKAGE, SPVW

# UP: Writes a String to a Stream element by element.
# write_string(&stream,string);
# > string: String
# > stream: Stream
# < stream: Stream
# can trigger GC
extern maygc void write_string (const gcv_object_t* stream_, object string);
# is used by PACKAGE, DEBUG

# UP: Writes an object to a Stream.
# prin1(&stream,obj);
# > obj: Objekt
# > stream: Stream
# < stream: Stream
# can trigger GC
extern maygc void prin1 (const gcv_object_t* stream_, object obj);
# is used by EVAL, DEBUG, PACKAGE, ERROR, SPVW

# UP: Writes a Newline to a Stream.
# terpri(&stream);
# > stream: Stream
# < stream: Stream
# can trigger GC
# extern maygc void terpri (const gcv_object_t* stream_);
#define terpri(stream_)  write_ascii_char(stream_,NL)
# is used by IO, DEBUG, PACKAGE, ERROR, SPVW

# ####################### LISTBIBL for LIST.D ############################## #

# UP: Copies a list
# copy_list(list)
# > list: list
# < result: copy of the list
# can trigger GC
extern maygc object copy_list (object list);
# is used by PACKAGE
%% puts("extern object copy_list (object old_list);");

# UP: Reverses a list constructively.
# reverse(list)
# > list: list (x1 ... xm)
# < result: reversed list (xm ... x1)
# can trigger GC
extern maygc object reverse (object list);
# is used by SEQUENCE, PACKAGE, PATHNAME

/* UP: get the length of a list and the last atom
 llength1(obj,last)
 > obj: object
 < uintL result: length of obj, interpreted as list
 < last: the last atom
 Does not test for circular lists. */
extern uintL llength1 (object obj, object* last);
/* used in SEQUENCE */
#define llength(obj)  llength1(obj,NULL)
/* used by CONTROL, EVAL, RECORD, IO, PACKAGE, HASHTABL, STREAM */

# UP: Makes a list with exactly len elements
# make_list(len)
# > STACK_0: Initial value for the elements
# > uintL len: desired list length
# < result: list with len elements
# can trigger GC
extern maygc object make_list (uintL len);
# is used by
%% #if notused
%% puts("extern object make_list (uintL len);");
%% #endif

# UP: reverses a list destructively.
# nreverse(list)
# > list: list (x1 ... xm)
# < result: list (xm ... x1), EQ to the old one
extern object nreverse (object list);
# is used by SEQUENCE, EVAL, CONTROL, IO, PATHNAME, ERROR, DEBUG, PACKAGE
%% puts("extern object nreverse (object list);");

# UP: A0 := (nreconc A0 A1)
# nreconc(list,obj)
# > list: list
# > obj: object
# < result: (nreconc A0 A1)
extern object nreconc (object list, object obj);
# is used by SEQUENCE, IO, PATHNAME, CONTROL, DEBUG

# UP: Build (delete obj (the list list) :test #'EQ)
# deleteq(list,obj)
# Remove all elements that are EQ to obj from the list.
# > obj: element to be removed
# > list: list
# < result: modified list
extern object deleteq (object list, object obj);
# is used by PACKAGE, STREAM
%% puts("extern object deleteq (object list, object obj);");

/* UP: check whether OBJ ends a proper list
 endp(obj)
 > obj: object
 < result: true if obj is the list end NIL,
           false if obj is a Cons.
           error otherwise */
extern bool endp (object obj);
/* used by CONTROL */
%% puts("extern bool endp (object obj);");

/* Finds the length of a possibly circular or dotted list.
 list_length(list,&dotted)
 > list: an object
 < result: the length (integer >= 0, or NIL for circular lists)
 < dotted: if non-circular, the last atom, i.e., the indicator whether the list
           is dotted
 can trigger GC */
extern maygc object list_length (object list, object *dottedp);
/* used by SEQUENCE */

/* proper_list_p(obj)
   returns true if obj is a proper list, i.e. a list which is neither dotted
   nor circular, i.e. a list which ends in NIL. */
extern bool proper_list_p (object obj);
/* used by PREDTYPE */

# UP: Creates a list with given elements.
# listof(len)
# > uintC len: desired list length
# > auf STACK: len objects, first one on top
# < result: list of those objects
# Increases STACK
# modifies STACK, can trigger GC
extern maygc object listof (uintC len);
/* used by STREAM, PATHNAME, PACKAGE, ARRAY, EVAL, PREDTYPE, ERROR, SPVW */
%% puts("extern object listof (uintC len);");

# UP: find OBJ in LIS: (MEMBER OBJ LIS :TEST #'EQ)
extern object memq (const object obj, const object lis);
# used by RECORD
%% puts("extern object memq (const object obj, const object lis);");

# ####################### MISCBIBL for MISC.D ############################## #

#ifdef HAVE_ENVIRONMENT
# Modify the environment variables.
# clisp_setenv(name,value) sets the value of the environment variable `name'
# to `value' and returns 0. Returns -1 if not enough memory.
  extern int clisp_setenv (const char * name, const char * value);

#endif

/* for modules */
typedef struct { long c_const; gcv_object_t *l_const; } c_lisp_pair_t;
typedef struct {
  const c_lisp_pair_t *table;         /* C <--> Lisp */
  const unsigned int size;            /* table length */
  const long default_value; /* what to use when Lisp value is missing */
  const bool have_default_value_p;   /* use default_value? */
  const bool use_default_function_p; /* use L_to_I when there is no map entry */
  const char *name;                  /* map name for error messages */
} c_lisp_map_t;
extern maygc long map_lisp_to_c (object obj, const c_lisp_map_t *map);
extern maygc object map_c_to_lisp (long val, const c_lisp_map_t *map);
extern maygc object map_c_to_list (long val, const c_lisp_map_t *map);
global maygc long map_list_to_c (object obj, const c_lisp_map_t *map);
%% emit_typedef("struct { long c_const; gcv_object_t *l_const; }","c_lisp_pair_t");
%% emit_typedef("struct { const c_lisp_pair_t *table; const unsigned int size; const long default_value; const bool have_default_value_p;  const bool use_default_function_p; const char *name; }","c_lisp_map_t");
%% puts("extern long map_lisp_to_c (object obj, const c_lisp_map_t *map);");
%% puts("extern object map_c_to_lisp (long val, const c_lisp_map_t *map);");
%% puts("extern object map_c_to_list (long val, const c_lisp_map_t *map);");
%% puts("extern long map_list_to_c (object obj, const c_lisp_map_t *map);");

# ####################### ERRBIBL for ERROR.D ############################## #

# Classification of the known condition-types:
# (More precisely, all these are the SIMPLE-... types.)
typedef enum {
  condition, # all kinds of conditions
    serious_condition, # conditions that require interactive intervention
      error, # serious conditions that occur deterministically
        program_error, # mostly statically detectable errors of a program
          source_program_error, # statically detectable errors of a program,
                                # source available
        control_error, # not statically detectable errors in program control
        arithmetic_error, # errors that occur while doing arithmetic operations
          division_by_zero, # eval a mathematical function at a singularity
          floating_point_overflow, # trying to get too close to infinity and...
          floating_point_underflow, # trying to get too close to zero
                                    # in the floating point domain
        cell_error, # trying to access a location which contains #<UNBOUND>
          unbound_variable, # trying to get the value of an unbound variable
          undefined_function, # trying to get the global function definition
                              # of an undefined function
          unbound_slot, # trying to get the value of an unbound slot
        type_error, # some datum does not belong to the expected type
          keyword_error, # a keyword is not one of the allowed keywords
          charset_type_error, # a character does not belong to a character set
          argument_list_dotted, # an argument list in APPLY is dotted
        package_error, # errors during operation on packages
        print_not_readable, # attempted violation of *PRINT-READABLY*
        parse_error, # errors related to parsing
        stream_error, # errors while doing stream I/O
          end_of_file, # unexpected end of stream
          reader_error, # parsing/tokenization error during READ
        file_error, # errors with pathnames, OS level errors with streams
        os_error, # general OS errors
      storage_condition, # "Virtual memory exhausted"
      interrupt_condition, # "User break"
    warning, # conditions for which user notification is appropriate
  # junk
  condition_for_broken_compilers_that_dont_like_trailing_commas
} condition_t;
%% printf("typedef enum { condition=%d, serious_condition=%d, error=%d, program_error=%d, source_program_error=%d, control_error=%d, arithmetic_error=%d, division_by_zero=%d, floating_point_overflow=%d, floating_point_underflow=%d, cell_error=%d, unbound_variable=%d, undefined_function=%d, unbound_slot=%d, type_error=%d, keyword_error=%d, charset_type_error=%d, package_error=%d, print_not_readable=%d, parse_error=%d, stream_error=%d, end_of_file=%d, reader_error=%d, file_error=%d, os_error=%d, storage_condition=%d, interrupt_condition=%d, warning=%d } condition_t;\n",condition, serious_condition, error, program_error, source_program_error, control_error, arithmetic_error, division_by_zero, floating_point_overflow, floating_point_underflow, cell_error, unbound_variable, undefined_function, unbound_slot, type_error, keyword_error, charset_type_error, package_error, print_not_readable, parse_error, stream_error, end_of_file, reader_error, file_error, os_error, storage_condition, interrupt_condition, warning);

/* Error with error-string. Does not return.
 fehler(errortype,errorstring);
 > errortype: condition-type
 > errorstring: constant ASCIZ-String, in UTF-8 Encoding.
   At every tilde-S, a LISP-object is taken from the STACK and printed
   instead of the tilde-S.
 > on the STACK: initial values for the Condition, depending on error-type */
nonreturning_function(extern, fehler, (condition_t errortype, const char * errorstring));
/* used by all modules */
%% puts("nonreturning_function(extern, fehler, (condition_t errortype, const char * errorstring));");

/* Report an error and try to recover by asking the user to supply a value.
 check_value(errortype,errorstring);
 > errortype: condition-type
 > errorstring: constant ASCIZ-String, in UTF-8 Encoding.
   At every tilde-S, a LISP-object is taken from the STACK and printed
   instead of the tilde-S.
 > on the STACK: PLACE (form to be shown to the user) or NIL, then
   the initial values for the Condition, depending on error-type
 < value1, value2: return values from CHECK-VALUE:
   value1 = value supplied by the user,
   value2 = indicates whether PLACE should be filled
 < STACK: cleaned up
 can trigger GC */
extern maygc void check_value (condition_t errortype, const char * errorstring);
/* used by all modules */
%% puts("extern void check_value (condition_t errortype, const char * errorstring);");

/* Report an error and try to recover by asking the user to choose among some
 alternatives.
 correctable_error(errortype,errorstring);
 > errortype: condition-type
 > errorstring: constant ASCIZ-String, in UTF-8 Encoding.
   At every tilde-S, a LISP-object is taken from the STACK and printed
   instead of the tilde-S.
 > on the STACK: list of alternatives
   ((restart-name restart-help-string . value-returned-by-the-restart)*), then
   the initial values for the Condition, depending on error-type
 < value1: return value from CORRECTABLE-ERROR, one of the CDDRs of the
   alternatives
 < STACK: cleaned up
 can trigger GC */
extern maygc void correctable_error (condition_t errortype, const char* errorstring);
/* use by PACKAGE */

# Just like OS_error, but signal a FILE-ERROR.
# OS_file_error(pathname);
# > pathname: Pathname
# > end_system_call() already called
nonreturning_function(extern, OS_file_error, (object pathname));
#if defined(DEBUG_OS_ERROR)
  # Show the file and line number of the caller of OS_file_error().
  # For debugging.
  #define OS_file_error(pathname)  \
    (fprintf(stderr,"\n[%s:%d] ",__FILE__,__LINE__), (OS_file_error)(pathname))
#endif
%% puts("nonreturning_function(extern, OS_file_error, (object pathname));");

# Just like OS_error, but takes a channel stream and signals a FILE-ERROR.
# OS_filestream_error(stream);
# > stream: a channel stream
# > end_system_call() already called
nonreturning_function(extern, OS_filestream_error, (object stream));
#if defined(DEBUG_OS_ERROR)
  # Show the file and line number of the caller of OS_filestream_error(). For debugging.
  #define OS_filestream_error(stream)  \
    (fprintf(stderr,"\n[%s:%d] ",__FILE__,__LINE__), (OS_filestream_error)(stream))
#endif
%% puts("nonreturning_function(extern, OS_filestream_error, (object stream));");

/* Prints error directly via the OS:  errno_out_low(errorcode,FILE,LINE);
 > errorcode: error code
 > FILE: Filename (with quotation marks) as constant ASCIZ-String
 > LINE: line number */
#if defined(UNIX)
extern void errno_out_low (int errorcode, const char* file, uintL line);
#endif
#if defined(WIN32_NATIVE)
extern void errno_out_low (DWORD errorcode, const char* file, uintL line);
#endif
/* Show the file and line number of the caller of errno_out(). */
#define errno_out(e)   errno_out_low(e,__FILE__,__LINE__)

# UP: Executes break-loop because of a keyboard-interrupt.
# > -(STACK) : calling funtion
# modifies STACK, can trigger GC
extern maygc void tast_break (void);
# is used by EVAL, IO, SPVW, STREAM

/* check_classname(obj,type)
 > obj: an object
 > classname: a symbol expected to name a class with "proper name" classname
 < result: an object of the given type, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_classname (object obj, object type);
%% puts("extern object check_classname (object obj, object classname);");

#ifdef FOREIGN
/* check_fpointer(obj,restart_p)
 > obj: an object
 > restart_p: flag whether to allow entering a replacement
 < result: a valid foreign pointer, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_fpointer_replacement (object obj, bool restart_p);
#ifndef COMPILE_STANDALONE
static inline maygc object check_fpointer (object obj, bool restart_p) {
  if (!(fpointerp(obj) && fp_validp(TheFpointer(obj))))
    obj = check_fpointer_replacement(obj,restart_p);
  return obj;
}
#endif
/* used by FOREIGN and REGEXP (amd maybe other non-FFI modules) */
#endif
%% #ifdef FOREIGN
%%   puts("extern object check_fpointer_replacement (object obj, bool restart_p);");
%%   puts("#ifndef COMPILE_STANDALONE");
%%   puts("static inline object check_fpointer (object obj, bool restart_p) {"
%%          " if (!(fpointerp(obj) && fp_validp(TheFpointer(obj))))"
%%            " obj = check_fpointer_replacement(obj,restart_p);"
%%          " return obj;"
%%        " }");
%%   puts("#endif");
%% #endif

# Error message, if an object isn't a list.
# fehler_list(obj);
# > obj: non-list
nonreturning_function(extern, fehler_list, (object obj));
/* used by LIST, EVAL, STREAM */
%% #if notused
%% puts("nonreturning_function(extern, fehler_list, (object obj));");
%% #endif

/* check_list(obj)
 > obj: an object
 < result: a list, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_list_replacement (object obj);

/* create check function for name */
#define MAKE_CHECK_LOW(name,test)                               \
  static inline maygc object check_##name (object obj) {        \
    if (!test)                                                  \
      obj = check_##name##_replacement(obj);                    \
    return obj;                                                 \
  }
#define MAKE_CHECK(name) MAKE_CHECK_LOW(name,name##p(obj))
#define MAKE_CHECK_(name) MAKE_CHECK_LOW(name,name##_p(obj))

#ifndef COMPILE_STANDALONE
MAKE_CHECK(list)
#endif
/* used by PATHNAME */
%% puts("extern object check_list_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK(list));
%% puts("#endif");

# Error message, if an object isn't a proper list because it is dotted.
# fehler_proper_list_dotted(caller,obj);
# > caller: caller (a symbol)
# > obj: End of list, non-list
nonreturning_function(extern, fehler_proper_list_dotted, (object caller, object obj));
# is used by LIST
%% puts("nonreturning_function(extern, fehler_proper_list_dotted, (object caller, object obj));");

# Error message, if an object isn't a proper list because it is circular.
# fehler_proper_list_circular(caller,obj);
# > caller: caller (a symbol)
# > obj: circular list
nonreturning_function(extern, fehler_proper_list_circular, (object caller, object obj));
# is used by LIST

/* check_symbol(obj)
 > obj: an object
 < result: a symbol, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_symbol_replacement (object obj);
#ifndef COMPILE_STANDALONE
static inline maygc object check_symbol (object obj) {
  if (!symbolp(obj))
    obj = check_symbol_replacement(obj);
  return obj;
}
#endif
/* used by CONTROL, EVAL, I18N, PACKAGE, PREDTYPE, RECORD, STREAM, SYMBOL */

/* check_symbol_non_constant(obj,caller)
 > obj: an object
 > caller: a symbol
 < result: a non-constant symbol, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_symbol_non_constant_replacement (object obj, object caller);
#ifndef COMPILE_STANDALONE
static inline maygc object check_symbol_non_constant (object obj, object caller) {
  if (!(symbolp(obj) && !constant_var_p(TheSymbol(obj))))
    obj = check_symbol_non_constant_replacement(obj,caller);
  return obj;
}
#endif
/* used by EVAL, CONTROL */

/* UP: signal an error if a non-symbol was declared special
 returns the symbol
 can trigger GC */
extern maygc object check_symbol_special (object obj, object caller);
/* used by EVAL, CONTROL */

/* Error message, if an object isn't a Simple-Vector.
 fehler_kein_svector(caller,obj);
 > caller: caller (a Symbol)
 > obj: non-Svector */
nonreturning_function(extern, fehler_kein_svector, (object caller, object obj));
/* is used by ARRAY, EVAL */
%% #if notused
%% puts("nonreturning_function(extern, fehler_kein_svector, (object caller, object obj));");
%% #endif

/* Error message, if an object isn't a vector.
 fehler_vector(obj);
 > obj: non-vector */
nonreturning_function(extern, fehler_vector, (object obj));
/* is used by ARRAY */
%% #if notused
%% puts("nonreturning_function(extern, fehler_vector, (object obj));");
%% #endif

/* check_array(obj)
 > obj: an object
 < result: an array, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_array_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK(array)
#endif
/* used by ARRAY, modules */
%% puts("extern object check_array_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK(array));
%% puts("#endif");

/* check_vector(obj)
 > obj: an object
 < result: an vector, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_vector_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK(vector)
#endif
/* used by ARRAY, ENCODING, modules */
%% puts("extern object check_vector_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK(vector));
%% puts("#endif");

/* check_byte_vector_replacement(obj)
 > obj: not an (ARRAY (UNSIGNED-BYTE 8) (*))
 < result: an (ARRAY (UNSIGNED-BYTE 8) (*)), a replacement
 can trigger GC */
extern maygc object check_byte_vector_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_LOW(byte_vector,bit_vector_p(Atype_8Bit,obj))
#endif
/* used by STREAM, modules */
%% puts("extern object check_byte_vector_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_LOW(byte_vector,bit_vector_p(Atype_8Bit,obj)));
%% puts("#endif");


/* error-message, if an object is not an environment.
 fehler_environment(obj);
 > obj: non-vector */
nonreturning_function(extern, fehler_environment, (object obj));
/* used by EVAL, CONTROL */

/* Error message, if an argument isn't a Fixnum >=0:
 > obj: the faulty argument */
nonreturning_function(extern, fehler_posfixnum, (object obj));
/* used by DEBUG, ENCODING, FOREIGN, IO, STREAM, TIME */

/* check_posfixnum(obj)
 > obj: an object
 < result: a fixnum >= 0, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_posfixnum_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK(posfixnum)
#endif
/* used by STREAM, LISPARIT */
%% puts("extern object check_posfixnum_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK(posfixnum));
%% puts("#endif");

/* check_integer(obj)
 > obj: an object
 < result: an integer, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_integer_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK(integer)
#endif
/* used by LISPARIT */

/* check_pos_integer(obj)
 > obj: an object
 < result: an integer >= 0, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_pos_integer_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_LOW(pos_integer,(integerp(obj)&&!R_minusp(obj)))
#endif
/* used by LISPARIT, LIST */
%% puts("extern object check_pos_integer_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_LOW(pos_integer,(integerp(obj)&&!R_minusp(obj))));
%% puts("#endif");

# Error message, if an argument isn't a Character:
# fehler_char(obj);
# > obj: the faulty argument
nonreturning_function(extern, fehler_char, (object obj));
/* used by IO, STREAM */

/* check_char(obj)
 > obj: an object
 < result: a character, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_char_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK(char)
#endif
/* used by CHARSTRG, ENCODING, IO */
%% #if notused
%% puts("extern object check_char_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK(char));
%% puts("#endif");
%% #endif

/* check_string(obj)
 > obj: an object
 < result: a string, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_string_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK(string)
#endif
/* used by CHARSTRG, FOREIGN, MISC, PACKAGE, PATHNAME, STREAM, SOCKET, I18N */
%% puts("extern object check_string_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK(string));
%% puts("#endif");

# Error message, if an argument isn't a Simple-String:
# fehler_sstring(obj);
# > obj: the faulty argument
nonreturning_function(extern, fehler_sstring, (object obj));
# is used by CHARSTRG
%% #if notused
%% puts("nonreturning_function(extern, fehler_sstring, (object obj));");
%% #endif

# Checks a simple-string for being mutable.
# check_sstring_mutable(string);
  #define check_sstring_mutable(obj)  \
    if (sstring_immutable(TheSstring(obj))) \
      fehler_sstring_immutable(obj);
  # Error message, if a Simple-String is immutable:
  # fehler_sstring_immutable(obj);
  # > obj: der String
  nonreturning_function(extern, fehler_sstring_immutable, (object obj));
  # is used by Macro check_sstring_mutable

# Error message, if an argument is not of type (OR STRING INTEGER).
# fehler_string_integer(obj);
nonreturning_function(extern, fehler_string_integer, (object obj));
%% puts("nonreturning_function(extern, fehler_string_integer, (object obj));");

# Error message, if a string size is too big.
# fehler_stringsize(size);
# > size: the desired string length
nonreturning_function(extern, fehler_stringsize, (uintV size));

# Check a string size, reporting an error when it's too big.
#define check_stringsize(size)  \
  if ((size) > stringsize_limit_1) \
    fehler_stringsize(size)/*;*/

/* error message if an argument is not a class.
 fehler_class(caller,obj);
 > obj: the erroneous argument */
nonreturning_function(extern, fehler_class, (object obj));

/* Report an error when the argument is not an encoding:
 check_encoding(obj,&default,keyword_p)
 > obj: the (possibly) bad argument
 > default: what to return for :DEFAULT
 > keyword_p: true if the object comes from the :EXTERNAL-FORMAT argument
 < result: an encoding
 can trigger GC */
extern maygc object check_encoding (object obj, const gcv_object_t* e_default,
                                    bool keyword_p);
/* used by ENCODING, FOREIGN */

/* Error when the property list has odd length
 fehler_plist_odd(caller,plist);
 > plist: bad plist */
nonreturning_function(extern, fehler_plist_odd, (object plist));
%% puts("nonreturning_function(extern, fehler_plist_odd, (object plist));");

/* error-message for non-paired keyword-arguments
 fehler_key_odd(argcount,caller);
 > argcount: the number of arguments on the STACK
 > caller: function */
nonreturning_function(extern, fehler_key_odd, (uintC argcount, object caller));
%% puts("nonreturning_function(extern, fehler_key_odd, (uintC argcount, object caller));");

/* error-message for flawed keyword
 fehler_key_notkw(kw);
 > key: Non-Symbol
 > caller: function */
nonreturning_function(extern, fehler_key_notkw, (object key, object caller));

/* error-message for flawed keyword
 fehler_key_badkw(fun,kw,kwlist);
 > fun: function
 > key: illegal keyword
 > val: its value
 > kwlist: list of legal keywords */
nonreturning_function(extern, fehler_key_badkw,
                      (object fun, object key, object val, object kwlist));
%% puts("nonreturning_function(extern, fehler_key_badkw, (object fun, object key, object val, object kwlist));");

/* check_function(obj)
 > obj: an object
 < result: a function, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_function_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK(function)
#endif
/* used by RECORD, EVAL, SEQUENCE, SYMBOL, FOREIGN */

/* error if funname does not have a function definition
 check_fdefinition(funname,caller)
 > funname: symbol or (setf symbol)
 > caller: symbol
 < a function object, possibly also installed as (FDEFINITION funname)
 can trigger GC */
extern maygc object check_fdefinition (object funname, object caller);
/* used by EVAL, CONTROL */

/* check_funname(obj)
 > errtype: type of condition to signal if obj is not a function name,
            either type_error or source_program_error
 > caller: a symbol
 > obj: an object
 < result: a function name, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_funname_replacement (condition_t errtype, object caller, object obj);
#ifndef COMPILE_STANDALONE
static inline maygc object check_funname (condition_t errtype, object caller, object obj) {
  if (!funnamep(obj))
    obj = check_funname_replacement(errtype,caller,obj);
  return obj;
}
#endif
/* used by EVAL, CONTROL */

# Error message, if an argument is a lambda-expression instead of a function:
# fehler_lambda_expression(caller,obj);
# caller: caller (a symbol)
# obj: the faulty argument
nonreturning_function(extern, fehler_lambda_expression, (object caller, object obj));
# is used by EVAL, SYMBOL

# too many/few arguments in a function call
# > caller : the function that is reporting the error (unbound == EVAL/APPLY)
# > func   : the function being incorrectly called
# > ngiven : the number of arguments given
# < nmax   : the maximum number of arguments accepted
# < nmin   : the minimum number of arguments required
nonreturning_function(extern, fehler_too_many_args,
                      (object caller, object func, uintL ngiven, uintL nmax));
nonreturning_function(extern, fehler_too_few_args,
                      (object caller, object func, uintL ngiven, uintL nmin));

# used by EVAL, FOREIGN

# Error message, if an argument isn't of a given elementary C type.
# fehler_<ctype>(obj);
# > obj: the faulty argument
nonreturning_function(extern, fehler_uint8, (object obj));
nonreturning_function(extern, fehler_sint8, (object obj));
nonreturning_function(extern, fehler_uint16, (object obj));
nonreturning_function(extern, fehler_sint16, (object obj));
nonreturning_function(extern, fehler_uint32, (object obj));
nonreturning_function(extern, fehler_sint32, (object obj));
nonreturning_function(extern, fehler_uint64, (object obj));
nonreturning_function(extern, fehler_sint64, (object obj));
# nonreturning_function(extern, fehler_uint, (object obj));
# nonreturning_function(extern, fehler_sint, (object obj));
#if (int_bitsize==16)
  #define fehler_uint  fehler_uint16
  #define fehler_sint  fehler_sint16
#else # (int_bitsize==32)
  #define fehler_uint  fehler_uint32
  #define fehler_sint  fehler_sint32
#endif
# nonreturning_function(extern, fehler_ulong, (object obj));
# nonreturning_function(extern, fehler_slong, (object obj));
#if (long_bitsize==32)
  #define fehler_ulong  fehler_uint32
  #define fehler_slong  fehler_sint32
#else # (long_bitsize==64)
  #define fehler_ulong  fehler_uint64
  #define fehler_slong  fehler_sint64
#endif
/* used by STREAM, ENCODING, modules */
%% puts("nonreturning_function(extern, fehler_uint8, (object obj));");
%% puts("nonreturning_function(extern, fehler_sint8, (object obj));");
%% puts("nonreturning_function(extern, fehler_uint16, (object obj));");
%% puts("nonreturning_function(extern, fehler_sint16, (object obj));");
%% puts("nonreturning_function(extern, fehler_uint32, (object obj));");
%% puts("nonreturning_function(extern, fehler_sint32, (object obj));");
%% puts("nonreturning_function(extern, fehler_uint64, (object obj));");
%% puts("nonreturning_function(extern, fehler_sint64, (object obj));");
%% #if (int_bitsize==16)
%%   emit_define("fehler_uint","fehler_uint16");
%%   emit_define("fehler_sint","fehler_sint16");
%% #else
%%   emit_define("fehler_uint","fehler_uint32");
%%   emit_define("fehler_sint","fehler_sint32");
%% #endif
%% #if (long_bitsize==32)
%%   emit_define("fehler_ulong","fehler_uint32");
%%   emit_define("fehler_slong","fehler_sint32");
%% #else
%%   emit_define("fehler_ulong","fehler_uint64");
%%   emit_define("fehler_slong","fehler_sint64");
%% #endif

/* Check whether an object can be converted to an elementary C type.
 check_<ctype>(obj)
 > obj: an object
 < result: an object that can be converted to the C type, either the same
           as obj or a replacement
 can trigger GC */
extern maygc object check_uint8_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(uint8)
#endif
extern maygc object check_sint8_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(sint8)
#endif
extern maygc object check_uint16_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(uint16)
#endif
extern maygc object check_sint16_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(sint16)
#endif
extern maygc object check_uint32_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(uint32)
#endif
extern maygc object check_sint32_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(sint32)
#endif
extern maygc object check_uint64_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(uint64)
#endif
extern maygc object check_sint64_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(sint64)
#endif
extern maygc object check_uint_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(uint)
#endif
extern maygc object check_sint_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(sint)
#endif
extern maygc object check_ulong_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(ulong)
#endif
extern maygc object check_slong_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_(slong)
#endif
extern maygc object check_ffloat_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_LOW(ffloat,single_float_p(obj))
#endif
extern maygc object check_dfloat_replacement (object obj);
#ifndef COMPILE_STANDALONE
MAKE_CHECK_LOW(dfloat,double_float_p(obj))
#endif
# is used by STREAM, FFI
%% puts("extern object check_uint8_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(uint8));
%% puts("#endif");
%% puts("extern object check_sint8_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(sint8));
%% puts("#endif");
%% puts("extern object check_uint16_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(uint16));
%% puts("#endif");
%% puts("extern object check_sint16_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(sint16));
%% puts("#endif");
%% puts("extern object check_uint32_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(uint32));
%% puts("#endif");
%% puts("extern object check_sint32_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(sint32));
%% puts("#endif");
%% puts("extern object check_uint64_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(uint64));
%% puts("#endif");
%% puts("extern object check_sint64_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(sint64));
%% puts("#endif");
%% puts("extern object check_uint_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(uint));
%% puts("#endif");
%% puts("extern object check_sint_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(sint));
%% puts("#endif");
%% puts("extern object check_ulong_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(ulong));
%% puts("#endif");
%% puts("extern object check_slong_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_(slong));
%% puts("#endif");
%% puts("extern object check_ffloat_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_LOW(ffloat,single_float_p(obj)));
%% puts("#endif");
%% puts("extern object check_dfloat_replacement (object obj);");
%% puts("#ifndef COMPILE_STANDALONE");
%% export_literal(MAKE_CHECK_LOW(dfloat,double_float_p(obj)));
%% puts("#endif");

# ##################### PACKBIBL for PACKAGE.D ############################# #

# UP: tests whether a symbol is accessible in a package and isn't hidden
# by a another symbol with the same name.
# accessiblep(sym,pack)
# > sym: Symbol
# > pack: Package
# < result: true if sym is accessible in pack and nod hidden,
#             else false
extern bool accessiblep (object sym, object pack);
# is used by IO

# UP: tests whether a symbol is accessible as external symbol in a package
# externalp(sym,pack)
# > sym: Symbol
# > pack: Package
# < result: true if sym is accessible as external symbol in pack ,
#           else false
extern bool externalp (object sym, object pack);
# is used by IO

/* UP: locates an external symbol with a given printname in a package.
 find_external_symbol(string,invert,pack,&sym)
 > string: string
 > invert: whether to implicitly case-invert the string
 > pack: package
 < result: true, if an external symbol with that printname has been found in pack.
 < sym: this symbol, if found. */
extern bool find_external_symbol (object string, bool invert, object pack, object* sym_);
# is used by IO

# UP: locates a package with a given name or nickname
# find_package(string)
# > string: String
# < result: Package with that name or NIL
extern object find_package (object string);
# is used by IO, EVAL
%% puts("extern object find_package (object string);");

# UP: Interns a symbol with a given printname in a package.
# intern(string,invert,pack,&sym)
# > string: String
# > invert: whether to implicitly case-invert the string
# > pack: Package
# < sym: Symbol
# < result: 0, if not found but newly created
#           1, if found as external symbol
#           2, if inherited through use-list
#           3, if exists as internal symbol
# can trigger GC
extern maygc uintBWL intern (object string, bool invert, object pack, object* sym_);
# is used by IO, SPVW
%% puts("extern uintBWL intern (object string, object pack, object* sym_);");

# UP: Interns a symbol with a given printname in the Keyword-Package.
# intern_keyword(string)
# > string: String
# < result: Symbol, a keyword
# can trigger GC
extern maygc object intern_keyword (object string);
# is used by IO, EVAL, GRAPH
%% puts("extern object intern_keyword (object string);");

# UP: Imports a symbol into a package
# import(&sym,&pack);
# > sym: Symbol (on STACK)
# > pack: Package (on STACK)
# < sym: Symbol, EQ with the old
# < pack: Package, EQ with the old
# can trigger GC
extern maygc void import (const gcv_object_t* sym_, const gcv_object_t* pack_);
# is used by SPVW

# UP: Exports a symbol from a package
# export(&sym,&pack);
# > sym: Symbol (on STACK)
# > pack: Package (on STACK)
# < sym: Symbol, EQ with the old
# < pack: Package, EQ with the old
# can trigger GC
extern maygc void export (const gcv_object_t* sym_, const gcv_object_t* pack_);
# is used by SPVW

/* UP: gets the current package
 get_current_package()
 < result: current Package
 can trigger GC */
extern maygc object get_current_package (void);
/* is used by IO */

/* check whether package lock prevents assignment to symbol
 can trigger GC */
extern maygc void symbol_value_check_lock (object caller, object symbol);
/* used by EVAL */

# UP: Initializes the package-management
# init_packages();
# can trigger GC
extern void init_packages (void);
# is used by SPVW

# ##################### PATHBIBL for PATHNAME.D ############################ #

/* Check that the namestring for path will be parsed into a similar object
 used by pr_orecord() in io.d
 can trigger GC */
extern maygc bool namestring_correctly_parseable_p (gcv_object_t *path_);

/* Converts an object into an absolute physical pathname and returns its
   namestring.
 physical_namestring(thing)
 > thing: an object
 < result: the namestring of the pathname denoted by thing
 can trigger GC */
extern maygc object physical_namestring (object thing);
%% puts("extern object physical_namestring (object obj);");

# UP: Gives the directory-namestring in OS-format of a halfway checked
#     pathname assuming that the directory of the pathname exists.
# assume_dir_exists()
# > STACK_0: absolute pathname, halfway checked
# < STACK_0: (poss. the same) pathname, better resolved
# < result:
#     if Name=NIL: Directory-Namestring (for the OS)
#     if Name/=NIL: Namestring (for the OS, with nullbyte at the end)
# can trigger GC
extern maygc object assume_dir_exists (void);
# is used by STREAM

/* Converts a directory pathname to an OS directory specification.
 > pathname: an object
 > use_default: whether to use the current default directory
 < result: a simple-bit-vector containing an ASCIZ string in OS format
 can trigger GC */
extern maygc object pathname_to_OSdir (object pathname, bool use_default);
/* used by modules (I18N) */
%% puts("extern object pathname_to_OSdir (object pathname, bool use_default);");

/* Converts an OS directory specification to a directory pathname.
 > path: a pathname referring to a directory
 < result: a pathname without name and type
 can trigger GC */
extern maygc object OSdir_to_pathname (const char* path);
/* used by modules (I18N) */
%% puts("extern object OSdir_to_pathname (const char* path);");

# UP: Initializes the pathname-system.
# init_pathnames();
# can trigger GC
extern maygc void init_pathnames (void);
# is used by SPVW

/* Duplicate an open file handle.
 handle_dup(oldfd)
 Similar to dup(oldfd), with error checking.
 To be called only inside begin/end_system_call(). */
extern Handle handle_dup (Handle old_handle);
%% puts("extern Handle handle_dup (Handle old_handle);");

/* Duplicate an open file handle.
 handle_dup2(oldfd,newfd)
 Similar to dup2(oldfd,newfd), with error checking. The result may or may not
 be equal to newfd.
 To be called only inside begin/end_system_call(). */
extern Handle handle_dup2 (Handle old_handle, Handle new_handle);
%% puts("extern Handle handle_dup2 (Handle old_handle, Handle new_handle);");

# Locates the executable program immediately after the program start.
# find_executable(argv[0])
extern int find_executable (const char * program_name);
# is used by SPVW

# check the :DIRECTION argument
# return one of the following:
typedef enum {
  # see READ_P, RO_P and WRITE_P in <stream.d> regarding the choice of values
  DIRECTION_PROBE           = 0,
  DIRECTION_INPUT           = 1,
  DIRECTION_INPUT_IMMUTABLE = 3,
  DIRECTION_OUTPUT          = 4,
  DIRECTION_IO              = 5,
  /* Work around a g++-3.4.0 bug, see
     http://gcc.gnu.org/bugzilla/show_bug.cgi?id=15069 */
  DIRECTION_DUMMY_TO_AVOID_GXX_BUG = 100
} direction_t;
extern direction_t check_direction (object dir);
%% printf("typedef enum { DIRECTION_PROBE=%d, DIRECTION_INPUT=%d, DIRECTION_INPUT_IMMUTABLE=%d, DIRECTION_OUTPUT=%d, DIRECTION_IO=%d} direction_t;\n",
%%        DIRECTION_PROBE, DIRECTION_INPUT, DIRECTION_INPUT_IMMUTABLE,
%%        DIRECTION_OUTPUT, DIRECTION_IO);
%% puts("extern direction_t check_direction (object dir);");

# check the :IF-DOES-NOT-EXIST argument
# check_if_does_not_exist(argument)
# return one of the following:
typedef enum {
  IF_DOES_NOT_EXIST_UNBOUND,
  IF_DOES_NOT_EXIST_ERROR,
  IF_DOES_NOT_EXIST_NIL,
  IF_DOES_NOT_EXIST_CREATE
} if_does_not_exist_t;
extern if_does_not_exist_t check_if_does_not_exist (object if_not_exist);
%% emit_typedef("enum { IF_DOES_NOT_EXIST_UNBOUND, IF_DOES_NOT_EXIST_ERROR, IF_DOES_NOT_EXIST_NIL, IF_DOES_NOT_EXIST_CREATE }","if_does_not_exist_t");
%% puts("extern if_does_not_exist_t check_if_does_not_exist (object if_not_exist);");

# Converts a :IF-DOES-NOT-EXIST enum item to a symbol.
# if_does_not_exist_symbol(item)
extern object if_does_not_exist_symbol (if_does_not_exist_t if_not_exist);
%% puts("extern object if_does_not_exist_symbol (if_does_not_exist_t if_not_exist);");

# check the :IF-EXISTS argument
# check_if_exists(argument)
# return one of the following:
typedef enum {
  IF_EXISTS_UNBOUND,
  IF_EXISTS_ERROR,
  IF_EXISTS_NIL,
  IF_EXISTS_RENAME,
  IF_EXISTS_RENAME_AND_DELETE,
  IF_EXISTS_SUPERSEDE,
  IF_EXISTS_APPEND,
  IF_EXISTS_OVERWRITE
} if_exists_t;
extern if_exists_t check_if_exists (object if_exists);
%% emit_typedef("enum { IF_EXISTS_UNBOUND, IF_EXISTS_ERROR, IF_EXISTS_NIL, IF_EXISTS_RENAME, IF_EXISTS_RENAME_AND_DELETE, IF_EXISTS_SUPERSEDE, IF_EXISTS_APPEND, IF_EXISTS_OVERWRITE }","if_exists_t");
%% puts("extern if_exists_t check_if_exists (object if_exists);");

# Converts a :IF-EXISTS enum item to a symbol.
# if_exists_symbol(item)
extern object if_exists_symbol (if_exists_t if_exists);
%% puts("extern object if_exists_symbol (if_exists_t if_exists);");

#if defined(WIN32_NATIVE)
/* ------------------- Functions defined in w32shell.c ------------------- */

/* shell_quote() - surround dangerous strings with double quotes.
 escape quotes and backslashes.
 dest should be twice as large as source
  + 2 (for quotes) + 1 for zero byte + 1 for possible endslash */
extern int shell_quote (char * dest, const char * source);
/* used by PATHNAME and the driver clisp.exe */

/* real_path() - the ultimate shortcut megaresolver
   style inspired by directory_search_scandir
 > namein: filename pointing to file or directory
            wildcards (only asterisk) may appear only as filename
 < nameout: filename with directory and file shortcuts resolved
             on failure holds filename resolved so far
 < result:  true if resolving succeeded */
extern BOOL real_path (LPCSTR namein, LPSTR nameout);
/* used by PATHNAME, SPVW [for loadmem()] and the driver clisp.exe */

#endif

# ##################### PREDBIBL for PREDTYPE.D ############################ #

# UP: test for atomic equality EQL
# eql(obj1,obj2)
# > obj1,obj2: Lisp-objects
# < result: true, if objects are equal
extern bool eql (object obj1, object obj2);
# is used by CONTROL, EVAL, HASHTABL, LISPARIT
%% #if notused
%% puts("extern bool eql (object obj1, object obj2);");
%% #endif

# UP: tests for equality EQUAL
# equal(obj1,obj2)
# > obj1,obj2: Lisp-objects
# < result: true, if objects are equal
extern bool equal (object obj1, object obj2);
# is used by EVAL, PATHNAME, HASHTABL, MISC
%% #if notused
%% puts("extern bool equal (object obj1, object obj2);");
%% #endif

# UP: tests for a more lax equality EQUALP
# equalp(obj1,obj2)
# > obj1,obj2: Lisp-objects
# < result: true, if objects are equal
extern bool equalp (object obj1, object obj2);
# is used by PATHNAME, HASHTABL
%% #if notused
%% puts("extern bool equalp (object obj1, object obj2);");
%% #endif

/* typep_class(obj,clas)
 > obj: an object
 > clas: a class object
 < true if the object is an instance of the class, false otherwise
 clobbers value1, mv_count */
extern bool typep_class (object obj, object clas);
%% puts("extern bool typep_class (object obj, object clazz);");

/* typep_classname(obj,classname)
 > obj: an object
 > classname: a symbol expected to name a class with "proper name" classname
 < true if the object is an instance of the class, false otherwise
 clobbers value1, mv_count */
extern bool typep_classname (object obj, object classname);
%% puts("extern bool typep_classname (object obj, object classname);");

# UP: expand all DEFTYPE definitions in the type spec
# (recursively, unless once_p is true)
# > type_spec: Lisp object
# < result: the expansion (when not a deftyped type, returns the argument)
# can trigger GC
extern maygc object expand_deftype (object type_spec, bool once_p);
# used by predtype.d, sequence.d

# UP: Makes a statistic about the action of a GC.
# with_gc_statistics(fun);
# > fun: Function that does a GC
typedef void gc_function_t (void);
extern void with_gc_statistics (gc_function_t* fun);
# is used by SPVW

# ####################### RECBIBL for RECORD.D ############################# #

/* check_structure(obj)
 > obj: an object
 < result: a structure object, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_structure_replacement (object obj);
#ifndef COMPILE_STANDALONE
static inline maygc object check_structure (object obj) {
  if (!structurep(obj))
    obj = check_structure_replacement(obj);
  return obj;
}
#endif
/* used by IO */

# instance_un_realloc(obj);
# walks over forward pointers left by instance reallocation (CHANGE-CLASS
# and/or redefined classes).
# > obj: a CLOS instance, possibly a forward pointer
# < obj: the same CLOS instance, not a forward pointer
# Note that the forwarded instance must not be leaked to "userland", because
# the forward pointer and the forwarded instance are not EQ.
#define instance_un_realloc(obj) \
  if (record_flags(TheInstance(obj)) & instflags_forwarded_B) {        \
    (obj) = TheInstance(obj)->inst_class_version;                      \
    /* We know that there is at most one indirection. */               \
    ASSERT(!(record_flags(TheInstance(obj)) & instflags_forwarded_B)); \
  }

# update_instance(user_obj,obj)
# updates a CLOS instance after its class or one of its superclasses has been
# redefined.
# > user_obj: a CLOS instance, possibly a forward pointer
# > obj: the same CLOS instance, not a forward pointer
# < result: the same CLOS instance, not a forward pointer
# can trigger GC
extern maygc object update_instance (object user_obj, object obj);

# instance_valid_p(obj)
# Tests whether a CLOS instance can be used without first updating it.
# > obj: a CLOS instance, not a forward pointer
# < result: false if its class was redefined since the instance was last used
#define instance_valid_p(obj) \
  nullp(TheClassVersion(TheInstance(obj)->inst_class_version)->cv_next)

# instance_update(user_obj,obj);
# performs necessary CLOS instance updates on obj.
# > user_obj: a CLOS instance, possibly a forward pointer
# > obj: the same CLOS instance, not a forward pointer
# < obj: the same CLOS instance, not a forward pointer
# can trigger GC
#define instance_update(user_obj,obj) \
  if (!instance_valid_p(obj)) \
    (obj) = update_instance(user_obj,obj); \
  GCTRIGGER1(obj)/*;*/

/* Test for CLOS instance of a given class
 > obj: a Lisp object
 > clas: a class that doesn't have obsolete instances */
#ifndef COMPILE_STANDALONE
static inline bool instanceof (object obj, object clas) {
  if (!instancep(obj)) return false;
  var object obj_forwarded = obj;
  instance_un_realloc(obj_forwarded);
  /*instance_update(obj,obj_forwarded); - not needed since we don't access a slot */
  var object cv = TheInstance(obj_forwarded)->inst_class_version;
  var object objclas = TheClassVersion(cv)->cv_newest_class;
  return !eq(gethash(clas,TheClass(objclas)->all_superclasses,false),nullobj);
}
#endif

# ###################### SEQBIBL for SEQUENCE.D ############################ #

# UP: Converts an object into a sequence of a given type.
# coerce_sequence(obj, result_type, error_p)
# > obj: Object, should be a sequence
# > result_type: identifier (symbol) of the sequence-type
# > error_p: when result_type does not name a sequence:
#              when true, signal an error; when false, return nullobj
# < value: Sequence of type result_type
# can trigger GC
extern maygc Values coerce_sequence (object sequence, object result_type,
                                     bool error_p);
# is used by PREDTYPE, EVAL

# UP:  Traverses a sequence and calls a function for every element.
# map_sequence(obj,fun,arg);
# > obj: Object, should be a sequence
# > fun: Function, fun(arg,element) may trigger GC
# > arg: arbitrary given argument
# can trigger GC
typedef maygc void map_sequence_function_t (void* arg, object element);
extern maygc void map_sequence (object obj, map_sequence_function_t* fun, void* arg);
/* used by ARRAY, modules */
/* do not use emit_typedef_f because g++ complains:
 error: invalid application of `sizeof' to a function type */
%% puts("typedef void map_sequence_function_t (void* arg, object element);");
%% puts("extern void map_sequence (object obj, map_sequence_function_t* fun, void* arg);");

# Error, if both :TEST, :TEST-NOT - argumente have been given.
# fehler_both_tests();
nonreturning_function(extern, fehler_both_tests, (void));
# is used by LIST, SEQUENCE, WEAK

# ###################### STRMBIBL for STREAM.D ############################# #

/* Error message, if an argument isn't a stream:
 check_stream_replacement(obj);
 > obj: not a stream
 < obj: a stream
 can trigger GC */
extern maygc object check_stream_replacement (object obj);
#ifndef COMPILE_STANDALONE
static inline maygc object check_stream (object obj) {
  if (!streamp(obj))
    obj = check_stream_replacement(obj);
  return obj;
}
#endif
/* is used by IO, STREAM, DEBUG */

/* parse timeout argument
 sec = posfixnum or (SEC . USEC) or (SEC USEC) or float or ratio or nil/unbound
 usec = posfixnum or nil/unbound
 can trigger GC */
extern maygc struct timeval * sec_usec (object sec, object usec, struct timeval *tv);
%% puts("extern struct timeval * sec_usec (object sec, object usec, struct timeval *tv);");

/* Convert C sec/usec (struct timeval et al) pair into Lisp number (of seconds)
 if abs_p is true, add UNIX_LISP_TIME_DIFF
 can trigger GC */
#if defined(SIZEOF_STRUCT_TIMEVAL) && SIZEOF_STRUCT_TIMEVAL == 16
global maygc object sec_usec_number (uint64 sec, uint64 usec, bool abs_p);
#else
global maygc object sec_usec_number (uint32 sec, uint32 usec, bool abs_p);
#endif
%% #if defined(SIZEOF_STRUCT_TIMEVAL) && SIZEOF_STRUCT_TIMEVAL == 16
%% puts("extern object sec_usec_number (uint64 sec, uint64 usec, bool abs_p);");
%% #else
%% puts("extern object sec_usec_number (uint32 sec, uint32 usec, bool abs_p);");
%% #endif

/* UP: Initializes the stream variables.
 init_streamvars(batch_p);
 > batch_p: Flag, whether *standard-input*, *standard-output*, *error-output*
            should be initialized to the C stdio handle-streams
            (deviates from the standard)
 can trigger GC */
extern maygc void init_streamvars (bool batch_p);
/* used by SPVW */

# Error-message, if a stream-operation is not permitted on a stream.
# fehler_illegal_streamop(caller,stream);
# > caller: Caller (a symbol)
# > stream: Stream
nonreturning_function(extern, fehler_illegal_streamop, (object caller, object stream));
# is used by IO

# Reads a byte from a stream.
# read_byte(stream)
# > stream: Stream
# < result: read Integer (eof_value at EOF)
# can trigger GC
extern maygc object read_byte (object stream);
# is used by SEQUENCE

# Writes a byte onto a stream.
# write_byte(stream,byte);
# > stream: Stream
# > byte: Integer to be written
# can trigger GC
extern maygc void write_byte(object stream, object byte);
# is used by SEQUENCE

# Reads a character from a stream.
# read_char(&stream)
# > stream: Stream
# < stream: Stream
# < result: read character (eof_value at EOF)
# can trigger GC
extern maygc object read_char (const gcv_object_t* stream_);
# is used by IO, DEBUG, SEQUENCE

# Pushes the last read character back onto a stream.
# unread_char(&stream,ch);
# > ch: last read character
# > stream: Stream
# < stream: Stream
# can trigger GC
extern maygc void unread_char (const gcv_object_t* stream_, object ch);
# is used by IO, DEBUG

# Reads a character from a stream without using it.
# peek_char(&stream)
# > stream: Stream
# < stream: Stream
# < result: read character (eof_value at EOF)
# can trigger GC
extern maygc object peek_char (const gcv_object_t* stream_);
# is used by IO

# Reads a line of characters from a stream.
# read_line(&stream,&buffer)
# > stream: stream
# > buffer: a semi-simple string
# < stream: stream
# < buffer: contains the read characters, excluding the terminating #\Newline
# < result: true is EOF was seen before newline, else false
# can trigger GC
extern maygc bool read_line (const gcv_object_t* stream_, const gcv_object_t* buffer_);
# used by IO

# Write a character onto a stream.
# write_char(&stream,ch);
# > ch: Character to be written
# > stream: Stream
# < stream: Stream
# can trigger GC
extern maygc void write_char (const gcv_object_t* stream_, object ch);
# is used by LISPARIT, IO, ERROR, SEQUENCE

# Writes a character onto a stream.
# write_code_char(&stream,ch);
# > ch: a character
# > stream: Stream
# < stream: Stream
# can trigger GC
# extern maygc void write_code_char (const gcv_object_t* stream_, chart ch);
#define write_code_char(stream_,ch)  write_char(stream_,code_char(ch))
# is used by LISPARIT, IO

# Writes a fixed standard-char onto a stream.
# write_ascii_char(&stream,ch);
# > ch: a standard char, in ASCII encoding
# > stream: Stream
# < stream: Stream
# can trigger GC
# extern maygc void write_ascii_char (const gcv_object_t* stream_, uintB ch);
#define write_ascii_char(stream_,ch)  write_char(stream_,code_char(as_chart(ch)))
# is used by LISPARIT, IO, DEBUG, Macro TERPRI

#ifdef UNICODE
# Changes a terminal stream's external format.
# > stream: a stream
# > encoding: an encoding
# can trigger GC
extern maygc void set_terminalstream_external_format (object stream, object encoding);
# used by ENCODING
#endif

# UP: Determines whether a stream is "interactive",
#     ie. whether the input from the stream
#     depends from a promt that has propably just been printed.
# interactive_stream_p(stream)
# > stream: Stream
extern bool interactive_stream_p (object stream);
# is used by DEBUG

/* UP: Closes a stream.
 builtin_stream_close(&stream,abort);
 > stream: Builtin-Stream
 > abort: flag: non-0 => ignore errors
 < stream: Builtin-Stream
 can trigger GC */
extern maygc void builtin_stream_close (const gcv_object_t* stream_, uintB abort);
/* is used by PATHNAME, SPVW, DEBUG, MISC */
%% puts("extern void builtin_stream_close (const gcv_object_t* stream_, uintB abort);");

# UP: Closes a list of open files.
# close_some_files(list);
# > list: List of open builtin-streams
# can trigger GC
extern maygc void close_some_files (object list);
# is used by SPVW

# UP: Closes all open files.
# close_all_files();
# can trigger GC
extern maygc void close_all_files (void);
# is used by SPVW

# UP: declares all open file-streams closed.
# closed_all_files();
extern void closed_all_files (void);
# is used by SPVW

# UP: determines whether a char is available in the Stream stream
# listen_char(stream)
# > stream: Stream
# < result: ls_avail if a character is available,
#             ls_eof   if EOF is reached,
#             ls_wait  if no character is available, but not because of EOF
# can trigger GC
extern maygc signean listen_char (object stream);
#define ls_avail  0
#define ls_eof   -1
#define ls_wait   1
#define ls_avail_p(x)  ((x) == 0)
#define ls_eof_p(x)  ((x) < 0)
#define ls_wait_p(x)  ((x) > 0)
# is used by IO, DEBUG

# UP: clears an already entered interactive input from a Stream stream.
# clear_input(stream)
# > stream: Stream
# < result: true if input has been deleted
# can trigger GC
extern maygc bool clear_input (object stream);
# is used by IO, DEBUG

# UP: Determines whether a stream has a byte immediately available.
# listen_byte(stream)
# > stream: a stream with element-type ([UN]SIGNED-BYTE 8)
# < result: ls_avail if a byte is available,
#           ls_eof   if EOF is reached,
#           ls_wait  if no byte is available, but not because of EOF
# can trigger GC
extern maygc signean listen_byte (object stream);
# is used by

# UP: Finishes waiting output of a Stream stream
# finish_output(stream);
# > stream: Stream
# can trigger GC
extern maygc void finish_output (object stream);
# is used by IO

# UP: Forces waiting output of a Stream stream
# force_output(stream);
# > stream: Stream
# can trigger GC
extern maygc void force_output (object stream);
# is used by IO, DEBUG

# UP: clear the waiting output of a stream.
# clear_output(stream);
# > stream: Stream
# can trigger GC
extern maygc void clear_output (object stream);
# is used by IO

# UP: Gives the line position of a stream:
# get_line_position(stream)
# > stream: Stream
# < result: Line-Position (Fixnum >=0 or NIL)
# can trigger GC
extern maygc object get_line_position (object stream);
# is used by IO, DEBUG

# Writes a newline on a stream, if it is not already positioned at column 0.
# fresh_line(&stream);
# > stream: Stream
# < stream: Stream
# < result: true if did output a newline
# can trigger GC
extern maygc bool fresh_line (const gcv_object_t* stream_);
# is used by IO

# Writes a newline on a stream, delayed and nullified if the next character
# written would be a newline anyway.
# elastic_newline(&stream);
# > stream: Stream
# < stream: Stream
# can trigger GC
extern maygc void elastic_newline (const gcv_object_t* stream_);
# is used by IO

/* UP: give away corresponding underlying handle
 making sure buffers were flushed. One can then use the
 handle outside of stream object as far as the latter
 is not used and not GCed.
 stream_lend_handle(stream, inputp, handletype)
 > stream_: stream for handle to extract
 > inputp: whether its input or output side is requested.
 < stream_: corrected stream (if the original argument was not a handle stream)
 < int * handletype 0:reserved, 1:file, 2:socket
 < Handle result - extracted handle
 can trigger GC */
extern maygc Handle stream_lend_handle (gcv_object_t *stream_, bool inputp, int * handletype);
/* used by STREAM */
%% puts("extern Handle stream_lend_handle (gcv_object_t *stream_, bool inputp, int * handletype);");

/* extract the OS file handle from the file stream
 > stream: open Lisp file stream
 < fd: OS file handle
 < result: either stream, or a corrected stream in case stream was invalid
 for syscall module
 can trigger GC */
extern maygc object open_file_stream_handle (object stream, Handle *fd);
%% puts("extern object open_file_stream_handle (object stream, Handle *fd);");

/* Function: Reads several bytes from a stream.
 read_byte_array(&stream,&bytearray,start,len,persev)
 > stream: stream (on the STACK)
 > object bytearray: simple-8bit-vector (on the STACK)
 > uintL start: start index of byte sequence to be filled
 > uintL len: length of byte sequence to be filled
 > perseverance_t persev: how to react on incomplete I/O
 < uintL result: number of bytes that have been filled
 can trigger GC */
extern maygc uintL read_byte_array (const gcv_object_t* stream_, const gcv_object_t* bytearray_, uintL start, uintL len, perseverance_t persev);
/* used by SEQUENCE, PATHNAME */
%% puts("extern uintL read_byte_array (const gcv_object_t* stream_, const gcv_object_t* bytearray_, uintL start, uintL len, perseverance_t persev);");

# Function: Writes several bytes to a stream.
# write_byte_array(&stream,&bytearray,start,len,no_hang)
# > stream: Stream (on the STACK)
# > object bytearray: simple-8bit-vector (on the STACK)
# > uintL start: start index of byte sequence to be written
# > uintL len: length of byte sequence to be written
# > perseverance_t persev: how to react on incomplete I/O
# < uintL result: number of bytes that have been written
# can trigger GC
extern maygc uintL write_byte_array (const gcv_object_t* stream_, const gcv_object_t* bytearray_, uintL start, uintL len, perseverance_t persev);
# is used by SEQUENCE
%% puts("extern uintL write_byte_array (const gcv_object_t* stream_, const gcv_object_t* bytearray_, uintL start, uintL len, perseverance_t persev);");

# Function: Reads several characters from a stream.
# read_char_array(&stream,&chararray,start,len)
# > stream: stream (on the STACK)
# > object chararray: a mutable string that is or was simple (on the STACK)
# > uintL start: start index of character sequence to be filled
# > uintL len: length of character sequence to be filled
# < uintL result: number of characters that have been filled
# can trigger GC
extern maygc uintL read_char_array (const gcv_object_t* stream_, const gcv_object_t* chararray_, uintL start, uintL len);
# is used by SEQUENCE

# Function: Writes several characters to a stream.
# write_char_array(&stream,&chararray,start,len)
# > stream: stream (on the STACK)
# > object chararray: not-reallocated simple-string (on the STACK)
# > uintL start: start index of character sequence to be written
# > uintL len: length of character sequence to be written
# can trigger GC
extern maygc void write_char_array (const gcv_object_t* stream_, const gcv_object_t* chararray_, uintL start, uintL len);
# is used by SEQUENCE

# UP: Gives the stream that is the value of a variable
# var_stream(sym,strmflags)
# > sym: Variable (symbol)
# > strmflags: Set of operations that should work on the stream
# < result: Stream
extern object var_stream (object sym, uintB strmflags);
# is used by IO, PACKAGE, ERROR, DEBUG, SPVW

/* return the file stream truename
 > s: file stream (open or closed) - no type check is done!
 < truename of the file associated with the stream
 for syscall module */
extern object file_stream_truename (object s);
%% puts("extern object file_stream_truename (object s);");

# UP: makes a file-stream
# make_file_stream(direction,append_flag,handle_fresh)
# > STACK_5: Filename, a Pathname or NIL
# > STACK_4: Truename, a Pathname or NIL
# > STACK_3: :BUFFERED argument
# > STACK_2: :EXTERNAL-FORMAT argument
# > STACK_1: :ELEMENT-TYPE argument
# > STACK_0: Handle of the open file
# > direction: Mode (0 = :PROBE, 1 = :INPUT, 4 = :OUTPUT, 5 = :IO,
#                    3 = :INPUT-IMMUTABLE)
# > append_flag: true if the stream should immediately be positioned at the end
#                ,else false
# > handle_fresh: whether the handle is freshly created.
#                 This means 1. that it is currently positioned at position 0,
#                 2. if (direction & bit(2)), it is opened for read/write, not
#                 only for write.
#                 If the handle refers to a regular file, this together means
#                 that it supports file_lseek, reading/repositioning/writing
#                 and close/reopen.
# If direction==5, handle_fresh must be true.
# < result: File-Stream (or evtl. File-Handle-Stream)
# < STACK: cleaned
# can trigger GC
extern maygc object make_file_stream (direction_t direction, bool append_flag, bool handle_at_pos_0);
# is used by PATHNAME
%% puts("extern object make_file_stream (direction_t direction, bool append_flag,bool handle_fresh);/*+6 arguments on the STACK!*/");

/* check whether the object is a handle stream or a socket-server
 and return its socket-like handle(s) */
extern void stream_handles (object obj, bool check_open, bool* char_p, SOCKET* in_sock, SOCKET* out_sock);
%% puts("extern void stream_handles (object obj, bool check_open, bool* char_p, SOCKET* in_sock, SOCKET* out_sock);");

#ifdef PIPES
/* mkops_from_handles(pipe,process_id)
   Make a PIPE-OUTPUT-STREAM from pipe handle and a process-id
   > STACK_0: buffered
   > STACK_1: element-type
   > STACK_2: encoding
   < STACK_0: result - a PIPE-OUTPUT-STREAM
   Used in LAUNCH
   Can trigger GC */
extern maygc void mkops_from_handles (Handle opipe, int process_id);
# is used by PATHNAME
#endif

#ifdef PIPES
/* mkips_from_handles(pipe,process_id)
   Make a PIPE-INPUT-STREAM from pipe handle and a process-id
   > STACK_0: buffered
   > STACK_1: element-type
   > STACK_2: encoding
   < STACK_0: result - a PIPE-INPUT-STREAM
   Used in LAUNCH
   Can trigger GC */
extern maygc void mkips_from_handles (Handle ipipe, int process_id);
# is used by PATHNAME
#endif

# Makes a Broadcast-Stream using a Stream stream.
# make_broadcast1_stream(stream)
# can trigger GC
extern maygc object make_broadcast1_stream (object stream);
# is used by IO

# Makes a Two-Way-stream using an Input-Stream and an Output-Stream.
# make_twoway_stream(input_stream,output_stream)
# > input_stream : Input-Stream
# > output_stream : Output-Stream
# < result : Two-Way-Stream
# can trigger GC
extern maygc object make_twoway_stream (object input_stream, object output_stream);
# is used by SPVW

# Makes a string-output-stream.
# make_string_output_stream()
# can trigger GC
extern maygc object make_string_output_stream (void);
# is used by IO, EVAL, DEBUG, ERROR

# UP: Returns the collected contents of a String-Output-Stream.
# get_output_stream_string(&stream)
# > stream: String-Output-Stream
# < stream: emptied Stream
# < result: the aggregation, a Simple-String
# can trigger GC
extern maygc object get_output_stream_string (const gcv_object_t* stream_);
# is used by IO, EVAL, DEBUG, ERROR

# UP: Makes a pretty-printer help stream
# make_pphelp_stream()
# can trigger GC
extern maygc object make_pphelp_stream (void);
# is used by IO

# UP: Tells whether a stream is buffered.
# stream_isbuffered(stream)
# > stream: a channel or socket stream
# < result: bit(1) set if input side is buffered,
#           bit(0) set if output side is buffered
extern uintB stream_isbuffered (object stream);
# is used by IO

# UP: Returns the current line number of a stream.
# stream_line_number(stream)
# > stream: a stream
# < result: an integer or NIL
# can trigger GC
extern maygc object stream_line_number (object stream);
# is used by IO

# Function: Returns the last character read (and not yet unread) from a stream.
# stream_get_lastchar(stream)
# > stream: a stream
# < result: the last character read, or NIL
# can trigger GC
extern maygc object stream_get_lastchar (object stream);
# is used by DEBUG

# Function: Returns true if a stream allows read-eval.
# stream_get_read_eval(stream)
# > stream: a stream
# < result: true if read-eval is allowed from the stream, else false
extern maygc bool stream_get_read_eval (object stream);
# used by IO

# Function: Changes the read-eval state of a stream.
# stream_set_read_eval(stream,value);
# > stream: a stream
# > value: true if read-eval shall be allowed from the stream, else false
extern maygc void stream_set_read_eval (object stream, bool value);
# used by

#if defined(UNIX) && !defined(NEXTAPP)
  # UP: return terminal to normal mode
  # terminal_sane();
  extern void terminal_sane (void);
  # is used by SPVW
#endif

# Function: test whether a stream is a terminal stream.
extern bool terminal_stream_p(object stream);

# check whether the charset is valid
# signal an error when code is invalid and charset is not nullobj
# return false otherwise
extern bool check_charset (const char * code, object charset);
# used in encoding.d

# ###################### SOCKBIBL for SOCKET.D ############################# #

#if defined(UNIX) || defined(WIN32_NATIVE)
/* Convert the IP address from C format to Lisp
 > type: address type (AF_INET..)
 > addr: whatever the address is for this type
 < lisp string representing the address in a human-readable format
 for syscalls & rawsock modules
 can trigger GC */
extern maygc object addr_to_string (short type, char *addr);
#endif
%% #if defined(UNIX) || defined(WIN32_NATIVE)
%%   puts("extern object addr_to_string (short type, char *addr);");
%% #endif

#if (defined(UNIX) || defined(WIN32_NATIVE)) && defined(HAVE_GETHOSTBYNAME)
/* A wrapper around the connect() function.
 To be used inside begin/end_system_call() only. */
extern int nonintr_connect (SOCKET fd, struct sockaddr * name, int namelen);
#endif

#if (defined(UNIX) || defined(WIN32_NATIVE)) && defined(HAVE_GETHOSTBYNAME) && defined(TCPCONN)
/* Convert the IP address from C format to Lisp
 > name: FQDN or dotted quad or IPv6 address
 < lisp string for FQDN or integer for IPv[46] numerics
 for syscalls & rawsock modules
 can trigger GC */
extern maygc object string_to_addr (const char* name);
#endif
%% #if (defined(UNIX) || defined(WIN32_NATIVE)) && defined(HAVE_GETHOSTBYNAME) && defined(TCPCONN)
%%   puts("extern object string_to_addr (const char *name);");
%% #endif

#if (defined(UNIX) || defined(WIN32_NATIVE)) && defined(HAVE_GETHOSTBYNAME) && defined(TCPCONN)
/* Return the hostent specified by the host designator
 > arg: host name designator:
        :DEFAULT - current host
        string/symbol: FQDN is resolved (gethostbyname)
        uint32: raw IPv4 address (gethostbyaddr)
        uint128: raw IPv6 address (gethostbyaddr)
        bit vector: raw IPv4[46] address (gethostbyaddr)
 < static hostent descriptor from LIBC
 for syscalls & rawsock modules */
extern struct hostent* resolve_host (object arg);
#endif
%% #if (defined(UNIX) || defined(WIN32_NATIVE)) && defined(HAVE_GETHOSTBYNAME) && defined(TCPCONN)
%%   puts("extern struct hostent* resolve_host (object arg);");
%% #endif

#if (defined(UNIX) || defined(WIN32_NATIVE)) && defined(HAVE_GETHOSTBYNAME)
/* connect_to_x_server(host,display)
 Attempts to connect to server, given host name and display number.
 Returns file descriptor (network socket). Returns -1 and sets errno
 if connection fails.
 An empty hostname is interpreted as the most efficient local connection to
 a server on the same machine (usually a UNIX domain socket).
 hostname="unix" is interpreted as a UNIX domain connection.
 To be used inside begin/end_system_call() only. */
extern SOCKET connect_to_x_server (const char* host, int display);
#endif

#if (defined(UNIX) || defined(WIN32_NATIVE)) && defined(HAVE_GETHOSTBYNAME) && defined(SOCKET_STREAMS)
/* socket_getlocalname(socket_handle,hd)
 Return the IP name of the localhost for the given socket.
 Fills all of *hd.
 To be used inside begin/end_system_call() only. */
extern host_data_t * socket_getlocalname (SOCKET socket_handle, host_data_t * hd, bool resolve_p);
#endif

#if (defined(UNIX) || defined(WIN32_NATIVE)) && defined(HAVE_GETHOSTBYNAME) && defined(SOCKET_STREAMS)
/* socket_getpeername(socket_handle,hd)
 Returns the name of the host to which IP socket fd is connected.
 Fills all of *hd.
 To be used inside begin/end_system_call() only. */
extern host_data_t * socket_getpeername (SOCKET socket_handle, host_data_t * hd, bool resolve_p);
#endif

#if (defined(UNIX) || defined(WIN32_NATIVE)) && defined(HAVE_GETHOSTBYNAME) && defined(SOCKET_STREAMS)
/* Creates a socket to which other processes can connect. */
extern SOCKET create_server_socket_by_socket (host_data_t *hd, SOCKET sock,
                                              unsigned int port, int backlog);

extern SOCKET create_server_socket_by_string (host_data_t *hd,
                                              const char *ip_interface,
                                              unsigned int port, int backlog);
#endif

#if (defined(UNIX) || defined(WIN32_NATIVE)) && defined(HAVE_GETHOSTBYNAME) && defined(SOCKET_STREAMS)
/* Waits for a connection to another process. */
extern SOCKET accept_connection (SOCKET socket_handle);
#endif

#if (defined(UNIX) || defined(WIN32_NATIVE)) && defined(HAVE_GETHOSTBYNAME) && defined(SOCKET_STREAMS)
/* Creates a connection to a server (which must be waiting
   on the specified host and port). */
extern SOCKET create_client_socket (const char* hostname, unsigned int port, void* timeout);
#endif

# ####################### SYMBIBL for SYMBOL.D ############################# #

# UP: Returns the gobal definition of a symbol's function,
# and tests, whether the symbol is a global function.
# Symbol_function_checked(symbol)
# > symbol: symbol
# < result: the global definition of the function
extern object Symbol_function_checked (object symbol);
# is used by

# UP: gets a property from a symbol's property list.
# get(symbol,key)
# > symbol: a symbold
# > key: a key that is comparable with EQ
# < value: corresponding value from the property list of 'symbol', or unbound.
extern object get (object symbol, object key);
# is used by IO, CONTROL, EVAL, PREDTYPE, SEQUENCE
%% #if notused
%% puts("extern object get (object symbol, object key);");
%% #endif

# ##################### ARITBIBL for LISTARIT.D ############################ #

/* check_real(obj)
 > obj: an object
 < result: a real number, either the same as obj or a replacement
 can trigger GC */
extern maygc object check_real_replacement (object obj);
#ifndef COMPILE_STANDALONE
static inline maygc object check_real (object obj) {
  if_realp(obj, ; , { obj = check_real_replacement(obj); });
  return obj;
}
#endif
/* used by IO */

# UP: Initializes the arithmetics.
# init_arith();
# can trigger GC
extern maygc void init_arith (void);
# is used by SPVW

# Converts a longword into an Integer.
# L_to_I(wert)
# > wert: value of the Integer, a signed 32-Bit-Integer.
# < result: Integer with that value.
# can trigger GC
extern maygc object L_to_I (sint32 wert);
/* used by TIME */
%% puts("extern object L_to_I (sint32 wert);");

# Converts an unsigned longword into an Integer >=0.
# UL_to_I(wert)
# > wert: value of the Integer, an unsigned 32-bit-Integer.
# < result: Integer with that value.
# can trigger GC
#if (intLsize<=oint_data_len)
  #ifdef DEBUG_GCSAFETY
    static inline maygc object UL_to_I (uintL wert) {
      return fixnum(wert);
    }
  #else
    #define UL_to_I(wert)  fixnum((uintL)(wert))
  #endif
#else
  extern maygc object UL_to_I (uintL wert);
#endif
# is used by MISC, TIME, STREAM, PATHNAME, HASHTABL, SPVW, ARRAY
%% #if (intLsize<=oint_data_len)
%%   export_def(UL_to_I(wert));
%% #else
%%   puts("extern object UL_to_I (uintL wert);");
%% #endif

# converts a double-longword into an Integer.
# L2_to_I(wert_hi,wert_lo)
# > wert_hi|wert_lo: value of the Integer, an signed 64-bit-Integer.
# < result: Integer with that value.
# can trigger GC
#if (intVsize>32)
  #define L2_to_I(wert_hi,wert_lo)  \
    Q_to_I(((sint64)(sint32)(wert_hi)<<32)|(sint64)(uint32)(wert_lo))
#else
  extern maygc object L2_to_I (sint32 wert_hi, uint32 wert_lo);
#endif
# is used by TIME, FOREIGN
%% #if (intVsize>32)
%%   export_def(L2_to_I(wert_hi,wert_lo));
%% #else
%%   puts("extern object L2_to_I (sint32 wert_hi, uint32 wert_lo);");
%% #endif

# Converts an unsigned double-longword into an Integer.
# UL2_to_I(wert_hi,wert_lo)
# > wert_hi|wert_lo: value of the Integer, an unsigned 64-bit-Integer.
# < result: Integer with that value.
# can trigger GC
#if (intVsize>32)
  #define UL2_to_I(wert_hi,wert_lo)  \
    UQ_to_I(((uint64)(uint32)(wert_hi)<<32)|(uint64)(uint32)(wert_lo))
#else
  extern maygc object UL2_to_I (uint32 wert_hi, uint32 wert_lo);
#endif
# is used by TIME, FOREIGN, and by the FFI
%% #if (intVsize>32)
%%   export_def(UL2_to_I(wert_hi,wert_lo));
%% #else
%%   puts("extern object UL2_to_I (uint32 wert_hi, uint32 wert_lo);");
%% #endif

#if defined(intQsize) || (intVsize>32)
  # Converts a quadword into an Integer.
  # Q_to_I(wert)
  # > wert: value of the Integer, a signed 64-bit-Integer.
  # < result: Integer with that value
  # can trigger GC
  extern maygc object Q_to_I (sint64 wert);
  # is used by the FFI
#endif
%% #if defined(intQsize) || (intVsize>32)
%%   puts("extern object Q_to_I (sint64 wert);");
%% #endif

#if defined(intQsize) || (intVsize>32) || defined(WIDE_HARD) || (SIZEOF_OFF_T > 4) || (SIZEOF_INO_T > 4)
  # Converts an unsigned quadword into an Integer >=0.
  # UQ_to_I(wert)
  # > wert: value of the Integer, an unsigned 64-bit-Integer.
  # < result: Integer with that value
  # can trigger GC
  extern maygc object UQ_to_I (uint64 wert);
  # is used by MISC, TIME, FFI
#endif
%% #if defined(intQsize) || (intVsize>32)
%%   puts("extern object UQ_to_I (uint64 wert);");
%% #endif

# Converts a sintV into an Integer.
# V_to_I(wert)
# > wert: value of the Integer, a signed intVsize-bit-Integer.
# < result: Integer with that value
# can trigger GC
# extern maygc object V_to_I (uintV wert);
#if (intVsize<=32)
  #define V_to_I(wert)  L_to_I(wert)
#else
  #define V_to_I(wert)  Q_to_I(wert)
#endif
# is used by LISPARIT
%% #if notused
%% #if (intVsize<=32)
%%   emit_define("V_to_I(wert)","L_to_I(wert)");
%% #else
%%   emit_define("V_to_I(wert)","Q_to_I(wert)");
%% #endif
%% #endif

# Converts an uintV into an Integer >=0.
# UV_to_I(wert)
# > wert: value of the Integer, an unsigned intVsize-bit-Integer.
# < result: Integer with that value
# can trigger GC
# extern maygc object UV_to_I (uintV wert);
#if (intVsize<=32)
  #define UV_to_I(wert)  UL_to_I(wert)
#else
  #define UV_to_I(wert)  UQ_to_I(wert)
#endif
# is used by LISPARIT
%% #if notused
%% #if (intVsize<=32)
%%   emit_define("UV_to_I(wert)","UL_to_I(wert)");
%% #else
%%   emit_define("UV_to_I(wert)","UQ_to_I(wert)");
%% #endif
%% #endif

# Converts a C-Integer of a given type into an Integer
# val should be a variable
#define uint8_to_I(val)  fixnum((uint8)(val))
#define sint8_to_I(val)  L_to_I((sint32)(sint8)(val))
#define uint16_to_I(val)  fixnum((uint16)(val))
#define sint16_to_I(val)  L_to_I((sint32)(sint16)(val))
#define uint32_to_I(val)  UL_to_I((uint32)(val))
#define sint32_to_I(val)  L_to_I((sint32)(val))
#if defined(intQsize) || (intVsize>32)
  #define uint64_to_I(val)  UQ_to_I((uint64)(val))
  #define sint64_to_I(val)  Q_to_I((sint64)(val))
#elif defined(HAVE_FFI)
  #define uint64_to_I(val)  UL2_to_I((uint32)((val)>>32),(uint32)(val))
  #define sint64_to_I(val)  L2_to_I((sint32)((val)>>32),(uint32)(val))
#endif
#if (int_bitsize==16)
  #define uint_to_I(val)  uint16_to_I(val)
  #define sint_to_I(val)  sint16_to_I(val)
#else # (int_bitsize==32)
  #define uint_to_I(val)  uint32_to_I(val)
  #define sint_to_I(val)  sint32_to_I(val)
#endif
#if (long_bitsize==32)
  #define ulong_to_I(val)  uint32_to_I(val)
  #define slong_to_I(val)  sint32_to_I(val)
#else # (long_bitsize==64)
  #define ulong_to_I(val)  uint64_to_I(val)
  #define slong_to_I(val)  sint64_to_I(val)
#endif
# is used by MISC, for FFI
%% export_def(uint8_to_I(val));
%% export_def(sint8_to_I(val));
%% export_def(uint16_to_I(val));
%% export_def(sint16_to_I(val));
%% export_def(uint32_to_I(val));
%% export_def(sint32_to_I(val));
%% export_def(uint64_to_I(val));
%% export_def(sint64_to_I(val));
%% export_def(uint_to_I(val));
%% export_def(sint_to_I(val));
%% export_def(ulong_to_I(val));
%% export_def(slong_to_I(val));

# Converts a uintM integer into an Integer.
#if intMsize <= intLsize
  #define uintM_to_I(val)  UL_to_I(val)
#else
  #define uintM_to_I(val)  UQ_to_I(val)
#endif

# Converts a sintM integer into an Integer.
#if intMsize <= intLsize
  #define sintM_to_I(val)  L_to_I(val)
#else
  #define sintM_to_I(val)  Q_to_I(val)
#endif

/* converts off_t to an Integer */
#if defined(WIN32_NATIVE)
  #define off_to_I(val)  L2_to_I((sint32)((val)>>32),(uint32)(val))
#elif SIZEOF_OFF_T > 4
  #define off_to_I  Q_to_I
#else
  #define off_to_I  L_to_I
#endif

# Converts an Integer >=0 into an unsigned longword.
# I_to_UL(obj)
# > obj: an object, should be an Integer >=0, <2^32
# < result: the Integer's value as unsigned longword
extern uintL I_to_UL (object obj);
# is used by TIME, ARRAY
%% puts("extern uintL I_to_UL (object obj);");

# Converts an Integer into a signed longword
# I_to_L(obj)
# > obj: an object, should be an Integer >=-2^31, <2^31
# < result: the Integer's value as signed longword
extern sintL I_to_L (object obj);
# is used by
%% puts("extern sintL I_to_L (object obj);");

#if defined(HAVE_LONGLONG)
  # Converts an Integer >=0 into an unsigned quadword.
  # I_to_UQ(obj)
  # > obj: an object, should be an Integer >=0, <2^64
  # < result: the Integer's vaulue as unsigned quadword
  extern uint64 I_to_UQ (object obj);
  /* used by FOREIGN, for FFI, and by modules */
#endif
%% #ifdef HAVE_LONGLONG
%%   puts("extern uint64 I_to_UQ (object obj);");
%% #endif

#if defined(HAVE_LONGLONG)
  # Converts an Integer into a signed quadword.
  # I_to_Q(obj)
  # > obj: an object, should be an Integer >=-2^63, <2^63
  # < result: the Integer's value as quadword.
  extern sint64 I_to_Q (object obj);
  /* used by FOREIGN, for FFI, and by modules */
#endif
%% #ifdef HAVE_LONGLONG
%%   puts("extern sint64 I_to_Q (object obj);");
%% #endif

# Converts an Integer into a C-Integer of a given type.
# I_to_xintyy(obj) assumes that xintyy_p(obj) has already been checked.
#define I_to_uint8(obj)  (uint8)(as_oint(obj) >> oint_data_shift)
#define I_to_sint8(obj)  (sint8)(as_oint(obj) >> oint_data_shift)
#define I_to_uint16(obj)  (uint16)(as_oint(obj) >> oint_data_shift)
#define I_to_sint16(obj)  (sint16)(as_oint(obj) >> oint_data_shift)
#if (oint_data_len>=32)
  #define I_to_uint32(obj)  (uint32)(as_oint(obj) >> oint_data_shift)
#else
  #define I_to_uint32(obj)  I_to_UL(obj)
#endif
#if (oint_data_len>=31)
  #define I_to_sint32(obj)  (sint32)(as_oint(obj) >> oint_data_shift)
#else
  #define I_to_sint32(obj)  I_to_L(obj)
#endif
#ifdef HAVE_LONGLONG
  #define I_to_uint64(obj)  I_to_UQ(obj)
  #define I_to_sint64(obj)  I_to_Q(obj)
#endif
#if (int_bitsize==16)
  #define I_to_uint  I_to_uint16
  #define I_to_sint  I_to_sint16
#else # (int_bitsize==32)
  #define I_to_uint  I_to_uint32
  #define I_to_sint  I_to_sint32
#endif
#if defined(HAVE_FFI) || defined(HAVE_AFFI)
  #if (long_bitsize==32)
    #define I_to_ulong  I_to_uint32
    #define I_to_slong  I_to_sint32
  #else # (long_bitsize==64)
    #define I_to_ulong  I_to_uint64
    #define I_to_slong  I_to_sint64
  #endif
#endif
/* used by FFI, STREAM, modules */
%% export_def(I_to_uint8(obj));
%% export_def(I_to_sint8(obj));
%% export_def(I_to_uint16(obj));
%% export_def(I_to_sint16(obj));
%% export_def(I_to_uint32(obj));
%% export_def(I_to_sint32(obj));
%% export_def(I_to_uint64(obj));
%% export_def(I_to_sint64(obj));
%% export_def(I_to_uint);
%% export_def(I_to_sint);
%% export_def(I_to_ulong);
%% export_def(I_to_slong);

/* Unsigned Digit Sequence to Integer
 UDS_to_I(MSDptr,len)
 convert UDS MSDptr/len/.. into Integer >=0 .
 MSDptr[0] is the most significant digit, MSDptr[len-1] the least significant.
 there must be room for 1 digit below of MSDptr.
 can trigger GC */
extern maygc object UDS_to_I (uintD* MSDptr, uintC len);
/* is used by modules */
%% puts("extern object UDS_to_I (uintD* MSDptr, uintC len);");

/* Digit Sequence to Integer
 DS_to_I(MSDptr,len)
 convert DS MSDptr/len/.. into Integer.
 MSDptr[0] is the most significant digit, MSDptr[len-1] the least significant.
 can trigger GC */
extern maygc object DS_to_I (const uintD* MSDptr, uintC len);
/* is used by modules */
%% puts("extern object DS_to_I (const uintD* MSDptr, uintC len);");

# I_I_comp(x,y) compares two Integers x and y.
# Result: 0 if x=y, +1 if x>y, -1 if x<y.
extern signean I_I_comp (object x, object y);
# is used by SEQUENCE

# (1+ x), where x is an Integer. Result Integer.
# I_1_plus_I(x)
# can trigger GC
extern maygc object I_1_plus_I (object x);
# is used by SEQUENCE, SPVW, SYMBOL
%% #if notused
%% puts("extern object I_1_plus_I (object x);");
%% #endif

# (1- x), where x is an Integer. Result Integer.
# I_minus1_plus_I(x)
# can trigger GC
extern maygc object I_minus1_plus_I (object x);
# is used by SEQUENCE
%% #if notused
%% puts("extern object I_minus1_plus_I (object x);");
%% #endif

# (+ x y), where x and y are Integers. Result Integer.
# I_I_plus_I(x,y)
# can trigger GC
extern maygc object I_I_plus_I (object x, object y);
# is used by SEQUENCE
%% #if notused
%% puts("extern object I_I_plus_I (object x, object y);");
%% #endif

# (- x y), where x and y are Integers. Result Integer.
# I_I_minus_I(x,y)
# can trigger GC
extern maygc object I_I_minus_I (object x, object y);
# is used by SEQUENCE
%% #if notused
%% puts("extern object I_I_minus_I (object x, object y);");
%% #endif

# (ASH x y), where x and y are Integers. Result Integer.
# I_I_ash_I(x,y)
# can trigger GC
extern maygc object I_I_ash_I (object x, object y);
# is used by SEQUENCE

# (INTEGER-LENGTH x), where x is an Integer. Result uintL.
# I_integer_length(x)
extern uintL I_integer_length (object x);
# is used by ARRAY
%% puts("extern uintL I_integer_length (object x);");

# Converts a little-endian byte sequence to an unsigned integer.
# > bytesize: number of given 8-bit bytes of the integer,
#             < intDsize/8*uintWC_max
# > bufferptr: address of bytesize bytes in GC-invariant memory
# < result: an integer >= 0 with I_integer_length(result) <= 8*bytesize
# can trigger GC
extern maygc object LEbytes_to_UI (uintL bytesize, const uintB* bufferptr);
%% puts("extern object LEbytes_to_UI (uintL bytesize, const uintB* bufferptr);");

# Converts a little-endian byte sequence to an unsigned integer.
# > bytesize: number of given 8-bit bytes of the integer,
#             < intDsize/8*uintWC_max
# > *buffer_: address of a simple-8bit-vector (or of a fake)
#             containing bytesize bytes of memory
# < result: an integer >= 0 with I_integer_length(result) <= 8*bytesize
# can trigger GC
extern maygc object LESbvector_to_UI (uintL bytesize, const gcv_object_t* buffer_);
# is used by STREAM

# Converts a little-endian byte sequence to an integer.
# > bytesize: number of given 8-bit bytes of the integer, > 0,
#             < intDsize/8*uintWC_max
# > bufferptr: address of bytesize bytes in GC-invariant memory
# < result: an integer with I_integer_length(result) < 8*bytesize
# can trigger GC
extern maygc object LEbytes_to_I (uintL bytesize, const uintB* bufferptr);
%% puts("extern object LEbytes_to_I (uintL bytesize, const uintB* bufferptr);");

# Converts a little-endian byte sequence to an integer.
# > bytesize: number of given 8-bit bytes of the integer, > 0,
#             < intDsize/8*uintWC_max
# > *buffer_: address of a simple-8bit-vector (or of a fake)
#             containing bytesize bytes of memory
# < result: an integer with I_integer_length(result) < 8*bytesize
# can trigger GC
extern maygc object LESbvector_to_I (uintL bytesize, const gcv_object_t* buffer_);
# is used by STREAM

# Converts an unsigned integer to a little-endian byte sequence.
# > obj: an integer
# > bitsize: maximum number of bits of the integer
# > bufferptr: pointer to bytesize = ceiling(bitsize,8) bytes of memory
# < false and bufferptr[0..bytesize-1] filled, if obj >= 0 and
#                                              I_integer_length(obj) <= bitsize;
#   true, if obj is out of range
extern bool UI_to_LEbytes (object obj, uintL bitsize, uintB* bufferptr);
# is used by STREAM
%% puts("extern bool UI_to_LEbytes (object obj, uintL bitsize, uintB* bufferptr);");

# Converts an integer to a little-endian byte sequence.
# > obj: an integer
# > bitsize: maximum number of bits of the integer, including the sign bit
# > bufferptr: pointer to bytesize = ceiling(bitsize,8) bytes of memory
# < false and bufferptr[0..bytesize-1] filled, if I_integer_length(obj) < bitsize;
#   true, if obj is out of range
extern bool I_to_LEbytes (object obj, uintL bitsize, uintB* bufferptr);
# is used by STREAM
%% puts("extern bool I_to_LEbytes (object obj, uintL bitsize, uintB* bufferptr);");

# c_float_to_FF(&val) converts an IEEE-single-float val into an single-float.
# can trigger GC
extern maygc object c_float_to_FF (const ffloatjanus* val_);
%% puts("extern object c_float_to_FF (const ffloatjanus* val_);");

# FF_to_c_float(obj,&val);
# converts single-float obj into an IEEE-single-float val.
extern void FF_to_c_float (object obj, ffloatjanus* val_);
%% puts("extern void FF_to_c_float (object obj, ffloatjanus* val_);");

# c_double_to_DF(&val) converts an IEEE-double-float val into a double-float.
# can trigger GC
extern maygc object c_double_to_DF (const dfloatjanus* val_);
%% puts("extern object c_double_to_DF (const dfloatjanus* val_);");

# DF_to_c_double(obj,&val);
# converts a double-float obj into an IEEE-double-float val.
extern void DF_to_c_double (object obj, dfloatjanus* val_);
%% puts("extern void DF_to_c_double (object obj, dfloatjanus* val_);");

/* hash-code of a Long-Float: mixture of exponent, length, first 32 bits */
extern uint32 hashcode_lfloat (object obj);

/* (complex x (float 0 x)) */
extern object F_complex_C (object x);

# UP: turns a string with Integer syntax into an Integer number
# Points will be ignored
# read_integer(base,sign,string,index1,index2)
# > base: read base(>=2, <=36)
# > sign: sign (/=0 if negative)
# > string: simple-string (contains digits with values
#   <base and eventually a point)
# > index1: Index of the first digit
# > index2: Index after the last digit
#   (thus index2-index1 digits, incl. a decimal point that can be at the end)
# < result: Integer
# can trigger GC
extern maygc object read_integer (uintWL base, signean sign, object string, uintL index1, uintL index2);
# is used by IO

# UP: turns a string with rational syntax into a rational number
# read_rational(base,sign,string,index1,index3,index2)
# > base: read base (>=2, <=36)
# > sign: sign (/=0 if negative)
# > string: Normal-Simple-String (contains digits with values
#           <base and fraction bar)
# > index1: Index of the first digit
# > index3: Index of '/'
# > index2: Index after the last digit
#   (thus index3-index1 digits of the numerator,
#    index2-index3-1 digits of the denominator)
# < result: rational number
# can trigger GC
extern maygc object read_rational (uintWL base, signean sign, object string, uintL index1, uintL index3, uintL index2);
# is used by IO

# UP: turns a string with float-syntax into a float
# read_float(base,sign,string,index1,index4,index2,index3)
# > base: read base (=10)
# > sign: Sign (/=0 if negative)
# > string: normal-simple-string (contains digits and eventually
#           point and exponent marker)
# > index1: Index of the beginning of the mantissa (without sign)
# > index4: Index after the end of the mantissa
# > index2: Index at the end of the character
# > index3: Index after the decimal point (=index4 if there is none)
#   (thus mantissa with index4-index1 characters: digits and max. 1 '.')
#   (thus index4-index3 digits after the decimal point)
#   (thus at index4<index2: index4 = index of the exponent marker,
#    index4+1 = index of the exponent's sign or the first digit
#    of the exponent)
# < result: Float
# can trigger GC
extern maygc object read_float (uintWL base, signean sign, object string, uintL index1, uintL index4, uintL index2, uintL index3);
# is used by IO

# UP: prints an Integer
# print_integer(z,base,&stream);
# > z: Integer
# > base: base (>=2, <=36)
# > stream: Stream
# < stream: Stream
# can trigger GC
extern maygc void print_integer (object z, uintWL base, const gcv_object_t* stream_);
# is used by IO

# UP: prints a float
# print_float(z,&stream);
# > z: Float
# > stream: Stream
# < stream: Stream
# can trigger GC
extern maygc void print_float (object z, const gcv_object_t* stream_);
# is used by IO

# UP: Multiply an Integer by 10 and add another digit
# mal_10_plus_x(y,x)
# > y: Integer Y (>=0)
# > x: digit value X (>=0,<10)
# < result: Integer Y*10+X (>=0)
# can trigger GC
extern maygc object mal_10_plus_x (object y, uintB x);
# is used by IO

# UP: decides whether two numbers are equal
# number_gleich(x,y)
# > x,y: two numbers
# < result: true, if (= x y) holds
extern bool number_gleich (object x, object y);
# is used by PREDTYPE

# UP: Converts an object into a float of a given type
# coerce_float(obj,type)
# > obj: Object
# > type: one of the symbols
#         FLOAT, SHORT-FLOAT, SINGLE-FLOAT, DOUBLE-FLOAT, LONG-FLOAT
# < result: (coerce obj type)
# can trigger GC
extern maygc object coerce_float (object obj, object type);
# is used by PREDTYPE

/* Converts a function's argument to a C 'double'.
 to_double(obj)
 > obj: an object, usually a real number
 < result: its value as a C 'double'
 can trigger GC */
extern maygc double to_double (object x);
%% puts("extern double to_double (object obj);");

/* Converts a function's argument to a C 'int'.
 to_int(obj)
 > obj: an object, usually an integer
 < result: its value as a C 'int'
 can trigger GC */
extern maygc int to_int (object x);
%% puts("extern int to_int (object obj);");

# UP: Returns the decimal string representation of an integer >= 0.
# decimal_string(x)
# > object x: an integer >= 0
# < object result: a normal-simple-string containing the digits
# can trigger GC
extern maygc object decimal_string (object x);
# is used by PATHNAME

# ###################### FOREIGNBIBL for FOREIGN.D ########################## #

#ifdef DYNAMIC_FFI
%% #ifdef DYNAMIC_FFI

# Return the pointer encoded by a Foreign-Pointer.
  #define Fpointer_value(obj) TheFpointer(obj)->fp_pointer

# Return the pointer encoded by a Foreign-Address. obj a variable
  #define Faddress_value(obj)  \
   ((void*)((uintP)Fpointer_value(TheFaddress(obj)->fa_base) + TheFaddress(obj)->fa_offset))

/* ensure that the Faddress is valid
 < fa: foreign address (not checked!)
 can trigger GC */
extern maygc object check_faddress_valid (object fa);
/* usd by FOREIGN, MISC */

# Registers a foreign variable.
# register_foreign_variable(address,name,flags,size);
# > address: address of a variable in memory
# > name: its name
# > flags: fv_readonly for read-only variables
# > size: its size in bytes
# can trigger GC
  extern maygc void register_foreign_variable (void* address, const char * name, uintBWL flags, uintL size);
# Specifies that the variable will not be written to.
#define fv_readonly  bit(0)
# Specifies that when the value is replaced and the variable contains pointers,
# the old storage will be free()d and new storage will be allocated via malloc().
#define fv_malloc    bit(1)
%%   puts("extern void register_foreign_variable (void* address, const char * name, uintBWL flags, uintL size);");

# Registers a foreign function.
# register_foreign_function(address,name,flags);
# > address: address of the function in memory
# > name: its name
# > flags: its language and parameter passing convention
# can trigger GC
  extern maygc void register_foreign_function (void* address, const char * name, uintWL flags);
# Flags for language:
#define ff_lang_asm       bit(8)  # no argument passing conventions
#define ff_lang_c         bit(9)  # K&R C, with argument type promotions
#define ff_lang_ansi_c    bit(10) # ANSI C, without argument type promotions
# define ff_lang_pascal   bit(11) # not yet supported
#define ff_lang_stdcall   bit(15) # `stdcall' calling convention
# Varargs functions are not supported.
# Set this if pointers within the arg should point to alloca()ed data, i.e.
# have dynamic extent: are valid for this call only.
#define ff_alloca         bit(0)
# Set this if pointers within the arg should point to malloc()ed data. The
# function takes over responsibility for that storage. For return values,
# set this if free() shall be called for pointers within the resulting value.
#define ff_malloc         bit(1)
# Set this if the arg should point to a place where a return value can be
# stored.
#define ff_out            bit(4)
# Set this if the arg is also treated as a return value.
#define ff_inout          bit(5)
%%   puts("extern void register_foreign_function (void* address, const char * name, uintWL flags);");

# Convert foreign data to Lisp data.
# can trigger GC
extern maygc object convert_from_foreign (object fvd, const void* data);
%% puts("extern object convert_from_foreign (object fvd, const void* data);");

/* Convert Lisp data to foreign data. */
typedef void* converter_malloc_t (void* old_data, uintL size, uintL alignment);
%% puts("typedef void* converter_malloc_t (void* old_data, uintL size, uintL alignment);");
%% puts("extern converter_malloc_t mallocing, nomalloc;");
/* Convert Lisp data to foreign data.
   Storage is allocated through converter_malloc().
 Only the toplevel storage must already exist; its address is given.
 can trigger GC  */
extern void convert_to_foreign (object fvd, object obj, void* data, converter_malloc_t *converter_malloc);
%% puts("extern void convert_to_foreign (object fvd, object obj, void* data, converter_malloc_t *converter_malloc);");

# Initialize the FFI.
  extern maygc void init_ffi (void);
# used by SPVW

# De-Initialize the FFI.
  extern void exit_ffi (void);
# used by SPVW

#endif
%% #endif

# ######################## THREADBIBL for THREAD.D ######################### #

#ifdef MULTITHREAD

# Structure containing all the per-thread global variables.
# (We could use a single instance of this structure also in the single-thread
# model, but it would make debugging less straightforward.)
  typedef struct {
    # Most often used:
      #if !defined(STACK_register)
        gcv_object_t* _STACK;
      #endif
      #if !defined(mv_count_register)
        uintC _mv_count;
      #endif
      #if !defined(value1_register)
        object _value1;
      #endif
      #if !defined(back_trace_register)
        p_backtrace_t _back_trace;
      #endif
    # Less often used:
      #ifndef NO_SP_CHECK
        void* _SP_bound;
      #endif
      void* _STACK_bound;
      unwind_protect_caller_t _unwind_protect_to_save;
      #ifdef NEED_temp_mv_count
        uintC _temp_mv_count;
      #endif
      #ifdef NEED_temp_value1
        object _temp_value1;
      #endif
      #ifdef HAVE_SAVED_STACK
        gcv_object_t* _saved_STACK;
      #endif
      #ifdef HAVE_SAVED_mv_count
        uintC _saved_mv_count;
      #endif
      #ifdef HAVE_SAVED_value1
        object _saved_value1;
      #endif
      #ifdef HAVE_SAVED_back_trace
        p_backtrace_t _saved_back_trace;
      #endif
      #if defined(HAVE_SAVED_REGISTERS)
        struct registers * _callback_saved_registers;
      #endif
      uintC _index; # this thread's index in allthreads[]
    # Used for exception handling only:
      handler_args_t _handler_args;
      stack_range_t* _inactive_handlers;
    # Big, rarely used arrays come last:
      object _mv_space [mv_limit-1];
    # Now the lisp objects (seen by the GC).
      # The Lisp object representing this thread:
      object _lthread;
      # The lexical environment:
      gcv_environment_t _aktenv;
      # The values of per-thread symbols:
      object _symvalues[unspecified];
  } clisp_thread_t;
  #define thread_size(nsymvalues)  \
    (offsetofa(clisp_thread_t,_symvalues)+nsymvalues*sizeof(gcv_object_t))
  #define thread_objects_offset(nsymvalues)  \
    (offsetof(clisp_thread_t,_lthread))
  #define thread_objects_anz(nsymvalues)  \
    ((offsetofa(clisp_thread_t,_symvalues)-offsetof(clisp_thread_t,_lthread))/sizeof(gcv_object_t)+(nsymvalues))

# Size of a single thread's stack region. Must be a power of 2.
  #define THREAD_SP_SHIFT  22  # 4 MB should be sufficient, and leaves room
                               # for about 128 threads.
  #define THREAD_SP_SIZE  bit(THREAD_SP_SHIFT)
# Returns the stack pointer, or some address near the stack pointer.
  # Important for efficiency: Multiple calls to this function within a single
  # function must be combined to a single, inlined call. To reach this, we
  # use __asm__, not __asm__ __volatile__, and we don't use a global register
  # variable.
  #if defined(ASM_get_SP_register)
    #define roughly_SP()  \
      ({ var aint __SP; __asm__ ASM_get_SP_register(__SP); __SP; })
  #else
    #define roughly_SP()  (aint)__builtin_frame_address(0)
    # Note: If (__GNUC__ == 2) && (__GNUC_MINOR__ >= 8) && (__GNUC_MINOR__ < 95)
    # one can write
    #   #define roughly_SP()  (aint)__builtin_sp()
    # but this isn't efficient because gcc somehow knows that the stack pointer
    # varies across the function (maybe because of our register declaration?).
  #endif
# Returns a pointer to the thread structure, given the thread's stack pointer.
  #ifdef SP_DOWN
    #ifndef MORRIS_GC
      #define sp_to_thread(sp)  \
        (clisp_thread_t*)((aint)(sp) & minus_bit(THREAD_SP_SHIFT))
    #else
      # Morris GC doesn't like the backpointers to have garcol_bit set.
      #define sp_to_thread(sp)  \
        (clisp_thread_t*)((aint)(sp) & (minus_bit(THREAD_SP_SHIFT) & ~wbit(garcol_bit_o)))
    #endif
  #endif
  #ifdef SP_UP
    #define sp_to_thread(sp)  \
      (clisp_thread_t*)(((aint)(sp) | (bit(THREAD_SP_SHIFT)-1)) - 0x1FFFF)
  #endif
/* Returns a pointer to the current thread structure. */
  typedef clisp_thread_t* current_thread_function_t (void);
  local inline const current_thread_function_t current_thread;
  local inline clisp_thread_t* current_thread (void)
  { return sp_to_thread(roughly_SP()); }

#endif

# ######################## BUILTBIBL for BUILT.D ######################### #

/* Returns a multiline string containing some info about the flags with which
   the executable was built. */
extern object built_flags (void);

# ####################### FOR DEBUGGING UNDER GDB ######################## #

#ifdef GENERATIONAL_GC
# Put a breakpoint here if you want to catch CLISP just before it dies.
extern void sigsegv_handler_failed (void* address);
#endif

/* For debugging: From within gdb, type: call ext_show_stack().
   Equivalent to (ext:show-stack) from the Lisp prompt. */
extern void ext_show_stack (void);

/*************************************************************************/

