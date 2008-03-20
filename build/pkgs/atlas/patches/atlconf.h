#ifndef ATLCONF_H
   #define ATLCONF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define NOS 12
static char *osnam[NOS] =
   {"UNKNOWN", "Linux", "SunOS", "SunOS4", "OSF1", "IRIX", "AIX",
    "Win9x", "WinNT", "HPUX", "FreeBSD", "OSX"};
enum OSTYPE {OSOther=0, OSLinux, OSSunOS, OSSunOS4, OSOSF1, OSIRIX, OSAIX,
             OSWin9x, OSWinNT, OSHPUX, OSFreeBSD, OSOSX};
#define OSIsWin(OS_) (((OS_) == OSWinNT) || ((OS_) == OSWin9x))

enum ARCHFAM {AFOther=0, AFPPC, AFSPARC, AFALPHA, AFX86, AFIA64, AFMIPS};

#define NMACH 33
static char *machnam[NMACH] =
   {"UNKNOWN", "POWER3", "POWER4", "POWER5", "PPCG4", "PPCG5",
    "P5", "P5MMX", "PPRO", "PII", "PIII", "PM", "CoreSolo",
    "CoreDuo", "Core2Solo", "Core2Duo", "P4", "P4D", "P4E", "Efficeon", "K7",
    "HAMMER", "AMD64K10h", "UNKNOWNx86", "IA64Itan", "IA64Itan2",
    "USI", "USII", "USIII", "USIV", "UnknownUS", "MIPSR1xK", "MIPSICE9"};
enum MACHTYPE {MACHOther, IbmPwr3, IbmPwr4, IbmPwr5, PPCG4, PPCG5,
               IntP5, IntP5MMX, IntPPRO, IntPII, IntPIII, IntPM, IntCoreS,
               IntCoreDuo, IntCore2Solo, IntCore2Duo, IntP4, IntP4D, IntP4E,
               TMEff, AmdAthlon, AmdHammer, Amd64K10h, x86X, IA64Itan, IA64Itan2,
               SunUSI, SunUSII, SunUSIII, SunUSIV, SunUSX,
               MIPSR1xK, /* includes R10K, R12K, R14K, R16K */
               MIPSICE9   /* SiCortex ICE9 -- like MIPS5K */
               };
#define MachIsX86(mach_) \
   ( (mach_) >= IntP5 && (mach_) <= x86X )
#define MachIsIA64(mach_) \
   ( (mach_) >= IA64Itan && (mach_) <= IA64Itan2 )
#define MachIsUS(mach_) \
   ( (mach_) >= SunUSI && (mach_) <= SunUSX )
#define MachIsMIPS(mach_) \
   ( (mach_) >= MIPSR1xK && (mach_) <= MIPSICE9 )
#define MachIsPPC(mach_) \
   ( (mach_) >= PPCG4 && (mach_) <= PPCG5 )

static char *f2c_namestr[5] = {"UNKNOWN","Add_", "Add__", "NoChange", "UpCase"};
static char *f2c_intstr[5] =
       {"UNKNOWN", "F77_INTEGER=int", "F77_INTEGER=long",
        "F77_INTEGER=\"long long\"", "F77_INTEGER=short"};
static char *f2c_strstr[5]=
       {"UNKNOWN", "SunStyle", "CrayStyle", "StructVal", "StructPtr"};

enum F2CNAME {f2c_NamErr=0, f2c_Add_, f2c_Add__, f2c_NoChange, f2c_UpCase};
enum F2CINT {f2c_IntErr=0, FintCint, FintClong, FintClonglong, FintCshort};
enum F2CSTRING {f2c_StrErr=0, fstrSun, fstrCray, fstrStructVal, fstrStructPtr};

#define NISA 6
static char *ISAXNAM[NISA] =
   {"", "AltiVec", "SSE3", "SSE2", "SSE1", "3DNow"};
enum ISAEXT {ISA_None=0, ISA_AV, ISA_SSE3, ISA_SSE2, ISA_SSE1, ISA_3DNow};

#define NASMD 7
enum ASMDIA
   {ASM_None=0, gas_x86_32, gas_x86_64, gas_sparc, gas_ppc, gas_parisc,
    gas_mips};
static char *ASMNAM[NASMD] =
   {"",     "GAS_x8632", "GAS_x8664", "GAS_SPARC", "GAS_PPC", "GAS_PARISC",
    "GAS_MIPS"};


/*
 * Used for archinfo probes (can pack in bitfield)
 */
enum WHATPROBE{Parch=1, P64=2, Pncpu=4, Pverb=8, Pncache=16, PCacheSize=32,
               PMhz=64, Pthrottle=128};

#define NARDEF 4
enum ARDEF{ADsk=0, ADdk, ADsm, ADdm};  /* m = matmul kernel, k = non-mm kern */
/*
 * Used for all the compilers ATLAS needs
 */
#define NCOMP 7
static char *COMPNAME[NCOMP]={"ICC","SMC","DMC","SKC","DKC","XCC","F77"};
#define ICC_ 0   /* Compiles non-computation routines, and all I/O */
#define SMC_ 1   /* single prec matmul compiler */
#define DMC_ 2   /* double prec matmul compiler */
#define SKC_ 3   /* single prec computation compiler (non-mm kernels) */
#define DKC_ 4   /* double prec computation compiler */
#define XCC_ 5   /* Compiler for frontend of cross-compilation */
#define F77_ 6   /* Valid fixed-format Fortran77 compiler */

typedef struct CompNode COMPNODE;
struct CompNode
{
   int priority;              /* priority of this definition */
   int comps[1];              /* bitfield: (1<<ICC)|...|(1<<F77) */
   int OS[(NOS+31)/32];       /* bitfield for OS */
   int arch[(NMACH+31)/32];   /* bitfields for architecture */
   char *comp, *flags;        /* compiler & flags as strings */
   COMPNODE *next;
};
#include "atlconf_misc.h"

#endif
