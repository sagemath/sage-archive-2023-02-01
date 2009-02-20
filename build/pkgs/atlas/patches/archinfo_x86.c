/*
 *             Automatically Tuned Linear Algebra Software v@(ver)
 *                    (C) Copyright 2006 R. Clint Whaley
 *
 * Code contributers : R. Clint Whaley, Dean Gaudet
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
/*
 * This code written for ATLAS use by R. Clint Whaley based on code and info
 * submitted by Dean Gaudet, with the later help of the following websites:
 *   http://www.sandpile.org/ia32/cpuid.htm
 *   http://en.wikipedia.org/wiki/CPUID
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "atlconf.h"

#define uint unsigned int

/*
 * This routine returns the contents of registers set by the cpuid instruction
 * in the array res:
 *   res[0] : eax
 *   res[1] : ebx
 *   res[2] : ecx
 *   res[3] : edx
 */
void do_cpuid(uint *res, uint level);
/* result defines */
#define EAX 0
#define EBX 1
#define ECX 2
#define EDX 3


/* My driver, based on Dean's */
int ProbeArch(char *vendor, unsigned *family, unsigned *model, int *x86_64)
/*
 * Returns 0 on success, non-zero on error
 */
{
   uint r[4];
   uint max_level;
   uint *vp = (uint*) vendor;

   *x86_64 = 0;
/*
 * In this call, we ask for max supported cpuid support, and return if
 * we can't get any usuable info.  Also sets ebx,edx and ecx (16 chars of data)
 * to vendor ID string
 */
   do_cpuid(r, 0);
   max_level = r[EAX];
   if (!max_level)
      return(1);
/*
 * Copy vendor string as 3 ints rather than 16 char, then null-term at 12
 */
   *vp = r[EBX];
   vp[1] = r[EDX];
   vp[2] = r[ECX];
   vendor[12] = '\0';

/*
 * Find processor family and model, ouput EAX
 * According to latest docs, extended family and model should always be
 * added in, not just in the cases shown in the commented-out if statements
 * below.  The original "only do it in certain cases" was from the official
 * IA32 ISA, but doing this causes problems on Xeons, so now we do like the
 * newer docs indicate and always add the extended values in
 */
   do_cpuid(r, 1);
   *family = (r[EAX] >> 8) & 0xf;      /* base family in bits 11-8 */
/*   if (*family == 0xf || *family == 0) */ /* extended family is added in */
       *family += ((r[EAX] >> 20) & 0xff);

   *model = (r[0] >> 4) & 0xf;         /* model in bits 7-4 */
/* if (*model == 0xf) */               /* extended model is concatenated */
      *model |= ((r[0] >> 12) & 0xf0);

/*
 * Find out if we have extended cpuid level, and if so, see if we've got
 * x86-64 capability or not
 */
   do_cpuid(r, 0x80000000);
   if (r[0] >= 0x80000001)
   {
      do_cpuid(r, 0x80000001);
      *x86_64 = (r[EDX] & (1<<29)) != 0;   /* x86-64 in bit 29 */
   }
   return(0);
}

/*
 * constants used to check family + extended family
 */
#define EF_486       4        /* also AMD 5x86 and Cyrix 5x86 */
#define EF_P5        5        /* P5, K5 and K6 */
#define EF_P6        6        /* P6, Core and K7 (athlon) */
#define EF_ITAN0     7        /* Itanium */
#define EF_K8_P4_EFF 0x00F    /* P4, Hammer, Efficien */
#define EF_K8_ITAN   0x01F    /* Hammer, Itanium */
#define EF_K8        0x02F    /* Hammer */
#define EF_ITAN      0x020    /* Itanium */
#define EF_K8b      16        /* 3rd gen opteron */

enum FAM {ERR,         /* cannot decipher */
          i486,        /* 486 & AMD 5x86 and Cyrix 5x86 */
          P5,          /* Original Pentium and AMD K5 & K6 */
          P6,          /* Intel PIII, Core and AMD K7 (orig athlon) */
          P7,          /* Intel P4, AMD hammer, Efficeon */
          P8B,         /* 3rd generation hammer */
          ITAN};       /* Intel Itanium */

enum FAM GetFamily(int efam)  /* efam = (family+ext fam) from cpuid */
/*
 * Translates CPUID (family+extended family) to FAM enum type
 */
{
   enum FAM iret;
   switch (efam)
   {
   case EF_486:               /* also AMD 5x86 and Cyrix 5x86 */
      iret = i486;
      break;
   case EF_P5:                /* P5, K5 and K6 */
      iret = P5;
      break;
   case EF_P6:                /* P6, Core and K7 (athlon) */
      iret = P6;
      break;
   case EF_K8_P4_EFF:         /* P4, Hammer, Efficien */
      iret = P7;
      break;
   case EF_K8_ITAN:           /* Hammer, Itanium */
   case EF_K8:                /* Hammer */
      iret = P7;
      break;
   case EF_K8b:
      iret = P8B;
      break;
   case EF_ITAN:              /* Itanium */
   case EF_ITAN0:             /* Itanium */
      iret = ITAN;
      break;
   default:
      iret = ERR;
   }
   return (iret);
}

enum VEND {VERR, Intel, AMD, TM};
enum VEND str2vend(char *vendor)
/*
 * Translates vendor string to enum type
 */
{
   enum VEND iret;
   if (strstr(vendor, "GenuineIntel") != NULL)
      iret = Intel;
   else if (strstr(vendor, "AuthenticAMD") != NULL)
      iret = AMD;
   else if (strstr(vendor, "GenuineTMx86") != NULL)
      iret = TM;
   else
      iret = VERR;
   return(iret);
}

/*
 * Specific chip (family, but disambiguated using vendor string
 */
enum CHIP {CERR, Pentium, IntP6, Pentium4, Itanium, K7, Hammer, HammerB,
           Crusoe, Efficeon};

enum CHIP Family2Chip(char *vendor, enum FAM family)
/*
 * Disambiguates family based on vendor string
 */
{
   enum CHIP iret=CERR;
   enum VEND ivend;

/*
 * Figure out the vendor
 */
   ivend = str2vend(vendor);
   if (ivend == VERR)
      return(CERR);

   switch(family)
   {
   case   i486:        /* 486 & AMD 5x86 and Cyrix 5x86; unsupported */
      break;
   case   P5:          /* Original Pentium and AMD K5 & K6 */
      if (ivend == Intel)
         iret = Pentium;
      break;
   case   P6:          /* Intel PIII, Core and AMD K7 (orig athlon) */
      if (ivend == Intel)
         iret = IntP6;
      else if (ivend == AMD)
         iret = K7;
      else if (ivend == TM)
         iret = Crusoe;
      break;
   case   P7:          /* Intel P4, AMD hammer, Efficeon */
      if (ivend == Intel)
         iret = Pentium4;
      else if (ivend == AMD)
         iret = Hammer;
      else if (ivend == TM)
         iret = Efficeon;
      break;
   case P8B:
      if (ivend == AMD)
         iret = HammerB;
      break;
   case   ITAN:        /* Intel Itanium */
      iret = Itanium;
      break;
   default:
      iret = CERR;
   }
   return(iret);
}

enum MACHTYPE Chip2Mach(enum CHIP chip, int model, int x8664)
/*
 * translates chip and cpuid's model to config's machine enum
 */
{
   enum MACHTYPE iret=MACHOther;

   switch(chip)
   {
   case Pentium:
      switch(model)
      {
      case 1:
         iret = IntP5;
         break;
      case 4:
      case 8:
         iret = IntP5MMX;
         break;
      default:
         iret = MACHOther;
      }
      break;
   case IntP6:  /* includes PPRO, PII, PIII, Core and Pentium-M */
      switch(model)
      {
      case 0:
      case 1:
         iret = IntPPRO;
         break;
      case 3:
      case 5:
      case 6:
         iret = IntPII;
         break;
      case 7:
      case 8:
      case 10:
      case 11:
         iret = IntPIII;
         break;
      case  9:
      case 13:
         iret = IntPM;
         break;
      case 14:
         iret = IntCoreDuo;
         break;
      case 15:
      case 23:
      case 29:
         iret = IntCore2;
         break;
      case 26:
         iret = IntCorei7;
         break;
      default:
         iret = MACHOther;
      }
      break;
   case Pentium4:
      switch(model)
      {
      case 0:
      case 1:
      case 2:
         iret = IntP4;
         break;
      case 3:
      case 4: ; case 6:
         iret = IntP4E;
         break;
      default:
         iret = MACHOther;
      }
      break;
   case Itanium:
      switch(model)
      {
      case 7:
         iret = IA64Itan;
         break;
      case 0x1F:
         iret = IA64Itan2;
         break;
      default:
         iret = MACHOther;
      }
      break;
   case K7:
      switch(model)
      {
      case 4:
      case 6:
      case 8:
      case 10:
         iret = AmdAthlon;
         break;
      default:
         iret = MACHOther;
      }
      break;
   case Hammer:
      iret = AmdHammer;
      break;
   case HammerB:
      iret = Amd64K10h;
      break;
   case Efficeon:
      iret = TMEff;
      break;
   case Crusoe:  /* unsupported */
   default:
      iret = MACHOther;
   }
   return(iret);
}

void PrintUsage(char *name, int i)
{
   fprintf(stderr, "USAGE: %s -v (verb) -b (@ bits) -a (arch) -n (ncpu) -c <ncache> -C <lvl> (cache size) -m (Mhz) -t (cpu throttling)\n", name);
   exit(i);
}

int GetFlags(int nargs, char **args, int *CacheLevel)
{
   int i, flag = 0;

   *CacheLevel = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0], i);
      switch(args[i][1])
      {
      case 'n':
         flag |= Pncpu;
         break;
      case 'c':
         flag |= Pncache;
         break;
      case 'C':
         if (++i > nargs)
            PrintUsage(args[0], i);
         *CacheLevel = atoi(args[i]);
         break;
      case 'v':
         flag |= Pverb;
         break;
      case 'm':
         flag |= PMhz;
         break;
      case 'a':
         flag |= Parch;
         break;
      case 'b':
         flag |= P64;
         break;
      case 't':
         flag |= Pthrottle;
         break;
      default:
         PrintUsage(args[0], i);
      }
   }
   if (!flag)
     flag = Parch | P64;
   return(flag);
}

main(int nargs, char **args)
{
   int ierr, x86_64, flags, CacheLevel;
   unsigned family, model;
   char *cpu="UNKNOWN", vendor[13];
   enum FAM fam;
   enum CHIP chip;
   enum MACHTYPE mach;

   flags = GetFlags(nargs, args, &CacheLevel);
   cpu = NULL;
   vendor[0] = '\0';
   ierr = ProbeArch(vendor, &family, &model, &x86_64);
/*
 * If ProbeArch worked, translate vendor+family+model to ATLAS config-name
 */
   if (!ierr)
   {
       fam = GetFamily(family);
       if (fam)
       {
          chip = Family2Chip(vendor, fam);
          if (chip)
          {
             mach = Chip2Mach(chip, model, x86_64);
             if (!mach) ierr = 300;
          }
          else ierr = 200;
       }
       else ierr = 100;
   }
   if (ierr)
   {
      fprintf(stderr, "ERROR: enum fam=%d, chip=%d, mach=%d\n",
              fam, chip, mach);
      printf("ERROR %d: vendor='%s', family=%d, model=%d, x86_64=%d\n",
             ierr, vendor, family, model, x86_64);
   }
   else
   {
/*
 *    If verbatim set, print strings as well as enums
 */
      if (flags & Parch)
      {
         if (flags & Pverb)
            printf("cpu: %s\n", machnam[mach]);
         printf("MACHTYPE=%d\n", mach);
      }
      if (flags & P64)
         printf("PTR BITS=%d\n", x86_64 ? 64 : 32);
/*
 *    Not sure how to detect this.  cpuid has some features that might work,
 *    will need to experiment later
 */
      if (flags & Pthrottle)
         printf("CPU THROTTLE=0\n");
/*
 *    These guys can't be supported by cpuid, AFAIK
 */
      if ((flags & PMhz) || (flags & Pncpu))
         printf("Mhz/ncpu=0\n");
/*
 *    Cache info could be returned, but I'm lazy, so don't
 */
      if ((flags & Pncache) || (flags & PCacheSize))
         printf("ncache/CacheSize=0\n");
      if ((flags & (~Pverb)) == 0)
         printf("family=%d, model=%d, cpu='%s', Ptr bits=%d, arch#=%d\n",
                family, model, machnam[mach], x86_64?64:32, mach);
   }
   exit(ierr);
}
