#include "atlconf.h"

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

enum MACHTYPE ProbeArch()
{
   enum ARCHFAM fam;
   enum MACHTYPE mach=MACHOther;
   int ierr, i;
   char res[1024];

   fam = ProbeArchFam(NULL);
   switch(fam)
   {
   case AFPPC:
      if ( !CmndOneLine(NULL, "cat /proc/cpuinfo | fgrep cpu", res) )
      {
#if 0
         if (strstr(res, "604e")) mach = PPC604e;
         else if (strstr(res, "604")) mach = PPC604;
         else
#endif
         if (strstr(res, "G4")) mach = PPCG4;
         else if (strstr(res, "7400")) mach = PPCG4;
         else if (strstr(res, "7410")) mach = PPCG4;
         else if (strstr(res, "7447")) mach = PPCG4;
         else if (strstr(res, "7455")) mach = PPCG4;
         else if (strstr(res, "PPC970FX")) mach = PPCG5;
         else if (strstr(res, "POWER5")) mach = IbmPwr5;
         else if (strstr(res, "POWER4")) mach = IbmPwr4;
      }
      break;
   case AFMIPS:
      res[0] = '\0';
      ierr = CmndOneLine(NULL, "fgrep 'cpu model' /proc/cpuinfo", res);
      if (!ierr && res[0] != '\0')
      {
         if (strstr(res, "ICE9"))
            mach = MIPSICE9;
/*
 *       I have no access to what cpuinfo on Linux does for this procs, so this
 *       is a WAG as to what it would say
 */
         else if (strstr(res, "R10000") || strstr(res, "R12000") ||
                  strstr(res, "R12000") || strstr(res, "R14000"))
            mach = MIPSR1xK;
      }
      break;
   case AFIA64:
      res[0] = '\0';
      ierr = CmndOneLine(NULL, "fgrep 'IA-64' /proc/cpuinfo", res);
      if (ierr || res[0] == '\0')
         ierr = CmndOneLine(NULL, "fgrep \"model name\" /proc/cpuinfo", res);
      if (!ierr && res[0] != '\0')
      {
         if (strstr(res, "IA-64") || strstr(res, "McKinley"))
            mach = IA64Itan2;
         else if (strstr(res, "Itanium")) mach = IA64Itan;
      }
      break;
   case AFX86:
      res[0] = '\0';
      ierr = CmndOneLine(NULL, "fgrep 'model name' /proc/cpuinfo", res);
      if (ierr || res[0] == '\0')
         ierr = CmndOneLine(NULL, "fgrep model /proc/cpuinfo", res);
      if (!ierr && res[0] != '\0')
      {
         if (strstr(res, "Pentium"))
         { /* Pentium of some flavor */
            if (strstr(res, " III ")) mach = IntPIII;
            else if (strstr(res, " II ")) mach = IntPII;
            else if (strstr(res, "Pro")) mach = IntPPRO;
            else if (strstr(res, "MMX")) mach = IntP5MMX;
            else if (strstr(res, " 4 "))
            {
               ierr = CmndOneLine(NULL,
                      "fgrep 'model' /proc/cpuinfo | fgrep -v 'name'", res);
               if (!ierr)
               {
                  i = GetLastInt(res);
                  if (i < 3) mach = IntP4;
                  else if (i == 3) mach = IntP4E;
               }
            }
         }
         else if (strstr(res, "Efficeon")) mach = TMEff;
         else if (strstr(res, "Athlon HX")) mach = AmdHammer;
         else if (strstr(res, "Opteron") || strstr(res, "Hammer") ||
                  strstr(res, "Athlon(tm) 64"))
            mach = AmdHammer;
         else if (strstr(res, "Athlon")) mach = AmdAthlon;
         else if (strstr(res, "AMD-K7")) mach = AmdAthlon;
      }
      break;
/*
 *    Add these back if we get machine access and can test
 */
   case AFSPARC:  /* don't know here anymore */
      #if 0
      if ( !CmndOneLine(NULL, "fgrep cpu /proc/cpuinfo", res) )
      {
         if (strstr(res, "UltraSparc II")) mach = SunUS2;
         else if (strstr(res, "UltraSparc I")) mach = SunUS1;
         else if (strstr(res, "UltraSparc")) mach = SunUSX;
      }
      #endif
      break;
   case AFALPHA:
      #if 0
      res[0] = '\0';
      ierr = CmndOneLine(NULL, "fgrep 'model name' /proc/cpuinfo", res);
      if (ierr || res[0] == '\0')
         ierr = CmndOneLine(NULL, "fgrep model /proc/cpuinfo", res);
      if (!ierr && res[0] != '\0')
      {
         if (strstr(res, "EV5")) mach = Dec21164;
         else if (strstr(res, "EV4")) mach = Dec21064;
         else if (strstr(res, "EV6")) mach = Dec21264;
      }
      #endif
      break;
   default:
#if 0
      if (!CmndOneLine(NULL, "fgrep 'cpu family' /proc/cpuinfo", res))
         if (strstr(res, "PA-RISC 2.0")) mach = HPPA20;
#else
     ;
#endif
   }
   return(mach);
}

int ProbeNCPU()
{
   int ncpu = 0;
   char *reslns, res[1024];

   #if 0
   if (mach == Dec21264 || mach == Dec21164 || mach == Dec21064)
   {
      if ( !CmndOneLine(NULL, "fgrep 'cpus detected' /proc/cpuinfo", res) )
         ncpu = GetLastInt(res);
   }
   #endif
   if (!ncpu)
   {
      reslns = CmndResults(NULL, "grep '^processor' /proc/cpuinfo");
      if (reslns) ncpu = fNumLines(reslns);
   }
   return(ncpu);
}

int ProbePointerBits(int *sure)
{
   int iret = 32;
   char *uname;
   char cmnd[1024], res[1024];

/*
 * Note this is a weak probe, archinfo_x86 much better . . .
 */
   *sure = 0;
   uname = FindUname(NULL);
   sprintf(cmnd, "%s -a", uname);
/*
 * This probe should be running on backend; if its ptr length is 8, we've
 * definitely got a 64 bit machine
 * NOTE: getting 4 could be a result of compiler flags on a 64-bit arch,
 * so reverse is not dispositive
 */
   if (sizeof(void*) == 8)
   {
      iret = 64;
      *sure = 1;
   }
   else if (!CmndOneLine(NULL, cmnd, res))
   {
/*
 *    If uname is a known 64-bit platform, we're sure we've got OS support
 *    for 64bits (may not have compiler support, but that's not our fault)
 */
      if (strstr(res, "x86_64") || strstr(res, "ppc64") || strstr(res, "ia64"))
      {
         iret = 64;
         *sure = 1;
      }
   }
   return(iret);
}

int ProbeMhz()
{
   int mhz=0;
   char res[1024];
   if (!CmndOneLine(NULL, "fgrep 'cpu MHz' /proc/cpuinfo", res))
      mhz = GetFirstInt(res);
   return(mhz);
}

int ProbeThrottle()
/*
 * RETURNS: 1 if cpu throttling is detected, 0 otherwise
 */
{
   int iret=0;
   int imax=0, imin=0, icur=0;
   char res[1024];

/*
 * If cpufreq directory doesn't exist, guess no throttling.  If
 * cpufreq exists, and cur Mhz < max, throttling is enabled,
 * throttling also enabled if governer is not "performance", and min freq
 * is less than max
 */
   if (!CmndOneLine(NULL,
       "cat /sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq", res) )
   {
      imax = GetFirstInt(res);
      assert(!CmndOneLine(NULL,
             "cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_min_freq", res));
      imin = GetFirstInt(res);
      assert(!CmndOneLine(NULL,
             "cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq", res));
      icur = GetFirstInt(res);
      assert(!CmndOneLine(NULL,
             "cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor", res));
      if (icur < imax)
         iret = 1;
      else if (!strstr(res, "performance") && imin < imax)
         iret = 1;
   }
   return(iret);
}

main(int nargs, char **args)
{
   int flags, CacheLevel, ncpu, mhz, bits, sure;
   enum MACHTYPE arch=MACHOther;

   flags = GetFlags(nargs, args, &CacheLevel);
   if (flags & Parch)
   {
      arch = ProbeArch();
      if (flags & Pverb)
         printf("Architecture detected as %s.\n", machnam[arch]);
      printf("MACHTYPE=%d\n", arch);
   }
   if (flags & Pncpu)
      printf("NCPU=%d\n", ProbeNCPU());
   if (flags & PMhz)
      printf("CPU MHZ=%d\n", ProbeMhz());
   if (flags & Pthrottle)
      printf("CPU THROTTLE=%d\n", ProbeThrottle());
   if (flags & P64)
   {
      bits = ProbePointerBits(&sure);
      printf("PTR BITS=%d, SURE=%d\n", bits, sure);
   }

/*
 * Here for future, presently unsupported
 */
   if (flags & Pncache)
      printf("NCACHES=0\n");
   if (flags & PCacheSize)
      printf("%d Cache size (kb) = 0\n", CacheLevel);
   exit(0);
}
