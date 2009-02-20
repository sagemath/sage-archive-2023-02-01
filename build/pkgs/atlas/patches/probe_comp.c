#include "atlconf.h"

COMPNODE *GetCompNode(void)
{
   COMPNODE *p;
   p = calloc(1, sizeof(COMPNODE));
   assert(p);
   return(p);
}
COMPNODE *KillCompNode(COMPNODE *die)
{
   COMPNODE *p=NULL;
   if (die)
   {
      p = die->next;
      if (die->comp)
         free(die->comp);
      if (die->flags)
         free(die->flags);
      free(die);
   }
   return(p);
}

void KillAllCompNodes(COMPNODE *kill)
{
   while(kill)
      kill = KillCompNode(kill);
}

COMPNODE *CloneCompNode(COMPNODE *orig)
{
   COMPNODE *new;

   new = GetCompNode();
/*
 * Copy everything but strings wt memcopy
 */
   memcpy(new, orig, sizeof(COMPNODE));
   if (orig->comp)
   {
      new->comp = malloc(sizeof(char)*(strlen(orig->comp)+1));
      assert(new->comp);
      strcpy(new->comp, orig->comp);
   }
   if (orig->flags)
   {
      new->flags = malloc(sizeof(char)*(strlen(orig->flags)+1));
      assert(new->flags);
      strcpy(new->flags, orig->flags);
   }
   return(new);
}

void PrintCompNodes(FILE *fpout, COMPNODE *q, char *id)
{
   int i;
   COMPNODE *p;

   fprintf(fpout, "\nCompiler nodes: %s\n", id);
   if (!q)
      fprintf(fpout, "***NONE***\n");
   for (i=0, p=q; p; p=p->next, i++)
   {
      fprintf(fpout,
         "%3d. priority=%d, comps,OS,arch[0]=(%d,%d,%d), comp='%s'\n",
              i, p->priority, p->comps[0], p->OS[0], p->arch[0], p->comp);
      fprintf(fpout, "     '%s'\n", p->flags);
   }
   fprintf(fpout, "\n");
}

COMPNODE *SortCompsByPriority(COMPNODE *q)
/*
 * Builds new queue out of q, by first adding largest priority to newq, etc.
 */
{
   COMPNODE *newq=NULL, *p, *prev, *maxprev;
   int maxpri;

   while(q)  /* BFI N^2 sort */
   {
/*
 *    Find remaining element wt largest priority
 */
      maxprev = NULL;
      prev = q;
      maxpri = q->priority;
      for (p=q->next; p; p = p->next)
      {
         if (p->priority > maxpri)
         {
            maxpri = p->priority;
            maxprev = prev;
         }
         prev = p;
      }
/*
 *    Take max node off of q, add to newq
 */
      if (maxprev)  /* max elt wasn't stop of stack (q) */
      {
         p = maxprev->next->next;
         maxprev->next->next = newq;
         newq = maxprev->next;
         maxprev->next = p;
      }
      else
      {
         p = q->next;
         q->next = newq;
         newq = q;
         q = p;
      }
   }
   return(newq);
}

void DivideCompsByComp(COMPNODE *q, COMPNODE **comps)
/*
 * Builds individual queues for each compiler camp, and then kills the original
 * queue.  Note that original q is sorted smallest-to-largest, and since we
 * add in order to a stack, we wind up with largest-to-smallest, as we want.
 * Note comps is really COMPNODE *comps[NCOMP], and can be indexed by ICC_, etc.
 */
{
   int i;
   COMPNODE *p, *new;

   for (i=0; i < NCOMP; i++)
      comps[i] = NULL;
   for (p=q; p; p = p->next)
   {
      for (i=0; i < NCOMP; i++)
      {
         if (IsBitSetInField(p->comps, i))
         {
            new = CloneCompNode(p);
            new->next = comps[i];
            comps[i] = new;
         }
      }
   }
   KillAllCompNodes(q);
}

static int OSMatches(enum OSTYPE OS, int *bits)
/*
 * RETURNS: 0 if no OS of arch in bitfield, nonzero otherwise
 */
{
   if (IsBitSetInField(bits, 0))   /* If 0 bit set, matches all OS */
      return(1);
   if (IsBitSetInField(bits, OS))  /* If OS bit set, matches this OS */
      return(2);
   return(0);
}

static int ArchMatches(enum MACHTYPE arch, int *bits)
/*
 * RETURNS: 0 if no match of arch in bitfield, nonzero otherwise
 */
{
   if (IsBitSetInField(bits, 0))     /* If 0 bit set, matches all archs */
      return(1);
   if (IsBitSetInField(bits, arch)) /* If arch bit set, matches this arch */
      return(2);
   return(0);
}

COMPNODE *KillBadArchOS(enum OSTYPE OS, enum MACHTYPE arch, COMPNODE *q)
/*
 * Deletes all non-matching OS/arch from queue; Note any node wt these
 * quantities set to 0 is a wildcard, and so stays
 */
{
   COMPNODE *prev, *next, *p;
   if (!OS && !arch)
      return(q);
/*
 * Delete all beginning nodes until we find one for this arch
 */
   while(q && (!ArchMatches(arch, q->arch) || !OSMatches(OS, q->OS)))
      q = KillCompNode(q);
/*
 * With good top of stack, delete trailing nodes that don't match
 */
   if (q)
   {
      prev = q;
      for (p=q->next; p; p = next)
      {
         next = p->next;
         if (!ArchMatches(arch, p->arch) || !OSMatches(OS, p->OS))
            prev->next = KillCompNode(p);
         else
            prev = p;
      }
   }
   return(q);
}

int OSNameToInt(char *name)
{
   int i;
   for (i=1; i < NOS; i++)
   {
      if (!strcmp(name, osnam[i]))
         return(i);
   }
   return(0);
}

int MachNameToInt(char *name)
{
   int i;
   for (i=1; i < NMACH; i++)
   {
      if (!strcmp(name, machnam[i]))
         return(i);
   }
   return(0);
}
void NamesToBitField(int MACH, char *str, int *bits)
/*
 * Takes a list of machine (MACH=1) or OS (MACH=0) names and translates them
 * to their enumerated type numbers, and sets the appropriate bit in the
 * bits field
 */
{
   char name[128];
   int i=0;
   while(*str)
   {
      if (*str == ',' || *str == ' ' || *str == '\t' || *str == '\n' ||
          *str == '\0')
      {  /* finished a name */
         name[i] = '\0';
         if (!strcmp(name, "all") || !strcmp(name, "ALL"))
            i = 0;
         else
         {
            if (MACH)
               i = MachNameToInt(name);
            else
               i = OSNameToInt(name);
            if (!i)
            {
               fprintf(stderr, "Nonsensical %s name in list: %s\n",
                       MACH ? "machine" : "OS", str);
               exit(1);
            }
         }
         SetBitInField(bits, i);
         if (*str != ',')  /* anything but ',' ends list */
            return;
         str++;
         i = 0;
      }
      else
         name[i++] = *str++;
   }
}

void NumListToBitfield(char *str, int *bits)
/*
 * Takes a list of number like : '5,3,2,8,0' and turns on those bits
 */
{
   char num[16], *sp;
   int i, j;

   while (*str)
   {
      for (sp=num,i=0; i < 15; i++)
      {
         if (!isdigit(str[i])) break;
         *sp++ = str[i];
      }
      assert(i != 15);
      *sp = '\0';
      j = atoi(num);
      SetBitInField(bits, j);
      if (str[i] != ',') break;
      str += i+1;
   }
}

COMPNODE *ParseNewCompLine(char *ln)
{
   COMPNODE *p;
   char *sp;
   p = GetCompNode();

   sp = strstr(ln, "MACH=");
   assert(sp);
   sp += 5;
/*   NumListToBitfield(sp, p->arch); */
   NamesToBitField(1, sp, p->arch);


   sp = strstr(ln, "OS=");
   assert(sp);
   sp += 3;
/*   NumListToBitfield(sp, p->OS); */
   NamesToBitField(0, sp, p->OS);

   sp = strstr(ln, "LVL=");
   assert(sp);
   sp += 4;
   p->priority = atoi(sp);
/*
 * Parse 'COMPS=[icc,smc,dmc,skc,dkc,xcc,f77]', at least one comp must exist
 */
   sp = strstr(ln, "COMPS=");
   assert(sp);
   sp += 6;
   while (*sp)
   {
      if (sp[0] == 'i' && sp[1] == 'c' && sp[2] == 'c')
         SetBitInField(p->comps, ICC_);
      else if (sp[0] == 's' && sp[1] == 'm' && sp[2] == 'c')
         SetBitInField(p->comps, SMC_);
      else if (sp[0] == 'd' && sp[1] == 'm' && sp[2] == 'c')
         SetBitInField(p->comps, DMC_);
      else if (sp[0] == 's' && sp[1] == 'k' && sp[2] == 'c')
         SetBitInField(p->comps, SKC_);
      else if (sp[0] == 'd' && sp[1] == 'k' && sp[2] == 'c')
         SetBitInField(p->comps, DKC_);
      else if (sp[0] == 'x' && sp[1] == 'c' && sp[2] == 'c')
         SetBitInField(p->comps, XCC_);
      else if (sp[0] == 'f' && sp[1] == '7' && sp[2] == '7')
         SetBitInField(p->comps, F77_);
      else
      {
         fprintf(stderr, "WTF(%d of %s): '%s'??\n", __LINE__, __FILE__, sp);
         exit(-1);
      }
      if (sp[3] != ',') break;
      sp += 4;
   }
   return(p);
}

char *CopySingleQuoteString(char *str, char *out)
/*
 * Finds the leading ' in str, and copies quoted text to out until closing '
 * is found or end of string.
 * RETURNS: pointer in str at closing ', NULL if closing ' not found
 */
{
/*
 * Skip text before staring quote, return if quote not found
 */
   *out = '\0';
   while (*str != '\'' && *str) str++;
   if (*str == '\0')
      return(NULL);
/*
 * Copy quoted text, return pointer to end quote if it exists
 */
   str++;
   while (*str != '\'' && *str) *out++ = *str++;
   *out++ = '\0';
   if (*str == '\'')
      return(str);
   return(NULL);
}

COMPNODE *ParseCompLine(char *ln)
{
   static int NewComp=1;
   static COMPNODE *p=NULL;
   int i;
   char *sp;
   char ln2[2048];

   if (NewComp)
      p = ParseNewCompLine(ln);
/*
 * Should be line of form "'compiler' 'flags'"
 */
   else
   {
      sp = CopySingleQuoteString(ln, ln2);
      assert(sp);
      i = strlen(ln2) + 1;
      p->comp = malloc(sizeof(char)*i);
      assert(p->comp);
      strcpy(p->comp, ln2);

      sp = CopySingleQuoteString(sp+1, ln2);
      i = strlen(ln2) + 1;
      p->flags = malloc(sizeof(char)*i);
      assert(p->flags);
      strcpy(p->flags, ln2);
   }
   NewComp = !NewComp;
   return(p);
}

COMPNODE *ReadComps(char *file)
/*
 * Reads in a file describing the compilers ATLAS knows about, and returns
 * a queue of them for later manipulation.
 */
{
   char ln[2048];
   FILE *fpin;
   COMPNODE *compq=NULL, *p;

   fpin = fopen(file, "r");
   while (fgets(ln, 2048, fpin))
   {
      if (ln[0] != '#')
      {
         KillUselessSpace(ln);
         if (ln[0] != '#' && ln[0] != '\0')
         {
            p = ParseCompLine(ln);
            if (p != compq)
            {
               p->next = compq;
               compq = p;
            }
         }
      }
   }
   fclose(fpin);
   return(compq);
}

COMPNODE **GetDefaultComps(enum OSTYPE OS, enum MACHTYPE arch, int verb)
/*
 * This routine reads the file atlcomp.txt, and returns them sorted by
 * order of priority for each compiler ATLAS needs.  This list can then
 * be matched with the user's input to give final compiler and flags.
 */
{
   COMPNODE *q, *p, **comps;

   comps = malloc(sizeof(COMPNODE)*NCOMP);
   q = ReadComps("atlcomp.txt");    /* get all compiler lines */
   if (verb > 1)
      PrintCompNodes(stderr, q, "Fresh Read");
   q = KillBadArchOS(OS, arch, q);  /* discard comps for other platforms */
   if (verb > 1)
      PrintCompNodes(stderr, q, "Targeted");
   q = SortCompsByPriority(q);      /* q is smallest, bottom is largest */
   if (verb > 1)
      PrintCompNodes(stderr, q, "Sorted");
   DivideCompsByComp(q, comps);     /* split into individual queues */
   return(comps);
}

int CompTest(int verb, char *targ, int icomp, char *comp, char *flag)
/*
 * Tries to build simple program and run it.
 * RETURNS: result of system call (0: success, else error code)
 */
{
   char ln[1024], res[1024];
   int i, iret;

   if (icomp == ICC_)
      targ = NULL;
   if (icomp == F77_)
      i = sprintf(ln, "make IRunF77Comp F77='%s' F77FLAGS='%s' ", comp, flag);
   else
      i = sprintf(ln, "make IRunCComp CC='%s' CCFLAGS='%s' ", comp, flag);
   if (targ)
      i += sprintf(ln+i, "atlrun=atlas_runX targ=%s ", targ);
   i += sprintf(ln+i, "| fgrep SUCCESS");
   iret = CmndOneLine(NULL, ln, res);
   if (verb > 1)
      fprintf(stderr, "cmnd=%s\n", ln);
   if (!iret)
      iret = !strstr(res, "SUCCESS");
   if (verb)
      fprintf(stderr, "   %s %s : %s!\n", comp, flag,
              iret ? "FAILURE":"SUCCESS");
   return(iret);
}

void CompError(int icomp)
/*
 * Prints out informative error message when we die because a compiler doesn't
 * work
 */
{
   fprintf(stderr, "\n\nUnable to find usable compiler for %s; aborting",
           COMPNAME[icomp]);
   fprintf(stderr, "Make sure compilers are in your path, and specify good compilers to configure\n");
   fprintf(stderr,
           "(see INSTALL.txt or 'configure --help' for details)");
   exit(icomp+1);
}

char *GetPtrbitsFlag(enum OSTYPE OS, enum MACHTYPE arch, int ptrbits,
                     char *comp)
/*
 * RETURNS: string forcing setting of ptrbits for gcc
 */
{
   char *sp = "";
   int i, j, k;

   if (MachIsIA64(arch))
      return(sp);
   if (MachIsMIPS(arch))
      return((ptrbits == 64) ? "-mabi=64" : "-mabi=n32");
   if (!CompIsGcc(comp))
   {
/*
 *    Add correct 64/32 bit flags for Sun workshop compilers
 */
      if (MachIsUS(arch) && CompIsSunWorkshop(comp))
      {
         if (ptrbits == 64)
            sp = (arch == SunUSI || arch == SunUSII) ?
                 "-xarch=v9" : "-xarch=v9b";
         else
            sp = (arch == SunUSI || arch == SunUSII) ?
                 "-xarch=v8plusa" : "-xarch=v8plusb";
      }
      else if (CompIsIBMXL(comp))  /* IBM xl compilers */
         sp = (ptrbits == 64) ? "-q64" : "-q32";
      return(sp);
   }
   GetGccVers(comp, &k, &j, &k, &k);
   if ( !(j >= 3 && (OS != OSOSX || j > 3 || !CompIsAppleGcc(comp))) )
      return(sp);
   else if (OS == OSAIX)
      sp = (ptrbits == 64) ? "-maix64" : "-maix32";
   else if (arch == IA64Itan2)
      printf("Itanium2 - not setting -m64"); // -m64 is not supported on RHEL 5/Itanium
   else if (ptrbits == 64)
     sp = "-m64";
   else if (ptrbits == 32)
     sp = "-m32";
   return(sp);
}
char *GetStandardCompName(char *comp)
{
   int i, j, k;
   char *ucomp;
/*
 * Recognize gnu compiler regardless of name string (eg. ev6-gcc-3.2)
 */
   if (CompIsGcc(comp))
   {
      GetGccVers(comp, &k, &j, &k, &k);
      if (j < 4)
      {
         if (i == F77_)
            ucomp = "g77";
         else
            ucomp = "gcc";
      }
      else if (i == F77_)
         ucomp = "gfortran";
      else
         ucomp = "gcc";
   }
   else
      ucomp = NameWithoutPath(comp);
   return(ucomp);
}
char *GetWinComp(enum OSTYPE OS, char *comp)
{
   char ln[1024];
   char *ucomp;
   if (!OSIsWin(OS))
      return(NULL);
   ucomp = GetStandardCompName(comp);
   if (!strcmp(ucomp, "icc") || !strcmp(ucomp, "icl"))
      ucomp = "ATLwin_icc";
   else if (!strcmp(ucomp, "ifort") || !strcmp(ucomp, "ivf"))
      ucomp = "ATLwin_ifort";
   else if (!strcmp(ucomp, "mvc") || !strcmp(ucomp, "cl"))
      ucomp = "ATLwin_cl";
   else /* not a recognized windows compiler that needs wrapping, done */
      return(NULL);
   sprintf(ln, "make wind=/usr/local/bin /usr/local/bin/%s.exe \n", ucomp);
   if (system(ln))
   {
      fprintf(stderr, "Unable to to build %s, quitting\n", ucomp);
      fprintf(stderr, "cmnd='%s'\n", ln);
      exit(-1);
   }
   return(ucomp);
}

void GetComps(enum OSTYPE OS, enum MACHTYPE arch, int verb, char *targ,
              int ptrbits, char **usrcomps, int nof77)
/*
 * This routine gives config a list of compilers to use.  The first NCOMP
 * entries in usrcomps indicate a user override of the default compiler,
 * and the next NCOMP entries indicate user override of flags.  The next
 * NCOMP entries indicate that those flags should be appended to prior flags.
 * A NULL in any entry says the user is happy to use the defaults (or no
 * appending).  Chosen compilers and flags are returned in usrcomps array.
 */
{
   COMPNODE **comps, *p;
   char *ucomp, *dcomp, *flg, *sp, *sp2;
   char cmnd[256], res[1024];
   int i, j, k, h;

   comps = GetDefaultComps(OS, arch, verb);
   if (verb > 1)
      fprintf(stdout, "Finding good compilers:\n");
   for (i=0; i < NCOMP; i++)
   {
      if (nof77 && i == F77_) continue;
/*
 *    If the user has not specified the compiler, look through all available
 *    compilers with one that works (with user flags, if specified)
 */
      if (!usrcomps[i])
      {
         for (p=comps[i]; p; p = p->next)
         {
            flg = NewStringCopy(usrcomps[NCOMP+i]?usrcomps[NCOMP+i]:p->flags);
            if (usrcomps[NCOMP*2+i])
               flg = NewAppendedString(flg, usrcomps[NCOMP*2+i]);
            sp = GetWinComp(OS, p->comp);
            if (!sp) sp = p->comp;
            if (ptrbits)
               flg = NewAppendedString(flg,
                                       GetPtrbitsFlag(OS, arch, ptrbits, sp));
            if (!CompTest(verb, targ, i, sp, flg))
               break;
            free(flg);
         }
         if (!p)
            CompError(i);
         else
            free(flg);
         usrcomps[i] = p->comp;
         p->comp = NULL;                /* so it isn't deleted by Kill */
         if (!usrcomps[NCOMP+i])
         {
            usrcomps[NCOMP+i] = p->flags;
            p->flags = NULL;            /* so it isn't deleted by Kill */
         }
      }
/*
 *    If user specified comp w/o flags, get default flags or error
 */
      else if (!usrcomps[NCOMP+i])
      {
         p = comps[i];
         ucomp = NameWithoutPath(usrcomps[i]);
/*
 *       Recognize gnu compiler regardless of name string (eg. ev6-gcc-3.2)
 */
         if (CompIsGcc(usrcomps[i]))
         {
            GetGccVers(usrcomps[i], &k, &j, &k, &k);
            if (j < 4)
            {
               if (i == F77_)
                  ucomp = "g77";
               else
                  ucomp = "gcc";
            }
            else if (i == F77_)
               ucomp = "gfortran";
            else
               ucomp = "gcc";
         }
         for (p=comps[i]; p; p = p->next)
         {
            dcomp = NameWithoutPath(p->comp);
            if (!strcmp(p->comp, ucomp))
               break;
            free(dcomp);
         }
         if (!p)
         {
            fprintf(stderr,
               "UNKNOWN COMPILER '%s' for %s: you must also supply flags!\n",
                    usrcomps[i], COMPNAME[i]);
            exit(i+1);
         }
         usrcomps[NCOMP+i] = p->flags;
         p->flags = NULL;
      } /* If user specifed both flags and compiler, accept them */
/*
 *    On windows, build compiler wrapper for MSVC++ or Intel compilers
 */
      sp = GetWinComp(OS, usrcomps[i]);
      if (sp)
      {
         free(usrcomps[i]);
         usrcomps[i] = NewStringCopy(sp);
      }
/*
 *    Test selected compiler and flags, and die if they don't work
 */
      flg = NewStringCopy(usrcomps[NCOMP+i]?usrcomps[NCOMP+i]:p->flags);
      if (usrcomps[NCOMP*2+i])
         flg = NewAppendedString(flg, usrcomps[NCOMP*2+i]);
      if (ptrbits)
         flg = NewAppendedString(flg,
                                 GetPtrbitsFlag(OS, arch, ptrbits,usrcomps[i]));
      if (CompTest(verb, targ, i, usrcomps[i], flg))
         CompError(i);
      free(flg);
   } /* end of loop over compilers */
/*
 * modify base flags by appending user flags
 */
   for (i=2*NCOMP; i < 3*NCOMP; i++)
   {
      if (usrcomps[i])  /* user has appended flags for compiler i-2*NCOMP */
         usrcomps[i-NCOMP] = NewAppendedString(usrcomps[i-NCOMP], usrcomps[i]);
   }
/*
 * If nof77, set fortran compiler & flags to ICC to avoid linking problems
 */
   if (nof77)
   {
      usrcomps[F77_] = NewStringCopy(usrcomps[ICC_]);
      usrcomps[F77_+NCOMP] = NewStringCopy(usrcomps[ICC_+NCOMP]);
   }
/*
 * If ptrbits is set to manual override, add -m32/64 to gnu compilers
 * but not on Itaniums or Apple's munged gcc 3 compiler!
 */
   if (ptrbits && arch != IA64Itan && arch != IA64Itan2)
   {
      for (i=0; i < NCOMP; i++)
      {
         sp = GetPtrbitsFlag(OS, arch, ptrbits, usrcomps[i]);
         usrcomps[i+NCOMP] = NewAppendedString(usrcomps[i+NCOMP], sp);
      }
   }
/*
 * Need to add target & bitwidth args for MIPSpro compilers on IRIX
 */
   if (OS == OSIRIX)
   {
      sprintf(cmnd, "%s -m", FindUname(targ));
      assert(!CmndOneLine(NULL, cmnd, res));
      sp = strstr(res, "IP");
      for (i=2; isdigit(sp[i]); i++);
      sp[i] = '\0';
      sprintf(cmnd, "-TARG:platform=%s", sp);
      if (ptrbits == 64 || !ptrbits)
         strcat(cmnd, " -64");
      else
         strcat(cmnd, " -32");
      for (i=0; i < NCOMP; i++)
      {
         if (CompIsMIPSpro(usrcomps[i]))
         {
            usrcomps[i+NCOMP] = NewAppendedString(usrcomps[i+NCOMP], cmnd);
         }
      }
      free(sp);
   }
}  /* end of routine GetComps */

void TestComps(enum OSTYPE OS, enum MACHTYPE arch, int verb, char *targ,
               char *targarg, char **comps, enum F2CNAME *f2cnam,
               enum F2CINT *f2cint, enum F2CSTRING *f2cstr, int nof77)
/*
 * This file tests that all C compilers work and interact w/o any changes,
 * and figure out how to have the fortran compiler call the C compiler
 */
{
   char cmnd[2048], res[1024];
   char *sp;
   int i, ierr;
   if (verb)
      fprintf(stdout, "C compiler interoperation probe unimplemented!\n\n");
/*
 * C interoperation checks
 */
   if (verb > 1)
      fprintf(stderr, "ICC interoperation tests:\n");
   for (i=0; i < NCOMP; i++)
   {
      if (i != XCC_ && i != F77_ && i != ICC_)
      {
         if (strcmp(comps[i], comps[ICC_])) /* only check if different */
         {
            sprintf(cmnd,
   "make IRunC2C CC='%s' CCFLAGS='%s' CC1='%s' CC1FLAGS='%s' | fgrep SUCCESS",
                    comps[ICC_], comps[NCOMP+ICC_], comps[i], comps[i+NCOMP]);
            if (verb > 1)
               fprintf(stderr, "cmnd='%s'\n", cmnd);
            ierr = CmndOneLine(NULL, cmnd, res);
            if (!ierr)
               ierr = !strstr(res, "SUCCESS");
            if (ierr)
            {
               fprintf(stderr, "Compiler %d (%s) does not interoperate with interface compiler (%s), aborting!\n", i, comps[i], comps[ICC_]);
               fprintf(stderr, "ierr=%d, res='%s'\n", ierr, res);
               exit(ierr);
            }
            if (verb > 1)
               fprintf(stderr,
                       "   C2C %s/%s -- SUCCESS\n", comps[ICC_], comps[i]);
         }
      }
   }
/*
 * F2c tests
 */
   if (nof77)
   {
      *f2cnam = f2c_NamErr;
      *f2cint = f2c_IntErr;
      *f2cstr = f2c_StrErr;
   }
   else
   {
      sprintf(cmnd, "make IRun_f2c args=\"%s -C ic '%s' -F ic '%s' -C if '%s' -F if '%s'\" | fgrep 'F2C=('",
              targarg, comps[ICC_], comps[ICC_+NCOMP],
              comps[F77_], comps[F77_+NCOMP]);
      *f2cnam = f2c_NamErr;
      *f2cint = f2c_IntErr;
      *f2cstr = f2c_StrErr;
      if (verb > 1)
         fprintf(stderr, "cmnd='%s'\n", cmnd);
      if (!CmndOneLine(NULL, cmnd, res))
      {
         if (verb > 1)
            fprintf(stderr, "res='%s'\n", res);
         i = sscanf(res, " F2C=(%d,%d,%d)", f2cnam, f2cint, f2cstr);
         if (verb > 1)
            fprintf(stderr, "nread=%d, f2cname=%d, f2cint=%d, f2cstr=%d\n",
                    i, *f2cnam, *f2cint, *f2cstr);
         if (i != 3)
           *f2cnam = *f2cint = *f2cstr = 0;
      }
      if (verb)
      {
         printf("F2C name = %s\n", f2c_namestr[*f2cnam]);
         printf("F2C int  = %s\n", f2c_intstr[*f2cint]);
         printf("F2C str  = %s\n", f2c_strstr[*f2cstr]);
      }
   }
}

void PrintCompResults(char *file, char **comps, enum F2CNAME f2cnam,
                      enum F2CINT f2cint, enum F2CSTRING f2cstr)
{
   FILE *fpout;
   int i;

   if (file)
      fpout = fopen(file, "w");
   else fpout = stdout;
   assert(fpout);

   for (i=0; i < NCOMP; i++)
   {
      if (comps[i])
         fprintf(fpout, "%d '%s' '%s'\n", i, comps[i], comps[i+NCOMP]);
   }
   if (comps[F77_])
      fprintf(fpout, "F2CNAME,F2CINT,F2CSTRING=(%d,%d,%d)\n",
              f2cnam, f2cint, f2cstr);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
}


void PrintUsage(char *name, int iarg, char *arg)
{
   fprintf(stderr, "\nERROR around arg %d (%s).\n", iarg,
           arg ? arg : "unknown");
   fprintf(stderr, "USAGE: %s [flags] where flags are:\n", name);
   fprintf(stderr, "   -v <verb> : verbosity level\n");
   fprintf(stderr, "   -O <enum OSTYPE #>  : set OS type\n");
   fprintf(stderr, "   -A <enum MACHTYPE #> : set machine/architecture\n");
   fprintf(stderr,
   "   -V #    # = ((1<<vecISA1) | (1<<vecISA2) | ... | (1<<vecISAN))\n");
   fprintf(stderr, "   -b <32/64> : set pointer bitwidth\n");
   fprintf(stderr, "   -o <outfile>\n");
   fprintf(stderr, "   -C [xc,ic,if,sk,dk,sm,dm,al,ac] <compiler>\n");
   fprintf(stderr, "   -F [xc,ic,if,sk,dk,sm,dm,al,ac,gc] '<comp flags>'\n");
   fprintf(stderr,    /* HERE */
           "   -Fa [xc,ic,if,sk,dk,sm,dm,al,ac,gc] '<comp flags to append>'\n");
   fprintf(stderr, "        al: append flags to all compilers\n");
   fprintf(stderr, "        ac: append flags to all C compilers\n");
   fprintf(stderr,
      "   -T <targ> : ssh target for cross-compilation (probably broken)\n");
   fprintf(stderr, "   -S[i/s] <handle> <val>  : special int/string arg\n");
      fprintf(stderr,
        "      -Si nof77 <0/1> : Have/don't have fortran compiler\n");
   fprintf(stderr,
      "NOTE: enum #s can be found by : make xprint_enums ; ./xprint_enums\n");
   exit(iarg);
}

void GetFlags(int nargs,                /* nargs as passed into main */
              char **args,              /* args as passed into main */
              int *verb,                /* verbosity setting */
              enum OSTYPE *OS,          /* OS to assume */
              int *vec,                 /* Vector ISA extension bitfield */
              enum MACHTYPE *mach,     /* machine/arch to assume */
              int *ptrbits             /* # of bits in ptr: 32/64 */,
              char **comps,
              char **outfile,
              int *NoF77,
              char **targ             /* mach to ssh to*/
             )
{
   int i, k, k0, kn, DoInt;
   char *sp, *sp0;

   *verb = 0;
   *outfile = NULL;
   *targ = NULL;
   for (k=0; k < NCOMP*3; k++)
      comps[k] = NULL;

   *ptrbits = 0;
   *mach = 0;
   *vec = 0;
   *OS = 0;
   *verb = 0;
   *NoF77 = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'b':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *ptrbits = atoi(args[i]);
         break;
      case 'A':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *mach = atoi(args[i]);
         break;
      case 'V':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *vec = atoi(args[i]);
         break;
      case 'O':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *OS = atoi(args[i]);
         break;
      case 'v':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *verb = atoi(args[i]);
         break;
      case 'T':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *targ = args[i];
         break;
      case 'S':
         if (args[i][2] != 'i' && args[i][2] != 's')
            PrintUsage(args[0], i, "-S needs i or s suffix!");
         DoInt = args[i][2] == 'i';
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         sp0 = args[i];
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         if (DoInt)
            k = atoi(args[i]);
         else
            sp = args[i];
          if (!strcmp(sp0, "nof77"))
            *NoF77 = k;
         else
            PrintUsage(args[0], i-1, sp0);
         break;
      case 'o':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *outfile = args[i];
         break;
      case 'C':
      case 'F':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         sp = args[i];
         k = -1;
         if (*sp == 'i' && sp[1] == 'c') k = ICC_;
         else if (*sp == 'i' && sp[1] == 'f') k = F77_;
         else if (*sp == 's' && sp[1] == 'k') k = SKC_;
         else if (*sp == 'd' && sp[1] == 'k') k = DKC_;
         else if (*sp == 's' && sp[1] == 'm') k = SMC_;
         else if (*sp == 'd' && sp[1] == 'm') k = DMC_;
         else if (*sp == 'x' && sp[1] == 'c') k = XCC_;
         if (*sp == 'a' && (sp[1] == 'l' || sp[1] == 'c'))
         {  /* only appended flags can be applied to all compilers */
            if (args[i-1][1] == 'F')
            {
               if (args[i-1][2] == 'a')
               {
                  k0 = NCOMP+NCOMP;
                  kn = k0 + NCOMP;
               }
               else
               {
                  k0 = NCOMP;
                  kn = NCOMP+NCOMP;
               }
            }
            else
            {
               k0 = 0;
               kn = NCOMP;
            }
            if (++i >= nargs)
               PrintUsage(args[0], i, "out of arguments");
            for (k=k0; k < kn; k++)
               if (sp[1] == 'l' || k-2*NCOMP != F77_)
                  comps[k] = args[i];
         }
         else
         {
            if (k < 0) PrintUsage(args[0], i, args[i]);
            if (args[i-1][1] == 'F')
            {
               k += NCOMP;
               if (args[i-1][2] == 'a')
                  k += NCOMP;
            }
            if (++i >= nargs)
               PrintUsage(args[0], i, "out of arguments");
            comps[k] = args[i];
         }
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
/*
 * allocate these strings ourselves so we can free them later if necessary
 */
   for (i=0; i < 3*NCOMP; i++)
   {
      if (comps[i])
      {
         if (!strcmp(comps[i], "default"))
            comps[i] = NULL;
         else
         {
            sp = malloc(sizeof(char)*(strlen(comps[i])+1));
            strcpy(sp, comps[i]);
            comps[i] = sp;
         }
      }
   }
   if (*ptrbits != 32 && *ptrbits != 64)
      *ptrbits = 0;
}

main(int nargs, char **args)
/*
 * probe_comp has the following responsibilities:
 * 1. Read in atlcomp.txt for recommended compiler and flags
 * 2. Accept user override of compiler/flags
 * 3. Append any user appended flags
 * 4. Ensure all non-xcc C compilers interoperate by calling probe_ccomps
 * 5. Figure out F77/C interoperating rules by calling probe_f772c
 * 6. Printing out results of these probes for later use
 */
{
   enum OSTYPE OS;
   enum MACHTYPE mach;
   int ptrbits, verb, vecexts, i, nof77;
   char *usrcomps[3*NCOMP];
   char *outfile, *targ, *targarg;
   enum F2CNAME f2cnam;
   enum F2CINT f2cint;
   enum F2CSTRING f2cstr;

   GetFlags(nargs, args, &verb, &OS, &vecexts, &mach, &ptrbits, usrcomps,
            &outfile, &nof77, &targ);

   if (verb > 1)
   {
      fprintf(stdout, "User Override Compilers:\n");
      for (i=0; i < NCOMP; i++)
         fprintf(stdout, "   '%s' : '%s' '%s'\n",
            usrcomps[i] ? usrcomps[i]:"none",
            usrcomps[i+NCOMP] ? usrcomps[i+NCOMP]:"none",
            usrcomps[i+2*NCOMP] ? usrcomps[i+2*NCOMP]:"none");
      fprintf(stdout, "\n");
   }
   if (targ)
   {
      targarg = malloc(sizeof(char)*(strlen(targ)+24));
      assert(targarg);
      sprintf(targarg, "-T '%s'", targ);
   }
   else
      targarg = "";
   GetComps(OS, mach, verb, targ, ptrbits, usrcomps, nof77);
   TestComps(OS, mach, verb, targ, targarg, usrcomps,
             &f2cnam, &f2cint, &f2cstr, nof77);
   if (verb)
   {
      fprintf(stdout, "Compilers:\n");
      for (i=0; i < NCOMP; i++)
         fprintf(stdout, "   '%s' : '%s'\n", usrcomps[i], usrcomps[NCOMP+i]);
      fprintf(stdout, "\n");
   }
   PrintCompResults(outfile, usrcomps,  f2cnam, f2cint, f2cstr);
   if (targ)
      free(targarg);
   exit(0);
}
