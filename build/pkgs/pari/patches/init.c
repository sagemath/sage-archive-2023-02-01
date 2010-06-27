/* $Id: init.c 12456 2010-06-23 15:25:47Z kb $

Copyright (C) 2000-2003  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/*******************************************************************/
/*                                                                 */
/*        INITIALIZING THE SYSTEM, ERRORS, STACK MANAGEMENT        */
/*                                                                 */
/*******************************************************************/
#include <string.h>
#include <stdlib.h>
#include "pari.h"
#include "paripriv.h"
#include "anal.h"
#ifdef _WIN32
#  ifndef WINCE
#    include <process.h>
#  endif
#endif

const double LOG2    = 0.6931471805599453; /* log(2) */
const double LOG10_2 = 0.3010299956639812; /* log_10(2) */
const double LOG2_10 = 3.321928094887362;  /* log_2(10) */

GEN     gnil, gen_0, gen_1, gen_m1, gen_2, gen_m2, ghalf;
THREAD GEN    bernzone;
GEN     primetab; /* private primetable */
byteptr diffptr;
FILE    *pari_outfile, *pari_errfile, *pari_logfile, *pari_infile;
char    *current_logfile, *current_psfile, *pari_datadir;
long    gp_colors[c_LAST];
int     disable_color;
ulong   DEBUGFILES, DEBUGLEVEL, DEBUGMEM;
ulong   compatible, precreal, precdl, logstyle;
gp_data *GP_DATA;

GEN pari_colormap, pari_graphcolors;

entree  **varentries;

THREAD pari_sp bot, top, avma;
THREAD size_t memused;

static void ** MODULES, ** OLDMODULES;
static pari_stack s_MODULES, s_OLDMODULES;
const long functions_tblsz = 135; /* size of functions_hash */
entree **functions_hash;

void *foreignHandler;                       /* Handler for foreign commands.   */
char foreignExprSwitch = 3;               /* Just some unprobable char.      */
GEN  (*foreignExprHandler)(char*);    /* Handler for foreign expressions.*/
entree* (*foreignAutoload)(const char*, long len); /* Autoloader         */
void (*foreignFuncFree)(entree *);    /* How to free external entree.    */

void (*cb_pari_ask_confirm)(const char *);
int  (*cb_pari_handle_exception)(long);
int  (*cb_pari_whatnow)(const char *, int);
void (*cb_pari_sigint)(void);
void (*cb_pari_err_recover)(long);

typedef struct {
  jmp_buf *penv;
  long flag;
} cell;

static THREAD pari_stack s_ERR_CATCH;
static THREAD cell *ERR_CATCH;
THREAD void *global_err_data;

const long CATCH_ALL = 1;

/*********************************************************************/
/*                                                                   */
/*                       BLOCKS & CLONES                             */
/*                                                                   */
/*********************************************************************/
/*#define DEBUG*/
static THREAD long next_block;
static THREAD GEN cur_block; /* current block in block list */
#ifdef DEBUG
static THREAD long NUM;
#endif

static void
pari_init_blocks(void)
{
  next_block = 0; cur_block = 0;
#ifdef DEBUG
  NUM = 0;
#endif
}

static void
pari_close_blocks(void)
{
  while (cur_block) killblock(cur_block);
}

/* Return x, where:
 * x[-4]: reference count
 * x[-3]: adress of next block
 * x[-2]: adress of preceding block.
 * x[-1]: number of allocated blocs.
 * x[0..n-1]: malloc-ed memory. */
GEN
newblock(size_t n)
{
  long *x = (long *) pari_malloc((n + BL_HEAD)*sizeof(long)) + BL_HEAD;

  bl_refc(x) = 1;
  bl_next(x) = NULL;
  bl_prev(x) = cur_block;
  bl_num(x)  = next_block++;
  if (cur_block) bl_next(cur_block) = x;
#ifdef DEBUG
  fprintferr("+ %ld\n", ++NUM);
#endif
  if (DEBUGMEM)
  {
    if (!n) pari_warn(warner,"mallocing NULL object in newblock");
    if (DEBUGMEM > 2)
      fprintferr("new block, size %6lu (no %ld): %08lx\n", n, next_block-1, x);
  }
  return cur_block = x;
}

void
gclone_refc(GEN x) { ++bl_refc(x); }

void
gunclone(GEN x)
{
  if (--bl_refc(x) > 0) return;
  BLOCK_SIGINT_START;
  if (bl_next(x)) bl_prev(bl_next(x)) = bl_prev(x);
  else
  {
    cur_block = bl_prev(x);
    next_block = bl_num(x);
  }
  if (bl_prev(x)) bl_next(bl_prev(x)) = bl_next(x);
  if (DEBUGMEM > 2)
    fprintferr("killing block (no %ld): %08lx\n", bl_num(x), x);
  free((void*)bl_base(x)); /* pari_free not needed: we already block */
  BLOCK_SIGINT_END;
#ifdef DEBUG
  fprintferr("- %ld\n", NUM--);
#endif
}

/* Recursively look for clones in the container and kill them. Then kill
 * container if clone. SIGINT could be blocked until it returns */
void
gunclone_deep(GEN x)
{
  long i, lx;
  GEN v;
  BLOCK_SIGINT_START;
  switch(typ(x))
  {
    case t_VEC: case t_COL: case t_MAT:
      lx = lg(x);
      for (i=1;i<lx;i++) gunclone_deep(gel(x,i));
      break;
    case t_LIST:
      v = list_data(x); lx = v? lg(v): 1;
      for (i=1;i<lx;i++) gunclone_deep(gel(v,i));
      pari_free(v); break;
  }
  if (isclone(x)) gunclone(x);
  BLOCK_SIGINT_END;
}

int
pop_entree_block(entree *ep, long loc)
{
  GEN x = (GEN)ep->value;
  if (bl_num(x) < loc) return 0; /* older */
  if (DEBUGMEM>2)
    fprintferr("popping %s (block no %ld)\n", ep->name, bl_num(x));
  gunclone_deep(x); return 1;
}

/*********************************************************************/
/*                                                                   */
/*                       C STACK SIZE CONTROL                        */
/*                                                                   */
/*********************************************************************/
/* Avoid core dump on deep recursion. Adapted Perl code by Dominic Dunlop */
THREAD void *PARI_stack_limit = NULL;

#ifdef STACK_CHECK

#  ifdef __EMX__                                /* Emulate */
void
pari_stackcheck_init(void *stack_base)
{
  (void) stack_base;
  if (!stack_base) { PARI_stack_limit = NULL; return; }
  PARI_stack_limit = get_stack(1./16, 32*1024);
}
#  else /* !__EMX__ */
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

/* Set PARI_stack_limit to (a little above) the lowest safe address that can be
 * used on the stack. Leave PARI_stack_limit at its initial value (NULL) to
 * show no check should be made [init failed]. Assume stack grows downward. */
void
pari_stackcheck_init(void *stack_base)
{
  struct rlimit rip;
  ulong size;
  if (!stack_base) { PARI_stack_limit = NULL; return; }
  if (getrlimit(RLIMIT_STACK, &rip)) return;
  size = rip.rlim_cur;
  if (size == (ulong)RLIM_INFINITY || size > (ulong)stack_base)
    PARI_stack_limit = (void*)(((ulong)stack_base) / 16);
  else
    PARI_stack_limit = (void*)((ulong)stack_base - (size/16)*15);
}
#  endif /* !__EMX__ */

#else
void
pari_stackcheck_init(void *stack_base)
{
  PARI_stack_limit = NULL;
}
#endif /* STACK_CHECK */

/*******************************************************************/
/*                         HEAP TRAVERSAL                          */
/*******************************************************************/
struct getheap_t { long n, l; };
static void
f_getheap(GEN x, void *D)
{
  struct getheap_t *T = (struct getheap_t*)D;
  T->n++;
  T->l += gsizeword(x);
}
GEN
getheap(void)
{
  struct getheap_t T = { 0, 0 };
  traverseheap(&f_getheap, &T);
  return mkvec2s(T.n, T.l + BL_HEAD * T.n);
}

void
traverseheap( void(*f)(GEN, void *), void *data )
{
  GEN x;
  for (x = cur_block; x; x = bl_prev(x)) f(x, data);
}

/*********************************************************************/
/*                          DAEMON / FORK                            */
/*********************************************************************/
#if defined(HAS_WAITPID) && defined(HAS_SETSID)
#  include <sys/wait.h>
/* Properly fork a process, detaching from main process group without creating
 * zombies on exit. Parent returns 1, son returns 0 */
int
pari_daemon(void)
{
  pid_t pid = fork();
  switch(pid) {
      case -1: return 1; /* father, fork failed */
      case 0:
        (void)setsid(); /* son becomes process group leader */
        if (fork()) exit(0); /* now son exits, also when fork fails */
        break; /* grandson: its father is the son, which exited,
                * hence father becomes 'init', that'll take care of it */
      default: /* father, fork succeeded */
        (void)waitpid(pid,NULL,0); /* wait for son to exit, immediate */
        return 1;
  }
  /* grandson */
  return 0;
}
#else
int
pari_daemon(void)
{
  pari_err(impl,"pari_daemon without waitpid & setsid");
  return 0;
}
#endif

/*********************************************************************/
/*                                                                   */
/*                       SYSTEM INITIALIZATION                       */
/*                                                                   */
/*********************************************************************/
static int try_to_recover = 0;
VOLATILE int PARI_SIGINT_block  = 0, PARI_SIGINT_pending = 0;
static void pari_sighandler(int sig);

/*********************************************************************/
/*                         SIGNAL HANDLERS                           */
/*********************************************************************/
static void
dflt_sigint_fun(void) { pari_err(talker, "user interrupt"); }

#if defined(_WIN32) || defined(__CYGWIN32__)
int win32ctrlc = 0;
void
dowin32ctrlc(void)
{
  win32ctrlc = 0;
  cb_pari_sigint();
}
#endif

static void
pari_handle_SIGINT(void)
{
#ifdef _WIN32
  if (++win32ctrlc >= 5) _exit(3);
#else
  cb_pari_sigint();
#endif
}

static void
pari_sighandler(int sig)
{
  const char *msg;
#ifndef HAS_SIGACTION
  /*SYSV reset the signal handler in the handler*/
  (void)os_signal(sig,pari_sighandler);
#endif
  switch(sig)
  {
#ifdef SIGBREAK
    case SIGBREAK:
      if (PARI_SIGINT_block) PARI_SIGINT_pending=SIGBREAK;
      else pari_handle_SIGINT();
      return;
#endif

#ifdef SIGINT
    case SIGINT:
      if (PARI_SIGINT_block) PARI_SIGINT_pending=SIGINT;
      else pari_handle_SIGINT();
      return;
#endif

#ifdef SIGSEGV
    case SIGSEGV:
      msg="PARI/GP (Segmentation Fault)"; break;
#endif
#ifdef SIGBUS
    case SIGBUS:
      msg="PARI/GP (Bus Error)"; break;
#endif
#ifdef SIGFPE
    case SIGFPE:
      msg="PARI/GP (Floating Point Exception)"; break;
#endif

#ifdef SIGPIPE
    case SIGPIPE:
    {
      pariFILE *f = GP_DATA->pp->file;
      if (f && pari_outfile == f->file)
      {
        pari_err(talker, "Broken Pipe, resetting file stack...");
        GP_DATA->pp->file = NULL; /* to avoid oo recursion on error */
        pari_outfile = stdout; pari_fclose(f);
      }
      /*Do not attempt to write to stdout in case it triggered the SIGPIPE*/
      return; /* not reached */
    }
#endif

    default: msg="signal handling"; break;
  }
  pari_err(bugparier,msg);
}

void
pari_sig_init(void (*f)(int))
{
#ifdef SIGBUS
  (void)os_signal(SIGBUS,f);
#endif
#ifdef SIGFPE
  (void)os_signal(SIGFPE,f);
#endif
#ifdef SIGINT
  (void)os_signal(SIGINT,f);
#endif
#ifdef SIGBREAK
  (void)os_signal(SIGBREAK,f);
#endif
#ifdef SIGPIPE
  (void)os_signal(SIGPIPE,f);
#endif
#ifdef SIGSEGV
  (void)os_signal(SIGSEGV,f);
#endif
}

/*********************************************************************/
/*                      STACK AND UNIVERSAL CONSTANTS                */
/*********************************************************************/
static void
init_universal_constants(void)
{
  /* 2 (gnil) + 2 (gen_0)
   * + 3 (gen_1) + 3 (gen_m1) + 3 (gen_2) + 3 (gen_m2)
   * + 3 (half) */
  static long universal_constants[19];
  GEN p = (GEN) universal_constants;
  gen_0 = p; p+=2; gnil = p; p+=2;
  gen_0[0] = gnil[0] = evaltyp(t_INT) | _evallg(2);
  gen_0[1] = gnil[1] = evallgefint(2);

  gen_1 = p; p+=3;
  gen_2 = p; p+=3;
  gen_1[0] = gen_2[0] = evaltyp(t_INT) | _evallg(3);
  gen_1[1] = gen_2[1] = evalsigne(1) | evallgefint(3);
  gen_1[2] = 1; gen_2[2]= 2;

  gen_m1 = p; p+=3;
  gen_m1[0] = evaltyp(t_INT) | _evallg(3);
  gen_m1[1] = evalsigne(-1) | evallgefint(3);
  gen_m1[2] = 1;

  gen_m2 = p; p+=3;
  gen_m2[0] = evaltyp(t_INT) | _evallg(3);
  gen_m2[1] = evalsigne(-1) | evallgefint(3);
  gen_m2[2] = 2;

  ghalf = p;
  ghalf[0] = evaltyp(t_FRAC) | _evallg(3);
  gel(ghalf,1) = gen_1;
  gel(ghalf,2) = gen_2;
}

static size_t
fix_size(size_t a)
{
  size_t b = (a / sizeof(long)) * sizeof(long); /* Align */
  if (b < 1024) b = 1024;
  return b;
}
/* old = current stack size (0 = unallocated), size = new size */
void
pari_init_stack(size_t size, size_t old)
{
  size_t s = fix_size(size);
  if (old != s) {
    if (old) pari_free((void*)bot);
    BLOCK_SIGINT_START;
    for (;; s>>=1)
    {
      char buf[128];
      bot = (pari_sp)malloc(s); /* NOT pari_malloc, memer would be deadly */
      if (bot) break;
      if (!s) pari_err(memer); /* no way out. Die */
      /* must use sprintf: pari stack is currently dead */
      sprintf(buf, "not enough memory, new stack %lu", (ulong)s);
      pari_warn(warner, buf, s);
    }
    BLOCK_SIGINT_END;
  }
  avma = top = bot+s;
  memused = 0;
}

/*********************************************************************/
/*                           INIT DEFAULTS                           */
/*********************************************************************/
void
pari_init_defaults(void)
{
  long i;
  initout(1);

#ifdef LONG_IS_64BIT
  precreal = 4;
#else
  precreal = 5;
#endif

  precdl = 16;
  compatible = NONE;
  DEBUGFILES = DEBUGLEVEL = DEBUGMEM = 0;
  disable_color = 1;
  logstyle = logstyle_none;

  current_psfile = pari_strdup("pari.ps");
  current_logfile= pari_strdup("pari.log");
  pari_logfile = NULL;

  pari_datadir = os_getenv("GP_DATA_DIR");
  if (!pari_datadir) pari_datadir = (char*)GPDATADIR;
  if (pari_datadir) pari_datadir = pari_strdup(pari_datadir);

  for (i=0; i<c_LAST; i++) gp_colors[i] = c_NONE;
  pari_colormap = NULL; pari_graphcolors = NULL;
  (void)sd_graphcolormap("[\"white\",\"black\",\"blue\",\"violetred\",\"red\",\"green\",\"grey\",\"gainsboro\"]", d_SILENT);
  (void)sd_graphcolors("[4, 5]", d_SILENT);
}

/*********************************************************************/
/*                   FUNCTION HASHTABLES, MODULES                    */
/*********************************************************************/

/* Initialize hashtable */
static void
init_hashtable(entree **table, long tblsz)
{
  long i;
  for (i = 0; i < tblsz; i++)
  {
    entree *last = NULL, *ep = table[i];
    table[i] = NULL;
    while (ep)
    {
      entree *EP = ep->next;
      switch(EpVALENCE(ep))
      {
        case EpVAR: case EpINSTALL:
        /* keep: attach it to last entree seen */
          if (last)
            last->next = ep;
          else
            table[i] = ep;
          ep->next = NULL; last = ep;
          break;
        default: freeep(ep);
      }
      ep = EP;
    }
  }
}
/* Load in hashtable hash the modules contained in A */
static int
gp_init_entrees(pari_stack *p_A, entree **hash)
{
  long i;
  entree **v = (entree **)*stack_base(p_A);
  init_hashtable(hash, functions_tblsz);
  for (i = 0; i < p_A->n; i++) pari_fill_hashtable(hash, v[i]);
  return (hash == functions_hash);
}
int
gp_init_functions(void)
{
  return gp_init_entrees(new_fun_set? &s_MODULES: &s_OLDMODULES, functions_hash);
}

static void
pari_init_functions(void)
{
  stack_init(&s_MODULES, sizeof(*MODULES),(void**)&MODULES);
  stack_pushp(&s_MODULES,functions_basic);
  stack_init(&s_OLDMODULES, sizeof(*OLDMODULES),(void**)&OLDMODULES);
  stack_pushp(&s_OLDMODULES,oldfonctions);
  functions_hash = (entree**) pari_calloc(sizeof(entree*)*functions_tblsz);
  pari_fill_hashtable(functions_hash,
                      new_fun_set? functions_basic:oldfonctions);
}

void
pari_add_module(entree *ep)
{
  if (new_fun_set)
    pari_fill_hashtable(functions_hash, ep);
  stack_pushp(&s_MODULES, ep);
}

void
pari_add_oldmodule(entree *ep)
{
  if (!new_fun_set)
    pari_fill_hashtable(functions_hash, ep);
  stack_pushp(&s_OLDMODULES, ep);
}

static void
pari_init_errcatch(void)
{
  stack_init(&s_ERR_CATCH, sizeof(cell), (void**)&ERR_CATCH);
  global_err_data = NULL;
}

/*********************************************************************/
/*                       PARI THREAD                                 */
/*********************************************************************/

static void
pari_stack_alloc(struct pari_mainstack *st, size_t s)
{
  st->bot = (pari_sp)pari_malloc(s);
  st->avma = st->top = st->bot+s;
  st->memused = 0;
}

static void
pari_stack_free(struct pari_mainstack *st)
{
  free((void*)st->bot);
  st->avma = st->top = st->bot = 0;
}

static void
pari_stack_use(struct pari_mainstack *st)
{
  bot = st->bot; top = st->top; avma = st->avma;
  memused = st->memused;
}

/* Initial PARI thread structure t with a stack of size s and
 * argument arg */

void
pari_thread_alloc(struct pari_thread *t, size_t s, GEN arg)
{
  pari_stack_alloc(&t->st,s);
  t->data = arg;
}

void
pari_thread_free(struct pari_thread *t)
{
  pari_stack_free(&t->st);
}

void
pari_thread_init(void)
{
  pari_init_blocks();
  pari_init_errcatch();
  pari_init_rand();
  pari_init_floats();
  pari_init_seadata();
  pari_init_parser();
  pari_init_compiler();
  pari_init_evaluator();
  pari_init_files();
}

void
pari_thread_close(void)
{
  pari_close_files();
  pari_close_evaluator();
  pari_close_compiler();
  pari_close_parser();
  pari_close_seadata();
  pari_close_floats();
  pari_close_blocks();
}

GEN
pari_thread_start(struct pari_thread *t)
{
  pari_stack_use(&t->st);
  pari_thread_init();
  return t->data;
}

/*********************************************************************/
/*                       LIBPARI INIT / CLOSE                        */
/*********************************************************************/

static void
pari_exit(void)
{

	abort();
	fprintferr("  ***   Error in the PARI system. End of program.\n");
  exit(1);
}

static void
dflt_err_recover(long errnum) { (void) errnum; pari_exit(); }

/* initialize PARI data. Initialize [new|old]fun to NULL for default set. */
void
pari_init_opts(size_t parisize, ulong maxprime, ulong init_opts)
{
  ulong u;

  cb_pari_whatnow = NULL;
  cb_pari_sigint = dflt_sigint_fun;
  cb_pari_handle_exception = NULL;
  cb_pari_err_recover = dflt_err_recover;

  pari_stackcheck_init(&u);
  if ((init_opts&INIT_DFTm)) {
    GP_DATA = default_gp_data();
    gp_expand_path(GP_DATA->path);
    pari_init_defaults();
  }

  if ((init_opts&INIT_SIGm)) pari_sig_init(pari_sighandler);
  pari_init_stack(parisize, 0);
  diffptr = initprimes(maxprime);
  init_universal_constants();
  if (pari_kernel_init()) pari_err(talker,"Cannot initialize kernel");

  primetab = cgetalloc(t_VEC, 1);
  varentries = (entree**) pari_calloc((MAXVARN+1)*sizeof(entree*));
  pari_thread_init();
  pari_init_functions();
  pari_var_init();
  try_to_recover = 1;
}

void
pari_init(size_t parisize, ulong maxprime)
{ pari_init_opts(parisize, maxprime, INIT_JMPm | INIT_SIGm | INIT_DFTm); }

void
pari_close_opts(ulong init_opts)
{
  long i;

  BLOCK_SIGINT_START;
  if ((init_opts&INIT_SIGm)) pari_sig_init(SIG_DFL);

  while (delete_var()) /* empty */;
  for (i = 0; i < functions_tblsz; i++)
  {
    entree *ep = functions_hash[i];
    while (ep) {
      entree *EP = ep->next;
      if (!EpSTATIC(ep)) { freeep(ep); free(ep); }
      ep = EP;
    }
  }
  free((void*)varentries);
  free((void*)primetab);
  pari_thread_close();

  free((void*)functions_hash);
  free((void*)bot);
  free((void*)diffptr);
  free(current_logfile);
  free(current_psfile);
  stack_delete(&s_MODULES);
  stack_delete(&s_OLDMODULES);
  if (pari_datadir) free(pari_datadir);
  if (init_opts&INIT_DFTm)
  { /* delete GP_DATA */
    if (GP_DATA->hist->res) free((void*)GP_DATA->hist->res);
    if (GP_DATA->pp->cmd) free((void*)GP_DATA->pp->cmd);
    if (GP_DATA->help) free((void*)GP_DATA->help);
    delete_dirs(GP_DATA->path);
    free((void*)GP_DATA->path->PATH);
  }
  BLOCK_SIGINT_END;
}

void
pari_close(void)
{ pari_close_opts(INIT_JMPm | INIT_SIGm | INIT_DFTm); }

/*******************************************************************/
/*                                                                 */
/*                         ERROR RECOVERY                          */
/*                                                                 */
/*******************************************************************/

void
gp_context_save(struct gp_context* rec)
{
  rec->file = pari_last_tmp_file();
  rec->listloc = next_block;
  evalstate_save(&rec->eval);
  parsestate_save(&rec->parse);
}

void
gp_context_restore(struct gp_context* rec)
{
  pariFILE *f;
  long i;

  if (!(GP_DATA->flags & RECOVER)) pari_exit();
  if (!try_to_recover) return;
  /* disable gp_context_restore() and SIGINT */
  try_to_recover = 0;
  BLOCK_SIGINT_START
  if (DEBUGMEM>2) fprintferr("entering recover(), loc = %ld\n", rec->listloc);
  evalstate_restore(&rec->eval);
  parsestate_restore(&rec->parse);

  /* delete all "temp" files open since last reference point */
  f = pari_last_tmp_file();
  while (f)
  {
    pariFILE *g = f->prev;
    if (f == rec->file) break;
    pari_fclose(f); f = g;
  }

  for (i = 0; i < functions_tblsz; i++)
  {
    entree *ep = functions_hash[i];
    while (ep)
    {
      entree *EP = ep->next;
      switch(EpVALENCE(ep))
      {
        case EpVAR:
          while (pop_val_if_newer(ep,rec->listloc)) /* empty */;
          break;
        case EpNEW: break;
      }
      ep = EP;
    }
  }
  if (DEBUGMEM>2) fprintferr("leaving recover()\n");
  BLOCK_SIGINT_END
  try_to_recover = 1;
}

long
err_catch(long errnum, jmp_buf *penv)
{
  long n;
  /* for fear of infinite recursion */
  if (errnum == memer) pari_err(talker, "can't trap memory errors");
  if (errnum == CATCH_ALL) errnum = noer;
  else if (errnum > noer) pari_err(talker, "no such error number: %ld", errnum);
  n = stack_new(&s_ERR_CATCH);
  ERR_CATCH[n].flag = errnum;
  ERR_CATCH[n].penv = penv; return n;
}

/* delete traps younger than n (included) */
void
err_leave(long n) { if (n >= 0) s_ERR_CATCH.n = n; }

/* Get last (most recent) handler for error n (or generic noer) killing all
 * more recent non-applicable handlers (now obsolete) */
static cell *
err_seek(long n)
{
  if (n <= bugparier) return NULL;
  for( ; s_ERR_CATCH.n; s_ERR_CATCH.n--)
  {
    cell *t = &ERR_CATCH[ s_ERR_CATCH.n-1 ];
    if (t->flag == n || t->flag == noer) return t;
  }
  return NULL;
}

void
err_recover(long numerr)
{
  evalstate_reset();
  parsestate_reset();
  initout(0);
  dbg_release();
  killallfiles();
  s_ERR_CATCH.n = 0; /* untrapped error: kill all error handlers */
  global_err_data = NULL;
  fprintferr("\n"); flusherr();

  cb_pari_err_recover(numerr);
}

static void
err_init(void)
{
  /* make sure pari_err msg starts at the beginning of line */
  if (!pari_last_was_newline()) pari_putc('\n');
  pariOut->flush();
  pariErr->flush();
  pariOut = pariErr;
  term_color(c_ERR);
}

static void
err_init_msg(int numerr)
{
  const char *gp_function_name;
  pari_puts("  *** ");
  if (numerr != user && (gp_function_name = closure_func_err()))
    pari_printf("%s: ", gp_function_name);
  else
    pari_puts("  ");
}

void
pari_warn(int numerr, ...)
{
  char *ch1;
  PariOUT *out = pariOut;
  va_list ap;

  va_start(ap,numerr);

  err_init();
  err_init_msg(numerr);
  switch (numerr)
  {
    case user:
      pari_puts("user warning: ");
      print0(va_arg(ap, GEN), f_RAW);
      break;

    case warnmem:
      pari_puts("collecting garbage in "); ch1=va_arg(ap, char*);
      pari_vprintf(ch1,ap); pari_putc('.');
      break;

    case warner:
      pari_puts("Warning: "); ch1=va_arg(ap, char*);
      pari_vprintf(ch1,ap); pari_putc('.');
      break;

    case warnprec:
      pari_vprintf("Warning: increasing prec in %s; new prec = %ld",ap);
      break;

    case warnfile:
      pari_puts("Warning: failed to "),
      ch1 = va_arg(ap, char*);
      pari_printf("%s: %s", ch1, va_arg(ap, char*));
      break;
  }
  term_color(c_NONE); va_end(ap);
  pari_putc('\n');
  pariOut = out;
  flusherr();
}
void
pari_sigint(const char *s)
{
  PariOUT *out = pariOut;
  err_init();
  closure_err();
  err_init_msg(talker);
  pari_puts(s); pari_putc('.');
  term_color(c_NONE);
  pariOut = out;
  flusherr();
  if (cb_pari_handle_exception &&
      cb_pari_handle_exception(-1)) return;
  err_recover(talker);
}

void
pari_err(int numerr, ...)
{
  PariOUT *out = pariOut;
  va_list ap;

  va_start(ap,numerr);

  global_err_data = NULL;
  if (s_ERR_CATCH.n)
  {
    cell *trapped;
    if ( (trapped = err_seek(numerr)) )
    {
      switch(numerr)
      {
        case invmoder:
          global_err_data = (void*)va_arg(ap, GEN);
          break;
        case alarmer:
          global_err_data = (char*)va_arg(ap, char*);
          break;
      }
      longjmp(*(trapped->penv), numerr);
    }
  }
  err_init();
  if (numerr == talker2)
  {
    const char *msg = va_arg(ap, char*);
    const char *s = va_arg(ap,char *);
    print_errcontext(msg,s,va_arg(ap,char *));
  }
  else
  {
    closure_err();
    err_init_msg(numerr); pari_puts(errmessage[numerr]);
    switch (numerr)
    {
      case talker: case alarmer: {
        const char *ch1 = va_arg(ap, char*);
        pari_vprintf(ch1,ap); pari_putc('.'); break;
      }
      case user:
        pari_puts("user error: ");
        print0(va_arg(ap, GEN), f_RAW);
        break;
      case invmoder:
        pari_printf("impossible inverse modulo: %Ps.", va_arg(ap, GEN));
        break;
      case openfiler: {
        const char *type = va_arg(ap, char*);
        pari_printf("error opening %s file: `%s'.", type, va_arg(ap,char*));
        break;
      }
      case overflower:
        pari_printf("overflow in %s.", va_arg(ap, char*));
        break;
      case notfuncer:
      {
        GEN fun = va_arg(ap, GEN);
        if (gcmpX(fun))
        {
          entree *ep = varentries[varn(fun)];
          const char *s = ep->name;
          if (cb_pari_whatnow) cb_pari_whatnow(s,1);
        }
        break;
      }

      case impl:
        pari_printf("sorry, %s is not yet implemented.", va_arg(ap, char*));
        break;
      case typeer: case mattype1: case negexper:
      case constpoler: case notpoler: case redpoler:
      case zeropoler: case consister: case flagerr: case precer:
        pari_printf(" in %s.",va_arg(ap, char*)); break;

      case bugparier:
        pari_printf("bug in %s, please report",va_arg(ap, char*)); break;

      case operi: case operf:
      {
        const char *f, *op = va_arg(ap, const char*);
        GEN x = va_arg(ap, GEN);
        GEN y = va_arg(ap, GEN);
        pari_puts(numerr == operi? "impossible": "forbidden");
        switch(*op)
        {
          case '+': f = "addition"; break;
          case '-':
            pari_printf(" negation - %s.",type_name(typ(x)));
            f = NULL; break;
          case '*': f = "multiplication"; break;
          case '/': case '%': case '\\': f = "division"; break;
          case 'g': op = ","; f = "gcd"; break;
          default: op = "-->"; f = "assignment"; break;
        }
        if (f)
          pari_printf(" %s %s %s %s.",f,type_name(typ(x)),op,type_name(typ(y)));
        break;
      }

      case primer1: {
        ulong c = va_arg(ap, ulong);
        if (c) pari_printf(", need primelimit ~ %lu.", c);
        break;
      }
    }
  }
  term_color(c_NONE); va_end(ap);
  if (numerr==errpile)
  {
    size_t d = top - bot;
    char buf[256];
    /* don't use pari_printf: it needs the PARI stack for %.3f conversion */
    sprintf(buf, "\n  current stack size: %lu (%.3f Mbytes)\n",
                 (ulong)d, (double)d/1048576.);
    pariErr->puts(buf);
    pariErr->puts("  [hint] you can increase GP stack with allocatemem()\n");
  }
  pariOut = out;
  flusherr();
  if (cb_pari_handle_exception &&
      cb_pari_handle_exception(numerr)) return;
  err_recover(numerr);
}

/* Try f (trapping error e), recover using r (break_loop, if NULL) */
GEN
trap0(const char *e, GEN r, GEN f)
{
  long numerr = CATCH_ALL;
  GEN x;
       if (!strcmp(e,"errpile")) numerr = errpile;
  else if (!strcmp(e,"typeer")) numerr = typeer;
  else if (!strcmp(e,"gdiver")) numerr = gdiver;
  else if (!strcmp(e,"impl")) numerr = impl;
  else if (!strcmp(e,"invmoder")) numerr = invmoder;
  else if (!strcmp(e,"archer")) numerr = archer;
  else if (!strcmp(e,"alarmer")) numerr = alarmer;
  else if (!strcmp(e,"talker")) numerr = talker;
  else if (!strcmp(e,"user")) numerr = user;
  else if (*e) pari_err(impl,"this trap keyword");
  /* TODO: complete the list */

  if (!f) {
    pari_warn(warner,"default handlers are no longer supported --> ignored");
    return gnil;
  }
  /* explicit recovery text */
  x = closure_trapgen(f, numerr);
  if (x == (GEN)1L) x = r? closure_evalgen(r): gnil;
  return x;
}

/*******************************************************************/
/*                                                                */
/*                       CLONING & COPY                            */
/*                  Replicate an existing GEN                      */
/*                                                                 */
/*******************************************************************/
/* lontyp[tx] = 0 (non recursive type) or number of codewords for type tx */
const  long lontyp[] = { 0,0,0,1,1,2,1,2,1,1, 2,2,0,1,1,1,1,1,1,1, 2,0,0,2 };

static GEN
list_internal_copy(GEN z, long nmax)
{
  long i, l;
  GEN a;
  if (!z) return NULL;
  l = lg(z);
  a = (GEN)pari_malloc((nmax+1) * sizeof(long));
  for (i = 1; i < l; i++) gel(a,i) = gclone( gel(z,i) );
  a[0] = z[0]; return a;
}

static void
listassign(GEN x, GEN y)
{
  long nmax = list_nmax(x);
  GEN L = list_data(x);
  if (!nmax && L) nmax = lg(L) + 32; /* not malloc'ed yet */
  list_nmax(y) = nmax;
  list_data(y) = list_internal_copy(L, nmax);
}

/* copy list on the PARI stack */
GEN
listcopy(GEN x)
{
  GEN y = listcreate(), L = list_data(x);
  if (L) list_data(y) = gcopy(L);
  return y;
}

GEN
gcopy(GEN x)
{
  long tx = typ(x), lx, i;
  GEN y;
  switch(tx)
  { /* non recursive types */
    case t_INT: return signe(x)? icopy(x): gen_0;
    case t_REAL:
    case t_STR:
    case t_VECSMALL: return leafcopy(x);
    /* one more special case */
    case t_LIST: return listcopy(x);
  }
  y = cgetg_copy(x, &lx);
  if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
  for (; i<lx; i++) gel(y,i) = gcopy(gel(x,i));
  return y;
}

/* as gcopy, but truncate to the first lx components if recursive type
 * [ leaves use their own lg ]. No checks. */
GEN
gcopy_lg(GEN x, long lx)
{
  long tx = typ(x), i;
  GEN y;
  switch(tx)
  { /* non recursive types */
    case t_INT: return signe(x)? icopy(x): gen_0;
    case t_REAL:
    case t_STR:
    case t_VECSMALL: return leafcopy(x);
    /* one more special case */
    case t_LIST: return listcopy(x);
  }
  y = cgetg(lx, tx);
  if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
  for (; i<lx; i++) gel(y,i) = gcopy(gel(x,i));
  return y;
}

/* cf cgetg_copy: "allocate" (by updating first codeword only) for subsequent
 * copy of x, as if avma = *AVMA */
INLINE GEN
cgetg_copy_avma(GEN x, long *plx, pari_sp *AVMA) {
  GEN z;
  *plx = lg(x);
  z = ((GEN)*AVMA) - *plx;
  z[0] = x[0] & (TYPBITS|LGBITS);
  *AVMA = (pari_sp)z; return z;
}
INLINE GEN
cgetlist_avma(pari_sp *AVMA)
{
  GEN y = ((GEN)*AVMA) - 3;
  y[0] = _evallg(3) | evaltyp(t_LIST);
  *AVMA = (pari_sp)y; return y;
}

/* copy x as if avma = *AVMA, update *AVMA */
GEN
gcopy_avma(GEN x, pari_sp *AVMA)
{
  long i, lx, tx = typ(x);
  GEN y;

  switch(typ(x))
  { /* non recursive types */
    case t_INT:
      *AVMA = (pari_sp)icopy_avma(x, *AVMA);
      return (GEN)*AVMA;
    case t_REAL:
    case t_STR:
    case t_VECSMALL:
      y = cgetg_copy_avma(x, &lx, AVMA);
      for (i=1; i<lx; i++) y[i] = x[i];
      break;

    /* one more special case */
    case t_LIST:
      y = cgetlist_avma(AVMA);
      listassign(x, y); return y;

    default:
      y = cgetg_copy_avma(x, &lx, AVMA);
      if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
      for (; i<lx; i++) gel(y,i) = gcopy_avma(gel(x,i), AVMA);
  }
  return y;
}

/* [copy_bin/bin_copy:] same as gcopy_avma but use NULL to code an exact 0, and
 * make shallow copies of t_LISTs */
static GEN
gcopy_av0(GEN x, pari_sp *AVMA)
{
  long i, lx, tx = typ(x);
  GEN y;

  switch(tx)
  { /* non recursive types */
    case t_INT:
      if (!signe(x)) return NULL; /* special marker */
      *AVMA = (pari_sp)icopy_avma(x, *AVMA);
      return (GEN)*AVMA;
    case t_REAL:
    case t_STR:
    case t_VECSMALL:
    /* one more special case */
    case t_LIST:
      y = cgetg_copy_avma(x, &lx, AVMA);
      for (i=1; i<lx; i++) y[i] = x[i];
      break;
    default:
      y = cgetg_copy_avma(x, &lx, AVMA);
      if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
      for (; i<lx; i++) gel(y,i) = gcopy_av0(gel(x,i), AVMA);
  }
  return y;
}

INLINE GEN
icopy_avma_canon(GEN x, pari_sp AVMA)
{
  long i, lx = lgefint(x);
  GEN y = ((GEN)AVMA) - lx;
  y[0] = evaltyp(t_INT)|evallg(lx); /* kills isclone */
  y[1] = x[1]; x = int_MSW(x);
  for (i=2; i<lx; i++, x = int_precW(x)) y[i] = *x;
  return y;
}

/* [copy_bin_canon/bin_copy_canon:] same as gcopy_av0, but copy integers in
 * canonical (native kernel) form and make a full copy of t_LISTs */
static GEN
gcopy_av0_canon(GEN x, pari_sp *AVMA)
{
  long i, lx, tx = typ(x);
  GEN y;

  switch(tx)
  { /* non recursive types */
    case t_INT:
      if (!signe(x)) return NULL; /* special marker */
      *AVMA = (pari_sp)icopy_avma_canon(x, *AVMA);
      return (GEN)*AVMA;
    case t_REAL:
    case t_STR:
    case t_VECSMALL:
      y = cgetg_copy_avma(x, &lx, AVMA);
      for (i=1; i<lx; i++) y[i] = x[i];
      break;

    /* one more special case */
    case t_LIST:
    {
      GEN y = cgetlist_avma(AVMA), z = list_data(x);
      if (z) {
        list_data(y) = gcopy_av0_canon(z, AVMA);
        list_nmax(y) = lg(z)-1;
      } else {
        list_data(y) = NULL;
        list_nmax(y) = 0;
      }
      return y;
    }
    default:
      y = cgetg_copy_avma(x, &lx, AVMA);
      if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
      for (; i<lx; i++) gel(y,i) = gcopy_av0_canon(gel(x,i), AVMA);
  }
return y;
}

/* [copy_bin/bin_copy:] size (number of words) required for gcopy_av0(x) */
static long
taille0(GEN x)
{
  long i,n,lx, tx = typ(x);
  switch(tx)
  { /* non recursive types */
    case t_INT: return signe(x)? lgefint(x): 0;
    case t_REAL:
    case t_STR:
    case t_VECSMALL: return lg(x);

    /* one more special case */
    case t_LIST:
    {
      GEN L = list_data(x);
      return L? 3 + taille0(L): 3;
    }
    default:
      n = lx = lg(x);
      for (i=lontyp[tx]; i<lx; i++) n += taille0(gel(x,i));
      return n;
  }
}

/* [copy_bin/bin_copy:] size (number of words) required for gcopy_av0(x) */
static long
taille0_nolist(GEN x)
{
  long i,n,lx, tx = typ(x);
  switch(tx)
  { /* non recursive types */
    case t_INT:
      lx = lgefint(x);
      return lx == 2? 0: lx;
    case t_REAL:
    case t_STR:
    case t_VECSMALL:
    case t_LIST:
      return lg(x);
    default:
      n = lx = lg(x);
      for (i=lontyp[tx]; i<lx; i++) n += taille0_nolist(gel(x,i));
      return n;
  }
}

long
gsizeword(GEN x)
{
  long i,n,lx, tx = typ(x);
  switch(tx)
  { /* non recursive types */
    case t_INT: return lgefint(x);
    case t_REAL:
    case t_STR:
    case t_VECSMALL: return lg(x);

    case t_LIST: return 3;
    default:
      n = lx = lg(x);
      for (i=lontyp[tx]; i<lx; i++) n += gsizeword(gel(x,i));
      return n;
  }
}
long
gsizebyte(GEN x) { return gsizeword(x) * sizeof(long); }

/* return a clone of x structured as a gcopy */
GENbin*
copy_bin(GEN x)
{
  long t = taille0_nolist(x);
  GENbin *p = (GENbin*)pari_malloc(sizeof(GENbin) + t*sizeof(long));
  pari_sp AVMA = (pari_sp)(GENbinbase(p) + t);
  p->canon = 0;
  p->len = t;
  p->x   = gcopy_av0(x, &AVMA);
  p->base= (GEN)AVMA; return p;
}

/* same, writing t_INT in canonical native form */
GENbin*
copy_bin_canon(GEN x)
{
  long t = taille0(x);
  GENbin *p = (GENbin*)pari_malloc(sizeof(GENbin) + t*sizeof(long));
  pari_sp AVMA = (pari_sp)(GENbinbase(p) + t);
  p->canon = 1;
  p->len = t;
  p->x   = gcopy_av0_canon(x, &AVMA);
  p->base= (GEN)AVMA; return p;
}

GEN
gclone(GEN x)
{
  long i,lx,tx = typ(x), t = gsizeword(x);
  GEN y = newblock(t);
  switch(tx)
  { /* non recursive types */
    case t_INT:
      lx = lgefint(x);
      y[0] = evaltyp(t_INT)|evallg(lx);
      for (i=1; i<lx; i++) y[i] = x[i];
      break;
    case t_REAL:
    case t_STR:
    case t_VECSMALL:
      lx = lg(x);
      for (i=0; i<lx; i++) y[i] = x[i];
      break;

    /* one more special case */
    case t_LIST:
      y[0] = evaltyp(t_LIST)|_evallg(3);
      listassign(x, y);
      break;
    default: {
      pari_sp AVMA = (pari_sp)(y + t);
      lx = lg(x);
      y[0] = x[0];
      if (lontyp[tx] == 1) i = 1; else { y[1] = x[1]; i = 2; }
      for (; i<lx; i++) gel(y,i) = gcopy_avma(gel(x,i), &AVMA);
    }
  }
  setisclone(y); return y;
}

void
shiftaddress(GEN x, long dec)
{
  long i, lx, tx = typ(x);
  if (is_recursive_t(tx) && tx != t_LIST)
  {
    lx = lg(x);
    for (i=lontyp[tx]; i<lx; i++) {
      if (!x[i]) gel(x,i) = gen_0;
      else
      {
        x[i] += dec;
        shiftaddress(gel(x,i), dec);
      }
    }
  }
}

void
shiftaddress_canon(GEN x, long dec)
{
  long i, lx, tx = typ(x);
  switch(tx)
  { /* non recursive types */
    case t_INT: {
      GEN y;
      lx = lgefint(x); if (lx <= 3) return;
      y = x + 2;
      x = int_MSW(x);  if (x == y) return;
      while (x > y) { lswap(*x, *y); x = int_precW(x); y++; }
      break;
    }
    case t_REAL:
    case t_STR:
    case t_VECSMALL:
      break;

    /* one more special case */
    case t_LIST: {
      GEN Lx = list_data(x);
      if (Lx) {
        pari_sp av = avma;
        GEN L = (GEN)((long)Lx+dec);
        shiftaddress_canon(L, dec);
        list_data(x) = list_internal_copy(L, lg(L)); avma = av;
      }
    }
    default:
      lx = lg(x);
      for (i=lontyp[tx]; i<lx; i++) {
        if (!x[i]) gel(x,i) = gen_0;
        else
        {
          x[i] += dec;
          shiftaddress_canon(gel(x,i), dec);
        }
      }
  }
}

/*******************************************************************/
/*                                                                 */
/*                         STACK MANAGEMENT                        */
/*                                                                 */
/*******************************************************************/
INLINE void
dec_gerepile(pari_sp *x, pari_sp av0, pari_sp av, pari_sp tetpil, size_t dec)
{
  if (*x < av && *x >= av0)
  { /* update address if in stack */
    if (*x < tetpil) *x += dec;
    else pari_err(talker, "significant pointers lost in gerepile! (please report)");
  }
}

void
gerepileallsp(pari_sp av, pari_sp tetpil, int n, ...)
{
  const pari_sp av0 = avma;
  const size_t dec = av-tetpil;
  int i;
  va_list a; va_start(a, n);
  (void)gerepile(av,tetpil,NULL);
  for (i=0; i<n; i++) dec_gerepile((pari_sp*)va_arg(a,GEN*), av0,av,tetpil,dec);
}

/* Takes an array of pointers to GENs, of length n.
 * Cleans up the stack between av and tetpil, updating those GENs. */
void
gerepilemanysp(pari_sp av, pari_sp tetpil, GEN* gptr[], int n)
{
  const pari_sp av0 = avma;
  const size_t dec = av-tetpil;
  int i;
  (void)gerepile(av,tetpil,NULL);
  for (i=0; i<n; i++) dec_gerepile((pari_sp*)gptr[i], av0, av, tetpil, dec);
}

/* Takes an array of GENs (cast to longs), of length n.
 * Cleans up the stack between av and tetpil, updating those GENs. */
void
gerepilecoeffssp(pari_sp av, pari_sp tetpil, long *g, int n)
{
  const pari_sp av0 = avma;
  const size_t dec = av-tetpil;
  int i;
  (void)gerepile(av,tetpil,NULL);
  for (i=0; i<n; i++,g++) dec_gerepile((pari_sp*)g, av0, av, tetpil, dec);
}

static int
dochk_gerepileupto(GEN av, GEN x)
{
  long i,lx,tx;
  if (!isonstack(x)) return 1;
  if (x > av)
  {
    pari_warn(warner,"bad object %Ps",x);
    return 0;
  }
  tx = typ(x);
  if (! is_recursive_t(tx)) return 1;

  lx = lg(x);
  for (i=lontyp[tx]; i<lx; i++)
    if (!dochk_gerepileupto(av, gel(x,i)))
    {
      pari_warn(warner,"bad component %ld in object %Ps",i,x);
      return 0;
    }
  return 1;
}
/* check that x and all its components are out of stack, or have been
 * created after av */
int
chk_gerepileupto(GEN x) { return dochk_gerepileupto(x, x); }

/* print stack between avma & av */
void
dbg_gerepile(pari_sp av)
{
  GEN x = (GEN)avma;
  while (x < (GEN)av)
  {
    const long tx = typ(x), lx = lg(x);
    GEN a;

    pari_printf(" [%ld] %Ps:", x - (GEN)avma, x);
    if (! is_recursive_t(tx)) { pari_putc('\n'); x += lx; continue; }
    a = x + lontyp[tx]; x += lx;
    for (  ; a < x; a++) pari_printf("  %Ps,", *a);
    pari_printf("\n");
  }
}
void
dbg_gerepileupto(GEN q)
{
  fprintferr("%Ps:\n", q);
  dbg_gerepile((pari_sp) (q+lg(q)));
}

GEN
gerepile(pari_sp av, pari_sp tetpil, GEN q)
{
  const size_t dec = av - tetpil;
  const pari_sp av0 = avma;
  GEN x, a;

  if (dec == 0) return q;
  if ((long)dec < 0) pari_err(talker,"lbot>ltop in gerepile");

  /* dec_gerepile(&q, av0, av, tetpil, dec), saving 1 comparison */
  if (q >= (GEN)av0 && q < (GEN)tetpil)
    q = (GEN) (((pari_sp)q) + dec);

  for (x = (GEN)av, a = (GEN)tetpil; a > (GEN)av0; ) *--x = *--a;
  avma = (pari_sp)x;
  while (x < (GEN)av)
  {
    const long tx = typ(x), lx = lg(x);

    if (! is_recursive_t(tx)) { x += lx; continue; }
    a = x + lontyp[tx]; x += lx;
    for (  ; a < x; a++) dec_gerepile((pari_sp*)a, av0, av, tetpil, dec);
  }
  return q;
}

long
allocatemoremem(size_t newsize)
{
  size_t s, old = top - bot;
  if (!newsize) newsize = old << 1;
  pari_init_stack(newsize, old);
  s = top - bot;
  pari_warn(warner,"new stack size = %lu (%.3f Mbytes)", s, s/1048576.);
  return s;
}

void
fill_stack(void)
{
  GEN x = ((GEN)bot);
  while (x < (GEN)avma) *x++ = 0xfefefefeUL;
}

void
debug_stack(void)
{
  GEN z;
  fprintferr("bot=0x%lx\ttop=0x%lx\n", bot, top);
  for (z = (GEN)top; z >= (GEN)avma; z--)
    fprintferr("%p:\t0x%lx\t%lu\n",z,*z,*z);
}

long
getstack(void) { return top-avma; }

/*******************************************************************/
/*                                                                 */
/*                               TIMER                             */
/*                                                                 */
/*******************************************************************/

#if !defined(USE_GETRUSAGE) && !defined(USE_FTIME)
static long
_get_time(pari_timer *T, long Ticks, long TickPerSecond)
{
  long s  = Ticks / TickPerSecond;
  long us = (long) ((Ticks % TickPerSecond) * (1000000. / TickPerSecond));
  long delay = 1000 * (s - T->s) + (us - T->us) / 1000;
  T->us = us;
  T->s  = s; return delay;
}
#endif

#ifdef USE_TIMES

# include <sys/times.h>
# include <sys/time.h>
# include <time.h>
long
TIMER(pari_timer *T)
{
  struct tms t; times(&t);
  return _get_time(T, t.tms_utime,
#ifdef _SC_CLK_TCK
                      sysconf(_SC_CLK_TCK)
#else
                      CLK_TCK
#endif
  );
}
#elif defined(USE_GETRUSAGE)

# include <sys/time.h>
# include <sys/resource.h>
long
TIMER(pari_timer *T)
{
  struct rusage r;
  struct timeval t;
  long delay;

  getrusage(RUSAGE_SELF,&r); t = r.ru_utime;
  delay = 1000 * (t.tv_sec - T->s) + (t.tv_usec - T->us) / 1000;
  T->us = t.tv_usec;
  T->s  = t.tv_sec; return delay;
}
#elif defined(USE_FTIME)

# include <sys/timeb.h>
long
TIMER(pari_timer *T)
{
  struct timeb t;
  long delay;

  ftime(&t);
  delay = 1000 * (t.time - T->s) + (t.millitm - T->us / 1000);
  T->us = t.millitm * 1000;
  T->s  = t.time; return delay;
}
#elif defined(WINCE)
long
TIMER(pari_timer *T)
{
  return _get_time(T, GetTickCount(), 1000);
}
#else

# include <time.h>
# ifndef CLOCKS_PER_SEC
#   define CLOCKS_PER_SEC 1000000 /* may be false on YOUR system */
# endif
long
TIMER(pari_timer *T)
{
  return _get_time(T, clock(), CLOCKS_PER_SEC);
}
#endif
void
TIMERstart(pari_timer *T) { T->s = 0; T->us = 0; (void)TIMER(T); }
long
TIMERread(pari_timer *T) {
  long s = T->s, us = T->us, delay = TIMER(T);
  T->s = s; T->us = us; return delay;
}

long
timer(void)   { static THREAD pari_timer T; return TIMER(&T);}
long
timer2(void)  { static THREAD pari_timer T; return TIMER(&T);}

void
msgTIMER(pari_timer *T, const char *format, ...)
{
  va_list args;
  PariOUT *out = pariOut; pariOut = pariErr;

  pari_puts("Time "); va_start(args, format);
  pari_vprintf(format,args); va_end(args);
  pari_printf(": %ld\n", TIMER(T)); pari_flush();
  pariOut = out;
}

void
msgtimer(const char *format, ...)
{
  va_list args;
  PariOUT *out = pariOut; pariOut = pariErr;

  pari_puts("Time "); va_start(args, format);
  pari_vprintf(format,args); va_end(args);
  pari_printf(": %ld\n", timer2()); pari_flush();
  pariOut = out;
}

long
gettime(void) { return timer(); }

/*******************************************************************/
/*                                                                 */
/*                   FUNCTIONS KNOWN TO THE ANALYZER               */
/*                                                                 */
/*******************************************************************/
GEN
pari_version(void) {
  GEN v = cgetg(4, t_VEC);
  const ulong mask = (1UL<<PARI_VERSION_SHIFT) - 1;
  ulong major, minor, patch, n = PARI_VERSION_CODE;
  patch = n & mask; n >>= PARI_VERSION_SHIFT;
  minor = n & mask; n >>= PARI_VERSION_SHIFT;
  major = n;
  gel(v,1) = utoi(major);
  gel(v,2) = utoi(minor);
  gel(v,3) = utoi(patch); return v;
}

/* List of GP functions:
 * ---------------------
 * Format (struct entree) :
 *   char *name    : name (under GP).
 *   ulong valence : used to form arg list, now often handled by code.
 *   void *value   : For PREDEFINED FUNCTIONS: C function to call.
 *                   For USER FUNCTIONS: pointer to defining data (block) =
 *                    entree*: NULL, list of entree (arguments), NULL
 *                    char*  : function text
 *   long menu     : which help section do we belong to (See below).
 *   char *code    : argument list (See below).
 *   entree *next  : next entree (init to NULL, used in hashing code).
 *   char *help    : short help text (init to NULL).
 *   void *args    : For USER FUNCTIONS: default arguments (NULL terminated).
 *                   For VARIABLES: (esp. loop indexes): push_val history.
 *                   (while processing a loop, ep->value may not be a block)
 * menu:
 * -----
 *  1: Standard monadic or dyadic OPERATORS
 *  2: CONVERSIONS and similar elementary functions
 *  3: TRANSCENDENTAL functions
 *  4: NUMBER THEORETICAL functions
 *  5: Functions related to ELLIPTIC CURVES
 *  6: Functions related to general NUMBER FIELDS
 *  7: POLYNOMIALS and power series
 *  8: Vectors, matrices, LINEAR ALGEBRA and sets
 *  9: SUMS, products, integrals and similar functions
 *  10: GRAPHIC functions
 *  11: PROGRAMMING under GP
 *
 * code: describe function prototype. NULL = use valence instead.
 * -----
 * Arguments:
 *  I  closure whose value is ignored, like in for() loop
 *  E  closure whose value is used, like in sum() loop
 *  G  GEN
 *  L  long
 *  S  symbol (i.e GP function name) as a entree *
 *  V  lexical variable
 *  C  lexical context
 *  n  variable number
 *  &  *GEN
 *  f  Fake *long (function requires it, but we don't use the resulting long)
 *  p  real precision (prec for the C library)
 *  P  series precision (precdl for the C library)
 *  r  raw input (treated as a string without quotes).
 *     Quoted args are copied as strings. Stops at first unquoted ')' or ','.
 *     Special chars can be quoted using '\'.  Ex : aa"b\n)"c => "aab\n)c".
 *  s  expanded string. Example: Pi"x"2 yields "3.142x2".
 *     The unquoted components can be of any pari type (converted according to
 *     the current output format)
 *  s* any number of strings (see s)
 *  M  Mnemonic or a flag (converted to a long); description follows
 *         after \n at the end of the argument description
 *  D  Has a default value. Format is "Dvalue,type," (the ending comma is
 *     mandatory). Ex: D0,L, (arg is long, 0 by default).
 *     Special syntax:
 *       if type = G, &, I or V:  D[G&IV] all send NULL.
 *       if type = n: Dn sends -1.
 *
 *     The user-given args are read first, then completed by the defaults
 *
 * Return type (first char or immediately after 'x'): GEN by default, otherwise
 *  l Return long
 *  i Return int
 *  v Return void
 *
 * Syntax requirements:
 *  = Separator '=' required.
 *
 * Origin:
 *  x Installed foreign function. Put the ep of the function as the
 *       first argument, fill the rest with PARI arguments,
 *       then call installedHandler with these arguments.
 *       Should be the first char in the code.
 *
 ****************************************************************************
 * If new codes are added, change identifier and skipidentifier.
 *
 * Currently the following functions have no code word, but a valence code.
 * 'O' 50, 'if' 80, 'until' 82, 'while' 81, 'global' 88,
 * Valences:
 * 0  for functions without mandatory args.
 * 1  for functions with mandatory args.
 * 50 'O'
 * 80 'if'
 * 82 'until'
 * 81 'while'
 * 88 'global'
 */
#include "init.h"
