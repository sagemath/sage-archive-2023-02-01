/* $Id: paripriv.h,v 1.106.2.2 2007/03/29 08:58:00 kb Exp $

Copyright (C) 2004  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

BEGINEXTERN

/* for qsort */
typedef int (*QSCOMP)(const void *, const void *);

/* swap */
#define lswap(x,y) {long _z=x; x=y; y=_z;}
#define pswap(x,y) {GEN *_z=x; x=y; y=_z;}
#define swap(x,y)  {GEN  _z=x; x=y; y=_z;}
#define dswap(x,y) { double _t=x; x=y; y=_t; }
#define pdswap(x,y) { double* _t=x; x=y; y=_t; }
#define swapspec(x,y, nx,ny) {swap(x,y); lswap(nx,ny);}

/* */
#define both_odd(x,y) ((x)&(y)&1)

/* unused */
GEN ellheightoo(GEN e, GEN z, long prec);
void ellprint(GEN e);
GEN mpinv(GEN b);

/* generic */
GEN arith_proto(long f(GEN), GEN x, int do_error);
GEN arith_proto2(long f(GEN,GEN), GEN x, GEN n);
GEN arith_proto2gs(long f(GEN,long), GEN x, long y);
GEN gassoc_proto(GEN f(GEN,GEN),GEN,GEN);
GEN garith_proto(GEN f(GEN), GEN x, int do_error);
GEN garith_proto2gs(GEN f(GEN,long), GEN x, long y);
GEN trans_fix_arg(long *prec, GEN *s0, GEN *sig, pari_sp *av, GEN *res);
GEN transc(GEN (*f) (GEN, long), GEN x, long prec);
GEN sort_factor(GEN y, int (*cmp)(GEN,GEN));
GEN sort_factor_gen(GEN y, int (*cmp)(GEN,GEN));
GEN sort_factor_gen_aux(GEN y, void *data, int (*cmp)(void *,GEN,GEN));
GEN sort_vecpol(GEN a, int (*cmp)(GEN,GEN));

/* loops */
GEN incloop(GEN a);
GEN resetloop(GEN a, GEN b);
GEN setloop(GEN a);

/* multiprecision */
GEN   icopy_spec(GEN x, long nx);
GEN   addrex01(GEN x);
GEN   addumului(ulong a, ulong b, GEN Y);
void  affr_fixlg(GEN z, GEN y);
GEN   cxnorm(GEN x);
int   lgcdii(ulong* d, ulong* d1, ulong* u, ulong* u1, ulong* v, ulong* v1, ulong vmax);
ulong rgcduu(ulong d, ulong d1, ulong vmax, ulong* u, ulong* u1, ulong* v, ulong* v1, long *s);
ulong xgcduu(ulong d, ulong d1, int f, ulong* v, ulong* v1, long *s);
GEN   quadnorm(GEN x);
ulong xxgcduu(ulong d, ulong d1, int f, ulong* u, ulong* u1, ulong* v, ulong* v1, long *s);
GEN   divgsns(GEN x, long i);
GEN   divrsns(GEN x, long i);
GEN   init_remiimul(GEN M);
GEN   ishiftr_lg(GEN x, long lx, long n);
GEN   logagmcx(GEN q, long prec);
GEN   muliispec(GEN x, GEN y, long nx, long ny);
GEN   padic_to_Fp(GEN x, GEN Y);
ulong padic_to_Fl(GEN x, ulong p);
GEN   red_montgomery(GEN T, GEN N, ulong inv);
GEN   remiimul(GEN x, GEN sy);
GEN   sqrispec(GEN x, long nx);
GEN   subrex01(GEN x);
GEN   mulcxI(GEN x);
GEN   mulcxmI(GEN x);

int approx_0(GEN x, GEN y);
GEN bernfrac_using_zeta(long n);
int OK_bern(long nb, long prec);

/* FIXME: adapt/use mpn_[lr]shift instead */
#define shift_left(z2,z1,imin,imax,f, sh) {\
  register const ulong _m = BITS_IN_LONG - (sh);\
  shift_left2((z2),(z1),(imin),(imax),(f),(sh),(_m)); }

#define shift_right(z2,z1,imin,imax,f, sh) {\
  register const ulong _m = BITS_IN_LONG - (sh);\
  shift_right2((z2),(z1),(imin),(imax),(f),(sh),(_m)); }

/* powers */
#define sqrs(b) mulss((b),(b))
#define sqru(b) muluu((b),(b))
GEN    rpowuu(ulong a, ulong n, long prec);
GEN    powrshalf(GEN x, long s);
GEN    powrfrac(GEN x, long n, long d);
ulong  u_pow10(int n);

/* floats */
double dabs(double s, double t);
long   dblexpo(double x);
ulong  dblmantissa(double x);
void   dcxlog(double s, double t, double *a, double *b);
double dnorm(double s, double t);
double dbllog2(GEN z);
ulong  usqrtsafe(ulong a);

/* "abs" routines for t_REAL ( disregard sign ) */
int absrnz_egal1(GEN x);
int absrnz_egal2n(GEN x);
GEN exp1r_abs(GEN x);
GEN logagmr_abs(GEN q);
GEN logr_abs(GEN x);
GEN sqrtr_abs(GEN x);

/* hnf */
GEN gauss_triangle_i(GEN A, GEN B,GEN t);
GEN hnfadd(GEN m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,GEN extramat,GEN extraC);
GEN hnfadd_i(GEN m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,GEN extramat,GEN extraC);
GEN hnf_gauss(GEN A, GEN B);
GEN hnf_invimage(GEN A, GEN b);
GEN hnfmerge_get_1(GEN A, GEN B);
GEN hnfperm_i(GEN A, GEN *ptU, GEN *ptperm);
GEN hnfspec_i(long** m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,long k0);
GEN hnfspec(long** m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,long k0);
GEN mathnfspec(GEN x, GEN *ptperm, GEN *ptdep, GEN *ptB, GEN *ptC);

/* LLL */
GEN lllint_fp_ip(GEN x, long D);
GEN lllfp_marked(long *M, GEN x, long D, long flag, long prec, int gram);
GEN lllint_marked(long *M, GEN x, long D, int g, GEN *h, GEN *f, GEN *B);
GEN LLL_check_progress(GEN Bnorm, long n0, GEN m, int final, long *ti_LLL);
GEN extendedgcd(GEN A);
GEN sqred1intern(GEN a);
GEN lllgramintern(GEN x, long alpha, long flag, long prec);
GEN lllintern(GEN x, long D, long flag, long prec);

/* miscellaneous linear algebra */
GEN  diagonal_i(GEN x);
GEN  F2V_red_ip(GEN v);
GEN  gauss_realimag(GEN x, GEN y);
GEN  imagecomplspec(GEN x, long *nlze);
GEN  R_from_QR(GEN x, long prec);
void rowselect_p(GEN A, GEN B, GEN p, long init);
GEN  split_realimag(GEN x, long r1, long r2);
GEN  sqred1_from_QR(GEN x, long prec);
GEN  supnorm(GEN L, long prec);
GEN  znstar_hnf_elts(GEN Z, GEN H);
GEN  ZV_lincomb(GEN u, GEN v, GEN X, GEN Y);
GEN  vec_setconst(GEN v, GEN x);

GEN  GS_norms(GEN B, long prec);
GEN  dim1proj(GEN prh);
GEN  detcyc(GEN cyc, long *L);
GEN  close_modinvertible(GEN x, GEN y);
GEN  colreducemodHNF(GEN x, GEN y, GEN *Q);
GEN  col_to_ff(GEN x, long v);

/* famat */
GEN factorback_i(GEN fa, GEN e, GEN nf, int red);
GEN factorbackprime(GEN nf, GEN L, GEN e);
GEN famat_inv(GEN f);
GEN famat_makecoprime(GEN nf, GEN g, GEN e, GEN pr, GEN prk, GEN EX);
GEN famat_mul(GEN f, GEN g);
GEN famat_pow(GEN f, GEN n);
GEN famat_reduce(GEN fa);
GEN famat_to_arch(GEN nf, GEN fa, long prec);
GEN famat_to_nf_modideal_coprime(GEN nf, GEN g, GEN e, GEN id, GEN EX);
GEN famat_to_nf_modidele(GEN nf, GEN g, GEN e, GEN bid);
GEN merge_factor_i(GEN f, GEN g);
GEN to_famat_all(GEN x, GEN y);
GEN to_famat(GEN g, GEN e);
GEN trivfact(void);

/* integer factorization / discrete log */
int   BSW_isprime(GEN x);
int   BSW_isprime_small(GEN x);
GEN   coprime_part(GEN x, GEN f);
GEN   Z_factor_limit(GEN n, GEN limit);
GEN   Fp_PHlog(GEN a, GEN g, GEN p, GEN ord);
GEN   Fp_shanks(GEN x,GEN g0,GEN p, GEN q);
ulong is_kth_power(GEN x, ulong p, GEN *pt, byteptr d);
long  ifac_decomp_break(GEN n, long (*B)(GEN,GEN,GEN,GEN), GEN s, long hint);
long  ifac_moebius(GEN n, long hint);
long  ifac_issquarefree(GEN n, long hint);
long  ifac_omega(GEN n, long hint);
long  ifac_bigomega(GEN n, long hint);
GEN   ifac_totient(GEN n, long hint);
GEN   ifac_numdiv(GEN n, long hint);
GEN   ifac_sumdiv(GEN n, long hint);
GEN   ifac_sumdivk(GEN n, long k, long hint);
int   miller(GEN n, long k);
GEN   mpqs(GEN N);
ulong ugcd(ulong x, ulong y);
long  Z_lvalrem_stop(GEN n, ulong p, int *stop);

/* quadratic forms, quadratic numbers */
long cornacchia(GEN d, GEN p, GEN *px, GEN *py);
long cornacchia2(GEN d, GEN p, GEN *px, GEN *py);
GEN  primeform_u(GEN x, ulong p);
GEN  qf_disc(GEN x);
void qfb_comp(GEN z,GEN x,GEN y);
GEN  qfr_to_qfr5(GEN x, long prec);
GEN  qfr3_comp(GEN x, GEN y, GEN D, GEN isqrtD);
GEN  qfr3_pow(GEN x, GEN n, GEN D, GEN isqrtD);
GEN  qfr3_red(GEN x, GEN D, GEN isqrtD);
GEN  qfr3_rho(GEN x, GEN D, GEN isqrtD);
GEN  qfr3_to_qfr(GEN x, GEN z);
GEN  qfr5_dist(GEN e, GEN d, long prec);
GEN  qfr5_comp(GEN x, GEN y, GEN D, GEN sqrtD, GEN isqrtD);
GEN  qfr5_pow(GEN x, GEN n, GEN D, GEN sqrtD, GEN isqrtD);
GEN  qfr5_red(GEN x, GEN D, GEN sqrtD, GEN isqrtD);
GEN  qfr5_rho(GEN x, GEN D, GEN sqrtD, GEN isqrtD);
GEN  qfr_pow(GEN x, GEN n);
GEN  qfr_unit(GEN x);
GEN  qfi_unit(GEN x);
GEN  quad_polmod_conj(GEN x, GEN y);
GEN  quad_polmod_norm(GEN x, GEN y);

/* Polynomials */
/* a) Arithmetic/conversions */
GEN  addmulXn(GEN x, GEN y, long d);
GEN  addshiftpol(GEN x, GEN y, long d);
GEN  lift_if_rational(GEN x);
GEN  monomial(GEN a, long degpol, long v);
GEN  monomialcopy(GEN a, long degpol, long v);
GEN  mulmat_pol(GEN A, GEN x);
long polegal_spec(GEN x, GEN y);
GEN  polrecip_i(GEN x);
GEN  pol_to_monic(GEN pol, GEN *lead);
GEN  revpol(GEN x);
GEN  ser2pol_i(GEN x, long lx);
GEN  ser2rfrac_i(GEN x);
GEN  shiftpol_i(GEN x, long v);
GEN  swap_vars(GEN b0, long v);
GEN  translate_pol(GEN P, GEN c);

/* b) Modular */
GEN  bezout_lift_fact(GEN T, GEN Tmod, GEN p, long e);
GEN  caractducos(GEN p, GEN x, long v);
GEN  FpXQX_from_Kronecker(GEN z, GEN pol, GEN p);
GEN  FpX_quad_root(GEN x, GEN p, int unknown);
long FpX_split_Berlekamp(GEN *t, GEN pp);
GEN  FqX_split_all(GEN z, GEN T, GEN p);
long FqX_split_by_degree(GEN *pz, GEN u, GEN q, GEN T, GEN p);
long FqX_split_deg1(GEN *pz, GEN u, GEN q, GEN T, GEN p);
GEN  FqX_split_roots(GEN z, GEN T, GEN p, GEN pol);
GEN  polratlift(GEN P, GEN mod, GEN amax, GEN bmax, GEN denom);
GEN  polsym_gen(GEN P, GEN y0, long n, GEN T, GEN N);
GEN  ZX_caract_sqf(GEN A, GEN B, long *lambda, long v);
GEN  ZX_disc_all(GEN,ulong);
long ZX_get_prec(GEN x);
GEN  ZX_resultant_all(GEN A, GEN B, GEN dB, ulong bound);
/*
Coment out the following two lines, as they cause problems on
Solaris. I was originally going to do this on all platforms, not just
Solaris, but I noticed there was already an altered partpriv-osx.h,
so this modified version will only be copied over if using Solaris.
GEN  ZY_ZXY_resultant_all(GEN A, GEN B0, long *lambda, GEN *LPRS);
GEN  ZY_ZXY_rnfequation(GEN A, GEN B0, long *lambda); */
GEN  RgXQ_u_pow(GEN x, ulong n, GEN T);
GEN  RgX_gcd_simple(GEN x, GEN y);
GEN  RgX_extgcd_simple(GEN a, GEN b, GEN *pu, GEN *pv);
GEN  RgXQ_minpoly_naive(GEN y, GEN P);
GEN lift_intern0(GEN x,long v);
#define lift_intern(x) (lift_intern0((x),-1))

/* b') Chinese Remainder Theorem */
GEN ZM_init_CRT(GEN Hp, ulong p);
int ZM_incremental_CRT(GEN H, GEN Hp, GEN q, GEN qp, ulong p);

/* c) factorization */
double cauchy_bound(GEN p);
GEN chk_factors_get(GEN lt, GEN famod, GEN c, GEN T, GEN N);
int cmbf_precs(GEN q, GEN A, GEN B, long *a, long *b, GEN *qa, GEN *qb);
GEN ZX_DDF(GEN x, long hint);
GEN fact_from_DDF(GEN fa, GEN e, long n);
GEN initgaloisborne(GEN T, GEN dn, long prec, GEN *pL, GEN *pprep, GEN *pdis);
GEN logmax_modulus_bound(GEN p);
GEN polint_i(GEN xa, GEN ya, GEN x, long n, GEN *ptdy);
GEN quicktrace(GEN x, GEN sym);
GEN roots_to_pol_intern(GEN L, GEN a, long v, int plus);
GEN roots_to_pol_r1r2(GEN a, long r1, long v);
GEN special_pivot(GEN x);
GEN vandermondeinversemod(GEN L, GEN T, GEN den, GEN mod);
GEN ZX_monic_factorpadic(GEN f, GEN p, long prec);

#include "parinf.h"

/* Allocation / gerepile */
void   debug_stack(void);
void   fill_stack(void);
void   init_dalloc();
double *dalloc(size_t n);
void   gerepilecoeffs2(pari_sp av, GEN x, int n, GEN y, int o);
void   minim_alloc(long n, double ***q, GEN *x, double **y,  double **z, double **v);
int    pop_entree_bloc(entree *ep, long loc);
int    pop_val_if_newer(entree *ep, long loc);
void   gclone_refc(GEN x);
void   free_ep_args(entree *ep);

/* naive grow-arrays */
typedef struct {
  void **v;
  long len; /* len cells allocated */
  long n; /* first n cells occupied */
} __pari_growarray;
typedef __pari_growarray growarray[1];

growarray *pari_get_modules();
growarray *pari_get_oldmodules();
void    grow_append(growarray A, void *e);
void    grow_copy(growarray A, growarray B);
void    grow_init(growarray A);
void    grow_kill(growarray A);

/* Interfaces (GP, etc.) */
void  errcontext(char *msg, char *s, char *entry);
GEN   geni(void);
void* get_stack(double fraction, long min);
GEN   gpreadseq(char *c, int strict);
void  initout(int initerr);
void  init80col(long n);
char* itostr(GEN x, int minus);
void  kill_from_hashlist(entree *ep, long n);
void  member_err(char *s);
int   pari_kernel_init(void);
int   pari_last_was_newline(void);
void  pari_set_last_newline(int last);
void  print_functions_hash(const char *s);
void  print_user_fun(entree *ep);
void  print_user_member(entree *ep);
void  print_all_user_fun(void);
void  print_all_user_member(void);
void  pop_val(entree *ep);
void  push_val(entree *ep, GEN a);
GEN   readbin(const char *name, FILE *f, int *vector);
void  recover(int flag);
int   term_height(void);
int   term_width(void);
void  var_make_safe();
void  whatnow_new_syntax(char *f, long n);

/* defaults */
#define is_default(s) setdefault((s),"",d_EXISTS) == gen_1
enum { d_SILENT, d_ACKNOWLEDGE, d_INITRC, d_RETURN, d_EXISTS };
extern ulong precreal;

GEN sd_TeXstyle(const char *v, long flag);
GEN sd_colors(char *v, long flag);
GEN sd_compatible(const char *v, long flag);
GEN sd_datadir(char *v, long flag);
GEN sd_debug(const char *v, long flag);
GEN sd_debugfiles(const char *v, long flag);
GEN sd_debugmem(const char *v, long flag);
GEN sd_echo(const char *v, long flag);
GEN sd_factor_add_primes(char *v, long flag);
GEN sd_filename(const char *v, long flag, char *s, char **f);
GEN sd_format(const char *v, long flag);
GEN sd_help(char *v, long flag);
GEN sd_histsize(const char *v, long flag);
GEN sd_lines(const char *v, long flag);
GEN sd_log(const char *v, long flag);
GEN sd_logfile(const char *v, long flag);
GEN sd_new_galois_format(char *v, long flag);
GEN sd_output(const char *v, long flag);
GEN sd_parisize(const char *v, long flag);
GEN sd_path(char *v, long flag);
GEN sd_prettyprinter(char *v, long flag);
GEN sd_primelimit(const char *v, long flag);
GEN sd_prompt(const char *v, long flag);
GEN sd_prompt_cont(const char *v, long flag);
GEN sd_psfile(const char *v, long flag);
GEN sd_realprecision(const char *v, long flag);
GEN sd_rl(const char *v, long flag);
GEN sd_secure(const char *v, long flag);
GEN sd_seriesprecision(const char *v, long flag);
GEN sd_simplify(const char *v, long flag);
GEN sd_strictmatch(const char *v, long flag);
GEN sd_timer(const char *v, long flag);
GEN setdefault(const char *s, const char *v, long flag);

/* history */
typedef struct {
  GEN *res;    /* array of previous results, FIFO */
  size_t size; /* # res */
  ulong total; /* # of results computed since big bang */
} gp_hist;

/* prettyprinter */
typedef struct {
  pariFILE *file;
  char *cmd;
} gp_pp;

/* path */
typedef struct {
  char *PATH;
  char **dirs;
} gp_path;

/* for output */
typedef struct {
  char format; /* e,f,g */
  long fieldw; /* 0 (ignored) or field width */
  long sigd;   /* -1 (all) or number of significant digits printed */
  int sp;      /* 0 = suppress whitespace from output */
  int prettyp; /* output style: raw, prettyprint, etc */
  int TeXstyle;
} pariout_t;

char *GENtostr0(GEN x, pariout_t *T, void (*do_out)(GEN, pariout_t*));
void lim_lines_output(GEN z, pariout_t *fmt, long n, long max);
void gen_output(GEN x, pariout_t *T);

/* GP_DATA */
typedef struct {
  ulong primelimit;
  jmp_buf env;
  gp_hist *hist;
  gp_pp *pp;
  gp_path *path;
  pariout_t *fmt;
  ulong flags, lim_lines;
  char *help, *prompt, *prompt_cont;
  pari_timer *T;
} gp_data;
  /* GP_DATA->flags */
/*

Coment out the following two enum, as it causes problems on
Solaris. I was originally going to do this on all platforms, not just
Solaris, but I noticed there was already an altered partpriv-osx.h,
so this modified version will only be copied over if using Solaris.
enum { QUIET=1, TEST=2, SIMPLIFY=4, CHRONO=8, ECHO=16, STRICTMATCH=32,
       USE_READLINE=64, SECURE=128, EMACS=256, TEXMACS=512};
*/

extern gp_data *GP_DATA;

typedef struct Buffer {
  char *buf;
  ulong len;
  jmp_buf env;
} Buffer;
Buffer *new_buffer(void);
void delete_buffer(Buffer *b);
void fix_buffer(Buffer *b, long newlbuf);

typedef struct {
  char *s, *t, *end; /* source, target, last char read */
  int in_string, in_comment, more_input, wait_for_brace, downcase;
  Buffer *buf;
} filtre_t;
void init_filtre(filtre_t *F, Buffer *buf);
char *filtre(char *s, int flag);
void check_filtre(filtre_t *F);

gp_data *default_gp_data(void);
GEN  gp_history(gp_hist *H, long p, char *old, char *entry);
GEN  set_hist_entry(gp_hist *H, GEN x);

void delete_dirs(gp_path *p);
void gp_expand_path(gp_path *p);
const char *pari_default_path();
char *expand_prompt(char *prompt, filtre_t *F);

typedef struct input_method {
/* mandatory */
  char * (*fgets)(char *,int,FILE*);
  char * (*getline)(char**, int f, struct input_method*, filtre_t *F);
  int free; /* boolean: must we free the output of getline() ? */
/* for interactive methods */
  char *prompt, *prompt_cont;
/* for non-interactive methods */
  FILE *file;
} input_method;

int input_loop(filtre_t *F, input_method *IM);
char *file_input(char **s0, int junk, input_method *IM, filtre_t *F);

/* By files */

/* Flx.c */

GEN     FlxX_recipspec(GEN x, long l, long n, long vs);
GEN     FlxX_sub(GEN x, GEN y, ulong p);
GEN     FlxX_subspec(GEN x, GEN y, ulong p, long lx, long ly);
GEN     FlxX_swap(GEN x, long n, long ws);
GEN     FlxX_to_Kronecker(GEN P, GEN Q);
GEN     FlxX_to_Kronecker_spec(GEN P, long lp, GEN Q);
GEN     Flx_addshift(GEN x, GEN y, ulong p, long d);
GEN     Flx_addspec(GEN x, GEN y, ulong p, long lx, long ly);
GEN     Flx_even_odd_comb(GEN P, ulong u, ulong v, ulong p);
GEN     Flx_mulspec(GEN a, GEN b, ulong p, long na, long nb);
GEN     Flx_negspec(GEN x, ulong p, long l);
GEN     Flx_recipspec(GEN x, long l, long n);
GEN     Flx_sqrspec(GEN a, ulong p, long na);
GEN     Flx_subspec(GEN x, GEN y, ulong p, long lx, long ly);
GEN     FlxqX_from_Kronecker(GEN z, GEN T, ulong p);
GEN     FlxqX_invmontgomery(GEN T, GEN Q, ulong p);
GEN     FlxqX_mulspec(GEN x, GEN y, GEN T, ulong p, long lx, long ly);
GEN     FlxqX_rem_montgomery(GEN x, GEN mg, GEN T, GEN Q, ulong p);

/* Qfb.c */

GEN     powimagraw(GEN x, long n);
GEN     redimagsl2(GEN q, GEN *U);
GEN     redrealsl2(GEN V);
GEN     redrealsl2step(GEN A);
GEN     sqcompimagraw(GEN x);
GEN     sqcomprealraw(GEN x);

/* RgX.c */

GEN     RgXQ_matrix_pow(GEN y, long n, long m, GEN P);
GEN     RgX_copy_spec(GEN x, long nx);
GEN     zx_copy_spec(GEN x, long nx);

/* alglin1.c */

GEN     Flm_Fl_mul_inplace(GEN y, ulong x, ulong p);
GEN     gauss_get_col(GEN a, GEN b, GEN p, long li);
GEN     gauss_intern(GEN a, GEN b);
GEN     matid_intern(long n,GEN myun,GEN myzero);
void    vecselect_p(GEN A, GEN B, GEN p, long init, long lB);

/* alglin2.c */

GEN     Frobeniusform(GEN V, long n);
GEN     ZM_copy(GEN x);
GEN     ZV_Z_mul(GEN X, GEN c);
GEN     ZV_copy(GEN x);
void    ZV_neg_ip(GEN M);
GEN     hnf0(GEN A, int remove);
GEN     hnf_special(GEN x, long remove);
GEN     hnfall0(GEN A, long remove);

/* anal.c */

long unsigned int  parse_option_string(char *arg, char *tmplate, long flag, char **failure, char **failure_arg);
GEN     readexpr(char *t);

/* aprcl.c */

GEN     aprcl(GEN N);

/* arith1.c */

GEN     Z_chinese_coprime(GEN a, GEN b, GEN A, GEN B, GEN C);
GEN     bestappr_mod(GEN x, GEN A, GEN B);
long    hilii(GEN x, GEN y, GEN p);
long    krouu_s(ulong x, ulong y, long s);

/* arith2.c */

GEN     ibitand(GEN x, GEN y);
GEN     ibitnegimply(GEN x, GEN y);
GEN     ibitor(GEN x, GEN y);
GEN     ibitxor(GEN x, GEN y);
byteptr initprimes0(ulong maxnum, long *lenp, ulong *lastp);
long    set_optimize(long what, GEN g);

/* base1.c */

void    nfbasic_init(GEN x, long flag, GEN fa, nfbasic_t *T);
GEN     nffromhnfbasis(GEN nf, GEN x);
GEN     nftohnfbasis(GEN nf, GEN x);
GEN     polgaloisnames(long a, long b);
GEN     primitive_pol_to_monic(GEN pol, GEN *ptlead);

/* base2.c */

long    FpX_val(GEN x0, GEN t, GEN p, GEN *py);
GEN     fast_respm(GEN f, GEN g, GEN p, long M);
GEN     gen_if_principal(GEN bnf, GEN x);
int     nfissquarefree(GEN nf, GEN x);
GEN     nfreducemodpr_i(GEN x, GEN prh);
GEN     polsymmodp(GEN g, GEN p);
GEN     respm(GEN x, GEN y, GEN pm);
GEN     rnfdet0(GEN nf, GEN x, GEN y);

/* base3.c */

GEN     FpXQ_gener(GEN T, GEN p);
void    check_nfelt(GEN x, GEN *den);
GEN     ff_PHlog(GEN a, GEN g, GEN T, GEN p);
GEN     nf_PHlog(GEN nf, GEN a, GEN g, GEN pr);

/* base4.c */

GEN     arch_inv(GEN x);
GEN     arch_pow(GEN x, GEN n);
void    check_listpr(GEN x);
GEN     factor_norm(GEN x);
GEN     famat_to_nf(GEN nf, GEN f);
long    val_norm(GEN x, GEN p, long *vz);

/* base5.c */

GEN     check_and_build_nfabs(GEN rnf);
GEN     check_and_build_norms(GEN rnf);
GEN     checkrnfeq(GEN x);
GEN     hnfcenter_ip(GEN M);

/* bibli1.c */

GEN     gscal(GEN x, GEN y);
GEN     lll_scaled(long MARKED, GEN X0, long D);
GEN     lllall(GEN x, long D, int gram, long flag);
GEN     lllintpartialall(GEN m, long flag);
GEN     plindep(GEN x);
GEN     pslq(GEN x);
GEN     pslqL2(GEN x);
GEN     sqscal(GEN x);

/* bibli2.c */

long    ZV_search(GEN x, GEN y);
GEN     ZV_sort_uniq(GEN L);
GEN     gen_vecsort(GEN x, GEN k, long flag);

/* buch1.c */

GEN     buchquad(GEN D, double cbach, double cbach2, long RELSUP, long prec);
GEN     form_to_ideal(GEN x);
GEN     getallforms(GEN D, long *pth, GEN *ptz);

/* buch2.c */

GEN     check_and_build_matal(GEN bnf);
GEN     extract_full_lattice(GEN x);
GEN     init_red_mod_units(GEN bnf, long prec);
GEN     isprincipalarch(GEN bnf, GEN col, GEN kNx, GEN e, GEN dx, long *pe);
GEN     red_mod_units(GEN col, GEN z, long prec);
int     trunc_error(GEN x);

/* buch3.c */

GEN     minkowski_bound(GEN D, long N, long r2, long prec);

/* elliptic.c */

GEN     CM_CardEFp(GEN E, GEN p);
GEN     CM_ellap(GEN E, GEN p);
GEN     apell1(GEN e, GEN p);
void    checkpt(GEN z);
GEN     multi_invmod(GEN x, GEN p);
GEN     ratroot(GEN p);
GEN     weipell0(GEN e, long prec, long PREC);

/* es.c */

GEN     Str0(GEN g, long flag);
void    bruti(GEN g, pariout_t *T, int addsign);
const char * eng_ord(long i);
char *  env_ok(char *s);
void    matbruti(GEN g, pariout_t *T);
int     pari_is_dir(char *name);
GEN     readobj(FILE *f, int *ptc);
void    sori(GEN g, pariout_t *T);
void    texi(GEN g, pariout_t *T, int addsign);
void    writeGEN(GEN x, FILE *f);
void    writenamedGEN(GEN x, char *s, FILE *f);

/* galconj.c */

GEN     fixedfieldfactmod(GEN Sp, GEN p, GEN Tmod);
GEN     fixedfieldfactor(GEN L, GEN O, GEN perm, GEN M, GEN den, GEN mod, long x, long y);
GEN     fixedfieldinclusion(GEN O, GEN PL);
GEN     fixedfieldorbits(GEN O, GEN L);
GEN     fixedfieldsympol(GEN O, GEN mod, GEN l, GEN p, long v);
GEN     galoisconj2pol(GEN x, long nbmax, long prec);
GEN     galoiscosets(GEN O, GEN perm);
long    intheadlong(GEN x, GEN mod);
long    isomborne(GEN P, GEN den, GEN p);
GEN     listznstarelts(long m, long p);
GEN     matheadlong(GEN W, GEN mod);
GEN     matrixnorm(GEN M, long prec);
GEN     monomorphismlift(GEN P, GEN S, GEN Q, GEN p, long e);
long    polheadlong(GEN P, long n, GEN mod);
GEN     sympol_aut_evalmod(GEN sym, long g, GEN sigma, GEN Tp, GEN p);
GEN     sympol_eval(GEN v, GEN NS);
GEN     vandermondeinverseprep(GEN L);

/* galois.c */

GEN     partitions(long n);
GEN     polgaloisnamesbig(long n, long k);

/* gen1.c */

int     ff_poltype(GEN *x, GEN *p, GEN *pol);
GEN     gred_frac2(GEN x1, GEN x2);
GEN     gred_rfrac2(GEN x1, GEN x2);
GEN     gred_rfrac_simple(GEN n, GEN d);

/* gen2.c */

void    gopsg2z(GEN (*f) (GEN, GEN), long s, GEN y, GEN z);

/* gen3.c */

GEN     gsubst_expr(GEN pol, GEN from, GEN to);
GEN     inv_ser(GEN b);
GEN     mul_real(GEN x, GEN y);
GEN     poltoser(GEN x, long v, long prec);
GEN     qfbeval(GEN q, GEN x, GEN y);
GEN     rfractoser(GEN x, long v, long prec);

/* groupid.c */

long    group_ident_i(GEN G, GEN S);
long    group_ident_trans(GEN G, GEN S);
long    groupelts_sumorders(GEN S);
GEN     vecgroup_idxlist(GEN L, long order);
long    vecgroup_sumorders(GEN L);

/* ifactor1.c */

GEN     ellfacteur(GEN n, int insist);
long    ifac_decomp(GEN n, long hint);
GEN     ifac_primary_factor(GEN *partial, long *exponent);
void    ifac_realloc(GEN *partial, GEN *where, long new_lg);
GEN     ifac_start(GEN n, long moebius, long hint);
GEN     pollardbrent(GEN n);
ulong   snextpr(ulong p, byteptr *d, long *rcn, long *q, long k);
GEN     squfof(GEN n);

/* init.c */

void    err_recover(long numerr);
GEN     gcopy_av(GEN x, GEN *AVMA);
int     ok_gerepileupto(GEN x);
void    pari_init_defaults(void);

/* nffactor.c */

long    ZM_get_prec(GEN x);
int     nfissplit(GEN nf, GEN x);

/* perm.c */

long    cosets_perm_search(GEN C, GEN p);
GEN     group_export_GAP(GEN G);
GEN     group_export_MAGMA(GEN G);
GEN     perm_conj(GEN s, GEN t);
GEN     perm_generate(GEN S, GEN H, long o);
long    perm_relorder(GEN p, GEN S);
GEN     perm_to_GAP(GEN p);
GEN     quotient_subgroup_lift(GEN C, GEN H, GEN S);

/* polarit1.c */

GEN     Flx_Berlekamp_ker(GEN u, ulong p);
GEN     FpX_Berlekamp_ker(GEN u, GEN p);
GEN     FpX_factcantor(GEN f, GEN pp, long flag);
GEN     FqX_Berlekamp_ker(GEN u, GEN T, GEN q, GEN p);
GEN     FqX_rand(long d1, long v, GEN T, GEN p);
long    FqX_split_Berlekamp(GEN *t, GEN q, GEN T, GEN p);
GEN     Zp_appr(GEN f, GEN a);
int     cmp_padic(GEN x, GEN y);
GEN     factcantor0(GEN f, GEN pp, long flag);
GEN     zrhqr(GEN a, long prec);

/* polarit2.c */

GEN     DDF_roots(GEN pol, GEN polp, GEN p);
GEN     Q_divmuli_to_int(GEN x, GEN d, GEN n);
long    checkdeflate(GEN x);
void    factor_quad(GEN x, GEN res, long *ptcnt);
GEN     factorback_aux(GEN fa, GEN e, GEN (*_mul) (void *, GEN, GEN), GEN (*_pow) (void *, GEN, GEN), void *data);
GEN     matratlift(GEN M, GEN mod, GEN amax, GEN bmax, GEN denom);
GEN     pseudodiv(GEN x, GEN y, GEN *ptr);
long    s_centermod(long x, ulong pp, ulong pps2);
GEN     sort_vecpol_gen(GEN a, int (*cmp) (GEN, GEN));
GEN     sylvestermatrix_i(GEN x, GEN y);

/* polarit3.c */

GEN     Flm_Frobenius_pow(GEN M, long d, GEN T, ulong p);
GEN     Flxq_matrix_pow(GEN y, long n, long m, GEN P, ulong l);
GEN     FpM_Frobenius_pow(GEN M, long d, GEN T, GEN p);
GEN     FpXQ_sqrtl(GEN a, GEN l, GEN T, GEN p, GEN q, long e, GEN r, GEN y, GEN m);
GEN     FpXYQQ_pow(GEN x, GEN n, GEN S, GEN T, GEN p);
GEN     FpX_compositum(GEN A, GEN B, GEN p);
GEN     FpX_direct_compositum(GEN A, GEN B, GEN p);
GEN     FpX_div_by_X_x(GEN a, GEN x, GEN p, GEN *r);
GEN     FpX_gcd_check(GEN x, GEN y, GEN p);
GEN     RgX_to_FpXQX(GEN x, GEN T, GEN p);
GEN     Rg_to_FpXQ(GEN x, GEN T, GEN p);
int     ZX_incremental_CRT(GEN *ptH, GEN Hp, GEN q, GEN qp, ulong p);
GEN     ZX_init_CRT(GEN Hp, ulong p, long v);
ulong   ZY_ZXY_ResBound(GEN A, GEN B, GEN dB);
int     Z_incremental_CRT(GEN *H, ulong Hp, GEN q, GEN qp, ulong p);
GEN     ffinit_Artin_Shreier(GEN ip, long l);
GEN     ffinit_rand(GEN p, long n);
byteptr init_modular(ulong *p);
GEN     polint_triv(GEN xa, GEN ya);

/* rootpol.c */

GEN     FFT(GEN x, GEN Omega);
GEN     FFTinit(long k, long prec);

/* subcyclo.c */

GEN     bnr_to_znstar(GEN bnr, long *complex);
GEN     galoiscyclo(long n, long v);
GEN     subcyclo_complex_bound(pari_sp ltop, GEN V, long prec);
GEN     subcyclo_complex_roots(long n, long real, long prec);
GEN     subcyclo_cyclic(long n, long d, long m, long z, long g, GEN powz, GEN le);
GEN     subcyclo_orbits(long n, GEN H, GEN O, GEN powz, GEN le);
GEN     subcyclo_roots(long n, GEN zl);
GEN     subcyclo_start(long n, long d, long o, GEN borne, long *ptr_val, long *ptr_l);
GEN     znstar_bits(long n, GEN H);
long    znstar_conductor(long n, GEN H);
GEN     znstar_coset_bits(long n, GEN H, long c);
void    znstar_coset_bits_inplace(long n, GEN H, GEN bits, long c);
void    znstar_coset_func(long n, GEN H, void (*func) (void *, long), void *data, long c);
GEN     znstar_cosets(long n, long phi_n, GEN H);
GEN     znstar_elts(long n, GEN H);
GEN     znstar_generate(long n, GEN V);
GEN     znstar_hnf(GEN Z, GEN M);
GEN     znstar_hnf_generators(GEN Z, GEN M);
GEN     znstar_partial_bits(long n, GEN H, long d);
GEN     znstar_partial_coset_bits(long n, GEN H, long d, long c);
void    znstar_partial_coset_bits_inplace(long n, GEN H, GEN bits, long d, long c);
void    znstar_partial_coset_func(long n, GEN H, void (*func) (void *, long), void *data, long d, long c);
GEN     znstar_reduce_modulus(GEN H, long n);

/* thue.c */

GEN     bnfisintnormabs(GEN bnf, GEN a);

/* trans1.c */

GEN     constlog2(long prec);
GEN     padic_sqrt(GEN x);
GEN     padic_sqrtn(GEN x, GEN n, GEN *zetan);
GEN     padic_sqrtn_ram(GEN x, long e);
GEN     padic_sqrtn_unram(GEN x, GEN n, GEN *zetan);
void    pari_init_floats(void);
GEN     rootsof1complex(GEN n, long prec);
GEN     rootsof1padic(GEN n, GEN y);

/* trans2.c */

GEN     cxpsi(GEN s0, long prec);
double  darg(double s, double t);

/* trans3.c */

GEN     bernreal_using_zeta(long n, GEN iz, long prec);
GEN     czeta(GEN s0, long prec);
GEN     inv_szeta_euler(long n, double lba, long prec);
GEN     kbesselnew(GEN n, GEN z, long prec);
GEN     polylogd0(long m, GEN x, long flag, long prec);
GEN     twistpartialzeta(GEN p, GEN q, long f, long c, GEN va, GEN cff);

ENDEXTERN
