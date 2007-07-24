/* $Id: paripriv.h,v 1.60 2005/07/21 00:18:31 kb Exp $

Copyright (C) 2004  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

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

/* generic */
GEN arith_proto(long f(GEN), GEN x, int do_error);
GEN arith_proto2(long f(GEN,GEN), GEN x, GEN n);
GEN arith_proto2gs(long f(GEN,long), GEN x, long y);
GEN gassoc_proto(GEN f(GEN,GEN),GEN,GEN);
GEN garith_proto(GEN f(GEN), GEN x, int do_error);
GEN garith_proto2gs(GEN f(GEN,long), GEN x, long y);
GEN trans_fix_arg(long *prec, GEN *s0, GEN *sig, pari_sp *av, GEN *res);
GEN transc(GEN (*f) (GEN, long), GEN x, long prec);

/* loops */
GEN incloop(GEN a);
GEN incpos(GEN a);
GEN resetloop(GEN a, GEN b);
GEN setloop(GEN a);

/* multiprecision */
GEN   addrex01(GEN x);
GEN   addumului(ulong a, ulong b, GEN Y);
void  affr_fixlg(GEN z, GEN y);
GEN   cxnorm(GEN x);
GEN   quadnorm(GEN x);
GEN   divgsns(GEN x, long i);
GEN   divrsns(GEN x, long i);
GEN   init_remiimul(GEN M);
ulong invrev(ulong b);
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
GEN padic_sqrtn(GEN x, GEN n, GEN *zetan);

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
double mylog2(GEN z);
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
GEN hnfall_i(GEN A, GEN *ptB, long remove);
GEN hnf_gauss(GEN A, GEN B);
GEN hnf_invimage(GEN A, GEN b);
GEN hnfmerge_get_1(GEN A, GEN B);
GEN hnfperm_i(GEN A, GEN *ptU, GEN *ptperm);
GEN hnfspec_i(long** m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,long k0);
GEN hnfspec(long** m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,long k0);
GEN idealhermite_aux(GEN nf, GEN x);
GEN mathnfspec(GEN x, GEN *ptperm, GEN *ptdep, GEN *ptB, GEN *ptC);

/* LLL */
GEN lllint_fp_ip(GEN x, long D);
GEN lllfp_marked(long *M, GEN x, long D, long flag, long prec, int gram);
GEN lllint_marked(long *M, GEN x, long D, int g, GEN *h, GEN *f, GEN *B);
GEN LLL_check_progress(GEN Bnorm, long n0, GEN m, int final, long *ti_LLL);
GEN extendedgcd(GEN A);

/* miscellaneous linear algebra */
GEN  diagonal_i(GEN x);
GEN  F2V_red_ip(GEN v);
GEN  gauss_realimag(GEN x, GEN y);
GEN  imagecomplspec(GEN x, long *nlze);
GEN  matqpascal(long n, GEN q);
GEN  R_from_QR(GEN x, long prec);
void rowselect_p(GEN A, GEN B, GEN p, long init);
GEN  split_realimag(GEN x, long r1, long r2);
GEN  sqred1_from_QR(GEN x, long prec);
GEN  supnorm(GEN L, long prec);
GEN  znstar_hnf_elts(GEN Z, GEN H);
GEN  ZV_lincomb(GEN u, GEN v, GEN X, GEN Y);
GEN  vec_setconst(GEN v, GEN x);
GEN  vec_const(long n, GEN x);
GEN  col_const(long n, GEN x);

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
GEN famat_mul(GEN f, GEN g);
GEN famat_pow(GEN f, GEN n);
GEN famat_reduce(GEN fa);
GEN famat_to_arch(GEN nf, GEN fa, long prec);
GEN famat_to_nf_modideal_coprime(GEN nf, GEN g, GEN e, GEN id, GEN EX);
GEN famat_to_nf_modidele(GEN nf, GEN g, GEN e, GEN bid);
GEN famat_to_nf_modidele(GEN nf, GEN g, GEN e, GEN bid);
GEN merge_factor_i(GEN f, GEN g);
GEN to_famat_all(GEN x, GEN y);
GEN to_famat(GEN g, GEN e);
GEN trivfact(void);
GEN vconcat(GEN Q1, GEN Q2);

/* integer factorization / discrete log */
int   BSW_isprime(GEN x);
int   BSW_isprime_small(GEN x);
GEN   coprime_part(GEN x, GEN f);
GEN   decomp_limit(GEN n, GEN limit);
GEN   Fp_PHlog(GEN a, GEN g, GEN p, GEN ord);
GEN   Fp_shanks(GEN x,GEN g0,GEN p, GEN q);
ulong is_kth_power(GEN x, ulong p, GEN *pt, byteptr d);
long  ifac_decomp_break(GEN n, long (*B)(GEN,GEN,GEN,GEN), GEN s, long hint);
int   miller(GEN n, long k);
GEN   mpqs(GEN N);
ulong ucarrecomplet(ulong A);
ulong ugcd(ulong x, ulong y);
long  Z_lvalrem_stop(GEN n, ulong p, int *stop);
long  Z_issquarefree(GEN x);

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
GEN  fix_rfrac_if_pol(GEN x, GEN y);
GEN  gmulXn(GEN x, long d);
GEN  lift_if_rational(GEN x);
GEN  monomial(GEN a, int degpol, int v);
GEN  mulmat_pol(GEN A, GEN x);
long polegal_spec(GEN x, GEN y);
GEN  polrecip_i(GEN x);
GEN  pol_to_monic(GEN pol, GEN *lead);
GEN  revpol(GEN x);
GEN  ser_to_pol_i(GEN x, long lx);
GEN  shiftpol_i(GEN x, long v);
GEN  swap_vars(GEN b0, long v);
GEN  to_polmod(GEN x, GEN mod);
GEN  TR_pol(GEN P, GEN c);

/* b) Modular */
GEN  bezout_lift_fact(GEN T, GEN Tmod, GEN p, long e);
GEN  caractducos(GEN p, GEN x, int v);
GEN  FpXQX_from_Kronecker(GEN z, GEN pol, GEN p);
GEN  FpX_quad_root(GEN x, GEN p, int unknown);
long FpX_split_Berlekamp(GEN *t, GEN pp);
long FqX_is_squarefree(GEN P, GEN T, GEN p);
GEN  FqX_split_all(GEN z, GEN T, GEN p);
long FqX_split_by_degree(GEN *pz, GEN u, GEN q, GEN T, GEN p);
long FqX_split_deg1(GEN *pz, GEN u, GEN q, GEN T, GEN p);
GEN  FqX_split_roots(GEN z, GEN T, GEN p, GEN pol);
GEN  nfgcd(GEN P, GEN Q, GEN nf, GEN den);
GEN  polratlift(GEN P, GEN mod, GEN amax, GEN bmax, GEN denom);
GEN  polsym_gen(GEN P, GEN y0, long n, GEN T, GEN N);
GEN  ZX_caract_sqf(GEN A, GEN B, long *lambda, long v);
GEN  ZX_disc_all(GEN,ulong);
long ZX_get_prec(GEN x);
GEN  ZX_resultant_all(GEN A, GEN B, GEN dB, ulong bound);
/* GEN  ZY_ZXY_resultant_all(GEN A, GEN B0, long *lambda, GEN *LPRS);
GEN  ZY_ZXY_resultant(GEN A, GEN B0, long *lambda);
*/
GEN  RgXQ_u_pow(GEN x, ulong n, GEN T);
GEN  RgX_gcd_simple(GEN x, GEN y);
GEN  RgX_extgcd_simple(GEN a, GEN b, GEN *pu, GEN *pv);

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
GEN vandermondeinverse(GEN L, GEN T, GEN den, GEN prep);
GEN vandermondeinversemod(GEN L, GEN T, GEN den, GEN mod);
GEN ZX_squff(GEN f, GEN *ex);
GEN ZX_monic_factorpadic(GEN f, GEN p, long prec);

/* Number fields */
GEN    arch_mul(GEN x, GEN y);
GEN    archstar_full_rk(GEN x, GEN bas, GEN v, GEN gen);
GEN    bnrGetSurj(GEN bnr1, GEN bnr2);
GEN    check_and_build_cycgen(GEN bnf);
double check_bach(double cbach, double B);
GEN    _checkbnf(GEN bnf);
GEN    _checknf(GEN nf);
void   check_ZKmodule(GEN x, char *s);
void   dbg_rel(long s, GEN col);
GEN    element_mulidid(GEN nf, long i, long j);
GEN    element_powid_mod_p(GEN nf, long I, GEN n, GEN p);
GEN    eltabstorel(GEN x, GEN T, GEN pol, GEN k);
GEN    eltmulid_get_table(GEN nf, long i);
GEN    eltreltoabs(GEN rnfeq, GEN x);
GEN    galoisbig(GEN x, long prec);
GEN    get_arch(GEN nf,GEN x,long prec);
GEN    get_arch_real(GEN nf,GEN x,GEN *emb,long prec);
GEN    get_bas_den(GEN bas);
GEN    get_hnfid(GEN nf, GEN x);
GEN    get_mul_table(GEN x,GEN bas,GEN invbas);
GEN    get_nfindex(GEN bas);
GEN    get_proj_modT(GEN basis, GEN T, GEN p);
GEN    get_roots(GEN x,long r1,long prec);
GEN    get_theta_abstorel(GEN T, GEN pol, GEN k);
GEN    idealaddtoone_i(GEN nf, GEN x, GEN y);
GEN    idealcoprime_fact(GEN nf, GEN x, GEN fy);
GEN    idealhermite_aux(GEN nf, GEN x);
GEN    idealsqrtn(GEN nf, GEN x, GEN gn, int strict);
GEN    init_unif_mod_fZ(GEN L);
GEN    init_units(GEN BNF);
long   int_elt_val(GEN nf, GEN x, GEN p, GEN bp, GEN *t);
GEN    isprincipalfact(GEN bnf,GEN P, GEN e, GEN C, long flag);
GEN    make_integral(GEN nf, GEN L0, GEN f, GEN *listpr);
GEN    maxord_i(GEN p, GEN f, long mf, GEN w, long flag);
GEN    modprM(GEN z, GEN nf,GEN modpr);
GEN    modprV(GEN z, GEN nf,GEN modpr);
GEN    modprX(GEN x, GEN nf,GEN modpr);
GEN    nfpol_to_Flx(GEN nf, GEN pol, ulong *ptp);
GEN    nfreducemodideal_i(GEN x0,GEN ideal);
GEN    nfrootsall_and_pr(GEN nf, GEN pol);
GEN    norm_by_embed(long r1, GEN x);
GEN    perm_to_arch(GEN nf, GEN archp);
GEN    pidealprimeinv(GEN nf, GEN x);
GEN    primedec_apply_kummer(GEN nf,GEN pol,long e,GEN p);
GEN    prodid(GEN nf, GEN I);
GEN    pr_norm(GEN pr);
GEN    quadhilbertreal(GEN D, long prec);
GEN    rnfallbase(GEN nf, GEN pol, GEN *pD, GEN *pd, GEN *pfi);
GEN    _rnfequation(GEN A, GEN B, long *pk, GEN *pLPRS);
GEN    special_anti_uniformizer(GEN nf, GEN pr);
GEN    sqr_by_tab(GEN tab, GEN x);
GEN    subgroupcondlist(GEN cyc, GEN bound, GEN listKer);
GEN    T2_from_embed_norm(GEN x, long r1);
void   testprimes(GEN bnf, ulong bound);
GEN    to_Fp_simple(GEN nf, GEN x, GEN ffproj);
GEN    unif_mod_fZ(GEN pr, GEN F);
GEN    unnf_minus_x(GEN x);
void   wr_rel(GEN col);
GEN    zideallog_sgn(GEN nf, GEN x, GEN sgn, GEN bid);
GEN    zlog_units(GEN nf, GEN U, GEN sgnU, GEN bid);
GEN    zlog_units_noarch(GEN nf, GEN U, GEN bid);
GEN    zsign_from_logarch(GEN Larch, GEN invpi, GEN archp);

/* Dedekind zeta */
GEN  zeta_get_limx(long r1, long r2, long bit);
long zeta_get_i0(long r1, long r2, long bit, GEN limx);
long zeta_get_N0(GEN C,  GEN limx);


/* Allocation / gerepile */
void   init_dalloc();
double *dalloc(size_t n);
void   gerepilecoeffs2(pari_sp av, GEN x, int n, GEN y, int o);
void   minim_alloc(long n, double ***q, GEN *x, double **y,  double **z, double **v);
int    pop_entree_bloc(entree *ep, long loc);
int    pop_val_if_newer(entree *ep, long loc);

/* Interfaces (GP, etc.) */
void  errcontext(char *msg, char *s, char *entry);
void  free_graph(void);
GEN   geni(void);
void* get_stack(double,int);
GEN   gpreadseq(char *c, int strict);
void  init_defaults(int force);
void  init_graph(void);
void  initout(int initerr);
char* itostr(GEN x, int minus);
void  kill_from_hashlist(entree *ep, long n);
void  member_err(char *s);
int   pari_kernel_init(void);
void  pari_sig_init(void (*f)(int));
void  print_functions_hash(const char *s);
void  print_fun_list(char **matches, int nbli);
void  print_user_fun(entree *ep);
void  print_user_member(entree *ep);
void  print_all_user_fun(void);
void  print_all_user_member(void);
int   term_width(void);
void  texmacs_completion(char *s, long pos);
void  var_make_safe();
int   whatnow(char *s, int flag);
void  whatnow_new_syntax(char *f, long n);

/* defaults */
#define is_default(s) setdefault((s),"",d_EXISTS) == gen_1
enum { d_SILENT, d_ACKNOWLEDGE, d_INITRC, d_RETURN, d_EXISTS };
extern ulong prec;

GEN sd_TeXstyle(const char *v, int flag);
GEN sd_colors(char *v, int flag);
GEN sd_compatible(const char *v, int flag);
GEN sd_datadir(char *v, int flag);
GEN sd_debug(const char *v, int flag);
GEN sd_debugfiles(const char *v, int flag);
GEN sd_debugmem(const char *v, int flag);
GEN sd_echo(const char *v, int flag);
GEN sd_factor_add_primes(char *v, int flag);
GEN sd_filename(const char *v, int flag, char *s, char **f);
GEN sd_format(const char *v, int flag);
GEN sd_help(char *v, int flag);
GEN sd_histsize(const char *v, int flag);
GEN sd_lines(const char *v, int flag);
GEN sd_log(const char *v, int flag);
GEN sd_logfile(const char *v, int flag);
GEN sd_new_galois_format(char *v, int flag);
GEN sd_output(const char *v, int flag);
GEN sd_parisize(const char *v, int flag);
GEN sd_path(char *v, int flag);
GEN sd_prettyprinter(char *v, int flag);
GEN sd_primelimit(const char *v, int flag);
GEN sd_prompt(const char *v, int flag);
GEN sd_prompt_cont(const char *v, int flag);
GEN sd_psfile(const char *v, int flag);
GEN sd_realprecision(const char *v, int flag);
GEN sd_rl(const char *v, int flag);
GEN sd_secure(const char *v, int flag);
GEN sd_seriesprecision(const char *v, int flag);
GEN sd_simplify(const char *v, int flag);
GEN sd_strictmatch(const char *v, int flag);
GEN sd_timer(const char *v, int flag);
GEN setdefault(const char *s, const char *v, int flag);

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
/* enum { QUIET=1, TEST=2, SIMPLIFY=4, CHRONO=8, ECHO=16, STRICTMATCH=32,
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
char *color_prompt(char *prompt);
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
