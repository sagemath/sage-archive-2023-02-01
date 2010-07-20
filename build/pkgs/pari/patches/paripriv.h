/* $Id: paripriv.h 12463 2010-06-23 21:57:27Z bill $

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
/* hashtables */
hashtable *hashstr_import_static(hashentry *e, ulong size);
void hashstr_dbg(hashtable *h);

/* for qsort */
typedef int (*QSCOMP)(const void *, const void *);

/* for (ulong*) GENs */
typedef ulong *uGEN;
#define ucoeff(a,i,j)  (((uGEN*)(a))[j][i])
#define ugel(a,i) ((uGEN*)(a))[i]

/* swap */
#define lswap(x,y) {long _z=x; x=y; y=_z;}
#define pswap(x,y) {GEN *_z=x; x=y; y=_z;}
#define swap(x,y)  {GEN  _z=x; x=y; y=_z;}
#define dswap(x,y) { double _t=x; x=y; y=_t; }
#define pdswap(x,y) { double* _t=x; x=y; y=_t; }
#define swapspec(x,y, nx,ny) {swap(x,y); lswap(nx,ny);}

/* unused */
GEN ellheightoo(GEN e, GEN z, long prec);
void ellprint(GEN e);

/* generic */
GEN trans_fix_arg(long *prec, GEN *s0, GEN *sig, pari_sp *av, GEN *res);
GEN transc(GEN (*f) (GEN, long), GEN x, long prec);
GEN sort_factor_pol(GEN y, int (*cmp)(GEN,GEN));

/* loops */
GEN incloop(GEN a);
GEN resetloop(GEN a, GEN b);
GEN setloop(GEN a);

/* parser */
void forpari(GEN a, GEN b, GEN node);
void untilpari(GEN a, GEN b);
void whilepari(GEN a, GEN b);
GEN  ifpari(GEN g, GEN a, GEN b);
GEN  andpari(GEN a, GEN b);
GEN  orpari(GEN a, GEN b);
void ifpari_void(GEN g, GEN a, GEN b);
GEN  geval_gp(GEN x, GEN t);

GEN  gadde(GEN *x, GEN y);
GEN  gadd1e(GEN *x);
GEN  gdive(GEN *x, GEN y);
GEN  gdivente(GEN *x, GEN y);
GEN  gdivrounde(GEN *x, GEN y);
GEN  gmode(GEN *x, GEN y);
GEN  gmule(GEN *x, GEN y);
GEN  gshiftle(GEN *x, long n);
GEN  gshiftre(GEN *x, long n);
GEN  gstore(GEN *x, GEN y);
GEN  gsube(GEN *x, GEN y);
GEN  gsub1e(GEN *x);
GEN  gshift_right(GEN x, long n);

GEN  derivnum0(GEN a, GEN code, long prec);
GEN  derivnum1(GEN code, GEN args, long prec);
GEN  direuler0(GEN a, GEN b, GEN code, GEN c);
GEN  divsum(GEN num, GEN code);
void fordiv(GEN a, GEN code);
void forell(long a, long b, GEN code);
void forprime(GEN a, GEN b, GEN code);
void forstep(GEN a, GEN b, GEN s, GEN code);
void forsubgroup(GEN cyc, GEN bound, GEN code);
void forvec(GEN x, GEN code, long flag);
GEN  intcirc0(GEN a, GEN R, GEN code, GEN tab, long prec);
GEN  intfourcos0(GEN a, GEN b, GEN x, GEN code, GEN tab, long prec);
GEN  intfourexp0(GEN a, GEN b, GEN x, GEN code, GEN tab, long prec);
GEN  intfoursin0(GEN a, GEN b, GEN x, GEN code, GEN tab, long prec);
GEN  intfuncinit0(GEN a, GEN b, GEN code, long flag, long m, long prec);
GEN  intlaplaceinv0(GEN sig, GEN x, GEN code, GEN tab, long prec);
GEN  intmellininv0(GEN sig, GEN x, GEN code, GEN tab, long prec);
GEN  intnum0(GEN a, GEN b, GEN code, GEN tab, long prec);
GEN  intnuminit0(GEN a, GEN b, GEN tab, long prec);
GEN  intnuminitgen0(GEN a, GEN b, GEN code, long m, long flag, long prec);
GEN  intnumromb0(GEN a, GEN b, GEN code, long flag, long prec);
GEN  matrice(GEN nlig, GEN ncol, GEN code);
GEN  prodeuler0(GEN a, GEN b, GEN code, long prec);
GEN  prodinf0(GEN a, GEN code, long flag, long prec);
GEN  produit(GEN a, GEN b, GEN code, GEN x);
GEN  somme(GEN a, GEN b, GEN code, GEN x);
GEN  sumalt0(GEN a, GEN code,long flag, long prec);
GEN  suminf0(GEN a, GEN code, long prec);
GEN  sumnum0(GEN a, GEN sig, GEN code, GEN tab, long flag, long prec);
GEN  sumnumalt0(GEN a, GEN sig, GEN code, GEN tab, long flag, long prec);
GEN  sumnuminit0(GEN a, GEN tab, long sgn, long prec);
GEN  sumpos0(GEN a, GEN code, long flag,long prec);
GEN  vecteursmall(GEN nmax, GEN code);
GEN  vecteur(GEN nmax, GEN n);
GEN  vvecteur(GEN nmax, GEN n);
GEN  zbrent0(GEN a, GEN b, GEN code, long prec);

long    loop_break(void);

/* multiprecision */
GEN   addrex01(GEN x);
int   lgcdii(ulong* d, ulong* d1, ulong* u, ulong* u1, ulong* v, ulong* v1, ulong vmax);
ulong rgcduu(ulong d, ulong d1, ulong vmax, ulong* u, ulong* u1, ulong* v, ulong* v1, long *s);
ulong xgcduu(ulong d, ulong d1, int f, ulong* v, ulong* v1, long *s);
ulong xxgcduu(ulong d, ulong d1, int f, ulong* u, ulong* u1, ulong* v, ulong* v1, long *s);
GEN   divgunu(GEN x, ulong i);
GEN   divrunu(GEN x, ulong i);
GEN   init_remiimul(GEN M);
GEN   logagmcx(GEN q, long prec);
GEN   muliispec(GEN x, GEN y, long nx, long ny);
GEN   red_montgomery(GEN T, GEN N, ulong inv);
GEN   remiimul(GEN x, GEN y, GEN invy);
GEN   sqrispec(GEN x, long nx);
GEN   subrex01(GEN x);
GEN   modr_safe(GEN x, GEN y);
ulong *convi(GEN x, long *l);

int approx_0(GEN x, GEN y);
GEN bernfrac_using_zeta(long n);
int OK_bern(long nb, long prec);

/* powers */
GEN    rpowuu(ulong a, ulong n, long prec);
ulong  u_pow10(int n);

/* floats */
double dabs(double s, double t);
void   dcxlog(double s, double t, double *a, double *b);
double dnorm(double s, double t);
double dbllog2(GEN z);
ulong  usqrtsafe(ulong a);

/* "abs" routines for t_REAL ( disregard sign ) */
int absrnz_egal1(GEN x);
int absrnz_egal2n(GEN x);

/* hnf */
GEN hnfadd(GEN m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,GEN extramat,GEN extraC);
GEN hnfadd_i(GEN m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,GEN extramat,GEN extraC);
GEN hnfspec_i(GEN m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,long k0);
GEN hnfspec(GEN m,GEN p,GEN* ptdep,GEN* ptA,GEN* ptC,long k0);
GEN mathnfspec(GEN x, GEN *ptperm, GEN *ptdep, GEN *ptB, GEN *ptC);

GEN LLL_check_progress(GEN Bnorm, long n0, GEN m, int final, long *ti_LLL);
GEN extendedgcd(GEN A);

/* miscellaneous linear algebra */
GEN  imagecomplspec(GEN x, long *nlze);
void rowselect_p(GEN A, GEN B, GEN p, long init);

GEN  dim1proj(GEN prh);
GEN  detcyc(GEN cyc, long *L);

GEN merge_factor_i(GEN f, GEN g);

/* integer factorization / discrete log */
GEN   coprime_part(GEN x, GEN f);
ulong ucoprime_part(ulong x, ulong f);
ulong is_kth_power(GEN x, ulong p, GEN *pt, byteptr d);
long  ifac_decomp_break(GEN n, long (*B)(GEN,GEN,GEN,GEN), GEN s, long hint);
long  ifac_moebius(GEN n, long hint);
long  ifac_issquarefree(GEN n, long hint);
long  ifac_omega(GEN n, long hint);
long  ifac_bigomega(GEN n, long hint);
GEN   ifac_totient(GEN n, long hint);
GEN   ifac_numdiv(GEN n, long hint);
GEN   ifac_sumdivk(GEN n, long k, long hint);
GEN   mpqs(GEN N);
ulong gcduodd(ulong x, ulong y);
long  Z_lvalrem_stop(GEN n, ulong p, int *stop);
long  u_lvalrem_stop(ulong *n, ulong p, int *stop);

/* Polynomials */
/* a) Arithmetic/conversions */
GEN  addmulXn(GEN x, GEN y, long d);
GEN  addshiftpol(GEN x, GEN y, long d);
GEN  lift_if_rational(GEN x);
GEN  monomial(GEN a, long degpol, long v);
GEN  monomialcopy(GEN a, long degpol, long v);
GEN  mulmat_pol(GEN A, GEN x);
GEN  ser2pol_i(GEN x, long lx);
GEN  ser2rfrac_i(GEN x);
GEN  shiftpol_i(GEN x, long v);
GEN  swap_vars(GEN b0, long v);
GEN  RgX_recipspec_shallow(GEN x, long l, long n);

/* b) Modular */
GEN  bezout_lift_fact(GEN T, GEN Tmod, GEN p, long e);
GEN  Kronecker_to_FpXQX(GEN z, GEN pol, GEN p);
GEN  FpX_quad_root(GEN x, GEN p, int unknown);
long FpX_split_Berlekamp(GEN *t, GEN pp);
GEN  FqX_split_all(GEN z, GEN T, GEN p);
long FqX_split_by_degree(GEN *pz, GEN u, GEN q, GEN T, GEN p);
long FqX_split_deg1(GEN *pz, GEN u, GEN q, GEN T, GEN p);
GEN  FqX_split_roots(GEN z, GEN T, GEN p, GEN pol);
GEN  polsym_gen(GEN P, GEN y0, long n, GEN T, GEN N);
GEN  ZXQ_charpoly_sqf(GEN A, GEN B, long *lambda, long v);
GEN  ZX_disc_all(GEN,ulong);
GEN  ZX_resultant_all(GEN A, GEN B, GEN dB, ulong bound);
GEN  ZX_ZXY_resultant_all(GEN A, GEN B, long *lambda, GEN *LPRS);
GEN  RgX_gcd_simple(GEN x, GEN y);
GEN  RgX_extgcd_simple(GEN a, GEN b, GEN *pu, GEN *pv);
GEN  RgXQ_minpoly_naive(GEN y, GEN P);
GEN lift_intern0(GEN x,long v);
#define lift_intern(x) (lift_intern0((x),-1))

/* c) factorization */
double cauchy_bound(GEN p);
GEN chk_factors_get(GEN lt, GEN famod, GEN c, GEN T, GEN N);
long cmbf_maxK(long nb);
GEN ZX_DDF(GEN x);
GEN fact_from_DDF(GEN fa, GEN e, long n);
GEN initgaloisborne(GEN T, GEN dn, long prec, GEN *pL, GEN *pprep, GEN *pdis);
GEN logmax_modulus_bound(GEN p);
GEN polint_i(GEN xa, GEN ya, GEN x, long n, GEN *ptdy);
GEN quicktrace(GEN x, GEN sym);
GEN special_pivot(GEN x);
GEN vandermondeinversemod(GEN L, GEN T, GEN den, GEN mod);
GEN ZX_monic_factorpadic(GEN f, GEN p, long prec);

/* Finite field */

enum { t_FF_FpXQ = 0, t_FF_Flxq = 1, t_FF_F2xq = 2 };

/* for Buchall_param */
enum { fupb_NONE, fupb_RELAT, fupb_LARGE, fupb_PRECI };

/* Allocation / gerepile */
void   debug_stack(void);
void   fill_stack(void);
void   init_dalloc(void);
double *dalloc(size_t n);
void   gerepilecoeffs2(pari_sp av, GEN x, int n, GEN y, int o);
void   minim_alloc(long n, double ***q, GEN *x, double **y,  double **z, double **v);
int    pop_entree_block(entree *ep, long loc);
int    pop_val_if_newer(entree *ep, long loc);
void   gclone_refc(GEN x);

/* Interfaces (GP, etc.) */
pariFILE *pari_last_tmp_file();
void  print_errcontext(const char *msg, const char *s, const char *entry);
void* get_stack(double fraction, long min);
void  init_graph(void);
void  free_graph(void);
void  initout(int initerr);
void  init80col(void);
int   pari_kernel_init(void);
int   pari_last_was_newline(void);
void  pari_set_last_newline(int last);
void  print_functions_hash(const char *s);
void  print_all_user_fun(int member);
GEN   readbin(const char *name, FILE *f, int *vector);
void  term_color(long c);
const char *term_get_color(long c);
int   term_height(void);
int   term_width(void);
void  whatnow_new_syntax(const char *f, long n);
/* gp_colors */
void decode_color(long n, long *c);
enum { c_ERR, c_HIST, c_PROMPT, c_INPUT, c_OUTPUT, c_HELP, c_TIME, c_LAST,
       c_NONE = 0xffffUL };
extern GEN pari_colormap, pari_graphcolors;

/* defaults */
extern ulong precreal;

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
  long sigd;   /* -1 (all) or number of significant digits printed */
  int sp;      /* 0 = suppress whitespace from output */
  int prettyp; /* output style: raw, prettyprint, etc */
  int TeXstyle;
} pariout_t;

void lim_lines_output(char *s, long n, long max);
void gen_output(GEN x, pariout_t *T);

void parsestate_reset(void);
void parsestate_save(struct pari_parsestate *state);
void parsestate_restore(struct pari_parsestate *state);

void compilestate_reset(void);
void compilestate_save(struct pari_compilestate *comp);
void compilestate_restore(struct pari_compilestate *comp);

void evalstate_clone(void);
void evalstate_reset(void);
void evalstate_restore(struct pari_evalstate *state);
void evalstate_save(struct pari_evalstate *state);

/* GP_DATA */
typedef struct {
  gp_hist *hist;
  gp_pp *pp;
  gp_path *path;
  pariout_t *fmt;
  ulong flags, lim_lines;
  char *help, *prompt, *prompt_cont;
  pari_timer *T;
} gp_data;
  /* GP_DATA->flags */

/* The ECHO symbol is already defined in Sage when building on OS X. Since this is paripriv.h,
   it's very unlikely that anything in here is used, so this should be safe. */

enum { QUIET=1, TEST=2, SIMPLIFY=4, CHRONO=8, xxECHOxx=16, STRICTMATCH=32,
       USE_READLINE=64, SECURE=128, EMACS=256, TEXMACS=512, BREAKLOOP=1024,
       RECOVER=2048};

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
  const char *s; /* source */
  char *t, *end; /* target, last char read */
  int in_string, in_comment, more_input, wait_for_brace, downcase;
  Buffer *buf;
} filtre_t;
void init_filtre(filtre_t *F, Buffer *buf);
char *filtre(const char *s, int flag);
void check_filtre(filtre_t *F);

gp_data *default_gp_data(void);
GEN  gp_history(gp_hist *H, long p, char *old, char *entry);

void delete_dirs(gp_path *p);
void gp_expand_path(gp_path *p);
const char *pari_default_path(void);
const char *expand_prompt(const char *prompt, filtre_t *F);

typedef struct input_method {
/* mandatory */
  char * (*fgets)(char *,int,FILE*);
  char * (*getline)(char**, int f, struct input_method*, filtre_t *F);
  int free; /* boolean: must we free the output of getline() ? */
/* for interactive methods */
  const char *prompt, *prompt_cont;
/* for non-interactive methods */
  FILE *file;
} input_method;

int input_loop(filtre_t *F, input_method *IM);
char *file_input(char **s0, int junk, input_method *IM, filtre_t *F);
char *file_getline(Buffer *b, char **s0, input_method *IM);

/* By files */

/* Flx.c */

GEN     FlxX_recipspec(GEN x, long l, long n, long vs);
GEN     FlxX_sub(GEN x, GEN y, ulong p);
GEN     FlxX_subspec(GEN x, GEN y, ulong p, long lx, long ly);
GEN     FlxX_swap(GEN x, long n, long ws);
GEN     zxX_to_Kronecker(GEN P, GEN Q);
GEN     zxX_to_Kronecker_spec(GEN P, long lp, GEN Q);
GEN     Flx_addshift(GEN x, GEN y, ulong p, long d);
GEN     Flx_addspec(GEN x, GEN y, ulong p, long lx, long ly);
GEN     Flx_even_odd_comb(GEN P, ulong u, ulong v, ulong p);
GEN     Flx_mulspec(GEN a, GEN b, ulong p, long na, long nb);
GEN     Flx_negspec(GEN x, ulong p, long l);
GEN     Flx_recipspec(GEN x, long l, long n);
GEN     Flx_sqrspec(GEN a, ulong p, long na);
GEN     Flx_subspec(GEN x, GEN y, ulong p, long lx, long ly);
GEN     Kronecker_to_FlxqX(GEN z, GEN T, ulong p);
GEN     FlxqX_invMontgomery(GEN T, GEN Q, ulong p);
GEN     FlxqX_mulspec(GEN x, GEN y, GEN T, ulong p, long lx, long ly);
GEN     FlxqX_rem_Montgomery(GEN x, GEN mg, GEN T, GEN Q, ulong p);

/* Qfb.c */

GEN     redimagsl2(GEN q, GEN *U);
GEN     redrealsl2(GEN V);
GEN     redrealsl2step(GEN A);

/* alglin1.c */
typedef long (*pivot_fun)(GEN,GEN,long,GEN);
GEN RgM_pivots(GEN x0, GEN data, long *rr, pivot_fun pivot);
void    vecselect_p(GEN A, GEN B, GEN p, long init, long lB);

/* arith1.c */

GEN     bestappr_mod(GEN x, GEN A, GEN B);
int     is_gener_Fp(GEN x, GEN p, GEN p_1, GEN L);
int     is_gener_Fl(ulong x, ulong p, ulong p_1, GEN L);

/* arith2.c */

byteptr initprimes0(ulong maxnum, long *lenp, ulong *lastp);
long    set_optimize(long what, GEN g);

/* base1.c */

void    nfbasic_add_disc(nfbasic_t *T);
void    nfbasic_init(GEN x, long flag, GEN fa, nfbasic_t *T);
GEN     nffromhnfbasis(GEN nf, GEN x);
GEN     nftohnfbasis(GEN nf, GEN x);

/* base2.c */

GEN     gen_if_principal(GEN bnf, GEN x);
int     nfissquarefree(GEN nf, GEN x);
GEN     nfreducemodpr_i(GEN x, GEN prh);
GEN     polsymmodp(GEN g, GEN p);

/* base3.c */

void    check_nfelt(GEN x, GEN *den);
GEN     zk_ei_mul(GEN nf, GEN x, long i);

/* base4.c */

void    check_listpr(GEN x);
GEN     extideal_HNF_mul(GEN nf, GEN x, GEN y);
GEN     factor_norm(GEN x);
GEN     factorbackprime(GEN nf, GEN L, GEN e);
long    val_norm(GEN x, GEN p, long *vz);

/* base5.c */

GEN     check_and_build_nfabs(GEN rnf);
GEN     check_and_build_norms(GEN rnf);
GEN     checkrnfeq(GEN x);

/* bibli1.c */

GEN     plindep(GEN x);
GEN     pslq(GEN x);
GEN     pslqL2(GEN x);

/* buch1.c */

GEN     form_to_ideal(GEN x);
GEN     qfbforms(GEN D);

/* buch2.c */

typedef struct GRHcheck_t { double cD, cN; } GRHcheck_t;
void    init_GRHcheck(GRHcheck_t *S, long N, long R1, double LOGD);
int     GRHok(GRHcheck_t *S, double L, double SA, double SB);
GEN     check_and_build_matal(GEN bnf);
GEN     extract_full_lattice(GEN x);
GEN     init_red_mod_units(GEN bnf, long prec);
GEN     isprincipalarch(GEN bnf, GEN col, GEN kNx, GEN e, GEN dx, long *pe);
GEN     red_mod_units(GEN col, GEN z);

/* buch3.c */

GEN     minkowski_bound(GEN D, long N, long r2, long prec);
int     subgroup_conductor_ok(GEN H, GEN L);
GEN     subgrouplist_cond_sub(GEN bnr, GEN C, GEN bound);

/* elliptic.c */

GEN     CM_CardEFp(GEN E, GEN p);
GEN     weipell0(GEN e, long prec, long PREC);

/* ellsea.c */

void    pari_close_seadata(void);
void    pari_init_seadata(void);

/* es.c */

const char * eng_ord(long i);
char *  env_ok(const char *s);
void    killallfiles(void);
pariFILE* newfile(FILE *f, const char *name, int type);
int     popinfile(void);
GEN     readobj(FILE *f, int *ptc);
pariFILE* try_pipe(const char *cmd, int flag);
void    writeGEN(GEN x, FILE *f);
void    writenamedGEN(GEN x, const char *s, FILE *f);

/* galconj.c */

GEN     galoiscosets(GEN O, GEN perm);
long    intheadlong(GEN x, GEN mod);
long    isomborne(GEN P, GEN den, GEN p);
GEN     listznstarelts(long m, long p);
GEN     matheadlong(GEN W, GEN mod);
GEN     matrixnorm(GEN M, long prec);
GEN     monomorphismlift(GEN P, GEN S, GEN Q, GEN p, long e);
long    polheadlong(GEN P, long n, GEN mod);
GEN     vandermondeinverseprep(GEN L);

/* galois.c */

GEN     polgaloisnamesbig(long n, long k);

/* gen1.c */

int     ff_poltype(GEN *x, GEN *p, GEN *pol);
GEN     gred_frac2(GEN x1, GEN x2);
GEN     gred_rfrac2(GEN x1, GEN x2);
GEN     gred_rfrac_simple(GEN n, GEN d);
GEN     mulcxI(GEN x);
GEN     mulcxmI(GEN x);
GEN     sqr_ser_part(GEN x, long l1, long l2);

/* gen3.c */

GEN     gsubst_expr(GEN pol, GEN from, GEN to);
GEN     poltoser(GEN x, long v, long prec);
GEN     rfractoser(GEN x, long v, long prec);

/* ifactor1.c */

GEN     ellfacteur(GEN n, int insist);
long    ifac_decomp(GEN n, long hint);
GEN     ifac_primary_factor(GEN *partial, long *exponent);
void    ifac_realloc(GEN *partial, GEN *where, long new_lg);
GEN     ifac_start(GEN n, long moebius, long hint);
GEN     pollardbrent(GEN n);
ulong   snextpr(ulong p, byteptr *d, long *rcn, long *q, long k);
GEN     squfof(GEN n);

/* prime.c */

long    BPSW_psp_nosmalldiv(GEN N);
int     Fl_MR_Jaeschke(ulong n, long k);
int     MR_Jaeschke(GEN n, long k);
int     uisprime_nosmalldiv(ulong n);

/* init.c */

void    err_recover(long numerr);
void    pari_init_defaults(void);
void    pari_init_stack(size_t size, size_t old);

/* nffactor.c */

int     nfissplit(GEN nf, GEN x);

/* perm.c */

long    cosets_perm_search(GEN C, GEN p);
GEN     group_export_GAP(GEN G);
GEN     group_export_MAGMA(GEN G);
GEN     perm_generate(GEN S, GEN H, long o);
long    perm_relorder(GEN p, GEN S);
GEN     perm_to_GAP(GEN p);
GEN     quotient_subgroup_lift(GEN C, GEN H, GEN S);

/* polarit1.c */

GEN     F2x_Berlekamp_ker(GEN u);
GEN     Flx_Berlekamp_ker(GEN u, ulong p);
GEN     FpX_Berlekamp_ker(GEN u, GEN p);
GEN     FpX_factcantor(GEN f, GEN pp, long flag);
GEN     FqX_Berlekamp_ker(GEN u, GEN T, GEN q, GEN p);
GEN     FqX_rand(long d1, long v, GEN T, GEN p);
long    FqX_split_Berlekamp(GEN *t, GEN q, GEN T, GEN p);
GEN     Zp_appr(GEN f, GEN a);
int     cmp_padic(GEN x, GEN y);
GEN     factcantor0(GEN f, GEN pp, long flag);
GEN     trivfact(void);

/* polarit2.c */

GEN     sylvestermatrix_i(GEN x, GEN y);

/* QX_factor */

GEN     DDF_roots(GEN pol, GEN polp, GEN p);
void    factor_quad(GEN x, GEN res, long *ptcnt);
long    logint(GEN B, GEN y, GEN *ptq);

/* FpX.c */

GEN     FpX_gcd_check(GEN x, GEN y, GEN p);

/* polarit3.c */

GEN     Flm_Frobenius_pow(GEN M, long d, GEN T, ulong p);
GEN     FpM_Frobenius_pow(GEN M, long d, GEN T, GEN p);
GEN     FpXYQQ_pow(GEN x, GEN n, GEN S, GEN T, GEN p);
GEN     FpX_compositum(GEN A, GEN B, GEN p);
GEN     FpX_direct_compositum(GEN A, GEN B, GEN p);
ulong   ZX_ZXY_ResBound(GEN A, GEN B, GEN dB);
GEN     ffinit_Artin_Shreier(GEN ip, long l);
GEN     ffinit_rand(GEN p, long n);
byteptr init_modular(ulong *p);
GEN     polint_triv(GEN xa, GEN ya);

/* random.c */

void    pari_init_rand(void);

/* rootpol.c */

GEN     FFT(GEN x, GEN Omega);
GEN     FFTinit(long k, long prec);

/* subcyclo.c */

GEN     bnr_to_znstar(GEN bnr, long *complex);
GEN     galoiscyclo(long n, long v);
GEN     polsubcyclo_complex_bound(pari_sp ltop, GEN V, long prec);
GEN     polsubcyclo_complex_roots(long n, long real, long prec);
GEN     polsubcyclo_cyclic(long n, long d, long m, long z, long g, GEN powz, GEN le);
GEN     polsubcyclo_orbits(long n, GEN H, GEN O, GEN powz, GEN le);
GEN     polsubcyclo_roots(long n, GEN zl);
GEN     polsubcyclo_start(long n, long d, long o, GEN borne, long *ptr_val, long *ptr_l);
GEN     znstar_bits(long n, GEN H);
long    znstar_conductor(long n, GEN H);
GEN     znstar_coset_bits(long n, GEN H, long c);
void    znstar_coset_bits_inplace(long n, GEN H, GEN bits, long c);
void    znstar_coset_func(long n, GEN H, void (*func) (void *, long), void *data, long c);
GEN     znstar_cosets(long n, long phi_n, GEN H);
GEN     znstar_elts(long n, GEN H);
GEN     znstar_generate(long n, GEN V);
GEN     znstar_hnf(GEN Z, GEN M);
GEN     znstar_hnf_elts(GEN Z, GEN H);
GEN     znstar_hnf_generators(GEN Z, GEN M);
GEN     znstar_partial_bits(long n, GEN H, long d);
GEN     znstar_partial_coset_bits(long n, GEN H, long d, long c);
void    znstar_partial_coset_bits_inplace(long n, GEN H, GEN bits, long d, long c);
void    znstar_partial_coset_func(long n, GEN H, void (*func) (void *, long), void *data, long d, long c);
GEN     znstar_reduce_modulus(GEN H, long n);
GEN     znstar_small(GEN zn);

/* trans1.c */

void    pari_init_floats(void);
void    pari_close_floats(void);
GEN     rootsof1complex(GEN n, long prec);
GEN     rootsof1padic(GEN n, GEN y);

/* trans2.c */

GEN     cxpsi(GEN s0, long prec);
double  darg(double s, double t);

/* trans3.c */

GEN     bernreal_using_zeta(long n, GEN iz, long prec);
GEN     czeta(GEN s0, long prec);
GEN     inv_szeta_euler(long n, double lba, long prec);
GEN     polylogd0(long m, GEN x, long flag, long prec);
GEN     twistpartialzeta(GEN q, long f, long c, GEN va, GEN cff);

ENDEXTERN
