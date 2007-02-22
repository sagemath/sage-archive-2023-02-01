"""nodoctest
Interface to PARI
"""

doc=r"""
Interface to PARI

This module provides an interface to the PARI C-library.
This interface is modeled on the GP interpreter, and compared
to GP it has advantages and disadvantages.  First
we list some advantages:
\begin{enumerate}
\item  Improved memory management and garbage collection:
\begin{enumerate}
           \item Automatic garbage collection:  Each SAGE/PARI object
              has its own piece of memory allocated from the Python
              heap.  When objects go out of scope the Python memory
              manager deletes them.  (Implementation: I grafted
              genuine memory management onto PARI by creatinging a
              mini ``PARI stack'' in the Python heap for each GEN that
              is wrapped by a Python object.  In some cases it is
              necessary to copy this GEN back to the main PARI stack
              for computations, which entails a performance penalty.)


            \item It is less likely that the ``PARI stack overflows''.
            This is because the main stack is {\em only} used for
            individual computations; as soon as an object is computed
            it is moved of the stack into its own chunk of memory in
            the Python heap.
            \end{enumerate}

      \item Exception handling: try, except, raise.  GP doesn't have
      any exception handling.

      \item Python is a genuine object-oriented language with user
        defined classes, multiple inheritence, polymorphism, etc.  The
        GP language was written initially for testing the C library,
        and is not as sophisticated a language as Python.

      \item Object serialization: This can make loading and saving
        objects to disk easier.

            sage: a = pari(5)
            sage: import pickle
            sage: s = pickle.dumps(a);
            sage: z = pickle.loads(s);
            sage: print z
            5
            sage: type(z)
            <type '_py_pari.gen'>
            sage: a == z
            True
            sage: a is z
            False

      \item Documentation: Each method has (or will have) example
      usage, so there will be more examples of usage than in the PARI
      documentation.

\end{enumerate}

There are also several disadvantages:

\begin{enumerate}
 \item \module{py_pari} currently does not export all the functionality
 of PARI.  For example, there is currently no support for numerical integration.

 \item There is overhead associated with doing more complicated memory
 management than GP, so SAGE/PARI can be slower than GP.

 \item PARI functions run at the C-level and cannot be interrupted from
        the Python interpreter with ctrl-C until they return from C level.
        This is in accord with how the Python interpreter works.

 \item Some of the PARI error and log messages print to stdout.  This
        can be annoying, since it might interfere with output the
        user is trying to send to a file.  This is fortunately rare.
\end{enumerate}



"""


#
"""
TODO
* Make it so when asking for, e.g., the n-th prime, if not enough primes have
  been precomputed, py_pari automatically precomputes enough primes and answers
  the question.  Make py_pari robust like MAGMA, and not lame like gp.

"""

#########################################################
## NOTES about memory management.

## The idea for the trick I arrived at for adding memory management to
## PARI is the following.  In some cases it has a performance penalty,
## but in practice (and benchmarks) so far it is barely noticeable and
## worth the tradeoff in my opinion.  The memory model in PARI is one
## big stack.  When this stack runs out of space, game over.  Or at
## least that's what the documentation suggest.  After much reading of
## the source code I figured out that when you type "allocatemem()" in
## the GP interpreter, it creates yet another stack -- so now there
## are TWO stacks.  And if you do it again, it creates a third, etc.
## There are some global variables (top, bot, and avma) that define
## what is the current stack for PARI library functions.  So you can
## take any chunk of memory at all, point top, bot and avma at it, and
## all the PARI library functions will view that chunck of memory as
## "the" PARI stack.

## So here's what I do.  First, when you create a PARI object (a GEN),
## you do this by using some PARI library function, which puts the
## result and all objects created as part of your object on the stack.
## Also the avma pointer is moved down (since the PARI stack is upside
## down).  So what I do is just forcecopy into the Python heap all the
## memory between where avma is now and where it was before making
## your object.  I then set avma back to where it was before all this
## happened, thus freeing that memory.  This copy, which is a little
## piece of memory, is managed by the Python memory manager, so when
## the associated Python object goes out of scope (no more references
## to it), then it is garbage collected.  This garbage collection
## doesn't affect the PARI stack at all, and doesn't require anything
## like calling grepilemany.

## And that's basically it.  Also, when a PARI computation fails for
## lack of memory, I just toss the main stack and make a new one.
## Since none of my objects are stored on the stack, this doesn't mess
## anything up.

## Also, for some reason, some PARI operations fail if the GEN they
## take as input is not on the main stack.  Thus in many cases I use
## the forcecopy() PARI library function to copy the object from the
## Python heap back to the PARI stack, do the operation, then copy
## back the result.  There is overhead with these object copies, but
## in practice it isn't too much, and it's worth it for such an easy
## way of adding better memory management to PARI.

## In summary: the PARI library does it's game on the stack, but all
## the user variables are stored elsewhere, and only put on the stack
## for a moment when taking part in a computation.  It's very simple,
## gives proper memory management, and doesn't require rewriting even
## a single line of PARI, with the only tradeoffs lots of extra
## memcpy's, which are extremely fast anyways.

#########################################

cdef int MEMERR
MEMERR = 25

import gc

include "interrupt.pxi"

cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void free(void *ptr)
    void *malloc(size_t size)
    void *realloc(void *ptr, size_t size)
    void exit (int __status)


cdef extern from "string.h":
    void *memmove(void *dest, void *src, size_t n)
    void *memcpy(void *dest, void *src, size_t n)

cdef extern from "Python.h":
    void PyMem_Free(void *p)
    void* PyMem_Malloc(int)
    void* PyMem_Realloc(void *p, size_t n)

cdef extern from 'setjmp.h':
    struct __jmp_buf_tag:
        pass
    ctypedef __jmp_buf_tag jmp_buf
    int setjmp (jmp_buf __env)
    int longjmp (jmp_buf __env, int val)

ctypedef unsigned long ulong

cdef extern from 'pari/paricfg.h':
    char* PARIVERSION


cdef extern from 'pari/pari.h':
    ctypedef long* GEN
    ctypedef void entree   # fake -- not used
    ctypedef void FILE   # fake -- not used
    ctypedef void pariFILE   # fake -- not used
    ctypedef void va_list   # fake -- not used
    ctypedef void pari_timer   # fake -- not used
    ctypedef void stackzone   # fake -- not used


    # misc ...
    GEN     gp_read_str(char *t)
    long    typ(GEN x)
    long    itos(GEN x)
    double  gtodouble(GEN x)
    GEN     stoi(long s)
    #GEN     dbltor(double s)
    # These types are actually an enum type, but I can't get Pyrex to "properly"
    # wrap enums.  It doesn't matter as long as they are treated as ints by pyrexc.
    extern int t_INT, t_REAL, t_INTMOD, t_FRAC, t_COMPLEX, t_PADIC, t_QUAD, \
           t_POLMOD, t_POL, t_SER, t_RFRAC, t_QFR, t_QFI, t_VEC, t_COL,  \
           t_MAT, t_LIST, t_STR, t_VECSMALL

    extern unsigned long precdl, prec
    extern char * diffptr

    extern int BYTES_IN_LONG

    ctypedef unsigned long pari_sp
    pari_sp avma, bot, top, zonetop

    GEN     sd_realprecision(char* n, int flag)


    # Flx.c
    GEN     Fl_to_Flx(ulong x, long sv)
    GEN     Flm_id(long n)
    GEN     Flm_to_FlxV(GEN x, long sv)
    GEN     Flm_to_FlxX(GEN x, long v,long w)
    GEN     Flm_to_ZM(GEN z)
    GEN     Flv_to_Flx(GEN x, long vs)
    GEN     Flv_to_ZC(GEN z)
    GEN     Flv_to_ZV(GEN z)
    GEN     Flv_polint(GEN xa, GEN ya, ulong p, long vs)
    GEN     Flv_roots_to_pol(GEN a, ulong p, long vs)
    GEN     Flx_Fl_mul(GEN y, ulong x, ulong p)
    GEN     Flx_to_Flv(GEN x, long N)
    GEN     Flx_to_ZX(GEN z)
    GEN     Flx_to_ZX_inplace(GEN z)
    GEN     Flx_add(GEN x, GEN y, ulong p)
    GEN     Flx_deriv(GEN z, ulong p)
    GEN     Flx_div_by_X_x(GEN a, ulong x, ulong p, ulong *rem)
    GEN     Flx_divrem(GEN x, GEN y, ulong p, GEN *pr)
    ulong   Flx_eval(GEN x, ulong y, ulong p)
    GEN     Flx_extgcd(GEN a, GEN b, ulong p, GEN *ptu, GEN *ptv)
    ulong   Flx_extresultant(GEN a, GEN b, ulong p, GEN *ptU, GEN *ptV)
    GEN     Flx_gcd(GEN a, GEN b, ulong p)
    GEN     Flx_gcd_i(GEN a, GEN b, ulong p)
    GEN     Flx_invmontgomery(GEN T, ulong p)
    int     Flx_is_squarefree(GEN z, ulong p)
    GEN     Flx_mul(GEN x, GEN y, ulong p)
    GEN     Flx_neg(GEN x, ulong p)
    GEN     Flx_neg_inplace(GEN x, ulong p)
    GEN     Flx_normalize(GEN z, ulong p)
    GEN     Flx_pow(GEN x, long n, ulong p)
    GEN     Flx_recip(GEN x)
    GEN     Flx_red(GEN z, ulong p)
    GEN     Flx_rem_montgomery(GEN x, GEN mg, GEN T, ulong p)
    GEN     Flx_rem(GEN x, GEN y, ulong p)
    GEN     Flx_renormalize(GEN x, long l)
    ulong   Flx_resultant(GEN a, GEN b, ulong p)
    GEN     Flx_shift(GEN a, long n)
    GEN     Flx_sqr(GEN x, ulong p)
    GEN     Flx_sub(GEN x, GEN y, ulong p)
    long    Flx_valuation(GEN x)
    GEN     FlxM_to_ZXM(GEN z)
    GEN     FlxV_Flv_innerprod(GEN V, GEN W, ulong p)
    GEN     FlxV_to_Flm(GEN v, long n)
    GEN     FlxV_to_ZXC(GEN x)
    GEN     FlxX_add(GEN P, GEN Q, ulong p)
    GEN     FlxX_shift(GEN a, long n)
    GEN     FlxX_to_Flm(GEN v, long n)
    GEN     FlxX_to_ZXX(GEN B)
    GEN     FlxYqQ_pow(GEN x, GEN n, GEN S, GEN T, ulong p)
    GEN     Flxq_inv(GEN x,GEN T,ulong p)
    GEN     Flxq_invsafe(GEN x, GEN T, ulong p)
    GEN     Flxq_mul(GEN y,GEN x,GEN T,ulong p)
    GEN     Flxq_pow(GEN x, GEN n, GEN T, ulong p)
    GEN     Flxq_powers(GEN x, long l, GEN T, ulong p)
    GEN     Flxq_sqr(GEN y,GEN T,ulong p)
    GEN     FlxqX_normalize(GEN z, GEN T, ulong p)
    GEN     FlxqX_Flxq_mul(GEN P, GEN U, GEN T, ulong p)
    GEN     FlxqX_red(GEN z, GEN T, ulong p)
    GEN     FlxqX_mul(GEN x, GEN y, GEN T, ulong p)
    GEN     FlxqX_safegcd(GEN P, GEN Q, GEN T, ulong p)
    GEN     FlxqX_sqr(GEN x, GEN T, ulong p)
    GEN     FlxqX_divrem(GEN x, GEN y, GEN T, ulong p, GEN *pr)
    GEN     FlxqXQ_pow(GEN x, GEN n, GEN S, GEN T, ulong p)
    GEN     Z_to_Flx(GEN x, ulong p, long v)
    GEN     ZM_to_Flm(GEN x, ulong p)
    GEN     ZV_to_Flv(GEN x, ulong p)
    GEN     ZX_to_Flx(GEN x, ulong p)
    GEN     ZXV_to_FlxV(GEN v, ulong p)
    GEN     ZXX_to_FlxX(GEN B, ulong p, long v)
    GEN     polx_Flx(long sv)
    GEN     zero_Flx(long sv)

     # alglin1.c

    GEN     Flm_Flv_mul(GEN x, GEN y, ulong p)
    GEN     Flm_deplin(GEN x, ulong p)
    GEN     Flm_indexrank(GEN x, ulong p)
    GEN     Flm_inv(GEN x, ulong p)
    GEN     Flm_ker(GEN x, ulong p)
    GEN     Flm_ker_sp(GEN x, ulong p, long deplin)
    GEN     Flm_mul(GEN x, GEN y, ulong p)
    GEN     FlxqM_ker(GEN x, GEN T, ulong p)
    GEN     FpC_FpV_mul(GEN x, GEN y, GEN p)
    GEN     FpM_FpV_mul(GEN x, GEN y, GEN p)
    GEN     FpM_deplin(GEN x, GEN p)
    GEN     FpM_image(GEN x, GEN p)
    GEN     FpM_intersect(GEN x, GEN y, GEN p)
    GEN     FpM_inv(GEN x, GEN p)
    GEN     FpM_invimage(GEN m, GEN v, GEN p)
    GEN     FpM_ker(GEN x, GEN p)
    GEN     FpM_mul(GEN x, GEN y, GEN p)
    long    FpM_rank(GEN x, GEN p)
    GEN     FpM_indexrank(GEN x, GEN p)
    GEN     FpM_suppl(GEN x, GEN p)
    GEN     FqM_gauss(GEN a, GEN b, GEN T, GEN p)
    GEN     FqM_ker(GEN x, GEN T, GEN p)
    GEN     FqM_suppl(GEN x, GEN T, GEN p)
    GEN     QM_inv(GEN M, GEN dM)
    GEN     ZM_inv(GEN M, GEN dM)
    void    appendL(GEN x, GEN t)
    GEN     cget1(long l, long t)
    GEN     concat(GEN x, GEN y)
    GEN     concatsp(GEN x, GEN y)
    GEN     concatsp3(GEN x, GEN y, GEN z)
    GEN     deplin(GEN x)
    GEN     det(GEN a)
    GEN     det0(GEN a,long flag)
    GEN     det2(GEN a)
    GEN     dethnf(GEN x)
    GEN     dethnf_i(GEN mat)
    GEN     detint(GEN x)
    GEN     diagonal(GEN x)
    GEN     eigen(GEN x, long prec)
    GEN     extract(GEN x, GEN l)
    GEN     extract0(GEN x, GEN l1, GEN l2)
    GEN     gaddmat(GEN x, GEN y)
    GEN     gaddmat_i(GEN x, GEN y)
    GEN     gauss(GEN a, GEN b)
    GEN     gaussmodulo(GEN M, GEN D, GEN Y)
    GEN     gaussmodulo2(GEN M, GEN D, GEN Y)
    GEN     gscalcol(GEN x, long n)
    GEN     gscalcol_i(GEN x, long n)
    GEN     gscalcol_proto(GEN z, GEN myzero, long n)
    GEN     gscalmat(GEN x, long n)
    GEN     gscalsmat(long x, long n)
    GEN     gtomat(GEN x)
    GEN     gtrans(GEN x)
    GEN     gtrans_i(GEN x)
    int     hnfdivide(GEN A, GEN B)
    GEN     idmat(long n)
    GEN     idmat_intern(long n,GEN myun,GEN myzero)
    GEN     image(GEN x)
    GEN     image2(GEN x)
    GEN     imagecompl(GEN x)
    GEN     indexrank(GEN x)
    GEN     inverseimage(GEN mat, GEN y)
    long    isdiagonal(GEN x)
    long    isscalarmat(GEN x, GEN s)
    GEN     ker(GEN x)
    GEN     keri(GEN x)
    GEN     matextract(GEN x, GEN l1, GEN l2)
    GEN     matimage0(GEN x,long flag)
    GEN     matker0(GEN x, long flag)
    GEN     matmuldiagonal(GEN x, GEN d)
    GEN     matmultodiagonal(GEN x, GEN y)
    GEN     matsolvemod0(GEN M, GEN D, GEN Y,long flag)
    GEN     mattodiagonal(GEN m)
    GEN     mattodiagonal_i(GEN m)
    long    rank(GEN x)
    GEN     row(GEN A, long x1)
    GEN     row_i(GEN A, long x0, long x1, long x2)
    GEN     rowextract_i(GEN A, long x1, long x2)
    GEN     rowextract_ip(GEN A, GEN p, long x1, long x2)
    GEN     rowextract_p(GEN A, GEN p)
    GEN     sindexrank(GEN x)
    GEN     sum(GEN v, long a, long b)
    GEN     suppl(GEN x)
    GEN     vconcat(GEN A, GEN B)
    GEN     vec_ei(long n, long i)
    GEN     vec_Cei(long n, long i, GEN c)
    GEN     vecextract_i(GEN A, long y1, long y2)
    GEN     vecextract_ip(GEN A, GEN p, long y1, long y2)
    GEN     vecextract_p(GEN A, GEN p)

    # alglin2.c

    GEN     QuickNormL1(GEN x,long prec)
    GEN     QuickNormL2(GEN x,long prec)
    GEN     ZM_to_zm(GEN z)
    GEN     adj(GEN x)
    GEN     assmat(GEN x)
    GEN     caract(GEN x, int v)
    GEN     caract2(GEN p, GEN x, int v)
    GEN     caradj(GEN x, long v, GEN *py)
    GEN     caradj0(GEN x, long v)
    GEN     carhess(GEN x, long v)
    GEN     charpoly0(GEN x, int v,long flag)
    GEN     conjvec(GEN x,long prec)
    GEN     gconj(GEN x)
    GEN     gnorm(GEN x)
    GEN     gnorml1(GEN x,long prec)
    GEN     gnorml2(GEN x)
    GEN     gtrace(GEN x)
    GEN     hess(GEN x)
    GEN     hnf(GEN x)
    GEN     hnfall(GEN x)
    GEN     hnflll(GEN x)
    GEN     hnflll_i(GEN A, GEN *ptB, int remove)
    GEN     hnfmod(GEN x, GEN detmat)
    GEN     hnfmodid(GEN x,GEN p)
    GEN     hnfmodidpart(GEN x, GEN p)
    GEN     hnfperm(GEN x)
    GEN     intersect(GEN x, GEN y)
    GEN     jacobi(GEN a, long prec)
    GEN     matfrobenius(GEN M, long flag)
    GEN     matrixqz(GEN x, GEN pp)
    GEN     matrixqz0(GEN x, GEN pp)
    GEN     matrixqz2(GEN x)
    GEN     matrixqz3(GEN x)
    GEN     signat(GEN a)
    GEN     smith(GEN x)
    GEN     smith2(GEN x)
    GEN     smithall(GEN x, GEN *ptU, GEN *ptV)
    GEN     smithclean(GEN z)
    GEN     sqred(GEN a)
    GEN     sqred1(GEN a)
    GEN     sqred1intern(GEN a)
    GEN     sqred3(GEN a)
    GEN     zm_to_ZM(GEN z)
    GEN     zx_to_ZX(GEN z)

    # anal.c

    void    addhelp(entree *ep, char *s)
    void    delete_named_var(entree *ep)
    long    delete_var()
    entree  *fetch_named_var(char *s, int doerr)
    long    fetch_user_var(char *s)
    long    fetch_var()
    GEN     flisexpr(char *t)
    GEN     flisseq(char *t)
    void    freeep(entree *ep)
    entree  *gp_variable(char *s)
    long    hashvalue(char *s)
    entree* install(void *f, char *name, char *code)
    entree  *is_entry(char *s)
    void    kill0(entree *ep)
    GEN     lisexpr(char *t)
    GEN     lisseq(char *t)
    long    manage_var(long n, entree *ep)
    void    name_var(long n, char *s)
    GEN     readseq(char *c, int strict)
    GEN     strtoGENstr(char *s)
    GEN     type0(GEN x)

    # arith1.c

    GEN     bestappr0(GEN x, GEN a, GEN b)
    GEN     bestappr(GEN x, GEN k)
    long    carrecomplet(GEN x, GEN *pt)
    long    cgcd(long a,long b)
    void    check_quaddisc(GEN x, long *s, long *r, char *f)
    GEN     chinese(GEN x, GEN y)
    GEN     chinois(GEN x, GEN y)
    GEN     classno2(GEN x)
    GEN     classno(GEN x)
    long    clcm(long a,long b)
    GEN     compimag(GEN x, GEN y)
    GEN     compimagraw(GEN x, GEN y)
    GEN     compraw(GEN x, GEN y)
    GEN     compreal(GEN x, GEN y)
    GEN     comprealraw(GEN x, GEN y)
    GEN     contfrac0(GEN x, GEN b, long flag)
    GEN     fibo(long n)
    GEN     Fp_gener_fact(GEN p, GEN fa)
    GEN     Fp_gener(GEN p)
    GEN     FpXQ_gener(GEN T, GEN p)
    GEN     fundunit(GEN x)
    GEN     gboundcf(GEN x, long k)
    GEN     gcarrecomplet(GEN x, GEN *pt)
    GEN     gcarreparfait(GEN x)
    GEN     gcf2(GEN b, GEN x)
    GEN     gcf(GEN x)
    GEN     gener(GEN m)
    GEN     gfundunit(GEN x)
    GEN     ggener(GEN m)
    long    gisanypower(GEN x, GEN *pty)
    GEN     gisfundamental(GEN x)
    GEN     gisprime(GEN x, long flag)
    GEN     gispseudoprime(GEN x, long flag)
    GEN     gispsp(GEN x)
    GEN     gkrogs(GEN x, long y)
    GEN     gkronecker(GEN x, GEN y)
    GEN     gmillerrabin(GEN n, long k)
    GEN     gnextprime(GEN n)
    GEN     gprecprime(GEN n)
    GEN     gracine(GEN a)
    GEN     gregula(GEN x, long prec)
    GEN     hclassno(GEN x)
    long    hil0(GEN x, GEN y, GEN p)
    long    hil(GEN x, GEN y, GEN p)
    long    isanypower(GEN x, GEN *y)
    long    isfundamental(GEN x)
    long    ispower(GEN x, GEN k, GEN *pty)
    long    isprimeAPRCL(GEN N)
    long    isprime(GEN x)
    long    isprimeSelfridge(GEN x)
    long    ispseudoprime(GEN x, long flag)
    long    ispsp(GEN x)
    long    krois(GEN x, long y)
    long    kronecker(GEN x, GEN y)
    long    krosi(long s, GEN x)
    long    kross(long x, long y)
    long    krouu(ulong x, ulong y)
    GEN     lcmii(GEN a, GEN b)
    GEN     mpfact(long n)
    GEN     mpfactr(long n, long prec)
    GEN     Fp_inv(GEN a, GEN m)
    GEN     Fp_invsafe(GEN a, GEN m)
    GEN     Fp_pow(GEN a, GEN n, GEN m)
    GEN     Fp_sqrt(GEN a, GEN p)
    GEN     Fp_sqrtn(GEN a, GEN n, GEN p, GEN *zetan)
    GEN     nucomp(GEN x, GEN y, GEN l)
    GEN     nudupl(GEN x, GEN l)
    GEN     nupow(GEN x, GEN n)
    GEN     order(GEN x)
    GEN     pnqn(GEN x)
    GEN     powraw(GEN x, long n)
    GEN     powrealraw(GEN x, long n)
    GEN     primeform(GEN x, GEN p, long prec)
    GEN     Qfb0(GEN x, GEN y, GEN z, GEN d, long prec)
    GEN     qfbclassno0(GEN x,long flag)
    GEN     qfbimagsolvep(GEN Q, GEN n)
    GEN     qfbred0(GEN x, long flag, GEN D, GEN isqrtD, GEN sqrtD)
    GEN     qfbsolve(GEN Q, GEN n)
    GEN     qfi(GEN x, GEN y, GEN z)
    GEN     qfr(GEN x, GEN y, GEN z, GEN d)
    GEN     quaddisc(GEN x)
    GEN     racine(GEN a)
    GEN     redimag(GEN x)
    GEN     redreal(GEN x)
    GEN     redrealnod(GEN x, GEN isqrtD)
    GEN     regula(GEN x, long prec)
    GEN     rhoreal(GEN x)
    GEN     rhorealnod(GEN x, GEN isqrtD)
    GEN     seq_umul(ulong a, ulong b)
    GEN     sqcompimag(GEN x)
    GEN     sqcompreal(GEN x)
    ulong   Fl_sqrt(ulong a, ulong p)
    ulong   Fl_gener_fact(ulong p, GEN fa)
    ulong   Fl_gener(ulong p)
    GEN     znstar(GEN x)

    # arith2.c

    GEN     addprimes(GEN primes)
    GEN     auxdecomp(GEN n, long all)
    long    bigomega(GEN n)
    GEN     binaire(GEN x)
    long    bittest(GEN x, long n)
    GEN     boundfact(GEN n, long lim)
    GEN     core(GEN n)
    GEN     corepartial(GEN n, long l)
    GEN     core0(GEN n,long flag)
    GEN     core2(GEN n)
    GEN     core2partial(GEN n, long l)
    GEN     coredisc(GEN n)
    GEN     coredisc0(GEN n,long flag)
    GEN     coredisc2(GEN n)
    GEN     decomp(GEN n)
    GEN     decomp_small(long n)
    GEN     divisors(GEN n)
    GEN     factorint(GEN n, long flag)
    GEN     gbigomega(GEN n)
    GEN     gbitand(GEN x, GEN y)
    GEN     gbitneg(GEN x, long n)
    GEN     gbitnegimply(GEN x, GEN y)
    GEN     gbitor(GEN x, GEN y)
    GEN     gbittest(GEN x, GEN n)
    GEN	    gbittest3(GEN x, GEN n, long c)
    GEN     gbitxor(GEN x, GEN y)
    GEN     gboundfact(GEN n, long lim)
    GEN     gissquarefree(GEN x)
    GEN     gmu(GEN n)
    GEN     gnumbdiv(GEN n)
    GEN     gomega(GEN n)
    GEN     gphi(GEN n)
    GEN     gsumdiv(GEN n)
    GEN     gsumdivk(GEN n,long k)
    char*   initprimes(ulong maxnum)
    long    issquarefree(GEN x)
    ulong   maxprime()
    void    maxprime_check(ulong c)
    long    mu(GEN n)
    GEN     numbdiv(GEN n)
    long    omega(GEN n)
    GEN     phi(GEN n)
    GEN     prime(long n)
    GEN     primepi(GEN x)
    GEN     primes(long n)
    GEN     removeprimes(GEN primes)
    GEN     smallfact(GEN n)
    GEN     sumdiv(GEN n)
    GEN     sumdivk(GEN n,long k)

    # base1.c

    GEN     bnfnewprec(GEN nf, long prec)
    GEN     bnrnewprec(GEN bnr, long prec)
    void    check_pol_int(GEN x, char *s)
    GEN     check_units(GEN x, char *f)
    void    checkbid(GEN bid)
    GEN     checkbnf(GEN bnf)
    GEN     checkbnf_discard(GEN bnf)
    void    checkbnr(GEN bnr)
    void    checkbnrgen(GEN bnr)
    void    checkid(GEN x, long N)
    GEN     checknf(GEN nf)
    GEN     checknfelt_mod(GEN nf, GEN x, char *s)
    void    checkprimeid(GEN bid)
    void    checkrnf(GEN rnf)
    GEN     galois(GEN x, long prec)
    GEN     galoisapply(GEN nf, GEN aut, GEN x)
    GEN     get_bnf(GEN x,int *t)
    GEN     get_bnfpol(GEN x, GEN *bnf, GEN *nf)
    GEN     get_nf(GEN x,int *t)
    GEN     get_nfpol(GEN x, GEN *nf)
    GEN     get_primeid(GEN x)
    GEN     glambdak(GEN nfz, GEN s, long prec)
    int     gpolcomp(GEN p1, GEN p2)
    GEN     gsmith(GEN x)
    GEN     gsmith2(GEN x)
    GEN     gzetak(GEN nfz, GEN s, long prec)
    GEN     gzetakall(GEN nfz, GEN s, long flag, long prec)
    GEN     initalg(GEN x, long prec)
    GEN     initalgred(GEN x, long prec)
    GEN     initalgred2(GEN x, long prec)
    GEN     initzeta(GEN pol, long prec)
    GEN     mathnf0(GEN x,long flag)
    GEN     matsnf0(GEN x,long flag)
    long    nf_get_r1(GEN nf)
    long    nf_get_r2(GEN nf)
    void    nf_get_sign(GEN nf, long *r1, long *r2)
    long    nfgetprec(GEN x)
    GEN     nfinit0(GEN x,long flag, long prec)
    GEN     nfnewprec(GEN nf, long prec)
    GEN     nfnewprec_i(GEN nf, long prec)
    GEN     rootsof1(GEN x)
    GEN     tschirnhaus(GEN x)

    # base2.c

    GEN     allbase(GEN f, int flag, GEN *dx, GEN *dK, GEN *index, GEN *ptw)
    GEN     base(GEN x, GEN *y)
    GEN     base2(GEN x, GEN *y)
    void    checkmodpr(GEN modpr)
    GEN     compositum(GEN pol1, GEN pol2)
    GEN     compositum2(GEN pol1, GEN pol2)
    GEN     discf(GEN x)
    GEN     discf2(GEN x)
    GEN     factoredbase(GEN x, GEN p, GEN *y)
    GEN     factoreddiscf(GEN x, GEN p)
    GEN     ff_to_nf(GEN x, GEN modpr)
    GEN     fix_relative_pol(GEN nf, GEN x, int chk_lead)
    GEN     gcdpm(GEN f1,GEN f2,GEN pm)
    long    idealval(GEN nf,GEN ix,GEN vp)
    GEN     idealprodprime(GEN nf, GEN L)
    GEN     indexpartial(GEN P, GEN DP)
    GEN     modprX(GEN x, GEN nf,GEN modpr)
    GEN     modprX_lift(GEN x, GEN modpr)
    GEN     modprM(GEN z, GEN nf,GEN modpr)
    GEN     modprM_lift(GEN z, GEN modpr)
    GEN     nf_to_ff_init(GEN nf, GEN *pr, GEN *T, GEN *p)
    GEN     nf_to_ff(GEN nf, GEN x, GEN modpr)
    GEN     nfbasis(GEN x, GEN *y,long flag,GEN p)
    GEN     nfbasis0(GEN x,long flag,GEN p)
    GEN     nfdiscf0(GEN x,long flag, GEN p)
    GEN     nfreducemodideal(GEN nf,GEN x,GEN ideal)
    GEN     nfreducemodpr(GEN nf, GEN x, GEN modpr)
    GEN     polcompositum0(GEN pol1, GEN pol2,long flag)
    GEN     primedec(GEN nf,GEN p)
    GEN     rnfbasis(GEN bnf, GEN order)
    GEN     rnfdet(GEN nf, GEN order)
    GEN     rnfdet2(GEN nf, GEN A, GEN I)
    GEN     rnfdiscf(GEN nf, GEN pol)
    GEN     rnfequation(GEN nf, GEN pol2)
    GEN     rnfequation0(GEN nf, GEN pol2, long flall)
    GEN     rnfequation2(GEN nf, GEN pol)
    GEN     rnfhermitebasis(GEN bnf, GEN order)
    long    rnfisfree(GEN bnf, GEN order)
    GEN     rnflllgram(GEN nf, GEN pol, GEN order,long prec)
    GEN     rnfpolred(GEN nf, GEN pol, long prec)
    GEN     rnfpolredabs(GEN nf, GEN pol, long flag)
    GEN     rnfpseudobasis(GEN nf, GEN pol)
    GEN     rnfsimplifybasis(GEN bnf, GEN order)
    GEN     rnfsteinitz(GEN nf, GEN order)
    GEN     smallbase(GEN x, GEN *y)
    GEN     smalldiscf(GEN x)
    long    val_fact(ulong n, ulong p)
    GEN     zk_to_ff_init(GEN nf, GEN *pr, GEN *T, GEN *p)
    GEN     zk_to_ff(GEN x, GEN modpr)
    GEN     zkmodprinit(GEN nf, GEN pr)

    # base3.c

    GEN     _algtobasis(GEN nf, GEN x)
    GEN     _algtobasis_cp(GEN nf, GEN x)
    GEN     _basistoalg(GEN nf, GEN x)
    GEN     algtobasis(GEN nf, GEN x)
    GEN     algtobasis_i(GEN nf,GEN x)
    GEN     arch_to_perm(GEN arch)
    GEN     basistoalg(GEN nf, GEN x)
    GEN     element_div(GEN nf, GEN x, GEN y)
    GEN     element_inv(GEN nf, GEN x)
    GEN     element_invmodideal(GEN nf, GEN x, GEN ideal)
    GEN     element_mul(GEN nf,GEN x,GEN y)
    GEN     element_muli(GEN nf,GEN x,GEN y)
    GEN     element_mulid(GEN nf, GEN x, long i)
    GEN     element_mulvec(GEN nf, GEN x, GEN v)
    GEN     element_pow(GEN nf,GEN x,GEN k)
    GEN     element_pow_mod_p(GEN nf, GEN x, GEN n, GEN p)
    GEN     element_powmodideal(GEN nf,GEN x,GEN k,GEN ideal)
    GEN     element_powmodidele(GEN nf,GEN x,GEN k,GEN idele,GEN structarch)
    GEN     element_sqr(GEN nf,GEN x)
    GEN     element_sqri(GEN nf, GEN x)
    long    element_val(GEN nf, GEN x, GEN vp)
    GEN     eltmul_get_table(GEN nf, GEN x)
    GEN     ideallist(GEN nf,long bound)
    GEN     ideallist0(GEN nf,long bound, long flag)
    GEN     ideallistarch(GEN nf, GEN list, GEN arch)
    GEN     ideallistarch0(GEN nf, GEN list, GEN arch,long flag)
    GEN     ideallistarchgen(GEN nf, GEN list, GEN arch)
    GEN     ideallistunit(GEN nf,long bound)
    GEN     ideallistunitarch(GEN bnf,GEN list,GEN arch)
    GEN     ideallistunitarchgen(GEN bnf,GEN list,GEN arch)
    GEN     ideallistunitgen(GEN nf,long bound)
    GEN     ideallistzstar(GEN nf,long bound)
    GEN     ideallistzstargen(GEN nf,long bound)
    GEN     idealstar0(GEN nf, GEN x,long flag)
    int     isnfscalar(GEN x)
    GEN     lllreducemodmatrix(GEN x,GEN y)
    GEN     nfdiveuc(GEN nf, GEN a, GEN b)
    GEN     nfdivrem(GEN nf, GEN a, GEN b)
    GEN     nfmod(GEN nf, GEN a, GEN b)
    GEN     nfreducemodidele(GEN nf,GEN g,GEN idele,GEN structarch)
    GEN     reducemodinvertible(GEN x, GEN y)
    GEN     reducemodmatrix(GEN x, GEN y)
    GEN     reducemodHNF(GEN x, GEN y, GEN *Q)
    GEN     set_sign_mod_idele(GEN nf, GEN x, GEN y, GEN idele, GEN sarch)
    GEN     smithrel(GEN H, GEN *newU, GEN *newUi)
    GEN     vecmodii(GEN a, GEN b)
    GEN     zarchstar(GEN nf,GEN x,GEN arch)
    GEN     zideallog(GEN nf,GEN x,GEN bigideal)
    GEN     zidealstar(GEN nf, GEN x)
    GEN     zidealstarinit(GEN nf, GEN x)
    GEN     zidealstarinitall(GEN nf, GEN x,long flun)
    GEN     zidealstarinitgen(GEN nf, GEN x)
    GEN     znlog(GEN x, GEN g)
    GEN     zsigne(GEN nf,GEN alpha,GEN arch)
    GEN     zsigns(GEN nf,GEN alpha)

    # base4.c

    int     Z_ishnfall(GEN x)
    GEN     element_divmodpr(GEN nf, GEN x, GEN y, GEN modpr)
    GEN     element_invmodpr(GEN nf, GEN y, GEN modpr)
    GEN     element_mulmodpr(GEN nf, GEN x, GEN y, GEN modpr)
    GEN     element_powmodpr(GEN nf, GEN x, GEN k, GEN modpr)
    GEN     element_reduce(GEN nf, GEN x, GEN ideal)
    GEN     ideal_two_elt(GEN nf, GEN ix)
    GEN     ideal_two_elt0(GEN nf, GEN ix, GEN a)
    GEN     ideal_two_elt2(GEN nf, GEN x, GEN a)
    GEN     idealadd(GEN nf, GEN x, GEN y)
    GEN     idealaddmultoone(GEN nf, GEN list)
    GEN     idealaddtoone(GEN nf, GEN x, GEN y)
    GEN     idealaddtoone0(GEN nf, GEN x, GEN y)
    GEN     idealappr(GEN nf, GEN x)
    GEN     idealappr0(GEN nf, GEN x, long fl)
    GEN     idealapprfact(GEN nf, GEN x)
    GEN     idealchinese(GEN nf, GEN x, GEN y)
    GEN     idealcoprime(GEN nf, GEN x, GEN y)
    GEN     idealdiv(GEN nf, GEN x, GEN y)
    GEN     idealdiv0(GEN nf, GEN x, GEN y,long flag)
    GEN     idealdivexact(GEN nf, GEN x, GEN y)
    GEN     idealdivpowprime(GEN nf, GEN x, GEN vp, GEN n)
    GEN     idealmulpowprime(GEN nf, GEN x, GEN vp, GEN n)
    GEN     idealfactor(GEN nf, GEN x)
    GEN     idealhermite(GEN nf, GEN x)
    GEN     idealhermite2(GEN nf, GEN a, GEN b)
    GEN     idealhnf0(GEN nf, GEN a, GEN b)
    GEN     idealintersect(GEN nf, GEN x, GEN y)
    GEN     idealinv(GEN nf, GEN ix)
    GEN     ideallllred(GEN nf,GEN ix,GEN vdir,long prec)
    GEN     idealred_elt(GEN nf, GEN I)
    GEN     ideallllred_elt(GEN nf, GEN I, GEN vdir)
    GEN     idealmul(GEN nf, GEN ix, GEN iy)
    GEN     idealmul0(GEN nf, GEN ix, GEN iy, long flag, long prec)
    GEN     idealmulh(GEN nf, GEN ix, GEN iy)
    GEN     idealmulprime(GEN nf,GEN ix,GEN vp)
    GEN     idealmulred(GEN nf, GEN ix, GEN iy, long prec)
    GEN     idealnorm(GEN nf, GEN x)
    GEN     idealpow(GEN nf, GEN ix, GEN n)
    GEN     idealpow0(GEN nf, GEN ix, GEN n, long flag, long prec)
    GEN     idealpowred(GEN nf, GEN ix, GEN n,long prec)
    GEN     idealpows(GEN nf, GEN ideal, long iexp)
    long    idealtyp(GEN *ideal, GEN *arch)
    GEN     ideleaddone(GEN nf, GEN x, GEN idele)
    int     ishnfall(GEN x)
    int     isidentity(GEN x)
    GEN     hnfall_i(GEN A, GEN *ptB, long remove)
    long    isideal(GEN nf,GEN x)
    long    isinvector(GEN v, GEN x, long n)
    GEN     minideal(GEN nf,GEN ix,GEN vdir,long prec)
    GEN     mul_content(GEN cx, GEN cy)
    GEN     nfdetint(GEN nf,GEN pseudo)
    GEN     nfhermite(GEN nf, GEN x)
    GEN     nfhermitemod(GEN nf, GEN x, GEN detmat)
    GEN     nfkermodpr(GEN nf, GEN x, GEN modpr)
    GEN     nfmodprinit(GEN nf, GEN pr)
    GEN     nfsmith(GEN nf, GEN x)
    GEN     nfsolvemodpr(GEN nf, GEN a, GEN b, GEN modpr)
    GEN     prime_to_ideal(GEN nf, GEN vp)
    GEN     principalideal(GEN nf, GEN a)
    GEN     principalidele(GEN nf, GEN a, long prec)
    GEN     vecdiv(GEN x, GEN y)
    GEN     vecinv(GEN x)
    GEN     vecmul(GEN x, GEN y)
    GEN     vecpow(GEN x, GEN n)

    # base5.c

    GEN     lift_to_pol(GEN x)
    GEN     matalgtobasis(GEN nf, GEN x)
    GEN     matbasistoalg(GEN nf, GEN x)
    GEN     rnfalgtobasis(GEN rnf, GEN x)
    GEN     rnfbasistoalg(GEN rnf, GEN x)
    GEN     rnfelementabstorel(GEN rnf, GEN x)
    GEN     rnfelementdown(GEN rnf, GEN x)
    GEN     rnfelementreltoabs(GEN rnf, GEN x)
    GEN     rnfelementup(GEN rnf, GEN x)
    GEN     rnfidealabstorel(GEN rnf, GEN x)
    GEN     rnfidealdown(GEN rnf, GEN x)
    GEN     rnfidealhermite(GEN rnf, GEN x)
    GEN     rnfidealmul(GEN rnf,GEN x,GEN y)
    GEN     rnfidealnormabs(GEN rnf, GEN x)
    GEN     rnfidealnormrel(GEN rnf, GEN x)
    GEN     rnfidealreltoabs(GEN rnf, GEN x)
    GEN     rnfidealtwoelement(GEN rnf,GEN x)
    GEN     rnfidealup(GEN rnf, GEN x)
    GEN     rnfinitalg(GEN nf,GEN pol,long prec)

    # bibli1.c

    GEN     ZM_zc_mul(GEN x, GEN y)
    GEN     ZM_zm_mul(GEN x, GEN y)
    GEN     T2_from_embed(GEN x, long r1)
    GEN     algdep(GEN x, long n, long prec)
    GEN     algdep0(GEN x, long n, long bit,long prec)
    GEN     algdep2(GEN x, long n, long bit)
    GEN     factoredpolred(GEN x, GEN p)
    GEN     factoredpolred2(GEN x, GEN p)
    GEN     kerint(GEN x)
    GEN     kerint1(GEN x)
    GEN     lindep(GEN x, long prec)
    GEN     lindep0(GEN x, long flag,long prec)
    GEN     lindep2(GEN x, long bit)
    GEN     lll(GEN x, long prec)
    GEN     lllgen(GEN x)
    GEN     lllgram(GEN x, long prec)
    GEN     lllgramall(GEN x, long alpha, long flag)
    GEN     lllgramgen(GEN x)
    GEN     lllgramint(GEN x)
    GEN     lllgramintern(GEN x, long alpha, long flag, long prec)
    GEN     lllgramkerim(GEN x)
    GEN     lllgramkerimgen(GEN x)
    GEN     lllint(GEN x)
    GEN     lllint_i(GEN x, long alpha, int gram, GEN *h, GEN *ptfl, GEN *ptB)
    GEN     lllint_ip(GEN x, long alpha)
    GEN     lllintern(GEN x, long D, long flag, long prec)
    GEN     lllintpartial(GEN mat)
    GEN     lllintpartial_ip(GEN mat)
    GEN     lllkerim(GEN x)
    GEN     lllkerimgen(GEN x)
    GEN     matkerint0(GEN x,long flag)
    GEN     minim(GEN a, GEN borne, GEN stockmax)
    GEN     nf_get_LLL(GEN nf)
    GEN     qfrep0(GEN a, GEN borne, long flag)
    GEN     qfminim0(GEN a, GEN borne, GEN stockmax,long flag, long prec)
    GEN     minim2(GEN a, GEN borne, GEN stockmax)
    GEN     ordred(GEN x)
    GEN     perf(GEN a)
    GEN     polred(GEN x)
    GEN     polred0(GEN x, long flag, GEN p)
    GEN     polred2(GEN x)
    GEN     polredabs(GEN x)
    GEN     polredabs0(GEN x, long flag)
    GEN     polredabs2(GEN x)
    GEN     polredabsall(GEN x, long flun)
    GEN     qflll0(GEN x, long flag, long prec)
    GEN     qflllgram0(GEN x, long flag, long prec)
    GEN     smallpolred(GEN x)
    GEN     smallpolred2(GEN x)
    char    *stackmalloc(size_t N)
    GEN     zncoppersmith(GEN P0, GEN N, GEN X, GEN B)

    # bibli2.c

    GEN     binome(GEN x, long k)
    int     cmp_pol(GEN x, GEN y)
    int     cmp_prime_ideal(GEN x, GEN y)
    int     cmp_prime_over_p(GEN x, GEN y)
    int     cmp_vecint(GEN x, GEN y)
    GEN     convol(GEN x, GEN y)
    GEN     cyclo(long n, long v)
    GEN     dirdiv(GEN x, GEN y)
    GEN     dirmul(GEN x, GEN y)
    GEN     dirzetak(GEN nf, GEN b)
    long    gen_search(GEN x, GEN y, int flag, int (*cmp)(GEN,GEN))
    GEN     gen_setminus(GEN set1, GEN set2, int (*cmp)(GEN,GEN))
    GEN     gen_sort(GEN x, int flag, int (*cmp)(GEN,GEN))
    GEN     genrand(GEN N)
    GEN     getheap()
    long    getrand()
    long    getstack()
    long    gettime()
    GEN     gprec(GEN x, long l)
    GEN     gprec_trunc(GEN x, long pr)
    GEN     gprec_w(GEN x, long pr)
    GEN     ggrando(GEN x, long n)
    GEN     gtor(GEN x, long l)
    GEN     gtoset(GEN x)
    GEN     indexlexsort(GEN x)
    GEN     indexsort(GEN x)
    GEN     laplace(GEN x)
    GEN     legendre(long n, long v)
    GEN     lexsort(GEN x)
    GEN     mathilbert(long n)
    GEN     matqpascal(long n, GEN q)
    GEN     modreverse_i(GEN a, GEN T)
    GEN     numtoperm(long n, GEN x)
    int     pari_compare_int(int *a,int *b)
    int     pari_compare_long(long *a,long *b)
    GEN     permtonum(GEN x)
    GEN     polint(GEN xa, GEN ya, GEN x, GEN *dy)
    GEN     polrecip(GEN x)
    GEN     polymodrecip(GEN x)
    GEN     roots_to_pol(GEN a, long v)
    GEN     setintersect(GEN x, GEN y)
    long    setisset(GEN x)
    GEN     setminus(GEN x, GEN y)
    long    setrand(long seed)
    long    setsearch(GEN x, GEN y, long flag)
    GEN     setunion(GEN x, GEN y)
    GEN     sindexlexsort(GEN x)
    GEN     sindexsort(GEN x)
    GEN     sort(GEN x)
    long    tablesearch(GEN T, GEN x, int (*cmp)(GEN,GEN))
    GEN     tayl(GEN x, long v, long precdl)
    GEN     tchebi(long n, long v)
    GEN     vecbinome(long n)
    GEN     vecsort(GEN x, GEN k)
    GEN     vecsort0(GEN x, GEN k, long flag)

    # buch1.c

    GEN     buchimag(GEN D, GEN gcbach, GEN gcbach2, GEN gCO)
    GEN     buchreal(GEN D, GEN gsens, GEN gcbach, GEN gcbach2, GEN gRELSUP, long prec)
    GEN     cgetalloc(long t, size_t l)
    GEN     quadclassunit0(GEN x, long flag,GEN data, long prec)
    GEN     quadhilbert(GEN D, GEN flag, long prec)
    GEN     quadray(GEN bnf, GEN f, GEN flag, long prec)


    # buch2.c

    GEN     bnfclassunit0(GEN P,long flag,GEN data,long prec)
    GEN     bnfinit0(GEN P,long flag,GEN data,long prec)
    GEN     bnfmake(GEN sbnf,long prec)
    GEN     buchall(GEN P, double bach, double bach2, long nbrelpid, long flun, long prec)
    GEN     buchfu(GEN bignf)
    GEN     check_and_build_obj(GEN S, int tag, GEN (*build)(GEN))
    GEN     classgrouponly(GEN P,GEN data,long prec)
    GEN     isprincipal(GEN bignf, GEN x)
    GEN     isprincipalall(GEN bignf, GEN x,long flall)
    GEN     isprincipalfact(GEN bnf,GEN P, GEN e, GEN C, long flag)
    GEN     isprincipalforce(GEN bignf,GEN x)
    GEN     isprincipalgen(GEN bignf, GEN x)
    GEN     isprincipalgenforce(GEN bignf,GEN x)
    GEN     isunit(GEN bignf, GEN x)
    GEN     quick_isprincipalgen(GEN bnf, GEN x)
    GEN     regulator(GEN P,GEN data,long prec)
    GEN     signunits(GEN bignf)
    GEN     smallbuchinit(GEN pol,double bach,double bach2,long nbrelpid,long prec)
    GEN     zsignunits(GEN bnf, GEN archp, int add_zu)

    # buch3.c

    GEN     bnrclass0(GEN bignf, GEN ideal, long flag)
    GEN     bnrconductor(GEN arg0,GEN arg1,GEN arg2,GEN flag)
    GEN     bnrconductorofchar(GEN bnr,GEN chi)
    GEN     bnrdisc0(GEN arg0, GEN arg1, GEN arg2, long flag)
    GEN     bnrdisclist0(GEN bnf,GEN borne, GEN arch, long all)
    GEN     bnrinit0(GEN bignf,GEN ideal,long flag)
    long    bnrisconductor(GEN arg0,GEN arg1,GEN arg2)
    GEN     buchnarrow(GEN bignf)
    GEN     buchray(GEN bignf,GEN ideal)
    GEN     buchrayinit(GEN bignf,GEN ideal)
    GEN     buchrayinitgen(GEN bignf,GEN ideal)
    long    certifybuchall(GEN bnf)
    GEN     conductor(GEN bnr,GEN subgroup,long all)
    GEN     decodemodule(GEN nf, GEN fa)
    GEN     discrayabs(GEN bnr,GEN subgroup)
    GEN     discrayabscond(GEN bnr,GEN subgroup)
    GEN     discrayabslist(GEN bnf,GEN listes)
    GEN     discrayabslistarch(GEN bnf, GEN arch, long bound)
    GEN     discrayabslistlong(GEN bnf, long bound)
    GEN     discrayrel(GEN bnr,GEN subgroup)
    GEN     discrayrelcond(GEN bnr,GEN subgroup)
    GEN     idealmodidele(GEN bnr, GEN x)
    GEN     isprincipalray(GEN bignf, GEN x)
    GEN     isprincipalrayall(GEN bignf, GEN x,long flall)
    GEN     isprincipalraygen(GEN bignf, GEN x)
    GEN     rayclassno(GEN bignf,GEN ideal)
    GEN     rayclassnolist(GEN bnf,GEN listes)
    GEN     rnfconductor(GEN bnf, GEN polrel, long flag)
    GEN     rnfkummer(GEN bnr, GEN subgroup, long all, long prec)
    GEN     rnfnormgroup(GEN bnr, GEN polrel)
    GEN     subgrouplist0(GEN bnr, GEN indexbound, long all)

    # buch4.c

    GEN     bnfisnorm(GEN bnf,GEN x,long flag,long PREC)
    GEN     rnfisnorm(GEN S, GEN x, long flag)
    GEN     rnfisnorminit(GEN bnf, GEN relpol, int galois)
    GEN     bnfissunit(GEN bnf,GEN suni,GEN x)
    GEN     bnfsunit(GEN bnf,GEN s,long PREC)
    long    nfhilbert(GEN bnf,GEN a,GEN b)
    long    nfhilbert0(GEN bnf,GEN a,GEN b,GEN p)
    long    nfhilbertp(GEN bnf,GEN a,GEN b,GEN p)
    long    qpsoluble(GEN pol,GEN p)
    long    qpsolublenf(GEN bnf,GEN pol,GEN p)
    long    zpsoluble(GEN pol,GEN p)
    long    zpsolublenf(GEN bnf,GEN pol,GEN p)

    # elliptic.c

    GEN     addell(GEN e, GEN z1, GEN z2)
    GEN     akell(GEN e, GEN n)
    GEN     anell(GEN e, long n)
    GEN     apell(GEN e, GEN p)
    GEN     apell2(GEN e, GEN p)
    GEN     bilhell(GEN e, GEN z1, GEN z2, long prec)
    GEN     coordch(GEN e, GEN ch)
    GEN     ellap0(GEN e, GEN p, long flag)
    GEN     elleisnum(GEN om, long k, long flag, long prec)
    GEN     elleta(GEN om, long prec)
    GEN     ellheight0(GEN e, GEN a, long flag,long prec)
    GEN     ellinit0(GEN x,long flag,long prec)
    GEN     ellminimalmodel(GEN E, GEN *ptv)
    long    ellrootno(GEN e, GEN p)
    GEN     ellsigma(GEN om, GEN z, long flag, long prec)
    GEN     elltors0(GEN e, long flag)
    GEN     ellwp0(GEN e, GEN z, long flag, long prec, long PREC)

    GEN     ellzeta(GEN om, GEN z, long prec)
    GEN     ghell(GEN e, GEN a, long prec)
    GEN     ghell2(GEN e, GEN a, long prec)
    GEN     globalreduction(GEN e1)
    GEN     initell(GEN x, long prec)
    GEN     elllocalred(GEN e, GEN p1)
    GEN     lseriesell(GEN e, GEN s, GEN A, long prec)
    GEN     mathell(GEN e, GEN x, long prec)
    int     oncurve(GEN e, GEN z)
    GEN     ordell(GEN e, GEN x, long prec)
    GEN     orderell(GEN e, GEN p)
    GEN     pointch(GEN x, GEN ch)
    GEN     pointell(GEN e, GEN z, long prec)
    GEN     powell(GEN e, GEN z, GEN n)
    GEN     smallinitell(GEN x)
    GEN     subell(GEN e, GEN z1, GEN z2)
    GEN     taniyama(GEN e)
    GEN     torsell(GEN e)
    GEN     weipell(GEN e, long precdl)
    GEN     zell(GEN e, GEN z, long prec)

    # es.c
    GEN     GENtocanonicalstr(GEN x)
    GEN     GENtoGENstr(GEN x)
    char*   GENtoTeXstr(GEN x)
    char*   GENtostr(GEN x)
    GEN     Str(GEN g)
    GEN     Strchr(GEN g)
    GEN     Strexpand(GEN g)
    GEN     Strtex(GEN g)
    void    brute(GEN g, char format, long dec)
    void    bruteall(GEN g, char f, long d, long sp)
    void    bruterr(GEN x,char format,long dec)
    void    error0(GEN g)
    void    etatpile(unsigned int n)
    char*   expand_tilde(char *s)
    int     file_is_binary(FILE *f)
    void    flusherr()
    void    fprintferr(char* pat, ...)
    void    killallfiles(int check)
    int     killfile(pariFILE *f)
    GEN     lisGEN(FILE *f)
    void    matbrute(GEN g, char format, long dec)
    pariFILE* newfile(FILE *f, char *name, int type)
    void    os_close(long fd)
    char*   os_getenv(char *s)
    long    os_open(char *s, int mode)
    void    os_read(long fd, char ch[], long s)
    void    (*os_signal(int sig, void (*f)(int)))(int)
    void    outbeaut(GEN x)
    void    outbeauterr(GEN x)
    void    outbrute(GEN x)
    void    outerr(GEN x)
    void    outmat(GEN x)
    void    output(GEN x)
    void    outsor(GEN x)
    void    outtex(GEN x)
    char*   pGENtostr(GEN g, long flag)
    void    pari_fclose(pariFILE *f)
    pariFILE*   pari_fopen(char *s, char *mode)
    pariFILE*   pari_safefopen(char *s, char *mode)
    char*   pari_strdup(char *s)
    char*   pari_strndup(char *s, long n)
    char*   pari_unique_filename(char *s)
    void    pari_unlink(char *s)
    void    pariflush()
    void    pariputc(char c)
    void    pariputs(char *s)
    void    pariputsf(char *format, ...)
    int     popinfile()
    #void    print(GEN g)   # syntax error
    void    print1(GEN g)
    void    printp(GEN g)
    void    printp1(GEN g)
    void    printtex(GEN g)
    GEN     readbin(char *name, FILE *f)
    void    sor(GEN g, char fo, long dd, long chmp)
    void    switchin(char *name)
    void    switchout(char *name)
    void    texe(GEN g, char format, long dec)
    pariFILE* try_pipe(char *cmd, int flag)
    char*   type_name(long t)
    void    voir(GEN x, long nb)
    #void    vpariputs(char* format, va_list args)
    void    write0(char *s, GEN g)
    void    write1(char *s, GEN g)
    void    writebin(char *name, GEN x)
    void    writetex(char *s, GEN g)

    # galconj.c

    GEN     checkgal(GEN gal)
    GEN     galoisconj(GEN nf)
    GEN     galoisconj0(GEN nf, long flag, GEN d, long prec)
    GEN     galoisconj2(GEN x, long nbmax, long prec)
    GEN     galoisconj4(GEN T, GEN den, long flag, long karma)
    GEN     galoisexport(GEN gal, long format)
    GEN     galoisfixedfield(GEN gal, GEN v, long flag, long y)
    GEN     galoisidentify(GEN gal)
    GEN     galoisinit(GEN nf, GEN den, long karma)
    GEN     galoisisabelian(GEN gal, long flag)
    GEN     galoispermtopol(GEN gal, GEN perm)
    GEN     galoissubgroups(GEN G)
    GEN     galoissubfields(GEN G, long flag, long v)
    long    numberofconjugates(GEN T, long pdepart)
    GEN     vandermondeinverse(GEN L, GEN T, GEN den, GEN prep)
    # gen1.c

    GEN     gadd(GEN x, GEN y)
    GEN     gaddsg(long x, GEN y)
    GEN     gdiv(GEN x, GEN y)
    GEN     gmul(GEN x, GEN y)
    GEN     gsqr(GEN x)
    GEN     gsub(GEN x, GEN y)

    # gen2.c
    GEN     gopgs2(GEN (*f)(GEN, GEN), GEN y, long s)
    GEN     gopsg2(GEN (*f)(GEN, GEN), long s, GEN y)
    void    gopsgz(GEN (*f)(GEN, GEN), long s, GEN y, GEN z)
    int     opgs2(int (*f)(GEN, GEN), GEN y, long s)

    long    ZX_valuation(GEN x, GEN *Z)
    GEN     cgetimag()
    GEN     cgetp(GEN x)
    GEN     cvtop(GEN x, GEN p, long l)
    GEN     cvtop2(GEN x, GEN y)
    GEN     from_Kronecker(GEN z, GEN pol)
    GEN     gabs(GEN x, long prec)
    void    gaffect(GEN x, GEN y)
    void    gaffsg(long s, GEN x)
    GEN     gclone(GEN x)
    int     gcmp(GEN x, GEN y)
    int     gcmpsg(long x, GEN y)
    int     gcmp0(GEN x)
    int     gcmp1(GEN x)
    int     gcmp_1(GEN x)
    GEN     gcvtop(GEN x, GEN p, long r)
    int     gegal(GEN x, GEN y)
    int     gegalsg(long s, GEN x)
    long    gexpo(GEN x)
    long    ggval(GEN x, GEN p)
    long    glength(GEN x)
    GEN     gmax(GEN x, GEN y)
    GEN     gmin(GEN x, GEN y)
    GEN     gneg(GEN x)
    GEN     gneg_i(GEN x)
    GEN     greffe(GEN x, long l, long use_stack)
    int     gsigne(GEN x)
    GEN     gtofp(GEN z, long prec)
    GEN     gtolist(GEN x)
    long    gtolong(GEN x)
    int     lexcmp(GEN x, GEN y)
    GEN     listconcat(GEN list1, GEN list2)
    GEN     listcreate(long n)
    GEN     listinsert(GEN list, GEN object, long index)
    void    listkill(GEN list)
    GEN     listput(GEN list, GEN object, long index)
    GEN     listsort(GEN list, long flag)
    GEN     matsize(GEN x)
    GEN     normalize(GEN x)
    GEN     normalizepol(GEN x)
    GEN     normalizepol_i(GEN x, long lx)
    long    polvaluation(GEN x, GEN *z)
    long    polvaluation_inexact(GEN x, GEN *Z)
    GEN     pureimag(GEN x)
    GEN     quadtoc(GEN x, long l)
    long    sizedigit(GEN x)
    long    u_lval(ulong x, ulong p)
    long    u_lvalrem(ulong x, ulong p, ulong *py)
    long    u_pvalrem(ulong x, GEN p, ulong *py)
    GEN     to_Kronecker(GEN P, GEN Q)
    GEN     vecmax(GEN x)
    GEN     vecmin(GEN x)
    long    Z_lval(GEN n, ulong p)
    long    Z_lvalrem(GEN n, ulong p, GEN *py)
    long    z_pval(long n, GEN p)
    long    Z_pval(GEN n, GEN p)
    long    Z_pvalrem(GEN x, GEN p, GEN *py)

    # gen3.c

    GEN     _toser(GEN x)
    GEN     Mod0(GEN x, GEN y,long flag)
    GEN     ceil_safe(GEN x)
    GEN     ceilr(GEN x)
    GEN     centerlift(GEN x)
    GEN     centerlift0(GEN x,long v)
    GEN     coefs_to_col(long n, ...)
    GEN     coefs_to_int(long n, ...)
    GEN     coefs_to_pol(long n, ...)
    GEN     coefs_to_vec(long n, ...)
    GEN     compo(GEN x, long n)
    GEN     deg1pol(GEN x1, GEN x0,long v)
    GEN     deg1pol_i(GEN x1, GEN x0,long v)
    long    degree(GEN x)
    GEN     denom(GEN x)
    GEN     deriv(GEN x, long v)
    GEN     derivpol(GEN x)
    GEN     derivser(GEN x)
    GEN     diviiround(GEN x, GEN y)
    GEN     divrem(GEN x, GEN y, long v)
    GEN     gand(GEN x, GEN y)
    GEN     gceil(GEN x)
    GEN     gcvtoi(GEN x, long *e)
    GEN     gdivent(GEN x, GEN y)
    GEN     gdiventres(GEN x, GEN y)
    GEN     gdivgs(GEN x, long s)
    GEN     gdivmod(GEN x, GEN y, GEN *pr)
    GEN     gdivround(GEN x, GEN y)
    GEN     geq(GEN x, GEN y)
    GEN     geval(GEN x)
    GEN     gfloor(GEN x)
    GEN     gfloor2n(GEN x, long s)
    GEN     gfrac(GEN x)
    GEN     gge(GEN x, GEN y)
    GEN     ggprecision(GEN x)
    GEN     ggt(GEN x, GEN y)
    GEN     gimag(GEN x)
    GEN     ginv(GEN x)
    GEN     gle(GEN x, GEN y)
    GEN     glt(GEN x, GEN y)
    GEN     gmod(GEN x, GEN y)
    GEN     gmodulcp(GEN x,GEN y)
    GEN     gmodulo(GEN x,GEN y)
    GEN     gmodulsg(long x, GEN y)
    GEN     gmodulss(long x, long y)
    GEN     gmul2n(GEN x, long n)
    GEN     gmulsg(long s, GEN y)
    GEN     gne(GEN x, GEN y)
    GEN     gnot(GEN x)
    GEN     gor(GEN x, GEN y)
    GEN     gpolvar(GEN y)
    long    gprecision(GEN x)
    GEN     gram_matrix(GEN M)
    GEN     greal(GEN x)
    GEN     grndtoi(GEN x, long *e)
    GEN     ground(GEN x)
    GEN     gshift(GEN x, long n)
    GEN     gsubst(GEN x, long v, GEN y)
    GEN     gsubstpol(GEN x, GEN v, GEN y)
    GEN     gtocol(GEN x)
    GEN     gtopoly(GEN x, long v)
    GEN     gtopolyrev(GEN x, long v)
    GEN     gtoser(GEN x, long v)
    GEN     gtovec(GEN x)
    GEN     gtovecsmall(GEN x)
    GEN     gtrunc(GEN x)
    int     gvar(GEN x)
    int     gvar2(GEN x)
    int     gvar9(GEN x)
    GEN     hqfeval(GEN q, GEN x)
    GEN     imag_i(GEN x)
    GEN     int2n(long n)
    GEN     integ(GEN x, long v)
    int     iscomplex(GEN x)
    int     isexactzero(GEN g)
    int     isexactzeroscalar(GEN g)
    int     isinexactreal(GEN x)
    long    isint(GEN n, long *ptk)
    int     ismonome(GEN x)
    GEN     lift(GEN x)
    GEN     lift0(GEN x,long v)
    GEN     lift_intern0(GEN x,long v)
    GEN     truncr(GEN x)
    GEN     mulmat_real(GEN x, GEN y)
    GEN     numer(GEN x)
    long    padicprec(GEN x, GEN p)
    GEN     polcoeff0(GEN x,long n,long v)
    GEN     polcoeff_i(GEN x, long n, long v)
    long    poldegree(GEN x,long v)
    GEN     poleval(GEN x, GEN y)
    GEN     pollead(GEN x,long v)
    long    precision(GEN x)
    GEN     precision0(GEN x,long n)
    GEN     qf_base_change(GEN q, GEN M, int flag)
    GEN     qfeval(GEN q, GEN x)
    GEN     real_i(GEN x)
    GEN     real2n(long n, long prec)
    GEN     recip(GEN x)
    GEN     round0(GEN x, GEN *pte)
    GEN     roundr(GEN x)
    GEN     scalarpol(GEN x, long v)
    GEN     scalarser(GEN x, long v, long prec)
    GEN     simplify(GEN x)
    GEN     simplify_i(GEN x)
    GEN     truecoeff(GEN x, long n)
    GEN     trunc0(GEN x, GEN *pte)
    GEN     u2toi(ulong a, ulong b)
    GEN     zerocol(long n)
    GEN     zeromat(long m, long n)
    GEN     zeropadic(GEN p, long e)
    GEN     zeropol(long v)
    GEN     zeroser(long v, long prec)
    GEN     zerovec(long n)

    # groupid.c

    long    group_ident(GEN G, GEN S)

    # ifactor1.c

    long    BSW_psp(GEN N)
    long    is_357_power(GEN x, GEN *pt, ulong *mask)
    long    is_odd_power(GEN x, GEN *pt, ulong *curexp, ulong cutoffbits)
    long    millerrabin(GEN n, long k)
    GEN     nextprime(GEN n)
    GEN     plisprime(GEN N, long flag)
    GEN     precprime(GEN n)

    # init.c

    long    TIMER(pari_timer *T)
    void    TIMERstart(pari_timer *T)
    long    allocatemoremem(size_t newsize)
    GEN     changevar(GEN x, GEN y)
    void    disable_dbg(long val)
    GEN     dummycopy(GEN x)
    void    err(long numerr, ...)
    #void   *err_catch(long errnum, jmp_buf *penv)
    void    err_leave(void **v)
    GEN     forcecopy(GEN x)
    void    freeall()
    GEN     gcopy(GEN x)
    GEN     gcopy_i(GEN x, long lx)
    GEN     gerepile(pari_sp ltop, pari_sp lbot, GEN q)
    void    gerepileall(pari_sp av, int n, ...)
    void    gerepileallsp(pari_sp av, pari_sp tetpil, int n, ...)
    void    gerepilecoeffs(pari_sp av, GEN x, int n)
    void    gerepilecoeffssp(pari_sp av, pari_sp tetpil, long *g, int n)
    GEN     gerepilecopy(pari_sp av, GEN x)
    void    gerepilemany(pari_sp av, GEN* g[], int n)
    void    gerepilemanysp(pari_sp av, pari_sp tetpil, GEN* g[], int n)
    GEN     gerepileupto(pari_sp av, GEN q)
    GEN     gerepileuptoint(pari_sp av, GEN q)
    GEN     gerepileuptoleaf(pari_sp av, GEN q)
    char*   gpmalloc(size_t bytes)
    char*   gprealloc(void *pointer,size_t size)
    void    gunclone(GEN x)
    void    killbloc(GEN x)
    void    msgTIMER(pari_timer *T, char *format, ...)
    void    msgtimer(char *format, ...)
    GEN     newbloc(long n)
    void    pari_init(size_t parisize, ulong maxprime)
    GEN     reorder(GEN x)
    void    stackdummy(GEN x, long l)
    GEN     stackify(GEN x)
    stackzone* switch_stack(stackzone *z, long n)
    long    taille(GEN x)
    long    taille2(GEN x)
    long    timer()
    long    timer2()

    # members.c

    GEN     member_a1(GEN x)
    GEN     member_a2(GEN x)
    GEN     member_a3(GEN x)
    GEN     member_a4(GEN x)
    GEN     member_a6(GEN x)
    GEN     member_area(GEN x)
    GEN     member_b2(GEN x)
    GEN     member_b4(GEN x)
    GEN     member_b6(GEN x)
    GEN     member_b8(GEN x)
    GEN     member_bnf(GEN x)
    GEN     member_c4(GEN x)
    GEN     member_c6(GEN x)
    GEN     member_clgp(GEN x)
    GEN     member_codiff(GEN x)
    GEN     member_cyc(GEN clg)
    GEN     member_diff(GEN x)
    GEN     member_disc(GEN x)
    GEN     member_e(GEN x)
    GEN     member_eta(GEN x)
    GEN     member_f(GEN x)
    GEN     member_fu(GEN x)
    GEN     member_futu(GEN x)
    GEN     member_gen(GEN x)
    GEN     member_group(GEN x)
    GEN     member_index(GEN x)
    GEN     member_j(GEN x)
    GEN     member_mod(GEN x)
    GEN     member_nf(GEN x)
    GEN     member_no(GEN clg)
    GEN     member_omega(GEN x)
    GEN     member_orders(GEN x)
    GEN     member_p(GEN x)
    GEN     member_pol(GEN x)
    GEN     member_reg(GEN x)
    GEN     member_roots(GEN x)
    GEN     member_sign(GEN x)
    GEN     member_t2(GEN x)
    GEN     member_tate(GEN x)
    GEN     member_tufu(GEN x)
    GEN     member_tu(GEN x)
    GEN     member_w(GEN x)
    GEN     member_zk(GEN x)
    GEN     member_zkst(GEN bid)

    # mp.c

    int     absi_cmp(GEN x, GEN y)
    int     absi_equal(GEN x, GEN y)
    int     absr_cmp(GEN x, GEN y)
    GEN     addii_sign(GEN x, long sx, GEN y, long sy)
    GEN     addir_sign(GEN x, long sx, GEN y, long sy)
    GEN     addrr_sign(GEN x, long sx, GEN y, long sy)
    GEN     addsi_sign(long x, GEN y, long sy)
    GEN     addsr(long x, GEN y)
    GEN     addss(long x, long y)
    void    affir(GEN x, GEN y)
    void    affrr(GEN x, GEN y)
    GEN     bezout(GEN a, GEN b, GEN *u, GEN *v)
    long    cbezout(long a,long b,long *uu,long *vv)
    void    cgiv(GEN x)
    int     cmpii(GEN x, GEN y)
    int     cmprr(GEN x, GEN y)
    int     cmpsi(long x, GEN y)
    int     cmpui(ulong x, GEN y)
    GEN     dbltor(double x)
    GEN     diviiexact(GEN x, GEN y)
    GEN     diviuexact(GEN x, ulong y)
    GEN     divir(GEN x, GEN y)
    GEN     divis(GEN y, long x)
    GEN     divis_rem(GEN x, long y, long *rem)
    GEN     diviu_rem(GEN y, ulong x, ulong *rem)
    GEN     divri(GEN x, GEN y)
    GEN     divrr(GEN x, GEN y)
    GEN     divrs(GEN x, long y)
    GEN     divsi(long x, GEN y)
    GEN     divsr(long x, GEN y)
    GEN     dvmdii(GEN x, GEN y, GEN *z)
    int     egalii(GEN x, GEN y)
    GEN     floorr(GEN x)
    GEN     gcdii(GEN x, GEN y)
    GEN     int_normalize(GEN x, long known_zero_words)
    int     invmod(GEN a, GEN b, GEN *res)
    ulong   invrev(ulong b)
    ulong   Fl_inv(ulong x, ulong p)
    GEN     ishiftr(GEN x, long n)
    GEN     modii(GEN x, GEN y)
    void    modiiz(GEN x, GEN y, GEN z)
    void    mpdivz(GEN x, GEN y, GEN z)
    GEN     mulii(GEN x, GEN y)
    GEN     mulir(GEN x, GEN y)
    GEN     mulrr(GEN x, GEN y)
    GEN     mulsi(long x, GEN y)
    GEN     mulsr(long x, GEN y)
    GEN     mulss(long x, long y)
    GEN     mului(ulong x, GEN y)
    GEN     mulur(ulong x, GEN y)
    GEN     muluu(ulong x, ulong y)
    long    pari_rand31()
    GEN     randomi(GEN x)
    int     ratlift(GEN x, GEN m, GEN *a, GEN *b, GEN amax, GEN bmax)
    GEN     resmod2n(GEN x, long n)
    double  rtodbl(GEN x)
    GEN     shifti(GEN x, long n)
    GEN     shifti3(GEN x, long n, long flag)
    GEN     sqri(GEN x)
    #define sqrti(x) sqrtremi((x),NULL)
    GEN     sqrtremi(GEN S, GEN *R)
    GEN     subsr(long x, GEN y)
    GEN     truedvmdii(GEN x, GEN y, GEN *z)
    ulong   umodiu(GEN y, ulong x)
    long    vals(ulong x)

    # nffactor.c

    GEN     nffactor(GEN nf,GEN x)
    GEN     nffactormod(GEN nf,GEN pol,GEN pr)
    int     nfisgalois(GEN nf, GEN x)
    GEN     nfroots(GEN nf,GEN pol)
    GEN     rnfcharpoly(GEN nf,GEN T,GEN alpha,int n)
    GEN     rnfdedekind(GEN nf,GEN T,GEN pr)
    GEN     unifpol(GEN nf,GEN pol,long flag)

    # part.c

    GEN     numbpart(GEN x)

    # perm.c

    GEN     abelian_group(GEN G)
    GEN     bitvec_alloc(long n)
    void    bitvec_clear(GEN bitvec, long b)
    void    bitvec_set(GEN bitvec, long b)
    GEN     bitvec_shorten(GEN bitvec, long n)
    long    bitvec_test(GEN bitvec, long b)
    GEN     cyclicgroup(GEN g, long s)
    GEN     cyclicperm(long l, long d)
    GEN     cyc_pow(GEN cyc, long exp)
    GEN     cyc_pow_perm(GEN cyc, long exp)
    GEN     dicyclicgroup(GEN g1, GEN g2, long s1, long s2)
    GEN     group_abelianHNF(GEN G, GEN L)
    GEN     group_abelianSNF(GEN G, GEN L)
    long    group_domain(GEN G)
    GEN     group_elts(GEN G, long n)
    GEN     group_export(GEN G, long format)
    long    group_isA4S4(GEN G)
    long    group_isabelian(GEN G)
    GEN     group_leftcoset(GEN G, GEN g)
    long    group_order(GEN G)
    long    group_perm_normalize(GEN N, GEN g)
    GEN     group_quotient(GEN G, GEN H)
    GEN     group_rightcoset(GEN G, GEN g)
    GEN     group_subgroups(GEN G)
    GEN     groupelts_center(GEN S)
    GEN     groupelts_abelian_group(GEN S)
    int     perm_commute(GEN p, GEN q)
    GEN     perm_cycles(GEN v)
    GEN     perm_identity(long l)
    GEN     perm_inv(GEN x)
    GEN     perm_mul(GEN s, GEN t)
    long    perm_order(GEN perm)
    GEN     perm_pow(GEN perm, long exp)
    GEN     quotient_group(GEN C, GEN G)
    GEN     quotient_perm(GEN C, GEN p)
    GEN     vec_to_vecsmall(GEN z)
    GEN     vecperm_orbits(GEN v, long n)
    GEN     vecsmall_append(GEN V, long s)
    long    vecsmall_coincidence(GEN u, GEN v)
    GEN     vecsmall_concat(GEN u, GEN v)
    GEN     vecsmall_const(long n, long c)
    GEN     vecsmall_copy(GEN x)
    int     vecsmall_lexcmp(GEN x, GEN y)
    long    vecsmall_pack(GEN V, long base, long mod)
    int     vecsmall_prefixcmp(GEN x, GEN y)
    GEN     vecsmall_prepend(GEN V, long s)
    GEN     vecsmall_shorten(GEN v, long n)
    void    vecsmall_sort(GEN V)
    GEN     vecsmall_to_col(GEN z)
    GEN     vecsmall_to_vec(GEN z)
    GEN     vecsmall_uniq(GEN V)
    GEN     vecvecsmall_indexsort(GEN x)
    GEN     vecvecsmall_sort(GEN x)
    long    vecvecsmall_search(GEN x, GEN y, long flag)

    # polarit1.c

    long    Flx_nbfact(GEN z, ulong p)
    long    Flx_nbroots(GEN f, ulong p)
    GEN     FpX_degfact(GEN f, GEN p)
    long    FpX_is_irred(GEN f, GEN p)
    long    FpX_is_squarefree(GEN f, GEN p)
    long    FpX_is_totally_split(GEN f, GEN p)
    GEN     FpX_factor(GEN f, GEN p)
    long    FpX_nbfact(GEN f, GEN p)
    long    FpX_nbroots(GEN f, GEN p)
    GEN     FpXQX_gcd(GEN P, GEN Q, GEN T, GEN p)
    GEN     FqX_factor(GEN x, GEN T, GEN p)
    GEN     FqX_gcd(GEN P, GEN Q, GEN T, GEN p)
    long    FqX_is_squarefree(GEN P, GEN T, GEN p)
    long    FqX_nbfact(GEN u, GEN T, GEN p)
    long    FqX_nbroots(GEN f, GEN T, GEN p)
    GEN     FpX_rand(long d, long v, GEN p)
    GEN     FpX_roots(GEN f, GEN p)
    GEN     apprgen(GEN f, GEN a)
    GEN     apprgen9(GEN f, GEN a)
    GEN     factcantor(GEN x, GEN p)
    GEN     factmod(GEN f, GEN p)
    GEN     factmod9(GEN f, GEN p, GEN a)
    GEN     factormod0(GEN f, GEN p,long flag)
    GEN     factorpadic0(GEN f,GEN p,long r,long flag)
    GEN     factorpadic2(GEN x, GEN p, long r)
    GEN     factorpadic4(GEN x, GEN p, long r)
    GEN     ffinit(GEN p,long n, long v)
    int     gdvd(GEN x, GEN y)
    long    hensel_lift_accel(long n, long *pmask)
    GEN     init_Fq(GEN p, long n, long v)
    GEN     padicsqrtnlift(GEN a, GEN n, GEN S, GEN p, long e)
    int     poldvd(GEN x, GEN y, GEN *z)
    GEN     poldivrem(GEN x, GEN y, GEN *pr)
    GEN     poldivrem_i(GEN x, GEN y, GEN *pr, long vx)
    GEN     rootmod(GEN f, GEN p)
    GEN     rootmod0(GEN f, GEN p,long flag)
    GEN     rootmod2(GEN f, GEN p)
    GEN     rootpadic(GEN f, GEN p, long r)
    GEN     rootpadicfast(GEN f, GEN p, long e)
    GEN     rootpadiclift(GEN T, GEN S, GEN q, long e)
    GEN     rootpadicliftroots(GEN f, GEN S, GEN q, long e)
    GEN     roots2(GEN pol,long PREC)
    GEN     rootsold(GEN x, long l)
    GEN     simplefactmod(GEN f, GEN p)

    # polarit2.c

    GEN     Newton_exponents(long e)
    GEN     Q_content(GEN x)
    GEN     Q_denom(GEN x)
    GEN     Q_div_to_int(GEN x, GEN c)
    GEN     Q_muli_to_int(GEN x, GEN d)
    GEN     Q_primitive_part(GEN x, GEN *ptc)
    GEN     Q_primpart(GEN x)
    GEN     Q_remove_denom(GEN x, GEN *ptd)
    GEN     RgX_extgcd(GEN x, GEN y, GEN *U, GEN *V)
    GEN     ZX_squff(GEN f, GEN *ex)
    GEN     centermod(GEN x, GEN p)
    GEN     centermod_i(GEN x, GEN p, GEN ps2)
    GEN     centermodii(GEN x, GEN p, GEN po2)
    GEN     concat_factor(GEN f, GEN g)
    GEN     content(GEN x)
    GEN     deg1_from_roots(GEN L, long v)
    GEN     discsr(GEN x)
    GEN     divide_conquer_prod(GEN x, GEN (*mul)(GEN,GEN))
    GEN     factor(GEN x)
    GEN     factor0(GEN x,long flag)
    GEN     factorback(GEN fa,GEN nf)
    GEN     factorback0(GEN fa,GEN e, GEN nf)
    GEN     factorbackelt(GEN fa, GEN e, GEN nf)
    GEN     factpol(GEN x, long hint)
    GEN     gbezout(GEN x, GEN y, GEN *u, GEN *v)
    GEN     gcd0(GEN x, GEN y,long flag)
    GEN     gdeflate(GEN x, long v, long d)
    GEN     gdivexact(GEN x, GEN y)
    GEN     ggcd(GEN x, GEN y)
    GEN     ginvmod(GEN x, GEN y)
    GEN     gisirreducible(GEN x)
    GEN     glcm(GEN x, GEN y)
    GEN     glcm0(GEN x, GEN y)
    GEN     hensel_lift_fact(GEN pol, GEN Q, GEN T, GEN p, GEN pe, long e)
    GEN     leftright_pow(GEN,GEN,void*,GEN (*sqr)(void*,GEN),GEN (*mul)(void*,GEN,GEN))
    GEN     leftright_pow_u(GEN x, ulong n, void *data, GEN (*sqr)(void*,GEN), GEN (*mul)(void*,GEN,GEN))
    long    logint(GEN B, GEN y, GEN *ptq)
    GEN     newtonpoly(GEN x, GEN p)
    GEN     nfgcd(GEN P, GEN Q, GEN nf, GEN den)
    GEN     nfisincl(GEN a, GEN b)
    GEN     nfisisom(GEN a, GEN b)
    GEN     nfrootsQ(GEN x)
    GEN     poldeflate(GEN x0, long *m)
    GEN     poldeflate_i(GEN x0, long d)
    GEN     poldisc0(GEN x, long v)
    GEN     polfnf(GEN a, GEN t)
    GEN     polhensellift(GEN pol, GEN fct, GEN p, long exp)
    GEN     polinflate(GEN x0, long d)
    GEN     polresultant0(GEN x, GEN y,long v,long flag)
    GEN     polsym(GEN x, long n)
    GEN     primitive_part(GEN x, GEN *c)
    GEN     primpart(GEN x)
    GEN     pseudorem(GEN x, GEN y)
    GEN     quadgen(GEN x)
    GEN     quadpoly(GEN x)
    GEN     quadpoly0(GEN x, long v)
    GEN     reduceddiscsmith(GEN pol)
    GEN     resultant2(GEN x, GEN y)
    GEN     resultantducos(GEN x, GEN y)
    GEN     roots_from_deg1(GEN x)
    GEN     sort_factor(GEN y, int (*cmp)(GEN,GEN))
    GEN     sort_factor_gen(GEN y, int (*cmp)(GEN,GEN))
    GEN     sort_vecpol(GEN a, int (*cmp)(GEN,GEN))
    GEN     srgcd(GEN x, GEN y)
    long    sturmpart(GEN x, GEN a, GEN b)
    GEN     subresall(GEN u, GEN v, GEN *sol)
    GEN     subresext(GEN x, GEN y, GEN *U, GEN *V)
    GEN     sylvestermatrix(GEN x,GEN y)
    GEN     vecbezout(GEN x, GEN y)
    GEN     vecbezoutres(GEN x, GEN y)

    # polarit3.c

    ulong   Fl_pow(ulong x, ulong n, ulong p)
    GEN     Fp_pows(GEN A, long k, GEN N)
    GEN     Fp_powu(GEN x, ulong k, GEN p)
    GEN     FpM_red(GEN z, GEN p)
    GEN     FpM_to_mod(GEN z, GEN p)
    GEN     FpV_polint(GEN xa, GEN ya, GEN p)
    GEN     FpV_red(GEN z, GEN p)
    GEN     FpV_roots_to_pol(GEN V, GEN p, long v)
    GEN     FpV_to_mod(GEN z, GEN p)
    GEN     FpX_FpXQ_compo(GEN f,GEN x,GEN T,GEN p)
    GEN     FpX_FpXQV_compo(GEN f,GEN x,GEN T,GEN p)
    GEN     FpX_Fp_add(GEN y,GEN x,GEN p)
    GEN     FpX_Fp_mul(GEN y,GEN x,GEN p)
    GEN     FpX_add(GEN x,GEN y,GEN p)
    GEN     FpX_center(GEN T,GEN mod)
    GEN     FpX_chinese_coprime(GEN x,GEN y,GEN Tx,GEN Ty,GEN Tz,GEN p)
    GEN     FpX_divrem(GEN x, GEN y, GEN p, GEN *pr)
    GEN     FpX_eval(GEN x,GEN y,GEN p)
    GEN     FpX_extgcd(GEN x, GEN y, GEN p, GEN *ptu, GEN *ptv)
    GEN     FpX_factorff_irred(GEN P, GEN Q, GEN l)
    void    FpX_ffintersect(GEN P,GEN Q,long n,GEN l,GEN *SP,GEN *SQ,GEN MA,GEN MB)
    GEN     FpX_ffisom(GEN P,GEN Q,GEN l)
    GEN     FpX_gcd(GEN x, GEN y, GEN p)
    GEN     FpX_mul(GEN x,GEN y,GEN p)
    GEN     FpX_neg(GEN x,GEN p)
    GEN     FpX_normalize(GEN z, GEN p)
    GEN     FpX_red(GEN z, GEN p)
    GEN     FpX_sqr(GEN x,GEN p)
    GEN     FpX_sub(GEN x,GEN y,GEN p)
    GEN     FpX_to_mod(GEN z, GEN p)
    GEN     FpXQ_charpoly(GEN x, GEN T, GEN p)
    GEN     FpXQ_div(GEN x,GEN y,GEN T,GEN p)
    GEN     FpXQ_ffisom_inv(GEN S,GEN Tp, GEN p)
    GEN     FpXQ_inv(GEN x,GEN T,GEN p)
    GEN     FpXQ_invsafe(GEN x, GEN T, GEN p)
    GEN     FpXQ_matrix_pow(long n, long m, GEN y, GEN P, GEN l)
    GEN     FpXQ_minpoly(GEN x, GEN T, GEN p)
    GEN     FpXQ_mul(GEN y,GEN x,GEN T,GEN p)
    GEN     FpXQ_pow(GEN x, GEN n, GEN T, GEN p)
    GEN     FpXQ_powers(GEN x, long l, GEN T, GEN p)
    GEN     FpXQ_sqr(GEN y, GEN T, GEN p)
    GEN     FpXQ_sqrtn(GEN a, GEN n, GEN T, GEN p, GEN *zetan)
    GEN     FpXQX_mul(GEN x, GEN y, GEN T, GEN p)
    GEN     FpXQX_red(GEN z, GEN T, GEN p)
    GEN     FpXQX_sqr(GEN x, GEN T, GEN p)
    GEN     FpXQX_extgcd(GEN x, GEN y, GEN T, GEN p, GEN *ptu, GEN *ptv)
    GEN     FpXQX_divrem(GEN x, GEN y, GEN T, GEN p, GEN *pr)
    GEN     FpXQX_safegcd(GEN P, GEN Q, GEN T, GEN p)
    GEN     FpXQXV_prod(GEN V, GEN Tp, GEN p)
    GEN     FpXQYQ_pow(GEN x, GEN n, GEN S, GEN T, GEN p)
    GEN     FpXV_FpV_innerprod(GEN V, GEN W, GEN p)
    GEN     FpXV_prod(GEN V, GEN p)
    GEN     FpXV_red(GEN z, GEN p)
    GEN     FpXX_red(GEN z, GEN p)
    GEN     FpX_rescale(GEN P, GEN h, GEN p)
    GEN     FpY_FpXY_resultant(GEN a, GEN b0, GEN p)
    GEN     Fq_inv(GEN x, GEN T, GEN p)
    GEN     Fq_invsafe(GEN x, GEN T, GEN p)
    GEN     Fq_add(GEN x, GEN y, GEN T, GEN p)
    GEN     Fq_mul(GEN x, GEN y, GEN T, GEN p)
    GEN     Fq_neg(GEN x, GEN T, GEN p)
    GEN     Fq_neg_inv(GEN x, GEN T, GEN p)
    GEN     Fq_pow(GEN x, GEN n, GEN T, GEN p)
    GEN     Fq_red(GEN x, GEN T, GEN p)
    GEN     Fq_sub(GEN x, GEN y, GEN T, GEN p)
    GEN     FqM_to_FlxM(GEN x, GEN T, GEN pp)
    GEN     FqV_roots_to_pol(GEN V, GEN T, GEN p, long v)
    GEN     FqV_red(GEN z, GEN T, GEN p)
    GEN     FqV_to_FlxC(GEN v, GEN T, GEN pp)
    GEN     FqX_Fq_mul(GEN P, GEN U, GEN T, GEN p)
    GEN     FqX_div(GEN x, GEN y, GEN T, GEN p)
    GEN     FqX_divrem(GEN x, GEN y, GEN T, GEN p, GEN *z)
    GEN     FqX_normalize(GEN z, GEN T, GEN p)
    GEN     FqX_red(GEN z, GEN T, GEN p)
    GEN     FqX_rem(GEN x, GEN y, GEN T, GEN p)
    GEN     FqX_mul(GEN x, GEN y, GEN T, GEN p)
    GEN     FqX_sqr(GEN x, GEN T, GEN p)
    GEN     QXQ_inv(GEN A, GEN B)
    GEN     ZX_caract(GEN A, GEN B, long v)
    GEN     ZX_disc(GEN x)
    int     ZX_is_squarefree(GEN x)
    GEN     ZX_resultant(GEN A, GEN B)
    GEN     ZX_QX_resultant(GEN A, GEN B)
    GEN     ZX_s_add(GEN y,long x)
    long    brent_kung_optpow(long d, long n)
    GEN     modulargcd(GEN a,GEN b)
    GEN     rescale_pol(GEN P, GEN h)
    GEN     unscale_pol(GEN P, GEN h)
    GEN     stopoly(ulong m, ulong p, long v)
    GEN     stopoly_gen(GEN m, GEN p, long v)
    int     u_pow(int p, int k)
    long    u_val(ulong n, ulong p)

    # RgX.c

    int     is_rational(GEN x)
    int     RgX_is_rational(GEN x)
    GEN     RgM_to_RgXV(GEN x, long v)
    GEN     RgM_to_RgXX(GEN x, long v,long w)
    GEN     RgM_zc_mul(GEN x, GEN y)
    GEN     RgM_zm_mul(GEN x, GEN y)
    GEN     RgV_to_RgX(GEN x, long v)
    GEN     RgV_zc_mul(GEN x, GEN y)
    GEN     RgV_zm_mul(GEN x, GEN y)
    GEN     RgX_divrem(GEN x,GEN y,GEN *r)
    GEN     RgX_mul(GEN x,GEN y)
    GEN     RgX_mulspec(GEN a, GEN b, long na, long nb)
    GEN     RgX_powers(GEN a, GEN T, long l)
    GEN     RgXQX_divrem(GEN x,GEN y,GEN T,GEN *r)
    GEN     RgXQX_mul(GEN x,GEN y,GEN T)
    GEN     RgXQX_red(GEN P, GEN T)
    GEN     RgXQX_RgXQ_mul(GEN x, GEN y, GEN T)
    GEN     RgX_Rg_mul(GEN y, GEN x)
    GEN     RgX_RgX_compo(GEN f, GEN x, GEN T)
    GEN     RgX_shift(GEN x, long n)
    GEN     RgX_sqr(GEN x)
    GEN     RgX_sqrspec(GEN a, long na)
    GEN     RgXV_to_RgM(GEN v, long n)
    GEN     RgX_to_RgV(GEN x, long N)
    GEN     RgXX_to_RgM(GEN v, long n)
    GEN     RgXY_swap(GEN x, long n, long w)

    # rootpol.c

    GEN     cleanroots(GEN x,long l)
    int     isrealappr(GEN x, long l)
    GEN     roots(GEN x,long l)
    GEN     roots0(GEN x,long flag,long l)

    #subcyclo.c

    GEN     galoissubcyclo(GEN N, GEN sg, long flag, long v)
    GEN     polsubcyclo(long n, long d, long v)
    GEN     subcyclo(long n, long d, long v)
    GEN     znstar_small(GEN zn)

    # subfields.c

    GEN     subfields(GEN nf,GEN d)
    GEN     subfields0(GEN nf,GEN d)

    # subgroup.c

    void    forsubgroup(entree *oep, GEN cyc, GEN bound, char *och)
    GEN     subgrouplist(GEN cyc, GEN bound)

    # stark.c

    GEN     bnrL1(GEN bnr, GEN sbgrp, long flag, long prec)
    GEN     bnrrootnumber(GEN bnr, GEN chi, long flag, long prec)
    GEN     bnrstark(GEN bnr, GEN subgroup, long prec)

    # sumiter.c

    GEN     direuler(void *E, GEN (*eval)(GEN,void*), GEN ga, GEN gb, GEN c)
    GEN     direuler0(entree *ep, GEN a, GEN b, char *ch, GEN c)
    GEN     divsum(GEN num,entree *ep, char *ch)
    void    fordiv(GEN a, entree *ep, char *ch)
    void    forpari(entree *ep, GEN a, GEN b, char *ch)
    void    forprime(entree *ep, GEN a, GEN b, char *ch)
    void    forstep(entree *ep, GEN a, GEN b, GEN s, char *ch)
    void    forvec(entree *ep, GEN x, char *ch, long flag)
    GEN     forvec_start(GEN x, long flag, GEN *d, GEN (**next)(GEN,GEN))
    GEN     intnum(void *E, GEN (*e)(GEN, void*), GEN a,GEN b, GEN tab, long prec)
    long    intnumstep(long prec)
    GEN     intnumromb0(entree *ep, GEN a, GEN b, char *ch, long flag, long prec)
    GEN     intnum0(entree *ep, GEN a, GEN b, char *ch, GEN tab, long prec)
    GEN     intcirc0(entree *ep, GEN a, GEN R, char *ch, GEN tab, long prec)
    GEN     intmellininv0(entree *ep, GEN sig, GEN x, char *ch, GEN tab, long prec)
    GEN     intmellininvshort(GEN sig, GEN x, GEN tab, long prec)
    GEN     intlaplaceinv0(entree *ep, GEN sig, GEN x, char *ch, GEN tab, long prec)
    GEN     intnuminit(GEN a, GEN b, long m, long prec)
    GEN     intnuminitgen0(entree *ep, GEN a, GEN b, char *ch, long m, long flag, long prec)
    GEN     intfuncinit0(entree *ep, GEN a, GEN b, char *ch, long flag, long m, long prec)
    # GEN     intnumdoub0(entree *epx, GEN a, GEN b, entree *epy, char *chc, char *chd, char *chf, GEN tabext, GEN tabint, long prec)
    GEN     intfoursin0(entree *ep, GEN a, GEN b, GEN x, char *ch, GEN tab, long prec)
    GEN     intfourcos0(entree *ep, GEN a, GEN b, GEN x, char *ch, GEN tab, long prec)
    GEN     intfourexp0(entree *ep, GEN a, GEN b, GEN x, char *ch, GEN tab, long prec)
    GEN     sumnum(void *E, GEN (*f)(GEN,void*), GEN a,GEN sig,GEN tab,long flag,long prec)
    GEN     sumnum0(entree *ep, GEN a, GEN sig, char *ch, GEN tab, long flag, long prec)
    GEN     sumnumalt(void *E, GEN (*f)(GEN,void*),GEN a,GEN s,GEN tab,long flag,long prec)
    GEN     sumnumalt0(entree *ep, GEN a, GEN sig, char *ch, GEN tab, long flag, long prec)
    GEN     sumnuminit(GEN sig, long m, long sgn, long prec)
    GEN     matrice(GEN nlig, GEN ncol,entree *ep1, entree *ep2, char *ch)
    GEN     polzag(long n, long m)
    GEN     polzagreel(long n, long m, long prec)
    GEN     prodeuler(void *E, GEN (*eval)(GEN,void*), GEN ga, GEN gb, long prec)
    GEN     prodeuler0(entree *ep, GEN a, GEN b, char *ch, long prec)
    GEN     prodinf(void *E, GEN (*eval)(GEN,void*), GEN a, long prec)
    GEN     prodinf0(entree *ep, GEN a, char *ch, long flag, long prec)
    GEN     prodinf1(void *E, GEN (*eval)(GEN,void*), GEN a, long prec)
    GEN     produit(entree *ep, GEN a, GEN b, char *ch, GEN x)
    GEN     somme(entree *ep, GEN a, GEN b, char *ch, GEN x)
    GEN     sumalt(void *E, GEN (*eval)(GEN,void*), GEN a, long prec)
    GEN     sumalt0(entree *ep, GEN a, char *ch,long flag, long prec)
    GEN     sumalt2(void *E, GEN (*eval)(GEN,void*), GEN a, long prec)
    GEN     sumpos(void *E, GEN (*eval)(GEN,void*), GEN a, long prec)
    GEN     sumpos2(void *E, GEN (*eval)(GEN,void*), GEN a, long prec)
    GEN     sumpos0(entree *ep, GEN a, char *ch, long flag,long prec)
    GEN     suminf(void *E, GEN (*eval)(GEN,void*), GEN a, long prec)
    GEN     suminf0(entree *ep, GEN a, char *ch, long prec)
    GEN     vecteur(GEN nmax, entree *ep, char *ch)
    GEN     vecteursmall(GEN nmax, entree *ep, char *ch)
    GEN     vvecteur(GEN nmax, entree *ep, char *ch)
    GEN     zbrent0(entree *ep, GEN a, GEN b, char *ch, long prec)
    GEN     zbrent(void *E, GEN (*eval)(GEN,void*), GEN a, GEN b, long prec)

    # thue.c

    GEN     bnfisintnorm(GEN x, GEN y)
    GEN     thue(GEN thueres, GEN rhs, GEN ne)
    GEN     thueinit(GEN pol, long flag, long prec)

    # trans1.c

    GEN     Pi2n(long n, long prec)
    GEN     PiI2(long prec)
    GEN     PiI2n(long n, long prec)
    void    consteuler(long prec)
    void    constpi(long prec)
    GEN     exp_Ir(GEN x)
    GEN     gcos(GEN x, long prec)
    GEN     gcotan(GEN x, long prec)
    GEN     gexp(GEN x, long prec)
    GEN     glog(GEN x, long prec)
    GEN     gpow(GEN x, GEN n, long prec)
    GEN     gpowgs(GEN x, long n)
    GEN     gsin(GEN x, long prec)
    void    gsincos(GEN x, GEN *s, GEN *c, long prec)
    GEN     gsqrt(GEN x, long prec)
    GEN     gsqrtn(GEN x, GEN n, GEN *zetan, long prec)
    GEN     gtan(GEN x, long prec)
    GEN     log0(GEN x,long flag, long prec)
    GEN     mpcos(GEN x)
    GEN     mpeuler(long prec)
    GEN     mpexp(GEN x)
    GEN     mpexp1(GEN x)
    GEN     mplog(GEN x)
    GEN     mplog2(long prec)
    GEN     mppi(long prec)
    GEN     mpsin(GEN x)
    void    mpsincos(GEN x, GEN *s, GEN *c)
    GEN     sqrtr(GEN x)
    GEN     sqrtnr(GEN x, long n)
    GEN     palog(GEN x)
    GEN     powgi(GEN x, GEN n)
    GEN     teich(GEN x)

    # trans2.c

    GEN     bernfrac(long n)
    GEN     bernreal(long n, long prec)
    GEN     bernvec(long nomb)
    GEN     gach(GEN x, long prec)
    GEN     gacos(GEN x, long prec)
    GEN     garg(GEN x, long prec)
    GEN     gash(GEN x, long prec)
    GEN     gasin(GEN x, long prec)
    GEN     gatan(GEN x, long prec)
    GEN     gath(GEN x, long prec)
    GEN     gch(GEN x, long prec)
    GEN     ggamd(GEN x, long prec)
    GEN     ggamma(GEN x, long prec)
    GEN     glngamma(GEN x, long prec)
    GEN     gpsi(GEN x, long prec)
    GEN     gsh(GEN x, long prec)
    GEN     gth(GEN x, long prec)
    void    mpbern(long nomb, long prec)

    # trans3.c

    GEN     agm(GEN x, GEN y, long prec)
    GEN     dilog(GEN x, long prec)
    GEN     eint1(GEN x, long prec)
    GEN     eta(GEN x, long prec)
    GEN     eta0(GEN x, long flag,long prec)
    GEN     gerfc(GEN x, long prec)
    GEN     gpolylog(long m, GEN x, long prec)
    void    gpolylogz(long m, GEN x, GEN y)
    GEN     gzeta(GEN x, long prec)
    GEN     hyperu(GEN a, GEN b, GEN gx, long prec)
    GEN     incgam(GEN a, GEN x, long prec)
    GEN     incgam0(GEN a, GEN x, GEN z,long prec)
    GEN     incgam1(GEN a, GEN x, long prec)
    GEN     incgam2(GEN a, GEN x, long prec)
    GEN     incgam3(GEN a, GEN x, long prec)
    GEN     incgamc(GEN a, GEN x, long prec)
    GEN     hbessel1(GEN n, GEN z, long prec)
    GEN     hbessel2(GEN n, GEN z, long prec)
    GEN     ibessel(GEN n, GEN z, long prec)
    GEN     jbessel(GEN n, GEN z, long prec)
    GEN     jbesselh(GEN n, GEN z, long prec)
    GEN     nbessel(GEN n, GEN z, long prec)
    GEN     jell(GEN x, long prec)
    GEN     kbessel(GEN nu, GEN gx, long prec)
    GEN     kbessel0(GEN nu, GEN gx, long flag,long prec)
    GEN     kbessel2(GEN nu, GEN x, long prec)
    GEN     polylog(long m, GEN x, long prec)
    GEN     polylog0(long m, GEN x, long flag, long prec)
    GEN     polylogd(long m, GEN x, long prec)
    GEN     polylogdold(long m, GEN x, long prec)
    GEN     polylogp(long m, GEN x, long prec)
    GEN     szeta(long x, long prec)
    GEN     theta(GEN q, GEN z, long prec)
    GEN     thetanullk(GEN q, long k, long prec)
    GEN     trueeta(GEN x, long prec)
    GEN     veceint1(GEN C, GEN nmax, long prec)
    GEN     vecthetanullk(GEN q, long k, long prec)
    GEN     weber0(GEN x, long flag,long prec)
    GEN     wf(GEN x, long prec)
    GEN     wf2(GEN x, long prec)

    GEN     padicfieldslist(GEN p, GEN m, GEN d, long flag)

cdef extern from 'pari/pari.h':
    extern GEN geuler
    extern GEN gpi


cdef extern from 'pari/paripriv.h':
    struct __x:
        char format  # e,f,g
        long fieldw  # 0 (ignored) or field width
        long sigd    # -1 (all) or number of significant digits printed */
        int sp       # 0 = suppress whitespace from output */
        int prettyp  # output style: raw, prettyprint, etc */
        int TeXstyle
    ctypedef __x pariout_t

    struct __z:
        jmp_buf env
        pariout_t *fmt
    ctypedef __z gp_data
    extern gp_data* GP_DATA

# Initialize the relevant part of GP for exception trapping.
cdef gp_data __GP_DATA
GP_DATA = &__GP_DATA

cdef int REAL_PREC
REAL_PREC = 5

cdef int is_int_or_real_type(GEN g):
    cdef int t
    t = typ(g)
    return t==t_INT or t==t_REAL

cdef size_t fix_size(size_t a):
    cdef size_t b
    b = a - (a & (BYTES_IN_LONG-1))     # sizeof(long) | b <= a
    if b < 1024:
        b = 1024
    return b

cdef int init_stack(size_t size) except -1:
    cdef size_t s

    global top, bot, avma, stack_avma

    # delete this if get core dumps and change the 2* to a 1* below.
    if bot:
        #print "Freeing the old stack."
        PyMem_Free(<void*>bot)

    if size == 0:
        size = 2*(top-bot)

    # if size == -1, then allocate the biggest chunk possible
    if size == -1:
        s = 4294967295
        while True:
            s = fix_size(s)
            bot = <pari_sp> PyMem_Malloc(s)
            if bot:
                break
            s = s/2
    else:
        # Decide on size
        s = fix_size(size)

        # Alocate memory for new stack using Python's memory allocator.
        # As explained in the python/C api reference, using this instead
        # of malloc is much better (and more platform independent, etc.)
        bot = <pari_sp> PyMem_Malloc(s)
        if not bot:
            raise MemoryError, "Unable to allocate %s bytes memory for PARI."%(<long>size)
    #endif

    top = bot + s
    avma = top
    stack_avma = avma


def _my_sigint(signum, frame):
    raise KeyboardInterrupt

#def _my_sigsegv(signum, frame):
#    exit(0)
    #raise RuntimeError

def _my_sigpipe(signum, frame):
    #raise IOError, (32,'Broken pipe')
    # If I do anything, it messes up Ipython's pager.
    pass

cdef pariout_t fmt
cdef unsigned long num_primes

def _init_fmt():
    GP_DATA.fmt.prettyp = 0

def pari_version():
    return str(PARIVERSION)

def _init(long size=200000000, unsigned long maxprime=500000):
#def _init(long size=-1, unsigned long maxprime=500000):
    """
    Initialize the PARI system.

    INPUT:
        size -- long, the number of bytes for the initial PARI stack
                (see note below)
        maxprime -- unsigned long, upper limit on a precomputed prime
                    number table  (default: 500000)

    NOTES:

        * In py_pari, the PARI stack is different than in gp or the
          PARI C library.  In Python, instead of the PARI stack
          holding the results of all computations, it *only* holds the
          results of an individual computation.  Each time a new
          Python/PARI object is computed, it it copied to its own
          space in the Python heap, and the memory it occupied on the
          PARI stack is freed.  Thus it is not necessary to make the
          stack very large.  Also, unlike in PARI, if the stack does
          overflow, in most cases the PARI stack is automatically
          increased and the relevant step of the computation rerun.

          This design obviously involves some performance penalties
          over the way PARI works, but it scales much better and is
          far more robus for large projects.

        * If you do not want prime numbers, put maxprime=2, but be
          careful because many PARI functions require this table.  If
          you get the error message "not enough precomputed primes",
          increase this parameter.

    """
    global initialized, num_primes, ZERO, ONE, TWO, avma, top, bot
    if bot:
        raise RuntimeError, "pari_init has already been called."
    #print "Initializing PARI (size=%s, maxprime=%s)"%(size,maxprime)
    pari_init(1024, maxprime)

    init_stack(size)
    _init_fmt()

    # Take control of SIGINT back from PARI.
    import signal
    signal.signal(signal.SIGINT, _my_sigint)
    signal.signal(signal.SIGSEGV, signal.SIG_DFL)

    # We do the following, since otherwise the IPython pager
    # causes sage to crash when it is exited early.    Again,
    # this is because when PARI was initialized it set a trap
    # for this signal.
    signal.signal(signal.SIGPIPE, _my_sigpipe)
    initialized = 1
    stack_avma = avma
    num_primes = maxprime
    ZERO = pari(0)    # todo: gen_0
    ONE = pari(1)
    TWO = pari(2)
    return 0


############################################################
# conversion between GEN and string types
# Note that string --> GEN evaluates the string in PARI,
# where GEN_to_str returns a Python string.
############################################################


cdef GEN str_to_GEN(char* s) except NULL:
    cdef int n
    n = setjmp(GP_DATA.env)
    if n: _error(n,s)
    return flisseq(s)

cdef object GEN_to_str(GEN g):
    cdef char* c
    cdef int n
    _sig_on
    c = GENtostr(g)
    _sig_off
    s = str(c)
    free(c)
    return s




############################################################
# Extension class that models the PARI GEN type.
############################################################

cdef class gen:
    """
    Python extension class that models the PARI GEN type.
    """
    cdef GEN g
    cdef object refers_to
    cdef pari_sp b

    cdef init(self, GEN g, pari_sp b):
        """
            g -- PARI GEN
            b -- pointer to memory chunk where PARI gen lives
                 (if nonzero then this memory is freed when the object
                  goes out of scope)
        """
        self.g = g
        self.b = b
        self.refers_to = {}

    def __dealloc__(self):
        if self.b:
            #print "Freeing a PARI/Python object."
            PyMem_Free(<void*> self.b)

    def __repr__(self):
        return(GEN_to_str(self.g))

    cdef GEN _gen(self):
        return self.g

    def __reduce__(self):
        z = str(self)
        import sage.pari.py_pari_py
        return sage.pari.py_pari_py.pari, (z,)

    ###########################################
    # ARITHMETIC
    ###########################################

    def _add(gen self, gen other):
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error adding PARI objects")
        _sig_on
        return new_gen(gadd(forcecopy(self.g), other.g))

    def _sub(gen self, gen other):
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error subtracting PARI objects")
        _sig_on
        return new_gen(gsub(forcecopy(self.g), other.g))

    def _mul(gen self, gen other):
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error multiplying PARI objects")
        _sig_on
        return new_gen(gmul(forcecopy(self.g), other.g))

    def _div(gen self, gen other):
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error dividing PARI objects")
        _sig_on
        return new_gen(gdiv(forcecopy(self.g), forcecopy(other.g)))

    def _mod(gen self, gen other):
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error quotienting PARI objects")
        _sig_on
        return new_gen(gmod(forcecopy(self.g), other.g))

    def __add__(self, other):
        self = pari(self); other = pari(other)
        return self._add(other)

    def __sub__(self, other):
        self = pari(self); other = pari(other)
        return self._sub(other)

    def __mul__(self, other):
        self = pari(self); other = pari(other)
        return self._mul(other)

    def __div__(self, other):
        self = pari(self); other = pari(other)
        return self._div(other)

    def __mod__(self, other):
        self = pari(self); other = pari(other)
        return self._mod(other)

    def __neg__(gen self):
        _sig_on
        return new_gen(gneg(forcecopy(self.g)))

    def __pow__(self, n, m):
        cdef gen _n, _self
        _n = pari(n)
        _self = pari(self)

        cdef int j
        j = setjmp(GP_DATA.env)
        if j:
            _error(j,"error exponentiating PARI object")

        _sig_on
        return new_gen(gpow(forcecopy(_self.g), _n.g, REAL_PREC))

    def __xor__(gen self, n):
        raise RuntimeError, "Use ** for exponentiation, not '^', which means xor\n"+\
              "in Python, and has the wrong precedence."

    def __rshift__(gen self, long n):
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error right shifting PARI object")
        _sig_on
        return new_gen(gshift(forcecopy(self.g), -n))

    def __lshift__(gen self, long n):
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error left shifting PARI object")
        _sig_on
        return new_gen(gshift(forcecopy(self.g), n))

    def __invert__(gen self):
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error inverting PARI object")
        _sig_on
        return new_gen(ginv(forcecopy(self.g)))

    ###########################################
    # ACCESS
    ###########################################
    def __getitem__(gen self, n):
        """
        Return a *copy* of the nth entry.

        The indexing is 0-based, like in Python.  However, *copying*
        the nth entry instead of returning reference is different
        than what one usually does in Python.  However, we do it
        this way for consistency with the PARI/GP interpreter, which
        does make copies.
        """
        if typ(self.g) == t_POL:
            return self.polcoeff(n)

        if isinstance(n, tuple):
            i, j = n[0], n[1]
            if i < 0 or i >= self.nrows():
                raise IndexError, "row i(=%s) must be between 0 and %s"%(i,self.nrows())
            if j < 0 or j >= self.ncols():
                raise IndexError, "column j(=%s) must be between 0 and %s"%(j,self.ncols())
            return new_ref( <GEN> (<GEN>(self.g)[j+1])[i+1], self )

        if n < 0 or n >= glength(self.g):
            raise IndexError, "index (%s) must be between 0 and %s"%(n,glength(self.g)-1)
        return new_ref(<GEN> (self.g[n+1]), self)


    def __getslice__(self,  Py_ssize_t i,  Py_ssize_t j):
        """
        EXAMPLES:
            sage: v = pari(xrange(20))
            sage: v[2:5]
            [2, 3, 4]
            sage: v[:]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[10:]
            [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[:5]
            [0, 1, 2, 3, 4]
            sage: v
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: v[-1]
            Traceback (most recent call last):
            ...
            IndexError: index (-1) must be between 0 and 19
            sage: v[:-3]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
            sage: v[5:]
            [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        """
        cdef long l, k
        l = glength(self.g)
        if j >= l: j = l
        if i < 0: i = 0
        if j <= i: return vector(0)
        v = vector(j-i)
        for k from i <= k < j:
            v[k-i] = self[k]
        return v

    def __setitem__(gen self, n, y):
        r"""
        Set the nth entry to a reference to y.

        \begin{notice}
        \begin{itemize}
            \item There is a known BUG: If v is a vector and entry i of v is a vector,
                       \code{v[i][j] = x}
               should set entry j of v[i] to x.  Instead it sets it to nonsense.
               I do not understand why this occurs.  The following is a safe way
               to do the same thing:
                        \code{tmp = v[i]; tmp[j] = x    }

            \item The indexing is 0-based, like everywhere else in Python, but
               {\em unlike} in GP/PARI.

            \item Assignment sets the nth entry to a reference to y,
               assuming y is an object of type gen.  This is the same
               as in Python, but {\em different} than what happens in the
               gp interpreter, where assignment makes a copy of y.

            \item Because setting creates references it is {\em possible} to
               make circular references, unlike in GP.  Do {\em not} do
               this (see the example below).  If you need circular
               references, work at the Python level (where they work
               well), not the PARI object level.
        \end{itemize}
        \end{notice}

        EXAMPLES:
            sage: v = pari(range(10))
            sage: v
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[0] = 10
            sage: w = pari([5,8,-20])
            sage: v
            [10, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: v[1] = w
            sage: v
            [10, [5, 8, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: w[0] = -30
            sage: v
            [10, [-30, 8, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: t = v[1]; t[1] = 10   # don't do v[1][1] !!! (because of mysterious BUG...)
            sage: v
            [10, [-30, 10, -20], 2, 3, 4, 5, 6, 7, 8, 9]
            sage: w
            [-30, 10, -20]

        Finally, we create a circular reference:
            sage: v = pari([0])
            sage: w = pari([v])
            sage: v
            [0]
            sage: w
            [[0]]
            sage: v[0] = w
            sage: # Now there is a circular reference.  Trying to access v[0] will crash SAGE.

        """
        cdef int e
        e = setjmp(GP_DATA.env)
        if e: _error(e,"__setitem__")

        cdef gen x
        x = pari(y)
        if isinstance(n, tuple):
            self.refers_to[n] = x
            i = n[0]
            j = n[1]
            if i < 0 or i >= self.nrows():
                raise IndexError, "row i(=%s) must be between 0 and %s"%(i,self.nrows())
            if j < 0 or j >= self.ncols():
                raise IndexError, "column j(=%s) must be between 0 and %s"%(j,self.ncols())
            (<GEN>(self.g)[j+1])[i+1] = <long> x.g
            return

        if n < 0 or n >= glength(self.g):
            raise IndexError, "index (%s) must be between 0 and %s"%(n,glength(self.g)-1)
        self.refers_to[n] = x       # so python memory manager will work correctly
                                    # and not free x if PARI part of self is the
                                    # only thing pointing to it.
        e = setjmp(GP_DATA.env)
        if e: _error(e,"__setitem__")
        (self.g)[n+1] = <long>(x.g)

    def __len__(gen self):
        return glength(self.g)

    def ncols(gen self):
        return glength(self.g)

    def nrows(gen self):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"nrows")
        return glength(<GEN>(self.g[1]))

    def polcoeff(gen self, long n, var=-1):
        """
        EXAMPLES:
            sage: f = pari("x^2 + y^3 + x*y")
            sage: f
            x^2 + y*x + y^3
            sage: f.polcoeff(1)
            y
            sage: f.polcoeff(3)
            0
            sage: f.polcoeff(3, "y")
            1
            sage: f.polcoeff(1, "y")
            x
        """
        return new_gen(polcoeff0(forcecopy(self.g), n, get_var(var)))


    ###########################################
    # ?
    ###########################################

    def __cmp__(gen self, gen other):
        """
        Comparisons.
        """
        if gegal(self.g, other.g):
            return 0
        cdef int tg, to
        tg = typ(self.g)
        to = typ(other.g)

        if is_int_or_real_type(self.g) and is_int_or_real_type(other.g):
            return gcmp(self.g, other.g)

        sg = str(self)
        so = str(other)
        if sg < so:
            return -1
        else:
            return 1


    def __richcmp__(gen self, other, int op):
        """
        EXAMPLES:
            sage: a = pari(5)
            sage: b = 10
            sage: a < b
            True
            sage: a <= b
            True
            sage: a <= 5
            True
            sage: a > b
            False
            sage: a >= b
            False
            sage: a >= pari(10)
            False
            sage: a == 5
            True
            sage: a is 5
            False
        """
        cdef int n
        cdef gen x
        x = pari(other)
        n = self.__cmp__(x)
        if op == 0:
            return bool(n < 0)
        elif op == 1:
            return bool(n <= 0)
        elif op == 2:
            return bool(n == 0)
        elif op == 3:
            return bool(n != 0)
        elif op == 4:
            return bool(n > 0)
        elif op == 5:
            return bool(n >= 0)

    def copy(gen self):
        return new_gen(forcecopy(self.g))


    ###########################################
    # Conversion --> Python
    # Try to convert to a meaningful python object
    # in various ways
    ###########################################

    def list_str(gen self):
        """
        Return str that might correctly evaluate to a Python-list.
        """
        s = str(self)
        if s[:4] == "Mat(":
            s = "[" + s[4:-1] + "]"
        s = s.replace("~","")
        if s.find(";") != -1:
            s = s.replace(";","], [")
            s = "[" + s + "]"
            return eval(s)
        else:
            return eval(s)

    def __int__(gen self):
        """
        Return Python int.  Very fast unless the number is too large too
        fit into a C int, in which case a string is used.
        """
        return int(str(self))
        # NOTE: Could use int_unsafe below, which would be much faster, but
        # the default PARI prints annoying stuff to the screen when
        # the number is large.

    def int_unsafe(gen self):
        """
        Returns int form of self.

        This is about 5 times faster than the usual int conversion.

        This function is fast, but if self does not fit in a C int, it
        will print an error and continue running.  This error message
        is printed by the PARI C library, so I can't turn it off.
        """
        cdef int err
        err = setjmp(GP_DATA.env)
        if err:
            print "(py_pari.int_unsafe(): The above PARI error message can be ignored.)"
            return int(str(self))
        return gtolong(self.g)


    def python_list(gen self):
        """
        Return a Python list of the PARI gens.  This object must be of
        type t_VEC

        INPUT: None
        OUTPUT:
           list -- Python list whose elements are the
                   elements of the input gen.
        EXAMPLES:
            sage: v=pari([1,2,3,10,102,10])
            sage: w = v.python_list()
            sage: w
            [1, 2, 3, 10, 102, 10]
            sage: type(w[0])
            <type '_py_pari.gen'>
        """
        if typ(self.g) != t_VEC:
            raise TypeError, "Object (=%s) must be of type t_VEC."%self
        cdef long n, m
        cdef gen t
        m = glength(self.g)
        V = []
        for n from 0 <= n < m:
            t = new_ref(<GEN> (self.g[n+1]), V)
            V.append(t)
        return V

    def python(gen self):
        """
        Return Python eval of self.
        """
        import sage.pari.py_pari_py
        return sage.pari.py_pari_py.python(self)

    def __long__(gen self):
        """
        Return Python long.
        """
        return long(str(self))
        #cdef int err
        #err = setjmp(GP_DATA.env)
        #if err:
        #    if err == 28: return long(str(self))
        #    _error(err,"error in __long__")
        #return long(gtolong(self.g))

    def __float__(gen self):
        """
        Return Python float.
        """
        #return float(str(self))
        cdef int err
        err = setjmp(GP_DATA.env)
        if err:
            if err == 28: return float(str(self))
            _error(err,"error in __float__")
        return gtodouble(self.g)

    def __bool__(gen self):
        cdef int err
        err = setjmp(GP_DATA.env)
        if err:
            _error(err,"error in __bool__")
        return bool(self.g != stoi(0))


    ###########################################
    # arith1.c
    ###########################################
    def isprime(gen self, flag=0):
        """
        isprime(x, flag=0): Returns True if x is a PROVEN prime
        number, and False otherwise.

        INPUT:
            flag -- int
                    0 (default): use a combination of algorithms.
                    1: certify primality using the Pocklington-Lehmer Test.
                    2: certify primality using the APRCL test.
        OUTPUT:
            bool -- True or False
        EXAMPLES:
            sage: pari(9).isprime()
            False
            sage: pari(17).isprime()
            True
            sage: n = pari(561)    # smallest Carmichael number
            sage: n.isprime()      # not just a pseudo-primality test!
            False
            sage: n.isprime(1)
            False
            sage: n.isprime(2)
            False
        """
        _sig_on
        t = bool(gisprime(self.g, flag) != stoi(0))
        _sig_off
        return t


    ###########################################
    # 1: Standard monadic or dyadic OPERATORS
    ###########################################
    def divrem(gen x, y, var=-1):
        """

        divrem(x, y, {v}): Euclidean division of x by y giving as a
            2-dimensional column vector the quotient and the
            remainder, with respect to v (to main variable if v is
            omitted).

        sage:
        """
        cdef gen _y
        _y = pari(y)
        _sig_on
        return new_gen(divrem(forcecopy(x.g), _y.g, get_var(var)))

    def lex(gen x, y):
        """

        lex(x,y): Compare x and y lexicographically (1 if x>y, 0 if
            x==y, -1 if x<y)

        """
        cdef gen _y
        _y = pari(y)
        return lexcmp(x.g, _y.g)

    def max(gen x, y):
        """
        max(x,y): Return the maximum of x and y.
        """
        cdef gen _y
        _y = pari(y)
        return new_gen(gmax(forcecopy(x.g), _y.g))

    def min(gen x, y):
        """
        min(x,y): Return the minumum of x and y.
        """
        cdef gen _y
        _y = pari(y)
        return new_gen(gmin(forcecopy(x.g), _y.g))

    def shift(gen x, long n):
        """
        shift(x,n): shift x left n bits if n>=0, right -n bits if n<0.
        """
        return new_gen(gshift(forcecopy(x.g), n))

    def shiftmul(gen x, long n):
        """
        shiftmul(x,n): Return the product of x by $2^n$.
        """
        return new_gen(gmul2n(forcecopy(x.g), n))

    def sign(gen x):
        """
        sign(x): Return the sign of x, where x is of type integer, real
            or fraction.
        """
        return gsigne(x.g)

    def vecmax(gen x):
        """
        vecmax(x): Return the maximum of the elements of the vector/matrix x,
        """
        _sig_on
        return new_gen(vecmax(forcecopy(x.g)))


    def vecmin(gen x):
        """
        vecmin(x): Return the maximum of the elements of the vector/matrix x,
        """
        _sig_on
        return new_gen(vecmin(forcecopy(x.g)))



    ###########################################
    # 2: CONVERSIONS and similar elementary functions
    ###########################################


    def Col(gen x):
        """
        Col(x): Transforms the object x into a column vector.

        The vector will have only one component, except in the
        following cases:

           * When x is a vector or a quadratic form, the resulting
             vector is the initial object considered as a column
             vector.

           * When x is a matrix, the resulting vector is the column of
             row vectors comprising the matrix.

           * When x is a character string, the result is a column of
             individual characters.

           * When x is a polynomial, the coefficients of the vector
             start with the leading coefficient of the polynomial.

           * When x is a power series, only the significant
             coefficients are taken into account, but this time by
             increasing order of degree.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari(1.5).Col()
            [1.500000000000000000000000000]~
            sage: pari([1,2,3,4]).Col()
            [1, 2, 3, 4]~
            sage: pari('[1,2; 3,4]').Col()
            [[1, 2], [3, 4]]~
            sage: pari('"SAGE"').Col()
            ["S", "A", "G", "E"]~
            sage: pari('3*x^3 + x').Col()
            [3, 0, 1, 0]~
            sage: pari('x + 3*x^3 + O(x^5)').Col()
            [1, 0, 3, 0]~
        """
        return new_gen(gtocol(x.g))

    def List(gen x):
        """
        List(x): transforms the PARI vector or list x into a list.

        EXAMPLES:
            sage: v = pari([1,2,3])
            sage: v
            [1, 2, 3]
            sage: v.type()
            't_VEC'
            sage: w = v.List()
            sage: w
            List([1, 2, 3])
            sage: w.type()
            't_LIST'
        """
        return new_gen(gtolist(forcecopy(x.g)))

    def Mat(gen x):
        """
        Mat(x): Returns the matrix defined by x.

           * If x is already a matrix, a copy of x is created and
             returned.

           * If x is not a vector or a matrix, this function returns a
             1x1 matrix.

           * If x is a row (resp. column) vector, this functions
             returns a 1-row (resp. 1-column) matrix, *unless* all
             elements are column (resp. row) vectors of the same
             length, in which case the vectors are concatenated
             sideways and the associated big matrix is returned.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- a PARI matrix

        EXAMPLES:
            sage: x = pari(5)
            sage: x.type()
            't_INT'
            sage: y = x.Mat()
            sage: y
            Mat(5)
            sage: y.type()
            't_MAT'
            sage: x = pari('[1,2;3,4]')
            sage: x.type()
            't_MAT'
            sage: x = pari('[1,2,3,4]')
            sage: x.type()
            't_VEC'
            sage: y = x.Mat()
            sage: y
            Mat([1, 2, 3, 4])
            sage: y.type()
            't_MAT'

            sage: v = pari('[1,2;3,4]').Vec(); v
            [[1, 3]~, [2, 4]~]
            sage: v.Mat()
            [1, 2; 3, 4]
            sage: v = pari('[1,2;3,4]').Col(); v
            [[1, 2], [3, 4]]~
            sage: v.Mat()
            [1, 2; 3, 4]
        """
        return new_gen(gtomat(forcecopy(x.g)))

    def Mod(gen x, y):
        """
        Mod(x, y): Returns the object x modulo y, denoted Mod(x, y).

        The input y must be a an integer or a polynomial:

           * If y is an INTEGER, x must also be an integer, a rational
             number, or a p-adic number compatible with the modulus y.

           * If y is a POLYNOMIAL, x must be a scalar (which is not a polmod),
             a polynomial, a rational function, or a power series.

        WARNING: This function is not the same as x % y, the result of
        which is an integer or a polynomial.

        INPUT:
            x -- gen
            y -- integer or polynomial

        OUTPUT:
            gen -- intmod or polmod

        EXAMPLES:
            sage: z = pari(3)
            sage: x = z.Mod(pari(7))
            sage: x
            Mod(3, 7)
            sage: x**2
            Mod(2, 7)
            sage: x**100
            Mod(4, 7)
            sage: x.type()
            't_INTMOD'

            sage: f = pari("x^2 + x + 1")
            sage: g = pari("x")
            sage: a = g.Mod(f)
            sage: a
            Mod(x, x^2 + x + 1)
            sage: a*a
            Mod(-x - 1, x^2 + x + 1)
            sage: a.type()
            't_POLMOD'
        """
        cdef gen _y
        _y = pari(y)
        return new_gen(gmodulcp(forcecopy(x.g),_y.g))

    def Pol(gen x, v=-1):
        """
        Pol(x, {v}): convert x into a polynomial with main variable v
        and return the result.

           * If x is a scalar, returns a constant polynomial.

           * If x is a power series, the effect is identical
              to \kbd{truncate}, i.e.~it chops off the $O(X^k)$.

           * If x is a vector, this function creates the polynomial
             whose coefficients are given in x, with x[0]
             being the leading coefficient (which can be zero).

        WARNING: This is *not* a substitution function. It will not
        transform an object containing variables of higher priority
        than v:
            sage.: pari('x+y').Pol('y')
            Traceback (most recent call last):
            ...
            RuntimeError: PARI error 8: error in Pol

        INPUT:
            x -- gen
            v -- (optional) which variable, defaults to 'x'
        OUTPUT:
            gen -- a polynomial
        EXAMPLES:
            sage: v = pari("[1,2,3,4]")
            sage: f = v.Pol()
            sage: f
            x^3 + 2*x^2 + 3*x + 4
            sage: f*f
            x^6 + 4*x^5 + 10*x^4 + 20*x^3 + 25*x^2 + 24*x + 16

            sage: v = pari("[1,2;3,4]")
            sage: v.Pol()
            [1, 3]~*x + [2, 4]~
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in Pol")
        return new_gen(gtopoly(forcecopy(x.g), get_var(v)))

    def Polrev(gen x, v=-1):
        """
        Polrev(x, {v}): Convert x into a polynomial with main variable
        v and return the result.  This is the reverse of Pol if x is a
        vector, otherwise it is identical to Pol.   By "reverse" we mean
        that the coefficients are reveresed.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- a polynomial
        EXAMPLES:
            sage: v = pari("[1,2,3,4]")
            sage: f = v.Polrev()
            sage: f
            4*x^3 + 3*x^2 + 2*x + 1
            sage: v.Pol()
            x^3 + 2*x^2 + 3*x + 4
            sage: v.Polrev('y')
            4*y^3 + 3*y^2 + 2*y + 1

        Note that Polrev does *not* reverse the coefficients of a polynomial!
            sage: f
            4*x^3 + 3*x^2 + 2*x + 1
            sage: f.Polrev()
            4*x^3 + 3*x^2 + 2*x + 1
            sage: v = pari("[1,2;3,4]")
            sage: v.Polrev()
            [2, 4]~*x + [1, 3]~
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in Polrev")
        return new_gen(gtopolyrev(forcecopy(x.g),get_var(v)))


    def Qfb(gen a, b, c, D=0):
        """
        Qfb(a,b,c,{D=0.}): Returns the binary quadratic form
        $$
                   ax^2 + bxy + cy^2.
        $$
        The optional D is 0 by default and initializes Shanks's
        distance if $b^2 - 4ac > 0$.

        NOTE: Negative definite forms are not implemented, so use their
        positive definitine counterparts instead.  (I.e., if f is a
        negative definite quadratic form, then -f is positive
        definite.)

        INPUT:
            a -- gen
            b -- gen
            c -- gen
            D -- gen (optional, defaults to 0)
        OUTPUT:
            gen -- binary quadratic form
        EXAMPLES:
            sage: pari(3).Qfb(7, 2)
            Qfb(3, 7, 2, 0.E-28)
        """
        cdef gen _b, _c, _D
        _b = pari(b)
        _c = pari(c)
        _D = pari(D)
        _sig_on
        return new_gen(Qfb0(forcecopy(a.g), _b.g, _c.g, _D.g, REAL_PREC))


    def Ser(gen x, v=-1):
        """
        Ser(x,{v=x}): Create a power series from x with main variable v
        and return the result.

           * If x is a scalar, this gives a constant power series with
             precision given by the default series precision, as returned
             by get_series_precision().

           * If x is a polynomial, the precision is the greatest of
             get_series_precision() and the degree of the polynomial.

           * If x is a vector, the precision is similarly given, and
             the coefficients of the vector are understood to be the
             coefficients of the power series starting from the
             constant term (i.e.~the reverse of the function Pol).

        WARNING: This is *not* a substitution function. It will not
        transform an object containing variables of higher priority than v.

        INPUT:
            x -- gen
            v -- PARI variable (default: x)
        OUTPUT:
            gen -- PARI object of PARI type t_SER
        EXAMPLES:
            sage: pari(2).Ser()
            2 + O(x^16)
            sage: x = pari([1,2,3,4,5])
            sage: x.Ser()
            1 + 2*x + 3*x^2 + 4*x^3 + 5*x^4 + O(x^5)
            sage: f = x.Ser('v'); print f
            1 + 2*v + 3*v^2 + 4*v^3 + 5*v^4 + O(v^5)
            sage: pari(1)/f
            1 - 2*v + v^2 + O(v^5)
            sage: pari(1).Ser()
            1 + O(x^16)
        """
        return new_gen(gtoser(x.g, get_var(v)))


    def Set(gen x):
        """
        Set(x): convert x into a set, i.e. a row vector of strings in
        increasing lexicographic order.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- a vector of strings in increasing lexicographic order.
        EXAMPLES:
            sage: pari([1,5,2]).Set()
            ["1", "2", "5"]
            sage: pari([]).Set()     # the empty set
            []
            sage: pari([1,1,-1,-1,3,3]).Set()
            ["-1", "1", "3"]
            sage: pari(1).Set()
            ["1"]
            sage: pari('1/(x*y)').Set()
            ["1/(y*x)"]
            sage: pari('["bc","ab","bc"]').Set()
            ["ab", "bc"]
        """
        _sig_on
        return new_gen(gtoset(x.g))


    def Str(gen x):
        """
        Str(x): Concatenate the entries of the vector x into a single
        string.  If x is not a t_VEC its print representation is
        returned.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- a PARI gen of type t_STR, i.e., a PARI string
        EXAMPLES:
            sage: pari([1,2,['abc',1]]).Str()
            12[abc, 1]
            sage: pari([1,1,1.54]).Str()
            111.540000000000000000000000000
            sage: pari([1,1.54,1]).Str()
            11.5400000000000000000000000001
            sage: pari(1).Str()       # 1 is automatically converted to string rep
            1
            sage: x = pari('x')       # PARI variable "x"
            sage: x.Str()             # is converted to string rep.
            x
            sage: x.Str().type()
            't_STR'
        """
        cdef char* c
        cdef int n
        if x.type() != 't_VEC':
            n = setjmp(GP_DATA.env)
            if n:
                _error(n,"converting a gen to a string")
            c = GENtostr(x.g)
            v = new_gen(strtoGENstr(c))
            free(c)
            return v
        _sig_on
        return new_gen(Str(x.g))

    def Strchr(gen x):
        """
        Strchr(x): converts x to a string, translating each integer
        into a character (in ASCII).

        NOTE: Vecsmall is (essentially) the inverse to Strchr().

        INPUT:
            x -- PARI vector of integers
        OUTPUT:
            gen -- a PARI string
        EXAMPLES:
            sage: pari([65,66,123]).Strchr()
            AB{
            sage: pari('"SAGE"').Vecsmall()   # pari('"SAGE"') --> PARI t_STR
            Vecsmall([83, 65, 71, 69])
            sage: _.Strchr()
            SAGE
            sage: pari([83, 65, 71, 69]).Strchr()
            SAGE
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            raise ValueError, "Cannot convert %s, interpreted as ASCII values, to a string."%x
        _sig_on
        return new_gen(Strchr(x.g))

    def Strexpand(gen x):
        """
        Strexpand(x): Concatenate the entries of the vector x into a
        single string, performing tilde expansion.

        NOTE: I have no clue what the point of this function is. -- William
        """
        if x.type() != 't_VEC':
            raise TypeError, "x must be of type t_VEC."
        _sig_on
        return new_gen(Strexpand(x.g))


    def Strtex(gen x):
        r"""
        Strtex(x): Translates the vector x of PARI gens to TeX format
        and returns the resulting concatenated strings as a PARI t_STR.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- PARI t_STR (string)
        EXAMPLES:
            sage: v=pari('x^2')
            sage: v.Strtex()
            x^2
            sage: v=pari(['1/x^2','x'])
            sage: v.Strtex()
            \frac{1}{x^2}x
            sage: v=pari(['1 + 1/x + 1/(y+1)','x-1'])
            sage: v.Strtex()
            \frac{ \left(\frac{y
             + 2}{y
             + 1}\right)  x
             + 1}{x}x
             - 1
        """
        if x.type() != 't_VEC':
            x = vector(1, [x])
        _sig_on
        return new_gen(Strtex(x.g))

    def printtex(gen x):
        return x.Strtex()

    def Vec(gen x):
        """
        Vec(x): Transforms the object x into a vector.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- of PARI type t_VEC
        EXAMPLES:
            sage: pari(1).Vec()
            [1]
            sage: pari('x^3').Vec()
            [1, 0, 0, 0]
            sage: pari('x^3 + 3*x - 2').Vec()
            [1, 0, 3, -2]
            sage: pari([1,2,3]).Vec()
            [1, 2, 3]
            sage: pari('ab').Vec()
            [1, 0]
        """
        _sig_on
        return new_gen(gtovec(forcecopy(x.g)))

    def Vecsmall(gen x):
        """
        Vecsmall(x): transforms the object x into a t_VECSMALL.

        INPUT:
            x -- gen
        OUTPUT:
            gen -- PARI t_VECSMALL
        EXAMPLES:
            sage: pari([1,2,3]).Vecsmall()
            Vecsmall([1, 2, 3])
            sage: pari('"SAGE"').Vecsmall()
            Vecsmall([83, 65, 71, 69])
            sage: pari(1234).Vecsmall()
            Vecsmall([1234])
        """
        n = setjmp(GP_DATA.env)
        if n:
            raise TypeError, "Cannot convert %s to Vecsmall."%x
        _sig_on
        return new_gen(gtovecsmall(x.g))

    def binary(gen x):
        """
        binary(x): gives the vector formed by the binary digits of
        abs(x), where x is of type t_INT.

        INPUT:
            x -- gen of type t_INT
        OUTPUT:
            gen -- of type t_VEC
        EXAMPLES:
            sage: pari(0).binary()
            [0]
            sage: pari(-5).binary()
            [1, 0, 1]
            sage: pari(5).binary()
            [1, 0, 1]
            sage: pari(2005).binary()
            [1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1]

            sage.: pari('"2"').binary()
            Traceback (most recent call last):
            ...
            TypeError: x (=2) must be of type t_INT, but is of type t_STR.
        """
        if typ(x.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(x,x.type())
        return new_gen(binaire(x.g))

    def bitand(gen x, y):
        """
        bitand(x,y): Bitwise and of two integers x and y. Negative
        numbers behave as if modulo some large power of 2.

        INPUT:
            x -- gen  (of type t_INT)
            y -- coercible to gen  (of type t_INT)
        OUTPUT:
            gen -- of type type t_INT
        EXAMPLES:
            sage: pari(8).bitand(4)
            0
            sage: pari(8).bitand(8)
            8
            sage: pari(6).binary()
            [1, 1, 0]
            sage: pari(7).binary()
            [1, 1, 1]
            sage: pari(6).bitand(7)
            6
            sage: pari(19).bitand(-1)
            19
            sage: pari(-1).bitand(-1)
            -1
        """
        cdef gen _y
        if typ(x.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(x,x.type())
        _y = pari(y)
        if typ(_y.g) != t_INT:
            raise TypeError, "y (=%s) must be of type t_INT, but is of type %s."%(_y,_y.type())
        return new_gen(gbitand(x.g, _y.g))


    def bitneg(gen x, long n=-1):
        r"""
        bitneg(x,{n=-1}): Bitwise negation of the integer x truncated
        to n bits.  n=-1 (the default) represents an infinite sequence
        of the bit 1.  Negative numbers behave as if modulo some large
        power of 2.

        With n=-1, this function returns -n-1.  With n >= 0, it returns
        a number a such that  $a\cong -n-1 \pmod{2^n}$.

        INPUT:
            x -- gen (t_INT)
            n -- long, default = -1
        OUTPUT:
            gen -- t_INT
        EXAMPLES:
            sage: pari(10).bitneg()
            -11
            sage: pari(1).bitneg()
            -2
            sage: pari(-2).bitneg()
            1
            sage: pari(-1).bitneg()
            0
            sage: pari(569).bitneg()
            -570
            sage: pari(569).bitneg(10)
            454
            sage: 454 % 2**10
            454
            sage: -570 % 2**10
            454
        """
        if typ(x.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(x,x.type())
        return new_gen(gbitneg(x.g,n))


    def bitnegimply(gen x, y):
        """
        bitnegimply(x,y): Bitwise negated imply of two integers x and
        y, in other words, x BITAND BITNEG(y). Negative numbers behave
        as if modulo big power of 2.

        INPUT:
            x -- gen  (of type t_INT)
            y -- coercible to gen  (of type t_INT)
        OUTPUT:
            gen -- of type type t_INT
        EXAMPLES:
            sage: pari(14).bitnegimply(0)
            14
            sage: pari(8).bitnegimply(8)
            0
            sage: pari(8+4).bitnegimply(8)
            4
        """
        cdef gen _y
        if typ(x.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(x,x.type())
        _y = pari(y)
        if typ(_y.g) != t_INT:
            raise TypeError, "y (=%s) must be of type t_INT, but is of type %s."%(_y,_y.type())
        return new_gen(gbitnegimply(x.g, _y.g))


    def bitor(gen x, y):
        """
        bitor(x,y): Bitwise or of two integers x and y. Negative
        numbers behave as if modulo big power of 2.

        INPUT:
            x -- gen  (of type t_INT)
            y -- coercible to gen  (of type t_INT)
        OUTPUT:
            gen -- of type type t_INT
        EXAMPLES:
            sage: pari(14).bitor(0)
            14
            sage: pari(8).bitor(4)
            12
            sage: pari(12).bitor(1)
            13
            sage: pari(13).bitor(1)
            13
        """
        cdef gen _y
        if typ(x.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(x,x.type())
        _y = pari(y)
        if typ(_y.g) != t_INT:
            raise TypeError, "y (=%s) must be of type t_INT, but is of type %s."%(_y,_y.type())
        return new_gen(gbitor(x.g, _y.g))


    def bittest(gen x, long n):
        """
        bittest(x, long n): Returns bit number n (coefficient of $2^n$ in binary)
        of the integer x. Negative numbers behave as if modulo a big power of 2.

        INPUT:
           x -- gen (pari integer)
        OUTPUT:
           bool -- a Python bool
        EXAMPLES:
            sage: x = pari(6)
            sage: x.bittest(0)
            False
            sage: x.bittest(1)
            True
            sage: x.bittest(2)
            True
            sage: x.bittest(3)
            False
            sage: pari(-3).bittest(0)
            True
            sage: pari(-3).bittest(1)
            False
            sage: [pari(-3).bittest(n) for n in range(10)]
            [True, False, True, True, True, True, True, True, True, True]
        """
        if typ(x.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(x,x.type())
        return bool(bittest(x.g, n))

    def bitxor(gen x, y):
        """
        bitxor(x,y): Bitwise exclusive or of two integers x and y.
        Negative numbers behave as if modulo big power of 2.

        INPUT:
            x -- gen  (of type t_INT)
            y -- coercible to gen  (of type t_INT)
        OUTPUT:
            gen -- of type type t_INT
        EXAMPLES:
            sage: pari(6).bitxor(4)
            2
            sage: pari(0).bitxor(4)
            4
            sage: pari(6).bitxor(0)
            6
        """
        cdef gen _y
        if typ(x.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(x,x.type())
        _y = pari(y)
        if typ(_y.g) != t_INT:
            raise TypeError, "y (=%s) must be of type t_INT, but is of type %s."%(_y,_y.type())
        return new_gen(gbitxor(x.g, _y.g))


    def ceil(gen x):
        """
        Return the smallest integer >= x.

        INPUT:
           x -- gen
        OUTPUT:
           gen -- integer
        EXAMPLES:
            sage: pari(1.4).ceil()
            2
            sage: pari(-1.4).ceil()
            -1
            sage: pari('x').ceil()
            x
            sage: pari('x^2+5*x+2.2').ceil()
            x^2 + 5*x + 2.200000000000000000000000000
            sage: pari('3/4').ceil()
            1
        """
        n = setjmp(GP_DATA.env)
        if n == 20:
            raise TypeError, "Ceiling of %s not defined."%x
        return new_gen(gceil(x.g))

    def centerlift(gen x, v=-1):
        """
        centerlift(x,{v}): Centered lift of x.  This function returns
        exactly the same thing as lift, except if x is an integer mod.

        INPUT:
            x -- gen
            v -- var (default: x)
        OUTPUT:
            gen
        EXAMPLES:
            sage: x = pari(-2).Mod(5)
            sage: x.centerlift()
            -2
            sage: x.lift()
            3
            sage: f = pari('x-1').Mod('x^2 + 1')
            sage: f.centerlift()
            x - 1
            sage: f.lift()
            x - 1
            sage: f = pari('x-y').Mod('x^2+1')
            sage: f
            Mod(x - y, x^2 + 1)
            sage: f.centerlift('x')
            x - y
            sage: f.centerlift('y')
            Mod(x - y, x^2 + 1)
        """
        return new_gen(centerlift0(x.g,get_var(v)))


    def changevar(gen x, y):
        """
        changevar(gen x, y): change variables of x according to the vector y.

        INPUT:
            x -- gen
            y -- gen (or coercible to gen)
        OUTPUT:
            gen
        EXAMPLES:
            sage.: pari('X^3+1').changevar(['Y'])
            Y^3 + 1
            sage.: pari('X^3 + Y^2').changevar(['Y','Z'])
            Y^3 + Z^2
            sage.: pari('X^3 + Y^2').changevar(['Y','Z+1'])
            Y^3 + (Z^2 + 2*Z + 1)
        """# The above tests are not run because for some reason
        # they fail when run with all the other tests.  In an
        # interactive session they work fine.  Very frustrating.
        cdef gen _y
        _y = pari(y)
        if typ(_y.g) != t_VEC:
            raise TypeError, "y (=%s) must be of type t_VEC, but is of type %s."%(_y,_y.type())
        _sig_on
        return new_gen(changevar(x.g, _y.g))

    def component(gen x, long n):
        """
        component(x, long n): Return n'th component of the internal
        representation of x.  This this function is 1-based
        instead of 0-based.

        NOTE: For vectors or matrices, it is simpler to use x[n-1]. For
        list objects such as is output by nfinit, it is easier to use
        member functions.

        INPUT:
            x -- gen
            n -- C long (coercible to)
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari([0,1,2,3,4]).component(1)
            0
            sage: pari([0,1,2,3,4]).component(2)
            1
            sage: pari([0,1,2,3,4]).component(4)
            3
            sage: pari('x^3 + 2').component(1)
            2
            sage: pari('x^3 + 2').component(2)
            0
            sage: pari('x^3 + 2').component(4)
            1
            sage: pari('x').component(0)
            Traceback (most recent call last):
            ...
            IndexError: component (=0) must be between 1 and 2, inclusive
        """
        cdef long m
        m = glength(x.g)
        if n <= 0 or n > m:
            raise IndexError, "component (=%s) must be between 1 and %s, inclusive"%(n, m)
        return new_gen(compo(x.g, n))

    def conj(gen x):
        """
        conj(x): Return the algebraic conjugate of x.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('x+1').conj()
            x + 1
            sage: pari('x+I').conj()
            x - I
            sage: pari('1/(2*x+3*I)').conj()
            1/(2*x - 3*I)
            sage: pari([1,2,'2-I','Mod(x,x^2+1)', 'Mod(x,x^2-2)']).conj()
            [1, 2, 2 + I, Mod(-x, x^2 + 1), Mod(-x, x^2 - 2)]
            sage: pari('Mod(x,x^2-2)').conj()
            Mod(-x, x^2 - 2)

            sage.: pari('Mod(x,x^3-3)').conj()
            Traceback (most recent call last):
            ...
            TypeError: x (=Mod(x, x^3 - 3)) has no conjugate.
        """
        n = setjmp(GP_DATA.env)
        if n:
            if n == 20:
                raise TypeError, "x (=%s) has no conjugate."%x
            else:
                _error(n,"error in conj")
        _sig_on
        return new_gen(gconj(x.g))

    def conjvec(gen x):
        """
        conjvec(x): Returns the vector of all conjugates of the
        algebraic number x.  An algebraic number is a polynomial over
        Q modulo an irreducible polynomial.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('Mod(1+x,x^2-2)').conjvec()
            [-0.4142135623730950488016887242, 2.414213562373095048801688724]~
            sage: pari('Mod(x,x^3-3)').conjvec()
            [1.442249570307408382321638311, -0.7211247851537041911608191554 + 1.249024766483406479413179544*I, -0.7211247851537041911608191554 - 1.249024766483406479413179544*I]~
        """
        n = setjmp(GP_DATA.env)
        if n:
            if n == 20:
                raise TypeError, "x (=%s) has no conjugate."%x
            else:
                _error(n,"error in conj")
        _sig_on
        return new_gen(conjvec(x.g, REAL_PREC))

    def denominator(gen x):
        """
        denominator(x): Return the denominator of x.  When x is a
        vector, this is the least common multiple of the denominators
        of the components of x.

        what about poly?
        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('5/9').denominator()
            9
            sage: pari('(x+1)/(x-2)').denominator()
            x - 2
            sage: pari('2/3 + 5/8*x + 7/3*x^2 + 1/5*y').denominator()
            1
            sage: pari('2/3*x').denominator()
            1
            sage: pari('[2/3, 5/8, 7/3, 1/5]').denominator()
            120
        """
        _sig_on
        return new_gen(denom(x.g))

    def floor(gen x):
        """
        floor(x): Return the floor of x, which is the largest integer <= x.
        This function also works component-wise on polynomials, vectors, etc.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('5/9').floor()
            0
            sage: pari('11/9').floor()
            1
            sage: pari('1.17').floor()
            1
            sage: pari('x').floor()
            x
            sage: pari('x+1.5').floor()
            x + 1.500000000000000000000000000
            sage: pari('[1.5,2.3,4.99]').floor()
            [1, 2, 4]
            sage: pari('[[1.1,2.2],[3.3,4.4]]').floor()
            [[1, 2], [3, 4]]

            sage.: pari('"hello world"').floor()
            Traceback (most recent call last):
            ...
            TypeError: x (=hello world) has no floor.
        """
        n = setjmp(GP_DATA.env)
        if n:
            if n == 20:
                raise TypeError, "x (=%s) has no floor."%x
            else:
                _error(n,"error in floor")
        return new_gen(gfloor(x.g))

    def frac(gen x):
        """
        frac(x): Return the fractional part of x, which is x - floor(x).

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('1.7').frac()
            0.7000000000000000000000000000
            sage: pari('sqrt(2)').frac()
            0.4142135623730950488016887242

            sage.: pari('sqrt(-2)').frac()
            Traceback (most recent call last):
            ...
            TypeError: Fractional part of 1.4142135623730950488016887242*I not defined.
        """
        n = setjmp(GP_DATA.env)
        if n:
            if n == 20:
                raise TypeError, "Fractional part of %s not defined."%x
            else:
                _error(n,"error in frac")
        return new_gen(gfrac(x.g))

    def imag(gen x):
        """
        imag(x): Return the imaginary part of x.  This function also
        works component-wise.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('1+2*I').imag()
            2
            sage: pari('sqrt(-2)').imag()
            1.414213562373095048801688724
            sage: pari('x+I').imag()
            1
            sage: pari('x+2*I').imag()
            2
            sage: pari('(1+I)*x^2+2*I').imag()
            x^2 + 2
            sage: pari('[1,2,3] + [4*I,5,6]').imag()
            [4, 0, 0]
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            if n == 20:
                raise TypeError, "Imaginary part of %s not defined."%x
            else:
                _error(n,"error in imag")
        return new_gen(gimag(forcecopy(forcecopy(x.g))))

    def length(gen x):
        """
        length(x): Return the number of non-code words in x.  If x
        is a string, this is the number of characters of x.

        ?? terminator ?? carriage return ??
        """
        return glength(x.g)

    def lift(gen x, v=-1):
        """
        lift(x,{v}): Returns the lift of an element of Z/nZ to Z or
        R[x]/(P) to R[x] for a type R if v is omitted.  If v is given,
        lift only polymods with main variable v.  If v does not occur
        in x, lift only intmods.

        INPUT:
            x -- gen
            v -- (optional) variable
        OUTPUT:
            gen
        EXAMPLES:
            sage: x = pari("x")
            sage: a = x.Mod(x**3 + 17*x + 3)
            sage: a
            Mod(x, x^3 + 17*x + 3)
            sage: b = a**4; b
            Mod(-17*x^2 - 3*x, x^3 + 17*x + 3)
            sage: b.lift()
            -17*x^2 - 3*x

        ??? more examples
        """
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in lift")
        if v == -1:
            return new_gen(lift(forcecopy(x.g)))
        return new_gen(lift0(forcecopy(x.g), get_var(v)))

    def numerator(gen x):
        """
        numerator(x): Returns the numerator of x.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:

        """
        raise new_gen(numer(x.g))


    def numtoperm(gen k, long n):
        """
        numtoperm(k, n): Return the permutation number k (mod n!) of n
        letters, where n is an integer.

        INPUT:
            k -- gen, integer
            n -- int
        OUTPUT:
            gen -- vector (permutation of {1,...,n})
        EXAMPLES:
        """
        _sig_on
        return new_gen(numtoperm(n, k.g))


    def padicprec(gen x, p):
        """
        padicprec(x,p): Return the absolute p-adic precision of the object x.

        INPUT:
            x -- gen
        OUTPUT:
            int
        EXAMPLES:
        """
        cdef gen _p
        _p = pari(p)
        if typ(_p.g) != t_INT:
            raise TypeError, "p (=%s) must be of type t_INT, but is of type %s."%(
                _p, _p.type())
        return padicprec(x.g, _p.g)

    def permtonum(gen x):
        """
        permtonum(x): Return the ordinal (between 1 and n!) of permutation vector x.
        ??? Huh ???  say more.  what is a perm vector.  0 to n-1 or 1-n.

        INPUT:
            x -- gen (vector of integers)
        OUTPUT:
            gen -- integer
        EXAMPLES:
        """
        if typ(x.g) != t_VEC:
            raise TypeError, "x (=%s) must be of type t_VEC, but is of type %s."%(x,x.type())
        raise new_gen(permtonum(x.g))

    def precision(gen x, long n=-1):
        """
        precision(x,{n}): Change the precision of x to be n, where n
        is a C-integer). If n is omitted, output the real precision of x.

        INPUT:
            x -- gen
            n -- (optional) int
        OUTPUT:
            nothing
          or
            gen if n is omitted
        EXAMPLES:
        """
        if n <= -1:
            return precision(x.g)
        return new_gen(precision0(x.g, n))

    def random(gen N):
        r"""
        \code{random(\{N=$2^31$\})}: Return a pseudo-random integer between 0 and $N-1$.

        INPUT:
            N -- gen, integer
        OUTPUT:
            gen -- integer
        EXAMPLES:
        """
        if typ(N.g) != t_INT:
            raise TypeError, "x (=%s) must be of type t_INT, but is of type %s."%(N,N.type())
        _sig_on
        return new_gen(genrand(N.g))

    def real(gen x):
        """
        real(x): Return the real part of x.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in real")
        return new_gen(greal(forcecopy(x.g)))

    def round(gen x, estimate=False):
        """
        round(x,estimat=False):  If x is a real number, returns x rounded
        to the nearest integer (rounding up).  If the optional argument
        estimate is True, also returns the binary exponent e of the difference
        between the original and the rounded value (the "fractional part")
        (this is the integer ceiling of log_2(error)).

        When x is a general PARI object, this function returns the result
        of rounding every coefficient at every level of PARI object.
        Note that this is different than what the truncate function
        does (see the example below).

        One use of round is to get exact results after a long
        approximate computation, when theory tells you that the
        coefficients must be integers.

        INPUT:
            x -- gen
            estimate -- (optional) bool, False by default
        OUTPUT:
            * if estimate == False, return a single gen.
            * if estimate == True, return rounded verison of x and
              error estimate in bits, both as gens.
        EXAMPLES:
            sage: pari('1.5').round()
            2
            sage: pari('1.5').round(True)
            (2, -1)
            sage: pari('1.5 + 2.1*I').round()
            2 + 2*I
            sage: pari('1.0001').round(True)
            (1, -14)
            sage: pari('(2.4*x^2 - 1.7)/x').round()
            (2*x^2 - 2)/x
            sage: pari('(2.4*x^2 - 1.7)/x').truncate()
            2.400000000000000000000000000*x
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in round")
        if not estimate:
            return new_gen(ground(x.g))
        cdef long e
        cdef gen y
        y = new_gen(grndtoi(x.g, &e))
        return y, e

    def simplify(gen x):
        """
        simplify(x): Simplify the object x as much as possible, and return
        the result.

        A complex or quadratic number whose imaginary part is an exact 0
        (i.e., not an approximate one such as O(3) or 0.E-28) is converted
        to its real part, and a a polynomial of degree 0 is converted to
        its constant term.  Simplification occurs recursively.

        This function is useful before using arithmetic functions, which
        expect integer arguments:

        EXAMPLES:
            sage: y = pari('y')
            sage: x = pari('9') + y - y
            sage: x
            9
            sage: x.type()
            't_POL'
            sage: x.factor()
            matrix(0,2)
            sage: pari('9').factor()
            Mat([3, 2])
            sage: x.simplify()
            9
            sage: x.simplify().factor()
            Mat([3, 2])
            sage: x = pari('1.5 + 0*I')
            sage: x.type()
            't_COMPLEX'
            sage: x.simplify()
            1.500000000000000000000000000
            sage: y = x.simplify()
            sage: y.type()
            't_REAL'
        """
        _sig_on
        return new_gen(simplify(x.g))

    def sizebyte(gen x):
        """
        sizebyte(x): Return the total number of bytes occupied by the
        complete tree of the object x.

        INPUT:
            x -- gen
        OUTPUT:
            int (a Python int)
        EXAMPLES:
            sage: pari('1').sizebyte()
            12
            sage: pari('10').sizebyte()
            12
            sage: pari('10000000000000').sizebyte()
            16
            sage: pari('10^100').sizebyte()
            52
            sage: pari('x').sizebyte()
            36
            sage: pari('x^20').sizebyte()
            264
            sage: pari('[x, 10^100]').sizebyte()
            100
        """
        return taille2(x.g)

    def sizedigit(gen x):
        """

        sizedigit(x): Return a quick estimate for the maximal number of
        decimal digits before the decimal point of any component of x.

        INPUT:
            x -- gen
        OUTPUT:
            int -- Python integer
        EXAMPLES:
            sage: x = pari('10^100')
            sage: x.Str().length()
            101
            sage: x.sizedigit()
            101

        Note that digits after the decimal point are ignored.
            sage: x = pari('1.234')
            sage: x
            1.234000000000000000000000000
            sage: x.sizedigit()
            1

        The estimate can be one too big:
            sage: pari('7234.1').sizedigit()
            4
            sage: pari('9234.1').sizedigit()
            5
        """
        return sizedigit(x.g)

    def truncate(gen x, estimate=False):
        """
        truncate(x,estimate=False):  Return the truncation of x.
        If estimate is True, also return the number of error bits.

        When x is in the real numbers, this means that the part
        after the decimal point is chopped away, e is the binary
        exponent of the difference between the original and truncated
        value (the "fractional part").    If x is a rational
        function, the result is the integer part (Euclidean
        quotient of numerator by denominator) and if requested
        the error estimate is 0.

        When truncate is applied to a power series (in X), it
        transforms it into a polynomial or a rational function with
        denominator a power of X, by chopping away the $O(X^k)$.
        Similarly, when applied to a p-adic number, it transforms it
        into an integer or a rational number by chopping away the
        $O(p^k)$.

        INPUT:
            x -- gen
            estimate -- (optional) bool, which is False by default
        OUTPUT:
            * if estimate == False, return a single gen.
            * if estimate == True, return rounded verison of x and
              error estimate in bits, both as gens.
        OUTPUT:
        EXAMPLES:
            sage: pari('(x^2+1)/x').round()
            (x^2 + 1)/x
            sage: pari('(x^2+1)/x').truncate()
            x
            sage: pari('1.043').truncate()
            1
            sage: pari('1.043').truncate(True)
            (1, -5)
            sage: pari('1.6').truncate()
            1
            sage: pari('1.6').round()
            2
            sage: pari('1/3 + 2 + 3^2 + O(3^3)').truncate()
            34/3
            sage: pari('sin(x+O(x^10))').truncate()
            1/362880*x^9 - 1/5040*x^7 + 1/120*x^5 - 1/6*x^3 + x
            sage: pari('sin(x+O(x^10))').round()   # each coefficient has abs < 1
            x + O(x^10)
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in truncate")
        if not estimate:
            return new_gen(gtrunc(x.g))
        cdef long e
        cdef gen y
        y = new_gen(gcvtoi(x.g, &e))
        return y, e

    def valuation(gen x, p):
        """
        valuation(x,p): Return the valuation of x with respect to p.

        The valuation is the highest exponent of p dividing x.

           * If p is an integer, x must be an integer, an intmod whose
             modulus is divisible by p, a rational number, a p-adic
             number, or a polynomial or power series in which case the
             valuation is the minimal of the valuations of the
             coefficients.

           * If p is a polynomial, x must be a polynomial or a
             rational fucntion.  If p is a monomial then x may also be
             a power series.

           * If x is a vector, complex or quadratic number, then the
             valuation is the minimum of the component valuations.

           * If x = 0, the result is $2^31-1$ on 32-bit machines or
             $2^63-1$ on 64-bit machines if x is an exact object.
             If x is a p-adic number or power series, the result
             is the exponent of the zero.

        INPUT:
            x -- gen
            p -- coercible to gen
        OUTPUT:
            gen -- integer
        EXAMPLES:
            sage: pari(9).valuation(3)
            2
            sage: pari(9).valuation(9)
            1
            sage: x = pari(9).Mod(27); x.valuation(3)
            2
            sage: pari('5/3').valuation(3)
            -1
            sage: pari('9 + 3*x + 15*x^2').valuation(3)
            1
            sage: pari([9,3,15]).valuation(3)
            1
            sage: pari('9 + 3*x + 15*x^2 + O(x^5)').valuation(3)
            1

            sage: pari('x^2*(x+1)^3').valuation(pari('x+1'))
            3
            sage: pari('x + O(x^5)').valuation('x')
            1
            sage: pari('2*x^2 + O(x^5)').valuation('x')
            2

            sage.: pari(0).valuation(3)   # on 32-bit machine
            2147483647
        """
        cdef gen _p
        _p = pari(p)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in valuation")
        return ggval(x.g, _p.g)

    def variable(gen x):
        """
        variable(x): Return the main variable of the object x, or p
        if x is a p-adic number.

        This function raises a TypeError exception on scalars, i.e.,
        on objects with no variable associated to them.

        INPUT:
            x -- gen
        OUTPUT:
            gen
        EXAMPLES:
            sage: pari('x^2 + x -2').variable()
            x
            sage: pari('1+2^3 + O(2^5)').variable()
            2
            sage: pari('x+y0').variable()
            x
            sage: pari('y0+z0').variable()
            y0
        """
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in valuation")
        return new_gen(gpolvar(x.g))


    ###########################################
    # 3: TRANSCENDENTAL functions
    ###########################################
    def sqrt(gen x, int prec=0):
        cdef int err
        err = setjmp(GP_DATA.env)
        if err: _error(err,"error in sqrt")
        if prec < 3:
            prec = REAL_PREC
        _sig_on
        return new_gen(gsqrt(forcecopy(x.g), prec))

    def gamma(gen s, int prec=0):
        """
        gamma(x): Gamma function at s.
        INPUT:
            s -- gen (real or complex number
        OUTPUT:
            gen -- value of zeta at s.
        EXAMPLES:
        """
        cdef int err
        err = setjmp(GP_DATA.env)
        if err: _error(err,"error evaluating gamma(%s)"%s)
        if prec < 3:
            prec = REAL_PREC
        _sig_on
        return new_gen(ggamma(forcecopy(s.g), prec))


    def zeta(gen s, int prec=0):
        """
        zeta(s): Riemann zeta function at s with s a complex
                 or a p-adic number.
        INPUT:
            s -- gen (real or complex number
        OUTPUT:
            gen -- value of zeta at s.
        EXAMPLES:
        """
        cdef int err
        err = setjmp(GP_DATA.env)
        if err: _error(err,"error evaluating zeta(%s)"%s)
        if prec < 3:
            prec = REAL_PREC
        _sig_on
        return new_gen(gzeta(forcecopy(s.g), prec))


    ###########################################
    # 4: NUMBER THEORETICAL functions
    ###########################################

    def issquare(gen x, find_root=False):
        cdef int err
        cdef GEN G, t
        cdef gen g
        err = setjmp(GP_DATA.env)
        if err:
            _error(err,"error in issquare")
        if find_root:
            _sig_on
            t = gcarrecomplet(forcecopy(x.g), &G)
            _sig_off
            v = bool(new_gen_noclear(t))
            if v:
                return v, new_gen(G)
            else:
                return v, None
        else:
            _sig_on
            return new_gen(gcarreparfait(x.g))


    def issquarefree(gen self):
        _sig_on
        t = bool(issquarefree(forcecopy(self.g)))
        _sig_off
        return t

    def lcm(gen x, y):
        cdef gen _y
        _y = pari(y)
        cdef int err
        err = setjmp(GP_DATA.env)
        if err: _error(err,"error in lcm")
        _sig_on
        return new_gen(glcm(forcecopy(x.g), forcecopy(_y.g)))

    def gcd(gen x, y):
        cdef gen _y
        _y = pari(y)
        cdef int err
        err = setjmp(GP_DATA.env)
        if err: _error(err,"error in gcd")
        _sig_on
        return new_gen(ggcd(forcecopy(x.g), forcecopy(_y.g)))

    def bezout(gen x, y):
        cdef gen _y, u, v, g
        cdef GEN U, V, G
        _y = pari(y)
        cdef int err
        err = setjmp(GP_DATA.env)
        if err: _error(err,"error in bezout")
        _sig_on
        G = gbezout(forcecopy(x.g), forcecopy(_y.g), &U, &V)
        _sig_off
        g = new_gen_noclear(G)
        u = new_gen_noclear(U)
        v = new_gen(V)
        return g, u, v

    def xgcd(gen x, y):
        return x.bezout(y)

    ###########################################
    # 5: Functions related to ELLIPTIC CURVES
    ###########################################

    def elladd(gen self, z1, z2):
        """
        e.elladd(z1, z2)

        Sum of the points z1 and z2 on the elliptic curve e.

        INPUT:
            e -- elliptic curve E
            z1 -- point on E
            z2 -- point on E

        OUTPUT:
            point on E

        EXAMPLES:
        First we create an elliptic curve:

            sage: e = pari([0, 1, 1, -2, 0]).ellinit()
            sage: str(e)[:65]   # first part of output
            '[0, 1, 1, -2, 0, 4, -4, 1, -3, 112, -856, 389, 1404928/389, [0.90'

        Next we add two points on the elliptic curve.  Notice that
        the Python lists are automatically converted to PARI objects so
        you don't have to do that explicitly in your code.

            sage: e.elladd([1,0,1], [-1,1,1])
            [-3/4, -15/8]
        """
        cdef gen _z1, _z2
        _z1 = pari(z1); _z2 = pari(z2)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in elladd")
        _sig_on
        return new_gen(addell(forcecopy(self.g), _z1.g, _z2.g))

    def ellak(gen self, n):
        r"""
        e.ellak(n): Returns the coefficient $a_n$ of the $L$-function of
        the elliptic curve e, i.e. the coefficient of a newform of
        weight 2 newform.

        \begin{notice}
        The curve $e$ {\em must} be a medium or long vector of the type given
        by ellinit. For this function to work for every n and not just
        those prime to the conductor, e must be a minimal Weierstrass
        equation. If this is not the case, use the function
        ellminimalmodel first before using ellak (or you will get
        INCORRECT RESULTS!)
        \end{notice}

        INPUT:
            e -- a PARI elliptic curve.
            n -- integer ..

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.ellak(6)
            2
            sage: e.ellak(2005)
            2
            sage: e.ellak(-1)
            0
            sage: e.ellak(0)
            0
        """
        cdef gen _n
        _n = pari(n)
        cdef int r
        r = setjmp(GP_DATA.env)
        if r: _error(n,"error in ellak")
        _sig_on
        return new_gen(akell(forcecopy(self.g), _n.g))

    def ellan(gen self, long n):
        """
        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.ellan(3)
            [1, -2, -1]
            sage: e.ellan(20)
            [1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1, -2, 4, 4, -1, -4, -2, 4, 0, 2]
            sage: e.ellan(-1)
            []
        """
        cdef int r
        r = setjmp(GP_DATA.env)
        if r == MEMERR: # out of memory -- increase and try again
            init_stack(0)
            return self.ellan(n)
        elif r:
            _error(n,"error in ellan")
        _sig_on
        return new_gen(anell(forcecopy(self.g), n))

    def ellap(gen self, p):
        r"""
        e.ellap(p): Returns the prime-indexed coefficient $a_p$ of the
        $L$-function of the elliptic curve $e$, i.e. the coefficient of a
        newform of weight 2 newform.

        \begin{notice}
        If p is not prime, this function will return an {\bf incorrect}
        answer.

        The curve e must be a medium or long vector of the type given
        by ellinit. For this function to work for every n and not just
        those prime to the conductor, e must be a minimal Weierstrass
        equation. If this is not the case, use the function
        ellminimalmodel first before using ellak (or you will get
        INCORRECT RESULTS!)
        \end{notice}

        INPUT:
            e -- a PARI elliptic curve.
            p -- prime integer ..

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.ellap(2)
            -2
            sage: e.ellap(2003)
            4
            sage: e.ellak(-1)
            0
        """
        cdef gen _p
        _p = pari(p)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in apell, p='%s' not a prime"%p)
        _sig_on
        return new_gen(apell(forcecopy(self.g), _p.g))


    def ellbil(gen self, z1, z2):
        """
        EXAMPLES:
            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: e.ellbil([1, 0, 1], [-1, 1, 1])
            0.4181889844988605856298894585
            sage: get_real_precision()
            5
        """
##         Increasing the precision does not increase the precision
##         result, since quantities related to the elliptic curve were
##         computed to low precision.
##             sage: set_real_precision(10)
##             sage: e.ellbil([1, 0, 1], [-1, 1, 1])
##             0.4181889844988605856298894585
##         However, if we recompute the elliptic curve after increasing
##         the precision, then the bilinear pairing will be computed to
##         higher precision as well.
##             sage: e = pari([0,1,1,-2,0]).ellinit()
##             sage: e.ellbil([1, 0, 1], [-1, 1, 1])
##             0.4181889844988605856298894582
##             sage: set_real_precision(5)
        cdef gen _z1, _z2
        _z1 = pari(z1)
        _z2 = pari(z2)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in ellbil")
        _sig_on
        return new_gen(bilhell(forcecopy(self.g), _z1.g, _z2.g, REAL_PREC))

    def ellchangecurve(gen e, ch):
        """
        EXAMPLES:
            sage: e = pari([1,2,3,4,5]).ellinit()
            sage: e.ellglobalred()
            [10351, [1, -1, 0, -1], 1]
            sage: f = e.ellchangecurve([1,-1,0,-1])
            sage: f[:5]
            [1, -1, 0, 4, 3]
        """
        cdef gen _ch
        _ch = pari(ch)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in ellchangecurve")
        _sig_on
        return new_gen(coordch(forcecopy(e.g), _ch.g))

    def ellchangepoint(gen x, y):
        """
        x.ellchangepoint(y): change data on point or vector of points x
                             on an elliptic curve according to y=[u,r,s,t]

        EXAMPLES:
            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: x = pari([1,0,1])
            sage: e.ellisoncurve([1,4,4])
            False
            sage: e.ellisoncurve(x)
            True
            sage: f = e.ellchangecurve([1,2,3,-1])
            sage: f[:5]   # show only first five entries
            [6, -2, -1, 17, 8]
            sage: x.ellchangepoint([1,2,3,-1])
            [-1, 4]
            sage: f.ellisoncurve([-1,4])
            True
        """
        cdef gen _y
        _y = pari(y)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in ellchangepoint")
        _sig_on
        return new_gen(pointch(forcecopy(x.g), _y.g))

    def elleisnum(gen om, long k, int flag=0):
        """
        om.elleisnum(k, {flag=0}, {prec}): om=[om1,om2] being a
            2-component vector giving a basis of a lattice L and k an
            even positive integer, computes the numerical value of the
            Eisenstein series of weight k. When flag is non-zero and
            k=4 or 6, this gives g2 or g3 with the correct
            normalization.

        INPUT:
            om -- gen, 2-component vector giving a basis of a lattice L
            k  -- int (even positive)
            flag -- int (default 0)
            pref -- precision

        OUTPUT:
            gen -- numerical value of E_k

        EXAMPLES:
            sage: e = pari([0,1,1,-2,0]).ellinit()
            sage: om = e.omega()
            sage: om
            [2.490212560855055075321357792, 1.971737701551648204422407698*I]
            sage: om.elleisnum(2)
            -5.288649339654257621791534695
            sage: om.elleisnum(4)
            112.0000000000000000000000000
            sage: om.elleisnum(100)
            2.153142485760776361127070332 E50
            sage: om.elleisnum(5)
            Traceback (most recent call last):
            ...
            ArithmeticError: k(=5) must be even
            sage: om.elleisnum(0)
            Traceback (most recent call last):
            ...
            ArithmeticError: k(=0) must be positive
        """
        if k <= 0:
            raise ArithmeticError, "k(=%s) must be positive"%k
        if k % 2:
            raise ArithmeticError, "k(=%s) must be even"%k
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in elleisnum: invalid input")
        _sig_on
        return new_gen(elleisnum(forcecopy(om.g), k, flag, REAL_PREC))

    def elleta(gen self):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in elleta")
        _sig_on
        return new_gen(elleta(forcecopy(self.g), REAL_PREC))

    def ellglobalred(gen e1):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in ellglobalred")
        _sig_on
        return new_gen(globalreduction(forcecopy(e1.g)))

    def ellheight(gen e, a, flag=0):
        cdef gen _a
        _a = pari(a)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in ellheight")
        _sig_on
        return new_gen(ellheight0(forcecopy(e.g), _a.g, flag, REAL_PREC))

    def ellheightmatrix(gen e, x):
        """
        ellheightmatrix(e,x)

        Returns the height matrix for vector of points x on elliptic curve e using
        theta functions.
        """
        cdef gen _x
        _x = pari(x)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in ellheightmatrix")
        _sig_on
        return new_gen(mathell(forcecopy(e.g), _x.g, REAL_PREC))

    def ellinit(gen self, int flag=0, prec=0):
        r"""
        ellinit(x,{flag=0}): x being the vector [a1,a2,a3,a4,a6], gives
        the vector:
        \begin{verbatim}
                   [a1,a2,a3,a4,a6,
                    b2,b4,b6,b8,
                    c4,c6,
                    delta,j,
                    [e1,e2,e3],
                    w1,w2,eta1,eta2,area].
        \end{verbatim}
        If the curve is defined over a p-adic field, the last six components
        are replaced by root, $u^2$, u, q, w, 0. If optional flag is 1, omit them
        altogether.
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in ellinit: invalid input '%s' (flag=%s) to the elliptic curve creation function"%(self, flag))
        if prec < 3:
            prec = REAL_PREC
        if flag < 0 or flag >= 1:
            raise ValueError, "flag must be 0 or 1"
        _sig_on
        return new_gen(ellinit0(forcecopy(self.g), flag, prec))

    def ellisoncurve(gen e, x):
        cdef gen _x
        _x = pari(x)

        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in elleisoncurve: invalid input e=%s, x=%s"%(e,x))
        _sig_on
        t = bool(oncurve(e.g, _x.g) == 1)
        _sig_off
        return t

    def ellj(gen e):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in ellj: invalid input '%s', argument must be element of the upper half plane"%e)
        _sig_on
        return new_gen(jell(forcecopy(e.g), REAL_PREC))

    def elllocalred(gen e, p):
        cdef gen _p
        _p = pari(p)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in ellj: invalid input '%s', argument must be element of the upper half plane"%e)
        #_sig_on
        return new_gen(elllocalred(forcecopy(e.g), _p.g))

    def elllseries(gen e, s, A=1, long prec=0):
        cdef gen _s, _A
        _s = pari(s)
        _A = pari(A)
        if prec < 3:
            prec = REAL_PREC

        # The following line is a hack to get around a bug
        # in PARI 2.2.9-alpha, where L(E,n) crashes for n<=0 an integer.
        if s <= ZERO and int(s) == float(s): return ZERO

        n = setjmp(GP_DATA.env)
        if n:
            # I think this problem occurs only at nonpositive integer
            # points.  At these points MAYBE(?) the L-series always
            # vanishes, so just returning 0 would be a trivial fix.
            _error(n,"in elllseries: Error evaluating L(E,%s)"%s)

        _sig_on
        return new_gen(lseriesell(forcecopy(e.g), _s.g, _A.g, prec))

    def ellminimalmodel(gen e):
        """

        ellminimalmodel(e): return the standard minimal integral model
        of the rational elliptic curve e and the corresponding change
        of variables.
        INPUT:
            e -- gen (that defines an elliptic curve)
        OUTPUT:
            gen -- minimal model
            gen -- change of coordinates
        EXAMPLES:
            sage: e = pari([1,2,3,4,5]).ellinit()
            sage: F, ch = e.ellminimalmodel()
            sage: F[:5]
            [1, -1, 0, 4, 3]
            sage: ch
            [1, -1, 0, -1]
            sage: e.ellchangecurve(ch)[:5]
            [1, -1, 0, 4, 3]
        """
        cdef GEN x, y
        cdef gen model, change
        cdef pari_sp t
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in ellminimalmodel")
        _sig_on
        x = ellminimalmodel(forcecopy(e.g), &y)
        model = new_gen_noclear(x)
        change = new_gen(y)
        return model, change

    def ellorder(gen e, x):
        cdef gen _x
        _x = pari(x)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in ellorder")
        _sig_on
        return new_gen(orderell(forcecopy(e.g), _x.g))

    def ellordinate(gen e, x):
        cdef gen _x
        _x = pari(x)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in ellordinate")
        _sig_on
        return new_gen(ordell(forcecopy(e.g), _x.g, REAL_PREC))

    def ellpointtoz(gen e, z):
        cdef gen _x
        _x = pari(z)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in ellpointtoz")
        _sig_on
        return new_gen(zell(forcecopy(e.g), _x.g, REAL_PREC))

    def ellpow(gen e, z, n):
        cdef gen _z, _n
        _z = pari(z); _n = pari(n)
        cdef int r
        r = setjmp(GP_DATA.env)
        if r: _error(r,"in ellpow")
        _sig_on
        return new_gen(powell(forcecopy(e.g), _z.g, _n.g))

    def ellrootno(gen e, p=1):
        cdef gen _p
        _p = pari(p)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in ellrootno")
        _sig_on
        return ellrootno(forcecopy(e.g), _p.g)

    def ellsigma(gen om, z, flag=0):
        cdef gen _x
        _x = pari(z)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in ellsigma")
        _sig_on
        return new_gen(ellsigma(forcecopy(om.g), _x.g, flag, REAL_PREC))

    def ellsub(gen e, z1, z2):
        cdef gen _z1, _z2
        _z1 = pari(z1); _z2 = pari(z2)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in ellsub")
        _sig_on
        return new_gen(subell(forcecopy(e.g), _z1.g, _z2.g))

    def elltaniyama(gen e):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in elltaniyama")
        _sig_on
        return new_gen(taniyama(forcecopy(e.g)))

    def elltors(gen e, flag=0):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            if n == 22:
                # precision too low error
                raise RuntimeError, "Elliptic curve was created to too low of precision.  Call ellinit with the real_precision higher."
            _error(n,"in elltors")
        _sig_on
        return new_gen(elltors0(forcecopy(e.g), flag))

    def ellwp(gen e, z='z', long n=20, long flag=0):
        """
        ellwp(E, z,{flag=0}): Return the complex value of the Weierstrass
        P-function at z on the lattice defined by e.

        INPUT:
            E -- list OR elliptic curve
                  list -- [om1, om2], which are Z-generators for a lattice
                  elliptic curve -- created using ellinit

            z -- (optional) complex number  OR string (default = "z")
                   complex number -- any number in the complex plane
                   string (or PARI variable) -- name of a variable.

            n -- int (optional: default 20) if z is a variable, compute up to at least o(z^n).

            flag -- int: 0 (default): compute only P(z)
                         1 compute [P(z),P'(z)]
                         2 consider om or E as an elliptic curve and use P-function to
                           compute the point on E (with the Weierstrass equation for E)
                           P(z)
                           for that curve (identical to ellztopoint in this case).

        OUTPUT:
            gen -- complex number or list of two complex numbers

        EXAMPLES:

        We first define the elliptic curve X_0(11):
            sage: E = pari([0,-1,1,-10,-20]).ellinit()

        Compute P(1).
            sage: E.ellwp(1)
            13.96586952574849779469497769 + 2.660659330 E-28*I

        Compute P(1+I), where I = sqrt(-1).
            sage: E.ellwp(pari("1+I"))
            -1.115106825655550879209066487 + 2.334190523074698836184798800*I
            sage: E.ellwp("1+I")
            -1.115106825655550879209066487 + 2.334190523074698836184798800*I

        The series expansion, to the default 20 precision:
            sage: E.ellwp()
            z^-2 + 31/15*z^2 + 2501/756*z^4 + 961/675*z^6 + 77531/41580*z^8 + 1202285717/928746000*z^10 + 2403461/2806650*z^12 + 30211462703/43418875500*z^14 + 3539374016033/7723451736000*z^16 + 413306031683977/1289540602350000*z^18 + O(z^20)

        Compute the series for wp to lower precision:
            sage: E.ellwp(n=4)
            z^-2 + 31/15*z^2 + O(z^4)

        Next we use the version where the input is generators for a lattice:
            sage: pari([1.2692, '0.63 + 1.45*I']).ellwp(1)
            13.96561469366894364802003736 + 0.0006448292728105361474541635979*I

        With flag 1 compute the pair P(z) and P'(z):
            sage: E.ellwp(1, flag=1)
            [13.96586952574849779469497769 + 2.660659330 E-28*I, 50.56193008800727525558465689 + 1.615587134 E-27*I]

        With flag=2, the computed pair is (x,y) on the curve instead of [P(z),P'(z)]:
            sage: E.ellwp(1, flag=2)
            [14.29920285908183112802831103 + 2.660659330 E-28*I, 50.06193008800727525558465689 + 1.615587134 E-27*I]
        """
        cdef gen _z
        _z = pari(z)
        cdef int err
        err = setjmp(GP_DATA.env)
        if err: _error(err,"in ellwp")
        if n < 0:
            n = 0
        if n%2 == 1:
            n = n + 1
        _sig_on
        return new_gen(ellwp0(forcecopy(e.g), _z.g, flag, REAL_PREC, (n+2)/2))

    def ellzeta(gen om, z):
        cdef gen _z
        _z = pari(z)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in ellzeta")
        _sig_on
        return new_gen(ellzeta(forcecopy(om.g), _z.g, REAL_PREC))

    def ellztopoint(gen e, gen z):
        cdef gen _z
        _z = pari(z)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"in ellztopoint")
        _sig_on
        return new_gen(pointell(forcecopy(e.g), _z.g, REAL_PREC))

    def omega(gen e):
        """
        e.omega(): return basis for the period lattice of the elliptic curve e.

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.omega()
            [1.269209304279553421688794617, 0.6346046521397767108443973084 + 1.458816616938495229330889612*I]
        """
        return e[14:16]

    def disc(gen e):
        """
        e.disc(): return the discriminant of the elliptic curve e.

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.disc()
            -161051
            sage: _.factor()
            [-1, 1; 11, 5]
        """
        return e[11]

    def j(gen e):
        """
        e.j(): return the j-invariant of the elliptic curve e.

        EXAMPLES:
            sage: e = pari([0, -1, 1, -10, -20]).ellinit()
            sage: e.j()
            -122023936/161051
            sage: _.factor()
            [-1, 1; 2, 12; 11, -5; 31, 3]
        """
        return e[12]

    ###########################################
    # 6: Functions related to general NUMBER FIELDS
    ###########################################
    def bnfunit(gen k):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in bnfunit")
        _sig_on
        return new_gen(buchfu(forcecopy(k.g)))

    def bnfinit(gen f, long flag=0):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in nfinit")

        _sig_on
        return new_gen(bnfinit0(forcecopy(f.g), flag, <GEN>0, REAL_PREC))

    def idealfactor(gen nf, x):
        cdef gen _x
        _x = pari(x)
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in idealfactor")
        _sig_on
        return new_gen(idealfactor(forcecopy(nf.g), _x.g))

    def nfbasis(gen f, long flag=0, p=0):
        cdef gen _p
        cdef GEN g
        if p != 0:
            _p = pari(p)
            g = _p.g
        else:
            g = <GEN>NULL

        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in nfbasis")

        _sig_on
        return new_gen(nfbasis0(forcecopy(f.g), flag, g))

    def nfinit(gen f, long flag=0):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in nfinit")
        _sig_on
        return new_gen(nfinit0(forcecopy(f.g), flag, REAL_PREC))

    def polcompositum(gen pol1, pol2, long flag=0):
        cdef gen _pol2
        _pol2 = pari(pol2)

        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in polcompositum")
        _sig_on
        return new_gen(polcompositum0(forcecopy(pol1.g), _pol2.g, flag))


    ###########################################
    # 7: POLYNOMIALS and power series
    ###########################################
    def poldegree(gen x, var=-1):
        """
        f.poldegree(var={x}): Return the degree of this polynomial.
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in poldegree")
        return poldegree(x.g, get_var(var))

    def poldisc(gen x, var=-1):
        """
        f.poldist(var={x}):  Return the discriminant of this polynomial.
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in poldisc")
        _sig_on
        return new_gen(poldisc0(forcecopy(x.g), get_var(var)))

    def polgalois(gen x):
        """
        f.polgalois(): Galois group of the polynomial f
        """
        cdef int err
        err = setjmp(GP_DATA.env)
        if err:
            _error(err,"error in polgalois")
        _sig_on
        return new_gen(galois(forcecopy(x.g), REAL_PREC))

    def polisirreducible(gen x):
        """
        f.polisirreducible(): Returns True if f is an irreducible
        non-constant polynomial, or False if f is reducible or
        constant.
        """
        cdef int err
        if typ(x.g) != t_POL:
            raise TypeError, "x (=%s) must be a polynomial."%x
        err = setjmp(GP_DATA.env)
        if err:
            _error(err,"error in polisirreducible")
        _sig_on
        return bool(new_gen(gisirreducible(forcecopy(x.g))))


    def polresultant(gen x, y, var=-1, flag=0):
        cdef gen _y
        _y = pari(y)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in polresultant")

        _sig_on
        return new_gen(polresultant0(forcecopy(x.g), _y.g, get_var(var), flag))

    def polroots(gen x, flag=0):
        """
        polroots(x,{flag=0}): complex roots of the polynomial x. flag
        is optional, and can be 0: default, uses Schonhage's method modified
        by Gourdon, or 1: uses a modified Newton method.
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in polroots")

        _sig_on
        return new_gen(roots0(forcecopy(x.g), flag, REAL_PREC))


    ###########################################
    # 8: Vectors, matrices, LINEAR ALGEBRA and sets
    ###########################################

    ###########################################
    # 9: SUMS, products, integrals and similar functions
    ###########################################

    ###########################################
    # 11. PROGRAMMING under GP
    ###########################################

    ###########################################
    # polarit2.c
    ###########################################
    def factor(gen self, limit=-1):
        """
        factor(x,{lim}): factorization of x. lim is optional and can
        be set whenever x is of (possibly recursive) rational
        type. If lim is set return partial factorization, using
        primes up to lim (up to primelimit if lim=0)
        """
        _sig_on

        if sigsetjmp(_env,1): raise KeyboardInterrupt

        cdef int n

        _sig_on
        return new_gen(factor0(forcecopy(self.g), limit))





    ###########################################
    # misc (classify when I know where they go)
    ###########################################

    def eint1(gen self, long n=0):
        """
        eint1(x,{n})

        Returns the exponential integral E1(x). If the optional
        argument n is given, computes the vector of the first n
        values of the exponential integral E1(x*n), for x > 0.

        REFERENCE: See page 262, Prop 5.6.12, of Cohen's book
        "A Course in Computational Algebraic Number Theory".
        """
        cdef int err
        err = setjmp(GP_DATA.env)
        if err: _error(err,"error in eint1")

        cdef int t
        t = typ(self.g)
        if t != t_INT and t != t_REAL:
            raise TypeError, "Input (x=%s) must be a real number."%self

        if n <= 0:
            _sig_on
            return new_gen(eint1(forcecopy(self.g), REAL_PREC))
        else:
            if gsigne(self.g) <= 0:
                raise ArithmeticError, "Input (x=%s) must be a positive real number."%self
            _sig_on
            return new_gen(veceint1(forcecopy(self.g), stoi(n), REAL_PREC))


    def order(self):
        cdef int err
        err = setjmp(GP_DATA.env)
        if err: _error(err,"error in order")
        _sig_on
        return new_gen(order(forcecopy(self.g)))

    def znprimroot(self):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in znprimroot")
        _sig_on
        return new_gen(ggener(forcecopy(self.g)))

    def __abs__(self):
        return self.abs()

    def abs(gen self):
        cdef int err
        err = setjmp(GP_DATA.env)
        if err: _error(err,"error in abs")
        _sig_on
        return new_gen(gabs(forcecopy(self.g), REAL_PREC))

    def mattranspose(gen self):
        """
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in mattranspose")
        _sig_on
        return new_gen(gtrans(forcecopy(self.g)))

    def matker(gen self, long flag=0):
        """
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in matker")
        _sig_on
        return new_gen(matker0(forcecopy(self.g), flag))

    def matkerint(gen self, long flag=0):
        """
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in matkerint")
        _sig_on
        return new_gen(matkerint0(forcecopy(self.g), flag))

    def trace(gen self):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in trace")
        _sig_on
        return new_gen(gtrace(forcecopy(self.g)))

    def norm(gen self):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in norm")
        _sig_on
        return new_gen(gnorm(forcecopy(self.g)))

    def charpoly(gen self, var=-1, flag=0):
        """
        charpoly(A,{v=x},{flag=0}): det(v*Id-A) = characteristic polynomial of the
        matrix A using the comatrix. flag is optional and may be set to 1 (use
        Lagrange interpolation) or 2 (use Hessenberg form), 0 being the default.
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in charpoly")
        _sig_on
        return new_gen(charpoly0(forcecopy(self.g), get_var(var), flag))


    def nextprime(gen self):
        """
        nextprime(x): smallest pseudoprime >= x
        """
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in nextprime")
        _sig_on
        return new_gen(gnextprime(forcecopy(self.g)))

    def subst(gen self, var, y):
        """
        EXAMPLES:
           sage: x = pari("x"); y = pari("y")
           sage: f = x**3 + 17*x + 3
           sage: f.subst(x, y)
           y^3 + 17*y + 3
           sage: f.subst(x, "z")
           z^3 + 17*z + 3
           sage: f.subst(x, "z")**2
           z^6 + 34*z^4 + 6*z^3 + 289*z^2 + 102*z + 9
           sage: f.subst(x, "x+1")
           x^3 + 3*x^2 + 20*x + 21
           sage: f.subst(x, "xyz")
           xyz^3 + 17*xyz + 3
           sage: f.subst(x, "xyz")**2
           xyz^6 + 34*xyz^4 + 6*xyz^3 + 289*xyz^2 + 102*xyz + 9
        """
        cdef gen _y
        _y = pari(y)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in subst")
        _sig_on
        return new_gen(gsubst(forcecopy(self.g), get_var(var), _y.g))


    def __call__(gen self, y):
        cdef gen _y
        _y = pari(y)
        cdef int n
        n = setjmp(GP_DATA.env)
        if n: _error(n,"error in __call__")
        _sig_on
        return new_gen(gsubst(forcecopy(self.g), get_var(0), _y.g))

    def kronecker(gen self, y):
        cdef gen _y
        _y = pari(y)

        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in imag")
        _sig_on

        return new_gen(gkronecker(forcecopy(self.g), _y.g))

    def exp(gen self):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in imag")
        _sig_on
        return new_gen(gexp(forcecopy(self.g), REAL_PREC))

    def log(gen self):
        cdef int n
        n = setjmp(GP_DATA.env)
        if n:
            _error(n,"error in imag")
        _sig_on
        return new_gen(glog(forcecopy(self.g), REAL_PREC))


    def type(gen self):
        return str(type_name(typ(self.g)))


###########################################
# Create a gen from a GEN object.
# This *steals* a reference to the GEN.
###########################################
cdef pari_sp stack_avma

cdef save_avma():
    global avma, stack_avma
    stack_avma = avma

cdef gen new_gen_noclear(GEN x):
    """
    Create a new gen, but don't free any memory on the stack.
    """
    global top

    return _new_gen(x,top)

## cdef gen new_gen_from(GEN x, pari_sp av):
##     """
##     Create a new gen, then free everything on the stack below av.
##     The gen *must* not reference anything on the stack above av.
##     """
##     global avma
##     cdef gen g
##     g = _new_gen(x, av)
##     avma = av
##     return g

cdef gen new_gen(GEN x):
    """
    Create a new gen, then free the *entire* stack.
    """
    cdef gen g
    global top, avma
    g = _new_gen(x, top)
    avma = top
    return g

cdef gen __xxx_new_gen(GEN x, pari_sp start):
    cdef size_t s
    cdef pari_sp tmp_bot, tmp_top, tmp_avma
    cdef GEN h
    cdef gen p

    global avma, bot, top

    if <pari_sp>x < bot or <pari_sp>x > top:
        p = gen()
        p.init(x, 0)
        return p

    tmp_top = top
    tmp_bot = bot
    tmp_avma = avma
    s = <size_t> (start - avma)
    #print "Allocating %s bytes for PARI/Python object"%(<long> s)
    bot = <pari_sp> PyMem_Malloc(s)
    top = bot + s
    avma = top
    h = forcecopy(x)
    p = gen(); p.init(h, bot)

    # Restore the stack to how it was before x was created.
    top = tmp_top
    bot = tmp_bot
    avma = tmp_avma
    return p

cdef gen _new_gen(GEN x, pari_sp start):
    cdef size_t s
    cdef pari_sp tmp_bot, tmp_top, tmp_avma
    cdef GEN h
    cdef gen p

    global avma, bot, top

    if <pari_sp>x < bot or <pari_sp>x > top:
        p = gen()
        p.init(x, 0)
        return p

    tmp_top = top
    tmp_bot = bot
    tmp_avma = avma
    h = forcecopy(x)
    s = <size_t> (tmp_avma - avma)
    #print "Allocating %s bytes for PARI/Python object"%(<long> s)
    bot = <pari_sp> PyMem_Malloc(s)
    top = bot + s
    avma = top
    h = forcecopy(x)
    p = gen(); p.init(h, bot)

    # Restore the stack to how it was before x was created.
    top = tmp_top
    bot = tmp_bot
    avma = tmp_avma
    return p

cdef gen new_ref(GEN x, g):
    cdef gen p
    p = gen(); p.init(x, 0)
    p.refers_to[-1] = g  # so underlying memory won't get deleted
                         # out from under us.
    return p

def pari(s):
    """
    Create the PARI object got by evaluating s using PARI.
    """
    cdef gen p
    cdef GEN g
    cdef pari_sp prev_avma
    global avma

    prev_avma = avma

    if isinstance(s, gen):
        return s
    elif isinstance(s, (list, xrange, tuple)):
        v = vector(len(s))
        for i, x in enumerate(s):
            v[i] = pari(x)
        return v
    elif isinstance(s, bool):
        if s:
            return ONE
        return ZERO
    else:
        try:
            return s.pari()
        except AttributeError:
            pass

    t = str(s)
    again = True
    while again:
        try:
            _sig_on
            g = str_to_GEN(t)
            _sig_off
            again = False
        except MemoryError:
            again = True

    #return new_gen_from(g, prev_avma)
    return new_gen(g)


cdef int get_var(v):
    """
    Converts a Python string into a PARI variable reference number.
    Or if v = -1, returns -1.
    """
    if v != -1:
        s = str(v)
        return fetch_user_var(s)
    return -1

def _stack_ptr():
    return avma

def _stack_top():
    return top

def _stack_bottom():
    return bot

def _stack_size():
    """
    Stack size, in bytes
    """
    #return int((top-bot)/1048576)
    return top-bot

def stack_usage():
    """
    Return how much memory in megabytes is currently being used by the
    PARI stack.
    """
    #return (top-avma)/1048576
    return top-avma

def stack_available():
    """
    Return how much memory in megabytes is currently available on the
    PARI stack.
    """
    #return (avma-bot)/1048576
    return avma-bot


######################################################
# Direct interface to PARI C-level interpreter  (avoid using this)
######################################################

def execute(s):
    """
    Excecute s using the PARI C library.
    """
    cdef int n
    n = setjmp(GP_DATA.env)
    if n:
        _error(n, s)
        _sig_on
        execute(s)   # if call returns, then was memory error so try again
        _sig_off

    # Parse s in PARI.
    s = s + "\n"
    _sig_on
    flisseq(s)
    _sig_off

def eval(s):
    """
    Parse s in PARI and evaluate the result, then return it as a string.
    """
    cdef int n
    n = setjmp(GP_DATA.env)
    if n:
        _error(n, s)
        return eval(s)   # if call returns, then was memory error so try again

    cdef GEN x
    cdef char *c

    # Parse s in pari and evaluate result.
    s = s + "\n"
    _sig_on
    x = flisseq(s)
    # Get result as a C-string
    c = GENtostr(x)
    _sig_off
    # De-callocate memory used by x:
    cgiv(x)
    # Convert c to a Python string.
    s = str(c)
    # Free memory used by c
    free(c)
    return s

def _error(n, s=""):
    _sig_off
    import sys
##     if n == 20:
##         raise TypeError, "Incorrect type (%s)"%s
##     if n == 23:
##         raise TypeError, "Unable to assign I-->S.  (%s)"%s
##     if n == 24:
##         raise TypeError, "Unable to convert to integer: %s"%s
##     if n == 65:
##         raise OverflowError, "Object too big to fit in a codeword: %s"%s
##     if n == 83: # not enough precomputed primes... (might be useful)
##         global num_primes
##         print "Doubling number of precomputed primes."
##         init_primes(2*num_primes)
##         return
##     if n == 97:
##         raise InputError, "Bad arguments to an elliptic curve function: %s"%s
    if n == MEMERR: # out of memory
        print "PARI ran out of memory, so we automatically increase the memory available to PARI."
        init_stack(0)  # double current stack size; all data on stack is erased.
        return
    raise RuntimeError, "PARI error %s: %s"%(n,s)



def min(x,y):
    """
    min(x,y): Return the minimum of x and y.
    """
    if x <= y:
        return x
    return y

def max(x,y):
    """
    max(x,y): Return the maximum of x and y.
    """
    if x >= y:
        return x
    return y

def prime_list(long n):
    """
    prime_list(n): returns list of the first n primes

    INPUT:
        n -- C long
    OUTPUT:
        gen -- PARI list of first n primes
    EXAMPLES:
        sage: prime_list(0)
        []
        sage: prime_list(-1)
        []
        sage: prime_list(3)
        [2, 3, 5]
        sage: prime_list(10)
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        sage: prime_list(20)
        [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
        sage: len(prime_list(1000))
        1000
    """

    cdef int e
    e = setjmp(GP_DATA.env)
    if e: _error(e,"prime_list")
    _sig_on
    return new_gen(primes(n))

def nth_prime(long n):
    """
    nth_prime(n): returns the n-th prime, where n is a C-int
    """
    global num_primes

    if n <= 0:
        raise ArithmeticError, "nth prime meaningless for negative n (=%s)"%n
    cdef GEN g
    z = setjmp(GP_DATA.env)
    if z:
        init_primes(max(2*num_primes,20*n))
        return nth_prime(n)
        #_error(z, "Not enough precomputed primes")
    _sig_on
    g = prime(n)
    return new_gen(g)

def euler(int prec=0):
    """
    Return Euler's constant.

    EXAMPLES:
        sage.: euler()
        0.57721566490153286060651209008
        sage.: get_real_precision()
        5
        sage.: euler(6)
        0.57721566490153286060651209008240243104
        sage.: set_real_precision(10)
        sage.: euler()
        0.57721566490153286060651209008240243104215933593992359880576723488486772677767

    We restore precision to the default.
        sage: set_real_precision(5)

    Note that euler is still saved to the largest precision to which it has been
    computed so far:
        sage.: euler()
        0.57721566490153286060651209008240243104215933593992359880576723488486772677767
    """
    if prec<=3:
        consteuler(REAL_PREC)
    else:
        consteuler(prec)
    cdef gen g
    g = gen()
    g.init(geuler,0)
    return g

def pi(int prec=0):
    """
    Return the value of the constant pi = 3.1415....

    EXAMPLES:
        sage.: set_real_precision(5)
        sage.: pi()
        3.1415926535897932384626433833
        sage.: pi(6)
        3.14159265358979323846264338327950288419
        sage.: get_real_precision()
        5
        sage.: set_real_precision(10)
        sage.: pi()
        3.1415926535897932384626433832795028841971693993751058209749445923078164062862

    We restore precision to the default.
        sage.: set_real_precision(5)

    Note that pi is still saved to the largest precision to which it has been
    computed so far:
        sage.: pi()
        3.1415926535897932384626433832795028841971693993751058209749445923078164062862
    """
    if prec <= 3:
        constpi(REAL_PREC)
    else:
        constpi(prec)

    cdef gen g
    g = gen()
    g.init(gpi,0)
    return g

def factorial(long n):
    """
    Return the factorial of the integer n as a PARI gen.
    """
    _sig_on
    return new_gen(mpfact(n))

def vector(long n, entries=None):
    """
    vector(long n, entries=None):
    Create and return the length n PARI vector with given list of entries.
    """
    cdef int e
    e = setjmp(GP_DATA.env)
    if e: _error(e,"vector")
    cdef gen v
    _sig_on
    v = new_gen(zerovec(n))
    if entries != None:
        if len(entries) != n:
            raise IndexError, "length of entries (=%s) must equal n (=%s)"%\
                  (len(entries), n)
        for i, x in enumerate(entries):
            v[i] = x
    return v

def matrix(long m, long n, entries=None):
    """
    matrix(long m, long n, entries=None):
    Create and return the m x n PARI matrix with given list of entries.
    """
    cdef int e, i, j, k
    e = setjmp(GP_DATA.env)
    if e: _error(e,"matrix")
    cdef gen A
    _sig_on
    A = new_gen(forcecopy(zeromat(m,n)))
    if entries != None:
        if len(entries) != m*n:
            raise IndexError, "len of entries (=%s) must be %s*%s=%s"%(len(entries),m,n,m*n)
        k = 0
        for i from 0 <= i < m:
            for j from 0 <= j < n:
                A[i,j] = entries[k]
                k = k + 1
    return A

def finitefield_init(p, long n, var=-1):
    """

    finitefield_init(p, long n, var="x"): Return a polynomial f(x) so
    that the extension of F_p of degree n is k = F_p[x]/(f).  Moreover,
    the element x (mod f) of k is a generator for the multiplicative
    group of k.

    INPUT:
        p -- int, a prime number
        n -- int, positive integer
        var -- str, default to "x", but could be any pari variable.
    OUTPUT:
        pari polynomial mod p -- defines field
    EXAMPLES:
        sage: finitefield_init(97,1)
        Mod(1, 97)*x + Mod(92, 97)

    The last entry in each of the following two lines is
    determined by a random algorithm.
        sage.: finitefield_init(7,2)
        Mod(1, 7)*x^2 + Mod(6, 7)*x + Mod(3, 7)
        sage.: finitefield_init(2,3)
        Mod(1, 2)*x^3 + Mod(1, 2)*x^2 + Mod(1, 2)
    """
    cdef gen _p, _f2, s
    cdef int err
    cdef long x
    cdef GEN v, g
    _p = pari(int(p))
    err = setjmp(GP_DATA.env)
    if err: _error(err,"error in finitefield_init")
    if n < 1:
        raise ArithmeticError, "Degree n (=%s) must be at least 1."%n
    if _p < 2 or not _p.isprime():
        raise ArithmeticError, "Prime p (=%s) must be prime."%_p
    x = get_var(var)
    if n == 1:
        return new_gen(ffinit(_p.g, n, x)) - _p.znprimroot()

    _sig_on
    f = new_gen(ffinit(_p.g, n, x))
    _f2 = f.lift()
    g = FpXQ_gener(_f2.g, _p.g)
    s = new_gen(g)*ONE.Mod(p)
    return s.Mod(f).charpoly(var)


#################################################################################
## Note:  The function finitefield_init above makes use of FpXQ_gener, which
## is not defined in any PARI headers or available from gp.  I learned about
## it from a post by Karim Belabas from 2002 at
##
##    http://pari.math.u-bordeaux.fr/archives/pari-users-0205/msg00013.html
##
## "There is something built-in, very carefully hidden [ used by idealstar() ].
## provided you have an up-to-date development version from the CVS server, and
## provided install() works on your system, you can use the following:
##
##   install(FpXQ_gener,GG)
##
##   ffinitprim(p,n) =
##   { local(ffp);
##
##     ffp = ffinit(p,n,x);
##     Mod(FpXQ_gener(lift(ffp), p) * Mod(1,p), ffp)
##   }
##
## (14:07) gp > ranffinitprim(101, 40);
## time = 1mn, 40,980 ms.
## (14:08) gp > ranffinitprim(2, 100)
##  ***   user interrupt after 13mn, 57,260 ms.
##
##
## (14:28) gp > ffinitprim(101, 40);
## time = 940 ms.
## (14:28) gp > ffinitprim(2, 100)
## time = 330 ms.
##
## Obviously, this is quite a useful routine, so I'll have to make it directly
## available to gp someday (with a decent name).
##
## Cheers,
##
##     Karim."
#################################################################################


def set_real_precision(long n):
    """
    Set the default precision to n words (word = 4 bytes on
    a 32-bit machine).

    EXAMPLES:
        sage: get_real_precision()
        5
        sage: set_real_precision(10)
        sage: get_real_precision()
        10
        sage: set_real_precision(5)
    """
    global REAL_PREC
    if n <= 2: n = 3
    REAL_PREC = n

def get_real_precision():
    """
    Return the default PARI real precision.  This is the precision
    to which many objects are computed by default.

    EXAMPLES:
        sage: get_real_precision()
        5
        sage: set_real_precision(10)
        sage: get_real_precision()
        10
        sage: set_real_precision(5)
    """
    global REAL_PREC
    return REAL_PREC

def set_series_precision(long n):
    """
    Set the default precision used when coercing objects into power series.

    EXAMPLES:
        sage: default = get_series_precision()
        sage: f = pari('1+x + x^2 + x^3')
        sage: set_series_precision(3)
        sage: f.exp()
        2.718281828459045235360287471 + 2.718281828459045235360287471*x + 4.077422742688567853040431207*x^2 + 5.889610628327931343280622855*x^3 + O(x^4)
        sage: set_series_precision(10)
        sage: f.exp()
        2.718281828459045235360287471 + 2.718281828459045235360287471*x + 4.077422742688567853040431207*x^2 + 5.889610628327931343280622855*x^3 + 5.549825399770550688860586921*x^4 + 5.912262976898423386908625250*x^5 + 5.780124276903886465745277943*x^6 + 4.893446632858912186592041315*x^7 + 4.273810514670244409851059119*x^8 + 3.489008512344414242252330620*x^9 + O(x^10)
        sage: set_series_precision(default)
    """
    global precdl
    if n <= 2: n = 3
    precdl = n

def get_series_precision():
    """
    Get the default precision used when coercing objects into power series.

    EXAMPLES:
        sage: get_series_precision()
        16
    """
    global precdl
    return int(precdl)

def series_precision(new_val=None):
    """
    Return the default series precision, which is the default precision
    to which series are computed.
    """
    if new_val != None:
        set_series_precision(new_val)
    return get_series_precision()

def init_primes(unsigned long M):
    """
    Recompute the primes table including at least all primes up to M.

    EXAMPLES:
        sage: init_primes(200000)
    """
    global diffptr, num_primes
    free(<void*> diffptr)
    num_primes = M
    diffptr = initprimes(M)

def __read_script(char* s):
    cdef int err
    cdef pari_sp av
    err = setjmp(GP_DATA.env)
    if err:
        _error(err,"error in import_pari_script")
    _sig_on
    gp_read_str(s)
    #flisseq(s)
    _sig_off
    global top, avma
    top = avma

    # new gp_read_str, gp_read_file


def read(filename):
    r"""
    Read a script from the named filename into the interpreter, where
    s is a string.  The functions defined in the script are then
    available for use from SAGE/PARI.

    EXAMPLE:

        If foo.gp is a script that contains
        \begin{verbatim}
            {foo(n) =
                n^2
            }
        \end{verbatim}
        and you type \code{read("foo.gp")}, then the command
        \code{pari("foo(12)")} will create the Python/PARI gen which
        is the integer 144.

    CONSTRAINTS:
        The PARI script must *not* contain the following function calls:

             print, default, ???    (please report any others that cause trouble)

        Also multiline functions should be written in the following form:
        \begin{verbatim}
            {foo(x) =
                code...
            }
        \end{verbatim}

        and *NOT* in the form
        \begin{verbatim}
            foo(x) =
            {
                    code...
            }
        \end{verbatim}
    """
    F = open(filename).read()
    __read_script(F)
    return
    while True:
        i = F.find("{")
        if i == -1:
            __read_script(F)
            break
        __read_script(F[:i])
        j = F[i:].find("}") + i
        __read_script(F[i:j+1])
        F = F[j+1:]


def pari_real_precision(n):
    """
    Sets the PARI default real precision, both for creation of
    new objects and for printing.

    TODO: For example, log(2) is computed differently, depending on the
    pari_real_precision.
    """
    n = int(n)
    s = str(n)
    err = setjmp(GP_DATA.env)
    if err:
        raise ValueError
    global REAL_PREC, prec
    # REAL_PREC = n
    sd_realprecision(s, 2)
    REAL_PREC = prec
    # TODO : series precision should match precdl also.

# Some useful globals (initialized in init())
cdef gen ZERO, ONE, TWO


def allocate_mem():
    print "Doubling the PARI stack."
    init_stack(0)

# Initialize the PARI system.  There's a bit of a hack here.  I've
# found by experiment that the global PARI C-library variable precdl
# (which is the default series precision) is 0 if and only if the PARI
# system has not been initialized.  A PARI expert could probably
# replace this test by something more sensible.

if bot == 0:
    _init()
