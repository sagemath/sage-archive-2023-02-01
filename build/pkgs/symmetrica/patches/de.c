/* de.c SYMMETRICA */
#include "def.h"
#include "macro.h"


#include <sys/types.h>


#ifdef unix
#undef MSDOS
#include <sys/times.h>
#endif /* unix */

#include <time.h>   /* for the routine clock,time */


#ifdef unix
#include <sys/param.h>
#endif /* unix */

OP cons_drei;   /* global INTEGER variable 3 */
OP cons_zwei;   /* global INTEGER variable 2 */
OP cons_eins;   /* global INTEGER variable 1 */
OP cons_negeins;/* global INTEGER variable -1 */
OP cons_null;   /* global INTEGER variable 0 */
FILE *texout;   /* global variable for texoutput */
INT no_banner = TRUE; /* AK 281293 */
INT no_mem_check=TRUE; /* AK 100893 */
INT english_tableau=FALSE; /* AK 290995 */

INT doffset=0L;  /* global for debugprint AK 160393 */

INT freeall_speichersize_max = (INT) 1000000;
int SYM_free(a) char *a;
{
    if (sym_timelimit > 0L) check_time();
    free(a);
    return 0;
}

char * SYM_malloc(a) int a;
{
    INT erg = OK;
    char *res;
    INT err;
    if (sym_timelimit > 0L) check_time();
    SYMCHECK( (a < 0) , "SYM_malloc: size < 0");
sca:
    res =  (char*)malloc(a);
    if (res == NULL)
        {
        err=error("SYM_malloc: no memory");
        if (err==ERROR_RETRY) goto sca;
        if (err==ERROR_EXPLAIN) {
            fprintf(stderr,"I wanted %d Byte of Memory", a); }
        }
    return res;
    ENDTYP("SYM_malloc",char *);
}

char * SYM_calloc(a,b) int a,b;
{
    char *erg;
    INT err;
    if (sym_timelimit > 0L) check_time();

    if ( a < 0 )
        {
        err = error("SYM_calloc: negative number of entries");
        if (err==ERROR_EXPLAIN) {
            fprintf(stderr,"I wanted %d pieces of size %d", a,b); }
        return NULL;
        }
    else if ( b < 0 )
        {
        err = error("SYM_calloc: negative size");
        if (err==ERROR_EXPLAIN) {
            fprintf(stderr,"I wanted %d pieces of size %d", a,b); }
        return NULL;
        }
sca:
    erg=(char*) calloc(a,b);
    if (erg == NULL)
        {
        err=error("SYM_calloc: no memory");
        if (err==ERROR_RETRY)
            {
            goto sca;
            }
        if (err==ERROR_EXPLAIN) {
            fprintf(stderr,"I wanted %d pieces of size %d", a,b);
            goto sca;
            }
        }
    return erg;
}

char * SYM_realloc(a,b) char *a; int b;
{
    char *erg;
    INT err= -1;
    if (sym_timelimit > 0L) check_time();
sca:
    erg = (char *)realloc(a,b);

    if (erg == NULL)
        {
        err=error("SYM_realloc: no memory");
        if (err == ERROR_RETRY)
            {
            goto sca;
            }
        if (err==ERROR_EXPLAIN)
            {
            fprintf(stderr,"I wanted %d Byte of Memory", b);
            goto sca;
            }
        }
    return erg;
}



INT anfang()
/* AK 070890 V1.1 */ /* AK 210891 V1.3 */
/* AK 260298 V2.0 */
/* AK 280705 V3.0 */
{
    time_t l;
    INT erg = OK;
    void srand();
    if (not no_banner)
        {
        printeingabe("SYMMETRICA VERSION 3.0 - STARTING");
        printeingabe(TITELTEXT);
        }

    time(&l);
    l = l * l * clock();
#ifdef unix
    l = l * getpid();
#endif
    srand((unsigned long)l);
    memcheck("anfang");
    fflush(stdout); fflush(stderr);

    erg += speicher_anfang();
    NEW_INTEGER(cons_drei,3);
    NEW_INTEGER(cons_zwei,2);
    NEW_INTEGER(cons_eins,1);
    NEW_INTEGER(cons_negeins,-1);
    NEW_INTEGER(cons_null,0); /* needed in start_longint */

    texmath_yn=0L; /* not in math mode */

#ifdef LONGINTTRUE
    start_longint();
#endif /* LONGINTTRUE */

    check_time_co = NULL; /* co routine called in check time,
                may be set by other programms */
    texout = stdout;

#ifdef NUMBERTRUE    /* 291091: TPMcD */
/* The third parameter is NULL or the name of a file with cyclotomic data */
    setup_numbers(STD_BASIS,TRUE, NULL);

#endif /* NUMBERTRUE */

#ifdef BRUCHTRUE
    bruch_anfang(); /* AK 100893 */
#endif /* BRUCHTRUE */

#ifdef VECTORTRUE
    vec_anfang(); /* AK 100893 */
#endif /* VECTORTRUE */
#ifdef PARTTRUE
    part_anfang(); /* AK 040903 */
#endif /* PARTTRUE */

#ifdef TABLEAUXTRUE
    tab_anfang(); /* AK 100893 */
#endif /* TABLEAUXTRUE */

#ifdef PERMTRUE
    perm_anfang(); /* AK 100893 */
#endif /* PERMTRUE */

#ifdef LISTTRUE
    list_anfang(); /* AK 100893 */
#endif /* LISTTRUE */

#ifdef POLYTRUE
    monom_anfang(); /* AK 100893 */
#endif /* POLYTRUE */
#ifdef FFTRUE
    ff_anfang(); /* AK 011204 */
#endif /* FFTRUE */
#ifdef GRTRUE
    galois_anfang(); /* AK 271106  */
#endif /* GRTRUE */
#ifdef LOCALTRUE
    local_anfang(); /* AK 280705 */
#endif

    /* checks on type of constants */
    CTO(INTEGER,"anfang(e1)",cons_zwei);
    CTO(INTEGER,"anfang(e2)",cons_eins);
    CTO(INTEGER,"anfang(e3)",cons_negeins);
    CTO(INTEGER,"anfang(e4)",cons_null);
    CTO(INTEGER,"anfang(e5)",cons_drei);
    ENDR("anfang");
}


INT ende()
/* AK 070890 V1.1 */ /* AK 210891 V1.3 */
{
    INT erg = OK;
    char t[100];


#ifdef SCHURTRUE
    schur_ende();
#endif /* SCHURTRUE */

#ifdef NUMBERTRUE    /* 29.10.91: TPMcD */
    release_numbers();
#endif /* NUMBERTRUE */

#ifdef POLYTRUE
    monom_release();
#endif /* POLYTRUE */


#ifdef TABLEAUXTRUE
    tab_ende(); /* AK 100893 */
#endif /* TABLEAUXTRUE */

    hash_ende();

#ifdef POLYTRUE
    monom_ende(); /* AK 100893 */ /* nach schur ende */
#endif /* POLYTRUE */

#ifdef BRUCHTRUE
    bruch_ende(); /* AK 100893 */
#endif /* BRUCHTRUE */

#ifdef PARTTRUE
    part_ende();
#endif /* PARTTRUE */

#ifdef LISTTRUE
    list_ende(); /* AK 100893 */
#endif /* LISTTRUE */

#ifdef PERMTRUE
    perm_ende(); /* AK 100893 */
#endif /* PERMTRUE */

#ifdef FFTRUE
    ff_ende();
#endif /* FFTRUE */
#ifdef GRTRUE
    galois_ende(); /* AK 271106  */
#endif /* GRTRUE */

#ifdef LOCALTRUE
    local_ende(); /* AK 280705 */
#endif
#ifdef NUMBERTRUE    /* AK 310893 */
    nb_ende();
#endif /* NUMBERTRUE */

#ifdef LONGINTTRUE
    longint_ende();
#endif /* LONGINTTRUE */

#ifdef VECTORTRUE
    vec_ende(); /* AK 100893 */
#endif /* VECTORTRUE */


    if  (  /* AK 190194 */
        (S_O_K(cons_drei) != INTEGER) ||
        (S_O_K(cons_null) != INTEGER) ||
        (S_O_K(cons_zwei) != INTEGER) ||
        (S_O_K(cons_eins) != INTEGER) ||
        (S_O_K(cons_negeins) != INTEGER) ||
        (S_I_I(cons_null) != (INT) 0) ||
        (S_I_I(cons_zwei) != (INT) 2) ||
        (S_I_I(cons_eins) != (INT) 1) ||
        (S_I_I(cons_negeins) != (INT) -1)
        )
    erg += error("ende: wrong constant values e.g. cons_null");
    erg += freeall(cons_null);
    erg += freeall(cons_zwei);
    erg += freeall(cons_drei);
    erg += freeall(cons_eins);
    erg += freeall(cons_negeins);

    erg += speicher_ende();


    memcheck("ende");

    if (not no_banner)
        {
        printeingabe("\nSYMMETRICA VERSION 3.0 - ENDING");
        sprintf(t,"last changed: %s",TITELTEXT); /* AK 181194 */
        printeingabe(t);
        }

    fflush(stdout);
    fflush(stderr);
    return erg;
}
INT runtime(l) long *l;
/* AK 270689 V1.0 */ /* AK 070890 V1.1 */ /* AK 210891 V1.3 */
{
#ifdef UNDEF
#ifdef unix
    struct tms buffer;
    times(&buffer);
    *l = (long) buffer.tms_utime;
#else /* clock ist POSIX */
    *l = (long) clock()/60;
#endif /* unix */
#endif
    *l = (long) clock()/CLOCKS_PER_SEC;
    return OK;
}

INT get_time(a) OP a;
/* AK 160890 V1.1 */ /* AK 210891 V1.3 */
/* AK 300998 V2.0 */
{
    long l;
    runtime(&l);
    return m_i_i((INT)l,a);
}


INT print_time()
/* AK 160890 V1.1 */ /* AK 210891 V1.3 */
{
    long l;
    runtime(&l);
    printf("zeit:%ld\n",l);return OK;
}


INT fusedmemory(fn,stelle) FILE *fn; char *stelle;
/* AK 270689 V1.0 */ /* AK 010290 V1.1 */ /* AK 130691 V1.2 */
/* AK 210891 V1.3 */
{
#ifdef unix
#ifndef linux
/*
    struct mallinfo mallinfo();
    struct mallinfo ergebnis;
    free(calloc(1,1));
    ergebnis = mallinfo();
    fprintf(fn,"%s: ",stelle);
    fprintf(fn,"%d ",ergebnis.uordblks);
    fprintf(fn,"%d\n",ergebnis.usmblks);
    return(OK);
*/
#endif /* linux */
#endif /* unix */
#ifdef TURBOC
/*
    fprintf(fn,"%s: ",stelle);
    fprintf(fn,"%ul\n",coreleft());
    return(OK);
*/
#endif /* TURBOC */
    return(OK);
}

INT mem_small()
/* anzahl small memory zurueck */
/* AK 270689 V1.0 */ /* AK 070890 V1.1 */ /* AK 210891 V1.3 */
{
#ifdef unix
#ifndef linux
/*
    struct mallinfo mallinfo();
    struct mallinfo ergebnis;
    ergebnis = mallinfo();
    return(ergebnis.usmblks);
*/
#endif /*linux */
#endif /* unix */
    return(0);
}


INT memcheck(stelle) char *stelle;
/* informationen ueber memory 31/10/86 */
/* AK 270689 V1.0 */ /* AK 010290 V1.1 */ /* AK 210891 V1.3 */
{
#ifdef unix
#ifndef linux
/*
    struct mallinfo mallinfo();
    struct mallinfo ergebnis;

    if (no_mem_check == TRUE) return OK;
    SYM_free(SYM_calloc(1,1));
    ergebnis = mallinfo();
    printf("memory information  %s\n",stelle);
    printf("total space     %d\n",ergebnis.arena);
    printf("block number    %d\n",ergebnis.ordblks);
    printf("small blocks    %d\n",ergebnis.smblks);
    printf("used blocks     %d\n",ergebnis.uordblks);
    printf("free blocks     %d\n",ergebnis.fordblks);
    printf("used sm. blocks %d\n",ergebnis.usmblks);
    printf("free sm. blocks %d\n",ergebnis.fsmblks);
    return(OK);
*/
#endif /*linux */
#endif /* unix */
    return(OK);
}

INT sym_background = 0L;
INT sym_www = 0L;
INT sym_timelimit = 0L;
INT fatal_error(fehlertext) char *fehlertext;
/* AK 270295 */
{
    fprintf(stderr,"fatal error in function %s\n",fehlertext);
    exit(11);
    return OK;
}
INT error(fehlertext) char *fehlertext;
/* if answer == a ==> abort
   if answer == e ==> explain
   if answer == g ==> go on
   if answer == r ==> retry
   if answer == s ==> go on supress error texts
   if answer == f ==> go on forever
   else               exit */
/* AK 270689 V1.0 */ /* AK 070890 V1.1 */
/* AK 070291 V1.2 explanation of possible input */
/* AK 210891 V1.3 */
{
    char antwort[2];
    static int forever=0;
    if (forever==2) return ERROR;
    if (sym_www) {
        printf("ERROR: %s?: ",fehlertext);
        exit(ERROR_BACKGROUND);
        }
    fflush(stdout);
    fflush(stderr);
    fprintf(stderr,
"\nenter a to abort with core dump, g to go, f to supress\n");
    fprintf(stderr,
"s to supress further error text, r to retry,  e to explain, else stop\n");
    fprintf(stderr,"ERROR: %s?: ",fehlertext);


    fflush(stderr);

    if (sym_background) {
        fprintf(stderr,"\nerror occured in background mode finishing SYMMETRICA\n");
        exit(ERROR_BACKGROUND);
        }

    if (forever==1) return ERROR;

    antwort[0]='X';
    scanf("%s",antwort);
    if (antwort[0] == 'a') abort();
    if (antwort[0] == 'f') {forever = 1; return ERROR;}
    if (antwort[0] == 's') {forever = 2; return ERROR;}
    if (antwort[0] == 'g') return ERROR;
    if (antwort[0] == 'r') return ERROR_RETRY;
    if (antwort[0] == 'e') return ERROR_EXPLAIN;
    exit(1); /* AK 121192 */
}


INT no_memory()
/* AK 090792 */
{
    return error("no memory left");
}

INT debugprint(a) OP a;
/* AK 260788 */ /* AK 030789 V1.0 */ /* AK 130690 V1.1 */ /* AK 210891 V1.3 */
{
    OBJECTKIND kind;
    INT i,j,k;
    char *text=NULL;
    for (i=0L;i<doffset;i++) fprintf(stderr," ");
    if (a==NULL) {
        fprintf(stderr,"NULL\n");
        return(OK);
    }
    kind = s_o_k(a);
    switch ((int)kind)
    /* abschluss immer mit newline */
    {
    case 0:
        fprintf(stderr,"kind:0=empty self=%ld\n",s_o_s(a).ob_INT);
        break;
    case 1:
        fprintf(stderr,"kind:1=integer value:");
        fprintf(stderr,"%ld\n",s_i_i(a));
        return(OK);
#ifdef VECTORTRUE
        case 120199: case 31:
        case 26: case 19: case 15:
        case 2:
        case 211106:
        if (kind == 2) text="vector";
        if (kind == 15) text="integervector";
        if (kind == 19) text="word";
        if (kind == 26) text="comp";
        if (kind == 31) text="kranz";
        if (kind == 47) text="subset";
        if (kind == 120199) text="hashtable";
        if (kind == 211106) text="galois ring";
        fprintf(stderr,"kind:%d=%s length:\n",(int)kind, text);
        doffset += 2L;
        debugprint(s_v_l(a));
        doffset -= 2L;
        for (i=0L;i<s_v_li(a);i++)
        {
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        fprintf(stderr,"%s %ld-komponente:\n",text,i);
        doffset += 2L;
        debugprint(s_v_i(a,i));
        doffset -= 2L;
        }
        return(OK);
#endif /* VECTORTRUE */
#ifdef PARTTRUE
    case 3:
    case 12:
        {
        if (kind == 12) text="augpartition";
        if (kind == 3) text="partition";
        fprintf(stderr,"kind:%d=%s kind:%d hash:%d\n",(int)kind,text,
                    (int)s_pa_k(a),
                                        (int)s_pa_hash(a));
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        fprintf(stderr,"%s self:\n",text);
        doffset += 2L;
        debugprint(s_pa_s(a));
        doffset -= 2L;
        return(OK);
        }
#endif /* PARTTRUE */
#ifdef BRUCHTRUE
    case 4:
        {
        fprintf(stderr,"kind:4=bruch gekuerzt=%ld oben:\n", s_b_i(a));
        doffset += 2L;
        debugprint(s_b_o(a));
        doffset -= 2L;
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        fprintf(stderr,"bruch unten:\n");
        doffset += 2L;
        debugprint(s_b_u(a));
        doffset -= 2L;
        return(OK);
        }
#endif /* BRUCHTRUE */
#ifdef PERMTRUE
    case 6:
        {
        fprintf(stderr,"kind:6=permutation kind:%d\n",(int)s_p_k(a));
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        fprintf(stderr,"permutation self:\n");
        doffset += 2L;
        debugprint(s_p_s(a));
        doffset -= 2L;
        return(OK);
        }
#endif /* PERMTRUE */
#ifdef SKEWPARTTRUE
    case 7:
        {
        fprintf(stderr,"kind:7=skewpartition gross:\n");
        doffset += 2L;
        debugprint(s_spa_g(a));
        doffset -= 2L;
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        fprintf(stderr,"skewpartition klein:\n");
        doffset += 2L;
        debugprint(s_spa_k(a));
        doffset -= 2L;
        return(OK);
        }
#endif /* SKEWPARTTRUE */
#ifdef TABLEAUXTRUE
    case 8:
        {
        fprintf(stderr,"kind:8=tableaux self:\n");
        doffset += 2L;
        debugprint(s_t_s(a));
        doffset -= 2L;
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        fprintf(stderr,"tableaux umriss:\n");
        doffset += 2L;
        debugprint(s_t_u(a));
        doffset -= 2L;
        return(OK);
        }
#endif /* TABLEAUXTRUE */
#ifdef LISTTRUE
    case 13: case 10: case 29: case 28: case 33: case 32:case 14:
    case 20: case 9: case 42:
        {
        if (kind == 9)    text="polynom";
        if (kind == 20)    text="list";
        if (kind == 14)    text="schubert";
        if (kind == 10) text="schur";
        if (kind == 13) text="homsym";
        if (kind == 28) text="powsym";
        if (kind == 29) text="monomial";
        if (kind == 32) text="groupalgebra";
        if (kind == 33) text="elmsym";
        if (kind == 42) text="monopoly";
        fprintf(stderr,"kind:%d=%s self:\n",(int)kind,text);
        doffset += 2L;
        debugprint(s_l_s(a));
        doffset -= 2L;
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        fprintf(stderr,"%s next:\n",text);
        doffset += 2L;
        debugprint(s_l_n(a));
        doffset -= 2L;
        return(OK);
        }
#endif /* LISTTRUE */
#ifdef MATRIXTRUE
    case 27: case 11: case 40:
        {
        if (kind==11) text = "matrix";
        if (kind==27) text = "kranztypus";
        if (kind==40) text = "integermatrix";

        fprintf(stderr,"kind:%d=%s height:\n",(int)kind,text);
        doffset += 2L;
        debugprint(s_m_h(a));
        doffset -= 2L;
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        fprintf(stderr,"%s length:\n",text);
        doffset += 2L;
        debugprint(s_m_l(a));
        doffset -= 2L;
        fprintf(stderr,"%s hash:%d\n",text,s_m_hash(a));
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        for (i=0L;i<s_m_hi(a);i++)
        for (j=0L;j<s_m_li(a);j++)
        {
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        fprintf(stderr,"%s %ld %ld-komponente:\n",text,i,j);
        doffset += 2L;
        debugprint(s_m_ij(a,i,j));
        doffset -= 2L;
        }
        return(OK);
        }
#endif /* MATRIXTRUE */
#ifdef MONOMTRUE
    case 21:
        {
        fprintf(stderr,"kind:21=monom koeff:\n");
        doffset += 2L;
        debugprint(s_mo_k(a));
        doffset -= 2L;
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        fprintf(stderr,"monom self:\n");
        doffset += 2L;
        debugprint(s_mo_s(a));
        doffset -= 2L;
        return(OK);
        }
#endif /* MONOMTRUE */
#ifdef CHARTRUE
    case 18:
        {
        fprintf(stderr,"kind:18=symchar dim:\n");
        doffset += 2L;
                debugprint(s_sc_d(a));
        doffset -= 2L;
        fprintf(stderr,"symchar partitionen:\n");
        doffset += 2L;
                debugprint(s_sc_p(a));
        doffset -= 2L;
        fprintf(stderr,"symchar werte:\n");
        doffset += 2L;
                debugprint(s_sc_w(a));
        doffset -= 2L;
        return OK;
        }
#endif /* CHARTRUE */
#ifdef LONGINTTRUE
    case 22: return(debugprint_longint(a));
#endif /* LONGINTTRUE */
#ifdef NUMBERTRUE
    case 41: case 43:
        {
        if (kind == 41) text = "cyclotomic";
        if (kind == 43) text = "squareradical";

        fprintf(stderr,"kind:%d=%s self:\n",(int)kind,text);
        doffset += 2L;
        debugprint(s_n_s(a));
        doffset -= 2L;
        return(OK);
        }
#endif /*NUMBERTRUE*/
#ifdef VECTORTRUE
    case 44:
        {
        if (kind == 44) text = "bitvector";
        fprintf(stderr,"kind:%d=%s self:\n",(int)kind,text);
        doffset += 2L;
        C_O_K(a,VECTOR);
        for (k=0L;k<doffset;k++) fprintf(stderr," ");
        fprintf(stderr,"length = number of bits = %ld\n",s_v_li(a));
        C_O_K(a,BITVECTOR);
        doffset -= 2L;
                return(OK);
        }
#endif /*VECTORTRUE */
#ifdef FFTRUE
    case 35: return debugprint_ff(a);
#endif /* FFTRUE */
#ifdef REIHETRUE
    case 36: return debugprint_reihe(a);
#endif /* REIHETRUE */
    default:
        fprintf(stderr,"kind:%ld unknown\n",s_o_k(a));
        break;
    }
    return OK;
}

int SYM_isdigit(a) char a; /* AK 040194 */
{ return ((a >= '0') && (a <= '9')); }

int SYM_strlen(a) char *a; /* AK 030294 */
{ int i=0; while (*a++) i++; return i; }

int SYM_memcmp(a,b,c) char *a,*b; /* AK 210294 */
{ return memcmp(a,b,c); }

int SYM_abs(a) INT a; /* AK 230695 */
{ return (a>0 ) ? a : -a; }


INT mem_size(a) OP a;
/* AK 150295 */
{
    INT erg = OK;
    if (a == NULL)
                return 0;
    switch(S_O_K(a))
        {
        case EMPTY:
        case INTEGER:    return sizeof(struct object);
        case MATRIX:
        case INTEGERMATRIX:
        case KOSTKA:    return mem_size_matrix(a);
        case LONGINT:   return mem_size_longint(a); /* AK 080903 */
        case COMPOSITION:
        case WORD:
        case SUBSET:
        case INTEGERVECTOR:
        case VECTOR:    return mem_size_vector(a);
        case HASHTABLE:    return mem_size_hashtable(a);
        default:
            erg += WTO("mem_size",a);goto endr_ende;
        }
    ENDR("mem_size");
}
