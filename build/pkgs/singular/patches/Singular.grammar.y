/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/* $Id: grammar.y,v 1.123 2007/04/04 12:38:07 Singular Exp $ */
/*
* ABSTRACT: SINGULAR shell grammatik
*/
%{

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "mod2.h"
#include <mylimits.h>
#include "omalloc.h"
#include "tok.h"
#include "stype.h"
#include "ipid.h"
#include "intvec.h"
#include "febase.h"
#include "matpol.h"
#include "ring.h"
#include "kstd1.h"
#include "subexpr.h"
#include "ipshell.h"
#include "ipconv.h"
#include "sdb.h"
#include "ideals.h"
#include "numbers.h"
#include "polys.h"
#include "stairc.h"
#include "timer.h"
#include "cntrlc.h"
#include "maps.h"
#include "syz.h"
#include "lists.h"
#include "libparse.h"

#if 0
void debug_list(leftv v)
{
  idhdl r=basePackHdl;
  idhdl h;
  BOOLEAN found=FALSE;
  const char *nn=v->name;
  h=IDROOT->get(nn,myynest);
  if (h!=NULL)
  {
     Print("Curr::%s, (%s)\n",nn,Tok2Cmdname((int)IDTYP(h)));
     found=TRUE;
  }
  else         Print("`%s` not found in IDROOT\n",nn);
  while (r!=NULL)
  {
    if ((IDTYP(r)==PACKAGE_CMD)
    || (IDTYP(r)==RING_CMD)
    || (IDTYP(r)==QRING_CMD))
    {
      h=IDPACKAGE(r)->idroot->get(nn,myynest);
      if (h!=NULL)
      {
        Print("%s::%s, (%s)\n",r->id,nn,Tok2Cmdname((int)IDTYP(h)));
        found=TRUE;
      }
      else         Print("%s::%s not found\n",r->id,nn);
    }
    if (r==basePackHdl) r=IDPACKAGE(r)->idroot;
    r=r->next;
   if (r==basePackHdl) break;
  }
  if (!found)
  {
    listall(TRUE);
  }
}
#endif

/* From the bison docu:

     By defining the macro `YYMAXDEPTH', you can control how deep the
parser stack can become before a stack overflow occurs.  Define the
macro with a value that is an integer.  This value is the maximum number
of tokens that can be shifted (and not reduced) before overflow.  It
must be a constant expression whose value is known at compile time.

   The stack space allowed is not necessarily allocated.  If you
specify a large value for `YYMAXDEPTH', the parser actually allocates a
small stack at first, and then makes it bigger by stages as needed.
This increasing allocation happens automatically and silently.
Therefore, you do not need to make `YYMAXDEPTH' painfully small merely
to save space for ordinary inputs that do not need much stack.

   The default value of `YYMAXDEPTH', if you do not define it, is 10000.
*/
#define YYMAXDEPTH INT_MAX

extern int   yylineno;
extern FILE* yyin;

char       my_yylinebuf[80];
char *     currid;
BOOLEAN    yyInRingConstruction=FALSE;
BOOLEAN    expected_parms;
int        cmdtok;
int        inerror = 0;

#define TESTSETINT(a,i)                                \
   if ((a).Typ() != INT_CMD)                           \
   {                                                   \
     WerrorS("no int expression");                     \
     YYERROR;                                          \
   }                                                   \
   (i) = (int)((long)(a).Data());

#define MYYERROR(a) { WerrorS(a); YYERROR; }

void yyerror(char * fmt)
{

  BOOLEAN old_errorreported=errorreported;
  errorreported = TRUE;
  if (currid!=NULL)
  {
    killid(currid,&IDROOT);
    currid = NULL;
  }
  if(inerror==0)
  {
    #ifdef HAVE_TCL
    if (tclmode)
    { /* omit output of line number if tclmode and stdin */
      const char *n=VoiceName();
      if (strcmp(n,"STDIN")==0)
        Werror( "error occurred in %s: `%s`"
               ,n, my_yylinebuf);
      else
        Werror( "error occurred in %s line %d: `%s`"
               ,n, yylineno, my_yylinebuf);
    }
    else
    #endif
    {
      if ((strlen(fmt)>1)
      && (strncmp(fmt,"parse",5)!=0)
      && (strncmp(fmt,"syntax",6)!=0))
        WerrorS(fmt);
      Werror( "error occurred in %s line %d: `%s`"
             ,VoiceName(), yylineno, my_yylinebuf);
    }
    if (cmdtok!=0)
    {
      char *s=Tok2Cmdname(cmdtok);
      if (expected_parms)
      {
        Werror("expected %s-expression. type \'help %s;\'",s,s);
      }
      else
      {
        Werror("wrong type declaration. type \'help %s;\'",s);
      }
    }
    if (!old_errorreported && (lastreserved!=NULL))
    {
      Werror("last reserved name was `%s`",lastreserved);
    }
    inerror=1;
  }
  if ((currentVoice!=NULL)
  && (currentVoice->prev!=NULL)
  && (myynest>0)
#ifdef HAVE_SDB
  && ((sdb_flags &1)==0)
#endif
  )
  {
    Werror("leaving %s",VoiceName());
  }
#ifdef HAVE_LIBFAC_P
  extern int libfac_interruptflag;
  libfac_interruptflag=0;
#endif
}

%}

/* %expect 22 */
%pure_parser

/* special symbols */
%token DOTDOT
%token EQUAL_EQUAL
%token GE
%token LE
%token MINUSMINUS
%token NOT
%token NOTEQUAL
%token PLUSPLUS
%token COLONCOLON

/* types, part 1 (ring indep.)*/
%token <i> GRING_CMD
%token <i> INTMAT_CMD
%token <i> PROC_CMD
%token <i> RING_CMD

/* valid when ring defined ! */
%token <i> BEGIN_RING
/* types, part 2 */
%token <i> IDEAL_CMD
%token <i> MAP_CMD
%token <i> MATRIX_CMD
%token <i> MODUL_CMD
%token <i> NUMBER_CMD
%token <i> POLY_CMD
%token <i> RESOLUTION_CMD
%token <i> VECTOR_CMD
/* end types */

/* ring dependent cmd:*/
%token <i> BETTI_CMD
%token <i> COEFFS_CMD
%token <i> COEF_CMD
%token <i> CONTRACT_CMD
%token <i> DEGREE_CMD
%token <i> DEG_CMD
%token <i> DIFF_CMD
%token <i> DIM_CMD
%token <i> DIVISION_CMD
%token <i> ELIMINATION_CMD
%token <i> E_CMD
%token <i> FETCH_CMD
%token <i> FREEMODULE_CMD
%token <i> KEEPRING_CMD
%token <i> HILBERT_CMD
%token <i> HOMOG_CMD
%token <i> IMAP_CMD
%token <i> INDEPSET_CMD
%token <i> INTERRED_CMD
%token <i> INTERSECT_CMD
%token <i> JACOB_CMD
%token <i> JET_CMD
%token <i> KBASE_CMD
%token <i> KOSZUL_CMD
%token <i> LEADCOEF_CMD
%token <i> LEADEXP_CMD
%token <i> LEAD_CMD
%token <i> LEADMONOM_CMD
%token <i> LIFTSTD_CMD
%token <i> LIFT_CMD
%token <i> MAXID_CMD
%token <i> MINBASE_CMD
%token <i> MINOR_CMD
%token <i> MINRES_CMD
%token <i> MODULO_CMD
%token <i> MRES_CMD
%token <i> MULTIPLICITY_CMD
%token <i> ORD_CMD
%token <i> PAR_CMD
%token <i> PARDEG_CMD
%token <i> PREIMAGE_CMD
%token <i> QUOTIENT_CMD
%token <i> QHWEIGHT_CMD
%token <i> REDUCE_CMD
%token <i> REGULARITY_CMD
%token <i> RES_CMD
%token <i> SIMPLIFY_CMD
%token <i> SORTVEC_CMD
%token <i> SRES_CMD
%token <i> STD_CMD
%token <i> SUBST_CMD
%token <i> SYZYGY_CMD
%token <i> VAR_CMD
%token <i> VDIM_CMD
%token <i> WEDGE_CMD
%token <i> WEIGHT_CMD

/*system variables in ring block*/
%token <i> VALTVARS
%token <i> VMAXDEG
%token <i> VMAXMULT
%token <i> VNOETHER
%token <i> VMINPOLY

%token <i> END_RING
/* end of ring definitions */

%token <i> CMD_1
%token <i> CMD_2
%token <i> CMD_3
%token <i> CMD_12
%token <i> CMD_13
%token <i> CMD_23
%token <i> CMD_123
%token <i> CMD_M
%token <i> ROOT_DECL
        /* put variables of this type into the idroot list */
%token <i> ROOT_DECL_LIST
        /* put variables of this type into the idroot list */
%token <i> RING_DECL
        /* put variables of this type into the currRing list */
%token <i> EXAMPLE_CMD
%token <i> EXPORT_CMD
%token <i> HELP_CMD
%token <i> KILL_CMD
%token <i> LIB_CMD
%token <i> LISTVAR_CMD
%token <i> SETRING_CMD
%token <i> TYPE_CMD

%token <name> STRINGTOK BLOCKTOK INT_CONST
%token <name> UNKNOWN_IDENT RINGVAR PROC_DEF

/* control */
%token <i> BREAK_CMD
%token <i> CONTINUE_CMD
%token <i> ELSE_CMD
%token <i> EVAL
%token <i> QUOTE
%token <i> FOR_CMD
%token <i> IF_CMD
%token <i> SYS_BREAK
%token <i> WHILE_CMD
%token <i> RETURN
%token <i> PARAMETER

/* system variables */
%token <i> SYSVAR

%type <name> extendedid
%type <lv>   rlist ordering OrderingList orderelem
%type <name> stringexpr
%type <lv>   expr elemexpr exprlist expr_arithmetic
%type <lv>   declare_ip_variable left_value
%type <i>    ordername
%type <i>    cmdeq
%type <i>    currring_lists
%type <i>    setrings
%type <i>    ringcmd1

%type <i>    '=' '<' '>' '+' '-' COLONCOLON
%type <i>    '/' '[' ']' '^' ',' ';'


/*%nonassoc '=' PLUSEQUAL DOTDOT*/
/*%nonassoc '=' DOTDOT COLONCOLON*/
%nonassoc '=' DOTDOT
%left ','
%left '&'
%left EQUAL_EQUAL NOTEQUAL
%left '<'
%left '+' '-' ':'
%left '/' '*'
%left UMINUS NOT
%left  '^'
%left '[' ']'
%left '(' ')'
%left PLUSPLUS MINUSMINUS
%left COLONCOLON

%%
lines:
        /**/
        | lines pprompt
          {
            if (timerv)
            {
              writeTime("used time:");
              startTimer();
            }
            #ifdef HAVE_RTIMER
            if (rtimerv)
            {
              writeRTime("used real time:");
              startRTimer();
            }
            #endif
            prompt_char = '>';
#ifdef HAVE_SDB
            if (sdb_flags & 2) { sdb_flags=1; YYERROR; }
#endif
            if(siCntrlc)
            {
              siCntrlc=FALSE;
              MYYERROR("abort...");
            }
            if (errorreported)
            {
              yyerror("");
            }
            if (inerror==2) PrintLn();
            errorreported = inerror = cmdtok = 0;
            lastreserved = currid = NULL;
            expected_parms = siCntrlc = FALSE;
          }
        ;

pprompt:
        flowctrl                       /* if, while, for, proc */
        | command ';'                  /* commands returning no value */
          {currentVoice->ifsw=0;}
        | declare_ip_variable ';'      /* default initialization */
          { $1.CleanUp(); currentVoice->ifsw=0;}
        | returncmd
          {
            YYACCEPT;
          }
        | SYS_BREAK
          {
            currentVoice->ifsw=0;
            iiDebug();
          }
        | ';'                    /* ignore empty statements */
          {currentVoice->ifsw=0;}
        | error ';'
          {
            #ifdef SIQ
            siq=0;
            #endif
            yyInRingConstruction = FALSE;
            currentVoice->ifsw=0;
            if (inerror)
            {
              if ((inerror!=3) && ($1.i<UMINUS) && ($1.i>' '))
              {
                // 1: yyerror called
                // 2: scanner put actual string
                // 3: error rule put token+\n
                inerror=3;
                Print(" error at token `%s`\n",iiTwoOps($1.i));
              }
            }
            if (!errorreported) WerrorS("...parse error");
            yyerror("");
            yyerrok;
#ifdef HAVE_SDB
            if ((sdb_flags & 1) && currentVoice->pi!=NULL)
            {
              currentVoice->pi->trace_flag |=1;
            }
            else
#endif
            if (myynest>0)
            {
              feBufferTypes t=currentVoice->Typ();
              //PrintS("leaving yyparse\n");
              exitBuffer(BT_proc);
              if (t==BT_example)
                YYACCEPT;
              else
                YYABORT;
            }
            else if (currentVoice->prev!=NULL)
            {
              exitVoice();
            }
#ifdef HAVE_SDB
            if (sdb_flags &2) sdb_flags=1;
#endif
          }
        ;

flowctrl: ifcmd
          | whilecmd
          | example_dummy
          | forcmd
          | proccmd
          | filecmd
          | helpcmd
          | examplecmd
            {if (currentVoice!=NULL) currentVoice->ifsw=0;}
        ;

example_dummy : EXAMPLE_CMD BLOCKTOK { omFree((ADDRESS)$2); }
        ;

command: assign
         | exportcmd
         | killcmd
         | listcmd
         | parametercmd
         | ringcmd
         | scriptcmd
         | setringcmd
         | typecmd
         ;

assign: left_value exprlist
          {
            if(iiAssign(&$1,&$2)) YYERROR;
          }
        ;

elemexpr:
        RINGVAR
          {
            if (currRing==NULL) MYYERROR("no ring active");
            syMake(&$$,omStrDup($1));
          }
        | extendedid
          {
            syMake(&$$,$1);
          }
        | elemexpr COLONCOLON elemexpr
          {
            if(iiExprArith2(&$$, &$1, COLONCOLON, &$3)) YYERROR;
          }
        | elemexpr '('  ')'
          {
            if(iiExprArith1(&$$,&$1,'(')) YYERROR;
          }
        | elemexpr '(' exprlist ')'
          {
            if ($1.rtyp==LIB_CMD)
            {
              if(iiExprArith1(&$$,&$3,LIB_CMD)) YYERROR;
            }
            else
            {
              if ($1.Typ()==UNKNOWN)
              {
                if(iiExprArith2(&$$,&$1,'(',&$3)) YYERROR;
              }
              else
              {
                $1.next=(leftv)omAllocBin(sleftv_bin);
                memcpy($1.next,&$3,sizeof(sleftv));
                if(iiExprArithM(&$$,&$1,'(')) YYERROR;
              }
            }
          }
        | '[' exprlist ']'
          {
            if (currRingHdl==NULL) MYYERROR("no ring active");
            int j = 0;
            memset(&$$,0,sizeof(sleftv));
            $$.rtyp=VECTOR_CMD;
            leftv v = &$2;
            while (v!=NULL)
            {
              int i,t;
              sleftv tmp;
              memset(&tmp,0,sizeof(tmp));
              i=iiTestConvert((t=v->Typ()),POLY_CMD);
              if((i==0) || (iiConvert(t /*v->Typ()*/,POLY_CMD,i,v,&tmp)))
              {
                pDelete((poly *)&$$.data);
                $2.CleanUp();
                MYYERROR("expected '[poly,...'");
              }
              poly p = (poly)tmp.CopyD(POLY_CMD);
              pSetCompP(p,++j);
              $$.data = (void *)pAdd((poly)$$.data,p);
              v->next=tmp.next;tmp.next=NULL;
              tmp.CleanUp();
              v=v->next;
            }
            $2.CleanUp();
          }
        | INT_CONST
          {
            memset(&$$,0,sizeof($$));
            int i = atoi($1);
            /*remember not to omFree($1)
            *because it is a part of the scanner buffer*/
            $$.rtyp  = INT_CMD;
            $$.data = (void *)(long)i;

            /* check: out of range input */
            int l = strlen($1)+2;
            if (l >= MAX_INT_LEN)
            {
              char tmp[MAX_INT_LEN+5];
              sprintf(tmp,"%d",i);
              if (strcmp(tmp,$1)!=0)
              {
                if (currRing==NULL)
                {
                  Werror("`%s` greater than %d(max. integer representation)"
                         ,$1,MAX_INT_VAL);
                  YYERROR;
                }
                char *t1=omStrDup($1);
                syMake(&$$,t1);
              }
            }
          }
        | SYSVAR
          {
            memset(&$$,0,sizeof($$));
            $$.rtyp = $1;
            $$.data = $$.Data();
          }
        | stringexpr
          {
            memset(&$$,0,sizeof($$));
            $$.rtyp  = STRING_CMD;
            $$.data = $1;
          }
        ;

exprlist:
        exprlist ',' expr
          {
            leftv v = &$1;
            while (v->next!=NULL)
            {
              v=v->next;
            }
            v->next = (leftv)omAllocBin(sleftv_bin);
            memcpy(v->next,&($3),sizeof(sleftv));
            $$ = $1;
          }
        | expr
          {
            $$ = $1;
          }
        ;

expr:   expr_arithmetic
          {
            /*if ($1.typ == eunknown) YYERROR;*/
            $$ = $1;
          }
        | elemexpr       { $$ = $1; }
        | '(' exprlist ')'    { $$ = $2; }
        | expr '[' expr ',' expr ']'
          {
            if(iiExprArith3(&$$,'[',&$1,&$3,&$5)) YYERROR;
          }
        | expr '[' expr ']'
          {
            if(iiExprArith2(&$$,&$1,'[',&$3)) YYERROR;
          }
        | ROOT_DECL '(' expr ')'
          {
            if(iiExprArith1(&$$,&$3,$1)) YYERROR;
          }
        | ROOT_DECL_LIST '(' exprlist ')'
          {
            if(iiExprArithM(&$$,&$3,$1)) YYERROR;
          }
        | ROOT_DECL_LIST '(' ')'
          {
            if(iiExprArithM(&$$,NULL,$1)) YYERROR;
          }
        | RING_DECL '(' expr ')'
          {
            if(iiExprArith1(&$$,&$3,$1)) YYERROR;
          }
        | currring_lists '(' exprlist ')'
          {
            if(iiExprArithM(&$$,&$3,$1)) YYERROR;
          }
        | currring_lists '(' ')'
          {
            if(iiExprArithM(&$$,NULL,$1)) YYERROR;
          }
        | CMD_1 '(' expr ')'
          {
            if(iiExprArith1(&$$,&$3,$1)) YYERROR;
          }
        | CMD_2 '(' expr ',' expr ')'
          {
            if(iiExprArith2(&$$,&$3,$1,&$5,TRUE)) YYERROR;
          }
        | CMD_3 '(' expr ',' expr ',' expr ')'
          {
            if(iiExprArith3(&$$,$1,&$3,&$5,&$7)) YYERROR;
          }
        | CMD_23 '(' expr ',' expr ')'
          {
            if(iiExprArith2(&$$,&$3,$1,&$5,TRUE)) YYERROR;
          }
        | CMD_23 '(' expr ',' expr ',' expr ')'
          {
            if(iiExprArith3(&$$,$1,&$3,&$5,&$7)) YYERROR;
          }
        | CMD_12 '(' expr ')'
          {
            if(iiExprArith1(&$$,&$3,$1)) YYERROR;
          }
        | CMD_13 '(' expr ')'
          {
            if(iiExprArith1(&$$,&$3,$1)) YYERROR;
          }
        | CMD_12 '(' expr ',' expr ')'
          {
            if(iiExprArith2(&$$,&$3,$1,&$5,TRUE)) YYERROR;
          }
        | CMD_123 '(' expr ')'
          {
            if(iiExprArith1(&$$,&$3,$1)) YYERROR;
          }
        | CMD_123 '(' expr ',' expr ')'
          {
            if(iiExprArith2(&$$,&$3,$1,&$5,TRUE)) YYERROR;
          }
        | CMD_13 '(' expr ',' expr ',' expr ')'
          {
            if(iiExprArith3(&$$,$1,&$3,&$5,&$7)) YYERROR;
          }
        | CMD_123 '(' expr ',' expr ',' expr ')'
          {
            if(iiExprArith3(&$$,$1,&$3,&$5,&$7)) YYERROR;
          }
        | CMD_M '(' ')'
          {
            if(iiExprArithM(&$$,NULL,$1)) YYERROR;
          }
        | CMD_M '(' exprlist ')'
          {
            if(iiExprArithM(&$$,&$3,$1)) YYERROR;
          }
        | MATRIX_CMD '(' expr ',' expr ',' expr ')'
          {
            if(iiExprArith3(&$$,MATRIX_CMD,&$3,&$5,&$7)) YYERROR;
          }
        | MATRIX_CMD '(' expr ')'
          {
            if(iiExprArith1(&$$,&$3,MATRIX_CMD)) YYERROR;
          }
        | INTMAT_CMD '(' expr ',' expr ',' expr ')'
          {
            if(iiExprArith3(&$$,INTMAT_CMD,&$3,&$5,&$7)) YYERROR;
          }
        | INTMAT_CMD '(' expr ')'
          {
            if(iiExprArith1(&$$,&$3,INTMAT_CMD)) YYERROR;
          }
        | RING_CMD '(' rlist ',' rlist ',' ordering ')'
          {
            if(iiExprArith3(&$$,RING_CMD,&$3,&$5,&$7)) YYERROR;
          }
        | RING_CMD '(' expr ')'
          {
            if(iiExprArith1(&$$,&$3,RING_CMD)) YYERROR;
          }
        | quote_start expr quote_end
          {
            $$=$2;
          }
        | quote_start expr '=' expr quote_end
          {
            #ifdef SIQ
            siq++;
            if (siq>0)
            { if (iiExprArith2(&$$,&$2,'=',&$4)) YYERROR; }
            else
            #endif
            {
              memset(&$$,0,sizeof($$));
              $$.rtyp=NONE;
              if (iiAssign(&$2,&$4)) YYERROR;
            }
            #ifdef SIQ
            siq--;
            #endif
          }
        | EVAL  '('
          {
            #ifdef SIQ
            siq--;
            #endif
          }
          expr ')'
          {
            #ifdef SIQ
            if (siq<=0) $4.Eval();
            #endif
            $$=$4;
            #ifdef SIQ
            siq++;
            #endif
          }
          ;

quote_start:    QUOTE  '('
          {
            #ifdef SIQ
            siq++;
            #endif
          }
          ;

quote_end: ')'
          {
            #ifdef SIQ
            siq--;
            #endif
          }
          ;

expr_arithmetic:
          expr PLUSPLUS     %prec PLUSPLUS
          {
            if(iiExprArith1(&$$,&$1,PLUSPLUS)) YYERROR;
          }
        | expr MINUSMINUS   %prec MINUSMINUS
          {
            if(iiExprArith1(&$$,&$1,MINUSMINUS)) YYERROR;
          }
        | expr '+' expr
          {
            if(iiExprArith2(&$$,&$1,'+',&$3)) YYERROR;
          }
        | expr '-' expr
          {
            if(iiExprArith2(&$$,&$1,'-',&$3)) YYERROR;
          }
        | expr '/' expr
          {
            if(iiExprArith2(&$$,&$1,$<i>2,&$3)) YYERROR;
          }
        | expr '^' expr
          {
            if(iiExprArith2(&$$,&$1,'^',&$3)) YYERROR;
          }
        | expr '<' expr
          {
            if(iiExprArith2(&$$,&$1,$<i>2,&$3)) YYERROR;
          }
        | expr '&' expr
          {
            if(iiExprArith2(&$$,&$1,$<i>2,&$3)) YYERROR;
          }
        | expr NOTEQUAL expr
          {
            if(iiExprArith2(&$$,&$1,NOTEQUAL,&$3)) YYERROR;
          }
        | expr EQUAL_EQUAL expr
          {
            if(iiExprArith2(&$$,&$1,EQUAL_EQUAL,&$3)) YYERROR;
          }
        | expr DOTDOT expr
          {
            if(iiExprArith2(&$$,&$1,DOTDOT,&$3)) YYERROR;
          }
        | expr ':' expr
          {
            if(iiExprArith2(&$$,&$1,':',&$3)) YYERROR;
          }
        | NOT expr
          {
            memset(&$$,0,sizeof($$));
            int i; TESTSETINT($2,i);
            $$.rtyp  = INT_CMD;
            $$.data = (void *)(long)(i == 0 ? 1 : 0);
          }
        | '-' expr %prec UMINUS
          {
            if(iiExprArith1(&$$,&$2,'-')) YYERROR;
          }
        ;

left_value:
        declare_ip_variable cmdeq  { $$ = $1; }
        | exprlist '='
          {
            if ($1.rtyp==0)
            {
              Werror("`%s` is undefined",$1.Fullname());
              YYERROR;
            }
            $$ = $1;
          }
        ;


extendedid:
        UNKNOWN_IDENT
        | '`' expr '`'
          {
            if ($2.Typ()!=STRING_CMD)
            {
              MYYERROR("string expression expected");
            }
            $$ = (char *)$2.CopyD(STRING_CMD);
            $2.CleanUp();
          }
        ;

currring_lists:
        IDEAL_CMD | MODUL_CMD
        /* put variables into the current ring */
        ;

declare_ip_variable:
        ROOT_DECL elemexpr
          {
            #ifdef HAVE_NS
            if (iiDeclCommand(&$$,&$2,myynest,$1,&($2.req_packhdl->idroot)))
              YYERROR;
            #else
            if (iiDeclCommand(&$$,&$2,myynest,$1,&IDROOT)) YYERROR;
            #endif
          }
        | ROOT_DECL_LIST elemexpr
          {
            #ifdef HAVE_NS
            if (iiDeclCommand(&$$,&$2,myynest,$1,&($2.req_packhdl->idroot)))
              YYERROR;
            #else
            if (iiDeclCommand(&$$,&$2,myynest,$1,&IDROOT)) YYERROR;
            #endif
          }
        | RING_DECL elemexpr
          {
            if (iiDeclCommand(&$$,&$2,myynest,$1,&(currRing->idroot), TRUE)) YYERROR;
          }
        | currring_lists elemexpr
          {
            if (iiDeclCommand(&$$,&$2,myynest,$1,&(currRing->idroot), TRUE)) YYERROR;
          }
        | MATRIX_CMD elemexpr '[' expr ']' '[' expr ']'
          {
            if (iiDeclCommand(&$$,&$2,myynest,$1,&(currRing->idroot), TRUE)) YYERROR;
            int r; TESTSETINT($4,r);
            int c; TESTSETINT($7,c);
            if (r < 1)
              MYYERROR("rows must be greater than 0");
            if (c < 0)
              MYYERROR("cols must be greater than -1");
            leftv v=&$$;
            //while (v->next!=NULL) { v=v->next; }
            idhdl h=(idhdl)v->data;
            idDelete(&IDIDEAL(h));
            IDMATRIX(h) = mpNew(r,c);
            if (IDMATRIX(h)==NULL) YYERROR;
          }
        | MATRIX_CMD elemexpr
          {
            if (iiDeclCommand(&$$,&$2,myynest,$1,&(currRing->idroot), TRUE)) YYERROR;
          }
        | INTMAT_CMD elemexpr '[' expr ']' '[' expr ']'
          {
            int r; TESTSETINT($4,r);
            int c; TESTSETINT($7,c);
            if (r < 1)
              MYYERROR("rows must be greater than 0");
            if (c < 0)
              MYYERROR("cols must be greater than -1");
            #ifdef HAVE_NS
            if (iiDeclCommand(&$$,&$2,myynest,$1,&($2.req_packhdl->idroot)))
              YYERROR;
            #else
            if (iiDeclCommand(&$$,&$2,myynest,$1,&IDROOT)) YYERROR;
            #endif
            leftv v=&$$;
            idhdl h=(idhdl)v->data;
            delete IDINTVEC(h);
            IDINTVEC(h) = new intvec(r,c,0);
            if (IDINTVEC(h)==NULL) YYERROR;
          }
        | INTMAT_CMD elemexpr
          {
            #ifdef HAVE_NS
            if (iiDeclCommand(&$$,&$2,myynest,$1,&($2.req_packhdl->idroot)))
              YYERROR;
            #else
            if (iiDeclCommand(&$$,&$2,myynest,$1,&IDROOT)) YYERROR;
            #endif
            leftv v=&$$;
            idhdl h;
            do
            {
               h=(idhdl)v->data;
               delete IDINTVEC(h);
               IDINTVEC(h) = new intvec(1,1,0);
               v=v->next;
            } while (v!=NULL);
          }
        | declare_ip_variable ',' elemexpr
          {
            int t=$1.Typ();
            sleftv r;
            memset(&r,0,sizeof(sleftv));
            if ((BEGIN_RING<t) && (t<END_RING))
            {
              if (iiDeclCommand(&r,&$3,myynest,t,&(currRing->idroot), TRUE))
                YYERROR;
            }
            else
            {
              #ifdef HAVE_NS
              if (iiDeclCommand(&r,&$3,myynest,t,&($3.req_packhdl->idroot)))
                YYERROR;
              #else
              if (iiDeclCommand(&r,&$3,myynest,t,&IDROOT)) YYERROR;
              #endif
            }
            leftv v=&$1;
            while (v->next!=NULL) v=v->next;
            v->next=(leftv)omAllocBin(sleftv_bin);
            memcpy(v->next,&r,sizeof(sleftv));
            $$=$1;
          }
        | PROC_CMD elemexpr
          {
            #ifdef HAVE_NS
            if (iiDeclCommand(&$$,&$2,myynest,$1,&($2.req_packhdl->idroot)))
              YYERROR;
            #else
            if (iiDeclCommand(&$$,&$2,myynest,$1,&IDROOT)) YYERROR;
            #endif
          }
        ;

stringexpr:
        STRINGTOK
        ;

rlist:
        expr
        | '(' expr ',' exprlist ')'
          {
            leftv v = &$2;
            while (v->next!=NULL)
            {
              v=v->next;
            }
            v->next = (leftv)omAllocBin(sleftv_bin);
            memcpy(v->next,&($4),sizeof(sleftv));
            $$ = $2;
          }
        ;

ordername:
        UNKNOWN_IDENT
        {
          // let rInit take care of any errors
          $$=rOrderName($1);
        }
        ;

orderelem:
        ordername
          {
            memset(&$$,0,sizeof($$));
            intvec *iv = new intvec(2);
            (*iv)[0] = 1;
            (*iv)[1] = $1;
            $$.rtyp = INTVEC_CMD;
            $$.data = (void *)iv;
          }
        | ordername '(' exprlist ')'
          {
            memset(&$$,0,sizeof($$));
            leftv sl = &$3;
            int slLength;
            {
              slLength =  exprlist_length(sl);
              int l = 2 +  slLength;
              intvec *iv = new intvec(l);
              (*iv)[0] = slLength;
              (*iv)[1] = $1;

              int i = 2;
              while ((i<l) && (sl!=NULL))
              {
                if (sl->Typ() == INT_CMD)
                {
                  (*iv)[i++] = (int)((long)(sl->Data()));
                }
                else if ((sl->Typ() == INTVEC_CMD)
                ||(sl->Typ() == INTMAT_CMD))
                {
                  intvec *ivv = (intvec *)(sl->Data());
                  int ll = 0,l = ivv->length();
                  for (; l>0; l--)
                  {
                    (*iv)[i++] = (*ivv)[ll++];
                  }
                }
                else
                {
                  delete iv;
                  $3.CleanUp();
                  MYYERROR("wrong type in ordering");
                }
                sl = sl->next;
              }
              $$.rtyp = INTVEC_CMD;
              $$.data = (void *)iv;
            }
            $3.CleanUp();
          }
        ;

OrderingList:
        orderelem
        |  orderelem ',' OrderingList
          {
            $$ = $1;
            $$.next = (sleftv *)omAllocBin(sleftv_bin);
            memcpy($$.next,&$3,sizeof(sleftv));
          }
        ;

ordering:
        orderelem
        | '(' OrderingList ')'
          {
            $$ = $2;
          }
        ;

cmdeq:  '='
          {
            expected_parms = TRUE;
          }
        ;


/* --------------------------------------------------------------------*/
/* section of pure commands                                            */
/* --------------------------------------------------------------------*/

filecmd:
        '<' stringexpr
          { if ($<i>1 != '<') YYERROR;
            if((feFilePending=feFopen($2,"r",NULL,TRUE))==NULL) YYERROR; }
        ';'
          { newFile($2,feFilePending); }
        ;

helpcmd:
        HELP_CMD STRINGTOK ';'
          {
            feHelp($2);
            omFree((ADDRESS)$2);
          }
        | HELP_CMD ';'
          {
            feHelp(NULL);
          }
        ;

examplecmd:
        EXAMPLE_CMD STRINGTOK ';'
          {
            singular_example($2);
            omFree((ADDRESS)$2);
          }
       ;

exportcmd:
        EXPORT_CMD exprlist
        {
#ifdef HAVE_NS
          if (basePack!=$2.req_packhdl)
          {
            if(iiExport(&$2,0,currPackHdl)) YYERROR;
          }
          else
#endif /* HAVE_NS */
            if (iiExport(&$2,0)) YYERROR;
        }
        ;

killcmd:
        KILL_CMD elemexpr
        {
          leftv v=&$2;
          if (v->rtyp!=IDHDL)
          {
            if (v->name!=NULL)
            {
               Werror("`%s` is undefined in kill",v->name);
            }
            else               WerrorS("kill what ?");
          }
          else
          {
            #ifdef HAVE_NS
            killhdl((idhdl)v->data,v->req_packhdl);
            #else
            killhdl((idhdl)v->data);
            #endif
          }
        }
        | killcmd ',' elemexpr
        {
          leftv v=&$3;
          if (v->rtyp!=IDHDL)
          {
            if (v->name!=NULL)
            {
               Werror("`%s` is undefined in kill",v->name);
            }
            else               WerrorS("kill what ?");
          }
          else
          {
            #ifdef HAVE_NS
            killhdl((idhdl)v->data,v->req_packhdl);
            #else
            killhdl((idhdl)v->data);
            #endif
          }
        }
        ;

listcmd:
        LISTVAR_CMD '(' ROOT_DECL ')'
          {
            list_cmd($3,NULL,"// ",TRUE);
          }
        | LISTVAR_CMD '(' ROOT_DECL_LIST ')'
          {
            list_cmd($3,NULL,"// ",TRUE);
          }
        | LISTVAR_CMD '(' RING_DECL ')'
          {
            if ($3==QRING_CMD) $3=RING_CMD;
            list_cmd($3,NULL,"// ",TRUE);
          }
        | LISTVAR_CMD '(' currring_lists ')'
          {
            list_cmd($3,NULL,"// ",TRUE);
          }
        | LISTVAR_CMD '(' RING_CMD ')'
          {
            list_cmd(RING_CMD,NULL,"// ",TRUE);
          }
        | LISTVAR_CMD '(' MATRIX_CMD ')'
          {
            list_cmd(MATRIX_CMD,NULL,"// ",TRUE);
           }
        | LISTVAR_CMD '(' INTMAT_CMD ')'
          {
            list_cmd(INTMAT_CMD,NULL,"// ",TRUE);
          }
        | LISTVAR_CMD '(' PROC_CMD ')'
          {
            list_cmd(PROC_CMD,NULL,"// ",TRUE);
          }
        | LISTVAR_CMD '(' elemexpr ')'
          {
            list_cmd(0,$3.Fullname(),"// ",TRUE);
            $3.CleanUp();
          }
        | LISTVAR_CMD '(' elemexpr ',' ROOT_DECL ')'
          {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if($3.Typ() == PACKAGE_CMD)
              list_cmd($5,NULL,"// ",TRUE);
#endif /* HAVE_NS */
            $3.CleanUp();
          }
        | LISTVAR_CMD '(' elemexpr ',' ROOT_DECL_LIST ')'
          {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if($3.Typ() == PACKAGE_CMD)
              list_cmd($5,NULL,"// ",TRUE);
#endif /* HAVE_NS */
            $3.CleanUp();
          }
        | LISTVAR_CMD '(' elemexpr ',' RING_DECL ')'
          {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if($3.Typ() == PACKAGE_CMD)
              list_cmd($5,NULL,"// ",TRUE);
#endif /* HAVE_NS */
            $3.CleanUp();
          }
        | LISTVAR_CMD '(' elemexpr ',' currring_lists ')'
          {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if($3.Typ() == PACKAGE_CMD)
              list_cmd($5,NULL,"// ",TRUE);
#endif /* HAVE_NS */
            $3.CleanUp();
          }
        | LISTVAR_CMD '(' elemexpr ',' RING_CMD ')'
          {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if($3.Typ() == PACKAGE_CMD)
              list_cmd($5,NULL,"// ",TRUE);
#endif /* HAVE_NS */
            $3.CleanUp();
          }
        | LISTVAR_CMD '(' elemexpr ',' MATRIX_CMD ')'
          {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if($3.Typ() == PACKAGE_CMD)
              list_cmd($5,NULL,"// ",TRUE);
#endif /* HAVE_NS */
            $3.CleanUp();
          }
        | LISTVAR_CMD '(' elemexpr ',' INTMAT_CMD ')'
          {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if($3.Typ() == PACKAGE_CMD)
              list_cmd($5,NULL,"// ",TRUE);
#endif /* HAVE_NS */
            $3.CleanUp();
          }
        | LISTVAR_CMD '(' elemexpr ',' PROC_CMD ')'
          {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if($3.Typ() == PACKAGE_CMD)
              list_cmd($5,NULL,"// ",TRUE);
#endif /* HAVE_NS */
            $3.CleanUp();
          }
        //| LISTVAR_CMD '(' elemexpr ',' elemexpr ')'
        //  {
        //#ifdef HAVE_NS
        //    //PrintS("?????\n");
        //    //if($3.Typ() == PACKAGE_CMD)
        //    //  list_cmd($5,NULL,"// ",TRUE);
        //#endif /* HAVE_NS */
        //    $3.CleanUp();
        //  }
        | LISTVAR_CMD '(' ')'
          {
            list_cmd(-1,NULL,"// ",TRUE);
          }
        ;

ringcmd1:
       RING_CMD { yyInRingConstruction = TRUE; }
       ;

ringcmd:
        ringcmd1
          elemexpr cmdeq
          rlist     ','      /* description of coeffs */
          rlist     ','      /* var names */
          ordering           /* list of (multiplier ordering (weight(s))) */
          {
            BOOLEAN do_pop = FALSE;
            char *ring_name = $2.name;
            ring b=
            rInit(&$4,            /* characteristik and list of parameters*/
                  &$6,            /* names of ringvariables */
                  &$8);            /* ordering */
            idhdl newRingHdl=NULL;

            if (b!=NULL)
            {
              #ifdef HAVE_NS
                newRingHdl=enterid(ring_name, myynest, RING_CMD,
                                   &($2.req_packhdl->idroot));
              #else
                newRingHdl=enterid(ring_name, myynest, RING_CMD, &IDROOT);
              #endif
              $2.CleanUp();
              if (newRingHdl!=NULL)
              {
                omFreeSize(IDRING(newRingHdl),sizeof(ip_sring));
                IDRING(newRingHdl)=b;
              }
              else
              {
                rKill(b);
              }
            }
            yyInRingConstruction = FALSE;
            if (newRingHdl==NULL)
            {
              MYYERROR("cannot make ring");
            }
            else
            {
              rSetHdl(newRingHdl);
            }
          }
        | ringcmd1 elemexpr
          {
            BOOLEAN do_pop = FALSE;
            char *ring_name = $2.name;
            if (!inerror) rDefault(ring_name);
            yyInRingConstruction = FALSE;
            $2.CleanUp();
          }
        ;

scriptcmd:
         SYSVAR stringexpr
          {
            if (($1!=LIB_CMD)||(iiLibCmd($2,TRUE,TRUE,TRUE)))
            //if ($1==LIB_CMD)
            //{
            //  sleftv tmp;
            //  if(iiExprArith1(&tmp,&$2,LIB_CMD)) YYERROR;
            //}
            //else
                YYERROR;
          }
        ;

setrings:  SETRING_CMD | KEEPRING_CMD ;

setringcmd:
        setrings expr
          {
            if (($1==KEEPRING_CMD) && (myynest==0))
               MYYERROR("only inside a proc allowed");
            const char * n=$2.Name();
            if ((($2.Typ()==RING_CMD)||($2.Typ()==QRING_CMD))
            && ($2.rtyp==IDHDL))
            {
              idhdl h=(idhdl)$2.data;
              if ($2.e!=NULL) h=rFindHdl((ring)$2.Data(),NULL, NULL);
              //Print("setring %s lev %d (ptr:%x)\n",IDID(h),IDLEV(h),IDRING(h));
              if ($1==KEEPRING_CMD)
              {
                if (h!=NULL)
                {
                  if (IDLEV(h)!=0)
                  {
                    if (iiExport(&$2,myynest-1)) YYERROR;
#if 1
                    //if (TEST_OPT_KEEPVARS)
                    //{
                      idhdl p=IDRING(h)->idroot;
                      idhdl root=p;
                      int prevlev=myynest-1;
                      while (p!=NULL)
                      {
                        if (IDLEV(p)==myynest)
                        {
                          idhdl old=root->get(IDID(p),prevlev);
                          if (old!=NULL)
                          {
                            if (BVERBOSE(V_REDEFINE))
                              Warn("redefining %s",IDID(p));
                            killhdl2(old,&root,IDRING(h));
                            IDRING(h)->idroot=root;
                          }
                          IDLEV(p)=prevlev;
                        }
                        p=IDNEXT(p);
                      }
                      //IDRING(h)->idroot=root;
#endif
                    //}
                  }
#ifdef USE_IILOCALRING
                  iiLocalRing[myynest-1]=IDRING(h);
#endif
                  procstack->cRing=IDRING(h);
                  procstack->cRingHdl=h;
                }
                else
                {
                  Werror("%s is no identifier",n);
                  $2.CleanUp();
                  YYERROR;
                }
              }
              if (h!=NULL) rSetHdl(h);
              else
              {
                Werror("cannot find the name of the basering %s",n);
                $2.CleanUp();
                YYERROR;
              }
              $2.CleanUp();
            }
            else
            {
              Werror("%s is no name of a ring/qring",n);
              $2.CleanUp();
              YYERROR;
            }
          }
        ;

typecmd:
        TYPE_CMD expr
          {
            if ($2.rtyp!=IDHDL) MYYERROR("identifier expected");
            idhdl h = (idhdl)$2.data;
            type_cmd(h);
          }
        | exprlist
          {
            //Print("typ is %d, rtyp:%d\n",$1.Typ(),$1.rtyp);
            #ifdef SIQ
            if ($1.rtyp!=COMMAND)
            {
            #endif
              if ($1.Typ()==UNKNOWN)
              {
                if ($1.name!=NULL)
                {
                  Werror("`%s` is undefined",$1.name);
                  omFree((ADDRESS)$1.name);
                }
                YYERROR;
              }
            #ifdef SIQ
            }
            #endif
            $1.Print(&sLastPrinted);
            $1.CleanUp(currRing);
            if (errorreported) YYERROR;
          }
        ;

/* --------------------------------------------------------------------*/
/* section of flow control                                             */
/* --------------------------------------------------------------------*/

ifcmd: IF_CMD '(' expr ')' BLOCKTOK
          {
            int i; TESTSETINT($3,i);
            if (i!=0)
            {
              newBuffer( $5, BT_if);
            }
            else
            {
              omFree((ADDRESS)$5);
              currentVoice->ifsw=1;
            }
          }
        | ELSE_CMD BLOCKTOK
          {
            if (currentVoice->ifsw==1)
            {
              currentVoice->ifsw=0;
              newBuffer( $2, BT_else);
            }
            else
            {
              if (currentVoice->ifsw!=2)
              {
                Warn("`else` without `if` in level %d",myynest);
              }
              omFree((ADDRESS)$2);
            }
            currentVoice->ifsw=0;
          }
        | IF_CMD '(' expr ')' BREAK_CMD
          {
            int i; TESTSETINT($3,i);
            if (i)
            {
              if (exitBuffer(BT_break)) YYERROR;
            }
            currentVoice->ifsw=0;
          }
        | BREAK_CMD
          {
            if (exitBuffer(BT_break)) YYERROR;
            currentVoice->ifsw=0;
          }
        | CONTINUE_CMD
          {
            if (contBuffer(BT_break)) YYERROR;
            currentVoice->ifsw=0;
          }
      ;

whilecmd:
        WHILE_CMD STRINGTOK BLOCKTOK
          {
            /* -> if(!$2) break; $3; continue;*/
            char * s = (char *)omAlloc( strlen($2) + strlen($3) + 36);
            sprintf(s,"whileif (!(%s)) break;\n%scontinue;\n " ,$2,$3);
            newBuffer(s,BT_break);
            omFree((ADDRESS)$2);
            omFree((ADDRESS)$3);
          }
        ;

forcmd:
        FOR_CMD STRINGTOK STRINGTOK STRINGTOK BLOCKTOK
          {
            /* $2 */
            /* if (!$3) break; $5; $4; continue; */
            char * s = (char *)omAlloc( strlen($3)+strlen($4)+strlen($5)+36);
            sprintf(s,"forif (!(%s)) break;\n%s%s;\ncontinue;\n "
                   ,$3,$5,$4);
            omFree((ADDRESS)$3);
            omFree((ADDRESS)$4);
            omFree((ADDRESS)$5);
            newBuffer(s,BT_break);
            s = (char *)omAlloc( strlen($2) + 3);
            sprintf(s,"%s;\n",$2);
            omFree((ADDRESS)$2);
            newBuffer(s,BT_if);
          }
        ;

proccmd:
        PROC_CMD extendedid BLOCKTOK
          {
            procinfov pi;
            idhdl h = enterid($2,myynest,PROC_CMD,&IDROOT,TRUE);
            if (h==NULL) {omFree((ADDRESS)$2);omFree((ADDRESS)$3); YYERROR;}
            iiInitSingularProcinfo(IDPROC(h),"", $2, 0, 0);
            IDPROC(h)->data.s.body = (char *)omAlloc(strlen($3)+31);;
            sprintf(IDPROC(h)->data.s.body,"parameter list #;\n%s;return();\n\n",$3);
            omFree((ADDRESS)$3);
            omFree((ADDRESS)$2);
          }
        | PROC_DEF STRINGTOK BLOCKTOK
          {
            idhdl h = enterid($1,myynest,PROC_CMD,&IDROOT,TRUE);
            if (h==NULL)
            {
              omFree((ADDRESS)$1);
              omFree((ADDRESS)$2);
              omFree((ADDRESS)$3);
              YYERROR;
            }
            char *args=iiProcArgs($2,FALSE);
            omFree((ADDRESS)$2);
            procinfov pi;
            iiInitSingularProcinfo(IDPROC(h),"", $1, 0, 0);
            IDPROC(h)->data.s.body = (char *)omAlloc(strlen($3)+strlen(args)+14);;
            sprintf(IDPROC(h)->data.s.body,"%s\n%s;return();\n\n",args,$3);
            omFree((ADDRESS)args);
            omFree((ADDRESS)$3);
            omFree((ADDRESS)$1);
          }
        | PROC_DEF STRINGTOK STRINGTOK BLOCKTOK
          {
            omFree((ADDRESS)$3);
            idhdl h = enterid($1,myynest,PROC_CMD,&IDROOT,TRUE);
            if (h==NULL)
            {
              omFree((ADDRESS)$1);
              omFree((ADDRESS)$2);
              omFree((ADDRESS)$4);
              YYERROR;
            }
            char *args=iiProcArgs($2,FALSE);
            omFree((ADDRESS)$2);
            procinfov pi;
            iiInitSingularProcinfo(IDPROC(h),"", $1, 0, 0);
            omFree((ADDRESS)$1);
            IDPROC(h)->data.s.body = (char *)omAlloc(strlen($4)+strlen(args)+14);;
            sprintf(IDPROC(h)->data.s.body,"%s\n%s;return();\n\n",args,$4);
            omFree((ADDRESS)args);
            omFree((ADDRESS)$4);
          }
        ;

parametercmd:
        PARAMETER declare_ip_variable
          {
            //Print("par:%s, %d\n",$2.Name(),$2.Typ());
            //yylineno--;
            if (iiParameter(&$2)) YYERROR;
          }
        | PARAMETER expr
          {
            //Print("par:%s, %d\n",$2.Name(),$2.Typ());
            sleftv tmp_expr;
            //yylineno--;
            if ((iiDeclCommand(&tmp_expr,&$2,myynest,DEF_CMD,&IDROOT))
            || (iiParameter(&tmp_expr)))
              YYERROR;
          }
        ;

returncmd:
        RETURN '(' exprlist ')'
          {
            if(iiRETURNEXPR==NULL) YYERROR;
            iiRETURNEXPR[myynest].Copy(&$3);
            $3.CleanUp();
            if (exitBuffer(BT_proc)) YYERROR;
          }
        | RETURN '(' ')'
          {
            if ($1==RETURN)
            {
              if(iiRETURNEXPR==NULL) YYERROR;
              iiRETURNEXPR[myynest].Init();
              iiRETURNEXPR[myynest].rtyp=NONE;
              if (exitBuffer(BT_proc)) YYERROR;
            }
          }
        ;
