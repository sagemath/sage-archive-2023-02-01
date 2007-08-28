

typedef int boolean;
typedef char* string;
typedef long entry; /* see also |MaxEntry| and |MinEntry| */
typedef long index;
typedef short Short;
typedef unsigned short digit; /* see also |MaxDigit| */
typedef unsigned long twodigits;

#include "memtype.h"
#include "nodetype.h"
#include "infolrn.h"

typedef struct {
  objtype type;
  objclass class;
  char *name;
  strtype formal;
  fobject f;
  long arglst;
  long next;
} *Symblst, Symbrec;

typedef struct { strtype name; objtype type;} nametype;
typedef struct { strtype p[N_PARENTS]; short n; } par_tp;

typedef int cmp_tp;
typedef cmp_tp (*cmpfn_tp) (entry*,entry*,index);


#include  <setjmp.h>
#include  <stdio.h>
/* #include  <string.h> */
#include  <ctype.h>
#include  <stdlib.h>
#include  <signal.h>
#include  <stdarg.h>
#include  <limits.h>
#include  <assert.h>
#include  "lexer.h"
#include  "getl.h"
#include  "node.h"
#include  "mem.h"
#include  "gettype.h"
#include  "getvalue.h"
#include  "init.h"
#include  "sym.h"
#include  "main.h"
#include  "onoff.h"
#include  "ansi.h"
#define false  0
#define true   1
#define MaxEntry 		LONG_MAX
#define MinEntry 		LONG_MIN
#define MaxDigit 		((1<<15)-1)
#define max_obj_size 		UINT_MAX
#define MAXPTRS_DFLT 		99999
#define GCCRIT 		1000 \

#define MAXNODES_DFLT 	9999
#define LINELENGTH   80
#define NPOLY  10
#define NAMESTACK_LEN    500
#define REPORT_LEN       200
#define LMARGIN  5
#define PROMPT 	  "> "
#define PROMPT2 	  "\\ "
#define TO_LOOK  1
#define TO_EDIT  0
#define EXTEND 	 8
#define RANKMAXSUB  8
#define readmode  "rb"
#define writemode  "wb"
#define DEFAULT_EDITOR 	"emacs"
#define DEFAULT_PAGER 	"less"
#define array_size(a) (int)(sizeof(a)/sizeof(a[0]))
#define Integer(o) \
 (type_of(o)==INTEGER?(o)->i.intval:bigint2entry(&(o)->b))
#define is_int(t) ((t)==INTEGER ||(t)==BIGINT)
#define Max(a,b)  ((a)>=(b)?(a):(b))
#define Min(a,b)  ((a)<=(b)?(a):(b))
#define mul1(num,d)	mul1add(num,(digit)(d),(digit)0)
#define private_pol(p)  ( isshared(p) ? copypoly(p) : p )
#define C0(name,f,res) {res,OPERATOR,name,1,f,0,0},
#define C1(name,f,res,arg1) {arg1,DUMMY,NULL,0,NULL,0,0}, {res,OPERATOR,name,2,f,-1,0},
#define C2(name,f,res,arg1,arg2) \
{arg1,DUMMY,NULL,0,NULL,0,1}, {arg2,DUMMY,NULL,0,NULL,0,0}, \
{res,OPERATOR,name,3,f,-2,0},
#define C3(name,f,res,arg1,arg2,arg3) \
{arg1,DUMMY,NULL,0,NULL,0,1}, {arg2,DUMMY,NULL,0,NULL,0,1}, \
{arg3,DUMMY,NULL,0,NULL,0,0}, {res,OPERATOR,name,4,f,-3,0},
#define C4(name,f,res,arg1,arg2,arg3,arg4) \
{arg1,DUMMY,NULL,0,NULL,0,1}, {arg2,DUMMY,NULL,0,NULL,0,1}, \
{arg3,DUMMY,NULL,0,NULL,0,1}, {arg4,DUMMY,NULL,0,NULL,0,0}, \
{res,OPERATOR,name,5,f,-4,0},
#define C5(name,f,res,arg1,arg2,arg3,arg4,arg5) \
{arg1,DUMMY,NULL,0,NULL,0,1}, {arg2,DUMMY,NULL,0,NULL,0,1}, \
{arg3,DUMMY,NULL,0,NULL,0,1}, {arg4,DUMMY,NULL,0,NULL,0,1}, \
{arg5,DUMMY,NULL,0,NULL,0,0}, {res,OPERATOR,name,6,f,-5,0},
#define M0(name,f,res) {res,MAP,name,1,f,0,0},
#define M1(name,f,res,arg1) \
{arg1,DUMMY,NULL,0,NULL,0,0}, {res,MAP,name,2,f,-1,0},
#define M2(name,f,res,arg1,arg2) \
{arg1,DUMMY,NULL,0,NULL,0,1}, {arg2,DUMMY,NULL,0,NULL,0,0}, \
{res,MAP,name,3,f,-2,0},
#define M3(name,f,res,arg1,arg2,arg3) \
{arg1,DUMMY,NULL,0,NULL,0,1}, {arg2,DUMMY,NULL,0,NULL,0,1}, \
{arg3,DUMMY,NULL,0,NULL,0,0}, {res,MAP,name,4,f,-3,0},
#define M4(name,f,res,arg1,arg2,arg3,arg4) \
{arg1,DUMMY,NULL,0,NULL,0,1}, {arg2,DUMMY,NULL,0,NULL,0,1}, \
{arg3,DUMMY,NULL,0,NULL,0,1}, {arg4,DUMMY,NULL,0,NULL,0,0}, \
{res,MAP,name,5,f,-4,0},
#define M5(name,f,res,arg1,arg2,arg3,arg4,arg5) \
{arg1,DUMMY,NULL,0,NULL,0,1}, {arg2,DUMMY,NULL,0,NULL,0,1}, \
{arg3,DUMMY,NULL,0,NULL,0,1}, {arg4,DUMMY,NULL,0,NULL,0,1}, \
{arg5,DUMMY,NULL,0,NULL,0,0}, {res,MAP,name,6,f,-5,0},
#define A0(name,f,res) \
 push(top_definitions,creatopsym(0,match(name,false),f,res));
#define A1(name,f,res,arg1) \
 push(top_definitions,creatopsym(1,match(name,false),f,res,arg1));
#define A2(name,f,res,arg1,arg2) \
 push(top_definitions,creatopsym(2,match(name,false),f,res,arg1,arg2));
#define A3(name,f,res,arg1,arg2,arg3) \
 push(top_definitions,creatopsym(3,match(name,false),f,res,arg1,arg2,arg3));
#define A4(name,f,res,arg1,arg2,arg3,arg4) \
 push(top_definitions \
 ,creatopsym(4,match(name,false),f,res,arg1,arg2,arg3,arg4));
#define A5(name,f,res,arg1,arg2,arg3,arg4,arg5) \
 push(top_definitions \
 ,creatopsym(5,match(name,false),f,res,arg1,arg2,arg3,arg4,arg5));
#define FINISH  {0,DUMMY,NULL,0,NULL,0,0}

extern unsigned long  maxnodes, maxptrs, gccrit, maxenters, maxlabels;

extern strtype seq_name, if_name, dollar_name,
assign_name, assign_loc_name, break_name, return_name, block_name,
setdefault_name;

extern bigint* (*int2bin) (intcel*);
extern intcel* (*bin2int) (bigint*);
extern matrix* (*pol2mat) (poly*);
extern poly* (*mat2pol) (matrix*)
          ,* (*vec2pol) (vector*)
          ,* (*bin2pol) (bigint*)
          ,* (*int2pol) (intcel*);
extern symblst top_definitions,topsym;

extern bigint *one, *null, *minus_one;
extern intcel *bool_false, *bool_true;

extern boolean isargument, check_return;
extern int nstatic1,nstatic2,nstatic3,nstatic4,nstatic5,
    nstatic6,nstatic7;
extern Symbrec static1[],static2[],static3[],static4[],static5[],
    static6[],static7[];
extern boolean am_monitor, prompt, verbose, runtime, gc_set, bigint_set,
       lprint, parsing;
extern FILE *stderr_out, *monfile;
extern boolean with_Pre_on;
extern int tree_pt, object_pt, block_depth, lex_depth, label_pt, lmargin;

extern labeltp label, label_null; /* needed for error messages . */

extern strtype fun_name;

extern int line; /* line number needed for error messages */

extern index nrefl;

extern boolean alloc_gc;
  /* whether to use |allocmem| rather than |mlalloc| in |creatsym| */



extern nametype var[];
extern objtype return_stack[];

extern object grp, defaultgrp,	repair_obj;

extern cmpfn_tp cmpfn;

int no_terminal(FILE* f);


