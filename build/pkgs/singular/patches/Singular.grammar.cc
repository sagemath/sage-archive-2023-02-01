/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     DOTDOT = 258,
     EQUAL_EQUAL = 259,
     GE = 260,
     LE = 261,
     MINUSMINUS = 262,
     NOT = 263,
     NOTEQUAL = 264,
     PLUSPLUS = 265,
     COLONCOLON = 266,
     GRING_CMD = 267,
     INTMAT_CMD = 268,
     PROC_CMD = 269,
     RING_CMD = 270,
     BEGIN_RING = 271,
     IDEAL_CMD = 272,
     MAP_CMD = 273,
     MATRIX_CMD = 274,
     MODUL_CMD = 275,
     NUMBER_CMD = 276,
     POLY_CMD = 277,
     RESOLUTION_CMD = 278,
     VECTOR_CMD = 279,
     BETTI_CMD = 280,
     COEFFS_CMD = 281,
     COEF_CMD = 282,
     CONTRACT_CMD = 283,
     DEGREE_CMD = 284,
     DEG_CMD = 285,
     DIFF_CMD = 286,
     DIM_CMD = 287,
     DIVISION_CMD = 288,
     ELIMINATION_CMD = 289,
     E_CMD = 290,
     FETCH_CMD = 291,
     FREEMODULE_CMD = 292,
     KEEPRING_CMD = 293,
     HILBERT_CMD = 294,
     HOMOG_CMD = 295,
     IMAP_CMD = 296,
     INDEPSET_CMD = 297,
     INTERRED_CMD = 298,
     INTERSECT_CMD = 299,
     JACOB_CMD = 300,
     JET_CMD = 301,
     KBASE_CMD = 302,
     KOSZUL_CMD = 303,
     LEADCOEF_CMD = 304,
     LEADEXP_CMD = 305,
     LEAD_CMD = 306,
     LEADMONOM_CMD = 307,
     LIFTSTD_CMD = 308,
     LIFT_CMD = 309,
     MAXID_CMD = 310,
     MINBASE_CMD = 311,
     MINOR_CMD = 312,
     MINRES_CMD = 313,
     MODULO_CMD = 314,
     MRES_CMD = 315,
     MULTIPLICITY_CMD = 316,
     ORD_CMD = 317,
     PAR_CMD = 318,
     PARDEG_CMD = 319,
     PREIMAGE_CMD = 320,
     QUOTIENT_CMD = 321,
     QHWEIGHT_CMD = 322,
     REDUCE_CMD = 323,
     REGULARITY_CMD = 324,
     RES_CMD = 325,
     SIMPLIFY_CMD = 326,
     SORTVEC_CMD = 327,
     SRES_CMD = 328,
     STD_CMD = 329,
     SUBST_CMD = 330,
     SYZYGY_CMD = 331,
     VAR_CMD = 332,
     VDIM_CMD = 333,
     WEDGE_CMD = 334,
     WEIGHT_CMD = 335,
     VALTVARS = 336,
     VMAXDEG = 337,
     VMAXMULT = 338,
     VNOETHER = 339,
     VMINPOLY = 340,
     END_RING = 341,
     CMD_1 = 342,
     CMD_2 = 343,
     CMD_3 = 344,
     CMD_12 = 345,
     CMD_13 = 346,
     CMD_23 = 347,
     CMD_123 = 348,
     CMD_M = 349,
     ROOT_DECL = 350,
     ROOT_DECL_LIST = 351,
     RING_DECL = 352,
     EXAMPLE_CMD = 353,
     EXPORT_CMD = 354,
     HELP_CMD = 355,
     KILL_CMD = 356,
     LIB_CMD = 357,
     LISTVAR_CMD = 358,
     SETRING_CMD = 359,
     TYPE_CMD = 360,
     STRINGTOK = 361,
     BLOCKTOK = 362,
     INT_CONST = 363,
     UNKNOWN_IDENT = 364,
     RINGVAR = 365,
     PROC_DEF = 366,
     BREAK_CMD = 367,
     CONTINUE_CMD = 368,
     ELSE_CMD = 369,
     EVAL = 370,
     QUOTE = 371,
     FOR_CMD = 372,
     IF_CMD = 373,
     SYS_BREAK = 374,
     WHILE_CMD = 375,
     RETURN = 376,
     PARAMETER = 377,
     SYSVAR = 378,
     UMINUS = 379
   };
#endif
/* Tokens.  */
#define DOTDOT 258
#define EQUAL_EQUAL 259
#define GE 260
#define LE 261
#define MINUSMINUS 262
#define NOT 263
#define NOTEQUAL 264
#define PLUSPLUS 265
#define COLONCOLON 266
#define GRING_CMD 267
#define INTMAT_CMD 268
#define PROC_CMD 269
#define RING_CMD 270
#define BEGIN_RING 271
#define IDEAL_CMD 272
#define MAP_CMD 273
#define MATRIX_CMD 274
#define MODUL_CMD 275
#define NUMBER_CMD 276
#define POLY_CMD 277
#define RESOLUTION_CMD 278
#define VECTOR_CMD 279
#define BETTI_CMD 280
#define COEFFS_CMD 281
#define COEF_CMD 282
#define CONTRACT_CMD 283
#define DEGREE_CMD 284
#define DEG_CMD 285
#define DIFF_CMD 286
#define DIM_CMD 287
#define DIVISION_CMD 288
#define ELIMINATION_CMD 289
#define E_CMD 290
#define FETCH_CMD 291
#define FREEMODULE_CMD 292
#define KEEPRING_CMD 293
#define HILBERT_CMD 294
#define HOMOG_CMD 295
#define IMAP_CMD 296
#define INDEPSET_CMD 297
#define INTERRED_CMD 298
#define INTERSECT_CMD 299
#define JACOB_CMD 300
#define JET_CMD 301
#define KBASE_CMD 302
#define KOSZUL_CMD 303
#define LEADCOEF_CMD 304
#define LEADEXP_CMD 305
#define LEAD_CMD 306
#define LEADMONOM_CMD 307
#define LIFTSTD_CMD 308
#define LIFT_CMD 309
#define MAXID_CMD 310
#define MINBASE_CMD 311
#define MINOR_CMD 312
#define MINRES_CMD 313
#define MODULO_CMD 314
#define MRES_CMD 315
#define MULTIPLICITY_CMD 316
#define ORD_CMD 317
#define PAR_CMD 318
#define PARDEG_CMD 319
#define PREIMAGE_CMD 320
#define QUOTIENT_CMD 321
#define QHWEIGHT_CMD 322
#define REDUCE_CMD 323
#define REGULARITY_CMD 324
#define RES_CMD 325
#define SIMPLIFY_CMD 326
#define SORTVEC_CMD 327
#define SRES_CMD 328
#define STD_CMD 329
#define SUBST_CMD 330
#define SYZYGY_CMD 331
#define VAR_CMD 332
#define VDIM_CMD 333
#define WEDGE_CMD 334
#define WEIGHT_CMD 335
#define VALTVARS 336
#define VMAXDEG 337
#define VMAXMULT 338
#define VNOETHER 339
#define VMINPOLY 340
#define END_RING 341
#define CMD_1 342
#define CMD_2 343
#define CMD_3 344
#define CMD_12 345
#define CMD_13 346
#define CMD_23 347
#define CMD_123 348
#define CMD_M 349
#define ROOT_DECL 350
#define ROOT_DECL_LIST 351
#define RING_DECL 352
#define EXAMPLE_CMD 353
#define EXPORT_CMD 354
#define HELP_CMD 355
#define KILL_CMD 356
#define LIB_CMD 357
#define LISTVAR_CMD 358
#define SETRING_CMD 359
#define TYPE_CMD 360
#define STRINGTOK 361
#define BLOCKTOK 362
#define INT_CONST 363
#define UNKNOWN_IDENT 364
#define RINGVAR 365
#define PROC_DEF 366
#define BREAK_CMD 367
#define CONTINUE_CMD 368
#define ELSE_CMD 369
#define EVAL 370
#define QUOTE 371
#define FOR_CMD 372
#define IF_CMD 373
#define SYS_BREAK 374
#define WHILE_CMD 375
#define RETURN 376
#define PARAMETER 377
#define SYSVAR 378
#define UMINUS 379




/* Copy the first part of user declarations.  */
#line 8 "grammar.y"


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



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 535 "grammar.cc"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   1848

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  142
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  44
/* YYNRULES -- Number of rules.  */
#define YYNRULES  167
/* YYNRULES -- Number of states.  */
#define YYNSTATES  369

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   379

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   135,     2,
     139,   140,   137,   127,   133,   128,     2,   129,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   136,   134,
     125,   124,   126,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   130,     2,   131,   132,     2,   141,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   138
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    12,    15,    17,    19,
      21,    24,    26,    28,    30,    32,    34,    36,    38,    40,
      43,    45,    47,    49,    51,    53,    55,    57,    59,    61,
      64,    66,    68,    72,    76,    81,    85,    87,    89,    91,
      95,    97,    99,   101,   105,   112,   117,   122,   127,   131,
     136,   141,   145,   150,   157,   166,   173,   182,   187,   192,
     199,   204,   211,   220,   229,   233,   238,   247,   252,   261,
     266,   275,   280,   284,   290,   291,   297,   300,   302,   305,
     308,   312,   316,   320,   324,   328,   332,   336,   340,   344,
     348,   351,   354,   357,   360,   362,   366,   368,   370,   373,
     376,   379,   382,   391,   394,   403,   406,   410,   413,   415,
     417,   423,   425,   427,   432,   434,   438,   440,   444,   446,
     447,   452,   456,   459,   463,   466,   469,   473,   478,   483,
     488,   493,   498,   503,   508,   513,   518,   525,   532,   539,
     546,   553,   560,   567,   574,   578,   580,   589,   592,   595,
     597,   599,   602,   605,   607,   613,   616,   622,   624,   626,
     630,   636,   640,   644,   649,   652,   655,   660
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     143,     0,    -1,    -1,   143,   144,    -1,   145,    -1,   147,
     134,    -1,   159,   134,    -1,   185,    -1,   119,    -1,   134,
      -1,     1,   134,    -1,   180,    -1,   181,    -1,   146,    -1,
     182,    -1,   183,    -1,   167,    -1,   169,    -1,   170,    -1,
      98,   107,    -1,   148,    -1,   171,    -1,   172,    -1,   173,
      -1,   184,    -1,   175,    -1,   176,    -1,   178,    -1,   179,
      -1,   156,   150,    -1,   110,    -1,   157,    -1,   149,    11,
     149,    -1,   149,   139,   140,    -1,   149,   139,   150,   140,
      -1,   130,   150,   131,    -1,   108,    -1,   123,    -1,   160,
      -1,   150,   133,   151,    -1,   151,    -1,   155,    -1,   149,
      -1,   139,   150,   140,    -1,   151,   130,   151,   133,   151,
     131,    -1,   151,   130,   151,   131,    -1,    95,   139,   151,
     140,    -1,    96,   139,   150,   140,    -1,    96,   139,   140,
      -1,    97,   139,   151,   140,    -1,   158,   139,   150,   140,
      -1,   158,   139,   140,    -1,    87,   139,   151,   140,    -1,
      88,   139,   151,   133,   151,   140,    -1,    89,   139,   151,
     133,   151,   133,   151,   140,    -1,    92,   139,   151,   133,
     151,   140,    -1,    92,   139,   151,   133,   151,   133,   151,
     140,    -1,    90,   139,   151,   140,    -1,    91,   139,   151,
     140,    -1,    90,   139,   151,   133,   151,   140,    -1,    93,
     139,   151,   140,    -1,    93,   139,   151,   133,   151,   140,
      -1,    91,   139,   151,   133,   151,   133,   151,   140,    -1,
      93,   139,   151,   133,   151,   133,   151,   140,    -1,    94,
     139,   140,    -1,    94,   139,   150,   140,    -1,    19,   139,
     151,   133,   151,   133,   151,   140,    -1,    19,   139,   151,
     140,    -1,    13,   139,   151,   133,   151,   133,   151,   140,
      -1,    13,   139,   151,   140,    -1,    15,   139,   161,   133,
     161,   133,   165,   140,    -1,    15,   139,   151,   140,    -1,
     153,   151,   154,    -1,   153,   151,   124,   151,   154,    -1,
      -1,   115,   139,   152,   151,   140,    -1,   116,   139,    -1,
     140,    -1,   151,    10,    -1,   151,     7,    -1,   151,   127,
     151,    -1,   151,   128,   151,    -1,   151,   129,   151,    -1,
     151,   132,   151,    -1,   151,   125,   151,    -1,   151,   135,
     151,    -1,   151,     9,   151,    -1,   151,     4,   151,    -1,
     151,     3,   151,    -1,   151,   136,   151,    -1,     8,   151,
      -1,   128,   151,    -1,   159,   166,    -1,   150,   124,    -1,
     109,    -1,   141,   151,   141,    -1,    17,    -1,    20,    -1,
      95,   149,    -1,    96,   149,    -1,    97,   149,    -1,   158,
     149,    -1,    19,   149,   130,   151,   131,   130,   151,   131,
      -1,    19,   149,    -1,    13,   149,   130,   151,   131,   130,
     151,   131,    -1,    13,   149,    -1,   159,   133,   149,    -1,
      14,   149,    -1,   106,    -1,   151,    -1,   139,   151,   133,
     150,   140,    -1,   109,    -1,   162,    -1,   162,   139,   150,
     140,    -1,   163,    -1,   163,   133,   164,    -1,   163,    -1,
     139,   164,   140,    -1,   124,    -1,    -1,   125,   160,   168,
     134,    -1,   100,   106,   134,    -1,   100,   134,    -1,    98,
     106,   134,    -1,    99,   150,    -1,   101,   149,    -1,   172,
     133,   149,    -1,   103,   139,    95,   140,    -1,   103,   139,
      96,   140,    -1,   103,   139,    97,   140,    -1,   103,   139,
     158,   140,    -1,   103,   139,    15,   140,    -1,   103,   139,
      19,   140,    -1,   103,   139,    13,   140,    -1,   103,   139,
      14,   140,    -1,   103,   139,   149,   140,    -1,   103,   139,
     149,   133,    95,   140,    -1,   103,   139,   149,   133,    96,
     140,    -1,   103,   139,   149,   133,    97,   140,    -1,   103,
     139,   149,   133,   158,   140,    -1,   103,   139,   149,   133,
      15,   140,    -1,   103,   139,   149,   133,    19,   140,    -1,
     103,   139,   149,   133,    13,   140,    -1,   103,   139,   149,
     133,    14,   140,    -1,   103,   139,   140,    -1,    15,    -1,
     174,   149,   166,   161,   133,   161,   133,   165,    -1,   174,
     149,    -1,   123,   160,    -1,   104,    -1,    38,    -1,   177,
     151,    -1,   105,   151,    -1,   150,    -1,   118,   139,   151,
     140,   107,    -1,   114,   107,    -1,   118,   139,   151,   140,
     112,    -1,   112,    -1,   113,    -1,   120,   106,   107,    -1,
     117,   106,   106,   106,   107,    -1,    14,   157,   107,    -1,
     111,   106,   107,    -1,   111,   106,   106,   107,    -1,   122,
     159,    -1,   122,   151,    -1,   121,   139,   150,   140,    -1,
     121,   139,   140,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   364,   364,   366,   401,   402,   404,   406,   410,   415,
     417,   465,   466,   467,   468,   469,   470,   471,   472,   476,
     479,   480,   481,   482,   483,   484,   485,   486,   487,   490,
     497,   502,   506,   510,   514,   534,   562,   590,   596,   605,
     616,   622,   627,   628,   629,   633,   637,   641,   645,   649,
     653,   657,   661,   665,   669,   673,   677,   681,   685,   689,
     693,   697,   701,   705,   709,   713,   717,   721,   725,   729,
     733,   737,   741,   745,   763,   762,   780,   788,   797,   801,
     805,   809,   813,   817,   821,   825,   829,   833,   837,   841,
     845,   852,   859,   860,   873,   874,   886,   886,   891,   900,
     909,   913,   917,   933,   937,   957,   975,  1000,  1012,  1016,
    1017,  1031,  1039,  1048,  1093,  1094,  1103,  1104,  1110,  1123,
    1122,  1130,  1135,  1142,  1150,  1164,  1184,  1207,  1211,  1215,
    1220,  1224,  1228,  1232,  1236,  1240,  1245,  1254,  1263,  1272,
    1281,  1290,  1299,  1308,  1326,  1333,  1337,  1380,  1391,  1404,
    1404,  1407,  1483,  1489,  1518,  1531,  1548,  1557,  1562,  1570,
    1582,  1601,  1612,  1632,  1656,  1662,  1674,  1681
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "DOTDOT", "EQUAL_EQUAL", "GE", "LE",
  "MINUSMINUS", "NOT", "NOTEQUAL", "PLUSPLUS", "COLONCOLON", "GRING_CMD",
  "INTMAT_CMD", "PROC_CMD", "RING_CMD", "BEGIN_RING", "IDEAL_CMD",
  "MAP_CMD", "MATRIX_CMD", "MODUL_CMD", "NUMBER_CMD", "POLY_CMD",
  "RESOLUTION_CMD", "VECTOR_CMD", "BETTI_CMD", "COEFFS_CMD", "COEF_CMD",
  "CONTRACT_CMD", "DEGREE_CMD", "DEG_CMD", "DIFF_CMD", "DIM_CMD",
  "DIVISION_CMD", "ELIMINATION_CMD", "E_CMD", "FETCH_CMD",
  "FREEMODULE_CMD", "KEEPRING_CMD", "HILBERT_CMD", "HOMOG_CMD", "IMAP_CMD",
  "INDEPSET_CMD", "INTERRED_CMD", "INTERSECT_CMD", "JACOB_CMD", "JET_CMD",
  "KBASE_CMD", "KOSZUL_CMD", "LEADCOEF_CMD", "LEADEXP_CMD", "LEAD_CMD",
  "LEADMONOM_CMD", "LIFTSTD_CMD", "LIFT_CMD", "MAXID_CMD", "MINBASE_CMD",
  "MINOR_CMD", "MINRES_CMD", "MODULO_CMD", "MRES_CMD", "MULTIPLICITY_CMD",
  "ORD_CMD", "PAR_CMD", "PARDEG_CMD", "PREIMAGE_CMD", "QUOTIENT_CMD",
  "QHWEIGHT_CMD", "REDUCE_CMD", "REGULARITY_CMD", "RES_CMD",
  "SIMPLIFY_CMD", "SORTVEC_CMD", "SRES_CMD", "STD_CMD", "SUBST_CMD",
  "SYZYGY_CMD", "VAR_CMD", "VDIM_CMD", "WEDGE_CMD", "WEIGHT_CMD",
  "VALTVARS", "VMAXDEG", "VMAXMULT", "VNOETHER", "VMINPOLY", "END_RING",
  "CMD_1", "CMD_2", "CMD_3", "CMD_12", "CMD_13", "CMD_23", "CMD_123",
  "CMD_M", "ROOT_DECL", "ROOT_DECL_LIST", "RING_DECL", "EXAMPLE_CMD",
  "EXPORT_CMD", "HELP_CMD", "KILL_CMD", "LIB_CMD", "LISTVAR_CMD",
  "SETRING_CMD", "TYPE_CMD", "STRINGTOK", "BLOCKTOK", "INT_CONST",
  "UNKNOWN_IDENT", "RINGVAR", "PROC_DEF", "BREAK_CMD", "CONTINUE_CMD",
  "ELSE_CMD", "EVAL", "QUOTE", "FOR_CMD", "IF_CMD", "SYS_BREAK",
  "WHILE_CMD", "RETURN", "PARAMETER", "SYSVAR", "'='", "'<'", "'>'", "'+'",
  "'-'", "'/'", "'['", "']'", "'^'", "','", "';'", "'&'", "':'", "'*'",
  "UMINUS", "'('", "')'", "'`'", "$accept", "lines", "pprompt", "flowctrl",
  "example_dummy", "command", "assign", "elemexpr", "exprlist", "expr",
  "@1", "quote_start", "quote_end", "expr_arithmetic", "left_value",
  "extendedid", "currring_lists", "declare_ip_variable", "stringexpr",
  "rlist", "ordername", "orderelem", "OrderingList", "ordering", "cmdeq",
  "filecmd", "@2", "helpcmd", "examplecmd", "exportcmd", "killcmd",
  "listcmd", "ringcmd1", "ringcmd", "scriptcmd", "setrings", "setringcmd",
  "typecmd", "ifcmd", "whilecmd", "forcmd", "proccmd", "parametercmd",
  "returncmd", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,    61,    60,    62,    43,    45,    47,
      91,    93,    94,    44,    59,    38,    58,    42,   379,    40,
      41,    96
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   142,   143,   143,   144,   144,   144,   144,   144,   144,
     144,   145,   145,   145,   145,   145,   145,   145,   145,   146,
     147,   147,   147,   147,   147,   147,   147,   147,   147,   148,
     149,   149,   149,   149,   149,   149,   149,   149,   149,   150,
     150,   151,   151,   151,   151,   151,   151,   151,   151,   151,
     151,   151,   151,   151,   151,   151,   151,   151,   151,   151,
     151,   151,   151,   151,   151,   151,   151,   151,   151,   151,
     151,   151,   151,   151,   152,   151,   153,   154,   155,   155,
     155,   155,   155,   155,   155,   155,   155,   155,   155,   155,
     155,   155,   156,   156,   157,   157,   158,   158,   159,   159,
     159,   159,   159,   159,   159,   159,   159,   159,   160,   161,
     161,   162,   163,   163,   164,   164,   165,   165,   166,   168,
     167,   169,   169,   170,   171,   172,   172,   173,   173,   173,
     173,   173,   173,   173,   173,   173,   173,   173,   173,   173,
     173,   173,   173,   173,   173,   174,   175,   175,   176,   177,
     177,   178,   179,   179,   180,   180,   180,   180,   180,   181,
     182,   183,   183,   183,   184,   184,   185,   185
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     2,     2,     1,     1,     1,
       2,     1,     1,     1,     1,     1,     1,     1,     1,     2,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     2,
       1,     1,     3,     3,     4,     3,     1,     1,     1,     3,
       1,     1,     1,     3,     6,     4,     4,     4,     3,     4,
       4,     3,     4,     6,     8,     6,     8,     4,     4,     6,
       4,     6,     8,     8,     3,     4,     8,     4,     8,     4,
       8,     4,     3,     5,     0,     5,     2,     1,     2,     2,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       2,     2,     2,     2,     1,     3,     1,     1,     2,     2,
       2,     2,     8,     2,     8,     2,     3,     2,     1,     1,
       5,     1,     1,     4,     1,     3,     1,     3,     1,     0,
       4,     3,     2,     3,     2,     2,     3,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     6,     6,     6,     6,
       6,     6,     6,     6,     3,     1,     8,     2,     2,     1,
       1,     2,     2,     1,     5,     2,     5,     1,     1,     3,
       5,     3,     3,     4,     2,     2,     4,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,     0,     0,     0,     0,   145,    96,     0,
      97,   150,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   149,     0,
     108,    36,    94,    30,     0,   157,   158,     0,     0,     0,
       0,     0,     8,     0,     0,     0,    37,     0,     0,     0,
       9,     0,     0,     3,     4,    13,     0,    20,    42,   153,
      40,     0,    41,     0,    31,     0,     0,    38,    16,    17,
      18,    21,    22,    23,     0,    25,    26,     0,    27,    28,
      11,    12,    14,    15,    24,     7,    10,     0,     0,     0,
       0,     0,     0,    37,    90,     0,     0,   105,   107,    31,
       0,     0,   103,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    98,     0,    99,     0,   100,     0,    19,   124,
       0,   122,   125,     0,   152,     0,   155,    74,    76,     0,
       0,     0,     0,     0,   165,   164,   148,   119,    91,     0,
       0,     0,     5,     0,     0,    93,     0,     0,     0,    79,
       0,    78,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    29,     0,   101,   118,     0,     6,    92,     0,   147,
     151,     0,     0,   161,     0,   109,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    64,     0,     0,    48,
       0,     0,   123,   121,     0,     0,     0,     0,     0,     0,
       0,   144,     0,     0,     0,   162,     0,     0,     0,   159,
     167,     0,     0,    35,    43,    95,    32,    33,     0,    39,
      88,    87,    86,    84,    80,    81,    82,     0,    83,    85,
      89,     0,    77,    72,    51,     0,   106,   126,     0,     0,
      69,     0,    40,    71,     0,     0,    67,     0,    52,     0,
       0,     0,    57,     0,    58,     0,     0,    60,    65,    46,
      47,    49,   133,   134,   131,   132,   127,   128,   129,     0,
     135,   130,   163,     0,     0,     0,   166,   120,    34,    45,
       0,     0,    50,   109,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    75,   160,   154,   156,     0,
      73,     0,     0,     0,     0,     0,     0,     0,    53,     0,
      59,     0,     0,    55,     0,    61,   142,   143,   140,   141,
     136,   137,   138,   139,    44,     0,     0,     0,   110,   111,
       0,   112,   116,     0,     0,     0,     0,     0,     0,     0,
       0,    68,   104,   114,     0,     0,    70,    66,   102,    54,
      62,    56,    63,   146,     0,   117,     0,   115,   113
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,    53,    54,    55,    56,    57,    58,   140,    60,
     206,    61,   233,    62,    63,    64,    95,    66,    67,   176,
     341,   342,   354,   343,   167,    68,   212,    69,    70,    71,
      72,    73,    74,    75,    76,    77,    78,    79,    80,    81,
      82,    83,    84,    85
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -339
static const yytype_int16 yypact[] =
{
    -339,   277,  -339,  -131,  1219,   632,   215,  -130,  -339,   801,
    -339,  -339,  -125,  -122,   -93,   -66,   -64,   -58,   -50,   -47,
     806,   877,   932,   -21,  1219,   -96,   215,   -39,  -339,  1219,
    -339,  -339,  -339,  -339,   -51,  -339,  -339,   -40,   -37,   -36,
      -5,   -33,  -339,     1,   -26,   883,    14,    14,  1219,  1219,
    -339,  1219,  1219,  -339,  -339,  -339,   -18,  -339,    -4,  -118,
    1633,  1219,  -339,  1219,  -339,   989,   -80,  -339,  -339,  -339,
    -339,  -339,    13,  -339,   215,  -339,  -339,  1219,  -339,  -339,
    -339,  -339,  -339,  -339,  -339,  -339,  -339,   -24,  -130,   -11,
       8,    28,    30,  -339,    61,    39,  1219,    10,    -4,    51,
    1274,  1219,    21,  1219,  1219,  1219,  1219,  1219,  1219,  1219,
     938,  1219,    -4,   995,    -4,  1219,    -4,    47,  -339,    62,
      60,  -339,    -4,  1707,  1633,   -19,  -339,  -339,  -339,    93,
    1219,    97,  1050,   215,  1633,    74,  -339,  -339,    61,   -72,
    -113,     9,  -339,   215,  1107,  -339,  1219,  1219,  1219,  -339,
    1219,  -339,  1219,  1219,  1219,  1219,  1219,  1219,  1219,  1219,
      73,    62,  1162,    -4,  -339,   215,  -339,  -339,   215,    -7,
    1633,    87,  1219,  -339,  1219,   444,    77,   101,  1219,   474,
    1329,  1341,   176,   310,  1369,   340,  -339,  -104,   488,  -339,
     -98,   502,  -339,  -339,    66,    71,    78,    81,    85,    92,
      99,  -339,    49,   100,   135,  -339,  1219,   137,   516,  -339,
    -339,   -90,   113,  -339,  -339,  -339,  -339,  -339,   -84,  1633,
    1645,   464,   464,    27,   228,   228,    61,  1314,    23,  1657,
     228,  1219,  -339,  -339,  -339,   -82,    -4,    -4,  1274,  1219,
    -339,  1384,  1396,  -339,  1274,  1219,  -339,  1413,  -339,  1219,
    1219,  1219,  -339,  1219,  -339,  1219,  1219,  -339,  -339,  -339,
    -339,  -339,  -339,  -339,  -339,  -339,  -339,  -339,  -339,   231,
    -339,  -339,  -339,   553,   145,   -38,  -339,  -339,  -339,  -339,
    1219,   578,  -339,  1633,   120,  1425,   124,  1219,   122,  1475,
     128,   650,  1559,   664,  1571,   405,   419,   119,   121,   125,
     126,   130,   132,   139,   140,  -339,  -339,  -339,  -339,  1583,
    -339,  1274,  1219,  1219,   -76,  -108,  1219,  1219,  -339,  1219,
    -339,  1219,  1219,  -339,  1219,  -339,  -339,  -339,  -339,  -339,
    -339,  -339,  -339,  -339,  -339,   127,   687,  1609,  -339,  -339,
     155,   136,  -339,   141,   712,  1621,   726,   740,   754,   798,
    -108,  -339,  -339,   149,   143,  1219,  -339,  -339,  -339,  -339,
    -339,  -339,  -339,  -339,   155,  -339,   -68,  -339,  -339
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -339,  -339,  -339,  -339,  -339,  -339,  -339,    19,    -1,    18,
    -339,  -339,     3,  -339,  -339,   281,     7,   243,    52,  -233,
    -339,  -338,   -75,   -57,   131,  -339,  -339,  -339,  -339,  -339,
    -339,  -339,  -339,  -339,  -339,  -339,  -339,  -339,  -339,  -339,
    -339,  -339,  -339,  -339
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_int16 yytable[] =
{
      59,   339,   353,    86,   143,   284,   145,   143,    65,   100,
     120,   288,   147,   148,   103,   146,   149,   104,   150,   151,
     146,   143,    94,   119,    97,    98,   353,   214,   102,   146,
     149,   340,   143,   151,   149,   146,   258,   151,   121,   112,
     114,   116,   260,   146,   164,   122,   105,   124,   139,   146,
     276,   146,    65,   165,   166,   125,   278,   146,   282,   213,
     143,   146,   161,   134,   338,   146,   138,   126,   149,   307,
     141,   151,   368,   106,   308,   107,   147,   148,   335,   160,
     149,   108,   150,   151,   163,   117,   118,   204,   205,   109,
     147,   148,   110,   169,   149,   170,   150,   151,   136,   137,
     123,   129,   127,   128,   147,   148,   130,   131,   149,   187,
     150,   151,   190,   132,   171,    96,   142,   164,   175,   177,
      30,   179,   180,   181,   182,   183,   184,   185,   101,   188,
     203,   211,   144,   191,   152,   144,   153,   154,   155,   156,
     172,   157,   202,   218,   158,   159,   168,   111,   208,   144,
     215,   178,    98,   156,   153,   154,   155,   156,   173,   157,
     144,   235,   216,   159,   219,   220,   221,   113,   222,   115,
     223,   224,   225,   226,   227,   228,   229,   230,   162,   147,
     148,   192,   269,   149,   236,   150,   151,   237,   144,   270,
     241,   156,   242,   157,   193,   146,   247,   231,   152,   207,
     153,   154,   155,   156,   209,   157,   262,   165,   158,   159,
     244,   263,   152,   232,   153,   154,   155,   156,   264,   157,
     239,   265,   158,   159,   273,   266,   152,   240,   153,   154,
     155,   156,   267,   157,   245,   149,   158,   159,   151,   268,
     271,   246,   272,   274,   297,   298,   299,   277,     8,   281,
     300,    10,   306,   311,   313,   315,   283,   285,   317,   326,
     350,   327,   283,   289,   339,   328,   329,   291,   292,   293,
     330,   294,   331,   295,   296,   355,   304,     2,     3,   332,
     333,   356,   364,   365,   310,     4,   314,    99,   135,   367,
       5,     6,     7,   363,     8,     0,     9,    10,   309,     0,
     238,   152,     0,   153,   154,   155,   156,     0,   157,   251,
       0,   158,   159,   147,   148,    11,   252,   149,     0,   150,
     151,    30,     0,    31,    32,    33,   301,   302,   303,   283,
     336,   337,     0,     0,   344,   345,     0,   346,    93,   347,
     348,     0,   349,   147,   148,    49,     0,   149,     0,   150,
     151,     0,     0,     0,   366,     0,    52,   155,   156,     0,
     157,     0,     0,     0,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,     0,
      27,    28,    29,    30,     0,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,     0,    47,     0,     0,    48,     0,    49,   147,   148,
       0,    50,   149,     0,   150,   151,    51,     0,    52,     0,
       0,     0,   147,   148,     0,     0,   149,     0,   150,   151,
       0,     0,     0,     0,     0,   152,     0,   153,   154,   155,
     156,     0,   157,   253,     0,   158,   159,   147,   148,     0,
     254,   149,     0,   150,   151,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   152,     0,   153,   154,   155,
     156,   149,   157,   256,   151,   158,   159,   147,   148,     0,
     257,   149,     0,   150,   151,     0,     0,     0,     0,     0,
       0,   147,   148,     0,     0,   149,     0,   150,   151,     0,
       0,     0,     0,     0,     0,   147,   148,     0,     0,   149,
       0,   150,   151,     0,     0,     0,     0,     0,     0,   147,
     148,     0,     0,   149,     0,   150,   151,     0,     0,     0,
     152,     0,   153,   154,   155,   156,     0,   157,   322,     0,
     158,   159,     0,     0,   152,   323,   153,   154,   155,   156,
       0,   157,   324,     0,   158,   159,   147,   148,     0,   325,
     149,     0,   150,   151,     0,     0,     0,     0,     0,   152,
       0,   153,   154,   155,   156,     0,   157,     0,     0,   158,
     159,   147,   148,     0,   243,   149,     0,   150,   151,   152,
       0,   153,   154,   155,   156,     0,   157,     0,     0,   152,
     159,   153,   154,   155,   156,     0,   157,     0,     0,   158,
     159,     0,     0,   152,   248,   153,   154,   155,   156,     0,
     157,     0,     0,   158,   159,     0,     0,   152,   259,   153,
     154,   155,   156,     0,   157,     0,     0,   158,   159,     0,
       0,   152,   261,   153,   154,   155,   156,     0,   157,     0,
       0,   158,   159,   147,   148,     0,   275,   149,     0,   150,
     151,     0,     0,     0,     0,     0,     0,   147,   148,     0,
       0,   149,     0,   150,   151,     0,     0,     0,   152,     0,
     153,   154,   155,   156,     0,   157,     0,     0,   158,   159,
     147,   148,     0,   305,   149,     0,   150,   151,     0,     0,
       0,     0,     0,   152,     0,   153,   154,   155,   156,     0,
     157,     0,     0,   158,   159,   147,   148,     0,   232,   149,
       0,   150,   151,     0,     0,     0,     0,     0,     0,   147,
     148,     0,     0,   149,     0,   150,   151,     0,    30,     0,
      31,    32,    33,   147,   148,     0,     0,   149,     0,   150,
     151,     0,     0,     0,     0,    93,     0,   147,   148,     0,
       0,   149,    49,   150,   151,     0,     0,     0,     0,     0,
       0,    96,     0,    52,     0,   152,     0,   153,   154,   155,
     156,     0,   157,     0,     0,   158,   159,     0,     0,   152,
     318,   153,   154,   155,   156,     0,   157,     0,     0,   158,
     159,   147,   148,     0,   320,   149,     0,   150,   151,     0,
       0,     0,   152,     0,   153,   154,   155,   156,     0,   157,
       0,     0,   158,   159,     0,     0,     0,   351,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   152,     0,   153,
     154,   155,   156,     0,   157,     0,     0,   158,   159,     0,
       0,   152,   357,   153,   154,   155,   156,     0,   157,     0,
       0,   158,   159,     0,     0,   152,   359,   153,   154,   155,
     156,     0,   157,     0,     0,   158,   159,     0,     0,   152,
     360,   153,   154,   155,   156,     0,   157,     0,     0,   158,
     159,     4,     0,     0,   361,     0,     5,   133,    88,     0,
       8,     0,     9,    10,     0,     0,     0,    30,     0,    31,
      32,    33,    30,     0,    31,    32,    33,     0,     0,     0,
       0,     0,     0,   152,    93,   153,   154,   155,   156,    93,
     157,    49,     0,   158,   159,     0,    49,     0,   362,     0,
     101,     0,    52,     0,     0,   111,     4,    52,     0,     0,
       0,    87,     0,    88,     0,     8,     0,    89,    10,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,     0,     0,    30,     0,    31,    32,    33,     0,    30,
       0,    31,    32,    33,     0,     0,     0,     0,    38,    39,
      93,     0,     0,     4,     0,     0,    93,    49,    87,     0,
      88,    48,     8,    49,    89,    10,   113,     0,    52,     0,
       0,     0,    51,     0,    52,    12,    13,    14,    15,    16,
      17,    18,    19,    90,    91,    92,     0,     0,    30,     0,
      31,    32,    33,     0,    30,     0,    31,    32,    33,     0,
       0,     0,     0,    38,    39,    93,     0,     0,     4,     0,
       0,    93,    49,    87,     0,    88,    48,     8,    49,    89,
      10,   115,     0,    52,     0,     0,     0,    51,   186,    52,
       0,     0,    12,    13,    14,    15,    16,    17,    18,    19,
      90,    91,    92,     0,     0,    30,     0,    31,    32,    33,
       0,    30,     0,    31,    32,    33,     0,     0,     0,     0,
      38,    39,    93,     0,     0,     4,     0,     0,    93,    49,
      87,     0,    88,    48,     8,    49,    89,    10,   162,     0,
      52,     0,     0,     0,    51,   189,    52,    12,    13,    14,
      15,    16,    17,    18,    19,    90,    91,    92,     0,     0,
       0,     0,     0,     0,     0,     0,    30,     0,    31,    32,
      33,     0,     0,     0,     0,    38,    39,     0,     0,     0,
       4,     0,     0,    93,     0,    87,     0,    88,    48,     8,
      49,    89,    10,     0,     0,     0,     0,     0,     0,    51,
     210,    52,     0,     0,    12,    13,    14,    15,    16,    17,
      18,    19,    90,    91,    92,     0,     0,     0,     0,     0,
       0,     0,     0,    30,     0,    31,    32,    33,     0,     0,
       0,     0,    38,    39,     0,     0,     0,     4,     0,     0,
      93,     0,    87,     0,    88,    48,     8,    49,    89,    10,
       0,     0,     0,     0,     0,     0,    51,   217,    52,    12,
      13,    14,    15,    16,    17,    18,    19,    90,    91,    92,
       0,     0,     0,     0,     0,     0,     0,     0,    30,     0,
      31,    32,    33,     0,     0,     0,     0,    38,    39,     0,
       0,     0,     4,     0,     0,    93,     0,    87,     0,    88,
      48,     8,    49,    89,    10,     0,     0,     0,     0,     0,
       0,    51,   234,    52,     0,     0,    12,    13,    14,    15,
      16,    17,    18,    19,    90,    91,    92,   147,   148,     0,
       0,   149,     0,   150,   151,    30,     0,    31,    32,    33,
       0,     0,   147,   148,    38,    39,   149,     0,   150,   151,
       0,     0,    93,     0,   147,   148,     0,    48,   149,    49,
     150,   151,     0,     0,     0,     0,     0,     0,    51,     0,
      52,    12,    13,    14,    15,    16,    17,    18,    19,    90,
      91,    92,   147,   148,     0,     0,   149,     0,   150,   151,
      30,     0,    31,    32,    33,     0,     0,   147,   148,    38,
      39,   149,     0,   150,   151,     0,     0,    93,     0,   147,
     148,     0,    48,   149,    49,   150,   151,     0,     0,     0,
       0,     0,     0,   174,     0,    52,   147,   148,     0,     0,
     149,     0,   150,   151,     0,     0,     0,     0,   147,   148,
       0,     0,   149,     0,   150,   151,     0,     0,     0,   152,
       0,   153,   154,   155,   156,   279,   157,   280,     0,   158,
     159,     0,     0,     0,   152,     0,   153,   154,   155,   156,
       0,   157,   249,     0,   158,   159,   152,     0,   153,   154,
     155,   156,     0,   157,   250,     0,   158,   159,   147,   148,
       0,     0,   149,     0,   150,   151,     0,     0,     0,     0,
       0,     0,     0,     0,   152,     0,   153,   154,   155,   156,
       0,   157,   255,     0,   158,   159,     0,     0,     0,   152,
       0,   153,   154,   155,   156,   286,   157,     0,     0,   158,
     159,   152,     0,   153,   154,   155,   156,     0,   157,   287,
       0,   158,   159,     0,     0,     0,     0,     0,   152,     0,
     153,   154,   155,   156,   290,   157,     0,     0,   158,   159,
     152,     0,   153,   154,   155,   156,     0,   157,   312,     0,
     158,   159,   147,   148,     0,     0,   149,     0,   150,   151,
       0,     0,     0,     0,   147,   148,     0,     0,   149,     0,
     150,   151,     0,     0,     0,     0,   147,   148,     0,     0,
     149,     0,   150,   151,     0,     0,     0,     0,     0,     0,
     152,     0,   153,   154,   155,   156,     0,   157,   316,     0,
     158,   159,   147,   148,     0,     0,   149,     0,   150,   151,
       0,     0,     0,     0,   147,   148,     0,     0,   149,     0,
     150,   151,     0,     0,     0,     0,   147,   148,     0,     0,
     149,     0,   150,   151,     0,     0,     0,     0,    -1,   148,
       0,     0,   149,     0,   150,   151,     0,     0,     0,     0,
       0,   148,     0,     0,   149,     0,   150,   151,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   152,     0,   153,   154,   155,   156,
       0,   157,   319,     0,   158,   159,   152,     0,   153,   154,
     155,   156,     0,   157,   321,     0,   158,   159,   152,     0,
     153,   154,   155,   156,   334,   157,     0,     0,   158,   159,
     194,   195,   196,     0,     8,     0,   197,    10,     0,     0,
       0,     0,     0,     0,   152,     0,   153,   154,   155,   156,
     352,   157,     0,     0,   158,   159,   152,     0,   153,   154,
     155,   156,   358,   157,     0,     0,   158,   159,   152,     0,
     153,   154,   155,   156,     0,   157,     0,     0,   158,   159,
     152,     0,   153,   154,   155,   156,     0,   157,     0,     0,
     158,   159,   152,     0,   153,   154,   155,   156,     0,   157,
       0,     0,     0,   159,     0,     0,     0,     0,     0,     0,
       0,     0,   198,   199,   200,     0,     0,     0,     0,     0,
       0,     0,     0,    30,     0,    31,    32,    33,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      93,     0,     0,     0,     0,     0,     0,    49,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   201,    52
};

static const yytype_int16 yycheck[] =
{
       1,   109,   340,   134,    11,   238,   124,    11,     1,   139,
     106,   244,     3,     4,   139,   133,     7,   139,     9,    10,
     133,    11,     4,    24,     5,     6,   364,   140,     9,   133,
       7,   139,    11,    10,     7,   133,   140,    10,   134,    20,
      21,    22,   140,   133,   124,    26,   139,    29,    49,   133,
     140,   133,    45,   133,   134,   106,   140,   133,   140,   131,
      11,   133,    63,    45,   140,   133,    48,   107,     7,   107,
      52,    10,   140,   139,   112,   139,     3,     4,   311,    61,
       7,   139,     9,    10,    65,   106,   107,   106,   107,   139,
       3,     4,   139,    74,     7,    77,     9,    10,    46,    47,
     139,   106,   139,   139,     3,     4,   139,   106,     7,   110,
       9,    10,   113,   139,    96,   139,   134,   124,   100,   101,
     106,   103,   104,   105,   106,   107,   108,   109,   139,   111,
     123,   132,   139,   115,   125,   139,   127,   128,   129,   130,
     130,   132,   123,   144,   135,   136,   133,   139,   130,   139,
     141,   130,   133,   130,   127,   128,   129,   130,   107,   132,
     139,   162,   143,   136,   146,   147,   148,   139,   150,   139,
     152,   153,   154,   155,   156,   157,   158,   159,   139,     3,
       4,   134,   133,     7,   165,     9,    10,   168,   139,   140,
     172,   130,   174,   132,   134,   133,   178,   124,   125,   106,
     127,   128,   129,   130,   107,   132,   140,   133,   135,   136,
     133,   140,   125,   140,   127,   128,   129,   130,   140,   132,
     133,   140,   135,   136,   206,   140,   125,   140,   127,   128,
     129,   130,   140,   132,   133,     7,   135,   136,    10,   140,
     140,   140,   107,   106,    13,    14,    15,   134,    17,   231,
      19,    20,   107,   133,   130,   133,   238,   239,   130,   140,
     133,   140,   244,   245,   109,   140,   140,   249,   250,   251,
     140,   253,   140,   255,   256,   139,   269,     0,     1,   140,
     140,   140,   133,   140,   281,     8,   287,     6,    45,   364,
      13,    14,    15,   350,    17,    -1,    19,    20,   280,    -1,
     169,   125,    -1,   127,   128,   129,   130,    -1,   132,   133,
      -1,   135,   136,     3,     4,    38,   140,     7,    -1,     9,
      10,   106,    -1,   108,   109,   110,    95,    96,    97,   311,
     312,   313,    -1,    -1,   316,   317,    -1,   319,   123,   321,
     322,    -1,   324,     3,     4,   130,    -1,     7,    -1,     9,
      10,    -1,    -1,    -1,   355,    -1,   141,   129,   130,    -1,
     132,    -1,    -1,    -1,    87,    88,    89,    90,    91,    92,
      93,    94,    95,    96,    97,    98,    99,   100,   101,    -1,
     103,   104,   105,   106,    -1,   108,   109,   110,   111,   112,
     113,   114,   115,   116,   117,   118,   119,   120,   121,   122,
     123,    -1,   125,    -1,    -1,   128,    -1,   130,     3,     4,
      -1,   134,     7,    -1,     9,    10,   139,    -1,   141,    -1,
      -1,    -1,     3,     4,    -1,    -1,     7,    -1,     9,    10,
      -1,    -1,    -1,    -1,    -1,   125,    -1,   127,   128,   129,
     130,    -1,   132,   133,    -1,   135,   136,     3,     4,    -1,
     140,     7,    -1,     9,    10,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   125,    -1,   127,   128,   129,
     130,     7,   132,   133,    10,   135,   136,     3,     4,    -1,
     140,     7,    -1,     9,    10,    -1,    -1,    -1,    -1,    -1,
      -1,     3,     4,    -1,    -1,     7,    -1,     9,    10,    -1,
      -1,    -1,    -1,    -1,    -1,     3,     4,    -1,    -1,     7,
      -1,     9,    10,    -1,    -1,    -1,    -1,    -1,    -1,     3,
       4,    -1,    -1,     7,    -1,     9,    10,    -1,    -1,    -1,
     125,    -1,   127,   128,   129,   130,    -1,   132,   133,    -1,
     135,   136,    -1,    -1,   125,   140,   127,   128,   129,   130,
      -1,   132,   133,    -1,   135,   136,     3,     4,    -1,   140,
       7,    -1,     9,    10,    -1,    -1,    -1,    -1,    -1,   125,
      -1,   127,   128,   129,   130,    -1,   132,    -1,    -1,   135,
     136,     3,     4,    -1,   140,     7,    -1,     9,    10,   125,
      -1,   127,   128,   129,   130,    -1,   132,    -1,    -1,   125,
     136,   127,   128,   129,   130,    -1,   132,    -1,    -1,   135,
     136,    -1,    -1,   125,   140,   127,   128,   129,   130,    -1,
     132,    -1,    -1,   135,   136,    -1,    -1,   125,   140,   127,
     128,   129,   130,    -1,   132,    -1,    -1,   135,   136,    -1,
      -1,   125,   140,   127,   128,   129,   130,    -1,   132,    -1,
      -1,   135,   136,     3,     4,    -1,   140,     7,    -1,     9,
      10,    -1,    -1,    -1,    -1,    -1,    -1,     3,     4,    -1,
      -1,     7,    -1,     9,    10,    -1,    -1,    -1,   125,    -1,
     127,   128,   129,   130,    -1,   132,    -1,    -1,   135,   136,
       3,     4,    -1,   140,     7,    -1,     9,    10,    -1,    -1,
      -1,    -1,    -1,   125,    -1,   127,   128,   129,   130,    -1,
     132,    -1,    -1,   135,   136,     3,     4,    -1,   140,     7,
      -1,     9,    10,    -1,    -1,    -1,    -1,    -1,    -1,     3,
       4,    -1,    -1,     7,    -1,     9,    10,    -1,   106,    -1,
     108,   109,   110,     3,     4,    -1,    -1,     7,    -1,     9,
      10,    -1,    -1,    -1,    -1,   123,    -1,     3,     4,    -1,
      -1,     7,   130,     9,    10,    -1,    -1,    -1,    -1,    -1,
      -1,   139,    -1,   141,    -1,   125,    -1,   127,   128,   129,
     130,    -1,   132,    -1,    -1,   135,   136,    -1,    -1,   125,
     140,   127,   128,   129,   130,    -1,   132,    -1,    -1,   135,
     136,     3,     4,    -1,   140,     7,    -1,     9,    10,    -1,
      -1,    -1,   125,    -1,   127,   128,   129,   130,    -1,   132,
      -1,    -1,   135,   136,    -1,    -1,    -1,   140,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   125,    -1,   127,
     128,   129,   130,    -1,   132,    -1,    -1,   135,   136,    -1,
      -1,   125,   140,   127,   128,   129,   130,    -1,   132,    -1,
      -1,   135,   136,    -1,    -1,   125,   140,   127,   128,   129,
     130,    -1,   132,    -1,    -1,   135,   136,    -1,    -1,   125,
     140,   127,   128,   129,   130,    -1,   132,    -1,    -1,   135,
     136,     8,    -1,    -1,   140,    -1,    13,    14,    15,    -1,
      17,    -1,    19,    20,    -1,    -1,    -1,   106,    -1,   108,
     109,   110,   106,    -1,   108,   109,   110,    -1,    -1,    -1,
      -1,    -1,    -1,   125,   123,   127,   128,   129,   130,   123,
     132,   130,    -1,   135,   136,    -1,   130,    -1,   140,    -1,
     139,    -1,   141,    -1,    -1,   139,     8,   141,    -1,    -1,
      -1,    13,    -1,    15,    -1,    17,    -1,    19,    20,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      87,    88,    89,    90,    91,    92,    93,    94,    95,    96,
      97,    -1,    -1,   106,    -1,   108,   109,   110,    -1,   106,
      -1,   108,   109,   110,    -1,    -1,    -1,    -1,   115,   116,
     123,    -1,    -1,     8,    -1,    -1,   123,   130,    13,    -1,
      15,   128,    17,   130,    19,    20,   139,    -1,   141,    -1,
      -1,    -1,   139,    -1,   141,    87,    88,    89,    90,    91,
      92,    93,    94,    95,    96,    97,    -1,    -1,   106,    -1,
     108,   109,   110,    -1,   106,    -1,   108,   109,   110,    -1,
      -1,    -1,    -1,   115,   116,   123,    -1,    -1,     8,    -1,
      -1,   123,   130,    13,    -1,    15,   128,    17,   130,    19,
      20,   139,    -1,   141,    -1,    -1,    -1,   139,   140,   141,
      -1,    -1,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    -1,    -1,   106,    -1,   108,   109,   110,
      -1,   106,    -1,   108,   109,   110,    -1,    -1,    -1,    -1,
     115,   116,   123,    -1,    -1,     8,    -1,    -1,   123,   130,
      13,    -1,    15,   128,    17,   130,    19,    20,   139,    -1,
     141,    -1,    -1,    -1,   139,   140,   141,    87,    88,    89,
      90,    91,    92,    93,    94,    95,    96,    97,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   106,    -1,   108,   109,
     110,    -1,    -1,    -1,    -1,   115,   116,    -1,    -1,    -1,
       8,    -1,    -1,   123,    -1,    13,    -1,    15,   128,    17,
     130,    19,    20,    -1,    -1,    -1,    -1,    -1,    -1,   139,
     140,   141,    -1,    -1,    87,    88,    89,    90,    91,    92,
      93,    94,    95,    96,    97,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   106,    -1,   108,   109,   110,    -1,    -1,
      -1,    -1,   115,   116,    -1,    -1,    -1,     8,    -1,    -1,
     123,    -1,    13,    -1,    15,   128,    17,   130,    19,    20,
      -1,    -1,    -1,    -1,    -1,    -1,   139,   140,   141,    87,
      88,    89,    90,    91,    92,    93,    94,    95,    96,    97,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   106,    -1,
     108,   109,   110,    -1,    -1,    -1,    -1,   115,   116,    -1,
      -1,    -1,     8,    -1,    -1,   123,    -1,    13,    -1,    15,
     128,    17,   130,    19,    20,    -1,    -1,    -1,    -1,    -1,
      -1,   139,   140,   141,    -1,    -1,    87,    88,    89,    90,
      91,    92,    93,    94,    95,    96,    97,     3,     4,    -1,
      -1,     7,    -1,     9,    10,   106,    -1,   108,   109,   110,
      -1,    -1,     3,     4,   115,   116,     7,    -1,     9,    10,
      -1,    -1,   123,    -1,     3,     4,    -1,   128,     7,   130,
       9,    10,    -1,    -1,    -1,    -1,    -1,    -1,   139,    -1,
     141,    87,    88,    89,    90,    91,    92,    93,    94,    95,
      96,    97,     3,     4,    -1,    -1,     7,    -1,     9,    10,
     106,    -1,   108,   109,   110,    -1,    -1,     3,     4,   115,
     116,     7,    -1,     9,    10,    -1,    -1,   123,    -1,     3,
       4,    -1,   128,     7,   130,     9,    10,    -1,    -1,    -1,
      -1,    -1,    -1,   139,    -1,   141,     3,     4,    -1,    -1,
       7,    -1,     9,    10,    -1,    -1,    -1,    -1,     3,     4,
      -1,    -1,     7,    -1,     9,    10,    -1,    -1,    -1,   125,
      -1,   127,   128,   129,   130,   131,   132,   133,    -1,   135,
     136,    -1,    -1,    -1,   125,    -1,   127,   128,   129,   130,
      -1,   132,   133,    -1,   135,   136,   125,    -1,   127,   128,
     129,   130,    -1,   132,   133,    -1,   135,   136,     3,     4,
      -1,    -1,     7,    -1,     9,    10,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   125,    -1,   127,   128,   129,   130,
      -1,   132,   133,    -1,   135,   136,    -1,    -1,    -1,   125,
      -1,   127,   128,   129,   130,   131,   132,    -1,    -1,   135,
     136,   125,    -1,   127,   128,   129,   130,    -1,   132,   133,
      -1,   135,   136,    -1,    -1,    -1,    -1,    -1,   125,    -1,
     127,   128,   129,   130,   131,   132,    -1,    -1,   135,   136,
     125,    -1,   127,   128,   129,   130,    -1,   132,   133,    -1,
     135,   136,     3,     4,    -1,    -1,     7,    -1,     9,    10,
      -1,    -1,    -1,    -1,     3,     4,    -1,    -1,     7,    -1,
       9,    10,    -1,    -1,    -1,    -1,     3,     4,    -1,    -1,
       7,    -1,     9,    10,    -1,    -1,    -1,    -1,    -1,    -1,
     125,    -1,   127,   128,   129,   130,    -1,   132,   133,    -1,
     135,   136,     3,     4,    -1,    -1,     7,    -1,     9,    10,
      -1,    -1,    -1,    -1,     3,     4,    -1,    -1,     7,    -1,
       9,    10,    -1,    -1,    -1,    -1,     3,     4,    -1,    -1,
       7,    -1,     9,    10,    -1,    -1,    -1,    -1,     3,     4,
      -1,    -1,     7,    -1,     9,    10,    -1,    -1,    -1,    -1,
      -1,     4,    -1,    -1,     7,    -1,     9,    10,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   125,    -1,   127,   128,   129,   130,
      -1,   132,   133,    -1,   135,   136,   125,    -1,   127,   128,
     129,   130,    -1,   132,   133,    -1,   135,   136,   125,    -1,
     127,   128,   129,   130,   131,   132,    -1,    -1,   135,   136,
      13,    14,    15,    -1,    17,    -1,    19,    20,    -1,    -1,
      -1,    -1,    -1,    -1,   125,    -1,   127,   128,   129,   130,
     131,   132,    -1,    -1,   135,   136,   125,    -1,   127,   128,
     129,   130,   131,   132,    -1,    -1,   135,   136,   125,    -1,
     127,   128,   129,   130,    -1,   132,    -1,    -1,   135,   136,
     125,    -1,   127,   128,   129,   130,    -1,   132,    -1,    -1,
     135,   136,   125,    -1,   127,   128,   129,   130,    -1,   132,
      -1,    -1,    -1,   136,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    95,    96,    97,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   106,    -1,   108,   109,   110,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     123,    -1,    -1,    -1,    -1,    -1,    -1,   130,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   140,   141
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,   143,     0,     1,     8,    13,    14,    15,    17,    19,
      20,    38,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   103,   104,   105,
     106,   108,   109,   110,   111,   112,   113,   114,   115,   116,
     117,   118,   119,   120,   121,   122,   123,   125,   128,   130,
     134,   139,   141,   144,   145,   146,   147,   148,   149,   150,
     151,   153,   155,   156,   157,   158,   159,   160,   167,   169,
     170,   171,   172,   173,   174,   175,   176,   177,   178,   179,
     180,   181,   182,   183,   184,   185,   134,    13,    15,    19,
      95,    96,    97,   123,   151,   158,   139,   149,   149,   157,
     139,   139,   149,   139,   139,   139,   139,   139,   139,   139,
     139,   139,   149,   139,   149,   139,   149,   106,   107,   150,
     106,   134,   149,   139,   151,   106,   107,   139,   139,   106,
     139,   106,   139,    14,   151,   159,   160,   160,   151,   150,
     150,   151,   134,    11,   139,   124,   133,     3,     4,     7,
       9,    10,   125,   127,   128,   129,   130,   132,   135,   136,
     151,   150,   139,   149,   124,   133,   134,   166,   133,   149,
     151,   151,   130,   107,   139,   151,   161,   151,   130,   151,
     151,   151,   151,   151,   151,   151,   140,   150,   151,   140,
     150,   151,   134,   134,    13,    14,    15,    19,    95,    96,
      97,   140,   149,   158,   106,   107,   152,   106,   151,   107,
     140,   150,   168,   131,   140,   141,   149,   140,   150,   151,
     151,   151,   151,   151,   151,   151,   151,   151,   151,   151,
     151,   124,   140,   154,   140,   150,   149,   149,   166,   133,
     140,   151,   151,   140,   133,   133,   140,   151,   140,   133,
     133,   133,   140,   133,   140,   133,   133,   140,   140,   140,
     140,   140,   140,   140,   140,   140,   140,   140,   140,   133,
     140,   140,   107,   151,   106,   140,   140,   134,   140,   131,
     133,   151,   140,   151,   161,   151,   131,   133,   161,   151,
     131,   151,   151,   151,   151,   151,   151,    13,    14,    15,
      19,    95,    96,    97,   158,   140,   107,   107,   112,   151,
     154,   133,   133,   130,   150,   133,   133,   130,   140,   133,
     140,   133,   133,   140,   133,   140,   140,   140,   140,   140,
     140,   140,   140,   140,   131,   161,   151,   151,   140,   109,
     139,   162,   163,   165,   151,   151,   151,   151,   151,   151,
     133,   140,   131,   163,   164,   139,   140,   140,   131,   140,
     140,   140,   140,   165,   133,   140,   150,   164,   140
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (&yylval, YYLEX_PARAM)
#else
# define YYLEX yylex (&yylval)
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */






/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  /* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;

  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 3:
#line 367 "grammar.y"
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
          ;}
    break;

  case 5:
#line 403 "grammar.y"
    {currentVoice->ifsw=0;;}
    break;

  case 6:
#line 405 "grammar.y"
    { (yyvsp[(1) - (2)].lv).CleanUp(); currentVoice->ifsw=0;;}
    break;

  case 7:
#line 407 "grammar.y"
    {
            YYACCEPT;
          ;}
    break;

  case 8:
#line 411 "grammar.y"
    {
            currentVoice->ifsw=0;
            iiDebug();
          ;}
    break;

  case 9:
#line 416 "grammar.y"
    {currentVoice->ifsw=0;;}
    break;

  case 10:
#line 418 "grammar.y"
    {
            #ifdef SIQ
            siq=0;
            #endif
            yyInRingConstruction = FALSE;
            currentVoice->ifsw=0;
            if (inerror)
            {
              if ((inerror!=3) && ((yyvsp[(1) - (2)]).i<UMINUS) && ((yyvsp[(1) - (2)]).i>' '))
              {
                // 1: yyerror called
                // 2: scanner put actual string
                // 3: error rule put token+\n
                inerror=3;
                Print(" error at token `%s`\n",iiTwoOps((yyvsp[(1) - (2)]).i));
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
          ;}
    break;

  case 18:
#line 473 "grammar.y"
    {if (currentVoice!=NULL) currentVoice->ifsw=0;;}
    break;

  case 19:
#line 476 "grammar.y"
    { omFree((ADDRESS)(yyvsp[(2) - (2)].name)); ;}
    break;

  case 29:
#line 491 "grammar.y"
    {
            if(iiAssign(&(yyvsp[(1) - (2)].lv),&(yyvsp[(2) - (2)].lv))) YYERROR;
          ;}
    break;

  case 30:
#line 498 "grammar.y"
    {
            if (currRing==NULL) MYYERROR("no ring active");
            syMake(&(yyval.lv),omStrDup((yyvsp[(1) - (1)].name)));
          ;}
    break;

  case 31:
#line 503 "grammar.y"
    {
            syMake(&(yyval.lv),(yyvsp[(1) - (1)].name));
          ;}
    break;

  case 32:
#line 507 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv), &(yyvsp[(1) - (3)].lv), COLONCOLON, &(yyvsp[(3) - (3)].lv))) YYERROR;
          ;}
    break;

  case 33:
#line 511 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(1) - (3)].lv),'(')) YYERROR;
          ;}
    break;

  case 34:
#line 515 "grammar.y"
    {
            if ((yyvsp[(1) - (4)].lv).rtyp==LIB_CMD)
            {
              if(iiExprArith1(&(yyval.lv),&(yyvsp[(3) - (4)].lv),LIB_CMD)) YYERROR;
            }
            else
            {
              if ((yyvsp[(1) - (4)].lv).Typ()==UNKNOWN)
              {
                if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (4)].lv),'(',&(yyvsp[(3) - (4)].lv))) YYERROR;
              }
              else
              {
                (yyvsp[(1) - (4)].lv).next=(leftv)omAllocBin(sleftv_bin);
                memcpy((yyvsp[(1) - (4)].lv).next,&(yyvsp[(3) - (4)].lv),sizeof(sleftv));
                if(iiExprArithM(&(yyval.lv),&(yyvsp[(1) - (4)].lv),'(')) YYERROR;
              }
            }
          ;}
    break;

  case 35:
#line 535 "grammar.y"
    {
            if (currRingHdl==NULL) MYYERROR("no ring active");
            int j = 0;
            memset(&(yyval.lv),0,sizeof(sleftv));
            (yyval.lv).rtyp=VECTOR_CMD;
            leftv v = &(yyvsp[(2) - (3)].lv);
            while (v!=NULL)
            {
              int i,t;
              sleftv tmp;
              memset(&tmp,0,sizeof(tmp));
              i=iiTestConvert((t=v->Typ()),POLY_CMD);
              if((i==0) || (iiConvert(t /*v->Typ()*/,POLY_CMD,i,v,&tmp)))
              {
                pDelete((poly *)&(yyval.lv).data);
                (yyvsp[(2) - (3)].lv).CleanUp();
                MYYERROR("expected '[poly,...'");
              }
              poly p = (poly)tmp.CopyD(POLY_CMD);
              pSetCompP(p,++j);
              (yyval.lv).data = (void *)pAdd((poly)(yyval.lv).data,p);
              v->next=tmp.next;tmp.next=NULL;
              tmp.CleanUp();
              v=v->next;
            }
            (yyvsp[(2) - (3)].lv).CleanUp();
          ;}
    break;

  case 36:
#line 563 "grammar.y"
    {
            memset(&(yyval.lv),0,sizeof((yyval.lv)));
            int i = atoi((yyvsp[(1) - (1)].name));
            /*remember not to omFree($1)
            *because it is a part of the scanner buffer*/
            (yyval.lv).rtyp  = INT_CMD;
            (yyval.lv).data = (void *)(long)i;

            /* check: out of range input */
            int l = strlen((yyvsp[(1) - (1)].name))+2;
            if (l >= MAX_INT_LEN)
            {
              char tmp[MAX_INT_LEN+5];
              sprintf(tmp,"%d",i);
              if (strcmp(tmp,(yyvsp[(1) - (1)].name))!=0)
              {
                if (currRing==NULL)
                {
                  Werror("`%s` greater than %d(max. integer representation)"
                         ,(yyvsp[(1) - (1)].name),MAX_INT_VAL);
                  YYERROR;
                }
                char *t1=omStrDup((yyvsp[(1) - (1)].name));
                syMake(&(yyval.lv),t1);
              }
            }
          ;}
    break;

  case 37:
#line 591 "grammar.y"
    {
            memset(&(yyval.lv),0,sizeof((yyval.lv)));
            (yyval.lv).rtyp = (yyvsp[(1) - (1)].i);
            (yyval.lv).data = (yyval.lv).Data();
          ;}
    break;

  case 38:
#line 597 "grammar.y"
    {
            memset(&(yyval.lv),0,sizeof((yyval.lv)));
            (yyval.lv).rtyp  = STRING_CMD;
            (yyval.lv).data = (yyvsp[(1) - (1)].name);
          ;}
    break;

  case 39:
#line 606 "grammar.y"
    {
            leftv v = &(yyvsp[(1) - (3)].lv);
            while (v->next!=NULL)
            {
              v=v->next;
            }
            v->next = (leftv)omAllocBin(sleftv_bin);
            memcpy(v->next,&((yyvsp[(3) - (3)].lv)),sizeof(sleftv));
            (yyval.lv) = (yyvsp[(1) - (3)].lv);
          ;}
    break;

  case 40:
#line 617 "grammar.y"
    {
            (yyval.lv) = (yyvsp[(1) - (1)].lv);
          ;}
    break;

  case 41:
#line 623 "grammar.y"
    {
            /*if ($1.typ == eunknown) YYERROR;*/
            (yyval.lv) = (yyvsp[(1) - (1)].lv);
          ;}
    break;

  case 42:
#line 627 "grammar.y"
    { (yyval.lv) = (yyvsp[(1) - (1)].lv); ;}
    break;

  case 43:
#line 628 "grammar.y"
    { (yyval.lv) = (yyvsp[(2) - (3)].lv); ;}
    break;

  case 44:
#line 630 "grammar.y"
    {
            if(iiExprArith3(&(yyval.lv),'[',&(yyvsp[(1) - (6)].lv),&(yyvsp[(3) - (6)].lv),&(yyvsp[(5) - (6)].lv))) YYERROR;
          ;}
    break;

  case 45:
#line 634 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (4)].lv),'[',&(yyvsp[(3) - (4)].lv))) YYERROR;
          ;}
    break;

  case 46:
#line 638 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(3) - (4)].lv),(yyvsp[(1) - (4)].i))) YYERROR;
          ;}
    break;

  case 47:
#line 642 "grammar.y"
    {
            if(iiExprArithM(&(yyval.lv),&(yyvsp[(3) - (4)].lv),(yyvsp[(1) - (4)].i))) YYERROR;
          ;}
    break;

  case 48:
#line 646 "grammar.y"
    {
            if(iiExprArithM(&(yyval.lv),NULL,(yyvsp[(1) - (3)].i))) YYERROR;
          ;}
    break;

  case 49:
#line 650 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(3) - (4)].lv),(yyvsp[(1) - (4)].i))) YYERROR;
          ;}
    break;

  case 50:
#line 654 "grammar.y"
    {
            if(iiExprArithM(&(yyval.lv),&(yyvsp[(3) - (4)].lv),(yyvsp[(1) - (4)].i))) YYERROR;
          ;}
    break;

  case 51:
#line 658 "grammar.y"
    {
            if(iiExprArithM(&(yyval.lv),NULL,(yyvsp[(1) - (3)].i))) YYERROR;
          ;}
    break;

  case 52:
#line 662 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(3) - (4)].lv),(yyvsp[(1) - (4)].i))) YYERROR;
          ;}
    break;

  case 53:
#line 666 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(3) - (6)].lv),(yyvsp[(1) - (6)].i),&(yyvsp[(5) - (6)].lv),TRUE)) YYERROR;
          ;}
    break;

  case 54:
#line 670 "grammar.y"
    {
            if(iiExprArith3(&(yyval.lv),(yyvsp[(1) - (8)].i),&(yyvsp[(3) - (8)].lv),&(yyvsp[(5) - (8)].lv),&(yyvsp[(7) - (8)].lv))) YYERROR;
          ;}
    break;

  case 55:
#line 674 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(3) - (6)].lv),(yyvsp[(1) - (6)].i),&(yyvsp[(5) - (6)].lv),TRUE)) YYERROR;
          ;}
    break;

  case 56:
#line 678 "grammar.y"
    {
            if(iiExprArith3(&(yyval.lv),(yyvsp[(1) - (8)].i),&(yyvsp[(3) - (8)].lv),&(yyvsp[(5) - (8)].lv),&(yyvsp[(7) - (8)].lv))) YYERROR;
          ;}
    break;

  case 57:
#line 682 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(3) - (4)].lv),(yyvsp[(1) - (4)].i))) YYERROR;
          ;}
    break;

  case 58:
#line 686 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(3) - (4)].lv),(yyvsp[(1) - (4)].i))) YYERROR;
          ;}
    break;

  case 59:
#line 690 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(3) - (6)].lv),(yyvsp[(1) - (6)].i),&(yyvsp[(5) - (6)].lv),TRUE)) YYERROR;
          ;}
    break;

  case 60:
#line 694 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(3) - (4)].lv),(yyvsp[(1) - (4)].i))) YYERROR;
          ;}
    break;

  case 61:
#line 698 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(3) - (6)].lv),(yyvsp[(1) - (6)].i),&(yyvsp[(5) - (6)].lv),TRUE)) YYERROR;
          ;}
    break;

  case 62:
#line 702 "grammar.y"
    {
            if(iiExprArith3(&(yyval.lv),(yyvsp[(1) - (8)].i),&(yyvsp[(3) - (8)].lv),&(yyvsp[(5) - (8)].lv),&(yyvsp[(7) - (8)].lv))) YYERROR;
          ;}
    break;

  case 63:
#line 706 "grammar.y"
    {
            if(iiExprArith3(&(yyval.lv),(yyvsp[(1) - (8)].i),&(yyvsp[(3) - (8)].lv),&(yyvsp[(5) - (8)].lv),&(yyvsp[(7) - (8)].lv))) YYERROR;
          ;}
    break;

  case 64:
#line 710 "grammar.y"
    {
            if(iiExprArithM(&(yyval.lv),NULL,(yyvsp[(1) - (3)].i))) YYERROR;
          ;}
    break;

  case 65:
#line 714 "grammar.y"
    {
            if(iiExprArithM(&(yyval.lv),&(yyvsp[(3) - (4)].lv),(yyvsp[(1) - (4)].i))) YYERROR;
          ;}
    break;

  case 66:
#line 718 "grammar.y"
    {
            if(iiExprArith3(&(yyval.lv),MATRIX_CMD,&(yyvsp[(3) - (8)].lv),&(yyvsp[(5) - (8)].lv),&(yyvsp[(7) - (8)].lv))) YYERROR;
          ;}
    break;

  case 67:
#line 722 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(3) - (4)].lv),MATRIX_CMD)) YYERROR;
          ;}
    break;

  case 68:
#line 726 "grammar.y"
    {
            if(iiExprArith3(&(yyval.lv),INTMAT_CMD,&(yyvsp[(3) - (8)].lv),&(yyvsp[(5) - (8)].lv),&(yyvsp[(7) - (8)].lv))) YYERROR;
          ;}
    break;

  case 69:
#line 730 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(3) - (4)].lv),INTMAT_CMD)) YYERROR;
          ;}
    break;

  case 70:
#line 734 "grammar.y"
    {
            if(iiExprArith3(&(yyval.lv),RING_CMD,&(yyvsp[(3) - (8)].lv),&(yyvsp[(5) - (8)].lv),&(yyvsp[(7) - (8)].lv))) YYERROR;
          ;}
    break;

  case 71:
#line 738 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(3) - (4)].lv),RING_CMD)) YYERROR;
          ;}
    break;

  case 72:
#line 742 "grammar.y"
    {
            (yyval.lv)=(yyvsp[(2) - (3)].lv);
          ;}
    break;

  case 73:
#line 746 "grammar.y"
    {
            #ifdef SIQ
            siq++;
            if (siq>0)
            { if (iiExprArith2(&(yyval.lv),&(yyvsp[(2) - (5)].lv),'=',&(yyvsp[(4) - (5)].lv))) YYERROR; }
            else
            #endif
            {
              memset(&(yyval.lv),0,sizeof((yyval.lv)));
              (yyval.lv).rtyp=NONE;
              if (iiAssign(&(yyvsp[(2) - (5)].lv),&(yyvsp[(4) - (5)].lv))) YYERROR;
            }
            #ifdef SIQ
            siq--;
            #endif
          ;}
    break;

  case 74:
#line 763 "grammar.y"
    {
            #ifdef SIQ
            siq--;
            #endif
          ;}
    break;

  case 75:
#line 769 "grammar.y"
    {
            #ifdef SIQ
            if (siq<=0) (yyvsp[(4) - (5)].lv).Eval();
            #endif
            (yyval.lv)=(yyvsp[(4) - (5)].lv);
            #ifdef SIQ
            siq++;
            #endif
          ;}
    break;

  case 76:
#line 781 "grammar.y"
    {
            #ifdef SIQ
            siq++;
            #endif
          ;}
    break;

  case 77:
#line 789 "grammar.y"
    {
            #ifdef SIQ
            siq--;
            #endif
          ;}
    break;

  case 78:
#line 798 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(1) - (2)].lv),PLUSPLUS)) YYERROR;
          ;}
    break;

  case 79:
#line 802 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(1) - (2)].lv),MINUSMINUS)) YYERROR;
          ;}
    break;

  case 80:
#line 806 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (3)].lv),'+',&(yyvsp[(3) - (3)].lv))) YYERROR;
          ;}
    break;

  case 81:
#line 810 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (3)].lv),'-',&(yyvsp[(3) - (3)].lv))) YYERROR;
          ;}
    break;

  case 82:
#line 814 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (3)].lv),(yyvsp[(2) - (3)].i),&(yyvsp[(3) - (3)].lv))) YYERROR;
          ;}
    break;

  case 83:
#line 818 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (3)].lv),'^',&(yyvsp[(3) - (3)].lv))) YYERROR;
          ;}
    break;

  case 84:
#line 822 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (3)].lv),(yyvsp[(2) - (3)].i),&(yyvsp[(3) - (3)].lv))) YYERROR;
          ;}
    break;

  case 85:
#line 826 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (3)].lv),(yyvsp[(2) - (3)].i),&(yyvsp[(3) - (3)].lv))) YYERROR;
          ;}
    break;

  case 86:
#line 830 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (3)].lv),NOTEQUAL,&(yyvsp[(3) - (3)].lv))) YYERROR;
          ;}
    break;

  case 87:
#line 834 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (3)].lv),EQUAL_EQUAL,&(yyvsp[(3) - (3)].lv))) YYERROR;
          ;}
    break;

  case 88:
#line 838 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (3)].lv),DOTDOT,&(yyvsp[(3) - (3)].lv))) YYERROR;
          ;}
    break;

  case 89:
#line 842 "grammar.y"
    {
            if(iiExprArith2(&(yyval.lv),&(yyvsp[(1) - (3)].lv),':',&(yyvsp[(3) - (3)].lv))) YYERROR;
          ;}
    break;

  case 90:
#line 846 "grammar.y"
    {
            memset(&(yyval.lv),0,sizeof((yyval.lv)));
            int i; TESTSETINT((yyvsp[(2) - (2)].lv),i);
            (yyval.lv).rtyp  = INT_CMD;
            (yyval.lv).data = (void *)(long)(i == 0 ? 1 : 0);
          ;}
    break;

  case 91:
#line 853 "grammar.y"
    {
            if(iiExprArith1(&(yyval.lv),&(yyvsp[(2) - (2)].lv),'-')) YYERROR;
          ;}
    break;

  case 92:
#line 859 "grammar.y"
    { (yyval.lv) = (yyvsp[(1) - (2)].lv); ;}
    break;

  case 93:
#line 861 "grammar.y"
    {
            if ((yyvsp[(1) - (2)].lv).rtyp==0)
            {
              Werror("`%s` is undefined",(yyvsp[(1) - (2)].lv).Fullname());
              YYERROR;
            }
            (yyval.lv) = (yyvsp[(1) - (2)].lv);
          ;}
    break;

  case 95:
#line 875 "grammar.y"
    {
            if ((yyvsp[(2) - (3)].lv).Typ()!=STRING_CMD)
            {
              MYYERROR("string expression expected");
            }
            (yyval.name) = (char *)(yyvsp[(2) - (3)].lv).CopyD(STRING_CMD);
            (yyvsp[(2) - (3)].lv).CleanUp();
          ;}
    break;

  case 98:
#line 892 "grammar.y"
    {
            #ifdef HAVE_NS
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (2)].lv),myynest,(yyvsp[(1) - (2)].i),&((yyvsp[(2) - (2)].lv).req_packhdl->idroot)))
              YYERROR;
            #else
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (2)].lv),myynest,(yyvsp[(1) - (2)].i),&IDROOT)) YYERROR;
            #endif
          ;}
    break;

  case 99:
#line 901 "grammar.y"
    {
            #ifdef HAVE_NS
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (2)].lv),myynest,(yyvsp[(1) - (2)].i),&((yyvsp[(2) - (2)].lv).req_packhdl->idroot)))
              YYERROR;
            #else
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (2)].lv),myynest,(yyvsp[(1) - (2)].i),&IDROOT)) YYERROR;
            #endif
          ;}
    break;

  case 100:
#line 910 "grammar.y"
    {
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (2)].lv),myynest,(yyvsp[(1) - (2)].i),&(currRing->idroot), TRUE)) YYERROR;
          ;}
    break;

  case 101:
#line 914 "grammar.y"
    {
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (2)].lv),myynest,(yyvsp[(1) - (2)].i),&(currRing->idroot), TRUE)) YYERROR;
          ;}
    break;

  case 102:
#line 918 "grammar.y"
    {
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (8)].lv),myynest,(yyvsp[(1) - (8)].i),&(currRing->idroot), TRUE)) YYERROR;
            int r; TESTSETINT((yyvsp[(4) - (8)].lv),r);
            int c; TESTSETINT((yyvsp[(7) - (8)].lv),c);
            if (r < 1)
              MYYERROR("rows must be greater than 0");
            if (c < 0)
              MYYERROR("cols must be greater than -1");
            leftv v=&(yyval.lv);
            //while (v->next!=NULL) { v=v->next; }
            idhdl h=(idhdl)v->data;
            idDelete(&IDIDEAL(h));
            IDMATRIX(h) = mpNew(r,c);
            if (IDMATRIX(h)==NULL) YYERROR;
          ;}
    break;

  case 103:
#line 934 "grammar.y"
    {
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (2)].lv),myynest,(yyvsp[(1) - (2)].i),&(currRing->idroot), TRUE)) YYERROR;
          ;}
    break;

  case 104:
#line 938 "grammar.y"
    {
            int r; TESTSETINT((yyvsp[(4) - (8)].lv),r);
            int c; TESTSETINT((yyvsp[(7) - (8)].lv),c);
            if (r < 1)
              MYYERROR("rows must be greater than 0");
            if (c < 0)
              MYYERROR("cols must be greater than -1");
            #ifdef HAVE_NS
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (8)].lv),myynest,(yyvsp[(1) - (8)].i),&((yyvsp[(2) - (8)].lv).req_packhdl->idroot)))
              YYERROR;
            #else
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (8)].lv),myynest,(yyvsp[(1) - (8)].i),&IDROOT)) YYERROR;
            #endif
            leftv v=&(yyval.lv);
            idhdl h=(idhdl)v->data;
            delete IDINTVEC(h);
            IDINTVEC(h) = new intvec(r,c,0);
            if (IDINTVEC(h)==NULL) YYERROR;
          ;}
    break;

  case 105:
#line 958 "grammar.y"
    {
            #ifdef HAVE_NS
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (2)].lv),myynest,(yyvsp[(1) - (2)].i),&((yyvsp[(2) - (2)].lv).req_packhdl->idroot)))
              YYERROR;
            #else
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (2)].lv),myynest,(yyvsp[(1) - (2)].i),&IDROOT)) YYERROR;
            #endif
            leftv v=&(yyval.lv);
            idhdl h;
            do
            {
               h=(idhdl)v->data;
               delete IDINTVEC(h);
               IDINTVEC(h) = new intvec(1,1,0);
               v=v->next;
            } while (v!=NULL);
          ;}
    break;

  case 106:
#line 976 "grammar.y"
    {
            int t=(yyvsp[(1) - (3)].lv).Typ();
            sleftv r;
            memset(&r,0,sizeof(sleftv));
            if ((BEGIN_RING<t) && (t<END_RING))
            {
              if (iiDeclCommand(&r,&(yyvsp[(3) - (3)].lv),myynest,t,&(currRing->idroot), TRUE))
                YYERROR;
            }
            else
            {
              #ifdef HAVE_NS
              if (iiDeclCommand(&r,&(yyvsp[(3) - (3)].lv),myynest,t,&((yyvsp[(3) - (3)].lv).req_packhdl->idroot)))
                YYERROR;
              #else
              if (iiDeclCommand(&r,&(yyvsp[(3) - (3)].lv),myynest,t,&IDROOT)) YYERROR;
              #endif
            }
            leftv v=&(yyvsp[(1) - (3)].lv);
            while (v->next!=NULL) v=v->next;
            v->next=(leftv)omAllocBin(sleftv_bin);
            memcpy(v->next,&r,sizeof(sleftv));
            (yyval.lv)=(yyvsp[(1) - (3)].lv);
          ;}
    break;

  case 107:
#line 1001 "grammar.y"
    {
            #ifdef HAVE_NS
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (2)].lv),myynest,(yyvsp[(1) - (2)].i),&((yyvsp[(2) - (2)].lv).req_packhdl->idroot)))
              YYERROR;
            #else
            if (iiDeclCommand(&(yyval.lv),&(yyvsp[(2) - (2)].lv),myynest,(yyvsp[(1) - (2)].i),&IDROOT)) YYERROR;
            #endif
          ;}
    break;

  case 110:
#line 1018 "grammar.y"
    {
            leftv v = &(yyvsp[(2) - (5)].lv);
            while (v->next!=NULL)
            {
              v=v->next;
            }
            v->next = (leftv)omAllocBin(sleftv_bin);
            memcpy(v->next,&((yyvsp[(4) - (5)].lv)),sizeof(sleftv));
            (yyval.lv) = (yyvsp[(2) - (5)].lv);
          ;}
    break;

  case 111:
#line 1032 "grammar.y"
    {
          // let rInit take care of any errors
          (yyval.i)=rOrderName((yyvsp[(1) - (1)].name));
        ;}
    break;

  case 112:
#line 1040 "grammar.y"
    {
            memset(&(yyval.lv),0,sizeof((yyval.lv)));
            intvec *iv = new intvec(2);
            (*iv)[0] = 1;
            (*iv)[1] = (yyvsp[(1) - (1)].i);
            (yyval.lv).rtyp = INTVEC_CMD;
            (yyval.lv).data = (void *)iv;
          ;}
    break;

  case 113:
#line 1049 "grammar.y"
    {
            memset(&(yyval.lv),0,sizeof((yyval.lv)));
            leftv sl = &(yyvsp[(3) - (4)].lv);
            int slLength;
            {
              slLength =  exprlist_length(sl);
              int l = 2 +  slLength;
              intvec *iv = new intvec(l);
              (*iv)[0] = slLength;
              (*iv)[1] = (yyvsp[(1) - (4)].i);

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
                  (yyvsp[(3) - (4)].lv).CleanUp();
                  MYYERROR("wrong type in ordering");
                }
                sl = sl->next;
              }
              (yyval.lv).rtyp = INTVEC_CMD;
              (yyval.lv).data = (void *)iv;
            }
            (yyvsp[(3) - (4)].lv).CleanUp();
          ;}
    break;

  case 115:
#line 1095 "grammar.y"
    {
            (yyval.lv) = (yyvsp[(1) - (3)].lv);
            (yyval.lv).next = (sleftv *)omAllocBin(sleftv_bin);
            memcpy((yyval.lv).next,&(yyvsp[(3) - (3)].lv),sizeof(sleftv));
          ;}
    break;

  case 117:
#line 1105 "grammar.y"
    {
            (yyval.lv) = (yyvsp[(2) - (3)].lv);
          ;}
    break;

  case 118:
#line 1111 "grammar.y"
    {
            expected_parms = TRUE;
          ;}
    break;

  case 119:
#line 1123 "grammar.y"
    { if ((yyvsp[(1) - (2)].i) != '<') YYERROR;
            if((feFilePending=feFopen((yyvsp[(2) - (2)].name),"r",NULL,TRUE))==NULL) YYERROR; ;}
    break;

  case 120:
#line 1126 "grammar.y"
    { newFile((yyvsp[(2) - (4)].name),feFilePending); ;}
    break;

  case 121:
#line 1131 "grammar.y"
    {
            feHelp((yyvsp[(2) - (3)].name));
            omFree((ADDRESS)(yyvsp[(2) - (3)].name));
          ;}
    break;

  case 122:
#line 1136 "grammar.y"
    {
            feHelp(NULL);
          ;}
    break;

  case 123:
#line 1143 "grammar.y"
    {
            singular_example((yyvsp[(2) - (3)].name));
            omFree((ADDRESS)(yyvsp[(2) - (3)].name));
          ;}
    break;

  case 124:
#line 1151 "grammar.y"
    {
#ifdef HAVE_NS
          if (basePack!=(yyvsp[(2) - (2)].lv).req_packhdl)
          {
            if(iiExport(&(yyvsp[(2) - (2)].lv),0,currPackHdl)) YYERROR;
          }
          else
#endif /* HAVE_NS */
            if (iiExport(&(yyvsp[(2) - (2)].lv),0)) YYERROR;
        ;}
    break;

  case 125:
#line 1165 "grammar.y"
    {
          leftv v=&(yyvsp[(2) - (2)].lv);
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
        ;}
    break;

  case 126:
#line 1185 "grammar.y"
    {
          leftv v=&(yyvsp[(3) - (3)].lv);
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
        ;}
    break;

  case 127:
#line 1208 "grammar.y"
    {
            list_cmd((yyvsp[(3) - (4)].i),NULL,"// ",TRUE);
          ;}
    break;

  case 128:
#line 1212 "grammar.y"
    {
            list_cmd((yyvsp[(3) - (4)].i),NULL,"// ",TRUE);
          ;}
    break;

  case 129:
#line 1216 "grammar.y"
    {
            if ((yyvsp[(3) - (4)].i)==QRING_CMD) (yyvsp[(3) - (4)].i)=RING_CMD;
            list_cmd((yyvsp[(3) - (4)].i),NULL,"// ",TRUE);
          ;}
    break;

  case 130:
#line 1221 "grammar.y"
    {
            list_cmd((yyvsp[(3) - (4)].i),NULL,"// ",TRUE);
          ;}
    break;

  case 131:
#line 1225 "grammar.y"
    {
            list_cmd(RING_CMD,NULL,"// ",TRUE);
          ;}
    break;

  case 132:
#line 1229 "grammar.y"
    {
            list_cmd(MATRIX_CMD,NULL,"// ",TRUE);
           ;}
    break;

  case 133:
#line 1233 "grammar.y"
    {
            list_cmd(INTMAT_CMD,NULL,"// ",TRUE);
          ;}
    break;

  case 134:
#line 1237 "grammar.y"
    {
            list_cmd(PROC_CMD,NULL,"// ",TRUE);
          ;}
    break;

  case 135:
#line 1241 "grammar.y"
    {
            list_cmd(0,(yyvsp[(3) - (4)].lv).Fullname(),"// ",TRUE);
            (yyvsp[(3) - (4)].lv).CleanUp();
          ;}
    break;

  case 136:
#line 1246 "grammar.y"
    {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if((yyvsp[(3) - (6)].lv).Typ() == PACKAGE_CMD)
              list_cmd((yyvsp[(5) - (6)].i),NULL,"// ",TRUE);
#endif /* HAVE_NS */
            (yyvsp[(3) - (6)].lv).CleanUp();
          ;}
    break;

  case 137:
#line 1255 "grammar.y"
    {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if((yyvsp[(3) - (6)].lv).Typ() == PACKAGE_CMD)
              list_cmd((yyvsp[(5) - (6)].i),NULL,"// ",TRUE);
#endif /* HAVE_NS */
            (yyvsp[(3) - (6)].lv).CleanUp();
          ;}
    break;

  case 138:
#line 1264 "grammar.y"
    {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if((yyvsp[(3) - (6)].lv).Typ() == PACKAGE_CMD)
              list_cmd((yyvsp[(5) - (6)].i),NULL,"// ",TRUE);
#endif /* HAVE_NS */
            (yyvsp[(3) - (6)].lv).CleanUp();
          ;}
    break;

  case 139:
#line 1273 "grammar.y"
    {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if((yyvsp[(3) - (6)].lv).Typ() == PACKAGE_CMD)
              list_cmd((yyvsp[(5) - (6)].i),NULL,"// ",TRUE);
#endif /* HAVE_NS */
            (yyvsp[(3) - (6)].lv).CleanUp();
          ;}
    break;

  case 140:
#line 1282 "grammar.y"
    {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if((yyvsp[(3) - (6)].lv).Typ() == PACKAGE_CMD)
              list_cmd((yyvsp[(5) - (6)].i),NULL,"// ",TRUE);
#endif /* HAVE_NS */
            (yyvsp[(3) - (6)].lv).CleanUp();
          ;}
    break;

  case 141:
#line 1291 "grammar.y"
    {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if((yyvsp[(3) - (6)].lv).Typ() == PACKAGE_CMD)
              list_cmd((yyvsp[(5) - (6)].i),NULL,"// ",TRUE);
#endif /* HAVE_NS */
            (yyvsp[(3) - (6)].lv).CleanUp();
          ;}
    break;

  case 142:
#line 1300 "grammar.y"
    {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if((yyvsp[(3) - (6)].lv).Typ() == PACKAGE_CMD)
              list_cmd((yyvsp[(5) - (6)].i),NULL,"// ",TRUE);
#endif /* HAVE_NS */
            (yyvsp[(3) - (6)].lv).CleanUp();
          ;}
    break;

  case 143:
#line 1309 "grammar.y"
    {
#ifdef HAVE_NS
            //PrintS("?????\n");
            if((yyvsp[(3) - (6)].lv).Typ() == PACKAGE_CMD)
              list_cmd((yyvsp[(5) - (6)].i),NULL,"// ",TRUE);
#endif /* HAVE_NS */
            (yyvsp[(3) - (6)].lv).CleanUp();
          ;}
    break;

  case 144:
#line 1327 "grammar.y"
    {
            list_cmd(-1,NULL,"// ",TRUE);
          ;}
    break;

  case 145:
#line 1333 "grammar.y"
    { yyInRingConstruction = TRUE; ;}
    break;

  case 146:
#line 1342 "grammar.y"
    {
            BOOLEAN do_pop = FALSE;
            char *ring_name = (yyvsp[(2) - (8)].lv).name;
            ring b=
            rInit(&(yyvsp[(4) - (8)].lv),            /* characteristik and list of parameters*/
                  &(yyvsp[(6) - (8)].lv),            /* names of ringvariables */
                  &(yyvsp[(8) - (8)].lv));            /* ordering */
            idhdl newRingHdl=NULL;

            if (b!=NULL)
            {
              #ifdef HAVE_NS
                newRingHdl=enterid(ring_name, myynest, RING_CMD,
                                   &((yyvsp[(2) - (8)].lv).req_packhdl->idroot));
              #else
                newRingHdl=enterid(ring_name, myynest, RING_CMD, &IDROOT);
              #endif
              (yyvsp[(2) - (8)].lv).CleanUp();
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
          ;}
    break;

  case 147:
#line 1381 "grammar.y"
    {
            BOOLEAN do_pop = FALSE;
            char *ring_name = (yyvsp[(2) - (2)].lv).name;
            if (!inerror) rDefault(ring_name);
            yyInRingConstruction = FALSE;
            (yyvsp[(2) - (2)].lv).CleanUp();
          ;}
    break;

  case 148:
#line 1392 "grammar.y"
    {
            if (((yyvsp[(1) - (2)].i)!=LIB_CMD)||(iiLibCmd((yyvsp[(2) - (2)].name),TRUE,TRUE,TRUE)))
            //if ($1==LIB_CMD)
            //{
            //  sleftv tmp;
            //  if(iiExprArith1(&tmp,&$2,LIB_CMD)) YYERROR;
            //}
            //else
                YYERROR;
          ;}
    break;

  case 151:
#line 1408 "grammar.y"
    {
            if (((yyvsp[(1) - (2)].i)==KEEPRING_CMD) && (myynest==0))
               MYYERROR("only inside a proc allowed");
            const char * n=(yyvsp[(2) - (2)].lv).Name();
            if ((((yyvsp[(2) - (2)].lv).Typ()==RING_CMD)||((yyvsp[(2) - (2)].lv).Typ()==QRING_CMD))
            && ((yyvsp[(2) - (2)].lv).rtyp==IDHDL))
            {
              idhdl h=(idhdl)(yyvsp[(2) - (2)].lv).data;
              if ((yyvsp[(2) - (2)].lv).e!=NULL) h=rFindHdl((ring)(yyvsp[(2) - (2)].lv).Data(),NULL, NULL);
              //Print("setring %s lev %d (ptr:%x)\n",IDID(h),IDLEV(h),IDRING(h));
              if ((yyvsp[(1) - (2)].i)==KEEPRING_CMD)
              {
                if (h!=NULL)
                {
                  if (IDLEV(h)!=0)
                  {
                    if (iiExport(&(yyvsp[(2) - (2)].lv),myynest-1)) YYERROR;
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
                  (yyvsp[(2) - (2)].lv).CleanUp();
                  YYERROR;
                }
              }
              if (h!=NULL) rSetHdl(h);
              else
              {
                Werror("cannot find the name of the basering %s",n);
                (yyvsp[(2) - (2)].lv).CleanUp();
                YYERROR;
              }
              (yyvsp[(2) - (2)].lv).CleanUp();
            }
            else
            {
              Werror("%s is no name of a ring/qring",n);
              (yyvsp[(2) - (2)].lv).CleanUp();
              YYERROR;
            }
          ;}
    break;

  case 152:
#line 1484 "grammar.y"
    {
            if ((yyvsp[(2) - (2)].lv).rtyp!=IDHDL) MYYERROR("identifier expected");
            idhdl h = (idhdl)(yyvsp[(2) - (2)].lv).data;
            type_cmd(h);
          ;}
    break;

  case 153:
#line 1490 "grammar.y"
    {
            //Print("typ is %d, rtyp:%d\n",$1.Typ(),$1.rtyp);
            #ifdef SIQ
            if ((yyvsp[(1) - (1)].lv).rtyp!=COMMAND)
            {
            #endif
              if ((yyvsp[(1) - (1)].lv).Typ()==UNKNOWN)
              {
                if ((yyvsp[(1) - (1)].lv).name!=NULL)
                {
                  Werror("`%s` is undefined",(yyvsp[(1) - (1)].lv).name);
                  omFree((ADDRESS)(yyvsp[(1) - (1)].lv).name);
                }
                YYERROR;
              }
            #ifdef SIQ
            }
            #endif
            (yyvsp[(1) - (1)].lv).Print(&sLastPrinted);
            (yyvsp[(1) - (1)].lv).CleanUp(currRing);
            if (errorreported) YYERROR;
          ;}
    break;

  case 154:
#line 1519 "grammar.y"
    {
            int i; TESTSETINT((yyvsp[(3) - (5)].lv),i);
            if (i!=0)
            {
              newBuffer( (yyvsp[(5) - (5)].name), BT_if);
            }
            else
            {
              omFree((ADDRESS)(yyvsp[(5) - (5)].name));
              currentVoice->ifsw=1;
            }
          ;}
    break;

  case 155:
#line 1532 "grammar.y"
    {
            if (currentVoice->ifsw==1)
            {
              currentVoice->ifsw=0;
              newBuffer( (yyvsp[(2) - (2)].name), BT_else);
            }
            else
            {
              if (currentVoice->ifsw!=2)
              {
                Warn("`else` without `if` in level %d",myynest);
              }
              omFree((ADDRESS)(yyvsp[(2) - (2)].name));
            }
            currentVoice->ifsw=0;
          ;}
    break;

  case 156:
#line 1549 "grammar.y"
    {
            int i; TESTSETINT((yyvsp[(3) - (5)].lv),i);
            if (i)
            {
              if (exitBuffer(BT_break)) YYERROR;
            }
            currentVoice->ifsw=0;
          ;}
    break;

  case 157:
#line 1558 "grammar.y"
    {
            if (exitBuffer(BT_break)) YYERROR;
            currentVoice->ifsw=0;
          ;}
    break;

  case 158:
#line 1563 "grammar.y"
    {
            if (contBuffer(BT_break)) YYERROR;
            currentVoice->ifsw=0;
          ;}
    break;

  case 159:
#line 1571 "grammar.y"
    {
            /* -> if(!$2) break; $3; continue;*/
            char * s = (char *)omAlloc( strlen((yyvsp[(2) - (3)].name)) + strlen((yyvsp[(3) - (3)].name)) + 36);
            sprintf(s,"whileif (!(%s)) break;\n%scontinue;\n " ,(yyvsp[(2) - (3)].name),(yyvsp[(3) - (3)].name));
            newBuffer(s,BT_break);
            omFree((ADDRESS)(yyvsp[(2) - (3)].name));
            omFree((ADDRESS)(yyvsp[(3) - (3)].name));
          ;}
    break;

  case 160:
#line 1583 "grammar.y"
    {
            /* $2 */
            /* if (!$3) break; $5; $4; continue; */
            char * s = (char *)omAlloc( strlen((yyvsp[(3) - (5)].name))+strlen((yyvsp[(4) - (5)].name))+strlen((yyvsp[(5) - (5)].name))+36);
            sprintf(s,"forif (!(%s)) break;\n%s%s;\ncontinue;\n "
                   ,(yyvsp[(3) - (5)].name),(yyvsp[(5) - (5)].name),(yyvsp[(4) - (5)].name));
            omFree((ADDRESS)(yyvsp[(3) - (5)].name));
            omFree((ADDRESS)(yyvsp[(4) - (5)].name));
            omFree((ADDRESS)(yyvsp[(5) - (5)].name));
            newBuffer(s,BT_break);
            s = (char *)omAlloc( strlen((yyvsp[(2) - (5)].name)) + 3);
            sprintf(s,"%s;\n",(yyvsp[(2) - (5)].name));
            omFree((ADDRESS)(yyvsp[(2) - (5)].name));
            newBuffer(s,BT_if);
          ;}
    break;

  case 161:
#line 1602 "grammar.y"
    {
            procinfov pi;
            idhdl h = enterid((yyvsp[(2) - (3)].name),myynest,PROC_CMD,&IDROOT,TRUE);
            if (h==NULL) {omFree((ADDRESS)(yyvsp[(2) - (3)].name));omFree((ADDRESS)(yyvsp[(3) - (3)].name)); YYERROR;}
            iiInitSingularProcinfo(IDPROC(h),"", (yyvsp[(2) - (3)].name), 0, 0);
            IDPROC(h)->data.s.body = (char *)omAlloc(strlen((yyvsp[(3) - (3)].name))+31);;
            sprintf(IDPROC(h)->data.s.body,"parameter list #;\n%s;return();\n\n",(yyvsp[(3) - (3)].name));
            omFree((ADDRESS)(yyvsp[(3) - (3)].name));
            omFree((ADDRESS)(yyvsp[(2) - (3)].name));
          ;}
    break;

  case 162:
#line 1613 "grammar.y"
    {
            idhdl h = enterid((yyvsp[(1) - (3)].name),myynest,PROC_CMD,&IDROOT,TRUE);
            if (h==NULL)
            {
              omFree((ADDRESS)(yyvsp[(1) - (3)].name));
              omFree((ADDRESS)(yyvsp[(2) - (3)].name));
              omFree((ADDRESS)(yyvsp[(3) - (3)].name));
              YYERROR;
            }
            char *args=iiProcArgs((yyvsp[(2) - (3)].name),FALSE);
            omFree((ADDRESS)(yyvsp[(2) - (3)].name));
            procinfov pi;
            iiInitSingularProcinfo(IDPROC(h),"", (yyvsp[(1) - (3)].name), 0, 0);
            IDPROC(h)->data.s.body = (char *)omAlloc(strlen((yyvsp[(3) - (3)].name))+strlen(args)+14);;
            sprintf(IDPROC(h)->data.s.body,"%s\n%s;return();\n\n",args,(yyvsp[(3) - (3)].name));
            omFree((ADDRESS)args);
            omFree((ADDRESS)(yyvsp[(3) - (3)].name));
            omFree((ADDRESS)(yyvsp[(1) - (3)].name));
          ;}
    break;

  case 163:
#line 1633 "grammar.y"
    {
            omFree((ADDRESS)(yyvsp[(3) - (4)].name));
            idhdl h = enterid((yyvsp[(1) - (4)].name),myynest,PROC_CMD,&IDROOT,TRUE);
            if (h==NULL)
            {
              omFree((ADDRESS)(yyvsp[(1) - (4)].name));
              omFree((ADDRESS)(yyvsp[(2) - (4)].name));
              omFree((ADDRESS)(yyvsp[(4) - (4)].name));
              YYERROR;
            }
            char *args=iiProcArgs((yyvsp[(2) - (4)].name),FALSE);
            omFree((ADDRESS)(yyvsp[(2) - (4)].name));
            procinfov pi;
            iiInitSingularProcinfo(IDPROC(h),"", (yyvsp[(1) - (4)].name), 0, 0);
            omFree((ADDRESS)(yyvsp[(1) - (4)].name));
            IDPROC(h)->data.s.body = (char *)omAlloc(strlen((yyvsp[(4) - (4)].name))+strlen(args)+14);;
            sprintf(IDPROC(h)->data.s.body,"%s\n%s;return();\n\n",args,(yyvsp[(4) - (4)].name));
            omFree((ADDRESS)args);
            omFree((ADDRESS)(yyvsp[(4) - (4)].name));
          ;}
    break;

  case 164:
#line 1657 "grammar.y"
    {
            //Print("par:%s, %d\n",$2.Name(),$2.Typ());
            //yylineno--;
            if (iiParameter(&(yyvsp[(2) - (2)].lv))) YYERROR;
          ;}
    break;

  case 165:
#line 1663 "grammar.y"
    {
            //Print("par:%s, %d\n",$2.Name(),$2.Typ());
            sleftv tmp_expr;
            //yylineno--;
            if ((iiDeclCommand(&tmp_expr,&(yyvsp[(2) - (2)].lv),myynest,DEF_CMD,&IDROOT))
            || (iiParameter(&tmp_expr)))
              YYERROR;
          ;}
    break;

  case 166:
#line 1675 "grammar.y"
    {
            if(iiRETURNEXPR==NULL) YYERROR;
            iiRETURNEXPR[myynest].Copy(&(yyvsp[(3) - (4)].lv));
            (yyvsp[(3) - (4)].lv).CleanUp();
            if (exitBuffer(BT_proc)) YYERROR;
          ;}
    break;

  case 167:
#line 1682 "grammar.y"
    {
            if ((yyvsp[(1) - (3)].i)==RETURN)
            {
              if(iiRETURNEXPR==NULL) YYERROR;
              iiRETURNEXPR[myynest].Init();
              iiRETURNEXPR[myynest].rtyp=NONE;
              if (exitBuffer(BT_proc)) YYERROR;
            }
          ;}
    break;


/* Line 1267 of yacc.c.  */
#line 3976 "grammar.cc"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



