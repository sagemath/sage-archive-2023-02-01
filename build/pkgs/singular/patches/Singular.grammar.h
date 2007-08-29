/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

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




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



