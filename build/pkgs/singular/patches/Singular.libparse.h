#ifndef LIBPARSE_H
#define LIBPARSE_H
/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/* $Id: libparse.h,v 1.15 2003/12/16 11:25:34 Singular Exp $ */
/*
* ABSTRACT: lib parsing
*/
#ifndef STANDALONE_PARSER
#  include "structs.h"
#  include "subexpr.h"
#endif
typedef enum { LOAD_LIB, GET_INFO } lp_modes;
typedef enum { OLD_LIBSTYLE, NEW_LIBSTYLE } lib_style_types;

procinfo *iiInitSingularProcinfo(procinfov pi, char *libname,
                                 char *procname, int line, long pos,
                                 BOOLEAN pstatic = FALSE);
#ifdef HAVE_NS
int yylplex(char *libname, char *libfile, lib_style_types *lib_style,
           idhdl pl, BOOLEAN autoexport=FALSE, lp_modes=LOAD_LIB);
#else
int yylplex(char *libname, char *libfile, lib_style_types *lib_style,
            lp_modes=LOAD_LIB);
#endif /* HAVE_NS */

void reinit_yylp();

extern char * text_buffer;

#  define YYLP_ERR_NONE    0
#  define YYLP_DEF_BR2     1
#  define YYLP_BODY_BR2    2
#  define YYLP_BODY_BR3    3
#  define YYLP_BODY_TMBR2  4
#  define YYLP_BODY_TMBR3  5
#  define YYLP_EX_BR2      6
#  define YYLP_EX_BR3      7
#  define YYLP_BAD_CHAR    8
#  define YYLP_MISSQUOT    9
#  define YYLP_MISS_BR1   10
#  define YYLP_MISS_BR2   11
#  define YYLP_MISS_BR3   12

#  ifdef STANDALONE_PARSER
#ifndef unix
extern FILE* myfopen(char *path, char *mode);
extern size_t myfread(void *ptr, size_t size, size_t nmemb, FILE *stream);
#else
#define myfopen fopen
#define myfread fread
#endif
#  endif

#endif /* LIBPARSE_H */


