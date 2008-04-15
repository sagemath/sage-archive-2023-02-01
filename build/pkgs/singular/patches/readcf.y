/* emacs edit mode for this file is -*- C++ -*- */
/* $Id: readcf.y,v 1.10 2006/05/15 09:03:07 Singular Exp $ */

%{

#if defined(WINNT) && ! defined(__GNUC__)
#include <malloc.h>
#include <memory.h>
#define alloca _alloca
#endif

#include <cstring> // we need this fot gcc 4.3
#include <config.h>

#include <ctype.h>
#ifdef HAVE_IOSTREAM
#include <iostream>
#define ISTREAM std::istream
#define CERR std::cerr
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#define ISTREAM istream
#define CERR cerr
#endif


#include "assert.h"

#include "canonicalform.h"
#include "cf_defs.h"
#include "gfops.h"
#include "parseutil.h"
#include "variable.h"

#ifndef BISONPP
#define YYSTYPE ParseUtil
#else
#define YY_parse_USE_GOTO 1
#define YY_parse_STYPE ParseUtil
#endif

static char* readString( ISTREAM& );

#ifndef BISONPP
void yyerror( char * s );
int yylex();
#endif

static ISTREAM * defaultin = 0;

static CanonicalForm * retvalue = 0;

%}

/* BISON Declarations */

%token NUM
%left '-' '+'
%left '*' '/'
%left NEG
%right '^'

/* Grammar follows */

%%
input:	/* empty string */
	| input line
;

line:	';'
	| exp ';' { *retvalue = $1.getval(); return 0; }
;

exp:	NUM			{ $$ = $1; }
	| exp '+' exp		{ $$ = $1.getval() + $3.getval(); }
	| exp '-' exp		{ $$ = $1.getval() - $3.getval(); }
	| exp '*' exp		{ $$ = $1.getval() * $3.getval(); }
	| exp '/' exp		{ $$ = $1.getval() / $3.getval(); }
	| '-' exp %prec NEG	{ $$ = -$2.getval(); }
	| '+' exp %prec NEG	{ $$ = $2.getval(); }
	| exp '^' NUM		{ $$ = power( $1.getval(), $3.getintval() ); }
	| '(' exp ')'		{ $$ = $2.getval(); }
;

%%

#ifdef BISONPP
void YY_parse_CLASS::yyerror( char * s )
#else
void yyerror( char * s )
#endif
{
    CERR << s << "\n";
}

#ifdef BISONPP
int YY_parse_CLASS::yylex()
#else
int yylex()
#endif
{
    int c;

    while ((c = defaultin->get()) == ' ' || c == '\t' || c == '\n' ) ;
    if ( isdigit( c ) ) {
	defaultin->putback( c );
	yylval = ParseUtil( readString( *defaultin ) );
	return NUM;
    }
    else if ( isalpha( c ) ) {
	// look for generators of GF(q)
	if ( getCharacteristic() > 0 && getGFDegree() > 1 && c == gf_name ) {
#ifdef BISONPP
	    this->yylval = getGFGenerator();
#else
	    yylval = getGFGenerator();
#endif
	}
	else if ( c == getDefaultVarName() ) {
	    int cc;
	    cc = defaultin->get();
	    if ( cc == '_' ) {
		ParseUtil index( readString( *defaultin ) );
#ifdef BISONPP
		this->yylval = Variable( index.getintval() );
#else
		yylval = Variable( index.getintval() );
#endif
	    }
	    else {
		defaultin->putback( cc );
#ifdef BISONPP
		this->yylval = Variable( (char)c );
#else
		yylval = Variable( (char)c );
#endif
	    }
	}
	else {
#ifdef BISONPP
	    this->yylval = Variable( (char)c );
#else
	    yylval = Variable( (char)c );
#endif
	}
	return NUM;
    }
    return c;
}

CanonicalForm readCF( ISTREAM& str )
{
    CanonicalForm theRetvalue;
    retvalue = new CanonicalForm();
#ifdef BISONPP
    YY_parse_CLASS my_parser;
    defaultin = &str;
    if ( my_parser.yyparse() == 0 ) {
	theRetvalue = *retvalue;
	delete retvalue;
	return theRetvalue;
    }
    else {
	delete retvalue;
	return 0;
    }
#else
    defaultin = &str;
    if ( yyparse() == 0 ) {
	theRetvalue = *retvalue;
	delete retvalue;
	return theRetvalue;
    }
    else {
	delete retvalue;
	return 0;
    }
#endif
}

char* readString( ISTREAM& s )
{
    static char * buffer = 0;
    static int bufsize = 0;

    if ( buffer == 0 ) {
	bufsize = 10000;
	buffer = new char[bufsize];
    }
    int i = 0, c, goon = 1;
    while ( goon ) {
	while ( isdigit( c = s.get() ) && i < bufsize - 2 ) {
	    buffer[i] = c;
	    i++;
	}
	if ( isdigit( c ) ) {
	    bufsize += 1000;
	    char * newbuffer = (char*)memcpy( new char[bufsize], buffer, bufsize - 1000 );
	    delete [] buffer;
	    buffer = newbuffer;
	    buffer[i] = c;
	    i++;
	}
	else {
	    goon = 0;
	    buffer[i] = '\0';
	    s.putback( c );
	}
    }
    return buffer;
}
