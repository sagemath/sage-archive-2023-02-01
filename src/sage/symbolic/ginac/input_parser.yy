/** @file input_parser.yy
 *
 *  Input grammar definition for reading expressions.
 *  This file must be processed with yacc/bison. */

/*
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


/*
 *  Definitions
 */

%{
#include <stdexcept>

#include "ex.h"
#include "input_lexer.h"
#include "relational.h"
#include "operators.h"
#include "symbol.h"
#include "lst.h"
#include "power.h"
#include "exprseq.h"
#include "idx.h"
#include "indexed.h"
#include "matrix.h"
#include "inifcns.h"

namespace GiNaC {

#define YYERROR_VERBOSE 1

// Parsed output expression
ex parsed_ex;

// Last error message returned by parser
static std::string parser_error;

// Prototypes
ex attach_index(const ex & base, ex i, bool covariant);
%}

/* Tokens (T_LITERAL means a literal value returned by the parser, but not
   of class numeric or symbol (e.g. a constant or the FAIL object)) */
%token T_EOF T_NUMBER T_SYMBOL T_LITERAL T_DIGITS T_EQUAL T_NOTEQ T_LESSEQ T_GREATEREQ

/* Operator precedence and associativity */
%right '='
%left T_EQUAL T_NOTEQ
%left '<' '>' T_LESSEQ T_GREATEREQ
%left '+' '-'
%left '*' '/' '%'
%nonassoc NEG
%right '^'
%left '.' '~'
%nonassoc '!'

%start input


/*
 *  Grammar rules
 */

%%
input	: exp T_EOF {
		try {
			parsed_ex = $1;
			YYACCEPT;
		} catch (std::exception &err) {
			parser_error = err.what();
			YYERROR;
		}
	}
	;

exp	: T_NUMBER		{$$ = $1;}
	| T_SYMBOL {
		if (is_lexer_symbol_predefined($1))
			$$ = $1.eval();
		else
			throw (std::runtime_error("unknown symbol '" + get_symbol_name($1) + "'"));
	}
	| T_LITERAL		{$$ = $1;}
	| T_DIGITS		{$$ = $1;}
	| T_SYMBOL '(' exprseq ')' {
		std::string n = get_symbol_name($1);
		if (n == "sqrt") {
			if ($3.nops() != 1)
				throw (std::runtime_error("too many arguments to sqrt()"));
			$$ = sqrt($3.op(0));
		} else if (n == "pow" || n == "power") {
		  if ($3.nops() != 2) 
			  throw std::invalid_argument("wrong number of arguments to pow()");
			$$ = power($3.op(0), $3.op(0));
		} else {
			unsigned i = function::find_function(n, $3.nops());
			$$ = function(i, ex_to<exprseq>($3)).eval(1);
		}
	}
	| exp T_EQUAL exp	{$$ = $1 == $3;}
	| exp T_NOTEQ exp	{$$ = $1 != $3;}
	| exp '<' exp		{$$ = $1 < $3;}
	| exp T_LESSEQ exp	{$$ = $1 <= $3;}
	| exp '>' exp		{$$ = $1 > $3;}
	| exp T_GREATEREQ exp	{$$ = $1 >= $3;}
	| exp '+' exp		{$$ = $1 + $3;}
	| exp '-' exp		{$$ = $1 - $3;}
	| exp '*' exp		{$$ = $1 * $3;}
	| exp '/' exp		{$$ = $1 / $3;}
	| '-' exp %prec NEG	{$$ = -$2;}
	| '+' exp %prec NEG	{$$ = $2;}
	| exp '^' exp		{$$ = pow($1, $3);}
	| exp '.' exp		{$$ = attach_index($1, $3, true);}
	| exp '~' exp		{$$ = attach_index($1, $3, false);}
	| exp '!'		{$$ = factorial($1);}
	| '(' exp ')'		{$$ = $2;}
	| '{' list_or_empty '}'	{$$ = $2;}
	| '[' matrix ']'	{$$ = lst_to_matrix(ex_to<lst>($2));}
	;

exprseq	: exp			{$$ = exprseq($1);}
	| exprseq ',' exp	{exprseq es(ex_to<exprseq>($1)); $$ = es.append($3);}
	;

list_or_empty: /* empty */	{$$ = *new lst;}
	| list			{$$ = $1;}
	;

list	: exp			{$$ = lst($1);}
	| list ',' exp		{lst l(ex_to<lst>($1)); $$ = l.append($3);}
	;

matrix	: '[' row ']'		{$$ = lst($2);}
	| matrix ',' '[' row ']' {lst l(ex_to<lst>($1)); $$ = l.append($4);}
	;

row	: exp			{$$ = lst($1);}
	| row ',' exp		{lst l(ex_to<lst>($1)); $$ = l.append($3);}
	;


/*
 *  Routines
 */

%%
// Attach index to expression
ex attach_index(const ex & base, ex i, bool covariant)
{
	// Toggle index variance if necessary
	if (is_a<varidx>(i)) {
		const varidx &vi = ex_to<varidx>(i);
		if (vi.is_covariant() != covariant)
			i = vi.toggle_variance();
	} else if (!covariant)
		throw (std::runtime_error("index '" + get_symbol_name(i) + "' is not a varidx and cannot be contravariant"));

	// Add index to an existing indexed object, or create a new indexed
	// object if there are no indices yet
	if (is_a<indexed>(base)) {
		const ex &b = base.op(0);
		exvector iv;
		for (unsigned n=1; n<base.nops(); n++)
			iv.push_back(base.op(n));
		iv.push_back(i);
		return indexed(b, iv);
	} else
		return indexed(base, i);
}

// Get last error encountered by parser
std::string get_parser_error(void)
{
	return parser_error;
}

} // namespace GiNaC

// Error print routine (store error string in parser_error)
int ginac_yyerror(char *s)
{
	GiNaC::parser_error = std::string(s) + " at " + std::string(ginac_yytext);
	return 0;
}
