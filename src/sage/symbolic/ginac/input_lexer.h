/** @file input_lexer.h
 *
 *  Lexical analyzer definition for reading expressions.
 *  This file must be processed with flex. */

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

#ifndef __GINAC_INPUT_LEXER_H__
#define __GINAC_INPUT_LEXER_H__

extern "C" {
#include <stdio.h>
}
	
#include "config.h"

// yacc stack type
#define YYSTYPE ex

// lex functions/variables
extern int ginac_yyerror(char *s);
extern int ginac_yylex(void);
extern void ginac_yyrestart(FILE *f);
extern char *ginac_yytext;

namespace GiNaC {

class ex;

/** Set the input string to be parsed by ginac_yyparse() (used internally). */
extern void set_lexer_string(const std::string &s);

/** Get name of symbol/index (used internally). */
extern std::string get_symbol_name(const ex & s);

/** Set the list of predefined symbols for the lexer (used internally). */
extern void set_lexer_symbols(ex l);

/** Check whether lexer symbol was predefined (vs. created by the lexer, e.g. function names). */
extern bool is_lexer_symbol_predefined(const ex &s);

/** The expression parser function (used internally). */
extern int ginac_yyparse(void);

/** The expression returned by the parser (used internally). */
extern ex parsed_ex;

/** Get error message from the parser. */
extern std::string get_parser_error(void);

} // namespace GiNaC

#endif // ndef __GINAC_INPUT_LEXER_H__
