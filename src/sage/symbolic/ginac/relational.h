/** @file relational.h
 *
 *  Interface to relations between expressions. */

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

#ifndef __GINAC_RELATIONAL_H__
#define __GINAC_RELATIONAL_H__

#include "basic.h"
#include "ex.h"

namespace GiNaC {

/** This class holds a relation consisting of two expressions and a logical
 *  relation between them. */
class relational : public basic
{
	GINAC_DECLARE_REGISTERED_CLASS(relational, basic)

// types
public:
	enum operators {
		equal,
		not_equal,
		less,
		less_or_equal,
		greater,
		greater_or_equal
	};
	
	enum result {
		True,
		False,
		undecidable,
		notimplemented
	};

	// other constructors
public:
	relational(ex  lhs, ex  rhs, operators oper=equal);
	
	// functions overriding virtual functions from base classes
public:
	unsigned precedence() const {return 20;}
	bool info(unsigned inf) const;
	size_t nops() const;
	ex op(size_t i) const;
	ex map(map_function & f) const;
	ex subs(const exmap & m, unsigned options = 0) const;
	ex eval(int level=0) const;

protected:
	ex eval_ncmul(const exvector & v) const;
	bool match_same_type(const basic & other) const;
	unsigned return_type() const;
	tinfo_t return_type_tinfo() const;
	int64_t calchash() const;

	// new virtual functions which can be overridden by derived classes
protected:
	void print_rel(const print_context & c, unsigned level, bool latex) const;
	void do_print_dflt(const print_context & c, unsigned level) const;
	void do_print_latex(const print_context & c, unsigned level) const;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const;

public:
	virtual ex lhs() const;
	virtual ex rhs() const;
	virtual operators the_operator() const;

	// non-virtual functions in this class
public:
	operator bool() const {
		return decide() == result::True;
	}
	result decide() const;

// member variables
	
protected:
	ex lh;
	ex rh;
	operators o;
};

} // namespace GiNaC

#endif // ndef __GINAC_RELATIONAL_H__
