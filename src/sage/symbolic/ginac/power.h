/** @file power.h
 *
 *  Interface to GiNaC's symbolic exponentiation (basis^exponent). */

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

#ifndef __GINAC_POWER_H__
#define __GINAC_POWER_H__

#include "basic.h"
#include "ex.h"

namespace GiNaC {

class numeric;
class add;
class mul;

/** This class holds a two-component object, a basis and and exponent
 *  representing exponentiation. */
class power : public basic
{
	GINAC_DECLARE_REGISTERED_CLASS(power, basic)
	
	friend class mul;
	
// member functions
	
	// other constructors
public:
	power(const ex & lh, const ex & rh) : inherited(&power::tinfo_static), basis(lh), exponent(rh) {}
	template<typename T> power(const ex & lh, const T & rh) : inherited(&power::tinfo_static), basis(lh), exponent(rh) {}
	
	// functions overriding virtual functions from base classes
public:
	unsigned precedence() const {return 60;}
	bool info(unsigned inf) const;
	size_t nops() const;
	ex op(size_t i) const;
	ex map(map_function & f) const;
	bool is_polynomial(const ex & var) const;
	int degree(const ex & s) const;
	int ldegree(const ex & s) const;
	ex coeff(const ex & s, int n = 1) const;
	ex eval(int level=0) const;
	ex evalf(int level=0) const;
	ex evalm() const;
	ex series(const relational & s, int order, unsigned options = 0) const;
	ex subs(const exmap & m, unsigned options = 0) const;
	bool has(const ex & other, unsigned options = 0) const;
	ex normal(exmap & repl, exmap & rev_lookup, int level = 0) const;
	ex to_rational(exmap & repl) const;
	ex to_polynomial(exmap & repl) const;
	ex conjugate() const;
	ex real_part() const;
	ex imag_part() const;
	int compare(const basic& other) const;
	int compare_symbol(const symbol& other) const;
protected:
	ex derivative(const symbol & s) const;
	ex eval_ncmul(const exvector & v) const;
	unsigned return_type() const;
	tinfo_t return_type_tinfo() const;
	ex expand(unsigned options = 0) const;
	unsigned calchash() const;

	// new virtual functions which can be overridden by derived classes
	// none
	
	// non-virtual functions in this class
protected:
	void print_power(const print_context & c, const char *powersymbol, const char *openbrace, const char *closebrace, unsigned level) const;
	void do_print_dflt(const print_dflt & c, unsigned level) const;
	void do_print_latex(const print_latex & c, unsigned level) const;
	void do_print_csrc(const print_csrc & c, unsigned level) const;
	void do_print_python(const print_python & c, unsigned level) const;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const;

	ex expand_add(const add & a, int n, unsigned options) const;
	ex expand_add_2(const add & a, unsigned options) const;
	ex expand_mul(const mul & m, const numeric & n, unsigned options, bool from_expand = false) const;
	
// member variables
	
protected:
	ex basis;
	ex exponent;
};

// wrapper functions

/** Symbolic exponentiation.  Returns a power-object as a new expression.
 *
 *  @param b the basis expression
 *  @param e the exponent expression */
inline ex pow(const ex & b, const ex & e)
{
	return power(b, e);
}
template<typename T1, typename T2>
inline ex pow(const T1 & b, const T2 & e)
{
	return power(ex(b), ex(e));
}

/** Square root expression.  Returns a power-object with exponent 1/2. */
inline ex sqrt(const ex & a)
{
	extern const ex _ex1_2;
	return power(a,_ex1_2);
}

} // namespace GiNaC

#endif // ndef __GINAC_POWER_H__
