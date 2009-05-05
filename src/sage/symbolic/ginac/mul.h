/** @file mul.h
 *
 *  Interface to GiNaC's products of expressions. */

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

#ifndef __GINAC_MUL_H__
#define __GINAC_MUL_H__

#include "expairseq.h"
#include "power.h"

namespace GiNaC {

/** Product of expressions. */
class mul : public expairseq
{
	GINAC_DECLARE_REGISTERED_CLASS(mul, expairseq)
	
	friend class add;
	friend class ncmul;
	friend class power;
	
	// other constructors
public:
	mul(const ex & lh, const ex & rh);
	mul(const exvector & v);
	mul(const epvector & v);
	mul(const epvector & v, const ex & oc, bool do_index_renaming = false);
	mul(std::auto_ptr<epvector> vp, const ex & oc, bool do_index_renaming = false);
	mul(const ex & lh, const ex & mh, const ex & rh);
	
	// functions overriding virtual functions from base classes
public:
	unsigned precedence() const {return 50;}
	bool info(unsigned inf) const;
	int degree(const ex & s) const;
	int ldegree(const ex & s) const;
	ex coeff(const ex & s, int n = 1) const;
	bool has(const ex & other, unsigned options = 0) const;
	ex eval(int level=0) const;
	ex evalf(int level=0, int prec=0) const;
	ex real_part() const;
	ex imag_part() const;
	ex evalm() const;
	ex series(const relational & s, int order, unsigned options = 0) const;
	ex normal(exmap & repl, exmap & rev_lookup, int level = 0) const;
	numeric integer_content() const;
	ex smod(const numeric &xi) const;
	numeric max_coefficient() const;
	exvector get_free_indices() const;
protected:
	ex derivative(const symbol & s) const;
	ex eval_ncmul(const exvector & v) const;
	unsigned return_type() const;
	tinfo_t return_type_tinfo() const;
	ex thisexpairseq(const epvector & v, const ex & oc, bool do_index_renaming = false) const;
	ex thisexpairseq(std::auto_ptr<epvector> vp, const ex & oc, bool do_index_renaming = false) const;
	expair split_ex_to_pair(const ex & e) const;
	expair combine_ex_with_coeff_to_pair(const ex & e, const ex & c) const;
	expair combine_pair_with_coeff_to_pair(const expair & p, const ex & c) const;
	ex recombine_pair_to_ex(const expair & p) const;
	bool expair_needs_further_processing(epp it);
	ex default_overall_coeff() const;
	void combine_overall_coeff(const ex & c);
	void combine_overall_coeff(const ex & c1, const ex & c2);
	bool can_make_flat(const expair & p) const;
	ex expand(unsigned options=0) const;
	void find_real_imag(ex&, ex&) const;
	int compare(const basic& other) const;
	
	// new virtual functions which can be overridden by derived classes
	// none
	
	// non-virtual functions in this class
public:
	ex algebraic_subs_mul(const exmap & m, unsigned options) const;
	double total_degree() const;
	int compare_symbol(const symbol &other) const;
	int compare_pow(const power &other) const;
protected:
	void print_overall_coeff(const print_context & c,
			const char *mul_sym, bool latex=false) const;
	void print_exvector(const exvector & v, const print_context & c,
		const char* sep) const;
	void do_print(const print_context & c, unsigned level) const;
	void do_print_latex(const print_latex & c, unsigned level) const;
	void do_print_rat_func(const print_context & c, unsigned level, 
			bool latex_tags) const;
	void do_print_csrc(const print_csrc & c, unsigned level) const;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const;
	static bool can_be_further_expanded(const ex & e);
	std::auto_ptr<epvector> expandchildren(unsigned options) const;

	mutable double tdegree;
};

} // namespace GiNaC

#endif // ndef __GINAC_MUL_H__
