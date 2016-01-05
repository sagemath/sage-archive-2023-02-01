/** @file add.cpp
 *
 *  Implementation of GiNaC's sums of expressions. */

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

#include "add.h"
#include "mul.h"
#include "archive.h"
#include "operators.h"
#include "matrix.h"
#include "utils.h"
#include "clifford.h"
#include "ncmul.h"
#include "constant.h"
#include "infinity.h"
#include "compiler.h"
#include "order.h"

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <string>

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(add, expairseq,
  print_func<print_context>(&add::do_print).
  print_func<print_latex>(&add::do_print_latex).
  print_func<print_csrc>(&add::do_print_csrc).
  print_func<print_tree>(&add::do_print_tree).
  print_func<print_python_repr>(&add::do_print_python_repr))

//////////
// default constructor
//////////

add::add()
{
	tinfo_key = &add::tinfo_static;
}

//////////
// other constructors
//////////

// public

add::add(const ex & lh, const ex & rh)
{
	tinfo_key = &add::tinfo_static;
	overall_coeff = _ex0;
	construct_from_2_ex(lh,rh);
	GINAC_ASSERT(is_canonical());
}

add::add(const exvector & v, bool do_hold)
{
	tinfo_key = &add::tinfo_static;
	overall_coeff = _ex0;
	construct_from_exvector(v, do_hold);
	GINAC_ASSERT(is_canonical());
}

add::add(const epvector & v)
{
	tinfo_key = &add::tinfo_static;
	overall_coeff = _ex0;
	construct_from_epvector(v);
	GINAC_ASSERT(is_canonical());
}

add::add(const epvector & v, const ex & oc)
{
	tinfo_key = &add::tinfo_static;
	overall_coeff = oc;
	construct_from_epvector(v);
	GINAC_ASSERT(is_canonical());
}

//add::add(std::unique_ptr<epvector> vp, const ex & oc)
//{
//	tinfo_key = &add::tinfo_static;
//	GINAC_ASSERT(vp.get()!=0);
//	overall_coeff = oc;
//	construct_from_epvector(*vp);
//	GINAC_ASSERT(is_canonical());
//}

//////////
// archiving
//////////

DEFAULT_ARCHIVING(add)

//////////
// functions overriding virtual functions from base classes
//////////

// public

void add::print_add(const print_context & c, unsigned level, bool latex) const
{
	if (precedence() <= level){
		if (latex)
			c.s << "{\\left(";
		else
			c.s << '(';
	}

	bool first = true;

	const epvector & sorted_seq = get_sorted_seq();
	// Then proceed with the remaining factors
	for (const auto & elem : sorted_seq) {
		std::stringstream tstream;
		std::unique_ptr<print_context> tcontext_p;
		if (latex) {
			tcontext_p.reset(new print_latex(tstream, c.options));
		} else {
			tcontext_p.reset(new print_dflt(tstream, c.options));
		}
		mul(elem.rest, elem.coeff).print(*tcontext_p, precedence());

		if (!first) {
			if (tstream.peek() == '-') {
				tstream.ignore();
				c.s << " - ";
			} else 
				c.s << " + ";
		} else {
			first = false;
		}
		tstream.get(*(c.s.rdbuf()));
	}

	// Finally print the "overall" numeric coefficient, if present.
	// This is just the constant coefficient. 
	if (!(ex_to<numeric>(overall_coeff)).is_zero()) {
		std::stringstream tstream;
		std::unique_ptr<print_context> tcontext_p;
		if (latex) {
			tcontext_p.reset(new print_latex(tstream, c.options));
		} else {
			tcontext_p.reset(new print_dflt(tstream, c.options));
		}
		overall_coeff.print(*tcontext_p, 0);
		if (!first) {
			if (tstream.peek() == '-') {
				c.s << " - ";
				tstream.ignore();
			} else
				c.s << " + ";
		}

		tstream.get(*(c.s.rdbuf()));
	}

	if (precedence() <= level) {
		if (latex)
			c.s << "\\right)}";
		else
			c.s << ')';
	}
}

void add::do_print(const print_context & c, unsigned level) const
{
	print_add(c, level, false);
}

void add::do_print_latex(const print_latex & c, unsigned level) const
{
	print_add(c, level, true);
}

void add::do_print_csrc(const print_csrc & c, unsigned level) const
{
	if (precedence() <= level)
		c.s << "(";
	
	// Print arguments, separated by "+" or "-"
	char separator = ' ';
        for (const auto & elem : seq) {
		
		// If the coefficient is negative, separator is "-"
		if (elem.coeff.is_equal(_ex_1) || 
			ex_to<numeric>(elem.coeff).numer().is_equal(*_num_1_p))
			separator = '-';
		c.s << separator;
		if (elem.coeff.is_equal(_ex1) || elem.coeff.is_equal(_ex_1)) {
			elem.rest.print(c, precedence());
		} else if (ex_to<numeric>(elem.coeff).numer().is_equal(*_num1_p) ||
				 ex_to<numeric>(elem.coeff).numer().is_equal(*_num_1_p))
		{
			elem.rest.print(c, precedence());
			c.s << '/';
			ex_to<numeric>(elem.coeff).denom().print(c, precedence());
		} else {
			elem.coeff.print(c, precedence());
			c.s << '*';
			elem.rest.print(c, precedence());
		}
		
		separator = '+';
	}
	
	if (!overall_coeff.is_zero()) {
		if (overall_coeff.info(info_flags::positive)
		 || is_a<print_csrc_cl_N>(c) || !overall_coeff.info(info_flags::real))  // sign inside ctor argument
			c.s << '+';
		overall_coeff.print(c, precedence());
	}
		
	if (precedence() <= level)
		c.s << ")";
}

void add::do_print_python_repr(const print_python_repr & c, unsigned) const
{
	c.s << class_name() << '(';
	op(0).print(c);
	for (size_t i=1; i<nops(); ++i) {
		c.s << ',';
		op(i).print(c);
	}
	c.s << ')';
}

bool add::info(unsigned inf) const
{
	switch (inf) {
		case info_flags::polynomial:
		case info_flags::integer_polynomial:
		case info_flags::cinteger_polynomial:
		case info_flags::rational_polynomial:
		case info_flags::real:
		case info_flags::rational:
		case info_flags::integer:
		case info_flags::crational:
		case info_flags::cinteger:
		case info_flags::positive:
		case info_flags::nonnegative:
		case info_flags::posint:
		case info_flags::nonnegint:
		case info_flags::even:
		case info_flags::crational_polynomial:
		case info_flags::rational_function: {
			for (const auto & elem : seq) {
				if (!(recombine_pair_to_ex(elem).info(inf)))
					return false;
			}
			if (overall_coeff.is_zero() && (inf == info_flags::positive || inf == info_flags::posint))
                                return true;
                        return overall_coeff.info(inf);
		}
		case info_flags::inexact:
		case info_flags::algebraic: {
                        if (overall_coeff.info(inf))
                                return true;
                        for (const auto & elem : seq) {
				if ((recombine_pair_to_ex(elem).info(inf)))
					return true;
			}
			return false;
		}
	}
	return inherited::info(inf);
}

bool add::is_polynomial(const ex & var) const
{
	for (const auto & elem : seq) {
		if (!(elem.rest).is_polynomial(var)) {
			return false;
		}
	}
	return true;
}

int add::degree(const ex & s) const
{
	int deg = std::numeric_limits<int>::min();
	if (!overall_coeff.is_zero())
		deg = 0;
	
	// Find maximum of degrees of individual terms
        for (const auto & elem : seq) {
		int cur_deg = elem.rest.degree(s);
		if (cur_deg > deg)
			deg = cur_deg;
	}
	return deg;
}

int add::ldegree(const ex & s) const
{
	int deg = std::numeric_limits<int>::max();
	if (!overall_coeff.is_zero())
		deg = 0;
	
	// Find minimum of degrees of individual terms
        for (const auto & elem : seq) {
		int cur_deg = elem.rest.ldegree(s);
		if (cur_deg < deg)
			deg = cur_deg;
	}
	return deg;
}

ex add::coeff(const ex & s, int n) const
{
	epvector coeffseq;
	epvector coeffseq_cliff;
	int rl = clifford_max_label(s);
	bool do_clifford = (rl != -1);
	bool nonscalar = false;

	// Calculate sum of coefficients in each term
        for (const auto & elem : seq) {
		ex restcoeff = elem.rest.coeff(s, n);
 		if (!restcoeff.is_zero()) {
 			if (do_clifford) {
 				if (clifford_max_label(restcoeff) == -1) {
 					coeffseq_cliff.push_back(combine_ex_with_coeff_to_pair(ncmul(restcoeff, dirac_ONE(rl)), elem.coeff));
				} else {
 					coeffseq_cliff.push_back(combine_ex_with_coeff_to_pair(restcoeff, elem.coeff));
					nonscalar = true;
 				}
			}
			coeffseq.push_back(combine_ex_with_coeff_to_pair(restcoeff, elem.coeff));
		}
	}

	return (new add(nonscalar ? coeffseq_cliff : coeffseq,
	                n==0 ? overall_coeff : _ex0))->setflag(status_flags::dynallocated);
}

/** Perform automatic term rewriting rules in this class.  In the following
 *  x stands for a symbolic variables of type ex and c stands for such
 *  an expression that contain a plain number.
 *  - +(;c) -> c
 *  - +(x;0) -> x
 *
 *  @param level cut-off in recursive evaluation */
ex add::eval(int level) const
{
	std::unique_ptr<epvector> evaled_seqp = evalchildren(level);
	if (unlikely(evaled_seqp.get() != nullptr)) {
		// do more evaluation later
		return (new add(*evaled_seqp, overall_coeff))->
		       setflag(status_flags::dynallocated);
	}
	
#ifdef DO_GINAC_ASSERT
        for (const auto & elem : seq) {
		GINAC_ASSERT(!is_exactly_a<add>(elem.rest));
		if (is_exactly_a<numeric>(elem.rest))
			dbgprint();
		GINAC_ASSERT(!is_exactly_a<numeric>(elem.rest));
	}
#endif // def DO_GINAC_ASSERT
	
	if (flags & status_flags::evaluated) {
		GINAC_ASSERT(seq.size()>0);
		GINAC_ASSERT(seq.size()>1 || !overall_coeff.is_zero());
		return *this;
	}
		
	// handle infinity
        for (auto i = seq.begin(); i != seq.end(); i++)
                if (unlikely(is_exactly_a<infinity>(i->rest)))
                        return eval_infinity(i);

	/** Perform automatic term rewriting rules */
	int seq_size = seq.size();
	if (seq_size == 0) {
		// +(;c) -> c
		return overall_coeff;
	} else if (seq_size == 1 && overall_coeff.is_zero()) {
		// +(x;0) -> x
		return recombine_pair_to_ex(*(seq.begin()));
	} else if (!overall_coeff.is_zero() && seq[0].rest.return_type() != return_types::commutative) {
		throw (std::logic_error("add::eval(): sum of non-commutative objects has non-zero numeric term"));
	}

	// if any terms in the sum still are purely numeric, then they are more
	// appropriately collected into the overall coefficient
	int terms_to_collect = 0;
	for (const auto & elem : seq)
		if (unlikely(is_exactly_a<numeric>(elem.rest)))
			++terms_to_collect;
	if (terms_to_collect) {
                epvector s;
		s.reserve(seq_size - terms_to_collect);
		numeric oc = *_num0_p;
                for (const auto & elem : seq)
			if (unlikely(is_exactly_a<numeric>(elem.rest)))
				oc = oc.add((ex_to<numeric>(elem.rest)).mul(ex_to<numeric>(elem.coeff)));
			else
				s.push_back(elem);
		return (new add(s, ex_to<numeric>(overall_coeff).add_dyn(oc)))
		        ->setflag(status_flags::dynallocated);
	}

	return this->hold();
}



namespace { // anonymous namespace
	infinity infinity_from_iter(epvector::const_iterator i)
	{
		GINAC_ASSERT(is_exactly_a<infinity>(i->rest));
		GINAC_ASSERT(is_a<numeric>(i->coeff));
		infinity result = ex_to<infinity>(i->rest);
		result *= i->coeff;
		return result;
	}
} // end anonymous namespace


ex add::eval_infinity(epvector::const_iterator infinity_iter) const
{
	GINAC_ASSERT(is_exactly_a<infinity>(infinity_iter->rest));
	infinity result = infinity_from_iter(infinity_iter);

        for (auto i = seq.begin(); i != seq.end(); i++) {
                if (not is_exactly_a<infinity>(i->rest)) continue;
                if (i == infinity_iter) continue;
		infinity i_infty = infinity_from_iter(i);
		result += i_infty;
        }
	return result;
}


ex add::evalm() const
{
	// Evaluate children first and add up all matrices. Stop if there's one
	// term that is not a matrix.
	epvector s;
	s.reserve(seq.size());

	bool all_matrices = true;
	bool first_term = true;
	matrix sum;

	for (const auto & elem : seq) {
		const ex &m = recombine_pair_to_ex(elem).evalm();
		s.push_back(split_ex_to_pair(m));
		if (is_exactly_a<matrix>(m)) {
			if (first_term) {
				sum = ex_to<matrix>(m);
				first_term = false;
			} else
				sum = sum.add(ex_to<matrix>(m));
		} else
			all_matrices = false;
	}

	if (all_matrices)
		return sum + overall_coeff;
	else
		return (new add(s, overall_coeff))->setflag(status_flags::dynallocated);
}

ex add::conjugate() const
{
	std::unique_ptr<exvector> v(nullptr);
	for (size_t i=0; i<nops(); ++i) {
		if (v) {
			v->push_back(op(i).conjugate());
			continue;
		}
		ex term = op(i);
		ex ccterm = term.conjugate();
		if (are_ex_trivially_equal(term, ccterm))
			continue;
		v.reset(new exvector);
		v->reserve(nops());
		for (size_t j=0; j<i; ++j)
			v->push_back(op(j));
		v->push_back(ccterm);
	}
	if (v) {
		ex result = add(*v);
		return result;
	}
	return *this;
}

ex add::real_part() const
{
	epvector v;
	v.reserve(seq.size());
	for (const auto & elem : seq)
		if ((elem.coeff).info(info_flags::real)) {
			ex rp = (elem.rest).real_part();
			if (!rp.is_zero())
				v.push_back(expair(rp, elem.coeff));
		} else {
			ex rp=recombine_pair_to_ex(elem).real_part();
			if (!rp.is_zero())
				v.push_back(split_ex_to_pair(rp));
		}
	return (new add(v, overall_coeff.real_part()))
		-> setflag(status_flags::dynallocated);
}

ex add::imag_part() const
{
	epvector v;
	v.reserve(seq.size());
	for (const auto & elem : seq)
		if ((elem.coeff).info(info_flags::real)) {
			ex ip = (elem.rest).imag_part();
			if (!ip.is_zero())
				v.push_back(expair(ip, elem.coeff));
		} else {
			ex ip=recombine_pair_to_ex(elem).imag_part();
			if (!ip.is_zero())
				v.push_back(split_ex_to_pair(ip));
		}
	return (new add(v, overall_coeff.imag_part()))
		-> setflag(status_flags::dynallocated);
}

ex add::eval_ncmul(const exvector & v) const
{
	if (seq.empty())
		return inherited::eval_ncmul(v);
	else
		return seq.begin()->rest.eval_ncmul(v);
}    

// protected

/** Implementation of ex::diff() for a sum. It differentiates each term.
 *  @see ex::diff */
ex add::derivative(const symbol & y) const
{
	epvector s;
	s.reserve(seq.size());
	
	// Only differentiate the "rest" parts of the expairs. This is faster
	// than the default implementation in basic::derivative() although
	// if performs the same function (differentiate each term).
	for (const auto & elem : seq)
		s.push_back(combine_ex_with_coeff_to_pair(elem.rest.diff(y), elem.coeff));
	return (new add(s, _ex0))->setflag(status_flags::dynallocated);
}

int add::compare_same_type(const basic & other) const
{
	return inherited::compare_same_type(other);
}

unsigned add::return_type() const
{
	if (seq.empty())
		return return_types::commutative;
	else
		return seq.begin()->rest.return_type();
}

tinfo_t add::return_type_tinfo() const
{
	if (seq.empty())
		return this;
	else
		return seq.begin()->rest.return_type_tinfo();
}

// Note: do_index_renaming is ignored because it makes no sense for an add.
ex add::thisexpairseq(const epvector & v, const ex & oc, bool) const
{
	return (new add(v,oc))->setflag(status_flags::dynallocated);
}

//// Note: do_index_renaming is ignored because it makes no sense for an add.
ex add::thisexpairseq(std::unique_ptr<epvector> vp, const ex & oc, bool) const
{
	return (new add(*vp,oc))->setflag(status_flags::dynallocated);
}

expair add::split_ex_to_pair(const ex & e) const
{
	if (is_exactly_a<mul>(e)) {
		const mul &mulref(ex_to<mul>(e));
		const ex &numfactor = mulref.overall_coeff;
		auto mulcopyp = new mul(mulref);
		mulcopyp->overall_coeff = _ex1;
		mulcopyp->clearflag(status_flags::evaluated);
		mulcopyp->clearflag(status_flags::hash_calculated);
		mulcopyp->setflag(status_flags::dynallocated);
		return expair(*mulcopyp,numfactor);
	}
	return expair(e,_ex1);
}

expair add::combine_ex_with_coeff_to_pair(const ex & e,
										  const ex & c) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(c));
	if (is_exactly_a<mul>(e)) {
		const mul &mulref(ex_to<mul>(e));
		const ex &numfactor = mulref.overall_coeff;
		auto mulcopyp = new mul(mulref);
		mulcopyp->overall_coeff = _ex1;
		mulcopyp->clearflag(status_flags::evaluated);
		mulcopyp->clearflag(status_flags::hash_calculated);
		mulcopyp->setflag(status_flags::dynallocated);
		if (c.is_equal(_ex1))
			return expair(*mulcopyp, numfactor);
		else if (numfactor.is_equal(_ex1))
			return expair(*mulcopyp, c);
		else
			return expair(*mulcopyp, ex_to<numeric>(numfactor).mul_dyn(ex_to<numeric>(c)));
	} else if (is_exactly_a<numeric>(e)) {
		if (c.is_equal(_ex1))
			return expair(e, _ex1);
		return expair(ex_to<numeric>(e).mul_dyn(ex_to<numeric>(c)), _ex1);
	}
	return expair(e, c);
}

expair add::combine_pair_with_coeff_to_pair(const expair & p,
											const ex & c) const
{
	GINAC_ASSERT(is_exactly_a<numeric>(p.coeff));
	GINAC_ASSERT(is_exactly_a<numeric>(c));

	if (is_exactly_a<numeric>(p.rest)) {
		GINAC_ASSERT(ex_to<numeric>(p.coeff).is_equal(*_num1_p)); // should be normalized
		return expair(ex_to<numeric>(p.rest).mul_dyn(ex_to<numeric>(c)),_ex1);
	}

	return expair(p.rest,ex_to<numeric>(p.coeff).mul_dyn(ex_to<numeric>(c)));
}

ex add::recombine_pair_to_ex(const expair & p) const
{
	if (ex_to<numeric>(p.coeff).is_equal(*_num1_p))
		return p.rest;
	else
		return (new mul(p.rest,p.coeff))->setflag(status_flags::dynallocated);
}

ex add::expand(unsigned options) const
{
	std::unique_ptr<epvector> vp = expandchildren(options);
	if (vp.get() == nullptr) {
		// the terms have not changed, so it is safe to declare this expanded
		return (options == 0) ? setflag(status_flags::expanded) : *this;
	}

	return (new add(*vp, overall_coeff))->setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0));
}

const epvector & add::get_sorted_seq() const
{
	if (seq_sorted.empty() and not seq.empty()) {
		seq_sorted = epvector(seq.size());
		partial_sort_copy(seq.begin(), seq.end(),
				  seq_sorted.begin(), seq_sorted.end(),
				  print_order_pair());
	}
	return expairseq::get_sorted_seq();
}

ex add::lead_coeff() const {
	// if sorted_seq was computed before, we don't have to search
	if (seq_sorted.empty() and not seq.empty()) {
		return min_element(seq.begin(), seq.end(),
				print_order_pair())->coeff;
	} else {
		return seq_sorted.begin()->coeff;
	}
}

} // namespace GiNaC
